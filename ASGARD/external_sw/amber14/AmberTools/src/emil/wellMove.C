#include "config.h"
#ifdef USE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <sstream>

#include <math.h>
#include <stdlib.h>
#include <string.h>

using namespace std;
#include "emil.h"

const double Pi = M_PI;

void hsc::tryWellMove(){
  
  int       lIndex, atomsInChain, i;
  Particle *p;
  Well     *w, *ww;
  Cell     *c, *newC;
  double    Rold[3];
  int       newCellW;
  double    oldE, newE;
  double    newDist, oldDist;
  double    balanceFactor, log_pAccept, invxCellW[3];

  //choose a random liquid chain
  p  =  pickLiquidChain( &lIndex, &atomsInChain );
  w  =  p->theWell;  

  // save old configuration
  oldDist = p->rijsqW;
  for ( i=0; i<3; i++){
       Rold[i]  = w->R[i];
  }
  
  //get old energy wrt wells, using saved distance
  oldE    = p->theWell->computePhiOfDist(oldDist);
  
  // move into sphere if outside it
  if (p->rijsqW > rcut2Liquid ) {
    balanceFactor =  ( 4.0 * M_PI * rcut3Liquid / 3.0 ) / box->V;
    relocateWell( w );
  }
  else {  //move out of the sphere if inside it
    balanceFactor = 1.0 - ( 4.0 * M_PI * rcut3Liquid / 3.0 ) / box->V;
    do{
      w->R[0] = (mRan() - 0.5)*box->x[0];
      w->R[1] = (mRan() - 0.5)*box->x[1];
      w->R[2] = (mRan() - 0.5)*box->x[2];
    }while( computeDistance( *p ) <= rcut2Liquid ); 
  }
  
  // compute new distance from wells and calculate energy in passing
  newDist = computeDistance( *p );
  newE    = p->theWell->computePhiOfDist(newDist);
  
  balanceFactor = log(balanceFactor);
  log_pAccept   = balanceFactor + mixFunc_abs( restLambda ) * (oldE- newE);
  
  if (log(mRan()) < log_pAccept ){
    
    //need this to assign cell positions
    for(i = 0; i < 3; i++){
      //invxCell[i]  = 1.0 / box->xCell[i];
      invxCellW[i] = 1.0 / box->xCellW[i];
    }
  
    //find new cell
    newCellW   =  selectCellW( w->R, invxCellW );
    if (newCellW < box->nCellsW)
        newC = &cellsW[newCellW];
    else 
        newC = &cellsW[newCellW-1];
        
    if( testWellOverlap( w, newC ) ){
        //don't drop wells into overlap with eachother... return rejected.
        w->R[0] = Rold[0];
        w->R[1] = Rold[1];
        w->R[2] = Rold[2];
        return;
    }
        
    //Clear to Accept://///////////
    //update the old cell's well list.  Could use a doubly-linked list data structure, but not important.
    c  = w->cell;
    ww = c->firstWell;
    if( ww == w ){
      c->firstWell = w->next;
    }else{
      while( ww ){
        if( ww->next == w ){
          ww->next = w->next;
          break;
        }
        ww = ww->next;
      }
    }
    
    //update the list of wells in the new cell
    w->cell         = newC;
    w->next         = newC->firstWell;
    newC->firstWell = w;
        
    //formally move the well.
    relocAccept++;
    p->rijsqW = newDist;
    return;
  }
  else { //failed MC acceptance test
    w->R[0] = Rold[0];
    w->R[1] = Rold[1];
    w->R[2] = Rold[2];
  }
  
  return;
}

//check that the well is not in "overlap" with any other wells, solid or liquid.
bool hsc::testWellOverlap( Well *w, Cell *c ){
  Well   *ww;
  Cell   *cc;
  double  r2;
  
  ww = c->firstWell;
  while( ww ){
    if( ww != w ){
      r2  = ( w->R[0] - ww->R[0] ) * ( w->R[0] - ww->R[0] );
      r2 += ( w->R[1] - ww->R[1] ) * ( w->R[1] - ww->R[1] );
      r2 += ( w->R[2] - ww->R[2] ) * ( w->R[2] - ww->R[2] );
    
      if( sqrt(r2) < w->req + ww->req ){
        return( true );
      }
    }
    ww = ww->next;
  }
  for(int j=0; j<26; j++){
      cc = c->neighbours[j];
      ww = cc->firstWell;
      while( ww ){
        
        if( ww != w ){
          r2  = ( w->R[0] - ww->R[0] ) * ( w->R[0] - ww->R[0] );
          r2 += ( w->R[1] - ww->R[1] ) * ( w->R[1] - ww->R[1] );
          r2 += ( w->R[2] - ww->R[2] ) * ( w->R[2] - ww->R[2] );
        
          if(  sqrt(r2) < w->req + ww->req ){
            return( true );
          }
        }
        
        ww = ww->next;
      }
  }
  return( false );
}


void hsc::relocateWell( Well *w ){

  Particle *p;
  bool      outside;
  double    A[3];
  
  p = w->theParticle;
  
  //find a point in the sphere r = sqrt(rcut2)
  outside = true;
  do {
        for (int i = 0; i < 3; i++) 
          A[i] = (mRan() - 0.5)*dCutLiquid;
        if ((A[0]*A[0] + A[1]*A[1] + A[2]*A[2]) < rcut2Liquid )
          outside = false;
  }while(outside);

  for (int i = 0; i < 3; i++){
      w->R[i] = p->R[i] + A[i];
  }
  
}


double Well::computePhi(Box *b) {

    double Phi, *R;

    R = theParticle->R;
    Phi = 0.0;
    for(int i = 0; i < num_bonds; i++){
      double dx, dR, dR2;

      dR2 = 0.0;
      for( int j = 0; j < 3; j++){
        dx = bondedWells[i]->theParticle->R[j] - R[j];
        while( fabs(dx) > b->halfx[j] )
          dx -= copysign(b->x[j], dx);
        dR2 += dx * dx;
      }
      dR   = sqrt(dR2);
      dR  -= bond_eq_l[i];
      Phi += 0.5 * dR * dR * bond_k;
    }
    return( Phi );
}

double Well::compute_dPhi(double *dPhi_dr, Box *b) {

    double Phi, *r, *R, invDr;
    double dx[3], dR2, dR, delta;

    r   = theParticle->R;
    Phi = 0.0;
    dR2 = 0.0;

    dPhi_dr[0] = 0.0;
    dPhi_dr[1] = 0.0;
    dPhi_dr[2] = 0.0;

    for(int neigh = 0; neigh < num_bonds; neigh++){
      R   = bondedWells[neigh]->theParticle->R;
      for( int crd = 0; crd < 3; crd++){
        dx[crd] = R[crd] - r[crd];
        while( fabs(dx[crd]) > b->halfx[crd] )
          dx[crd] -= copysign(b->x[crd], dx[crd]);
      }
      dR2   = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];
      dR    = sqrt(dR2);
      delta = dR - bond_eq_l[neigh];
      invDr = 1.0 / dR;

      dPhi_dr[0] += bond_k * delta * dx[0] * invDr;
      dPhi_dr[1] += bond_k * delta * dx[1] * invDr;
      dPhi_dr[2] += bond_k * delta * dx[2] * invDr;
    }
    return( Phi );
  }


