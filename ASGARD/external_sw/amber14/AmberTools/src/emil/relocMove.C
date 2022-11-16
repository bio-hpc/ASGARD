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


double hsc::nonHastingsReloc(){//we can only do this non-Hastings reloc move when the amber hamiltonian is mixed out to 0.
  
  int       lIndex, atomsInChain, i, j, ii, jj;
  Particle *p, *q;
  bool      overlap = false;
  double   *Rold;
  int       newCell, newCellW;
  int      *oldCell, *oldCellW;
  double    oldE, newE;
  double    newDist, *oldDist;
  double    balanceFactor, log_pAccept, invxCell[3], invxCellW[3];
  double    translate[3];
  double    deSoft;

  //choose a random liquid chain
  p  =  pickLiquidChain( &lIndex, &atomsInChain );
  
  //need this to assign cell positions
  for(i = 0; i < 3; i++){
     invxCell[i]  = 1.0 / box->xCell[i];
     invxCellW[i] = 1.0 / box->xCellW[i];
  }
  
  // save old configuration
  Rold     = new double[atomsInChain * 3];
  oldDist  = new double[atomsInChain];
  oldCell  = new int[atomsInChain];
  oldCellW = new int[atomsInChain];
  
  q  = p;
  jj = 0;
  j  = 0;
  do{
    for ( ii=0; ii<3; ii++){
       Rold[jj++]  = q->R[ii];
    }
    oldCell[j]   = q->cell->myId;
    oldCellW[j]  = q->WCell->myId;
    oldDist[j]   = q->rijsqW;
    q = q->strandAtom;
    j = j + 1;
  }while( q );
  
  deSoft    = getSoftEnergyAll( p );

  // move into sphere if outside it
  if (p->rijsqW > rcut2Liquid ) {
    balanceFactor =  ( 4.0 * M_PI * rcut3Liquid / 3.0 ) / box->V;
    relocateChain( p );
  }
  else {  //move out of the sphere if inside it
    balanceFactor = 1.0 - ( 4.0 * M_PI * rcut3Liquid / 3.0 ) / box->V;
    do{
      p->R[0] = (mRan() - 0.5)*box->x[0];
      p->R[1] = (mRan() - 0.5)*box->x[1];
      p->R[2] = (mRan() - 0.5)*box->x[2];
    }while( computeDistance( *p ) <= rcut2Liquid ); 
      
    translate[0] = p->R[0] - Rold[0];
    translate[1] = p->R[1] - Rold[1];
    translate[2] = p->R[2] - Rold[2];
  
    q = p->strandAtom;
    while( q ){
      q->R[0] = q->R[0] + translate[0];//doesn't matter if the qs go outside the box.
      q->R[1] = q->R[1] + translate[1];
      q->R[2] = q->R[2] + translate[2];
      q       = q->strandAtom;
    }
  }

  //update pbc cell.
  q  = p;
  do {
    newCell = selectCell( q->R, invxCell );
    if ( &cells[newCell] != q->cell ){
      q->moveBetweenCells( *(q->cell), cells[newCell]);
    }
    q = q->strandAtom;
  } while( q );
  
  //for stability, we avoid dropping it within ~1 Angstrom of another particle...if the particles are going 
  //to move into eachother, the dynamics can handle this.
  //overlap = testChainOverlap( p );
  overlap = false;
  if (!overlap){
    
    //get old chain energy wrt wells, using vector of saved distances
    oldE    = p->theWell->computePhiOfDist(oldDist[0]);
  
    // compute new distance from wells and calculate energy in passing
    newDist = computeDistance( *p );
    newE    = p->theWell->computePhiOfDist(newDist);

    
    balanceFactor = log(balanceFactor);
    
    deSoft       -= getSoftEnergyAll( p );
    log_pAccept   =   balanceFactor  + mixFunc_abs( restLambda ) * (oldE- newE) 
                                     + mixFunc_soft( restLambda ) * deSoft;
   
    if (log(mRan()) < log_pAccept ){

    
      // move also with respect to well-list
      newCellW = selectCellW( p->R, invxCellW );
      if( newCellW != oldCellW[0] ) {
          p->moveBetweenCellsW( cellsW[oldCellW[0]], cellsW[newCellW]);
      }
      p->rijsqW = computeDistance( *p );
     
      
      return( 1 );
    }
    else { //failed MC acceptance test
      replaceTrialReloc( p, invxCell, Rold, oldCell );
    }
  }
  else { //trial move had overlap
      replaceTrialReloc( p, invxCell, Rold, oldCell );
  }
  
  delete[] Rold;
  delete[] oldDist;
  delete[] oldCell;
  delete[] oldCellW;
  return( 0 );
}




double hsc::tryRelocMove(){
  
  int       lIndex, atomsInChain, i, j, ii, jj;
  Particle *p, *q;
  bool      overlap = false;
  double   *Rold;
  int       newCell, newCellW;
  int      *oldCell, *oldCellW;
  double    oldE, newE, oldEnb, newEnb, oldElec, newElec;
  double    newDist, *oldDist;
  double    balanceFactor, log_pAccept, invxCell[3], invxCellW[3];
  double    translate[3];
  double    deSoft;

  //choose a random liquid chain
  p  =  pickLiquidChain( &lIndex, &atomsInChain );
  
  
  //for balance, if we do not drop particles into 'overlap'
  //then we shouldn't pick them up from overlap either
  if( testChainOverlap( p ) ){
    return( 0.0 );
  }
  
  //need this to assign cell positions
  for(i = 0; i < 3; i++){
     invxCell[i]  = 1.0 / box->xCell[i];
     invxCellW[i] = 1.0 / box->xCellW[i];
  }
  
  // save old configuration
  Rold     = new double[atomsInChain * 3];
  oldDist  = new double[atomsInChain];
  oldCell  = new int[atomsInChain];
  oldCellW = new int[atomsInChain];
  
  q  = p;
  jj = 0;
  j  = 0;
  do{//don't need this loop I think.  Come back later and get rid of it when sure.
    for ( ii=0; ii<3; ii++){
       Rold[jj++]  = q->R[ii];
    }
    oldCell[j]   = q->cell->myId;
    oldCellW[j]  = q->WCell->myId;
    oldDist[j]   = q->rijsqW;
    q = q->strandAtom;
    j = j + 1;
  }while( q );
  
  //save intermolecular energy
  oldEnb  = chainEnb( p, &oldElec ); //multiplied by 1/kT because up to this point we are using units from the external program

  deSoft  = getSoftEnergyAll( p );//A possible optimisation is to calculate this later, not now.

  //temporary: hack in a system-wide ee evaluation
  oldEnb -= oldElec;
  oldElec = wholeSystemEE( p );
  oldEnb += oldElec;
  
  // move into sphere if outside it
  if (p->rijsqW > rcut2Liquid ) {
    balanceFactor =  ( 4.0 * M_PI * rcut3Liquid / 3.0 ) / box->V;
    relocateChain( p );
  }
  else {  //move out of the sphere if inside it
    balanceFactor = 1.0 - ( 4.0 * M_PI * rcut3Liquid / 3.0 ) / box->V;
    do{
      p->R[0] = (mRan() - 0.5)*box->x[0];
      p->R[1] = (mRan() - 0.5)*box->x[1];
      p->R[2] = (mRan() - 0.5)*box->x[2];
    }while( computeDistance( *p ) <= rcut2Liquid ); 
      
    translate[0] = p->R[0] - Rold[0];
    translate[1] = p->R[1] - Rold[1];
    translate[2] = p->R[2] - Rold[2];
  
    q = p->strandAtom;
    while( q ){
      q->R[0] = q->R[0] + translate[0];//doesn't matter if the qs go outside the box.
      q->R[1] = q->R[1] + translate[1];
      q->R[2] = q->R[2] + translate[2];
      q       = q->strandAtom;
    }
  }

  //update pbc cell.
  q  = p;
  do {
    newCell = selectCell( q->R, invxCell );
    if ( &cells[newCell] != q->cell ){
      q->moveBetweenCells( *(q->cell), cells[newCell]);
    }
    q = q->strandAtom;
  } while( q );
  
  //for stability, we avoid dropping it within ~1 Angstrom of another particle...if the particles are going 
  //to move into eachother, the dynamics can handle this.
  overlap = testChainOverlap( p );
  //overlap = false;
  if (!overlap){
    
    //get old chain energy wrt wells, using vector of saved distances
    //oldE     = chainEnergy( p, oldDist );
    oldE     = p->theWell->computePhiOfDist(oldDist[0]);
    
    //get new chain vdw energy
    newEnb   =  chainEnb( p, &newElec ); 
    
    //temporary: hack in a system-wide ee evaluation
    newEnb  -= newElec;
    newElec  = wholeSystemEE( p );
    newEnb  += newElec;
  
    // compute new distance from wells and calculate energy in passing
    newDist = computeDistance( *p );
    newE    = p->theWell->computePhiOfDist(newDist);
    deSoft -= getSoftEnergyAll( p );
    
    balanceFactor = log(balanceFactor);
    log_pAccept   = balanceFactor + mixFunc_abs( restLambda ) * (oldE- newE) 
                                  + mixFunc_mol( 1.0 - amberLambda ) * beta * (oldEnb - newEnb) 
                                  + mixFunc_soft( restLambda ) * deSoft;
    
  
    if (log(mRan()) < log_pAccept ){

    
     newCellW = selectCellW( p->R, invxCellW );
     if( newCellW != oldCellW[0] ) {
          p->moveBetweenCellsW( cellsW[oldCellW[0]], cellsW[newCellW]);
     }
     p->rijsqW = computeDistance( *p );
  
      
      //save the deltaE hastings
      deHastings_estimated = newEnb - oldEnb;
      deAbstract           = (newE - oldE);
      deVdwHastings        = (newEnb - newElec) - (oldEnb - oldElec);
      deElecHastings       = newElec - oldElec;
      
     *logFile << "logging hastings move, de_abs: " << (newE - oldE) * mixFunc_abs(restLambda)  << endl;
      
     if( myTaskId == 0 ) {cerr << "logging hastings move, de_abs" << endl;}
     
      //save the system state before the move, in case it is rejected by AMBER.
      saveRelocAttempt( p, Rold, log_pAccept );
      
      delete[] Rold;
      delete[] oldDist;
      delete[] oldCell;
      delete[] oldCellW;
      
      return( log_pAccept - balanceFactor );
    }
    else { //failed MC acceptance test
      replaceTrialReloc( p, invxCell, Rold, oldCell );
    }
  }
  else { //trial move had overlap
      replaceTrialReloc( p, invxCell, Rold, oldCell );
  }
  
  delete[] Rold;
  delete[] oldDist;
  delete[] oldCell;
  delete[] oldCellW;
  return( 0 );
}



