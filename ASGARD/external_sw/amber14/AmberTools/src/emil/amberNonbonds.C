#include "config.h"
#ifdef USE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "emil.h"
using namespace std;

//Most of this code is retired for the moment... it was 
//supposed to provide for Hastings Metropolis moves

//find the softcoring onset radius
double hsc::sc_ljA(int pairType){
  
    double a;

      a = LJ_aCoeff[pairType];

    return(a);
    
}
  
//find the softcoring onset radius
double hsc::scEpsilon(int pairType){
  
    double epsilon;

  //load the epsilon used in softcoring the interaction between this atom pair
    if( LJ_bCoeff[pairType] > 0.0 ){
      //define epsilon as the zero-point of the LJ-force
      epsilon  = pow( 2.0 * LJ_aCoeff[pairType] /  LJ_bCoeff[pairType], 0.16666666 );
    }else{
       epsilon  = 0.5;
    }
    
    return(epsilon);
    
}
  
double hsc::sigmoid(double r, double intercept, double intersect){
 
    double gradient;
  
    if( r > intersect * amberLambda ) {
        return( r );
    }
    intercept = intersect * 0.99;

    //blend a quadradic with y intercept "intercept" to merge at "intersect"
    gradient = (intersect - intercept) /  ( intersect * intersect  * amberLambda);
    
    return( intercept + r * r * gradient ); 
}


// loop over all N for the electrostatic energy between particles.
double hsc::wholeSystemEE( Particle *p ){
  Particle *q;
  double    ee;
#ifdef USE_MPI
  double    eeTmp;
#endif
  int       qIndex, pId;
  double    rij[3], rijsq, r, rprime;
  int       pair, pairIndex;
  double    epsilon;
  
  ee  = 0.0;
  
  
  do{
  
    pId = p->myId;
  
    for(int ii = 0; ii < myNatoms; ii++ ){
      qIndex = myAtomIds[ii];
      
      q      = &part[qIndex];
      if( q->chainIndex == p->chainIndex )
        continue;
    
      //load the epsilon used in softcoring the interaction between this atom pair
      pairIndex = (LJ_atomTypes[p->myId] - 1) * nLJ_atomTypes + LJ_atomTypes[q->myId] - 1;
      pair      = LJ_pairTypes[pairIndex];
      epsilon   = softcore_scale * scEpsilon( pair );
      
      //measure particle-particle distance
      for (int j = 0; j < 3; j++){
        rij[j] = p->R[j] - q->R[j];
	if( periodic )
	  while (fabs(rij[j]) > box->halfx[j]) 
	    rij[j] -= copysign(box->x[j], rij[j]);
      } 
      rijsq   = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
      r         = sqrt(rijsq);
      
      
      //softcoring:
      if( r < epsilon ){
          //rprime = r + amberLambda * epsilon * exp(-1.0 * r / epsilon);
          rprime = sigmoid( r, 0.8 * epsilon, epsilon );
      }else{
          rprime = r;
      }
        //coulomb energy
      ee       +=  ( atomCharge[pId] * atomCharge[qIndex] / ( rprime ) );
    }
    
    p = p->strandAtom;
    
  } while( p );
  
  
#ifdef USE_MPI//reduce and sum data
  MPI_Allreduce ( &ee, &eeTmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ourComm );
  ee = eeTmp;
#endif
  
  return( ee );
  
}


// *ESTIMATE* pairwise  (intermolecular only) NB energy between particles.
// This value is used to then guide candidate MC moves, 
// and the more expensive full Hamiltonian is calculated
// for those candidates which pass initial selection.  
double hsc::chainEnb( Particle *p, double *eelec ){
  
  Cell     *curCell;
  Particle *curPart;
  double    evdw, ee;
  int       neighbourCount;

  neighbourCount = 0;
  ee        = 0.0;
  evdw      = 0.0;   
 
  do{
  
    //loop over particles in the same cell
    curCell = p->cell;
    curPart = curCell->firstParticle;
    do {    
      if( curPart->chainIndex != p->chainIndex ){
        neighbourCount++;
        pairEnb( p, curPart, &ee, &evdw );
      }
      curPart = curPart->next;
    } while (curPart);
    
    //loop over neighbour cells
    for(int j=0; j<26; j++){
      curCell = p->cell->neighbours[j];
      curPart = curCell->firstParticle;
      while (curPart){
        if( curPart->chainIndex != p->chainIndex ){
          neighbourCount++;
          pairEnb( p, curPart, &ee, &evdw );
        }
        curPart = curPart->next;
      }
    }
    
    p = p->strandAtom;
  } while( p );
  
 *eelec    = ee;
  
  return( ee + evdw );
  
}


// Pairwise EVDW according to AMBER Hamiltonian
double hsc::pairEnb( Particle *p, Particle *q, double *eelec, double *evdw ){
  
  double ev, ee;
  double rij[3], rijsq, r;
  double inv_r4, inv_r6, inv_r12, epsilon;
  int    pair, pairIndex;
   
  
  //load the epsilon used in softcoring the interaction between this atom pair
  pairIndex = (LJ_atomTypes[p->myId] - 1) * nLJ_atomTypes + LJ_atomTypes[q->myId] - 1;
  pair      = LJ_pairTypes[pairIndex];
  epsilon   = softcore_scale * scEpsilon( pair );
  
  //measure particle-particle distance
  for (int j=0; j<3; j++){
      rij[j] = p->R[j] - q->R[j];
      if( periodic )
	while (fabs(rij[j]) > box->halfx[j]) 
	  rij[j] -= copysign(box->x[j], rij[j]);
  } 
  rijsq   = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
  if( rijsq > ( rCutVdw * rCutVdw) ){
    return( 0.0 );
  }
  
  //fade rather than cut the electrostatics
  //fade =  tanh(rCutVdw - r);
  
  //softcoring:
  r = sqrt(rijsq);
  if( r < epsilon ){
     r = sigmoid(r, 0.8 * epsilon, epsilon );
  }
  rijsq = r * r;
  
  //coulomb energy
  //fade =  (1.0 - exp(-1+r/rcut))/(1.0 - exp(-1.0));
  ee   =  atomCharge[p->myId] * atomCharge[q->myId] / r ;
    
  //vdw energy
  rijsq   = 1.0 / rijsq;
  inv_r4  = rijsq * rijsq;
  inv_r6  = inv_r4 * rijsq;
  
  
  if( pair > 0 ){
    
    inv_r12 = inv_r6 * inv_r6;
  
    pair -= 1;
    ev  =  LJ_aCoeff[pair] * inv_r12; 
    ev -=  LJ_bCoeff[pair] * inv_r6; 
  } else{
    ev  = 0.0;
  }
  
//*logFile<< "atpair: " << p->myId << " " << q->myId << " evdw: " << ev  << " eelec: " << ee <<  " " << (*eelec) + ee << endl;
  
 *eelec = (*eelec) + ee;
 *evdw  = (*evdw) + ev;
  return( ev + ee );
  
}









