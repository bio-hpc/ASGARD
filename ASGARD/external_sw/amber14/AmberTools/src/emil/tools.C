#include "config.h"
#ifdef USE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <sstream>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

using namespace std;
#include "emil.h"


static FILE *logWob  = NULL;
static FILE *vtfTraj = NULL;
void hsc::writeVtf(){
  int i,j,k;

  //open file and write header info
  if( !vtfTraj ){
      vtfTraj = fopen("wobTraj.vtf", "w");
      //for( i = 0; i < N; i++){
      //if( part[i].wellType == HSC_POTENTIAL_WOBBLIUM ){
      fprintf(vtfTraj, "atom %d:%d radius 1.0 name O type 0\n", 0,N-1);
          //  }
          //}
      for( i = 0; i < N; i++){
        if( part[i].wellType == HSC_POTENTIAL_WOBBLIUM ){
          for( j = 0; j <  part[i].theWell->num_bonds; j++ ){
            if( part[i].theWell->bondedWells[j] ){
              k = part[i].theWell->bondedWells[j]->theParticle->myId;
              if( i < k ){
                fprintf(vtfTraj, "bond %d:%d\n", i,k);
              }
            }
          }
        }
      }
      fprintf(vtfTraj, "\n");
  }
  fprintf(vtfTraj, "timestep ordered\n");
  for( i = 0; i < N; i++){     
    fprintf(vtfTraj, "%f %f %f\n", part[i].R[0],part[i].R[1],part[i].R[2]);
  }
  fprintf(vtfTraj, "\n");
}


void hsc::saveRelocAttempt( Particle *chainP, double *chainCoords, double log_pAccept ){

  int       i, bufIndex, ccIndex;

  bufIndex = 3 * nPart_hastings;
  ccIndex  = 0;
  
  movedP[nHastings] = chainP;
  
  //save the coordinates of each particle of the chain
  do{
    for(i = 0; i < 3; i++){
 //       if(fabs (chainCoords[ccIndex]) > 2 * box->x[i] ){
 //         cerr << "Error, attempted reloc a long way outside box." << endl;
 //         cerr << i << " " << chainCoords[ccIndex] << " " << box->halfx[i] << endl;
 //         exit(8);
 //       }
        oldConfBuf[bufIndex++] = chainCoords[ccIndex++];
    } 
    nPart_hastings++; //number of particles saved
    chainP = chainP->strandAtom;
  } while( chainP );
  
  //"proposal weight" for metropolis-hastings move
  logHastingsRatio[nHastings]  = log_pAccept;
  nHastings++;
}


Particle *hsc::pickLiquidChain( int *lIndex, int *nAtoms ) {

  Particle *p, *q;
  int       pNr, l, n;

  pNr = int((numLiquidChains)*mRan());
    
  //identify the particular liquid it is in
  l = 0;
  while( pNr >= chainsInLiquid[l]){
    pNr -= chainsInLiquid[l];
    l++;
  }
  
  //find the root particle of the chain
  p  =  liquidLists[l][pNr];
  n  = 1;
  q  = p->strandAtom;
  while( q ){
    n++;
    q = q->strandAtom;
  }
 *lIndex = l;
 *nAtoms = n;
  return( p );

}


  //temp func for debug
void hsc::debugWrap(int i, int newCell){
  
    if ( i >= N || i < 0 ){
      cerr << "error: moving particle " << i << endl;
      exit( 8 );
    }
    if ( newCell >= box->nCells || newCell < 0 ) {
      cerr << "error: moving particle to cell " << newCell << endl;
      exit( 8 );
    }
    
    part[i].moveBetweenCells( *(part[i].cell), cells[newCell]);

}

//track particles which have been assigned new coordinates
void hsc::updateParticles() {

  double invxCell[3];
  int    newCell;

  
  for(int i = 0; i < 3; i++){
     invxCell[i] = 1.0 / box->xCell[i];
  }
  
  
  for (int i = 0 ; i < N; i++){
    
    //track the r2 for particle-well attraction
    part[i].rijsqW = computeDistance( part[i] );
  
    //check if the particle has moved to a new cell
    newCell = selectCell(part[i].R, invxCell );
  
    if ( &cells[newCell] != part[i].cell ){
      part[i].moveBetweenCells( *(part[i].cell), cells[newCell]);
    }

  }

}

bool hsc::testChainOverlap(Particle *p){

  do{
    if( checkOverlapAll( *p ) ) {
      return( true );
    }
    p = p->strandAtom;
  }while( p );

  return( false );
  
}

void hsc::displaceParticle(double newP[3], double oldP[3]) {

  double displace;

  for (int i=0; i<3; i++){
 
    displace = (mRan() - 0.5)*maxDisplace;

    newP[i] = -int((oldP[i]+displace)*box->halfxInv[i])*box->halfx[i] 
      + fmod((oldP[i]+displace),box->halfx[i]);

  }

}

double hsc::chainEnergy(Particle *p, double *dist){

  double E;
  int    i;
  
  i = 0;
  E = 0.0;
  
  E = E + p->theWell->computePhiOfDist(dist[i++]);
  p = p->strandAtom;
  while( p ){
    E  = E + dist[i++]; //phi stored instead of dist for non-root
    p  = p->strandAtom;
  } 

  return( E );
}

// move the whole chain such that the root particle is somewhere inside its well
void hsc::relocateChain(Particle *p){

  bool      outside;
  double    A[3], translate[3];
  Particle *q;
  
  //find a point in the sphere r = sqrt(rcut2)
  outside = true;
  do {
        for (int i = 0; i < 3; i++){ 
          A[i] = (mRan() - 0.5)*dCutLiquid;
        }
        if ((A[0]*A[0] + A[1]*A[1] + A[2]*A[2]) < rcut2Liquid )
          outside = false;
  }while(outside);

  for (int i = 0; i < 3; i++){
      translate[i] = -p->R[i] - 
                   int((p->theWell->R[i]+A[i])*box->halfxInv[i])*box->halfx[i] 
                         + fmod((p->theWell->R[i]+A[i]),box->halfx[i]);
  }
  
  q = p;
  do{
    q->R[0] = q->R[0] + translate[0];
    q->R[1] = q->R[1] + translate[1];
    q->R[2] = q->R[2] + translate[2];
    q = q->strandAtom;
  }while( q );
  
}

// places particle inside potential well
void hsc::relocateParticle(Particle &p){

  bool outside = true;
  double A[3];

  while(outside){
    for (int i=0; i<3; i++){ 
      A[i] = (mRan() - 0.5)*dCutLiquid;
    }
    if ((A[0]*A[0] + A[1]*A[1] + A[2]*A[2]) < rcut2Liquid)
      outside = false;
  }

  for (int i=0; i<3; i++){
    p.R[i] = -int((p.theWell->R[i]+A[i])*box->halfxInv[i])*box->halfx[i] 
      + fmod((p.theWell->R[i]+A[i]),box->halfx[i]);
  }
  
}

void hsc::testLoop( Particle *p ){

  Particle *q;
  
  q = p->next;
  
  while( q ){
     if( q == p ){
     *logFile<< "looped the list with particle " << q << endl;
     *logFile<< "id " << q->myId;
       exit( 8 );
     }
    q = q->next;
  }

}

// checks for overlap of particle pNr with all particles from list
bool hsc::checkOverlapAll(Particle &p) {

  bool      overlap = false;
  Cell     *curCell;
  Particle *curPart;
  
  curCell = p.cell;
  curPart = curCell->firstParticle;

  
  // Local cell
  while (curPart){    
    overlap = checkOverlapTwo(p, *curPart);
    if (overlap) {
      if (curPart == &p) overlap = false;
      else break;
    }
    curPart = curPart->next;
  }

  
  // Neighbour cells
  if (!overlap){
    for(int j=0; j<26; j++){
      curCell = p.cell->neighbours[j];
      curPart = curCell->firstParticle;
   
      while (curPart){
    
        overlap = checkOverlapTwo(p, *curPart);
        if (overlap) {
          if (curPart == &p) {
            overlap = false;
          } else {
            break;
          }
        }
        curPart = curPart->next;
      }
    
      if( overlap ){
        break;
      }
    }
  }
      
  return overlap;
} 



bool hsc::checkOverlapTwo(Particle &p1, Particle &p2) {

    double rij[3], rijsq;

    for (int j=0; j<3; j++){
      rij[j] = p1.R[j] - p2.R[j];

      if( periodic )
         while (fabs(rij[j]) > box->halfx[j]) 
             rij[j] -= copysign(box->x[j], rij[j]);
    } 
    rijsq = rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
    
    //define "overlap" as approach closer than 0.25 angstroms.
    if (rijsq < 0.25) return true;
    else return false;
}


double hsc::computePhiOfConfig() {

  double p, phi;
#ifdef USE_MPI
  double phiTmp;
#endif
  
  phi                 = 0.0;
  abstractSolidEnergy = 0.0;
  
  for (int ii=0; ii < myNatoms; ii++){
    
    int i = myAtomIds[ii];
    
    //get potential energy
    if( part[i].wellType != HSC_POTENTIAL_WOBBLIUM ){
        p    = part[i].theWell->computePhiOfDist(  part[i].rijsqW );
    }else{
        p    = part[i].theWell->computePhi( box );
    }
    phi += p;
    
    if( part[i].wellType == HSC_POTENTIAL_EINSTEIN_CRYSTAL 
     || part[i].wellType == HSC_POTENTIAL_WOBBLIUM ){
        abstractSolidEnergy += p;
    }
  }

#ifdef USE_MPI //reduce and sum data
  MPI_Allreduce ( &phi, &phiTmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ourComm );
  phi  = phiTmp + abstractEnergyZero;
  MPI_Allreduce ( &abstractSolidEnergy, &phiTmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ourComm );
  abstractSolidEnergy = phiTmp;
#endif

  phi +=  abstractEnergyZero;
 
  return phi;

}

void hsc::adjustAttemptRates(){
  
  //adjust number of MC attempts with softening parameter for stability.
  nSwap  = int( nSwap  +  0.1 * nSwap  * ( targetSwapRate  * nSwap  - swapAccept  ) ) + 1;
  nReloc = int( nReloc +  0.1 * nReloc * ( targetRelocRate * nReloc - nHastings ) ) + 1;
  
}

//find appropriate negihbour-list cell
int hsc::selectCell(double R[3], double invxCell[3]) {

    double rWrap[3];
    int    newCell;
    
    //check if the particle has moved to a new cell
    for(int j = 0; j < 3; j++){ 
      rWrap[j] =  box->imageR( R[j], j);
    }
    
    newCell  =  int( ( rWrap[0] + box->halfx[0] ) * invxCell[0] );
    newCell +=  int( ( rWrap[1] + box->halfx[1] ) * invxCell[1] ) * box->nxCell[0];
    newCell +=  int( ( rWrap[2] + box->halfx[2] ) * invxCell[2] ) * box->nxCell[0] * box->nxCell[1];
    

    return newCell;
 
}     
//find appropriate neighbour-list cell
int hsc::selectCellW(double R[3], double invxCellW[3]) {

    double rWrap[3];
    int    newCell;
    
    //check if the particle has moved to a new cell
    for(int j = 0; j < 3; j++){ 
      rWrap[j] =  box->imageR( R[j], j);
    }
    
    newCell  =  int( ( rWrap[0] + box->halfx[0] ) * invxCellW[0] );
    newCell +=  int( ( rWrap[1] + box->halfx[1] ) * invxCellW[1] ) * box->nxCellW[0];
    newCell +=  int( ( rWrap[2] + box->halfx[2] ) * invxCellW[2] ) * box->nxCellW[0] * box->nxCellW[1];
    

    return newCell;
 
}     
//find (squared) distance, and vector, between particle and well accounting for PBCs
double hsc::computeDistance( Particle &p ) {

  double rijsq;
  Well  *w;
  
  w = p.theWell;
    
  if( p.isRoot ){
    p.rij[0] = p.R[0] - w->R[0];
    p.rij[1] = p.R[1] - w->R[1];
    p.rij[2] = p.R[2] - w->R[2];
  } else{
    p.rij[0] = p.R[0] - (w->rootP->R[0] + w->R[0]);
    p.rij[1] = p.R[1] - (w->rootP->R[1] + w->R[1]);
    p.rij[2] = p.R[2] - (w->rootP->R[2] + w->R[2]);
  }

  if ( periodic ){
	while( fabs(p.rij[0]) > box->halfx[0] )  p.rij[0] -= copysign(box->x[0], p.rij[0]);
	while( fabs(p.rij[1]) > box->halfx[1] )  p.rij[1] -= copysign(box->x[1], p.rij[1]);
	while( fabs(p.rij[2]) > box->halfx[2] )  p.rij[2] -= copysign(box->x[2], p.rij[2]);
  }
  rijsq = p.rij[0]*p.rij[0]+p.rij[1]*p.rij[1]+p.rij[2]*p.rij[2];

  return rijsq;
}

void Well::clean(){
    while (firstParticle){
      Particle *p;
      p             = firstParticle;
      firstParticle = firstParticle->nextW;
      p->prevW = 0;
      p->nextW = 0;
    }
    npart = 0;    
}
  
//find (squared) distance, and vector, between particle and well accounting for PBCs
double hsc::computeDistance( Particle &p, Well *w ) {

  double rijsq;

  if( p.isRoot ){
    p.rij[0] = p.R[0] - w->R[0];
    p.rij[1] = p.R[1] - w->R[1];
    p.rij[2] = p.R[2] - w->R[2];
  } else{
    p.rij[0] = p.R[0] - (w->rootP->R[0] + w->R[0]);
    p.rij[1] = p.R[1] - (w->rootP->R[1] + w->R[1]);
    p.rij[2] = p.R[2] - (w->rootP->R[2] + w->R[2]);
  }

  if ( periodic ){
	while( fabs(p.rij[0]) > box->halfx[0] )  p.rij[0] -= copysign(box->x[0], p.rij[0]);
	while( fabs(p.rij[1]) > box->halfx[1] )  p.rij[1] -= copysign(box->x[1], p.rij[1]);
	while( fabs(p.rij[2]) > box->halfx[2] )  p.rij[2] -= copysign(box->x[2], p.rij[2]);
  }
  rijsq = p.rij[0]*p.rij[0]+p.rij[1]*p.rij[1]+p.rij[2]*p.rij[2];

  return rijsq;
}

double hsc::addBonds_spcfW( Particle *oxy, double *forceOxy, double *forceH1, double *forceH2, double mix ){

  Particle *h1, *h2;
  double    deltaE;
  double    vecOH1[3], vecOH2[3], vecHH[3]; 
  double    dvdr[3], r1p[3], r2p[3];
  double    dr1, dr2, r1, r2, theta, dTheta;
  double    r1r2, invr1, invr2, invr1p, invr2p;
  
  //Parameter set for SPCFw of Wu, Tepper & Voth, J. Chem. Phys. 124:024503, 2006.
  
  r1    = 0.0;
  r2    = 0.0;
  r1r2  = 0.0;
  
  //Molecule has root-heavy atom listed first.
  h1 = oxy->strandAtom;
  h2 = h1->strandAtom;
  
  //measure displacement vectors
  for (int j = 0; j < 3; j++) {
      vecOH1[j] = h1->R[j] - oxy->R[j];
      vecOH2[j] = h2->R[j] - oxy->R[j];
      vecHH[j]  = h2->R[j] - h1->R[j];
      while( fabs(vecOH1[j]) > box->halfx[j] ) 
        vecOH1[j] -= copysign(box->x[j], vecOH1[j]);
      while( fabs(vecOH2[j]) > box->halfx[j] ) 
        vecOH2[j] -= copysign(box->x[j], vecOH2[j]);
      while( fabs(vecHH[j]) > box->halfx[j] ) 
        vecHH[j] -= copysign(box->x[j], vecHH[j]);
      
      r1r2  += vecOH1[j] * vecOH2[j];
      r1    += vecOH1[j] * vecOH1[j];
      r2    += vecOH2[j] * vecOH2[j];
  } 
  r1    = sqrt(r1);
  r2    = sqrt(r2);
  invr1 = 1.0 / r1;
  invr2 = 1.0 / r2;
  theta = acos( r1r2 * invr1 * invr2 );
  
  
  dr1    = r1    - R_HWATO_SPCF;
  dr2    = r2    - R_HWATO_SPCF;
  dTheta = theta - THETA_0_SPCF;
  
  //equation 1, Wu/Tepper/Voth 2006. 
  deltaE  =  0.5 * K_B_SPCF * (dr1*dr1 + dr2*dr2);
  deltaE +=  0.5 * K_A_SPCF * dTheta * dTheta;
  
  //dvdr for each of r1,r2,theta:
  dvdr[0] = K_B_SPCF * dr1;
  dvdr[1] = K_B_SPCF * dr2;
  dvdr[2] = K_A_SPCF * dTheta;
  
  //all forces linear in dvdr, so bring in the mixing here:
  dvdr[0] *= mix;
  dvdr[1] *= mix;
  dvdr[2] *= mix;
  
  
  //record bond-stretch forces on each of the three atoms
  forceOxy[0] += ( vecOH1[0] * invr1 * dvdr[0] + vecOH2[0] * invr2 * dvdr[1] );
  forceOxy[1] += ( vecOH1[1] * invr1 * dvdr[0] + vecOH2[1] * invr2 * dvdr[1] );
  forceOxy[2] += ( vecOH1[2] * invr1 * dvdr[0] + vecOH2[2] * invr2 * dvdr[1] );
   
  forceH1[0] -= vecOH1[0] * invr1 * dvdr[0];
  forceH1[1] -= vecOH1[1] * invr1 * dvdr[0];
  forceH1[2] -= vecOH1[2] * invr1 * dvdr[0];
  
  forceH2[0] -= vecOH2[0] * invr2 * dvdr[1];
  forceH2[1] -= vecOH2[1] * invr2 * dvdr[1];
  forceH2[2] -= vecOH2[2] * invr2 * dvdr[1];
  
  
  //find vectors perpendicular to each bond
  r1p[0] = vecOH2[0] - r1r2 * invr1 * vecOH1[0];
  r1p[1] = vecOH2[1] - r1r2 * invr1 * vecOH1[1];
  r1p[2] = vecOH2[2] - r1r2 * invr1 * vecOH1[2];
  r2p[0] = vecOH1[0] - r1r2 * invr2 * vecOH2[0];
  r2p[1] = vecOH1[1] - r1r2 * invr2 * vecOH2[1];
  r2p[2] = vecOH1[2] - r1r2 * invr2 * vecOH2[2];
  
  //normalise them and scale by the moment r1 * dvdr[2]
  invr1p  = r1p[0] * r1p[0];
  invr1p += r1p[1] * r1p[1];
  invr1p += r1p[2] * r1p[2];
  invr1p  = r1 * dvdr[2] / (sqrt(invr1p));
  invr2p  = r2p[0] * r2p[0];
  invr2p += r2p[1] * r2p[1];
  invr2p += r2p[2] * r2p[2];
  invr2p  = r2 * dvdr[2] / (sqrt(invr2p));
  for(int j = 0; j < 3; j++ ){
    r1p[j] *= invr1p;
    r2p[j] *= invr2p;
  }
  
  //add the resulting force vectors
  forceOxy[0] -= ( r1p[0]  + r2p[0] );
  forceOxy[1] -= ( r1p[1]  + r2p[1] );
  forceOxy[1] -= ( r1p[2]  + r2p[2] );
  
  forceH1[0] += r1p[0];
  forceH1[1] += r1p[1];
  forceH1[2] += r1p[2];
  
  forceH2[0] += r2p[0];
  forceH2[1] += r2p[1];
  forceH2[2] += r2p[2];
  
  
  return( deltaE );
}

double hsc::addBonds_spcfWater( Particle *oxy, double *forceOxy, double *forceH1, double *forceH2, double mix ){

  Particle *h1, *h2;
  double    deltaE;
  double    vecOH1[3], vecOH2[3], vecHH[3]; 
  double    dvdr[3];
  double    dr1, dr2, dr3, r1, r2, r3;
  
  //SHOULD use:
  //Parameter set for SPC/Fw of Wu, Tepper & Voth, J. Chem. Phys. 124:024503, 2006.
  
  //params for SPCF from Toukan and Rahman physRevB 1985... 
  //in "mdyn per angstrom" == millidynes? WTFFF?
  //1 dyne is 10 micronewtons, 1 dyn / angstrom is..
  //10^-8 N/A
  //1md/A is therefore 10^-11 N/A
  //1md/A is therefore 10^-11 J/A metre
  //1md/A is therefore 10^-1  J/A^2
  //1md/A is therefore 2.39005736 x 10^-5 kcal/A^2
  //1md/A is therefore 1.43932636 × 10^19 kcal/A^2 mol
  const double mdynToKcalPerAngstromMol = 239.005736;
  const double a =  9.331 * mdynToKcalPerAngstromMol;
  const double b =  2.283 * mdynToKcalPerAngstromMol;
  const double c = -1.469 * mdynToKcalPerAngstromMol;
  const double d =  0.776 * mdynToKcalPerAngstromMol;
  r1 = 0.0;
  r2 = 0.0;
  r3 = 0.0;
  
  //Molecule has root-heavy atom listed first.
  h1 = oxy->strandAtom;
  h2 = h1->strandAtom;
  
  //measure displacement vectors
  for (int j = 0; j < 3; j++) {
      vecOH1[j] = h1->R[j] - oxy->R[j];
      vecOH2[j] = h2->R[j] - oxy->R[j];
      vecHH[j]  = h2->R[j]  - h1->R[j];
      while( vecOH1[j] > box->halfx[j] ) 
        vecOH1[j] -= copysign(box->x[j], vecOH1[j]);
      while( vecOH2[j] > box->halfx[j] ) 
        vecOH2[j] -= copysign(box->x[j], vecOH2[j]);
      while( vecHH[j] > box->halfx[j] ) 
        vecHH[j] -= copysign(box->x[j], vecHH[j]);
      
      r1 += vecOH1[j] * vecOH1[j];
      r2 += vecOH2[j] * vecOH2[j];
      r3 += vecHH[j]  * vecHH[j];
  } 
  r1 = sqrt(r1);
  r2 = sqrt(r2);
  r3 = sqrt(r3);
  
  dr1 = r1 - R_HWATO_SPCF;
  dr2 = r2 - R_HWATO_SPCF;
  dr3 = r3 - R_HWATH_SPCF;
  
  //equation 2, Toukan and Rahman 1985. 
  deltaE  =  0.5 * a * (dr1*dr1 + dr2*dr2) + 0.5 * b * dr3*dr3;
  deltaE += ( c * (dr1 + dr2) * dr3 + d * dr1 * dr2 );
  
  //dvdr for each of r1,r2,r3:
  dvdr[0] = a * dr1 + c * dr3 + d * dr2;
  dvdr[1] = a * dr2 + c * dr3 + d * dr1;
  dvdr[2] = b * dr3 + c * (dr1 + dr2);
  
  //all forces linear in dvdr, so bring in the mixing here:
  dvdr[0] *= mix;
  dvdr[1] *= mix;
  dvdr[2] *= mix;
  
  //record forces on each of the three atoms
  forceOxy[0] += ( vecOH1[0] * dvdr[0] + vecOH2[0] * dvdr[1] );
  forceOxy[1] += ( vecOH1[1] * dvdr[0] + vecOH2[1] * dvdr[1] );
  forceOxy[2] += ( vecOH1[2] * dvdr[0] + vecOH2[2] * dvdr[1] );
  
  forceH1[0] += ( -1.0 * vecOH1[0] * dvdr[0] + vecHH[0] * dvdr[2] );
  forceH1[1] += ( -1.0 * vecOH1[1] * dvdr[0] + vecHH[1] * dvdr[2] );
  forceH1[2] += ( -1.0 * vecOH1[2] * dvdr[0] + vecHH[2] * dvdr[2] );
  
  forceH2[0] += ( -1.0 * vecOH2[0] * dvdr[1] - vecHH[0] * dvdr[2] );
  forceH2[1] += ( -1.0 * vecOH2[1] * dvdr[1] - vecHH[1] * dvdr[2] );
  forceH2[2] += ( -1.0 * vecOH2[2] * dvdr[1] - vecHH[2] * dvdr[2] );
  
  
  
  return( deltaE );
}


double hsc::addBond( double *forceVec, Particle *p, Particle *bondP, double rMin, double epsilon ){

  double bondLength, deltaE, deltaF, r[3], dFdr;
  
  for (int j = 0; j < 3; j++) {
      r[j] = p->R[j] - bondP->R[j];
      while( r[j] > box->halfx[j] ) 
        r[j] -= copysign(box->x[j], r[j]);
  } 
  bondLength = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
  
  deltaE = rMin - bondLength;
  deltaF = deltaE * epsilon;
  deltaE = deltaE * deltaE * epsilon * 0.5;
  
  dFdr = deltaF / bondLength; //normalise, so that r[] is a unit vector
  
  if( bondLength > rMin ){//if bond is stretched
    forceVec[0] += r[0] * dFdr;
    forceVec[1] += r[1] * dFdr;
    forceVec[2] += r[2] * dFdr;
  }
  return( deltaE );
}



void hsc::acceptHastings(){

   Particle *p;
  
  //*logFile<< "Accepting " << nHastings << " hastings moves of rootP ids:\n";
    
     //log the number of move acceptances
     relocAccept   += nHastings;
     
     do {
       nHastings--;
       p = movedP[nHastings];
  
       //recalculate forces on rootP for new position... or just do everything, outside the loop
      // dPhiDr = p->theWell->compute_dPhiOfDist( p->rijsqW );
        
       do {
         //*logFile << " " << p->myId << " " << p->R[0] << " " << p->R[1] << " " << p->R[2] << "\n";
          p = p->strandAtom;
       } while( p );
      //*logFile << endl;
     }while( nHastings > 0 );
  
     nHastings      = 0;
     nPart_hastings = 0;

}

int hsc::rejectHastings( double log_pAccept ){

    double    invxCell[3], invxCellW[3];
    Particle *q;
    int       i, ii, newCell, revCount;

    nHastings--;
    
  *logFile<< "Rejecting hastings moves" << endl;
    
    //because this move was tested and rejected, this alters the proposal distribution
    //for the next move that we test:
    //P(next_prop, this_rejected) = P(this_rejected) * Prior(next_prop)
    //therefore, because we will subtract this from the next acceptance rate,
    // for each attempt rejected, the chance of a final acceptance gets... larger
      
    //this whole thing reminds me of the Monty Hall paradox... consider that this section of the code is retired for now.
    if( nHastings >= 1 ){
    *logFile<< nHastings << " hastings moves still outstanding, so modifying their probabilities to posterior" << endl;
      logHastingsRatio[nHastings - 1] += log( 1.0 - exp(log_pAccept));
    }
    
    //need this to assign cell positions
    for(i = 0; i < 3; i++){
      invxCell[i] = 1.0 / box->xCell[i];
      invxCellW[i] = 1.0 / box->xCellW[i];
    }
    
    //decrement the counter
    q = movedP[nHastings];
    revCount = 0;
    do{
    revCount++; 
    q = q->strandAtom;
    } while( q );
    nPart_hastings = nPart_hastings - revCount;
    
    //reset coordinates and neighbour lists
    ii = 0;
    q = movedP[nHastings];
    do {
      for( i = 0; i < 3; i++ ){ //reset the coordinates
        q->R[i] = oldConfBuf[nPart_hastings*3 + ii]; 
        ii++;
      }
      
      //reset neighbourlist cell
      newCell = selectCell( q->R, invxCell );
      if (&(cells[newCell]) != q->cell){
        q->moveBetweenCells( *q->cell, cells[newCell]);
      }
      
      //reset well cell
      newCell = selectCellW( q->R, invxCellW );
      if (&(cellsW[newCell]) != q->WCell){
        q->moveBetweenCellsW( *q->WCell, cellsW[newCell]);
      }
      
      q = q->strandAtom;
    } while( q );
    
    return( revCount );

}

void hsc::replaceTrialReloc( Particle *q, double *invxCell, double *rOld, int *oldCellIndex ){

      int i, j, jj, newCell;
      
      j  = 0;
      jj = 0;
      
      //loop over the chain
      do {
        newCell = selectCell( q->R, invxCell );
        for ( i = 0; i < 3; i++){
          q->R[i] = rOld[jj++];
        }      
        if ( newCell != oldCellIndex[j] ){
          q->moveBetweenCells( cells[newCell], cells[oldCellIndex[j]] );
        }
        q = q->strandAtom;
        j++;
      }while( q );
  
    return;

}

//debug function, useful for setting water parameters
void hsc::printAverageWaterBonds(){

  Particle *q, *qq;
  double    r, rOH_mean, rHH_mean, sd[2];
  int       count;
  
  rOH_mean = 0.0;
  rHH_mean = 0.0;
  sd[0]    = 0.0;
  sd[1]    = 0.0;
  count    = 0;
  for( int i = 0; i < N; i++ ){
     if( part[i].isRoot && part[i].wellType != HSC_POTENTIAL_EINSTEIN_CRYSTAL ){
       
       if( part[i].strandAtom && part[i].strandAtom->strandAtom ){
        count++;
       
        q = part[i].strandAtom;
        r  = (q->R[0] - part[i].R[0])*(q->R[0] - part[i].R[0]) ;
        r += (q->R[1] - part[i].R[1])*(q->R[1] - part[i].R[1]) ;
        r += (q->R[2] - part[i].R[2])*(q->R[2] - part[i].R[2]) ;
        rOH_mean += sqrt( r );
       
        qq = q->strandAtom;
        r  = (qq->R[0] - part[i].R[0])*(qq->R[0] - part[i].R[0]) ;
        r += (qq->R[1] - part[i].R[1])*(qq->R[1] - part[i].R[1]) ;
        r += (qq->R[2] - part[i].R[2])*(qq->R[2] - part[i].R[2]) ;
        rOH_mean += sqrt( r );
       
        r  = (qq->R[0] - q->R[0])*(qq->R[0] - q->R[0]) ;
        r += (qq->R[1] - q->R[1])*(qq->R[1] - q->R[1]) ;
        r += (qq->R[2] - q->R[2])*(qq->R[2] - q->R[2]) ;
        rHH_mean += sqrt( r );
       }
     }
  } 
 
  rOH_mean /= double(2.0 * count);
  rHH_mean /= double(count);
  
  for( int i = 0; i < N; i++ ){
     if( part[i].isRoot && part[i].wellType != HSC_POTENTIAL_EINSTEIN_CRYSTAL ){
       
       if( part[i].strandAtom && part[i].strandAtom->strandAtom ){
       
        q = part[i].strandAtom;
        r  = (q->R[0] - part[i].R[0])*(q->R[0] - part[i].R[0]) ;
        r += (q->R[1] - part[i].R[1])*(q->R[1] - part[i].R[1]) ;
        r += (q->R[2] - part[i].R[2])*(q->R[2] - part[i].R[2]) ;
        sd[0] += ( sqrt( r ) - rOH_mean ) * ( sqrt( r ) - rOH_mean );
       
        qq = q->strandAtom;
        r  = (qq->R[0] - part[i].R[0])*(qq->R[0] - part[i].R[0]) ;
        r += (qq->R[1] - part[i].R[1])*(qq->R[1] - part[i].R[1]) ;
        r += (qq->R[2] - part[i].R[2])*(qq->R[2] - part[i].R[2]) ;
        sd[0] += ( sqrt( r ) - rOH_mean ) * ( sqrt( r ) - rOH_mean );
       
        r  = (qq->R[0] - q->R[0])*(qq->R[0] - q->R[0]) ;
        r += (qq->R[1] - q->R[1])*(qq->R[1] - q->R[1]) ;
        r += (qq->R[2] - q->R[2])*(qq->R[2] - q->R[2]) ;
        sd[1] += ( sqrt( r ) - rHH_mean ) * ( sqrt( r ) - rHH_mean );
       }
     }
  } 
  
*logFile<< "mean(sd) water OH,HH bond lengths: " 
        << rOH_mean << " ( "     << (sqrt(sd[0]/double(2.0*count))) << " ) " 
        << rHH_mean  << " ( "    << (sqrt(sd[1]/double(count))) << " ) " << endl;
  
}


double hsc::computeAverageDistanceFromWell(){
  double rijAbs = 0.0; 
  for (int i=0; i<N; i++){
    rijAbs += sqrt(part[i].rijsqW);
  }
  return (rijAbs*invNheavy_solid);
}

double hsc::computeRatioOutside(){

  int n=0;

  for (int i=0; i<N; i++){
    if (part[i].isRoot && part[i].rijsqW > rcut2Liquid) n++;
  }
  return (n*invNheavy_solid);
}

void hsc::reportWellPart(){
    double rms, x;
    int    c, j;
    c   = 0;
    rms = 0.0;
    
    if( !logWob ){
      logWob = fopen("wobbleWells.dat", "w");
    }

    for(int i = 0; i < N; i++){
      
      if( part[i].wellType == HSC_POTENTIAL_EINSTEIN_CRYSTAL ){
        for(j = 0; j < 3; j++ ){
          x =  part[i].R[j] - part[i].theWell->R[j];
          if ( periodic )
	         while( fabs(x) >= box->halfx[j] ) {
              x -= copysign(box->x[j], x);
             }        
          rms += x * x;
        }
        c++;
      } else if(  part[i].wellType == HSC_POTENTIAL_WOBBLIUM ){

        double f[3]; 
        part[i].theWell->compute_dPhi(f, box);
        for( int j = 0; j < 4; j++ ){
        
          double dx[3], *R, *r;
          R = part[i].theWell->bondedWells[j]->theParticle->R;
          r = part[i].R;
          for( int crd = 0; crd < 3; crd++){
            dx[crd] = R[crd] - r[crd];
	        if ( periodic )
	           while( fabs(dx[crd]) > box->halfx[crd] ) 
		           dx[crd] -= copysign(box->x[crd], dx[crd]);
          }
          fprintf(logWob, "%f %f %f %f %f %f %f %f %f\n",
                  r[0],r[1],r[2],
                  f[0],f[1],f[2],
                  dx[0],dx[1],dx[2]);
        }

      }
    }
    fprintf(logWob,"\n\n");
  *logFile<< "rms heavy/root-well: " << sqrt(rms/double(c)) << " over " << c << " heavies/roots" << endl;
    
}
void hsc::randomVector(double v[3]) {

    double rana, ranb, ransq, factor;

    ransq=2;
    while(ransq >= 1){
        rana = 1-2*mRan();
        ranb = 1-2*mRan();
        ransq=rana*rana+ranb*ranb;
    }
    factor = 2*sqrt(1-ransq);
    v[0] = rana*factor;
    v[1] = ranb*factor;
    v[2] = 1-2*ransq;
    
}
