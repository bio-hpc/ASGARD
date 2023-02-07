#include "config.h"
#ifdef USE_MPI 
#include "mpi.h"
#endif

#include <iostream>
#include <iomanip>
#include <sstream>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;
#include "emil.h"

#define _SSTI_FORCE_WARN 130.0
#define _FORCE_CEIL      145.0


void hsc::computeForces_mix(){
  
  double    dPhiDr, normFac, amberKt;
  double    scaleMolec, scaleAbs, scaleSoft;

  amberKt    = 1.0 / beta; // inside SSTI, forces are in units kT/Angstrom
                           // in amber, forces are kCal/mol / Angstrom
                           //... the amber kT is already in these units.
  
  
  //if( myTaskId == 0 )reportWellPart();
  //if( myTaskId == 0 )writeVtf();

  //scale down the amber forces
  if( amberSoftcoring ){
    scaleMolec = 1.0;
  }else{
    scaleMolec = mixFunc_mol( 1.0 - amberLambda );
  }
  if( emilSoftForce ){
    scaleSoft  = mixFunc_soft( restLambda );
  }else{
    scaleSoft  = 0.0;
  }
  scaleAbs   = mixFunc_abs( restLambda );
 *logFile << scientific;

  
  //scale down the amber forces
  if( scaleMolec != 1.0 ){
     for (int ii = 0; ii < myNatoms; ii++){

         int i = myAtomIds[ii];

         forces[3*i]     *= scaleMolec;
         forces[3*i + 1] *= scaleMolec;
         forces[3*i + 2] *= scaleMolec;
      }
   }

   //add in the soft force
   softEnergy = 0.0;
   if( emilSoftForce ){
     for (int ii = 0; ii < myNatoms; ii++){
         int i = myAtomIds[ii];

       if( part[i].isRoot
        || part[i].wellType == HSC_POTENTIAL_EINSTEIN_CRYSTAL
        || part[i].wellType == HSC_POTENTIAL_WOBBLIUM ){
         softEnergy +=  softForceAll( &part[i], scaleSoft, &forces[3 * i] );
        }
     }
#ifdef USE_MPI //collect the total soft energy
     {
         double softEnergyTmp;
  MPI_Allreduce ( &softEnergy, &softEnergyTmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ourComm );
        softEnergy = softEnergyTmp;
     }
#endif
   }


  //do the EMIL Hamiltonian
  for (int ii = 0; ii < myNatoms; ii++){

      int i = myAtomIds[ii];
      
      //get force depending on well type (units of kT/Angstrom)
      dPhiDr = part[i].theWell->compute_dPhiOfDist( part[i].rijsqW );
      
      //add in the restraint-based force
      if( part[i].rijsqW > 0.0 ){//need to avoid an instability at r==0.
        
        //scale the abstract force, and
        //apply force as a vector to the well, accounting for PBCs
        normFac = amberKt / sqrt( part[i].rijsqW );
        dPhiDr *= normFac;
        
        //add the abstract forces
        forces[3 * i]     += dPhiDr * part[i].rij[0] * scaleAbs;
        forces[3 * i + 1] += dPhiDr * part[i].rij[1] * scaleAbs;
        forces[3 * i + 2] += dPhiDr * part[i].rij[2] * scaleAbs;
        
      }  
  }
}
 

void hsc::computeForces(){
  
  double dPhiDr, normFac, externalKt;
  int    atIndex;
  
  externalKt = 1.0 / beta;
  
  atIndex = 0;
  for (int i = 0; i < N; i++){
   
    //get force depending on well type (units of kT/Angstrom)
    if( part[i].isRoot ){
      dPhiDr = part[i].theWell->compute_dPhiOfDist( part[i].rijsqW );
    
      //apply force as a vector to the well, accounting for PBCs
      if( dPhiDr != 0.0 ){
        normFac = externalKt / sqrt( part[i].rijsqW );
        dPhiDr *= normFac;
        forces[atIndex++] = dPhiDr * part[i].rij[0];
        forces[atIndex++] = dPhiDr * part[i].rij[1];
        forces[atIndex++] = dPhiDr * part[i].rij[2]; 
      }
      else{
        forces[atIndex++] = 0.0;
        forces[atIndex++] = 0.0;
        forces[atIndex++] = 0.0; 
      }
    }
    else{//non-root particles have the force saved in rij already
        forces[atIndex++] = part[i].rij[0];
        forces[atIndex++] = part[i].rij[1];
        forces[atIndex++] = part[i].rij[2]; 
    }
  }
}
        
