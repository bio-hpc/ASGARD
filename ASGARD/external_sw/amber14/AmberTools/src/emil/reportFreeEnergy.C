#include "config.h"
#ifdef USE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <iomanip>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;
#include "emil.h"


////////////////////////////////////////////////////////
//
// These functions estimate free energy based on the 
// formulae in the supplementary material to 
// Berryman & Schilling, J Chem Theory Comput, 2013.
//
////////////////////////////////////////////////////////



//useful to be able to run some of these funcs standalone
#ifdef STANDALONE_FREE_ENERGY
int main(int argc, char* argv[]){
  
  int nArgs = 4;

  if( argc != nArgs ){
    fprintf(stderr, "Require %i args, got %i\n", nArgs, argc);
    fprintf(stderr, "   int:      N atoms\n");
    fprintf(stderr, "   float:    Temperature (Kelvin)\n");
    fprintf(stderr, "   filename: masses vector (1 per atom)\n");
    fprintf(stderr, "   filename: well type vector (1 per atom)\n");
    exit( 8 );
  }
}
#endif //STANDALONE_FREE_ENERGY


//return KE in AMBER units
double hsc::reportParticleKE(int i){

  double myKE;

  myKE  = velocities[3*i]     * velocities[3*i];
  myKE += velocities[3*i + 1] * velocities[3*i + 1];
  myKE += velocities[3*i + 2] * velocities[3*i + 2];
  
  myKE *= parms->atMasses[i];
  myKE *= 0.5;
  
  return( myKE );
  
}

//This function is currently a stub.
double hsc::setAbstractPotentialOffset( double eMol, double eAbs ){
  
  //abstractEnergyZero = eMol - ( eAbs + 1.5 * N );//add 3NkT/2 to the eabs term, it is not equil'd yet.
  
  abstractEnergyZero = 0.0;
  
  return( abstractEnergyZero );

}

double hsc::reportSystemEnthalpy(){
    *logFile << "#\n";
    *logFile << "# (mixed)  ABS_ENTH(kT) " << abstractEnergy * mixFunc_abs( restLambda );
    *logFile << " MOL_ENTH "            << emolec * mixFunc_mol( 1.0 - amberLambda ) * beta;
    *logFile << " TOT      "            << emolec * mixFunc_mol( 1.0 - amberLambda ) * beta + abstractEnergy * mixFunc_abs( restLambda ) 
 << endl;
   
    *logFile << "# (native) ABS_ENTH(kT) " << abstractEnergy;
    *logFile << " MOL_ENTH "            << emolec * beta;
    *logFile << " TOT      "            << emolec * beta + abstractEnergy 
<< endl;
    *logFile << "#\n";
    
    return(emolec * mixFunc_mol( 1.0 - amberLambda )  * beta + abstractEnergy * mixFunc_abs( restLambda ) 
    );
    
}
double hsc::reportSystemFreeEnergy(){

  std::ios  ioState(NULL);
  double    aSolid, aLiquid, liquidDegen;
  int       i, liquidIndex;
  char      out_cstring[256];  

  aSolid  = 0.0;
  aLiquid = 0.0;

 *logFile << "#\n# Reporting EMIL abstract Helmholtz free energy,\n# according to equations in supp. data of Berryman & Schilling, JCTC 2013: " << endl;

 
  for (i = 0; i < N; i++){
    if( part[i].wellType == HSC_POTENTIAL_EINSTEIN_CRYSTAL){
        aSolid += part[i].theWell->reportFreeEnergy();
    }else if( part[i].isRoot ){
        if( part[i].wellType == HSC_POTENTIAL_HARMONIC_LIQUID){
            aLiquid += part[i].theWell->reportFreeEnergy();
        }
        else{ 
            aLiquid += part[i].theWell->reportFreeEnergy( box->V );
        }
    }else{
        aLiquid += part[i].theWell->reportFreeEnergy();
    }
  }


  //find the extra term in the entropy due to exchange degeneracy of liquid particles.
  liquidDegen = 0.0;
  liquidIndex = 0;
  for( i = 0; i < numMaterials; i++ ){
     if( materialType[i] == HSC_POTENTIAL_HARMONIC_LIQUID
      || materialType[i] == HSC_POTENTIAL_CONEWELL_LIQUID
      || materialType[i] == HSC_POTENTIAL_WINGWELL_LIQUID ){
         
       liquidDegen += logFac( chainsInLiquid[liquidIndex] ); //log(N!)

       sprintf(out_cstring, "%6.6e", logFac( chainsInLiquid[liquidIndex] ));
       *logFile << "# liquid " << i << " has nchains: " <<  chainsInLiquid[liquidIndex] 
                              << " and log(N!): "   << out_cstring << endl; 
       liquidIndex++;
     }
  }
  sprintf(out_cstring, "%6.6e", (aSolid + aLiquid + liquidDegen));
 *logFile << "# Total free energy A of abstract model (units kT): " << out_cstring << endl;
  sprintf(out_cstring, "%6.6e", (aSolid));
 *logFile << "# A-Solid: " << out_cstring;
  sprintf(out_cstring, "%6.6e", (aLiquid));
 *logFile << " A-Liquid: " << out_cstring;
  sprintf(out_cstring, "%6.6e", (liquidDegen));
 *logFile  << " A-degen: " << out_cstring << "\n#" << endl;


  sprintf(out_cstring, "%6.6e", (aSolid + aLiquid + liquidDegen)/beta);
 *logFile << "# Total free energy A of abstract model (units of calling program): " << out_cstring << endl;
  sprintf(out_cstring, "%6.6e", (aSolid/beta));
 *logFile << "# A-Solid: " << out_cstring;
  sprintf(out_cstring, "%6.6e", (aLiquid/beta));
 *logFile << " A-Liquid: " << out_cstring;
  sprintf(out_cstring, "%6.6e", (liquidDegen/beta));
 *logFile  << " A-degen: " << out_cstring << "\n#" << endl;

*logFile  << "#\n# Generalised forces will be reported in units of calling program, not in kT/A\n#" << endl; 
  

  return(aSolid + aLiquid + liquidDegen);

}
