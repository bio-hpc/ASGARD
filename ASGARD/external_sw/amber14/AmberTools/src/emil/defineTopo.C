
#include "config.h"
#ifdef USE_MPI
#include <mpi.h>
#endif

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <string>
#include <iostream>
#include <fstream>
#include <sstream> 

#include <stdio.h>

#include "emil.h"

using namespace std;
const double Pi = M_PI;


char *trim(char *str){ //string trimmer
  char *end;
  while(isspace(*str)) str++;
  if( *str==0 ){
    return(str);
  }
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;
  *(end+1) = 0;
  return(str);
}

//This is just a driver routine for testing a few things prior to integration with an MD program.
void hsc::run(){

#ifdef NAB_PROTO
  char *tName, *cName, *cIName;
 
  //nab has a global variable for output pointer, 
  //which has to be set before anything else is done. 
  nabout = stdout;
  
  // read new simulation parameters from ParameterFile
  readParamFile();

  if (newSetUp)setUpNew( (double *)0, (double *)0, (double *)0 );
  else setUpFromFile();

  writeOutConfig("StartConfig");
  writeOutWells("Wells");

  //equil();

  //convert filenames to nab string type
  tName  = (char *)trajOutFileName.c_str();
  cIName = (char *)crdFileName.c_str();
  cName  = (char *)crdOutFileName.c_str();
  
  //set mixing parameter
  lambdaMix = lambdaMix_start;
  emil_setLambdaMix( lambdaMix );
  scaleAmberParms( m, mRef, lambdaMix );
  //call nab setup - must re-call after any change to molecule parms
  //or lambdaMix
  if( !openWriteTraj( &tName ) ){
    MERR << "Couldn't open traj file " << tName << endl;
    exit( 8 );
  }
  initMM( &m, &cIName, &cName, workspace );
  while( lambdaMix < 1.0 ){
    
    for( int frame = 0; frame < lambdaEqFrames; frame++ ) {
      
        //compute forces for abstract model system, eg Einstein crystal.
        computeForces();
        
        //compute molecular forces and do one frame of MD 
        emilLoop( &m, particleCoords, forces, velocities );
        
        //get the molecular component of the generalised force for TI 
        generalisedForce  = emil_getDvDl();
        
        //book-keep changes in particle-well interaction
        updateParticles();
        
        //get the abstract component of the energy - also equal to 
        //the abstract component of the generalised force: dv/dl = v.
        abstractEnergy = computePhiOfConfig();
        
        printFrame( particleCoords, &N );
        
      *logFile<< "lambda mix: "      << lambdaMix << " ";
      *logFile<< "abstract force: "  << abstractEnergy << " ";
      *logFile<< "molecular force: " << generalisedForce << endl;
        
     } 
      
     //attempt exchange & relocation jump moves 
     MC();   
               
     //set mixing parameter
     lambdaMix += lambdaStep;
     emil_setLambdaMix( lambdaMix );
     scaleAmberParms( m, mRef, lambdaMix );
    
     //call nab setup - must re-call after any change to molecule parms
     //or lambdaMix
     initMM( &m, &cIName, &cName, workspace  );
    
  }  
    
#endif     
    
  writeOutConfig((char *)"Config");
  writeOutBindingToWells((char *)"BindingToWells");
  writeOutSimData((char *)"Stats");

}



/* 
 * Function to set up the potential types of the different attractive wells
 */
void hsc::setWellTypes(){
  
  int      *matchMask, i, ii, j, maskIndex, liquidIndex, chainIndex;
  int       matchCount, total, pIndex, nStrands;
  int       startNextRes, resCount;
  char     *maskString;
  int      *heavyMask;
  int      *strandIds, *strandRoots;
  Particle *p;
  
  liquidLists     = new Particle**[numMaterials];
  chainsInLiquid  = new int[numMaterials];
  liquidIndex     = 0;
  chainIndex      = 0;
  numLiquidChains = 0;
  
  strandRoots = new int  [N];
  strandIds   = new int  [N];
  heavyMask   = new int  [N];
  matchMask   = new int  [N];
  total       = 0;
  
  //identify all heavy atoms.
  //  markHeavyAtoms( &m, heavyMask );
  matchCount = 0;
  for( i = 0; i < parms->nAtoms; i++ ){

    if( parms->atMasses[i] > 1.9 ){//heavier than a quite heavy hydrogen
      heavyMask[i] = 1;
      matchCount++;
    }else{
      heavyMask[i] = 0;
    }

  }

 *logFile << "# marking " << matchCount << " of " << N << " atoms as heavy" << endl;
  if( matchCount == 0 ){
    MERR    << "# Warning! Marking " << matchCount << " of " << N << " atoms as heavy" << endl;
  }
  //identify each strand
 *logFile<< "# matching " << parms->nMolecules << " atoms as chain roots" << endl; 
  if( parms->nMolecules == 0 ){
    MERR << "# Warning! Matching " << parms->nMolecules << " atoms as chain roots" << endl; 
  }
  
  nStrands = parms->nMolecules;
  ii = 0;
  for( j = 0; j < nStrands; j++ ){
    strandRoots[j] = ii; 
    for( i = 0; i < parms->molSizes[j]; i++){
	    strandIds[ii++] = j;
    }
  }


  MERR << "matching masks for: " <<  numMaterials << " material types/phases." <<  endl; 
  
  //Set the type of each potential well according to atommask regular expressions.
  for( maskIndex = 0; maskIndex < numMaterials; maskIndex++){
    //!!!HACK ALERT: haven't implemented proper atom regular expression masks yet;
    //!!!temporary solution is to look for substring matches in a simple
    //!!!list of residue names.  


    //identify all matches to this well type
    maskString = (char *)materialAtomMask[maskIndex].c_str();

    matchCount = 0; 
    i          = 0;
    resCount   = 1;
    startNextRes = parms->resFirstAtom[1] - 1;
      
    for( j = 0; j < parms->nMolecules; j++ ){

      char resName[5], *r;
      resName[4] = '\0';  

      //check if resname is a substring of the mask
      memcpy(resName, &(parms->resNames[resCount*4]), 4 * sizeof(char));

      r=trim(resName);
      
      //check if resname is a substring of the mask
      if(strstr( maskString, r )){


	    for( ii = i; ii < i + parms->molSizes[j]; ii++ ){ 
	        matchMask[ii] = 1;
		
		if( EMIL_RESTRAIN_HATOMS != 0 || heavyMask[ii] != 0 ){
		  part[ii].wellType = materialType[maskIndex];
		  matchCount++;
		}else{
		  part[ii].wellType = HSC_POTENTIAL_NULL;
		}

            if( ii == startNextRes ){
                resCount++;
                if(resCount < parms->nRes ){
                    startNextRes =  parms->resFirstAtom[resCount] - 1;
                }else{
                    startNextRes = parms->nAtoms;
                }
            }
	    }
      }else{

	    for( ii = i; ii < i + parms->molSizes[j]; ii++ ){ 
	        matchMask[ii] = 0;
            if( ii == startNextRes ){
                resCount++;
                if(resCount < parms->nRes ){
                    startNextRes =  parms->resFirstAtom[resCount] - 1;
                }else{
                    startNextRes = parms->nAtoms;
                }
            }
	    }   
      }

      i += parms->molSizes[j];
    }

    resCount = 1;
    startNextRes = parms->resFirstAtom[1] - 1;
    for( ii = 0; ii < nStrands; ii++ ) {
        char resName[5];
        resName[4] = '\0';  
        memcpy(resName, &(parms->resNames[resCount*4]), 4 * sizeof(char));  
            if( ii == startNextRes ){
                resCount++;
                if(resCount < parms->nRes ){
                    startNextRes =  parms->resFirstAtom[resCount];
                }else{
                    startNextRes = parms->nAtoms;
                }
            }
    }


    //Feedback success or otherwise of matching each mask.
  *logFile<< "# mask " << materialAtomMask[maskIndex] << " sets type of " << matchCount << " atom wells as " 
                     << materialType[maskIndex] << " of " << HSC_POTENTIAL_NUMPOTS << endl;
   if( matchCount == 0 ){
      MERR << "EMIL ERROR! mask: " << materialAtomMask[maskIndex] << " not recognised, sets type of 0 atom wells" << endl;
      exit(8);
   }
   
    //did we populate a new liquid?
    if( matchCount > 0 
             && materialType[maskIndex] != HSC_POTENTIAL_EINSTEIN_CRYSTAL
             && materialType[maskIndex] != HSC_POTENTIAL_WOBBLIUM )
    {
      
      int countLightChains;
      countLightChains = 0;

      liquidLists[liquidIndex] = new Particle*[matchCount];
      pIndex     = 0;
      
      for( ii = 0; ii < nStrands; ii++ ) {
        
        //find root heavy atom i of strand ii
        i = strandRoots[ii];
        if( matchMask[i] == 1 ) {
          

          //root HEAVY atom, not just root atom
          if( ii + 1 == nStrands ){
            while( i < N && heavyMask[i] != 1 ) {
                i++;
            }
            if(i >= N) {
                countLightChains++;
                i = strandRoots[ii];
                heavyMask[i] = 1;
            }
          }
          else{
            while( i < strandRoots[ii+1] && heavyMask[i] != 1 ) {
                i++;
            }
            if(i >= strandRoots[ii+1]) {
                countLightChains++;
                i = strandRoots[ii];
                heavyMask[i] = 1;
            }
          }

              
          //save it as representative of the whole strand
          liquidLists[liquidIndex][pIndex++] = &part[i];
          part[i].liquidIndex                = liquidIndex;
            
          //build a list of particles which are in the same chain
          p             = &part[i];
          p->isRoot     =  true;
          p->chainIndex =  chainIndex;
          for( j = strandRoots[ii]; j < N && strandIds[j] == ii; j++ ){
               if( j != i ){
                   p->strandAtom  = &part[j];
                   p              = p->strandAtom;
                   p->chainIndex  = (pIndex-1);
                   p->liquidIndex = liquidIndex;
                   p->isRoot      = false;
                   //MERR << " leaf: " << j << endl;
               }
          }
          
          {
            int ccount;
            ccount = 0;
            p             = &part[i];
            while(p){
              p = p->strandAtom;
              ccount++;
            }
            if( ccount != 3 && ccount != 1 )
              MERR << "Warning: created a liquid chain of " << ccount << " atoms.  The code is only tested so far for water and salt.\n";
          }
          chainIndex++;
        }
      }

    *logFile<< "#    ...created liquid with exchange-symmetry between " <<  pIndex << " chains.\n";

     if(countLightChains > 0){
            MERR << "Warning: " << countLightChains 
                 << " chains had no heavy atoms. Free hydrogen, or mistake?\n"
                 << "Assigning first atom of each strand as root & making it an honorary heavy.\n";
     }

      chainsInLiquid[liquidIndex++] = pIndex;
      numLiquidChains              += pIndex;
    }
    
    total += matchCount;
  }
  
  numLiquids = liquidIndex;
  
  delete [] heavyMask;
  delete [] strandRoots;
  delete [] strandIds;
  delete [] matchMask;
  
}

#if 0
//align the restraints onto the (solute heavy-atom) coordinates, so that we restrain shape 
//but not position and orientation:
void hsc::alignRestraints(){

  MATRIX_T mat;
  double   delta2;
  
  //this nab-alike function should rotate the coords of all restraints based on best alignment of the solute restraints.
 delta2 = alignCoords( topSoluteAtom, N, soluteMask, restraintCoords, particleCoords, box->halfx, mat );
 *logFile<< "squared displacement of template to solute: " << delta2 << " was removed and template was rotationally aligned." << endl;
  
 
}
#endif
  

