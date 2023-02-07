#include "config.h"
#ifdef USE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

using namespace std;
#include "emil.h"
#include "linearAssignment.h"

#define _TEST_WELLS

void hsc::initialize(){

  double rij[3], rmsd;
  int    myCountHeavy, liquidIndex;
      
  dCutLiquid  = rcutLiquid  * 2.0;
  rcut2Liquid = rcutLiquid  * rcutLiquid;
  rcut3Liquid = rcut2Liquid * rcutLiquid;
  

  //seed the RNG
  mRan.seed(*iran);
  
  //set counters
  transAccept = 0; transMoves=0;
  swapAccept = 0; swapMoves = 0;
  relocAccept = 0; relocMoves = 0;

  //more counters
  hastingsPropRate     = new ErrorState; 
  hastingsFinalAccRate = new ErrorState;
  swapAccRate          = new ErrorState;
  
  //allocate some workspaces... note that this is used even in serial case; numtasks is 1.
  allSwaps = new int*[numTasks];
  numSwaps = new int[numTasks];
  memset( numSwaps, 0, numTasks * sizeof(int));
  for( int i = 0; i < numTasks; i++ ){
    allSwaps[i] = new int[RESOLVE_SWAPS_EVERY * 2];//have to allocate enough mem for all accepted swaps to be on the same task.
  }
  
  integrationFinished  = false;
  if( fixLambda == false ){
    if( stepLambdaBackwards == true ){
      restLambda  = 1.0;
      amberLambda = 1.0;
    }else{
      restLambda  = 0.0;
      amberLambda = 0.0;
    }
  }
  
  //work done against thermostat.... this is just for curiosity at the moment.
  workInThermostat     = 0.0;
  
  //work in particle-particle repulsion.
  softEnergy           = 0.0;
  
  meanPhi = 0;
  meanDistanceFromWell = 0;
  meanRatioOutside = 0;
  setupList();

  //cout << "# <phi(0)> = " << -4.0/9.0*M_PI*M_PI*rcut3*rcut3/box->V << endl;
  *logFile << "#\n# Setting well types\n";
  setWellTypes();
  
  
  *logFile << "#\n# Initialising wells\n";
  initWells();

  //count number of heavy or solid atoms
  invNheavy_solid = 0.0;
  for (int i=0; i<N; i++){
    if( part[i].isRoot ){
      invNheavy_solid = invNheavy_solid + 1.0;
    }
  }
*logFile<<"# " << invNheavy_solid << " heavy or solute atoms." << endl;
  invNheavy_solid = 1.0 / invNheavy_solid;
  
  
  //initial number of MC attempts - consider autoadapt of this for a steady 
  //number of moves accepted
  nSwap    = int(  numLiquidChains * swapTriesPerChain );
  nReloc   = int(  numLiquidChains * relocTriesPerChain );
  nMcMoves = nSwap + nReloc;
  
  for (int i=0; i<N; i++){
    for(int j=i+1; j<N; j++){
      for (int k=0; k<3; k++){
	    rij[k] = wells[i]->R[k] - wells[j]->R[k];
	    if( periodic )
	      while(fabs(rij[k]) > box->halfx[k]) 
	        rij[k] -= copysign(box->x[k], rij[k]);
      }
    }
  }
  
*logFile<<"# Building well neighbour list. " << endl;
  setupListW();

  
  //assign wells to particles
*logFile<<"# Linking particles to wells. " << endl;
  rmsd = 0.0;
  myCountHeavy = 0;
  for (int i=0; i<N; i++){
      double myR2;
      
      part[i].theWell       = wells[i];
      wells[i]->theParticle = &part[i];
      part[i].rijsqW        = computeDistance( part[i] );
      if(  part[i].isRoot || part[i].wellType == HSC_POTENTIAL_EINSTEIN_CRYSTAL ){
        rmsd += part[i].rijsqW;
        myCountHeavy++;
      }
      
      myR2  = pow(part[i].R[0] - wells[i]->R[0], 2.0);
      myR2 += pow(part[i].R[1] - wells[i]->R[1], 2.0);
      myR2 += pow(part[i].R[2] - wells[i]->R[2], 2.0);
      
  }
*logFile<<"# Linked particles to wells, rms distance (heavies): " << sqrt(rmsd/double(myCountHeavy)) << endl;
*logFile<<"# N heavies or roots: " << myCountHeavy << endl;
 
  /* set up a particle-well assignment solver for each liquid of
     indistinguishable particles */
  linearAssignment = new LinearAssignment *[numLiquids];
  liquidIndex = 0;
  for( int i = 0; i < numMaterials; i++ ){
    if( materialType[i] != HSC_POTENTIAL_EINSTEIN_CRYSTAL 
     && materialType[i] != HSC_POTENTIAL_WOBBLIUM
     && materialType[i] != HSC_POTENTIAL_NULL ){

 *logFile<< "# Molecule-Well assignment for indistinguishable fluid: " << i
         << " will be updated using:\n";

       if( materialType[i] == HSC_POTENTIAL_HARMONIC_LIQUID ) {
          linearAssignment[liquidIndex] 
                 =  new LinearAssignment( chainsInLiquid[liquidIndex], 
                                          liquidLists[liquidIndex], 
                                          box, 
                                          periodic, assignmentMethod );
 *logFile << "# " << HSC_ASSIGNMENT_NAME(assignmentMethod)  << "\n";
          
       }else{
          linearAssignment[liquidIndex] = NULL;
 *logFile << "# " << HSC_ASSIGNMENT_NAME(HSC_ASSIGNMENT_MC) << "\n";
       }
       liquidIndex++;
    }
  }  
 
}



//function to scatter the well centres so that the system is pre-equilibrated away from the 
//potential minimum
void hsc::randomiseWells() {

   double delta, rVec[3];
   int    i, ii;
   
   for( i = 0; i < N; i++) {
    
      //find the displacement equal to kT/2
      switch( part[i].wellType ){
        case HSC_POTENTIAL_EINSTEIN_CRYSTAL:
          delta = sqrt( wells[i]->rcut2 / (2.0 * wells[i]->epsilon) );
        break;
        case HSC_POTENTIAL_CONEWELL_LIQUID:
          delta = sqrt( wells[i]->rcut2) / (2.0 * wells[i]->epsilon);
        break;
        case HSC_POTENTIAL_HARMONIC_LIQUID:
          delta = sqrt(  wells[i]->rcut2 / (2.0 * wells[i]->epsilon) );
        break;
        case HSC_POTENTIAL_WINGWELL_LIQUID:
          delta = wells[i]->req / wells[i]->scale;
        break;
    
        default:
          cerr << "Error, no well type was set for particle " << i << endl;
          exit(8);
      }
      
      randomVector( rVec );
      for( ii = 0; ii < 3; ii++ ){
            wells[i]->R[ii] += ( delta * rVec[ii] );
      }
   }  
}
  
void hsc::initWells() {

  int    i;
  int   *typeCount;
   
  typeCount = new int[HSC_POTENTIAL_NUMPOTS];
  memset( typeCount, '\0', HSC_POTENTIAL_NUMPOTS * sizeof(int));

  wells = new Well*[N];
  #define _TEST_WELLS
  #ifdef TEST_WELLS
  //test the well distributions
  {
    ofstream outfile( "testWells.dat", ios::out );
    Well    *w;
    w = new harmonicWell;
    w->rCut      =  rcutLiquid;
    w->rcut2     =  rcutLiquid * rcutLiquid;
    w->rcut2Inv  =  1.0 / w->rcut2;
    w->epsilon   =  epsilonSolid;
    outfile << "#harmonic well\n";
    w->testMe( &outfile );
    delete w;
    
    outfile << "\n\n";
    
    w = new wingWell( req, rcutLiquid, epsilonTrap, wingForce );
    outfile << "#wing well\n";
    w->testMe( &outfile );
    delete w;
    
    outfile.close();
  }
  #endif
  
  //allocate the well type according to pattern matching done in setWellTypes()
  //...it would be nice to do this in blocks, to ensure contiguous memory,
  // but it probably falls out that way anyhow.
  for( i = 0; i < N; i++) {
    

#ifdef EMIL_REFUSE_WATER
      if( part[i].wellType != HSC_POTENTIAL_EINSTEIN_CRYSTAL){
           *logFile << "Error! EMIL does not currently support potential types other than a simple harmonic well." << endl;
            flushLog();

            MERR << "Error! EMIL does not currently support potential types other than a simple harmonic well." << endl;
            exit(8);
      }
#endif

      switch( part[i].wellType ){
        
        case HSC_POTENTIAL_EINSTEIN_CRYSTAL:


          wells[i] = new harmonicWellNoCut;
          wells[i]->rCut     =  rcutSolid;
          wells[i]->rcut2    =  rcutSolid * rcutSolid;
          wells[i]->rcut2Inv =  1.0 / wells[i]->rcut2;
          wells[i]->epsilon  =  epsilonSolid;
          wells[i]->req      =  5.0; //just set this nice and big so that no liquid wells are relocced into the solid.
          wells[i]->req2     =  wells[i]->req * wells[i]->req; 


        break;
      
        case HSC_POTENTIAL_CONEWELL_LIQUID:
          if( part[i].isRoot ){
          wells[i] = new coneWell; 
          wells[i]->rCut     =  rcutLiquid;
          wells[i]->rcut2    =  rcutLiquid * rcutLiquid;
          wells[i]->rcut2Inv =  1.0 / wells[i]->rcut2;
          wells[i]->epsilon  =  epsilonTrap;
          wells[i]->req      =  1.0 / epsilonTrap;
          wells[i]->req2     =  wells[i]->req * wells[i]->req; 
          }else{
            part[i].wellType = HSC_POTENTIAL_RELATIVE_POS;
            wells[i]         = new harmonicBondWell; //non-root particles are attracted to their root, not to an absolute well.
            wells[i]->rootP  = liquidLists[part[i].liquidIndex][part[i].chainIndex];
            wells[i]->rCut     =  rcutSolid;
            wells[i]->rcut2    =  rcutSolid * rcutSolid;
            wells[i]->rcut2Inv =  1.0 / wells[i]->rcut2;
            wells[i]->epsilon  =  epsilonSolid;
            wells[i]->req      =  5.0;
            wells[i]->req2     =  wells[i]->req * wells[i]->req;
          }

        break;
        
        case HSC_POTENTIAL_HARMONIC_LIQUID:

          if( part[i].isRoot ){

          wells[i] = new harmonicWellNoCut;
          wells[i]->rCut     =  rcutSolid;
          wells[i]->rcut2    =  rcutSolid * rcutSolid;
          wells[i]->rcut2Inv =  1.0 / wells[i]->rcut2;
          wells[i]->epsilon  =  epsilonSolid;
          wells[i]->req      =  5.0;
          wells[i]->req2     =  wells[i]->req * wells[i]->req;

          }else{
            part[i].wellType = HSC_POTENTIAL_RELATIVE_POS;
            wells[i]         = new harmonicBondWell; //non-root particles are attracted to their root, not to an absolute well.
            wells[i]->rootP  = liquidLists[part[i].liquidIndex][part[i].chainIndex];
            wells[i]->rCut     =  rcutSolid;
            wells[i]->rcut2    =  rcutSolid * rcutSolid;
            wells[i]->rcut2Inv =  1.0 / wells[i]->rcut2;
            wells[i]->epsilon  =  epsilonSolid;
            wells[i]->req      =  5.0;
            wells[i]->req2     =  wells[i]->req * wells[i]->req;
          }

        break;
        
        case HSC_POTENTIAL_WINGWELL_LIQUID:
          if( part[i].isRoot ){
            wells[i] = new wingWell( req, rcutLiquid, epsilonTrap, wingForce  ); 
          }else{
            part[i].wellType = HSC_POTENTIAL_RELATIVE_POS;
            wells[i]         = new harmonicBondWell; //non-root particles are attracted to their root, not to an absolute well.
            wells[i]->rootP  = liquidLists[part[i].liquidIndex][part[i].chainIndex];
            wells[i]->rCut     =  rcutSolid;
            wells[i]->rcut2    =  rcutSolid * rcutSolid;
            wells[i]->rcut2Inv =  1.0 / wells[i]->rcut2;
            wells[i]->epsilon  =  epsilonSolid;
            wells[i]->req      =  5.0;
            wells[i]->req2     =  wells[i]->req * wells[i]->req;
          }
        break;
    
        case HSC_POTENTIAL_WOBBLIUM:
          wells[i] = new wobbliumWell();
        break;

        case HSC_POTENTIAL_NULL:
	      wells[i] = new nullWell();
        break;
      
        default:
          cerr << "Error, no well type for particle " << i << endl;
          exit(8);
      }
      
      /* count the number of different particle types, just for interest */
      typeCount[part[i].wellType]++;

      /* set up this pointer: nice to have the wells R in a contiguous memory block */
      wells[i]->R        = &restraintCoords[3*i];

      /* log the parameters for this well type, if we have not mentioned it already */
      if( typeCount[part[i].wellType] == 1 ){
            *logFile << "# Well type "
                     << HSC_POTENTIAL_NAME(part[i].wellType) 
                     << " has parameters:\n";

            if( part[i].wellType == HSC_POTENTIAL_EINSTEIN_CRYSTAL
             || part[i].wellType == HSC_POTENTIAL_HARMONIC_LIQUID
             || part[i].wellType == HSC_POTENTIAL_CONEWELL_LIQUID
             || part[i].wellType == HSC_POTENTIAL_WOBBLIUM
             || part[i].wellType == HSC_POTENTIAL_RELATIVE_POS ) {

              *logFile << "#   r s.t. E=0: "<< wells[i]->rCut   << "\n";
              *logFile << "#   epsilon:    "<< wells[i]->epsilon<< "\n";

            }
            else if ( part[i].wellType == HSC_POTENTIAL_WINGWELL_LIQUID ) {
               ((wingWell*)wells[i])->report(logFile);
            }

            if( part[i].wellType == HSC_POTENTIAL_WINGWELL_LIQUID 
             || part[i].wellType == HSC_POTENTIAL_CONEWELL_LIQUID ){
   
             /* wells in which the potential is flat for some region
              * have a volume-dependent free energy */               
              *logFile << "#   A/well:    " << wells[i]->reportFreeEnergy(box->V) << "\n";

            }else{

              *logFile << "#   A/well:    " << wells[i]->reportFreeEnergy() << "\n";

            }
            *logFile << "#\n";
      }
      
  }

  if (readWells){
    
    ifstream infile;
    char line[2000];
    char quantity[2000];
    float number;
    int check;

    int allSet = 0;
    int count = 0;
    double Rin[3], msd;

    msd = 0.0;
    
   *logFile << "#\n# Reading well positions from " << wellFileName << endl;
      
    infile.open( wellFileName.c_str() );
    if (infile.bad()) {
      cerr << "Error " << wellFileName << " not found\n";
      exit(8);
    }
 
    while (infile.peek() != EOF) {
      infile.getline(line,sizeof(line));
      if(allSet == 5){

       check = sscanf(line, "%le %le %le", &Rin[0], &Rin[1], &Rin[2]);

       if (check == 3) { 
          if (count < N){
            for (int j=0; j<3; j++){
              wells[count]->R[j] = Rin[j];
              msd = msd + pow(wells[count]->R[j] - part[count].R[j], 2.0 ); 
            }
            count++;
          }
          else {
            cerr << "Error: Too many lines in " << wellFileName << endl;
            exit(8);
          }
        }
        else { 
          cerr << "Error: wrong line format in " << wellFileName << endl;
          cerr << "line was:\n" << line << endl;
          exit(8);
        }
      }
      else {
        check = sscanf(line, "%s%f", quantity, &number);
        if (check != 2) {
          cerr << "Error: wrong line format in " << wellFileName << endl;
          cerr << "line was:\n" << line << endl;
          exit(8);
        }
        if (strcmp(quantity,"#Number")==0) {
          if(N!=int(number)){
            cerr << "Wrong number of wells in " << wellFileName << endl;
            exit(8);
          }
          allSet++;
        }
        else if (strcmp(quantity,"#Density")==0) {allSet++;} 
        else if (strcmp(quantity,"#Boxx")==0) {
                //box->x[0]=number; 
                allSet++;
        }
        else if (strcmp(quantity,"#Boxy")==0) {
                //box->x[1]=number; 
                allSet++;
        }
        else if (strcmp(quantity,"#Boxz")==0) {
            //box->x[2]=number; 
            allSet++;
        }
      }
    }
  
    infile.close();

    if (count==0) {
      MERR << "\n\nEMIL Error reading wells position input file: \"" << wellFileName;
      MERR << "\"\n...note: if this is the first EMIL run (not a restart),";
      MERR << "  \n...then you do not neccesarily need a wells input file anyway!" << endl;
      exit(8);
    }
    
   *logFile<< "# msd part-well is: " << (msd/double(count)) << endl;
    
  }
  
  // use particle positions as positions of wells
  else {
    
    for (int i=0; i<N; i++){
      
      if( part[i].wellType == HSC_POTENTIAL_EINSTEIN_CRYSTAL || part[i].isRoot ){
        //set well position to particle position, but inside box
        for (int j=0; j < 3; j++){
          wells[i]->R[j] = part[i].R[j];
        }
      }else{
        //well has no real `position', instead it is a bond length restriction relative to the root particle.
        Particle *rootP;
        
        rootP = liquidLists[part[i].liquidIndex][part[i].chainIndex];
        
        for (int j=0; j < 3; j++){
          wells[i]->R[j] = part[i].R[j] - rootP->R[j];
        }
      }
    }
   *logFile<< "# didn't read well positions: positioned wells at particle coords\n";
   *logFile<< "# ...scaling velocities so that particles don't freeze because they are now at the bottom of a well.\n";
   
    //compensate for all the energy we have just dumped (into/out of) the system
    //by rescaling kinetic energies
    boostKEOnWellPlace();
    
  //  randomiseWells();
  // *logFile<< "# randomised well positions by a small delta.\n"; 
  }
  
  //image the wells... only image liquid wells.
  if( periodic )
  for (int i=0; i<N; i++){
    if( part[i].wellType != HSC_POTENTIAL_EINSTEIN_CRYSTAL 
     && part[i].wellType != HSC_POTENTIAL_WOBBLIUM 
     && part[i].isRoot ){
      for (int j=0; j < 3; j++){
          while( fabs(wells[i]->R[j]) > box->halfx[j] ){
            wells[i]->R[j] -= copysign( box->x[j], wells[i]->R[j] );
          }
      }
    }
  }
  
  //set up Elastic-network wells
  for( i = 0; i < N; i++) {

    if( part[i].wellType == HSC_POTENTIAL_WOBBLIUM ){
      double dR2;

      for(int j = i + 1; j < N; j++) {
        if( part[j].wellType != HSC_POTENTIAL_WOBBLIUM ) continue;

        dR2 = wells[i]->pairDist2(wells[j], box);
        if( dR2 < WOBBLIUM_BOND_CUT * WOBBLIUM_BOND_CUT ){
          if(  wells[i]->num_bonds >= WOBBLIUM_MAX_BONDS 
            || wells[j]->num_bonds >= WOBBLIUM_MAX_BONDS ){
            MERR << "EMIL Error, too many bonds at " << i << " , consider reducing distance cutoff" << endl;
            exit( 8 );
          }
          wells[i]->bond_eq_l2[wells[i]->num_bonds] = dR2;
          wells[i]->bondedWells[wells[i]->num_bonds] = wells[j];
          wells[i]->num_bonds++;
          wells[j]->bond_eq_l2[wells[j]->num_bonds] = dR2;
          wells[j]->bondedWells[wells[j]->num_bonds] = wells[i];
          wells[j]->num_bonds++;  
        }
      }                  
      MERR << "well " << i << " has " << wells[i]->num_bonds << " bonds " << endl;
      wells[i]->wobbliumAddBonds(wells[i]->bondedWells, wells[i]->num_bonds, 2.4, box);
    }
  }

  
 *logFile<< "# Finished setting particle well types:\n";
#ifdef EMIL_REFUSE_WATER
  for( i = 0; i <= HSC_POTENTIAL_EINSTEIN_CRYSTAL; i++ ){
     *logFile << "# " << HSC_POTENTIAL_NAME(i) << " " << typeCount[i] << "\n";
  }
#else
  for( i = 0; i < HSC_POTENTIAL_NUMPOTS; i++ ){
     *logFile << "# " << HSC_POTENTIAL_NAME(i) << " " << typeCount[i] << "\n";
  }
#endif



}

//compensate for all the energy we have just dumped (into/out of) the system
//by rescaling kinetic energies
void hsc::boostKEOnWellPlace(){
      
  Particle *q;
  double    myU, myKE, mulV;
  double    delta[3];
  int       myId;
      
   for (int i=0; i<N; i++){
     
      if( part[i].wellType == HSC_POTENTIAL_EINSTEIN_CRYSTAL || 
        ( part[i].wellType != HSC_POTENTIAL_WOBBLIUM && part[i].isRoot ) ){
        
        //get the PE which has just been taken out of the system
        myU  =  wells[i]->computePhiOfDist( 0.0 ) * mixFunc_abs( restLambda );
        if( myU < -1.5 ){
          myU = -1.5; //average PE of radial harmonic well should be some kT/2
          //... from trial and error it seems we need 6kT/2 which is pretty weird.
          //In case this is showing up a bug from somewhere else, I'm going to leave it as 3kT/2.
        }
        //..... maybe because of the leapfrom Verlet algorithm; our increased vels are getting
        //averaged down with the new vels?
        
        //If using a Langevin thermostat, then this will only "nearly" work anyway, 
        //because the added KE will bleed out through damping before the oscillation has finished a cycle. 
        
        
        //convert to AMBER units
        myU = myU / beta; 
       
        //get KE
        q    = &part[i];
        myKE = 0.0;
        do{
          myId  = q->myId;
          myKE += reportParticleKE( i );
          q     = q->strandAtom;
        }while( q );
        
        //compensate for the missing PE by increasing velocities
        mulV    = sqrt(1.0 - myU / myKE );
        
        delta[0]  =  velocities[i*3]   * (mulV - 1.0);
        delta[1]  =  velocities[i*3+1] * (mulV - 1.0);
        delta[2]  =  velocities[i*3+2] * (mulV - 1.0);
        
        q = &part[i];
        do{
          myId = q->myId;
          velocities[myId*3]  += delta[0];
          velocities[myId*3+1]+= delta[1];
          velocities[myId*3+2]+= delta[2];
          q = q->strandAtom;
        }while( q );
      }
    }
}




