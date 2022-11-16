#ifdef USE_MPI
#include <mpi.h>
#endif
#include "config.h"

#include <iostream>
#include <sstream>

#include <math.h>
#include <stdlib.h>
#include <string.h>

using namespace std;
#include "emil.h"
#include "linearAssignment.h"

const double Pi = M_PI;

#ifndef PARALLEL_MC
void hsc::MC(){


  for( int i = 0; i < numLiquids; i++ ){
    if( linearAssignment[i] != NULL && sim_nStep % assignmentEvery == 0 ){

        int    didJv;
        double dcost;
 
        double oldSumR2, checkSumR2;


        oldSumR2 = 0.0;
        for( int j = 0; j < chainsInLiquid[i]; j++){
              oldSumR2 += liquidLists[i][j]->rijsqW;
        }

        //MERR <<  "enteringAssignment: " <<  didJv 
        //         << " cost:   "    <<  oldSumR2    <<  endl ;

        didJv = linearAssignment[i]->refreshAssignment(&dcost);

        checkSumR2 = 0.0;
        for( int j = 0; j < chainsInLiquid[i]; j++){
            checkSumR2 += liquidLists[i][j]->rijsqW;
        }

        //MERR <<  "reassigned: " <<  didJv  
        //         << " newcost:  "  <<  checkSumR2 <<  endl ;

        //log the number of wells reassigned and the change in energy/k. */
        linearAssignment[i]->assignmentRate->accumulate( didJv / (double)chainsInLiquid[i] );
        linearAssignment[i]->assignmentDeRate->accumulate( checkSumR2-oldSumR2 );

    }
  }


  //do swap moves
  while( swapMoves < nSwap ){ 
    swapMove();  
    swapMoves++;
  }

  //only call relocs when the amber hamiltonian is fully off.
  if( mixFunc_mol(1.0 - amberLambda) == 0.0 && nReloc > 0){ 
      while( relocMoves < nReloc ){ //try turning off the reloc moves gradually as we turn on the molecular potential
              //Attempt to teleport a chain.  
              relocAccept += nonHastingsReloc(); //"long jump" - useful to fill voids in the system
              relocMoves++;
      }
      hastingsFinalAccRate->accumulate( relocAccept / (double)nReloc );
  }

  if( nSwap > 0 ) {

    //log the acceptance rates
    swapAccRate->accumulate( swapAccept / (double)nSwap );
    swapAccept  = 0;
    swapMoves   = 0;
    relocAccept = 0;
    relocMoves  = 0;

    //*logFile <<  "SwapAccRate:   " <<  swapAccRate->mean;
    //if( mixFunc_mol(1.0 - amberLambda) != 0.0 ){
    //     *logFile << "\n";
    //}
   }
   //if( mixFunc_mol(1.0 - amberLambda) == 0.0 ){
   //      *logFile << " RelocAccRate:  " <<  hastingsFinalAccRate->mean << "\n"; 
   //}
}
#endif


//The Metropolis-Hastings parallel MC code is a dog's breakfast.
//Re-do from scratch at some point. 
#ifdef PARALLEL_MC
void hsc::MC(){

  int    ourSwaps, mySwaps;

#ifdef PARALLEL_MC
#ifdef USE_MPI //parallelism requires distinct RNG streams here
  cerr << "forking have: " << myTaskId << " " << mRan() << endl;
  mRan.fork(myTaskId);
#endif          
#endif          
  
  while( swapMoves < nSwap ){ 
    
      ourSwaps = nSwap - swapMoves;
#ifdef PARALLEL_MC
      if( ourSwaps > RESOLVE_SWAPS_EVERY){
        ourSwaps = RESOLVE_SWAPS_EVERY;
      }
      mySwaps  =  ourSwaps / numTasks;
      if( myTaskId < ourSwaps % numTasks ){
        mySwaps++;
      }
#else
      mySwaps  = ourSwaps;
#endif
      while( mySwaps > 0 ) {
    
               //Change the identity of a chain... no net effect except on 
               //the molecular configuration
              swapMove();  
              mySwaps--;

              if(mySwaps % 100 == 0){
MERR << "swapmoves: " << mySwaps << " of " << nSwap <<  " acc: " << swapAccept<< endl;
              }

              
      }
      //do a blocking reduce every few hundred attempts
      reduceSwaps();
      swapMoves += ourSwaps;
      
  }
#ifdef PARALLEL_MC
#ifdef USE_MPI
  cerr << "synching, have: " << mRan() << endl;
  mRan.synch(ourComm);
//#ifdef TEST_RNG_SYNCH
  {
    double values[2];
    values[myTaskId] = mRan();
    MPI_Bcast( &values[0], 1, MPI_DOUBLE, 0, ourComm );
    MPI_Bcast( &values[1], 1, MPI_DOUBLE, 1, ourComm );
    MERR << "synched, have: " << values[0] << " vs " << values[1] << endl;
  }  
//#endif
#endif 
#endif          
  
#ifdef WELL_MOVES 
  while( relocMoves < nReloc ){ 
      
      tryWellMove(); ///never call tryWellMove() and then follow it with tryRelocMove()... tryWellMove corrupts the well-particle lists.  
      relocMoves++;
  }
#else
   #ifdef HASTINGS
    while( relocMoves < nReloc ){ //try turning off the reloc moves gradually as we turn on the molecular potential
    //while( relocMoves < nReloc ){
              //Attempt to teleport a chain.  Cannot definitively accept/reject
              //with only emil information, so log the changes made and be prepared
              //to reverse them.
              tryRelocMove(); //"long jump" - useful to fill voids in the system
              relocMoves++;
  
              
              if( nHastings > 0 ) { //jump back out to test change in amber epot
                  hastingsPropRate->accumulate( 1.0 );
                  if( nHastings > 1 ){
                    cerr << "Error: currently set up to accumulate only one Hastings move at a time" << endl;
                    exit( 8 );
                  }
                  return;
              }else{
                  hastingsPropRate->accumulate( 0.0 );
              }
              
    }
   #else
   if( mixFunc_mol(1.0 - amberLambda) == 0.0 ){ //only call relocs when the amber hamiltonian is fully off.
      
      while( relocMoves < nReloc ){ //try turning off the reloc moves gradually as we turn on the molecular potential
              //Attempt to teleport a chain.  
              relocAccept += nonHastingsReloc(); //"long jump" - useful to fill voids in the system
              relocMoves++;
      }
   }
   #endif
#endif

  *logFile<< "Accepted " <<  swapAccept << "  of "  << swapMoves << " swap attempts and " 
          << relocAccept << " of " << nReloc << " relocs\n"; 

  *logFile <<  "SwapAccRate:   " <<  swapAccRate->mean 
           << " RelocPropRate: " <<  hastingsPropRate->mean 
           << " RelocAccRate:  " <<  hastingsFinalAccRate->mean << "\n"; `
  //reset for next frame


  MERR << "Reset? swapMoves: " << swapMoves << " of nSwap: " << nSwap << endl; 
  MERR << "relocMoves: " << relocMoves << " of scaled nReloc: " << int(nReloc * mixFunc_abs(restLambda)) << endl; 
  if( ( relocMoves >= int(nReloc * mixFunc_abs(restLambda)) || mixFunc_mol(1.0 - amberLambda) == 0.0 ) && swapMoves == nSwap ){
    
    
    //log the acceptance rates
    hastingsFinalAccRate->accumulate( relocAccept / (double)nReloc );
    swapAccRate->accumulate( swapAccept / (double)nSwap );
     
    relocAccept = 0;
    swapAccept  = 0;
    relocMoves  = 0;
    swapMoves   = 0;
  } 
}
#endif
