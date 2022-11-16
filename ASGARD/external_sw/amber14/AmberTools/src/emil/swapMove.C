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

void hsc::swapMove(){
  
  int       pNr, lIndex, i, j;
  Particle *p, *q;
  Well     *w;
  bool      inside  = false;   // particle inside its own well? 
  bool      insideQ = false;   // other particle inside well?
  Cell     *curCell;
  Particle *curPart;
  double    rijsq;
  double    PhiOldP, PhiOldQ, PhiNewP, PhiNewQ;
  double    distpq, distqp;
  double    OldE, NewE, balanceFactor;
  int       chainsInW;

  //choose a random liquid chain
  pNr = int( numLiquidChains * mRan() );
  
  //identify the particular liquid it is in
  lIndex = 0;
  while( pNr >= chainsInLiquid[lIndex]){
    pNr -= chainsInLiquid[lIndex];
    lIndex++;
  }
    
  //find the root particle of the chain
  p  =  liquidLists[lIndex][pNr];
  w  =  p->theWell;

  // create lists of particles that are within rcut of w
  // refresh...
  w->clean();
  chainsInW = 0;

  //cerr << "Looking for particles in well w "  << endl;
  curCell = w->cell;
  curPart = curCell->firstParticle;
  while ( curPart ) { //loop over particles in the same cell as p's well
    
    if( curPart->isRoot && curPart->liquidIndex == lIndex ){
      rijsq = computeDistance(*curPart, w);
      if ( rijsq < rcut2Liquid  ) { //assign them to the well if they are inside rcut2 of it.
        curPart->insertToWell( w ); 
        chainsInW++;
        if (curPart == p) inside = true;
      }
    }
    curPart = curPart->nextWC;
  }

  //loop over particles in the 26 neighbour cells to p's well
  for( j = 0; j < 26; j++ ){
    curCell = w->cell->neighbours[j];
    curPart = curCell->firstParticle;
    while (curPart){
      if( curPart->isRoot && curPart->liquidIndex == lIndex  ){
        rijsq = computeDistance(*curPart, w);
        if (rijsq < rcut2Liquid ) {
          curPart->insertToWell( w );
          chainsInW++;
          if (curPart == p) inside = true;
        }
      }
      curPart = curPart->nextWC;
    }
  }
 
  // if well w of p is empty, move is rejected
  if (chainsInW == 0) {
    //cerr << "Rejected a swap move, no chain roots in region of target well\n";
    return;
  }
  else {
    if (inside){ //if p is already inside w.
      
      //select another random chain, with root particle q
      pNr = int( chainsInLiquid[lIndex] * mRan() );
      q   = liquidLists[lIndex][pNr];

      //easy to swap a chain with itself...
      if ( q == p ){ return; }

      //test if q is already in the same well as p, 
      curPart   = w->firstParticle;
      while( curPart ){
         if ( q == curPart ) {
            insideQ = true;
            break;
         }
         curPart=curPart->nextW;
      }
      if( insideQ ){
        // attempting to swap two particles inside well w
        balanceFactor = 1.0;
      }
      else {
        // attempt to swap p inside with q outside
        balanceFactor = float( chainsInLiquid[lIndex] ) / float( chainsInW );
      }
    }
    //p is outside its home well w
    else{
      
      // Pick a chain, with root particle q, from inside well w.
      pNr = int( chainsInW * mRan() );
      i   = 0;
      
      //get first candidate q
      q   = w->firstParticle;
      while( i < pNr ){
        q = q->nextW;
        i++;
      }
      
      balanceFactor =  float( chainsInW ) /  float( chainsInLiquid[lIndex] );    
    }

    //calculate the change in energy after the proposed swap
    PhiOldP = p->theWell->computePhiOfDist(p->rijsqW);
    PhiOldQ = q->theWell->computePhiOfDist(q->rijsqW);
    
    //swap wells so we can see what the new energy would be
    distpq  = computeDistance( *p, q->theWell );
    distqp  = computeDistance( *q, p->theWell );
    
    PhiNewP = p->theWell->computePhiOfDist(distqp);
    PhiNewQ = q->theWell->computePhiOfDist(distpq);
 
#if 0
    //include the effects of swapping the rest of the chain also
    pp = p->strandAtom;
    qq = q->strandAtom;
    while( pp ){
    
      PhiOldP += pp->theWell->computePhiOfDist(pp->rijsqW);
      PhiOldQ += qq->theWell->computePhiOfDist(qq->rijsqW);
      distqp   = computeDistance(*qq, pp->theWell );
      distpq   = computeDistance(*pp, qq->theWell );
      PhiNewP += pp->theWell->computePhiOfDist(distqp);
      PhiNewQ += qq->theWell->computePhiOfDist(distpq);
    
      pp = pp->strandAtom;
      qq = qq->strandAtom;
      
    }
#endif
         
    OldE = PhiOldP + PhiOldQ;
    NewE = PhiNewP + PhiNewQ;

    balanceFactor = log(balanceFactor);

    
    if ( log(mRan()) < balanceFactor + (OldE - NewE) * mixFunc_abs(restLambda) ){
    
#ifdef PARALLEL_MC
      //save the swap to check for collisions later
      allSwaps[myTaskId][numSwaps[myTaskId]++] = p->myId;
      allSwaps[myTaskId][numSwaps[myTaskId]++] = q->myId;
#else
      p->theWell->theParticle = q;
      q->theWell->theParticle = p;
      
      {
        Well *dummy = q->theWell;
        q->theWell = p->theWell;
        p->theWell = dummy;
      }
        
      q->rijsqW = computeDistance( *q );
      p->rijsqW = computeDistance( *p );
      swapAccept++;
#endif
    }   
      
      
    return;
  }
}

void hsc::reduceSwaps(){
  
  int         task, lastTask, i, j;
  Particle   *p, *q;

  for( task = 0; task < numTasks; task++ ){
    
#ifdef USE_MPI
    //share out the size of the list of swaps
    MPI_Bcast( &numSwaps[task], 1, MPI_INT, task, ourComm );
    
    //share out the list of accepted swaps from this task
    if( numSwaps[task] > 0 ) {
      MPI_Bcast( allSwaps[task], numSwaps[task], MPI_INT, task, ourComm );
    }else{
      continue;
    }  
#endif
    swapAccept += (numSwaps[task] / 2); 
    
    //resolve collisions.  
    //scaling is getting bad here, but remember we are only 
    //comparing *accepted* moves, which is not usually very many.
    for( i = 0; i < numSwaps[task]; i+= 2 ){
      for( lastTask = 0; lastTask <= task; lastTask++ ){
        for( j = 0; j < numSwaps[lastTask]; j+= 2 ){
          
          //test if we have a duplicate swap
          if(  allSwaps[lastTask][j]   == allSwaps[task][i] 
            || allSwaps[lastTask][j+1] == allSwaps[task][i]
            || allSwaps[lastTask][j]   == allSwaps[task][i + 1] 
            || allSwaps[lastTask][j+1] == allSwaps[task][i + 1] ){
           
            if(  ! ( lastTask == task && i == j ) ){
              allSwaps[task][i]   = -1; //cancel the duplicate swap
              swapAccept--;
              break;
            }
          }
          
        }
        if( allSwaps[task][i] == -1 ) break;//if we jumped out, keep jumping.
      }
      
      if( allSwaps[task][i] != -1 ){//if this swap is OK, then save it.
         Well *dummy; 
        
         p = &(part[allSwaps[task][i]]);
         q = &(part[allSwaps[task][i + 1]]);
         
         dummy = q->theWell;
         q->theWell = p->theWell;
         p->theWell = dummy;
         
         q->theWell->theParticle = q;
         p->theWell->theParticle = p;
         
         q->rijsqW = computeDistance( *q );
         p->rijsqW = computeDistance( *p );
      
      }
    }//end loop over proposed swaps  
  }//end loop over tasks
    
  //refresh the counter
  numSwaps[myTaskId] = 0;
  
}






