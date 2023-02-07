#ifndef HAVE_LINEAR_ASSIGNMENT_H
#define HAVE_LINEAR_ASSIGNMENT_H

#ifdef USE_MPI
#include <mpi.h>
#endif

#define _LOGGING
#define PRETRIM_OPTIMAL_PARTICLES

#include <float.h>

/* require emil.h for Well and Particle classes */
#include "emil.h"


/* Hack to get execution time of a code block*/
#include <time.h>
#ifndef SYSOUT_F
#define SYSOUT_F(msg,time)  fprintf(stderr,"%s %.9fs\n",msg,time) 
#endif
#ifndef speedtest__             
#define speedtest__(data)   for (long blockTime = NULL; (blockTime == NULL ? (blockTime = clock()) != NULL : false); SYSOUT_F(data, (double) (clock() - blockTime) / CLOCKS_PER_SEC))
#endif

/* Usage: */
//speedtest__("Block Speed: ")
//{
    // The code goes here
//}
/*...end speed test hack*/

/* Algorithm to solve particle-to well mapping exactly */
class LinearAssignment {


  /* problem size */
  int     total_size, N;
  double   lowerLimit; /* best cost */

  /* store costs of particle-well assignment pairs */
  float   *costMatrix;
  int     *intCostMatrix;

  /* current forward and inverse part-well mappings */
  int     *partToWell;
  int     *wellToPart;  

  bool      periodic;  

  Well     **well;
  Particle **part, **savepart;
  float      box[3];
  float      halfBox[3];
  double     invxCell[3];

  /* return costs */
  //inline float cm_lazy(int partId, int wellId){
  //  return( costMatrix[partId*N+wellId] );
  //}
  inline int   cm_lazy(int partId, int wellId){
    return( intCostMatrix[partId*N+wellId] );
  }

  /* get pair cost for given particle-well assignment */
  float getPairCost( int partId, int wellId );

  /* solve linear assignment problem using
     (serial) Jonker Volgenant algorithm */
  double jonkerVolgenantAssign();
  int    integerJonkerVolgenantAssign();
  float  originalLap();

public:
  /* assignment algorithm to use */
  HSC_ASSIGNMENT_T method;
 
  /* some logging information */
  ErrorState      *assignmentRate, *assignmentDeRate;

  /* partial re-initialise on change of particle positions */
 int refreshAssignment(double *dcost){

    int    countChanged;  
    int    ii;


    for (int i = 0; i < total_size; i++ ){ 
        part[i]       = savepart[i];
        well[i]       = part[i]->theWell;
        partToWell[i] = i;
        wellToPart[i] = i;
    }

#ifdef PRETRIM_OPTIMAL_PARTICLES
   /* don't include in the search particles which are already optimal
           ...because we know for this system that each particle has
           its own unique optimal well */
    ii = 0;
    for (int i = 0; i < total_size; i++ ){ 

        int pairCost;
        pairCost = getPairCost(i,i); 

        if( pairCost == 0 ){
             continue;
        }else{
             part[ii]       = part[i];
             well[ii]       = well[i];
             partToWell[ii] = partToWell[i];
             wellToPart[ii] = wellToPart[i];
             ii++;
        }
    }
    N = ii;     

    //MERR << "Cutting problem size to " << N << " from " << total_size << endl;    


    lowerLimit = 0.0;
    for (int i = 0; i < N; i++ ){ 
        for( int j = 0; j < N; j++ ){
            intCostMatrix[i*N+j] = int(getPairCost(i,j));
        }
        lowerLimit += intCostMatrix[i*N+i];
    }

#else
    lowerLimit = 0.0;
    for (int i = 0; i < N; i++ ){ 
        for( int j = 0; j < N; j++ ){
            costMatrix[i*N+j] = getPairCost(i,j);
        }
        lowerLimit += costMatrix[i*N+i];
    }
#endif


   *dcost  = lowerLimit;

//speedtest__("Linear assignment time: ")
{
    /* do the magic */
     switch(method){
      case HSC_ASSIGNMENT_JONKER_VOLGENANT:
            *dcost = (double)integerJonkerVolgenantAssign() - (*dcost);
      break;
      default:
      break;
     }
}

//   *dcost = (double)originalLap() - (*dcost);
     


    /* reset the particle-well assignments for the caller */
    countChanged = 0;
    for (int i = 0; i < N; i++ ){ 
      
        if( part[i]->theWell != well[partToWell[i]] ){

#if 0
           cerr << "part " << i << " @ " 
                << part[i]->R[0] << " " << part[i]->R[1] << " " << part[i]->R[2] << " " 
                << " from well " << part[i]->theWell << " @ " 
                << part[i]->theWell->R[0] << " " 
                << part[i]->theWell->R[1] << " " 
                << part[i]->theWell->R[2] << " "
                << " to "        << well[partToWell[i]] <<  " @ " 
                << well[partToWell[i]]->R[0] << " " 
                << well[partToWell[i]]->R[1] << " " 
                << well[partToWell[i]]->R[2] << " " << endl;

           cerr << " rij2 from " <<  part[i]->rijsqW << " to: " <<  cm_lazy(i, partToWell[i])<<endl;
#endif

           /*update pointers for the sake of caller routine*/
           part[i]->theWell  = well[partToWell[i]];
           part[i]->theWell->theParticle = part[i];

           /* also save a bit of changed state */
           part[i]->rijsqW      = cm_lazy(i, partToWell[i]);
            
           countChanged++;
        }
    }

    return( countChanged );

  }





  /* constructor */
  LinearAssignment( int        problem_size, 
                    Particle **parts_in, 
                    Box       *box_in, 
                    bool       is_periodic, 
                    HSC_ASSIGNMENT_T method_in ){


    total_size = problem_size;
    N          = problem_size;
    lowerLimit = 0.0;
    method     = method_in;
    //costMatrix = new  float[N*N];
    intCostMatrix = new int[N*N];
    partToWell = new  int[N];
    wellToPart = new  int[N];
    part       = new Particle *[N];
    savepart   = new Particle *[N];
    well       = new Well *[N];
    periodic   = is_periodic;

    assignmentRate       = new ErrorState;
    assignmentDeRate     = new ErrorState;

    /* save repeated dereferencing and casting 
       by keeping the box data locally */
    for( int i = 0; i < 3; i++ ){
       box[i]      = (float)box_in->x[i];
       halfBox[i]  = (float)box_in->halfx[i];
       invxCell[i] = 1.0 / box_in->xCell[i];
    }

    for (int i = 0; i < N; i++ ){ 

        savepart[i]   = parts_in[i];
        well[i]       = savepart[i]->theWell;
        partToWell[i] = i;
        wellToPart[i] = i;

    }
  }

};
#endif //HAVE_LINEAR_ASSIGNMENT_H


