#include "config.h"
#ifdef USE_MPI 
#include "mpi.h"
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

#define MIX_ORDER_MOL  4
#define MIX_ORDER_ABS  2

//assume 2fs timestep,    100000 is therefore 200ps
#define SWITCH_SOFT_RATIO          0.5

#define INIT_TIMESTEPS_FWD       0 //10000
#define INIT_TIMESTEPS_REV       0 //10000

#define STEPL_EVERY                1

double hsc::mixFunc_mol(double lambda){

  //experimental mixing funtion: intended to ease the transition at high lambda, 
  //for use in conjunction with softcoring.
  double scale;
 
  if ( amberSoftcoring ) return( 1.0 );

  //scale lambda onto a smaller interval
  if( lambda <= SWITCH_SOFT_RATIO  ){
    lambda = 0.0;
  }else{
    lambda = (lambda - SWITCH_SOFT_RATIO ) / ( 1.0 - SWITCH_SOFT_RATIO );
  }
  
  scale = pow( lambda, MIX_ORDER_MOL );
  
  return( scale );
  

}
double hsc::mixFunc_soft(double lambda){

  //experimental mixing funtion: intended to ease the transition at high lambda, 
  //for use in conjunction with softcoring.
  double scale;
 
  if ( !emilSoftForce ) return( 0.0 );

  //`hump' function with zeroes at 0,1.
  scale  = lambda * lambda * pow(1.0 - lambda, 4.0);
  scale *= (729.0 / 16.0) * 0.01;
  

  return( scale );
}
double hsc::mixFunc_abs(double lambda){
  
  double scale;
  
  scale = pow( lambda, MIX_ORDER_ABS );
  
  return( scale );
  

}
double hsc::genForce_mol(double dvdl, double emodel, double lambda ){

  //experimental mixing funtion: intended to ease the transition at high lambda, 
  //for use in conjunction with softcoring.
  double genForce;
  
  // much easier just to let amber handle this
  if ( amberSoftcoring ) return( dvdl );

  //scale lambda onto a smaller interval
  if( lambda <= SWITCH_SOFT_RATIO ){
    lambda = 0.0;
  }else{
    lambda = (lambda - SWITCH_SOFT_RATIO) / (1.0 - SWITCH_SOFT_RATIO);
  }
  
  genForce = -MIX_ORDER_MOL * pow(lambda, MIX_ORDER_MOL - 1) * emodel;
  // + pow(lambda, MIX_ORDER_MOL ) * dvdl;
  
  return( genForce );
  

}

double hsc::genForce_abs(double emodel,  double lambda ){

  double genForce;
  
  //generalised force due to abstract model, which switches on with increasing lambda on [0..1]
  genForce = MIX_ORDER_ABS * pow( lambda, MIX_ORDER_ABS - 1) * emodel;
  
  return( genForce / beta );
}


double hsc::genForce_soft(double esoft, double lambda ){

  double genForce;
  
  if ( !emilSoftForce ) return( 0.0 );

  genForce  = esoft * ( 2.0 * lambda * pow(1.0 - lambda, 4.0) - 4.0 * lambda * lambda * pow(1.0 - lambda, 3.0) );
  genForce *= (729.0 / 16.0);
  
  return( genForce / beta );
}


//this code is for doing off-equilibrium TI: initial trials were not promising, so retired for the moment.
int hsc::stepLambda( double dHdL ){

  double            eGuess;
  time_t            stepEnd;
  double            stepLength;
  unsigned long int ETA;
  
  if( integrationFinished == true ){
    flushLog();
    return( EXIT_SUCCESS );
  }
  
  //(under)-estimate the standard error of the mean dHdL. 
  eGuess = eState_dHdL.accumulate( dHdL );
  
  *(logFile) << "L: " << restLambda << " dHdL running average: " << eState_dHdL.mean << " +/- " << eGuess;
  *(logFile) << " target error: " << 1.0/getLStep( eState_dHdL.mean ) << endl;
  
  
  if(  ( ( !stepLambdaBackwards && ( sim_nStep >= INIT_TIMESTEPS_FWD ) ) 
      || ( stepLambdaBackwards  && ( sim_nStep >= INIT_TIMESTEPS_REV ) ) ) ){                           
  
    
    //save the pair dHdL, lambda
    *logFile << " mean_at_step: "  << sim_nStep;
    *logFile << " lambda:       "  << restLambda;
    *logFile << " dHdL:         "  << eState_dHdL.mean;
    *logFile << " dHdL_n:       "  << eState_dHdL.n;
    *logFile << " dHdL_SD:      "  << eState_dHdL.SD;
    *logFile << " therm         "  << workInThermostat * beta << endl;

   
    
    //no need to make a printout every md step
    if( integrationStepCount % 10 == 0 ){
   
      //estimate time to convergence
      stepEnd       = time( NULL );
      stepLength    = (double)(stepEnd - stepStartTime);// / 10.0;
      stepStartTime = stepEnd;
   
      ETA           = (unsigned long int)stepLength * (1.0 - amberLambda)/ getLStep( eState_dHdL.mean );
    *(logFile) << "Duration was " << stepLength << " seconds" << endl;
    *(logFile) << "ETA for calculation is " 
             << ETA/3600 << " h " 
             << (ETA%3600)/60  << " ' " 
             << (ETA%3600)%60  << " ''" << endl;
    }
             
             
    integrationStepCount++;
    
    //test for finish
    if( ( ( amberLambda == 1.0 || restLambda == 1.0 ) && stepLambdaBackwards == false )
     || ( ( amberLambda <= 0.0 || restLambda <= 0.0 ) && stepLambdaBackwards == true  ) ){
      *(logFile) << "INTEGRATION STAGE COMPLETE" << endl;
      *(logFile) << "You must do the numerics yourself to get the work done in the cycle." << endl;
      
      
      integrationFinished = true;
      reportSystemFreeEnergy();
      flushLog();
#ifdef USE_MPI
      cerr << "Task: " << myTaskId << " exiting EMIL gracefully." << endl;
#endif      
      return( EXIT_SUCCESS );
    }
    
    //advance the integration.
    amberLambda += getLStep( eState_dHdL.mean );
    restLambda  += getLStep( eState_dHdL.mean );
  
    
    //test if we are starting the final datapoint
    if( amberLambda > 1.0 || restLambda > 1.0 ){
      amberLambda = 1.0;
      restLambda  = 1.0;
    }else if( amberLambda < 0.0 || restLambda < 0.0 ){
      amberLambda = 0.0;
      restLambda  = 0.0;
    }
    
      
    
    //save current value and refresh the estimator state
    dHdL_last  = eState_dHdL.mean;
    eState_dHdL.clean();
  
  } //end if making step
  return( 1 );
}

double hsc::getLStep( double dHdL_mean ){

  double dL;
  
  dL = 1.0 / (double)nStepsLimit;
  
  if( stepLambdaBackwards == true ){
    dL *= -1.0;
  }
  
  return( dL ); 
  
  
}
