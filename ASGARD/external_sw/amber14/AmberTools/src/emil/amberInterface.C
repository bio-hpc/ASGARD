//This file is an attempt to encapsulate the SSTI code to be called from inside pmemd

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
#include "linearAssignment.h"
const double Pi = M_PI;

//declare a static object of class Hsc to hold all of the (considerable) context for the emil calculation
static hsc Hsc;

static bool pmemd_called;

//Local function to parse incoming strings
void stringAssignIfNonNull(string *toSet, char *in_cstr, const char *default_cstr){

  //Use default filenames for i/o if none were passed in
  if( strlen(in_cstr) > 0 ) {
      toSet->assign(in_cstr);
  }
  else { 
      toSet->assign(default_cstr);
  }

}

#ifdef USE_MPI
//Local function to broadcast incoming strings
void stringBroadcast_Assign(string     *toSet, 
                            char       *in_cstr, 
                            const char *default_cstr, 
                            int         taskId, 
                            MPI_Comm    comm){

    int   tmp_strlen;
    char *tmp_cstr;

    /* find the length of the string at master task */
    if( taskId == 0 ){
        tmp_strlen = strlen(in_cstr) + 1;
    }
    MPI_Bcast( &tmp_strlen, 1, MPI_INT, 0, Hsc.ourComm);

    /* assign default if null */
    if( tmp_strlen <= 1 ){
        toSet->assign(default_cstr);
        return;
    }

    tmp_cstr = (char *)calloc(tmp_strlen, sizeof(char));
    if( taskId == 0 ){memcpy(tmp_cstr, in_cstr, tmp_strlen*sizeof(char));}
 
    /* share the c-string and allocate it as a C++ string */
    MPI_Bcast( tmp_cstr, tmp_strlen, MPI_CHAR, 0, comm);
    toSet->assign(tmp_cstr);

    free(tmp_cstr);
}
#endif


//declare C++ functions for external interface to emil
void emil_init(int    *nAtoms,
               int    *nRes,
               int    *nMolecules,
               int    *resFirstAtom,
               int    *molSizes,
               double *atMasses,
               char   *atNames,
               char   *resNames,
               int    *nAtom_types,
               double *rCutVdw,
               int    *atomType,
               int    *vdwPairType,
               double *charge,
               double *crd, 
               double *frc, 
               double *vel,
               double *box,
               double *LJ_aCoeff,
               double *LJ_bCoeff,
               double *amber_lambda, 
               double *restraint_lambda, 
               double *emil_softcore_scale, 
               double *beta,
               int      *myTaskId,
               int      *numTasks,
#ifdef USE_MPI
               MPI_Fint *commHandle,
#endif
               int      *nsteps_limit,
               int      *numex,
               int      *natex,
               char     *emil_paramfile,
               char     *emil_logfile,
               char     *emil_modelfile,
               char     *emil_model_outfile,
               int      *amber_softcoring
               );

void emil_forces( double *crd, 
                  double *frc, 
                  double *box, 
                  double *dvdl, 
                  double *molec_epot, 
                  double *ke, 
                  double *beta, 
                  int    *frameCount, 
                  int    *istart, 
                  int    *iend, 
                  int    *doPme );

void emil_mcAttempt(double *crd, 
                    double *frc, 
                    double *vel, 
                    double *box, 
                    double *epot, 
                    int    *emil_mc_done, 
                    double *beta, 
                    double *evdw_old, 
                    double *eelec_old);
                    
void emil_mcAcceptReject(double *crd, 
                          double *frc, 
                          double *vel, 
                          double *box, 
                          double *epot, 
                          double *beta, 
                          int    *emil_atoms_reverted, 
                          int    *emil_n_atoms_reverted, 
                          double *evdw, 
                          double *eelec);

//declare C functions which call the C++ interface functions... 
//...because fortran can call C but cannot (easily) call C++.
extern "C" { 
  void c_emil_init_(int    *nAtoms,  //atm_cnt
           int    *nRes,             //nres
           int    *nMolecules,       //gbl_mol_cnt
           int    *resFirstAtom,     //gbl_res_atoms
           int    *molSizes,         //atm_nsp
           double *atMasses,         //mass
           char   *atNames,          //atm_igraph
           char   *resNames,         //gbl_labres
           int    *nAtom_types,      //ntypes
           double *rCutVdw,          //vdw_cutoff
           int    *atomType,         //atm_iac
           int    *vdwPairType,      //typ_ico
           double *charge,           //atm_qterm
           double *crd,               
           double *frc, 
           double *vel, 
           double *box, 
           double *LJ_aCoeff,           //gbl_cn1 
           double *LJ_bCoeff,           //gbl_cn2
           double *amber_lambda,        //0
           double *restraint_lambda,    //0
           double *emil_softcore_scale, //emil_softcore_scale
           double *beta,                //beta
           int     *myTaskId,           //0
           int     *numTasks,           //1
#ifdef USE_MPI           
MPI_Fint *commHandle,
#endif
           int      *nsteps_limit,      //nstlim
           int      *numex,             //array of number of exclusions per atom
           int      *natex,             //array of atoms excluded,
           char     *emil_paramfile,    //filename to read params (can be null)
           char     *emil_logfile,      //filename for output log (can be null)     
           char     *emil_modelfile,    //filename for model coordinates (can be null),      
           char     *emil_model_outfile,//filename to save model coordinates (can be null)
           int      *amber_softcoring   //is mixing-out of AMBER forces handled outside?  
           ){ 


              emil_init( nAtoms,
                         nRes,
                         nMolecules,
                         resFirstAtom,
                         molSizes,
                         atMasses,
                         atNames,
                         resNames,
                         nAtom_types,
                         rCutVdw,
                         atomType,
                         vdwPairType,
                         charge,
                         crd, 
                         frc, 
                         vel, 
                         box, 
                         LJ_aCoeff,
                         LJ_bCoeff,
                         amber_lambda, 
                         restraint_lambda, 
                         emil_softcore_scale, 
                         beta,
                         myTaskId,
                         numTasks,
#ifdef USE_MPI                         
commHandle,
#endif 
                         nsteps_limit,
                         numex,
                         natex,
                         emil_paramfile,
                         emil_logfile,
                         emil_modelfile,
                         emil_model_outfile,
                         amber_softcoring
       );
  } 
  
  
  //Functions with a trailing underscore are called directly from external fortran program.
  void c_emil_forces_(double *crd, 
                      double *frc, 
                      double *box, 
                      double *dvdl, 
                      double *molec_epot, 
                      double *ke, 
                      double *beta, 
                      int    *frameCount, 
                      int    *istart, 
                      int    *iend, 
                      int    *doPme){ 
              pmemd_called = false;
              emil_forces(crd, frc, box, dvdl, molec_epot, ke, beta, frameCount, istart, iend, doPme); 
  }
  void c_emil_forces_pme_(double *crd, double *frc, double *box, double *dvdl, double *molec_epot, double *ke, double *beta, int *frameCount, int *my_atm_lst, int *my_at_cnt, int *doPme ){ 
      //just call the same version as from sander... leaving the pme_called flag as "true".       
      emil_forces(crd, frc, box, dvdl, molec_epot, ke, beta, frameCount, my_atm_lst, my_at_cnt, doPme); 
  }
  
  
  void c_emil_mcattempt_(double *crd, 
                         double *frc, 
                         double *vel, 
                         double *box, 
                         double *epot, 
                         int    *emil_mc_done, 
                         double *beta, 
                         double *evdw_old, 
                         double *eelec_old){
              emil_mcAttempt(crd, frc, vel, box, epot, emil_mc_done, beta, evdw_old, eelec_old);
  }
  void c_emil_mcacceptreject_(double *crd, 
                               double *frc, 
                               double *vel, 
                               double *box, 
                               double *epot, 
                               double *beta, 
                               int    *emil_atoms_reverted, 
                               int    *emil_n_atoms_reverted, 
                               double *evdw, 
                               double *eelec){
              emil_mcAcceptReject(crd, frc, vel, box, epot, beta, emil_atoms_reverted, emil_n_atoms_reverted, evdw, eelec );
  }
}





//function to setup thermodynamic integration to abstract material
void emil_init( int    *nAtoms,
                int    *nRes,
                int    *nMolecules,
                int    *resFirstAtom,
                int    *molSizes,
                double *atMasses,
                char   *atNames,
                char   *resNames,
                int    *nAtom_types,
                double *rCutVdw,
                int    *atomType,
                int    *vdwPairType,
                double *charge,
                double *crd, 
                double *frc, 
                double *vel, 
                double *box, 
                double *LJ_aCoeff,
                double *LJ_bCoeff,
                double *amber_lambda, 
                double *restraint_lambda, 
                double *emil_softcore_scale, 
                double *beta,
                int      *myTaskId,
                int      *numTasks,
#ifdef USE_MPI                
MPI_Fint *commHandle,
#endif
                int      *nsteps_limit,     
                int      *numex,             //array of number of exclusions per atom
                int      *natex,             //array of atoms excluded,
                char     *emil_paramfile,
                char     *emil_logfile,
                char     *emil_modelfile,
                char     *emil_model_outfile,
                int      *amber_softcoring
                ){


  Parms *parms;

  ofstream dbfile_atomTypes, dbfile_pairTypes, dbfile_atPairVDW;

  pmemd_called = true;//default.



#ifdef USE_MPI 
  //Link up the MPI comm and do a sanity check.
  Hsc.ourComm         =   MPI_Comm_f2c( *commHandle );
  
  //test the MPI comm
  cerr << "MPI testing comm.\n" << endl;
  MPI_Comm_size( Hsc.ourComm, &Hsc.numTasks);
  cerr << "task: " << *myTaskId << " of: " << Hsc.numTasks << endl;
  MPI_Barrier(Hsc.ourComm);
  Hsc.myTaskId        =  *myTaskId;
  if( Hsc.numTasks  !=  *numTasks ){
    cerr << "MPI fail here.\n" << endl;
    exit( 8 );
  }
#else
  Hsc.myTaskId        =  0;
  Hsc.numTasks        =  1;
#endif


  Hsc.parms           =  new Parms;
  parms               =  Hsc.parms;
#ifdef USE_MPI
  if(*myTaskId == 0 ){
#endif
  parms->nAtoms       = *nAtoms;
  parms->nAtomTypes   = *nAtom_types;
  parms->nRes         = *nRes;
  parms->nMolecules   = *nMolecules;
  parms->resFirstAtom =  resFirstAtom;
  parms->atNames      =  atNames;
  parms->resNames     =  resNames;
  parms->molSizes     =  molSizes;
  parms->atMasses     =  atMasses;
  parms->atCharges    =  charge;
  parms->atomTypes    =  atomType;
  parms->pairTypes    =  vdwPairType;
#ifdef USE_MPI
  }
#endif


#ifdef USE_MPI
  if(*myTaskId == 0 ){
#endif
    //do a quick sanity check of the inputs
    if( parms->nAtoms <= 0 ){
      cerr << *myTaskId << ": Error setting up emil: nAtoms = " << parms->nAtoms << endl;
      exit(8);
    }
    if( parms->nMolecules <= 0 ){
      cerr << *myTaskId << ": Error setting up emil: nMolecules = " << parms->nMolecules << endl;
      exit(8);
    }
    if( *myTaskId == 0 && ( !atNames || !resNames ) ){
      cerr << *myTaskId << ": Error setting up emil: null input data from AMBER"<< endl;
      cerr << *myTaskId << ": atNames:  " <<endl; cerr << atNames << endl;
      cerr << *myTaskId << ": resNames: " <<endl; cerr << resNames << endl;
      exit( 8 );
    }
    if( !resFirstAtom || !atMasses || !molSizes ){
      cerr << *myTaskId << ": Error setting up emil: null input data from AMBER"<< endl;
      cerr << *myTaskId << ": resFirstAtom " << resFirstAtom << endl;
      cerr << *myTaskId << ": atMasses " << atMasses << endl;
      cerr << *myTaskId << ": molSizes " << molSizes << endl;
      exit(8);
    }
#ifdef USE_MPI
  }
#endif

#ifndef USE_MPI
  if(*myTaskId != 0 ){
     cerr << "Warning!! taskId is :" <<  *myTaskId << " but USE_MPI was not defined at compile time!!!"<< endl;
  }
#endif

#ifdef USE_MPI
  //some codes (sander, I'm looking at you) don't store
  //resnames etc. for the non-master tasks.
  MPI_Bcast(&parms->nAtoms,     1, MPI_INT, 0, Hsc.ourComm );
  MPI_Bcast(&parms->nAtomTypes, 1, MPI_INT, 0, Hsc.ourComm );
  MPI_Bcast(&parms->nRes,       1, MPI_INT, 0, Hsc.ourComm );
  MPI_Bcast(&parms->nMolecules, 1, MPI_INT, 0, Hsc.ourComm );
   

  if(*myTaskId != 0 ){
       parms->resFirstAtom = new int[parms->nRes];
       parms->atomTypes    = new int[parms->nAtoms];
       parms->pairTypes    = new int[parms->nAtomTypes*parms->nAtomTypes];
       parms->atNames      = new char[parms->nAtoms * 4];
       parms->resNames     = new char[parms->nRes * 4];
       parms->atMasses     = new double[parms->nAtoms];
       parms->atCharges    = new double[parms->nAtoms];
       parms->molSizes     = new int[parms->nMolecules];
  }
     
  MPI_Bcast(parms->resFirstAtom,  parms->nRes,  MPI_INT, 0, Hsc.ourComm );
  MPI_Bcast(parms->atomTypes, parms->nAtoms,    MPI_INT, 0, Hsc.ourComm );
  MPI_Bcast(parms->pairTypes, parms->nAtomTypes*parms->nAtomTypes, MPI_INT, 0, Hsc.ourComm );
  MPI_Bcast(parms->atNames,  parms->nAtoms * 4, MPI_CHAR, 0, Hsc.ourComm );
  MPI_Bcast(parms->resNames,   parms->nRes * 4, MPI_CHAR, 0, Hsc.ourComm );
  MPI_Bcast(parms->atMasses, parms->nAtoms,     MPI_DOUBLE, 0, Hsc.ourComm );
  MPI_Bcast(parms->atCharges, parms->nAtoms,    MPI_DOUBLE, 0, Hsc.ourComm );
  MPI_Bcast(parms->molSizes, parms->nMolecules, MPI_INT, 0, Hsc.ourComm );
#endif


#ifdef USE_MPI
  //Share out filenames from master.  Not actually used in all cases, but simplest to send all.
  stringBroadcast_Assign(&Hsc.parameterFile, emil_paramfile,
                                      "emilParameters.in", *myTaskId, Hsc.ourComm);
  stringBroadcast_Assign(&Hsc.logFileName,   emil_logfile, 
                                      "emil.log", *myTaskId, Hsc.ourComm);
  stringBroadcast_Assign(&Hsc.wellFileName,  emil_modelfile, 
                                      "", *myTaskId, Hsc.ourComm);
  stringBroadcast_Assign(&Hsc.wellSaveFileName, emil_model_outfile, 
                                      "", *myTaskId, Hsc.ourComm);
#else
  //Use default filenames for i/o if none were passed in
  stringAssignIfNonNull(&Hsc.parameterFile, emil_paramfile, "emilParameters.in");
  stringAssignIfNonNull(&Hsc.logFileName,   emil_logfile,   "emil.log");
  stringAssignIfNonNull(&Hsc.wellFileName,  emil_modelfile, "");
  stringAssignIfNonNull(&Hsc.wellSaveFileName,  emil_model_outfile,  "");
#endif

  Hsc.readWells  = true;
  Hsc.writeWells = true;
  if(Hsc.wellFileName.length() == 0 ){
    MERR << "No input well filename supplied: placing wells at atom positions" << endl;
    Hsc.readWells = false;
  }
  if(Hsc.wellSaveFileName.length() == 0 ){
    MERR << "No output well filename supplied, so not writing well positions." << endl;
    Hsc.writeWells = false;
  }

  Hsc.N =  parms->nAtoms;
  Hsc.newSetUp        = true;
  Hsc.nLJ_atomTypes   =  parms->nAtomTypes;
  Hsc.LJ_atomTypes    =  parms->atomTypes;
  Hsc.atomCharge      =  parms->atCharges;
  Hsc.LJ_pairTypes    =  parms->pairTypes;
  Hsc.LJ_aCoeff       =  LJ_aCoeff;
  Hsc.LJ_bCoeff       =  LJ_bCoeff;
  Hsc.rCutVdw         = *rCutVdw;
  Hsc.stepStartTime   = time( NULL );
  Hsc.eState_dHdL.clean();
  Hsc.dHdL_last       = 0.0;
  Hsc.integrationStepCount = 0;
  Hsc.nStepsLimit          = *nsteps_limit;

  /* decide what is do be done locally and what in AMBER */
  Hsc.amberSoftcoring =  bool(*amber_softcoring);  
  if( Hsc.amberSoftcoring ) Hsc.emilSoftForce = false;
  else                      Hsc.emilSoftForce = true;

  Hsc.gfAccum_mol   = 0.0;
  Hsc.gfAccum_soft  = 0.0;
  Hsc.gfAccum_abs   = 0.0;
  Hsc.gfAccum_count = 0;


  //initialise the exclusion lists to a more friendly format than 
  //that provided by amber: natex needs to be randomly addressable
  int ex_accum;
  ex_accum = 0;
  Hsc.numex          = new int[Hsc.N];

#ifdef USE_MPI
  if( Hsc.myTaskId == 0 )
#endif
  for(int i=0; i < Hsc.N; i++){
        Hsc.numex[i] = numex[i] + ex_accum;
        ex_accum    += numex[i];
  }
#ifdef USE_MPI
  MPI_Bcast(Hsc.numex, Hsc.N, MPI_INT, 0, Hsc.ourComm );   
  MPI_Bcast(&ex_accum, 1, MPI_INT, 0, Hsc.ourComm );   
#endif

  Hsc.natex          = new int[ex_accum];
#ifdef USE_MPI
  if( Hsc.myTaskId == 0 )
#endif
  for(int i=0; i < ex_accum; i++){
    Hsc.natex[i] = natex[i] - 1; //Fortran array indexing is off by one. 
  }
#ifdef USE_MPI
  MPI_Bcast(Hsc.natex, ex_accum, MPI_INT, 0, Hsc.ourComm );  
#endif

  //set mixing parameters
  //Hsc.amberLambda = *amber_lambda;
  //Hsc.restLambda  = *restraint_lambda;
  //Hsc.dL_last         = 0.0;
  

  //init/setup array of atoms managed by this task:
  // this list will probably be modified later,
  // in emil_frc() for MPI version of the code.
  Hsc.myAtomIds = new int[1 + *nAtoms / Hsc.numTasks];
  {int i, ii, atsPerTask, atsLeft, myMin, myMax;
    ii         =  0;
    atsPerTask = *nAtoms / Hsc.numTasks;
    atsLeft    = *nAtoms % Hsc.numTasks;
  
    if( Hsc.myTaskId <= atsLeft ) {
        myMin = Hsc.myTaskId * (atsPerTask + 1);
    }else{
        myMin = atsLeft * (atsPerTask + 1) + (Hsc.myTaskId - atsLeft) * atsPerTask;
    }
    if( Hsc.myTaskId < atsLeft ) {
        myMax = myMin + atsPerTask + 1;
    }else{
        myMax = myMin + atsPerTask;
    }

    for(i = myMin; i < myMax; i++ ){
        Hsc.myAtomIds[ii++] = i;
    }
    Hsc.myNatoms=ii;
  }

  //log the simulation temperature
  Hsc.beta = *beta;

  // read simulation parameters from ParameterFile (can override passed-in lambdas)
  Hsc.readParamFile();

#ifdef USE_MPI
  {
   double bbox[3];
   if( *myTaskId == 0 ) memcpy( bbox, box, 3*sizeof(double));
   MPI_Bcast(bbox, 3, MPI_DOUBLE, 0, Hsc.ourComm);
   Hsc.setUpNew( crd, frc, vel, bbox );
  } 
#else
  // do the initialisation
  Hsc.setUpNew( crd, frc, vel, box );
#endif
  

  //move the restraints onto the particles, so that shape but 
  //not position or orientation are restrained.
#if 0
  Hsc.alignRestraints();//currently this function does nothing.
#endif

  *(Hsc.logFile) << "# Energy units: initial beta  " 
                 << Hsc.beta << " == 1.0 / " << 1.0/Hsc.beta << endl;
 *(Hsc.logFile)  << "# Scaling forces using polynomial mixing functions, at lambda = " 
                 << Hsc.restLambda << " :" << endl;

  if( *amber_softcoring ){

 *(Hsc.logFile) << "# accepting amber forces pre-scaled and softcored."  << endl; 
 *(Hsc.logFile) << "# not adding local 'soft force'."  << endl; 

  }else{

 *(Hsc.logFile) << "# scaling amber forces by " 
                << Hsc.mixFunc_mol( 1.0 - Hsc.amberLambda ) << endl; 
 *(Hsc.logFile) << "# scaling soft force by   " 
                << Hsc.mixFunc_soft( Hsc.restLambda )  << endl; 

  }
 *(Hsc.logFile) << "# scaling EMIL forces by  " 
                << Hsc.mixFunc_abs( Hsc.restLambda )   << endl; 


  Hsc.flushLog();

  Hsc.abstractEnergyZero = 0.0; //the energy zero is arbitrary
  Hsc.abstractEnergy     = Hsc.computePhiOfConfig();
 *(Hsc.logFile) << "# Uncorrected EABS is: " << Hsc.abstractEnergy  << endl;

  Hsc.flushLog();

  Hsc.setAbstractPotentialOffset( 0.0, Hsc.abstractEnergy );
  Hsc.abstractEnergy     = Hsc.computePhiOfConfig();
 *(Hsc.logFile) << "# With offset of " << Hsc.abstractEnergyZero << " this is: " << Hsc.abstractEnergy  << endl;
  
  //log the energy of the abstract system
  Hsc.reportSystemFreeEnergy();
  
  //share a tuning parameter for the softcore interaction
 *emil_softcore_scale =  Hsc.softcore_scale;
  
  //(over-generously) allocate storage for Hastings Moves.
  Hsc.oldConfBuf        = new double[3*Hsc.N];
  Hsc.logHastingsRatio  = new double[Hsc.N];
  Hsc.movedP            = new Particle*[Hsc.N];
  Hsc.nHastings         = 0;
  Hsc.nPart_hastings    = 0;
  Hsc.swapAccept        = 0;
  Hsc.swapMoves         = 0;
  Hsc.relocAccept       = 0;
  Hsc.relocMoves        = 0;
 
  
 *(Hsc.logFile) << "# emil mol_lambda: "  << Hsc.amberLambda;
 *(Hsc.logFile) << " abstract_lambda: " << Hsc.restLambda << endl;
  
  Hsc.reportSystemEnthalpy();
  Hsc.flushLog();
  
  
#ifdef SSTI_DEBUG_NONBONDS
  //debug
  Hsc.dbBuf  = new double[3*Hsc.N];
  
  
  
  {
    int    i,j, ii, pair;
    double R[3], r;
    
    dbfile_atomTypes.open("emil_debugAtomTypes.dat");
    dbfile_pairTypes.open("emil_debugPairTypes.dat");
    dbfile_atPairVDW.open("emil_debugPairVDW.dat");
 
    
    dbfile_pairTypes << "##i aCoeff[i] bCoeff[i] LJmin LJmin_SC" << endl;
    
    //save the atom VDW types to a logfile.
    for( i = 0; i < Hsc.N; i++ ){
      dbfile_atomTypes << i << " " <<  Hsc.LJ_atomTypes[i] << endl;
    }
 
    
     //save data about each pairwise interaction
     nPairTypes = 0;
     for( i = 0; i <  Hsc.N; i++ ){
       for(  j = 0; j < Hsc.N; j++ ){
          pair = (Hsc.LJ_atomTypes[i] - 1) * Hsc.nLJ_atomTypes + Hsc.LJ_atomTypes[j] - 1;
          if( Hsc.LJ_pairTypes[pair] > nPairTypes ){
            nPairTypes = Hsc.LJ_pairTypes[pair];
          }
       }
     }  
      
     cerr << nPairTypes << " distinct atomType pairs == " << sqrt(nPairTypes) << "**2" << endl;
     
     for(pair = 0; pair < nPairTypes; pair++ ){
        double rPrime, epsilon;
        double LJMin, LJMin_sc, wArg; 
        
        LJMin    = pow(2.0 * Hsc.LJ_bCoeff[pair]/ Hsc.LJ_aCoeff[pair], 0.166666666 );
        
        //epsilon =  softcore_scale * Hsc.scEpsilon( pair + 1 ) / ;
        epsilon  = LJMin * 0.5;
        
        
        wArg     = -1.0 * Hsc.amberLambda * exp(-LJMin/epsilon);
        LJMin_sc = epsilon * LambertW<0>( wArg ) + LJMin;
        
        dbfile_pairTypes << pair 
                         << " " << Hsc.LJ_aCoeff[pair] 
                         << " " << Hsc.LJ_bCoeff[pair] 
                         << " " << LJMin 
                         << " " << LJMin_sc
                         << " " << epsilon
                         << " " << wArg << endl;
                         
     } 
        
    // Log a lot of atom parameters to make sure that they are all being passed in correctly. 
    for( i = 0; i < 100; i++ ){
      for( j = 0; j < 100; j++ ){
  
        double enb, ee, evdw;
  
        r = 0.0;
        for( ii = 0; ii < 3; ii++ ){
          R[ii]  = ( Hsc.part[i].R[ii] - Hsc.part[j].R[ii] ) * ( Hsc.part[i].R[ii] - Hsc.part[j].R[ii] );
          if( Hsc.periodic ){
              while( fabs(R[ii]) > Hsc.box->halfx[ii] ){
	               R[ii] -= copysign(Hsc.box->x[ii], R[ii]);
	          }
	      }
          r += ( R[ii] * R[ii] );
        }
        r = sqrt( r );
        
                         
        ee   = 0.0;
        evdw = 0.0;
        enb  = Hsc.pairEnb( &(Hsc.part[i]), &(Hsc.part[j]), &ee, &evdw );
                        
        dbfile_atPairVDW << i << " " << j << " " <<  r << " "  <<  charge[i] << " " <<  charge[j] 
                         << " " <<  enb << " " << evdw << " " << ee << endl;
      }
    }
    
    dbfile_atomTypes.close();
    dbfile_pairTypes.close();
    dbfile_atPairVDW.close();
    
   }
#endif
  
  cerr << "# Task " << *myTaskId << " of " << *numTasks << " initted" << endl;
  
}

//function to do the abstract material forces, call every timestep
//of MD
void emil_forces(double *crd, double *frc, double *box, double *dvdl, 
                 double *molec_epot, double *ke, 
                 double *beta, int *frameCount, int *myAtomIds, 
                 int *myNatoms, int *doPme ){

  
  //just in case there has been a change in temperature
  Hsc.beta   = *beta;
  //Hsc.ekin   = *ke;
  Hsc.emolec = *molec_epot;
 
#ifdef USE_MPI
  if( !( *myNatoms == 0 && *myAtomIds == 0 ) ){ 
    if( pmemd_called == false ){
      //sander sends in the base atom as *myAtomIds, counting from 1. 
      Hsc.myNatoms  = *myNatoms - *myAtomIds + 1;
    
      //sanitize?
      if( Hsc.myNatoms <= 0 ){
        cerr << "Warning: task " << Hsc.myTaskId << " has atoms " << *myAtomIds - 1 << " to " << (Hsc.myNatoms - 1 + *myAtomIds - 1) << " inclusive, of " << Hsc.N << endl;
      }
        
  //sander uses a continuous block of atoms for each task
      delete[] Hsc.myAtomIds;
      Hsc.myAtomIds = new int[Hsc.myNatoms];
      for(int i = 0; i < Hsc.myNatoms; i++ ){
        Hsc.myAtomIds[i] = i + *myAtomIds - 1;//sander sends in the base atom as *myAtomIds, counting from 1. 
      }
    }else{
      
      //reassigning atom ids
      Hsc.myNatoms  = *myNatoms;
  
      delete[] Hsc.myAtomIds;
      Hsc.myAtomIds = new int[*myNatoms];
      for(int i = 0; i < *myNatoms; i++ ){
        Hsc.myAtomIds[i] = myAtomIds[i] - 1;//pmemd uses an array with a list of atom ids for each proc. 
      }
    
    }
  }
  
  if( pmemd_called )
      if( Hsc.myNatoms <= 0 ){
        cerr << "Warning: task " << Hsc.myTaskId << " has " << Hsc.myNatoms 
             << " ats, first @ " << Hsc.myAtomIds[0]  << " and last @ " << Hsc.myAtomIds[Hsc.myNatoms - 1] << endl;
        cerr << "at step " << *frameCount << endl; 
      }
#endif
 
   //refresh this local pointer in case the forces array has moved
   Hsc.forces = frc;
   
#if 0
   //move the restraints onto the particles, so that shape but 
   //not position or orientation are restrained.
   Hsc.alignRestraints();
#endif

   //book-keep changes in particle-well interaction
   Hsc.updateParticles();

   //compute the abstract forces locally
   Hsc.computeForces_mix();
  
   //get the abstract component of the energy... only call after calling computeForces_mix().
   Hsc.abstractEnergy = Hsc.computePhiOfConfig();
 
   //track the time
   Hsc.sim_nStep = *frameCount;

   //accumulate the generalised force for purposes of averaging
   //double precision == approx 10 significant figures.
   Hsc.gfAccum_soft += Hsc.genForce_soft( Hsc.softEnergy, Hsc.restLambda );
   Hsc.gfAccum_mol  += Hsc.genForce_mol( *dvdl, *molec_epot, 1.0 - Hsc.amberLambda );
   Hsc.gfAccum_abs  += Hsc.genForce_abs( Hsc.abstractEnergy, Hsc.restLambda );
   Hsc.gfAccum_count++;

   if( Hsc.sim_nStep % Hsc.printEvery == 0 ){
      if( Hsc.myTaskId == 0 ){

      char tmp_cstr[128];

  *(Hsc.logFile) << "nstep:  "  << Hsc.sim_nStep;

   /* some contortions here to get the format that I want without having to 
    * mess around with the stringstream class precision switches */
   sprintf(tmp_cstr,"%.8e",Hsc.gfAccum_soft/(double)Hsc.gfAccum_count);
  *(Hsc.logFile) << " soft_dHdL:  "  << tmp_cstr;

   sprintf(tmp_cstr,"%.8e",Hsc.gfAccum_mol/(double)Hsc.gfAccum_count);
  *(Hsc.logFile) << " molec_dHdL:  "  << tmp_cstr;

   sprintf(tmp_cstr,"%.8e",Hsc.gfAccum_abs/(double)Hsc.gfAccum_count);
  *(Hsc.logFile) << " abstr_dHdL:  "  << tmp_cstr << endl;

        //Hsc.reportSystemEnthalpy();
        Hsc.flushLog();
      }

     Hsc.gfAccum_soft  = 0.0;
     Hsc.gfAccum_mol   = 0.0;
     Hsc.gfAccum_abs   = 0.0;
     Hsc.gfAccum_count = 0;
   }   
    
   if( Hsc.fixLambda == false ){
      //check if is is time to increment lambda
      if( EXIT_SUCCESS ==   Hsc.stepLambda(Hsc.genForce_soft( Hsc.softEnergy, Hsc.restLambda )
                                        +  Hsc.genForce_abs( Hsc.abstractEnergy, Hsc.restLambda )  
                                        +  Hsc.genForce_soft( Hsc.softEnergy, Hsc.restLambda ) ) ){
        
        // save the well state
        if( Hsc.myTaskId == 0 ){
            Hsc.reportSystemEnthalpy();
            Hsc.writeOutWells();
            Hsc.flushLog();
            cerr << "Should be finishing on this step " << *frameCount << endl;
        }
       *doPme = 0;
      }
   }
    
    
   //let amber know if we currently need pme calc or not.
   if( Hsc.amberSoftcoring ){
      *doPme = 1;
   }else{
      if( Hsc.mixFunc_mol( 1.0 - Hsc.amberLambda ) != 0.0 ){
        *doPme = 1;
      }else{
        *doPme = 0;
      }
   }
  
    
   
   // save the well state
   if( Hsc.writeWells && (*frameCount % Hsc.saveWellsEvery == 0 && Hsc.myTaskId == 0) ){
      //*(Hsc.logFile) << "#writing wells at frame "<< *frameCount << " to " << Hsc.wellSaveFileName << endl;
      Hsc.writeOutWells();
   }
   //printFrame( Hsc.restraintCoords, &Hsc.N );

   
}

//function to attempt MC moves, and store partially-evaluated attempts in hastings buffers.
void emil_mcAttempt(double *crd, double *frc, double *vel, double *box, double *epot, int *emil_mc_done, double *beta, double *evdw_old, double *eelec_old){

  
  //if this is the first of several calls
  if( Hsc.relocMoves  == 0 ){
    Hsc.eMolecular_old  = *epot; //save the PE in units of kT
    Hsc.eVdwMolec_old   = *evdw_old; //save these for debug/feedback purposes
    Hsc.eElecMolec_old  = *eelec_old;
  }
  
  
  //For now assume that the coordinate array HAS NOT moved in memory
  //Hsc.particleCoords    = crd; //its too complicated to do an update properly unless absolutely required.
  //for( int i = 0; i < Hsc.N; i++ ){
  //  Hsc.part[i].R = &(crd[3*i]);
  //}
  
  //book-keep changes in particle-well interaction
  Hsc.updateParticles();
  
  //get the abstract component of the energy 
  Hsc.abstractEnergy     = Hsc.computePhiOfConfig();
  
  //attempt reassignment, and also exchange & relocation MC moves 
  Hsc.MC();
    
  if( Hsc.nHastings > 0 ){
    *emil_mc_done = 0;
  }else{
    *emil_mc_done = 1;
  }
  
  if( Hsc.nReloc > 0 && 
      Hsc.sim_nStep % Hsc.printEvery == 0 && 
      Hsc.mixFunc_mol(1.0 - Hsc.amberLambda) == 0.0)
    *(Hsc.logFile) << "Reloc acceptance rate: " <<  Hsc.hastingsFinalAccRate->mean << " ";

  if( Hsc.nSwap  > 0 && Hsc.sim_nStep % Hsc.printEvery == 0 )
    *(Hsc.logFile) << "Swap acceptance rate: " <<  Hsc.swapAccRate->mean  << "\n";

  if( Hsc.sim_nStep % Hsc.printEvery == 0 ){

    char tmp_cstr[128];

    for( int i = 0; i < Hsc.numLiquids; i++){
       if(Hsc.linearAssignment[i] != NULL){

          sprintf(tmp_cstr, "Reassignment rate %2i is: %.3f for d(sum(r^2)): %2.1e\n", i,
                      Hsc.linearAssignment[i]->assignmentRate->mean,
                      Hsc.linearAssignment[i]->assignmentDeRate->mean);
          *(Hsc.logFile) << tmp_cstr;

          Hsc.linearAssignment[i]->assignmentRate->clean();
          Hsc.linearAssignment[i]->assignmentDeRate->clean();
          
       }
    }
  }
  
}
 
//function which decides finally whether to accept the metropolis-hastings provisional moves
void emil_mcAcceptReject(double *crd, 
                          double *frc, 
                          double *vel, 
                          double *box, 
                          double *epot, 
                          double *beta, 
                          int    *emil_atoms_reverted, 
                          int    *emil_n_atoms_reverted, 
                          double *evdw, 
                          double *eelec){
  
  double    log_pAccept, delta_ePot;
  int       revCount;
  
  Hsc.beta = *beta;
  revCount = 0;
  
  delta_ePot   = *epot - Hsc.eMolecular_old;
  
  //in Metropolis-Hastings, final acceptance rate is premultiplied by the proposal distribution
  log_pAccept  = -1.0 * Hsc.logHastingsRatio[Hsc.nHastings - 1];//high proposal rate earlier == low acceptance rate now
  
  //this term represents the true Hamiltonian, including the expensive amber forcefield.
  log_pAccept -=  delta_ePot * Hsc.mixFunc_mol(1.0 - Hsc.amberLambda) * Hsc.beta;
  
  //this term represents the parts of the initial hastings ratio which were not guesses of the amber Hamiltonian,
  //but instead had a physical meaning of their own.
  log_pAccept -=  Hsc.deAbstract * Hsc.mixFunc_abs(Hsc.restLambda);
    
  if( log(Hsc.mRan())  <  log_pAccept || Hsc.nHastings == 0 ) {
     
     *(Hsc.logFile) << "accepted a " << Hsc.nHastings << "-ple reloc move with prob of " << exp( log_pAccept ) << endl;
     *(Hsc.logFile) << "DE_AMBER:    " << (*epot - Hsc.eMolecular_old)  << " DE_VDW_AMBER:    " << (*evdw - Hsc.eVdwMolec_old) 
          << " DE_ELEC_AMBER:    " << (*eelec - Hsc.eElecMolec_old) << endl; 
     *(Hsc.logFile) << "DE_HASTINGS: " << Hsc.deHastings_estimated      << " DE_VDW_HASTINGS: " << Hsc.deVdwHastings  
          << " DE_ELEC_HASTINGS: " << Hsc.deElecHastings << endl;
     *(Hsc.logFile) << "DE_ABSTRACT: " << (Hsc.deAbstract / Hsc.beta)   << " PREV_HASTINGS:   " << exp(Hsc.logHastingsRatio[Hsc.nHastings - 1])   
          << " LxDE_AMBER+DE_ABSTRACT: " << delta_ePot* Hsc.mixFunc_mol(1.0 - Hsc.amberLambda) + Hsc.deAbstract * Hsc.mixFunc_abs(Hsc.restLambda) / Hsc.beta << endl;
    
     Hsc.eMolecular_old  = *epot;
     Hsc.eVdwMolec_old   = *evdw;
     Hsc.eElecMolec_old  = *eelec;
     Hsc.acceptHastings();
   
     
  }else{ //delete one of the move attempts, ready to try again.
    
     *(Hsc.logFile) << "rejected a reloc move with prob of 1 - " << exp( log_pAccept ) << endl;
     *(Hsc.logFile) << "DE_AMBER:    " << (*epot - Hsc.eMolecular_old)  << " DE_VDW_AMBER:    " << (*evdw - Hsc.eVdwMolec_old) 
          << " DE_ELEC_AMBER:    " << (*eelec - Hsc.eElecMolec_old) << endl; 
     *(Hsc.logFile) << "DE_HASTINGS: " << Hsc.deHastings_estimated      << " DE_VDW_HASTINGS: " << Hsc.deVdwHastings  
          << " DE_ELEC_HASTINGS: " << Hsc.deElecHastings << endl;
     *(Hsc.logFile) << "DE_ABSTRACT: " << (Hsc.deAbstract / Hsc.beta)   << " PREV_HASTINGS:   " << exp(Hsc.logHastingsRatio[Hsc.nHastings - 1])
          << " LxDE_AMBER+DE_ABSTRACT: " << delta_ePot* Hsc.mixFunc_mol(1.0 - Hsc.amberLambda) + Hsc.deAbstract * Hsc.mixFunc_abs(Hsc.restLambda) / Hsc.beta << endl;
     
     revCount = Hsc.rejectHastings( log_pAccept );
     
  
  }
  
  *(Hsc.logFile) << "returning " << revCount << " atoms reverted" << endl;
  *emil_n_atoms_reverted = revCount;
  
}






