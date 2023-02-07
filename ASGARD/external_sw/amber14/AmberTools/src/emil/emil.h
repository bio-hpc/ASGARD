/*************************************************************************
 * TI code
 *************************************************************************/

#ifndef HAVE_EMIL_H
#define HAVE_EMIL_H

#include "mtrand.h"
#include "time.h"

//////////////////////////////////////////////////////////////
//Turn off this flag if you want to use explicit water, 
//as in Berryman & Schilling J Comput Chem 2013.
//
//The expense of this type of calculation is such,
//and the difficulty of tuning parameters such that 
//the simulation conserves energy,
//that explicit-water EMIL is not supported for the current 
//release of AMBER.
//
//If you really want to do this and you find that you would 
//like some help, then contact josh.berryman@uni.lu.
//
#define _EMIL_REFUSE_WATER
//////////////////////////////////////////////////////////////

class Particle;
class Well;
class coneWell;
class harmonicWell;
class Cell;
class Box;
class Liquid;
class LinearAssignment;

//some tuning parameters and compile switches
#define RESOLVE_SWAPS_EVERY (1024)
#define SSTI_BUFFER_LOG        
#define SSTI_SOFT_REPULSION   
#define EMIL_RESTRAIN_HATOMS 1
#ifdef SSTI_BUFFER_LOG
#  include <cstdio> // Certain compilers need this for FILE
#endif
//These hash-deffs define a fixed range SF_MAX_R for the 
//soft force between *all* pairs of particles.
#define SF_MAX_R2     16.0
#define SF_MAX_R       4.0

//define a macro for stderr output from master only.
//bit of a dirty hack really.
#ifdef USE_MPI
#define MERR if(myTaskId == 0) cerr
#else
#define MERR                   cerr 
#endif

//define a macro for debugWrap
#define PCHECK(p,m) if(p==0){MERR << "pointer is zero " << m << endl;}

////default water parameters: wu, tepper & voth "SPCFw"
//#define R_HWATO_SPCF (1.012)      //from wu tepper voth
//#define R_HWATO_SPCF (1.03398)    //observed value
#define R_HWATO_SPCF (1.0327)   //this value gives roughly the observed geometry

//#define THETA_0_SPCF (1.97641084) //radians:= pi * 113.24 degrees/180 
#define THETA_0_SPCF (1.873) //again, this value roughly reproduces the observed geometry


#define K_A_SPCF     (75.90)     //kcal/mol/rad^2
#define K_B_SPCF     (1059.162)   //kcal/mol/Angstrom^2

//in case I ever work out a triangular implementation
#define R_HWATH_SPCF (1.633)

#define DOT( a, b ) (a[0]*b[0]+a[1]*b[1]+a[2]*b[2])


/* define a flag to select different restraining potentials */
typedef enum{ HSC_POTENTIAL_NULL             = 0,
              HSC_POTENTIAL_EINSTEIN_CRYSTAL = 1,
              HSC_POTENTIAL_CONEWELL_LIQUID  = 2,
              HSC_POTENTIAL_HARMONIC_LIQUID  = 3,
              HSC_POTENTIAL_WINGWELL_LIQUID  = 4,
              HSC_POTENTIAL_WOBBLIUM         = 5,
              HSC_POTENTIAL_RELATIVE_POS     = 6,
              HSC_POTENTIAL_NUMPOTS          = 7
} HSC_POTENTIAL_T;
#define HSC_POTENTIAL_NAME(x) (x==HSC_POTENTIAL_NULL?            "Null             ":\
                               x==HSC_POTENTIAL_EINSTEIN_CRYSTAL?"Einstein Molecule":\
                               x==HSC_POTENTIAL_CONEWELL_LIQUID? "Conewell Liquid  ":\
                               x==HSC_POTENTIAL_HARMONIC_LIQUID? "Harmonic Liquid  ":\
                               x==HSC_POTENTIAL_WINGWELL_LIQUID? "Wingwell Liquid  ":\
                               x==HSC_POTENTIAL_WOBBLIUM?        "Wobblium         ":\
                               x==HSC_POTENTIAL_RELATIVE_POS?    "Relative position":\
                                                               "Error: unrecognised")

/* define a flag to select different means of solving/sampling the 
   particle-well assignment space, for liquids where the well-particle 
   relationship is not fixed.  */
typedef enum{ HSC_ASSIGNMENT_MC    = 0,
              HSC_ASSIGNMENT_JONKER_VOLGENANT  = 1,
              HSC_ASSIGNMENT_NUMASSIGNMENTS    = 2
} HSC_ASSIGNMENT_T;

#define HSC_ASSIGNMENT_NAME(x) \
           (x==HSC_ASSIGNMENT_MC?               "Monte Carlo      ":\
            x==HSC_ASSIGNMENT_JONKER_VOLGENANT? "Jonker Volgenant ":\
                                                  "Error: unrecognised")





class Box {

 public:

  double x[3], halfx[3], halfxInv[3];  // dimensions
  double V;      // Volume
  double xCell[3]; //Cell dimensions
  int    nxCell[3]; // number of cells
  int    nCells;
  double xCellW[3]; //Cell dimensions wrt rcut
  int    nxCellW[3]; // number of cells wrt rcut
  int    nCellsW;

  void update(){
      for (int i=0; i<3; i++){
	  halfx[i] = 0.5*x[i];
	  halfxInv[i] = 1.0/halfx[i];
      }
      V = x[0]*x[1]*x[2];
  }

  //return coordinate inside box
  inline double imageR( double r, int i ){
    while( r >= halfx[i] ){
      r -= x[i];
    }
    while( r < -1.0 * halfx[i] ){
      r += x[i];
    }
    return( r );
  }

};

#include "wells.h" //wells depend on box, so must come later in this header

using namespace std;


class Parms {
 public:
  int     nAtoms;
  int     nAtomTypes;
  int     nRes;
  int     nMolecules;
  int    *resFirstAtom;
  int    *atomTypes;
  int    *pairTypes;
  char   *atNames;
  char   *resNames;
  int    *molSizes;
  double *atMasses;
  double *atCharges;
};

class ErrorState {

  public:
    void clean(){
      mean = 0.0;
      SD   = 0.0;
      ESE  = 0.0;
      n    = 0;
      sumX = 0.0; sumXX = 0.0; nSn = 0.0;
    }
    ErrorState(void){
      clean();
    }
    
    
    double accumulate( double x ){
      n++;
      sumX   += x;
      sumXX  += x * x;
      mean   += (x - mean) / double(n);
      nSn     = n * mean * mean + sumXX - 2.0 * mean * sumX;
      SD      = sqrt(nSn / double(n));
      if( n > 1 ){
        ESE     = sqrt( nSn / double(n*(n-1)) );
      }else{
        ESE     = 0.0;//no valid error at n=1
      }
      return( ESE );
    };
    
    double             mean, SD, ESE;
    unsigned long int  n;
    
  private:
    double             sumX, sumXX, nSn;
    
};




class hsc {

 public:
  
  //Mersenne twister RNG.
  MTRand mRan;

  //Accumulators to prepare averages of generalised force.
  //Print every hundred or so and there should normally be no overflow.
  unsigned int gfAccum_count;
  double       gfAccum_soft, gfAccum_mol, gfAccum_abs;
  
  //some flags
  bool newSetUp, stepLambdaBackwards, fixLambda;
  bool readWells, writeWells;
  bool amberSoftcoring, emilSoftForce;

  string      topFileName, crdFileName, pdbFileName; //names of molecular description files
  string      trajOutFileName, crdOutFileName; //names of molecular description files
  string      wellFileName;                    //load start well positions from this file
  string      wellSaveFileName;                //save  well positions to this file 
  string      logFileName;                     //save progress to this logfile
  string      parameterFile, configFile;       // Input filenames
 #ifdef SSTI_BUFFER_LOG
  stringstream *logFile;
  FILE         *logFile_fp;
 #else
  ostream      *logFile;
 #endif
  
  int         saveWellsEvery, printEvery;
  
  int         myTaskId;
  int         numTasks;
  int        *myAtomIds, myNatoms;
 #ifdef USE_MPI
  MPI_Comm    ourComm;
  int        *blockDisplace, *blockLength;
 #endif

  void run();

  double *particleCoords;
  double *restraintCoords;
  double *forces;
  double *velocities;
  double *workspace;

  //data describing the integration
  ErrorState eState_dHdL;
  double     dHdL_last;
  int        integrationStepCount;
  int        nStepsLimit;
  bool       integrationFinished;
  bool       periodic;  

  //for benchmarking/progress estimation
  time_t  stepStartTime;
  
  //data describing external simulation program
  double  amberLambda, restLambda;
  double  beta; // beta = 1/kT
  int    *LJ_atomTypes, *LJ_pairTypes, nLJ_atomTypes;
  double *atomCharge, *LJ_aCoeff, *LJ_bCoeff; //defined pairwise over atoms: arrays of size nLJ_atomTypes**2.
  double  permittivity;         //permittivity of the medium; used for Coulomb force.   
  double  softcore_scale;       //parameter for extent of softcoring. 1.0 is maximum.   

  double  scEpsilon(int pairType); //get the softcoring parameter for this atom pair

  Particle *part;
  Well    **wells;
  Cell     *cells;           // cells for overlap detection
  Cell     *cellsW;          // cells for detection of particles in attr. wells
  Box      *box;

  //part-well assignment solvers (one per degenerate liquid)
  LinearAssignment **linearAssignment; 

  // Simulation Data


  unsigned long int sim_nStep;
  ErrorState       *hastingsPropRate, *hastingsFinalAccRate, *swapAccRate;

  HSC_ASSIGNMENT_T  assignmentMethod; //how to match liquid particles to their wells.

  int transAccept;
  int transMoves;
  int swapAccept;
  int swapMoves;
  int relocAccept;
  int relocMoves;
 
  //staging-space for resolving sap collisions
  int **allSwaps;
  int  *numSwaps;
  
  int   runCount;
  int   measureCount;
  long  seed;          //seed for ran3
  long *iran;         //pointer to seed

  int eqSteps;          
  int eqEvalSteps;      // # trajectory averages taken during equilibration.
  int eqSnapSteps;      // # configuration snapshots taken during equil.
  int mcSteps;
  int mcEvalSteps;      // # measurements for averages
  int mcSnapSteps;      // # snapshots during mc

  double maxDisplace;     

  double bx;         // box dimensions
  int    N;          // particle number
  double invNheavy_solid;

  int               numMaterials;     // number of different types of material eg Einstein Crysta, liquid 1 (water), liquid 2 (salt)
  HSC_POTENTIAL_T  *materialType;     // potential applied to each material present
  string           *materialAtomMask; // regular-expressions which identify which atoms belong to each material for nab.
  
  Particle       ***liquidLists;      // array of molecules of each type of liquid.
  int              *chainsInLiquid;

  int               numLiquids, numLiquidChains;
  
  int              *numex, *natex;    //array of number of excluded nb interactions, and of specific atoms to exclude

  //bookkeeping for Hastings-Metropolis
  double          *oldConfBuf;
  double          *logHastingsRatio; 
  Particle       **movedP;
  int              nHastings, nPart_hastings;
  double           eMolecular_old, deHastings_estimated, deAbstract;
  double           deVdwHastings, deElecHastings; //save these for debug/feedback purposes.
  double           eVdwMolec_old, eElecMolec_old; //save these for debug/feedback purposes.
  int              nMcMoves, nSwap, nReloc;
  double           swapTriesPerChain, relocTriesPerChain;
  double           targetSwapRate, targetRelocRate;
  int              assignmentEvery;  

  Parms       *parms;

  double generalisedForce; //AMBER generalised force dV/dl
  double abstractEnergy, abstractSolidEnergy, abstractEnergyZero, softEnergy;   //energy of the abstract system
  double ekin, emolec, workInThermostat;

  double epsilonBond, epsilonTriangle, epsilonTrap, epsilonSolid, wingForce; //coupling strength
  double rcutLiquid, rcut2Liquid, rcut3Liquid, dCutLiquid, 
         rcutSolid, req, rcut_longest, rcut_longest2;  //cutoff of range of attraction to wells
   
  // observables
  // potential energy with respect to attractive wells
  double meanPhi;  
  double meanDistanceFromWell;
  double meanRatioOutside;
  
  //routines
  //set up
  void setUpNew(double *crd, double *frc, double *vel, double *box);
  void setUpFromFile();
  void setUpWells();
  void initialize();
  void initWells();
  void randomiseWells();
  void setWellTypes();
  void readParamFile();
  

  void setupList(); 
  void setupListW();
  void setUpNeighbours();
  void setUpNeighboursW();

  // simulation

  void equil();
  void MC();

  void   swapMove();
  void   reduceSwaps();
  void   tryWellMove();
  double tryRelocMove();
  double nonHastingsReloc();
  void   saveRelocAttempt( Particle *chainP, double *chainCoords, double log_pAccept );
  void   adjustAttemptRates();
  void   acceptHastings();
  int    rejectHastings(double log_pAccept);
  void   alignRestraints();
  void   boostKEOnWellPlace();
  
  // manage AMBER molecule
  void   updateParticles(); //reset to take account of changes introduced by atomistic MD
  string solidMask, liquidMask;
  int   *soluteMask, topSoluteAtom;
  
  //temp funcs for debug
  void debugWrap(int i, int newCell);
  void reportWellPart();
  void writeVtf();

  //energies and forces
  double computePhiOfConfig();
  void   computeForces();
  void   computeForces_mix();
  int    stepLambda( double dHdL );
  double getLStep( double dHdL_mean );
  double mixFunc_mol( double lambda );
  double mixFunc_abs( double lambda );
  double mixFunc_soft( double lambda );
  double genForce_mol( double dvdl, double emodel, double lambdaMix);
  double genForce_abs( double emodel, double lambdaMix);
  double genForce_soft( double soft, double lambdaMix);
  double computeDistance(Particle &p);
  double computeDistance(Particle &p, Well *w);
  int    selectCell(double R[3], double xcellInv[3]); //which cell is point R in?
  int    selectCellW(double R[3], double xcellWInv[3]); //which cell is point R in?
  Particle *pickLiquidChain( int *lIndex ); //root of random liquid molecule, also saves index of the liquid to *lIndex.
  double setAbstractPotentialOffset( double eMol, double eAbs );
  double computeAverageDistanceFromWell();
  double computeRatioOutside();
  double reportSystemFreeEnergy();
  double reportSystemEnthalpy();
  double reportParticleKE(int i);
  double logFac(int x){ //calculate log(x!)
    double lnFact;
    lnFact = 0.0;
    while( x > 1 ){
      lnFact += log( (double) x );
      x--;
    }
    return( lnFact );
  };
  double    softForce(Particle *p1, Particle *p2, double *fij, double scaleSoft );
  double    softForceAll(Particle *p, double scaleSoft, double *fij );
  double    getSoftEnergy(Particle *p1, Particle *p2 );
  double    getSoftEnergyAll(Particle *p );
  double    sc_ljA(int pairType);
  
  //molecular mechanics
  double sigmoid(double r, double a, double b);
  double chainEnergy(Particle *p, double *dist);
  double wholeSystemEE( Particle *p );
  double chainEnb( Particle *p, double *eelec );  //intermolecular van der Waals + Coulomb energy of the chain.
  double pairEnb( Particle *p, Particle *q, double *eelec, double *evdw );  //pairwise evdw + coulomb.
  double addBond( double *forceVec, Particle *p, Particle *bondP, double rMin, double epsilon );
  double addBonds_spcfW( Particle *oxy, double *forceOxy, double *forceH1, double *forceH2, double mix );
  double addBonds_spcfWater( Particle *oxy, double *forceOxy, double *forceH1, double *forceH2, double mix );
  
  double rCutVdw;
  
  //MC
  Particle *pickLiquidChain(int *liq, int *size);
  void      displaceParticle(double newP[3], double oldP[3]);
  void      relocateParticle(Particle &p);
  void      relocateWell(Well *w);
  void      relocateChain(Particle *p);
  void      replaceTrialReloc( Particle *q, double *invxCell, double *rOld, int *oldCellIndex );
  bool      checkOverlapAll(Particle &p);
  bool      checkOverlapTwo(Particle &p1, Particle &p2);
  bool      checkOverlapSoft(Particle &p1, Particle &p2);
  bool      testChainOverlap(Particle *p);
  bool      testWellOverlap( Well *w, Cell *c );
  
  // io
  void writeOutWells();
  void writeOutConfig(char name[16]);
  void writeOutSimData(char name[16]);
  void writeOutBinding(char name[16]);
  void writeOutBindingToWells(char name[16]);
  void flushLog();
  
  //debug
  double *dbBuf;
  void testLoop( Particle *p );
  void printAverageWaterBonds();
  
  //util
  void randomVector(double v[3]);
};


class Cell {

 public:

  Particle *firstParticle;
  Well     *firstWell;
  Cell     *neighbours[26];      // neighbouring cells 
  int       myId;
  Cell(void){
    firstParticle = 0;
    firstWell     = 0;
    for (int i=0; i<26; i++){
      neighbours[i] = 0;
    }
  }
  

};





class Particle {

 public:
  
  double           *R;           // center of mass
  Well             *theWell;     // well it is attracted to
  HSC_POTENTIAL_T   wellType;    // the type of well which that is
  int               liquidIndex; // the type of liquid (if any) which the particle is
  int               chainIndex;  // the number of the chain (of all liquid chains) which the particle is in
  double            rijsqW;      // squared distance from its well of attraction
                                 // (updated after each translation and swap).... for non-root liquid particles, is instead the particle's potential energy.
  double            rij[3];        // vector to well, accounting for PBCs
  
  
  Particle         *next, *prev; // next, previous in cell list
  Cell             *cell;        // cell with respect to overlap radius
  Cell             *WCell;       // cell with respect to rcut
  Particle         *nextWC, *prevWC; //next, previous in cell of rcut
  Particle         *nextW, *prevW;   //linked list of well which is checked for swap
  Particle         *strandAtom;      //list of particles in same covalently bonded molecule - only set up for liquid mols.
  bool              isRoot;          //flag particles which are root of a chain
  int               myId;

  Particle(){
      next        = 0;
      prev        = 0;
      cell        = 0;
      nextW       = 0;
      prevW       = 0;
      nextWC      = 0;
      prevWC      = 0;
      WCell       = 0;
      theWell     = 0;
      wellType    = HSC_POTENTIAL_NUMPOTS;
      liquidIndex = -1;
      chainIndex  = -1;
      strandAtom  = 0;
      isRoot      = true;
      myId        = 0;
  }
  
  void insertToCell(Cell &c){

    if( c.firstParticle == this ){
      cerr << "Error, double-insertion" << endl;
      exit( 8 );
    }
    
    
    prev = 0;               // initialize pointer
    next = c.firstParticle; // -> firstParticle must be initialized to 0!
    if (next) next->prev=this;
    c.firstParticle=this;
    cell=&c;
    
  };

  void moveBetweenCells(Cell &from, Cell &to){
    
    // first remove the particle from the old cell
    if (prev) prev->next=next;
    else from.firstParticle=next;
    if (next) next->prev=prev;

    // now add particle
    prev=0;
    next=to.firstParticle;
    if (next) next->prev=this;
    to.firstParticle=this;
    cell=&to;
    
  };
  void moveBetweenCells(Cell *from, Cell *to){

    // first remove the particle from the old cell
    if (prev) prev->next=next;
    else from->firstParticle=next;
    if (next) next->prev=prev;

    // now add particle
    prev=0;
    next=to->firstParticle;
    if (next) next->prev=this;
    to->firstParticle=this;
    cell=to;
  };

  void insertToCellW(Cell &c){

    prevWC = 0;               // initialize pointer
    nextWC = c.firstParticle; // -> firstParticle must be initialized to 0!
    if (nextWC) nextWC->prevWC=this;
    c.firstParticle=this;
    WCell=&c;
    
  };

  void moveBetweenCellsW(Cell &from, Cell &to){

    // first remove the particle from the old cell
    if (prevWC) prevWC->nextWC=nextWC;
    else from.firstParticle=nextWC;
    if (nextWC) nextWC->prevWC=prevWC;

    // now add particle
    prevWC=0;
    nextWC=to.firstParticle;
    if (nextWC) nextWC->prevWC=this;
    to.firstParticle=this;
    WCell=&to;
  };


  void insertToWell(Well *w){

    prevW = 0;
    nextW = w->firstParticle; 
    if (nextW) nextW->prevW=this;
    w->firstParticle=this;
    w->npart++;
    
  };

};

//each item of class liquid is a collection of interchangeable molecules
class Liquid {

  //this is an array of pointers to the "base" heavy atom of each liquid molecule 
  Particle       **heavies;
  int              nHeavies;
  
  //this identifies the potential type of the liquid
  HSC_POTENTIAL_T  wellType;
  
  //this identifies atoms which belong in it
  string           atomMask;
  
};


#endif //HAVE_EMIL_H
























