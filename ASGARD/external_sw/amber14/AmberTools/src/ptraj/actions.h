/*  _______________________________________________________________________
 *
 *                        RDPARM/PTRAJ: 2008
 *  _______________________________________________________________________
 *
 *  This file is part of rdparm/ptraj.
 *
 *  rdparm/ptraj is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  rdparm/ptraj is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You can receive a copy of the GNU General Public License from
 *  http://www.gnu.org or by writing to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *  ________________________________________________________________________
 *
 *  CVS tracking:
 *
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/actions.h,v 10.2 2009/12/28 20:24:04 isjoung Exp $
 *
 *  Revision: $Revision: 10.2 $
 *  Date: $Date: 2009/12/28 20:24:04 $
 *  Last checked in by $Author: isjoung $
 *  ________________________________________________________________________
 *
 *
 *  CONTACT INFO: To learn who the code developers are, who to contact for
 *  more information, and to know what version of the code this is, refer
 *  to the CVS information and the include files (contributors.h && version.h)
 *
 */

#include "contributors.h"
#include "version.h"

/*  ________________________________________________________________________
 */



/*
 *  This is the header file for action.c which contains the basic structures
 *  necessary for implementing the ACTIONS called by ptraj().  In this file,
 *  the definitions are as follows:
 *
 * (1) EXTERNALLY VISIBLE DEFINITIONS
 * (2) GLOBAL VARIABLES
 * (3) EXTERNALLY VISIBLE FUNCTION PROTOTYPES
 * (4) LOCAL STRUCTURES
 */


/*
 * (1) EXTERNALLY VISIBLE DEFINITIONS
 */


   /*
    *  possible ptraj ACTIONS 
    */

typedef enum _actionType { 
  TRANSFORM_NOOP, 
  TRANSFORM_ACCEPTOR,
  TRANSFORM_ANALYZE,
  TRANSFORM_ANGLE,
  TRANSFORM_ATOMICFLUCT,
  TRANSFORM_ATOMICFLUCT3D,
  TRANSFORM_AVERAGE,
  TRANSFORM_BENCHMARK,
  TRANSFORM_BOX,
  TRANSFORM_CENTER,
  TRANSFORM_CHECKDNA,
  TRANSFORM_CHECKOVERLAP,
  TRANSFORM_CLOSESTWATERS,
  TRANSFORM_CLUSTER,
  TRANSFORM_CLUSTERATTRIBUTE,
  TRANSFORM_CONTACTS,
  TRANSFORM_DIFFUSION,
  TRANSFORM_DIHEDRAL,
  TRANSFORM_DIHEDRALCLUSTER,
  TRANSFORM_DIPOLE,
  TRANSFORM_DISTANCE,
  TRANSFORM_DONOR,
  TRANSFORM_DNAIONTRACKER,
  TRANSFORM_ECHO,
  TRANSFORM_ENERGY,
  TRANSFORM_GRID,
  TRANSFORM_HBOND,
  TRANSFORM_IMAGE,
  TRANSFORM_MATRIX,
  TRANSFORM_PRINCIPAL,
  TRANSFORM_PRNLEV,
  TRANSFORM_PROJECTION,
  TRANSFORM_PUCKER,
  TRANSFORM_RADIAL,
  TRANSFORM_RADIUSOFGYRATION,
  TRANSFORM_RANDOMIZEIONS,
  TRANSFORM_REFERENCE,
  TRANSFORM_RMS,
  TRANSFORM_RUNNINGAVERAGE,
  TRANSFORM_SCALE,
  TRANSFORM_SECONDARYSTRUCT,
  TRANSFORM_STRIP,
  TRANSFORM_SMARTIMAGE,
  TRANSFORM_SOLVENT,
  TRANSFORM_TRAJIN, 
  TRANSFORM_TRAJOUT, 
  TRANSFORM_TRANSFORM,
  TRANSFORM_TRANSLATE,
  TRANSFORM_TRUNCOCT,
  TRANSFORM_TEST,
  TRANSFORM_UNWRAP,
  TRANSFORM_VECTOR,
  TRANSFORM_CORRELATION,
  TRANSFORM_WATERSHELL,
  TRANSFORM_2DRMS
} actionType;


   /*
    *  ACTION FUNCTION TYPE DEFINITION
    *
    *    action->fxn(action, X, Y, Z, box, ptrajMode)
    */


#  ifdef __STDC__
typedef int (*actionFunction)( void *, double *, double *, double *, 
			       double *, int);
#  else
typedef int (*actionFunction)();
#  endif


   /*
    *  ACTION INFORMATION (passed in as the first argument 
    *  to the action function)
    */

typedef struct _actionInformation { 
  actionFunction fxn; 
  actionType type;
  ptrajState *state;
  int suppressProcessing;
  int performSecondPass;
  int *mask;
  int iarg1;
  int iarg2;
  int iarg3;
  int iarg4;
  int iarg5;
  int iarg6;
  int iarg7;
  int iarg8;
  double darg1;
  double darg2;
  double darg3;
  double darg4;
  double darg5;
  double darg6;
  double darg7;
  double darg8;
  void *carg1;
  void *carg2;
  void *carg3;
  void *carg4;
  void *carg5;
  void *carg6;
  void *carg7;
  void *carg8;
} actionInformation;

#define INITIALIZE_actionInformation(_p_) \
  _p_->fxn = NULL;            \
  _p_->type = TRANSFORM_NOOP; \
  _p_->state = NULL;          \
  _p_->suppressProcessing = 0;\
  _p_->performSecondPass = 0; \
  _p_->mask = NULL;           \
  _p_->iarg1 = 0;             \
  _p_->iarg2 = 0;             \
  _p_->iarg3 = 0;             \
  _p_->iarg4 = 0;             \
  _p_->iarg5 = 0;             \
  _p_->iarg6 = 0;             \
  _p_->iarg7 = 0;             \
  _p_->iarg8 = 0;             \
  _p_->darg1 = 0.0;           \
  _p_->darg2 = 0.0;           \
  _p_->darg3 = 0.0;           \
  _p_->darg4 = 0.0;           \
  _p_->darg5 = 0.0;           \
  _p_->darg6 = 0.0;           \
  _p_->darg7 = 0.0;           \
  _p_->darg8 = 0.0;           \
  _p_->carg1 = NULL;          \
  _p_->carg2 = NULL;          \
  _p_->carg3 = NULL;          \
  _p_->carg4 = NULL;          \
  _p_->carg5 = NULL;          \
  _p_->carg6 = NULL;          \
  _p_->carg7 = NULL;          \
  _p_->carg8 = NULL

typedef struct _trajectoryInfo {
  ptrajState *state;
  int atoms;
  int current;
  int allocated;
  int rollover;
  float *x;
  float *y;
  float *z;
} trajectoryInfo;

#define INITIALIZE_trajectoryInfo(_p_) \
  _p_->state = NULL; \
  _p_->atoms = 0; \
  _p_->current = 0; \
  _p_->allocated = 0; \
  _p_->rollover = 0; \
  _p_->x = NULL; \
  _p_->y = NULL; \
  _p_->z = NULL


typedef struct _windowInfoType {
  int windowWidth;
  int windowOffset;
  int windows;
  int *windowStart;
  int *windowStop;
  int *windowVisits;
  double *x;
  double *y;
  double *z;
} windowInfoType;


/*
 *  ACTION: TRANSFORM_CONTACTS  -- this structure needs to be global since it is
 *                                 accessed by ptrajSetupIO
 */

typedef enum _transformContactsType {
  CONTACTS_NULL = 0,
  CONTACTS_FIRST,
  CONTACTS_REFERENCE
} transformContactsType;

typedef struct _transformContactsInfo {
  double *refx, *refy, *refz;
  char *filename;
  void *outFile;
  void *outFile2;
  int currentFrame;
  int byResidue;
} transformContactsInfo;

#define INITIALIZE_transformContactsInfo(_p_) \
  _p_->refx      = NULL; \
  _p_->refy      = NULL; \
  _p_->refz      = NULL; \
  _p_->filename  = NULL; \
  _p_->outFile   = NULL; \
  _p_->outFile2  = NULL; \
  _p_->currentFrame = 0; \
  _p_->byResidue    = 0;


typedef struct _contactList *ptcontactList;

typedef struct _contactList {
  int index;
  char *name;
  struct _contactList *next;
} contactList;

#define INITIALIZE_contactList(_p_) \
  _p_->index = -1; \
  _p_->name  = NULL; \
  _p_->next  = NULL;


/*
 *  ACTION: TRANSFORM_SECONDARYSTRUCTURE
 *     -- this structure needs to be global since it is accessed by ptrajSetupIO
 */

typedef enum _transformSecStructType {
  SECSTRUCT_NULL = 0,
  SECSTRUCT_PARA,
  SECSTRUCT_ANTI,
  SECSTRUCT_3_10,
  SECSTRUCT_ALPHA,
  SECSTRUCT_PI,
  SECSTRUCT_TURN
} transformSecStructType;

typedef struct _transformSecondaryStructInfo {
  char *filename;
  FILE *outFile;
  int **hbondMatrix;
  int **turns;
  int **bridges;
  transformSecStructType *secStruct;
} transformSecondaryStructInfo;

#define INITIALIZE_transformSecondaryStructInfo(_p_) \
  _p_->filename  = NULL; \
  _p_->outFile   = NULL; \
  _p_->hbondMatrix= NULL; \
  _p_->turns     = NULL; \
  _p_->bridges   = NULL; \
  _p_->secStruct = NULL;

typedef struct _ResidueType {
  int C, N, H, O;
  int alphaHelixProb, threeHelixProb, piHelixProb;
  int parallelSheetProb, antiparallelSheetProb;
  int turnProb;
  int isActive;
} ResidueType;

#define INITIALIZE_ResidueType(_p_) \
  _p_.C = -1; \
  _p_.N = -1; \
  _p_.H = -1; \
  _p_.O = -1; \
  _p_.alphaHelixProb = 0; \
  _p_.piHelixProb = 0; \
  _p_.threeHelixProb = 0; \
  _p_.parallelSheetProb = 0; \
  _p_.antiparallelSheetProb = 0; \
  _p_.turnProb = 0; \
  _p_.isActive = 0;

/*
 *  ACTION: TRANSFORM_RADIUSOFGYRATION  
 *    -- this structure needs to be global since it is accessed by ptrajSetupIO
 */

typedef struct _transformRadiusOfGyrationInfo {
  char *filename;
  FILE *outFile;
  int currentFrame;
} transformRadiusOfGyrationInfo;

#define INITIALIZE_transformRadiusOfGyration(_p_) \
  _p_->filename  = NULL; \
  _p_->outFile   = NULL; \
  _p_->currentFrame = 0;



/*
 *  ACTION: TRANSFORM_RMS  -- this structure needs to be global since it is
 *                            accessed by ptrajSetupIO
 */

typedef struct _transformRmsInfo {
  double *refx, *refy, *refz;
  char *filename;
  char *name;
  double *rmsValues;
  int currentFrame;
} transformRmsInfo;

enum _rmsActions { RMS_NOOP, RMS_PREVIOUS, RMS_FIRST, 
		     RMS_REFERENCE };

#define INITIALIZE_transformRmsInfo(_p_) \
  _p_->refx      = NULL; \
  _p_->refy      = NULL; \
  _p_->refz      = NULL; \
  _p_->filename  = NULL; \
  _p_->name      = NULL; \
  _p_->rmsValues = NULL; \
  _p_->currentFrame = 0


/*
 *  ACTION: TRANSFORM_HBOND -- this structure needs to be global since it is
 *                             accessed by analyze.c functions
 */

typedef struct _transformHBondInfo {
  ptrajState *state;
  char *name;
  int numdonor;
  int *donor;
  int numacceptor;
  int *acceptor;
  int *acceptorH;
  int solventNeighbor;
  int numSolventDonor;
  int *solventDonor;
  int numSolventAcceptor;
  int *solventAcceptor;
  int *solventAcceptorH;
  int window;
  int series;
  int visit;
  float timeinterval;
  float distanceCutoff;
  float angleCutoff;
  float *occupied;
  float *distance;
  float *distance2;
  float *angle;
  float *angle2;
  float *seriesOccupied;
  float *seriesDistance;
  float *seriesAngle;
  FILE* hbfile;
    /* Integer array to hold how many times a solvent residue has been
     * hydrogen bonded. For example if a solvent molecule has
     * hbonded to the solute more than once (bridging) this is useful
     * to know.
     */
  int* resIsBonded; 
    /* arrays to calculate average # hbonds per residue. Since
     * the total # hbonds a residue makes over the course of a traj
     * may be large, these are long ints.
     */
  long int* avgHbonds;
  long int* avgInterSolute;
  long int* avgSoluteSolvent;
  long int* avgSoluteBridge;
    /* Contains information about what solute molecule(s) the solvent is bound
     * to this frame. E.g.: solvent2solute[x][0] holds the first solute residue #
     * that solvent x is bound to, [x][1] would be the second etc.
     */
  int** solvent2solute;
    /* Store number of frames in long int. NOTE: Is this already stored
     * somewhere else? In ptraj_local.h maxframes is an int.
     */
  long int numframes;
  /*store number of residues for freeing solvent2solute at the end*/
  int numres;

} transformHBondInfo;

#define INITIALIZE_transformHBondInfo(_p_) \
  _p_->state = NULL; \
  _p_->name = NULL; \
  _p_->numdonor = 0; \
  _p_->donor = NULL; \
  _p_->numacceptor = 0; \
  _p_->acceptor = NULL; \
  _p_->acceptorH = NULL; \
  _p_->solventNeighbor = 6; \
  _p_->numSolventDonor = 0; \
  _p_->solventDonor = NULL; \
  _p_->numSolventAcceptor = 0; \
  _p_->solventAcceptor = NULL; \
  _p_->solventAcceptorH = NULL; \
  _p_->window = 0; \
  _p_->series = 0; \
  _p_->visit = 0; \
  _p_->timeinterval = 1.0; \
  _p_->distanceCutoff = 3.0; \
  _p_->angleCutoff = 120.0; \
  _p_->occupied = NULL; \
  _p_->distance = NULL; \
  _p_->distance2 = NULL; \
  _p_->angle = NULL; \
  _p_->angle2 = NULL; \
  _p_->seriesOccupied = NULL; \
  _p_->seriesDistance = NULL; \
  _p_->hbfile = NULL; \
  _p_->resIsBonded = NULL; \
  _p_->avgHbonds = NULL; \
  _p_->avgInterSolute = NULL; \
  _p_->avgSoluteSolvent = NULL; \
  _p_->avgSoluteBridge = NULL; \
  _p_->solvent2solute = NULL; \
  _p_->numframes = 0; \
  _p_->numres = 0

/*
 *  ACTION: TRANSFORM_DIHEDRALCLUSTER
 */
typedef struct _DCnodetype {
  struct _DCnodetype **branch;
  int *bin;
  long int numBranch;
  long int *count;
  double **frames;
} DCnodetype;

#define INITIALIZE_transformDihedralCluster(_p_) \
  _p_->branch = NULL; \
  _p_->bin = NULL; \
  _p_->numBranch = 0; \
  _p_->count = NULL; \
  _p_->frames = NULL ;

typedef struct _DCarray {
  int* Bins;
  long int count;
  double* frames;
} DCarray;




/*
 *  ACTION: TRANSFORM_MATRIX_INFO -- this structure needs to be global since it is
 *                                   accessed by analyze.c functions
 */

typedef enum _transformMatrixType {
  MATRIX_NULL,
  MATRIX_DIST,
  MATRIX_COVAR,
  MATRIX_MWCOVAR,
  MATRIX_CORREL,
  MATRIX_DISTCOVAR,
  MATRIX_IDEA,
  MATRIX_IRED
} transformMatrixType;

typedef struct _transformMatrixInfo {
  ptrajState *state;
  char *name;
  transformMatrixType type;
  double *vect;
  double *vect2;
  double *mat;
  int vectsize;
  int matsize;
  int *mask1;
  int *mask2;
  int mask1tot;
  int mask2tot;
  int snap;
} transformMatrixInfo;

#define INITIALIZE_transformMatrixInfo(_p_) \
  _p_->state         = NULL;   \
  _p_->name          = NULL;   \
  _p_->type          = MATRIX_NULL; \
  _p_->vect          = NULL;   \
  _p_->vect2         = NULL;   \
  _p_->mat           = NULL;   \
  _p_->vectsize      = 0;      \
  _p_->matsize       = 0;      \
  _p_->mask1         = NULL;   \
  _p_->mask2         = NULL;   \
  _p_->mask1tot      = 0;      \
  _p_->mask2tot      = 0;      \
  _p_->snap          = 0;


/*
 *  ACTION: TRANSFORM_PROJECTION_INFO
 */

typedef struct _transformProjectionInfo {
  char *outfile;
  FILE *fp;
  int ibeg;
  int iend;
  int start;
  int stop;
  int offset;
  double *sqrtmasses;
} transformProjectionInfo;

#define INITIALIZE_transformProjectionInfo(_p_) \
  _p_->outfile    = NULL; \
  _p_->fp         = NULL; \
  _p_->ibeg       = 1;    \
  _p_->iend       = 2;    \
  _p_->start      = 1;    \
  _p_->stop       = -1;   \
  _p_->offset     = 1;    \
  _p_->sqrtmasses = NULL;  


/*
 *  ACTION: TRANSFORM_VECTOR
 */

typedef enum _vectorMode { 
  VECTOR_NOOP, 
  VECTOR_PRINCIPAL_X, 
  VECTOR_PRINCIPAL_Y, 
  VECTOR_PRINCIPAL_Z, 
  VECTOR_DIPOLE,
  VECTOR_BOX,
  VECTOR_MASK,
  VECTOR_IRED,
  VECTOR_CORRPLANE,
  VECTOR_CORR,
  VECTOR_CORRIRED
} vectorMode;

typedef struct _transformVectorInfo {
  char *name;
  char *filename;
  int totalFrames;
  int frame;
  vectorMode mode;
  int *mask;
  int *mask2;
  double *cx, *cy, *cz;
  double *vx, *vy, *vz;

  int master;
  modesInfo *modinfo;
  int ibeg;
  int iend;
  int order;
  int npair;
  double *avgcrd;
  double rave;
  double r3iave;
  double r6iave;
  double *cftmp;
  double *p2cftmp;
  double *rcftmp;
} transformVectorInfo;

#define INITIALIZE_transformVectorInfo(_p_) \
  _p_->name        = NULL;                \
  _p_->filename    = NULL;                \
  _p_->totalFrames = 0;                   \
  _p_->frame       = 0;                   \
  _p_->mode        = VECTOR_NOOP;         \
  _p_->mask        = NULL;                \
  _p_->mask2       = NULL;                \
  _p_->cx          = NULL;                \
  _p_->cy          = NULL;                \
  _p_->cz          = NULL;                \
  _p_->vx          = NULL;                \
  _p_->vy          = NULL;                \
  _p_->vz          = NULL;                \
  _p_->master      = 0;                   \
  _p_->modinfo     = NULL;                \
  _p_->ibeg        = 1;                   \
  _p_->iend        = 50;                  \
  _p_->order       = 2;                   \
  _p_->npair       = -1;                  \
  _p_->avgcrd      = NULL;                \
  _p_->rave        = 0.0;                 \
  _p_->r3iave      = 0.0;                 \
  _p_->r6iave      = 0.0;                 \
  _p_->cftmp       = NULL;                \
  _p_->p2cftmp     = NULL;                \
  _p_->rcftmp      = NULL;


/*
 * (2) GLOBAL VARIABLES
 */

#ifdef ACTION_MODULE

int closestWaters = 0;
int *closestWatersMask = NULL;

stackType *hbondDonorStack = NULL;
stackType *hbondAcceptorStack = NULL;
stackType *hbondAcceptorH1Stack = NULL;
stackType *hbondAcceptorH2Stack = NULL;
stackType *hbondAcceptorH3Stack = NULL;
int *hbondDonor = NULL;
int *hbondAcceptor = NULL;
int *hbondAcceptorH1 = NULL;
int *hbondAcceptorH2 = NULL;
int *hbondAcceptorH3 = NULL;

   /*
    * to keep track of the reference configurations for RMS calculations 
    */

#else

extern int closestWaters;
extern int *closestWatersMask;

extern stackType *hbondDonorStack;
extern stackType *hbondAcceptorStack;
extern stackType *hbondAcceptorH1Stack;
extern stackType *hbondAcceptorH2Stack;
extern stackType *hbondAcceptorH3Stack;
extern int *hbondDonor;
extern int *hbondAcceptor;
extern int *hbondAcceptorH1;
extern int *hbondAcceptorH2;
extern int *hbondAcceptorH3;

#endif



/*
 * (3) EXTERNALLY VISIBLE FUNCTION PROTOTYPES
 */

/* PF - multiptraj - new methods to standardize warnings */
void printError(char *, char *, ...);
void printParallelError(char *);

/*
 *  PF - multiptraj
 *  New methods to make it easy for actions to write to
 *  files, regardless of whether using MPI or not
 */
void *ptrajOpenW(char *);
void ptrajfprintf(void *, char *, ...);
void ptrajfprintfone(void *, char *, ...);
void ptrajCloseFile(void *);

#ifndef ACTION_MODULE


#  ifdef __STDC__

extern int distindex(int, int, int);

extern int actionTest(actionInformation *, 
		      double *, double *, double *, double *, int);
extern int transformAngle(actionInformation *, 
			  double *, double *, double *, double *, int);
extern int transformAtomicFluct(actionInformation *, 
				double *, double *, double *, double *, int);
extern int transformAtomicFluct3D(actionInformation *, 
				  double *, double *, double *, double *, int);
extern int transformAverage(actionInformation *, 
			    double *, double *, double *, double *, int);
extern int transformCenter(actionInformation *, 
			   double *, double *, double *, double *, int);
extern int transformCheckOverlap(actionInformation *, 
				 double *, double *, double *, double *, int);
extern int transformClosestWaters(actionInformation *, 
				  double *, double *, double *, double *, int);
extern int transformCluster(actionInformation *, 
				  double *, double *, double *, double *, int);
extern int transformClusterAttribute(actionInformation *, 
				  double *, double *, double *, double *, int);
extern int transformContacts(actionInformation *,
			     double *, double *, double *, double *, int);
extern int transformCorr(actionInformation *, 
			 double *, double *, double *, double *, int);
extern int transformDiffusion(actionInformation *, 
			      double *, double *, double *, double *, int);
extern int transformDihedral(actionInformation *, 
			     double *, double *, double *, double *, int);
extern int transformDihedralCluster(actionInformation *, 
				    double *, double *, double *, double *, int);
extern int transformDipole(actionInformation *, 
			   double *, double *, double *, double *, int);
extern int transformDistance(actionInformation *, 
			     double *, double *, double *, double *, int);
extern int transformDNAiontracker(actionInformation *, 
				  double *, double *, double *, double *, int);
extern int transformEcho(actionInformation *, 
			   double *, double *, double *, double *, int);
extern int transformEnergy(actionInformation *, 
			   double *, double *, double *, double *, int);
extern int transformGrid(actionInformation *, 
			 double *, double *, double *, double *, int);
extern int transformHBond(actionInformation *, 
			  double *, double *, double *, double *, int);
extern int transformImage(actionInformation *, 
			  double *, double *, double *, double *, int);
extern int transformMatrix(actionInformation *, 
			   double *, double *, double *, double *, int);
extern int transformPrincipal(actionInformation *, 
			      double *, double *, double *, double *, int);
extern int transformProjection(actionInformation *, 
			       double *, double *, double *, double *, int);
extern int transformPucker(actionInformation *, 
			   double *, double *, double *, double *, int);
extern int transformRadial(actionInformation *, 
			   double *, double *, double *, double *, int);
extern int transformRadiusOfGyration(actionInformation *,
			     double *, double *, double *, double *, int);
extern int transformRandomizeIons(actionInformation *, 
				  double *, double *, double *, double *, int);
extern int transformRMS(actionInformation *, 
			double *, double *, double *, double *, int);
extern int transformRunningAverage(actionInformation *, 
				   double *, double *, double *, double *, int);
extern int transformScale(actionInformation *, 
			  double *, double *, double *, double *, int);
extern int transformSecondaryStruct(actionInformation *,
			     double *, double *, double *, double *, int);
extern int transformStrip(actionInformation *, 
			  double *, double *, double *, double *, int);
extern int transformTranslate(actionInformation *, 
			      double *, double *, double *, double *, int);
extern int transformTruncOct(actionInformation *, 
			     double *, double *, double *, double *, int);
extern int transformUnwrap(actionInformation *,
                           double *, double *, double *, double *, int);
extern int transformVector(actionInformation *, 
			   double *, double *, double *, double *, int);
extern int transformWatershell(actionInformation *, 
			       double *, double *, double *, double *, int);
extern int transform2dRMS(actionInformation *, 
			  double *, double *, double *, double *, int);
extern actionInformation* ptrajCopyAction(actionInformation**);

#  else  /* __STDC__ */

extern int distindex();
extern int actionTest();
extern int transformAngle();
extern int transformAtomicFluct();
extern int transformAtomicFluct3D();
extern int transformAverage();
extern int transformCenter();
extern int transformCheckOverlap();
extern int transformClosestWaters();
extern int transformCluster();
extern int transformClusterAttribute();
extern int transformContacts();
extern int transformCorr();
extern int transformDihedral();
extern int transformDihedralCluster();
extern int transformDiffusion();
extern int transformDipole();
extern int transformDistance();
extern int transformDNAiontracker();
extern int transformEcho();
extern int transformEnergy();
extern int transformGrid();
extern int transformHBond();
extern int transformImage();
extern int transformMatrix();
extern int transformPrincipal();
extern int transformProjection();
extern int transformPucker();
extern int transformRadial();
extern int transformRadiusOfGyration();
extern int transformRandomizeIons();
extern int transformRMS();
extern int transformRunningAverage();
extern int transformScale();
extern int transformSecondaryStruct();
extern int transformStrip();
extern int transformTranslate();
extern int transformTruncOct();
extern int transformVector();
extern int transformWatershell();
extern int transform2dRMS();
extern actionInformation* ptrajCopyAction();

#  endif /* __STDC__ */


#endif /* ACTION_MODULE */



/*
 * (4) LOCAL STRUCTURES
 */


#ifdef ACTION_MODULE

/*
 *  ACTION: TRANSFORM_CORRELATION
 */

typedef struct _transformCorrInfo {
  char *name;
  char *filename;
  int totalFrames;
  int frame;
  int mode;
  int *mask;
  int *mask2;
  double *cx, *cy, *cz;
  double *vx, *vy, *vz;
  int tmin;
  int tcorr;
  int tmax;
} transformCorrInfo;

enum _corrInfo { CORR_BOX, CORR_MASK };

#define INITIALIZE_transformCorrInfo(_p_) \
  _p_->name        = NULL;                \
  _p_->filename    = NULL;                \
  _p_->totalFrames = 0;                   \
  _p_->frame       = 0;                   \
  _p_->mode        = 0;                   \
  _p_->mask        = NULL;                \
  _p_->mask2       = NULL;                \
  _p_->cx          = NULL;                \
  _p_->cy          = NULL;                \
  _p_->cz          = NULL;                \
  _p_->vx          = NULL;                \
  _p_->vy          = NULL;                \
  _p_->vz          = NULL;                \
  _p_->tmin        = 0;                   \
  _p_->tcorr       = 0;                   \
  _p_->tmax        = 0


/*
 *  ACTION: TRANSFORM_DIFFUSION
 */

typedef struct _transformDiffusionInfo {
  double *dx;           /* reference coordinates: MUST BE SET */
  double *dy;           /*  TO NULL IN transformSetup         */
  double *dz;           
  double timePerFrame;  /*  in picoseconds */
  double *prevx;        /*  the previous x for active atoms */
  double *prevy;        /*  the previous y                  */
  double *prevz;        /*  the previous z                  */
  double *distancex;    /*  the active atoms distances      */
  double *distancey;    /*  the active atoms distances      */
  double *distancez;    /*  the active atoms distances      */
  double *distance;     /*  the active atoms distances      */
  double *deltax;
  double *deltay;
  double *deltaz;
  int activeAtoms;      /*  the total number of active atoms in the mask */
  int elapsedFrames;    /*  number of frames so far... */
  char *outputFilenameRoot;
  FILE *outputx;
  FILE *outputy;
  FILE *outputz;
  FILE *outputr;
  FILE *outputa;
  FILE *outputxyz;
} transformDiffusionInfo;

#define INITIALIZE_transformDiffusionInfo(_p_) \
  _p_->dx             = NULL; \
  _p_->dy             = NULL; \
  _p_->dz             = NULL; \
  _p_->timePerFrame   = 0.0;  \
  _p_->prevx          = NULL; \
  _p_->prevy          = NULL; \
  _p_->prevz          = NULL; \
  _p_->distancex      = NULL; \
  _p_->distancey      = NULL; \
  _p_->distancez      = NULL; \
  _p_->distance       = NULL; \
  _p_->deltax         = NULL; \
  _p_->deltay         = NULL; \
  _p_->deltaz         = NULL; \
  _p_->activeAtoms    = 0; \
  _p_->elapsedFrames  = 0; \
  _p_->outputFilenameRoot = NULL; \
  _p_->outputx        = NULL; \
  _p_->outputy        = NULL; \
  _p_->outputz        = NULL; \
  _p_->outputr        = NULL; \
  _p_->outputxyz      = NULL



/*
 *  ACTION: TRANSFORM_GRID, TRANSFORM_DIPOLE
 */

typedef struct _transformGridInfo {
  double dx;
  double dy;
  double dz;
  int nx;
  int ny;
  int nz;
  int frames;
  float *grid;
  float *dipolex;
  float *dipoley;
  float *dipolez;
  char *filename;
} transformGridInfo;

#define INITIALIZE_transformGridInfo(_p_) \
  _p_->dx        = 0.0;  \
  _p_->dy        = 0.0;  \
  _p_->dz        = 0.0;  \
  _p_->nx        = 0;    \
  _p_->ny        = 0;    \
  _p_->nz        = 0;    \
  _p_->frames    = 0;    \
  _p_->grid      = NULL; \
  _p_->dipolex   = NULL; \
  _p_->dipoley   = NULL; \
  _p_->dipolez   = NULL; \
  _p_->filename  = NULL

/*
 *  ACTION: TRANSFORM_RADIAL_INFO
 */

typedef struct _transformRadialInfo {
  double spacing;
  double maximum;
  double maximum_obs;
  double density;
  double minimum;
  double volume;
  double *distances;
  int *solventMask;
  int solventCount;
  int *soluteMask;
  int soluteCount;
  int *histogram;
  int visits;
  int closest;
  int measurements;
  char *fileroot;
} transformRadialInfo;

#define INITIALIZE_transformRadialInfo(_p_) \
  _p_->spacing       = 0.0;   \
  _p_->maximum       = 0.0;   \
  _p_->maximum_obs   = 0.0;   \
  _p_->density       = 0.0;   \
  _p_->minimum       = 0.0;   \
  _p_->volume        = 0.0;   \
  _p_->distances     = NULL;  \
  _p_->solventMask   = NULL;  \
  _p_->solventCount  = 0;     \
  _p_->soluteMask    = NULL;  \
  _p_->soluteCount   = 0;     \
  _p_->histogram     = NULL;  \
  _p_->visits        = 0;     \
  _p_->closest       = 0;     \
  _p_->measurements  = 0;    \
  _p_->fileroot      = NULL



/*
 *  ACTION: TRANSFORM_WATERSHELL
 */

typedef struct _transformShellInfo {
  int *soluteMask;
  int *solventMask;
  int *activeResidues;
  int *lower;
  int *upper;
  int visits;
  double lowerCutoff;
  double upperCutoff;
  char *filename;

} transformShellInfo;

#define INITIALIZE_transformShellInfo(_p_) \
  _p_->soluteMask     = NULL; \
  _p_->solventMask    = NULL; \
  _p_->activeResidues = NULL; \
  _p_->lower          = NULL; \
  _p_->upper          = NULL; \
  _p_->visits         = 0;    \
  _p_->lowerCutoff    = 0.0;  \
  _p_->upperCutoff    = 0.0;  \
  _p_->filename       = NULL

#endif /* ACTION_MODULE */

/*
 *  ACTION: TRANSFORM_UNWRAP
 */

typedef struct _transformUnwrapInfo {
  double *refx;
  double *refy;
  double *refz;
} transformUnwrapInfo;
