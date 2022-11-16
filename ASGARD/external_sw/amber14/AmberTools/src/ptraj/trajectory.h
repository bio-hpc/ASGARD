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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/trajectory.h,v 10.6 2010/01/27 22:03:13 droe Exp $
 *
 *  Revision: $Revision: 10.6 $
 *  Date: $Date: 2010/01/27 22:03:13 $
 *  Last checked in by $Author: droe $
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

#ifdef MPI
#include <mpi.h>
#endif

/*  ________________________________________________________________________
 */



   /*
    *  possible coordinate types
    */
typedef enum _coordType { 
  COORD_UNKNOWN,
  COORD_AMBER_TRAJECTORY,
  COORD_AMBER_RESTART,
  COORD_AMBER_NETCDF,
  COORD_AMBER_REMD,
  COORD_PDB,
  COORD_BINPOS,
  COORD_CHARMM_TRAJECTORY
} coordType;

typedef enum _LesAction {
    LES_NONE,
    LES_SPLIT,
    LES_COMBINE,
    LES_AVERAGE
} LesAction;

typedef enum _LesStatus {
    LES_EMPTY,
    LES_DONE,
    LES_READY
} LesStatus;

// ---------- Information for a NETCDF file ----------
typedef struct _netcdfTrajectoryInfo {
  int ncid;
  int currentFrame;
  int frameDID;
  int spatialDID;
  int atomDID;
  int cell_spatialDID;
  int cell_angularDID;
  int labelDID;
  int spatialVID;
  int timeVID;
  int coordinateVID;
  int cell_spatialVID;
  int cell_angularVID;
  int cellLengthVID;
  int cellAngleVID;
  int velocityVID;
  int TempVarID;
  double velocityScale;
  char *timeUnits;
  char *coordinateUnits;
  char *cellLengthUnits;
  char *cellAngleUnits;
  char *velocityUnits;
  char *Conventions;
  char *ConventionVersion;
  float *R;

} netcdfTrajectoryInfo;

#define INITIALIZE_netcdfTrajectoryInfo( _p_ ) \
  _p_->ncid = -1; \
  _p_->currentFrame = 0; \
  _p_->frameDID = 0; \
  _p_->spatialDID = 0; \
  _p_->atomDID = 0; \
  _p_->cell_spatialDID = 0; \
  _p_->cell_angularDID = 0; \
  _p_->labelDID = 0; \
  _p_->spatialVID = 0; \
  _p_->timeVID = 0; \
  _p_->coordinateVID = 0; \
  _p_->cell_spatialVID = 0; \
  _p_->cell_angularVID = 0; \
  _p_->cellLengthVID = 0; \
  _p_->cellAngleVID = 0; \
  _p_->velocityVID = 0; \
  _p_->TempVarID = 0; \
  _p_->timeUnits = NULL; \
  _p_->coordinateUnits = NULL; \
  _p_->cellLengthUnits = NULL; \
  _p_->cellAngleUnits = NULL; \
  _p_->velocityUnits = NULL; \
  _p_->velocityScale = 0.0; \
  _p_->Conventions = NULL; \
  _p_->ConventionVersion = NULL; \
  _p_->R = NULL

#define FREE_netcdfTrajectoryInfo( _p_ ) \
  safe_free(_p_->timeUnits); \
  safe_free(_p_->coordinateUnits); \
  safe_free(_p_->cellLengthUnits); \
  safe_free(_p_->cellAngleUnits); \
  safe_free(_p_->velocityUnits); \
  safe_free(_p_->Conventions); \
  safe_free(_p_->ConventionVersion); \
  safe_free(_p_->R); \
  safe_free(_p_)

#define AMBER_NETCDF_FRAME "frame"
#define AMBER_NETCDF_SPATIAL "spatial"
#define AMBER_NETCDF_ATOM "atom"
#define AMBER_NETCDF_CELL_SPATIAL "cell_spatial"
#define AMBER_NETCDF_CELL_ANGULAR "cell_angular"
#define AMBER_NETCDF_COORDS "coordinates"
#define AMBER_NETCDF_TIME "time"
#define AMBER_NETCDF_LABEL "label"
#define AMBER_NETCDF_LABELLEN 5

    
// ---------- Information for trajectory  files ----------
/*
 *  information relating to a coordinate file
 */
typedef struct _coordinateInfo {
  FILE *file;      // File pointer
#ifdef MPI
  MPI_File *mfp;
#endif
  char *filename;  // File name
  int start;       // Frame to start processing
  int stop;        // Frame to end processing
  int Nframes;     // Total number of frames in the file.
  int offset;      // # of frames to skip
  int append;      // File will be appended to
  int isBox;       // File has box information
  int isVelocity;  
  int option1;
  int option2;
  int *mask;
  double *x;
  double *y;
  double *z;
  double *time;
  double *vx;
  double *vy;
  double *vz;
  char *title;
  char *application;
  char *program;
  char *version;
  void *info;          // Holds NETCDF or CHARMM trajectory info
  netcdfTrajectoryInfo *NCInfo;  // Holds NETCDF trajectory info
  coordType type;      // Identify coordinate type
  int accessMode;      // 0 for read, 1 for write, 2 for append
  int compressType;    // 1 gzip 2 bzip 3 zip 
    /*
     *  LES information
     */
  LesAction les_action;
  int nlescopy;
  LesStatus les_status;
    /*
     *  REMD Trajectory info
     */
  struct _coordinateInfo **REMDtraj; // Hold information of other replica trajectories
  int isREMDTRAJ;       /* 0 - Process as normal trajectory,
                         * 1 - Process as replica trajectory
                         */
  char *compressEXT;    /* If not NULL, contains compressed traj ext.*/
  int numREMDTRAJ;      /* How many replica trajectories are present */
  int firstREMDTRAJ;    /* Index of first replica                    */
  int EXTwidth;         /* Length of replica extension               */
  char* baseFilename;   /* To hold replica traj base filename        */
  int linesperset;      /* for fast scanthrough of other REMD files  */
  double remdtrajtemp;  /* Target temperature                        */
  /* Write file with MPI IO? */
  int isMPI;
  int isNetcdf;
  /*
   *  AMBER trajectory file information
   */
  int seekable;
  int titleSize;
  int frameSize;
  int numBox;
  char *buffer;
    
} coordinateInfo;

#define INITIALIZE_coordinateInfo(_p_) \
  _p_->file = NULL; \
  _p_->filename = NULL; \
  _p_->start = 1; \
  _p_->stop = -1; \
  _p_->Nframes = 0; \
  _p_->offset = 1; \
  _p_->append = 0; \
  _p_->isBox = 0; \
  _p_->isVelocity = 0; \
  _p_->option1 = 0; \
  _p_->option2 = 0; \
  _p_->mask = NULL; \
  _p_->frameSize = 0; \
  _p_->x = NULL; \
  _p_->y = NULL; \
  _p_->z = NULL; \
  _p_->time = NULL; \
  _p_->vx = NULL; \
  _p_->vy = NULL; \
  _p_->vz = NULL; \
  _p_->title = NULL; \
  _p_->application = NULL; \
  _p_->program = NULL; \
  _p_->version = NULL; \
  _p_->info = NULL; \
  _p_->NCInfo = NULL; \
  _p_->type = COORD_UNKNOWN; \
  _p_->accessMode = 0; \
  _p_->compressType = 0; \
  _p_->REMDtraj = NULL; \
  _p_->isREMDTRAJ = 0; \
  _p_->compressEXT = NULL; \
  _p_->numREMDTRAJ = 0; \
  _p_->firstREMDTRAJ = 0; \
  _p_->EXTwidth = 0; \
  _p_->baseFilename = NULL; \
  _p_->linesperset = 0; \
  _p_->remdtrajtemp = 0.0; \
  _p_->isMPI = 0; \
  _p_->isNetcdf = 0; \
  _p_->seekable = 0; \
  _p_->titleSize = 0; \
  _p_->frameSize = 0; \
  _p_->numBox = 0; \
  _p_->buffer = NULL; \

   /*
    *  union structures for binary I/O, i.e. read in as bytes which can be 
    *  interpreted as float or integer
    */
typedef union _byte {
  char   c[4];
  int    i;
  float  f;
} byte;


typedef union _bytebyte {
  char   c[8];
  double d;
} bytebyte;


typedef struct _charmmTrajectoryInfo {
  int byteorder;         /* endian-ness, 1:LITTLE or [3][2][1][0]; 0:BIG or [0][1][2][3] */
  byte magic;            /* the magic header, should be "CORD" */
  int icntrl[20];        /* the control variables */
  int ntitle;            /* the number of title lines */
  stackType *titleStack; /* a stack for the title */
  int natrec;            /* the number of written atom records */
  int nfreat;            /* the number of free atoms */
  int *freeat;           /* a list of the free atoms */
  double xtlabc[6];      /* the crystal information in XTLABC form */
  double *fixedx;        /* fixed atom coords */
  double *fixedy;        /* fixed atom coords */
  double *fixedz;        /* fixed atom coords */
} charmmTrajectoryInfo;

#define INITIALIZE_charmmTrajectoryInfo( _p_ ) \
  _p_->byteorder = 0; \
  _p_->magic.i = 0; \
  _p_->titleStack = NULL; \
  _p_->natrec = 0; \
  _p_->nfreat = 0; \
  _p_->freeat = NULL; \
  _p_->fixedx = NULL; \
  _p_->fixedy = NULL; \
  _p_->fixedz = NULL



#ifndef TRAJECTORY_MODULE
#  ifdef __STDC__

extern Restart * readAmberRestart(int, char *, FILE *);
extern int getCoordinatesFromRestart(FILE *, double *, double *, double *,
                                    double *, double *, double *,
				     double *, double *, double *);
extern int getCoordinatesFromPdb(pdb_record *, double *, double *, double *);
extern void dumpAmberRestart(FILE *, int, double *, double *, double *, 
			     double *, double *, double *, double *, char *);
extern void writeAmberRestart(Restart *, char *, FILE *);
extern int loadPdb(FILE *, pdb_record **);
extern void savePdb(FILE *, pdb_record *);
extern pdb_record *parmToPdb(Parm *, int *);
extern pdb_record *ptrajStateToPdb(ptrajState *, int *, int);
extern void putCoordinatesInPdb(pdb_record *, int, 
				double *, double *, double *);
extern int readAmberTrajectory(FILE *, int, double *, double *, double *, 
			       double *, int, coordinateInfo *);
extern int readAmberTrajectory_nobuffer(FILE *, int, double *, double *, double *, 
                                        double *, int, int);

extern void dumpAmberTrajectory(coordinateInfo *, int, double *, double *, double *, double *);

extern void dumpAmberBox(FILE *, double *);
extern void writebinpos(FILE *, int, double *, double *, double *);
extern int openbinpos();
extern int readbinpos(FILE *, int *, float *, int *);
extern int readCharmmTrajectory(FILE *, charmmTrajectoryInfo **, 
				double *, double *, double *, double *, int);
extern void dumpCharmmTrajectory(FILE *, charmmTrajectoryInfo *, int, 
				 double *, double *, double *, double *, int);
extern void printString(void *);

extern int openTraj(coordinateInfo *);
extern int closeTraj(coordinateInfo *);
extern void cleanTraj(coordinateInfo *);
extern char *trajFile_fgets(char *, int, coordinateInfo *);
extern int trajFile_fseek(coordinateInfo *, int);
extern int trajFile_rewind(coordinateInfo *);
extern int trajFile_fseek_end(coordinateInfo *); 
extern int trajFile_get_position(coordinateInfo *, long int *);
extern int trajFile_fread(coordinateInfo *);
extern int bufferToXYZ(char *, int, double *, double *, double *, int, double *);

#  else

extern Restart * readAmberRestart();
extern int getCoordinatesFromRestart();
extern int getCoordinatesFromPdb();
extern void dumpAmberRestart();
extern void writeAmberRestart();
extern int loadPdb();
extern void savePdb();
extern pdb_record *parmToPdb();
extern pdb_record *ptrajStateToPdb();
extern putCoordinatesInPdb();
extern int readAmberTrajectory();
extern int readAmberTrajectory_nobuffer();
extern void dumpAmberTrajectory();
extern void dumpAmberBox();
extern void writebinpos();
extern int openbinpos();
extern int readbinpos();
extern int readCharmmTrajectory();
extern void dumpCharmmTrajectory();
extern void printString();

extern int openTraj();
extern int closeTraj(); 
extern void cleanTraj();
extern char *trajFile_fgets();
extern int trajFile_fseek();
extern int trajFile_rewind();
extern int trajFile_fseek_end();
extern int trajFile_get_position();
extern int trajFile_fread();
extern int bufferToXYZ();

#  endif
#endif





