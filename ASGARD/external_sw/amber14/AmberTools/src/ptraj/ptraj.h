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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/ptraj.h,v 10.3 2010/01/27 22:05:31 droe Exp $
 *
 *  Revision: $Revision: 10.3 $
 *  Date: $Date: 2010/01/27 22:05:31 $
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

/*  ________________________________________________________________________
 */


/*
 *        This is the main header file for the ptraj and rdparm code.         
 *
 *   NOTE: All of the prototypes and structure definitions reside in 
 *   individual files based on the name of the source code file.     
 *   Despite this logical separation, all the individual include     
 *   files are included indirectly and automatically by including    
 *   this single include file.  To separate local and global definitions
 *   within the header files, environment variables are set in each of 
 *   the header files (#define X_MODULE where X is the name of the source 
 *   file) to differentiate local information needed for the source file
 *   and globabally accessible (and prototyped) definitions.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <memory.h>
#include <time.h>
                                                                   
#include "constants.h"
#include "ptraj_local.h"
#include "utility.h"
#include "pdb/pdb.h"
#include "evec.h"
#include "dispatch.h"
#include "actions.h"
#include "analyze.h"
#include "rdparm.h"
#include "io.h"
#include "help.h"
#include "vector.h"
#include "trajectory.h"
#include "netcdf_ptraj.h"
#include "parallel_ptraj.h"
#include "torsion.h"
#include "rms.h"
#include "display.h"
#include "interface.h"
#include "energy.h"
#include "mask.h"
#include "experimental.h"
#include "clusterLib.h"
#include "cluster.h"

#ifdef BINTRAJ
#ifdef MPI
#define nc_strerror ncmpi_strerror
#include "../../include/pnetcdf.h"
#else
#include "netcdf.h"
#endif
#endif

/*
 *   MACHINE DEPENDENCY
 *
 *   Although most of this code is written in C, FORTRAN routines may
 *   also be called although none currently is anymore...  
 *   Various systems handle calling FORTRAN routines from C differently.
 *   We can overcome this the use of the C preprocess (CPP).  Most systems 
 *   call FORTRAN function "foobar" from C via the "foobar_()" call.
 *   Exceptions:
 *     HPUX       -- no appended underscore (i.e. "foobar()")
 *     UNICOS     -- no appended underscore, all capitals! (i.e. "FOOBAR()")
 */


/*
 * CRAY special case to avoid expf call...
 */
#ifdef CRAY
#define expf exp
#endif


/*
 *  GLOBAL DEFINES and ENUMERATED TYPES local to ptraj.c
 */

#ifdef PTRAJ_MODULE

ptrajState *ptrajCurrentStatePointer = NULL;  
  /*
   *  this holds the current state which may evolve as actions are setup.
   *  Each action has a local copy (setup as the action is initialized...
   */

stackType *transformFileStack = NULL;
coordinateInfo *globalOutInfo = NULL;
stackType *transformReferenceStack = NULL;
coordinateInfo *referenceInfo = NULL;
stackType *transformActionStack = NULL;
stackType *transformAnalyzeStack = NULL;
stackType *vectorStack = NULL;
stackType *scalarStack = NULL;
stackType *hbondStack = NULL;
stackType *matrixStack = NULL;
stackType *modesStack = NULL;

stackType *attributeStack = NULL;
arrayType *attributeArray = NULL;
arrayType *attributeArrayTorsion = NULL;

int prnlev = 0;

/* MULTIPTRAJ */
int worldsize;   /* Total # of processors in run */
int worldrank;   /* Rank of current processor    */
/* MULTIPTRAJ */

int bench, benchshort;
char *benchfile;
double checkInputTime;
double inputTime;
double outputTime;
double actionTime;
double t1, t2;

#else

extern stackType *vectorStack;
extern stackType *scalarStack;
extern stackType *hbondStack;
extern stackType *matrixStack;
extern stackType *modesStack;
extern coordinateInfo *globalOutInfo;
extern stackType *transformReferenceStack;
extern coordinateInfo *referenceInfo;

extern stackType *attributeStack;
extern arrayType *attributeArray;
extern arrayType *attributeArrayTorsion;

extern int prnlev;

/* MULTIPTRAJ */
extern int worldsize;   /* Total # of processors in run */
extern int worldrank;   /* Rank of current processor    */
/* MULTIPTRAJ */

extern double checkInputTime;
extern double inputTime;
extern double outputTime;
extern double actionTime;
extern double t1, t2;

#endif /* ifdef PTRAJ_MODULE */

/*
 *   FUNCTION PROTOTYPES
 *
 *   The following functions are defined in ptraj.c and are externally
 *   accessible.  In order to hide this prototyping/referencing from the 
 *   ptraj.c source (which like the other source files includes 
 *   this header file) the C preprocessor (CPP) token "PTRAJ_MODULE" is 
 *   defined at the top of the ptraj.c source to prevent inclusion...  
 *   This general strategy has been used throughout the source
 *   where the defined CPP token is equivalent to the capitalized 
 *   file name without the ".c" appended with "_MODULE"
 */

#ifndef PTRAJ_MODULE
#  ifdef __STDC__

extern void ptraj(char *);
extern ptrajState** ptrajCurrentState();
extern ptrajState* ptrajCopyState(ptrajState **);
extern void checkCoordinatesWrap(char *);
extern coordinateInfo *checkCoordinates(char *, int);
extern int *processAtomMask(char *, ptrajState *);
extern void checkAtomMask(char *);
extern int *returnAtomMask(char *);
extern void atomMaskIsActive(int *, ptrajState *, int *, int *);
extern ptrajState **ptrajCurrentState();
extern void ptrajSetup(stackType *, int);
extern void ptrajSetupAnalyze(stackType *, int);
extern void ptrajSetupIO(stackType *, int);
extern int  checkDivisibility(int *, int, int);
extern void ptrajInitializeState(ptrajState **, char *);
extern int  atomToResidue(int, int, int *);
extern int  atomToSolventMolecule(int, int, int *, int *);
extern int  atomToMolecule(int, int, int *);
extern void modifyStateByMask(ptrajState **, ptrajState **, int *, int);
extern scalarInfo *scalarStackGetName(stackType **, char *);
extern transformHBondInfo *hbondInfoStackGetName(stackType **, char *);
extern transformMatrixInfo *matrixInfoStackGetName(stackType **, char *);
extern modesInfo *modesInfoStackGetName(stackType **, char *);
extern void boxToRecip(double box[6], double ucell[9], double recip[9]);
extern double calculateDistance2(int, int, double *, double *, double *,
				 double *, double *, double *, double, int);
extern double calculateMinImagedDistance2(double *, double *, double *, 
					  double, double, double, double, double, double,
					  int *, int *, int *, int);

extern void ptrajOutputCoordinates(coordinateInfo *, ptrajState *, int, int, int, int, int,
				   double *, double *, double *, double *);
#  else

extern void ptraj();
extern ptrajState** ptrajCurrentState();
extern ptrajState* ptrajCopyState();
extern void checkCoordinatesWrap();
extern coordinateInfo *checkCoordinates();
extern int *processAtomMask();
extern void checkAtomMask();
extern int *returnAtomMask();
extern void atomMaskIsActive();
extern ptrajState **ptrajCurrentState();
extern void ptrajSetup();
extern void ptrajSetupAnalyze();
extern void ptrajSetupIO();
extern void ptrajInitializeState();
extern int  atomToResidue();
extern int  atomToSolventMolecule();
extern int  atomToMolecule();
extern void modifyStateByMask();
extern scalarInfo *scalarStackGetName();
extern transformHBondInfo *hbondInfoStackGetName();
extern transformMatrixInfo *matrixInfoStackGetName();
extern modesInfo *modesInfoStackGetName();
extern void boxToRecip();
extern double calculateDistance2();
extern double calculateMinImagedDistance2();
extern void ptrajOutputCoordinates();

#  endif

#endif /* PTRAJ_MODULE */

void printBench(FILE *, int);

