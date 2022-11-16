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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/ptraj_local.h,v 10.3 2009/09/25 15:08:56 case Exp $
 *
 *  Revision: $Revision: 10.3 $
 *  Date: $Date: 2009/09/25 15:08:56 $
 *  Last checked in by $Author: case $
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




/* The parameter file assumes the atom, residue, symbol, etc. names to be *
 * four characters (we will store them as strings, requiring a newline,   *
 * hence the size is 5).                                                  */
#define NAME_SIZE 5
#define NAME_DEFAULT "    "
typedef char Name[NAME_SIZE];



/*
 *  state information necessary for the transformation routines
 */
typedef struct _ptrajState { 
  double box[6];             /* box lengths and angles */
  double *masses;            /* atom masses */
  double *charges;           /* atom charges */
  int atoms;                 /* number of atoms */
  int residues;              /* number of residues */
  int *ipres;                /* atoms in each residue */
  int IFBOX;                 /* is there box information? */
  int boxfixed[6];           /* equal to 1 if this box coordinate is fixed */
  int molecules;             /* total number of molecules */
  int *moleculeInfo;         /* number of atoms in each molecule */
  int *solventMask;          /* atoms in the solvent */
  int solventMolecules;      /* number of solvent molecules */
  int *solventMoleculeStart; /* pointer into solventMask for first atom of each solvent */
  int *solventMoleculeStop;  /* pointer into solventMask for last atom of each solvent */
  int solventAtoms;          /* number of solvent atoms */
  Name *atomName;            /* atom names */
  Name *residueName;         /* residue names */
  int maxFrames;             /* number of useful frames in 1 pass of trajin's */
  double temp0;              /* DAN TEST: for writing out netcdf temp trajs */
} ptrajState;

#define INITIALIZE_ptrajState( _p_ ) \
  _p_->box[0] = 0.0; _p_->box[1] = 0.0; _p_->box[2] = 0.0; \
  _p_->box[3] = 90.0; _p_->box[4] = 90.0; _p_->box[5] = 90.0; \
  _p_->masses = NULL; \
  _p_->charges = NULL; \
  _p_->atoms = 0; \
  _p_->residues = 0; \
  _p_->ipres = NULL; \
  _p_->IFBOX = 0; \
  _p_->boxfixed[0] = 0; _p_->boxfixed[1] = 0; _p_->boxfixed[2] = 0; \
  _p_->boxfixed[3] = 0; _p_->boxfixed[4] = 0; _p_->boxfixed[5] = 0; \
  _p_->molecules = 0; \
  _p_->moleculeInfo = NULL; \
  _p_->solventMask = NULL; \
  _p_->solventMolecules = 0; \
  _p_->solventMoleculeStart = NULL; \
  _p_->solventMoleculeStop = NULL; \
  _p_->solventAtoms = 0; \
  _p_->atomName = NULL; \
  _p_->residueName = NULL; \
  _p_->maxFrames = 0; \
  _p_->temp0 = 0.0

  

/*
 *  structures used in the processing of the atom specification
 */

enum _parseOperators { 
  PARSE_NOOP,
  PARSE_CONTINUATION,
  PARSE_ATOM,
  PARSE_RESIDUE
};

typedef struct _parseEntry {
  char *token;
  int operator;
  int isnumber;
} parseEntry;


typedef enum _ptrajMode {
  PTRAJ_NOOP,
  PTRAJ_ACTION,
  PTRAJ_FIRSTPASS,
  PTRAJ_SECONDPASS,
  PTRAJ_SETUP, 
  PTRAJ_STATUS,
  PTRAJ_PRINT,
  PTRAJ_CLEANUP
} ptrajMode;

