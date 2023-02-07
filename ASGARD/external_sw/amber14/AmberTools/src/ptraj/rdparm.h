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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/rdparm.h,v 10.0 2008/04/15 23:24:11 case Exp $
 *
 *  Revision: $Revision: 10.0 $
 *  Date: $Date: 2008/04/15 23:24:11 $
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




/* Titles in parameter files are 20A4 or 80 characters; we add an extra   *
 * character to store newline                                             */
#define TITLE_LENGTH 81


/* Defining the atoms.                                                    */
typedef struct _Atom {
  double chrg;                       /* the charge on the atom            */
  double amass;                      /* the mass of the atom              */
  double radii;                      /* continuum radii                   */
  double screen;                     /* GB screening factor               */
  int join;                          /* the tree joining info             */
  int irotat;                        /* last at to move if cur at moved   */
  int iac;                           /* atom types involved in L-J        */
  int numex;                         /* index into excluded atom list     */
  int res;                           /* the residue number                */
  Name igraph;                       /* the true atom name                */
  Name isymbl;                       /* the atom type                     */
  Name itree;                        /* the atom tree symbol              */
} Atom;

#define FREE_Atom(_at_) \
  safe_free(_at_); _at_ = NULL


/* Defining the residues.                                                 */
typedef struct _Residue {
  Name labres;                        /* the residue name                 */
  int ipres;                          /* the pntr list or all residues    */
} Residue;

#define FREE_Residue(_res_) \
  safe_free(_res_); _res_ = NULL


/* Defining the bonds, angles, and dihedrals.  Information is duplicated  *
 * in two forms.  The parm form (to allow re-creation of parm file) and   *
 * the more obvious form, where each ``bond'' has associated parameters   *
 * and optional scale factors.  (This is used in LES spasms version)      *
 */
typedef struct _Bond {
  int atom[2];                        /* the atoms in the bond (2)        */
  int active;                         /* whether this bond is active      */
  double rk;                          /* the force constant               */
  double req;                         /* the minimum bond length          */
  double scale;                       /* scale factor                     */
} Bond;

#define FREE_Bond(_bond_) \
  safe_free(_bond_); _bond_ = NULL


typedef struct _ParmBond {
  int ib;                             /* first atom in bond               */
  int jb;                             /* second atom in bond              */
  int icb;                            /* pointer to parameters            */
} ParmBond;

#define FREE_ParmBond(_pbond_) \
  safe_free(_pbond_); _pbond_ = NULL



typedef struct _Angle {
  int atom[3];                        /* the atoms in the angle (3)       */
  double tk;                          /* the force constant               */
  double teq;                         /* the minimum angle                */
  double scale;                       /* scale factor                     */
} Angle;

#define FREE_Angle(_ang_) \
  safe_free(_ang_); _ang_ = NULL


typedef struct _ParmAngle {
  int it;                             /* first atom in angle              */
  int jt;                             /* second atom in angle             */
  int kt;                             /* third atom in angle              */
  int ict;                            /* pointer to parameters            */
} ParmAngle;

#define FREE_ParmAngle(_pang_) \
  safe_free(_pang_); _pang_ = NULL



typedef struct _Dihedral {
  int atom[4];                        /* the atoms in the dihedral (4)    */
  double pk;                          /* peak height                      */
  double pn;                          /* the periodicity                  */
  double phase;                       /* the phase                        */
  double scale;                       /* scale factor                     */
} Dihedral;

#define FREE_dihedral(_dih_) \
  safe_free(_dih_); _dih_ = NULL

typedef struct _ParmDihedral {
  int ip;                             /* first atom in angle              */
  int jp;			      /* second atom in angle             */
  int kp;			      /* third atom in angle              */
  int lp;			      /* fourth atom in angle             */
  int icp;                            /* pointer to parameters            */
} ParmDihedral;

#define FREE_ParmDihedral(_pdih_) \
  safe_free(_pdih_); _pdih_ = NULL


typedef struct _Box {
  int iptres;                         /* 1st res that is part of solvent  */
  int nspm;                           /* final res that is part of solv.  */
  int nspsol;                         /* atom in solv used for per imag   */
  int *nsp;                           /* num atoms in i'th solv. molecule */
  double beta;                        /* the box angle                    */
  double box[3];                      /* box length in x,y,z              */
} Box;

#define FREE_Box(_box_)     \
  if (_box_ != NULL)        \
     safe_free(_box_->nsp); \
  safe_free(_box_);         \
  _box_ = NULL


typedef struct _Cap {
  int natcap;                         /* number of atoms in cap           */
  double cutcap;                      /* distance to edge of cap          */
  double xcap;                        /* x coord of center                */
  double ycap;                        /* y coord of center                */
  double zcap;                        /* z coord of center                */
} Cap;

#define FREE_Cap(_cap_) \
  safe_free(_cap_); _cap_ = NULL

typedef struct _Pert {
  int *ibper, *jbper;                 /* the bonds to be perturbed        */
  int *icbper;                        /* pntr to bond parameter arrays    */
  int *itper, *jtper, *ktper;         /* angles to be perturbed           */
  int *ictper;                        /* pntr to angle parameter arrays   */
  int *ipper, *jpper, *kpper, *lpper; /* dihedrals to be perturbed        */
  int *icpper;                        /* pntr to dihedral parameter array */
  Name *labper;                       /* residue names at lambda = 1      */
  Name *igrper;                       /* atom names at lambda = 1         */
  Name *ismper;                       /* atomic symbols at lambda = 1     */
  double *almper;                     /* value of lambda for each atom    */
  int *iaper;                         /* = 1 if atom is being perturbed   */
  int *iacper;                        /* atom types at lambda = 1         */
  double *cgper;                      /* atom charges at lambda = 1       */
} Pert;

#define FREE_Pert(_pert_)       \
  if ( _pert_ != NULL ) {       \
    safe_free(_pert_->ibper);   \
    safe_free(_pert_->jbper);   \
    safe_free(_pert_->icbper);  \
    safe_free(_pert_->itper);   \
    safe_free(_pert_->jtper);   \
    safe_free(_pert_->ktper);   \
    safe_free(_pert_->ictper);  \
    safe_free(_pert_->ipper);   \
    safe_free(_pert_->jpper);   \
    safe_free(_pert_->kpper);   \
    safe_free(_pert_->icpper);  \
    safe_free(_pert_->labper);  \
    safe_free(_pert_->igrper);  \
    safe_free(_pert_->ismper);  \
    safe_free(_pert_->almper);  \
    safe_free(_pert_->iaper);   \
    safe_free(_pert_->iacper);  \
    safe_free(_pert_->cgper);   \
  } safe_free(_pert_); _pert_ = NULL


typedef struct _Restart {
  char *title;               /* the restart file title line               */
  int natoms;                /* the number of atoms in the restart file   */
  double time;               /* the current time in the simulation        */
  double *x, *y, *z;         /* the coordinates                           */
  int restart;               /* flag for the presence of velocities       */
  double *vx, *vy, *vz;      /* the velocities, only in restart files     */
  int isbox;                 /* flag for the presence of box values       */
  double box[6];             /* the box lengths and angles                */
} Restart;

#define FREE_Restart(_rst_)    \
  if (_rst_ != NULL ) {        \
    safe_free(_rst_->title);   \
    safe_free(_rst_->x);       \
    safe_free(_rst_->y);       \
    safe_free(_rst_->z);       \
    safe_free(_rst_->vx);      \
    safe_free(_rst_->vy);      \
    safe_free(_rst_->vz);      \
  } safe_free(_rst_); _rst_ = NULL

  
typedef struct _Parm {
  char *filename;            /* the filename represented                  */
  FILE *fp;                  /* FILE descriptor                           */
  char title[TITLE_LENGTH];  /* the title                                 */
                 /* The integer control variables                         */
  int  NTOTAT;   /* total number of atoms in the system                   */
  int  NTYPES;   /* number of AMBER atom types used, max is 60            */
  int  NBONH;    /* number of bonds containing hydrogen                   */
  int  NBONA;    /* number of bonds without hydrogen                      */
  int  NTHETH;   /* number of angles containing hydrogen                  */
  int  NTHETA;   /* number of angles not containing hydrogen              */
  int  NPHIH;    /* number of dihedrals containing hydrogen               */
  int  NPHIA;    /* number of dihedrals not containing hydrogen           */
  int  JHPARM;   /* NOT USED                                              */
  int  JPARM;    /* NOT USED                                              */
  int  NEXT;     /* total number of excluded atoms                        */
  int  NTOTRS;   /* total number of residues                              */
  int  MBONA;    /* NBONA + number of constraint bonds                    */
  int  MTHETS;   /* NTHETS (sic) + number of constraint angles            */
  int  MPHIA;    /* NPHIA + number of constraint dihedral angles          */
  int  MUMBND;   /* total number of unique bond types                     */
  int  MUMANG;   /* total number of unique angle types                    */
  int  MPTRA;    /* total number of unique dihedral types                 */
  int  NATYP;    /* number of "atoms" defined in parameter file           */
  int  NHB;      /* number of types of hydrogen bonded pair interactions  */
  int  IFPERT;   /* =1 if perturbation info is to be read =0 otherwise    */
  int  NBPER;    /* number of bonds to be perturbed                       */
  int  NGPER;    /* number of angles to be perturbed                      */
  int  NDPER;    /* number of dihedrals to be perturbed                   */
  int  MBPER;    /* num of pert bonds across boundary to non-pert groups  */
  int  MGPER;    /* num of pert angles across boundary to non-pert groups */
  int  MDPER;    /* num of pert dihedrals across bndry to non-pert groups */
  int  IFBOX;    /* =1 if periodic box info to be read =0 otherwise       */
  int  NMXRS;    /* number of atoms in the largest residue                */
  int  IFCAP;    /* =1 if CAP option was used in edit, =0 otherwise       */
  int  NUMEXTRA; /* number of extra points (aka lone pairs)               */
  int  RADII;    /* =1 if radii present                                   */
  int  SCREEN;   /* =1 if screen present                                  */
  char RADIUS_SET[81];                /* the radius set                   */
  Atom *atom;                         /* the atoms and related info       */
  int *nno;                           /* index for non-bond of each type  */
  Residue *residue;                   /* the residues and related info    */
  double *rk, *req;                   /* the bond parameters              */
  double *tk, *teq;                   /* the angle parameters             */
  double *pk, *pn, *phase;            /* the dihedral parameters          */
  double *solty;                      /* NOT USED                         */
  double *cn1;                        /* L-J r**12 and r**6 for all pos...*/
  double *cn2;                        /* ...atom type interactions        */  
  ParmBond *pbondH;                   /* bonds with hydrogen              */
  ParmBond *pbond;                    /* bonds without hydrogen           */
  ParmAngle *pangleH;                 /* angles with hydrogen             */
  ParmAngle *pangle;                  /* angles without hydrogen          */
  ParmDihedral *pdihedralH;           /* dihedrals with hydrogen          */
  ParmDihedral *pdihedral;            /* dihedrals without hydrogen       */
  int *natex;                         /* excluded atom list               */
  double *ag;                         /* H-bond r**12 and r**10...        */
  double *bg;                         /* ...                              */
  double *hbcut;                      /* NO LONGER USED                   */
  Box *box;                           /* box info IFBOX > 1               */
  Cap *cap;                           /* cap info IFCAP > 1               */
  Pert *pert;                         /* pert info IFPERT > 1             */
  Bond *bond;                         /* the bonds...                     */
  Angle *angle;                       /* the angles...                    */
  Dihedral *dihedral;                 /* the dihedrals...                 */
  Restart *coords;                    /* optional AMBER coordinates       */
  int trajSize;                       /* actual # of coords in trajectory */
  int trajMaxSize;                    /* max # of coords in trajectory    */
  int nlestyp;                        /* number of LES types              */
  int *lestyp;                        /* the specific type                */
  double *lesfac;                     /* the LES factors                  */
  int *lescnum;                          
  int *lessubsp;
} Parm;

#define FREE_Parm(_prm_)                   \
  if ( _prm_ != NULL ) {                   \
    safe_free(_prm_->fp);                  \
    safe_free(_prm_->filename);            \
    FREE_Atom(_prm_->atom);                \
    safe_free(_prm_->nno);                 \
    safe_free(_prm_->residue);             \
    safe_free(_prm_->rk);                  \
    safe_free(_prm_->req);                 \
    safe_free(_prm_->tk);                  \
    safe_free(_prm_->teq);                 \
    safe_free(_prm_->pk);                  \
    safe_free(_prm_->pn);                  \
    safe_free(_prm_->phase);               \
    safe_free(_prm_->solty);               \
    safe_free(_prm_->cn1);                 \
    safe_free(_prm_->cn2);                 \
    FREE_ParmBond(_prm_->pbondH);          \
    FREE_ParmBond(_prm_->pbond);           \
    FREE_ParmAngle(_prm_->pangleH);        \
    FREE_ParmAngle(_prm_->pangle);         \
    FREE_ParmDihedral(_prm_->pdihedralH);  \
    FREE_ParmDihedral(_prm_->pdihedral);   \
    safe_free(_prm_->natex);               \
    safe_free(_prm_->ag);                  \
    safe_free(_prm_->bg);                  \
    safe_free(_prm_->hbcut);               \
    FREE_Box(_prm_->box);                  \
    FREE_Cap(_prm_->cap);                  \
    FREE_Pert(_prm_->pert);                \
    FREE_Bond(_prm_->bond);                \
    FREE_Angle(_prm_->angle);              \
    FREE_Dihedral(_prm_->dihedral);        \
    FREE_Restart(_prm_->restart);          \
    safe_free(_prm_->lestyp);              \
    safe_free(_prm_->lesfac);              \
    safe_free(_prm_->lescnum);             \
    safe_free(_prm_->lessubsp);            \
  } safe_free(_prm_); _prm_ = NULL



#ifndef RDPARM_MODULE
#  ifdef __STDC__

extern void     clearParm(Parm *);
extern Parm *   openParm( char * );
extern void     parmInfo( );
extern void     printLJ( );
extern void     printAtomTypes( );
extern void     printExcluded( );
extern void     printAtomInfo(char *);
extern void     countAtoms(char *);
extern void     chargeOnAtoms(char *);
extern void     printAngles(char *);
extern void     printBonds(char *);
extern void     printPerturbedBonds(char *);
extern void     printPerturbedAngles(char *);
extern void     printPerturbedDihedrals(char *);
extern void     printDelphiCharge(char *);
extern void     printDihedrals(char *);
extern void     PrintDihedral(int, FILE *, int *);
extern void     dumpDihedrals();
extern void     readParm();
extern void     translateBox(char *);
extern void     translateRestart(char *);
extern void     writeParm( FILE *, int );
extern void     verbosity( );
extern void     initializeParm(Parm *);
extern void     getParm(char *);
extern void     putParm(char *);
extern void     quit();
extern void     doMardi2Sander(char *);
extern void     modifyBoxInfo();
extern void     modifyMolInfo();
extern int      getAtomTypes(Name, double *, double *, int);
extern int      getAtomCharge(Name, int, double *, int);
extern int      getAtomRadii(Name, int, double *, int);
extern void     testit();
extern void     checkVelocity();
extern void     deleteBAD(char *, int);
extern void     deletePerturbedBAD(char *, int);
extern void     restrainBAD(char *);
extern void     prnlevSet(char *);

#  else

extern void     clearParm();
extern Parm *   openParm();
extern void     parmInfo();
extern void     printLJ();
extern void     printAtomTypes( );
extern void     countAtoms( );
extern void     chargeOnAtoms( );
extern void     printExcluded( );
extern void     printAtomInfo();
extern void     printAngles();
extern void     printBonds();
extern void     printPerturbedBonds();
extern void     printPerturbedAngles();
extern void     printPerturbedDihedrals();
extern void     printDelphiCharge();
extern void     printDihedrals();
extern void     PrintDihedral();
extern void     dumpDihedrals();
extern void     readParm();
extern void     translateBox();
extern void     translateRestart();
extern void     writeParm();
extern void     verbosity();
extern void     getParm();
extern void     initializeParm();
extern void     putParm();
extern void     quit();
extern void     doMardi2Sander();
extern void     modifyBoxInfo();
extern void     modifyMolInfo();
extern int      getAtomTypes();
extern int      getAtomCharge();
extern int      getAtomRadii();
extern void     testit();
extern void     checkVelocity();
extern void     deleteBAD();
extern void     deletePerturbedBAD();
extern void     restrainBAD();
extern void     prnlevSet();

#  endif

/* global variables defined in rdparm */
extern int verboseStatus;
extern int isModifiedParm;
extern int isWrittenParm;
extern int isErrorParm;
extern Parm *parm;

#endif /* RDPARM_MODULE */

