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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/energy.h,v 10.0 2008/04/15 23:24:11 case Exp $
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



/*
 *  This is the header file for energy.c
 */

#ifdef ENERGY_MODULE



#else


#endif



/*
 *  EXTERNALLY VISIBLE FUNCTION PROTOTYPES
 */

#define MAX_BONDS_PER_ATOM 6
#define MAX_ANGLES_PER_ATOM 6
#define MAX_DIHEDRALS_PER_ATOM 36



typedef struct _atomInfo {

  double energy;
  double force;
  int bondP[MAX_BONDS_PER_ATOM];
  int bonds;
  int angleP[MAX_ANGLES_PER_ATOM];
  int angles;
  int dihedralP[MAX_DIHEDRALS_PER_ATOM];
  int dihedrals;

} atomInfo;



typedef struct _rtfInfo {

  int atoms;
  atomInfo *atom;

  int bonds;
  Bond *bond;
  double *Ebond;


  double *aK;
  int    *ai;
  int    *aj;
  int    *ak;
  double *E_angle;

  double *dPhase;
  double *dPeriod;
  double *dPeak;
  int    *di;
  int    *dj;
  int    *dk;
  int    *dl;

} rtfInfo;


#define INITIALIZE_rtfInfo(_p_) \
   _p_->bonds = 0; \
   _p_->bond = NULL; \
   _p_->Ebond = NULL






#ifndef ENERGY_MODULE


#  ifdef __STDC__

rtfInfo *loadEnergyInfoFromPrmtop(char *);
double calculateBondEnergy(rtfInfo *, double *, double *, double *);

#  else  /* __STDC__ */

rtfInfo *loadEnergyInfoFromPrmtop();
double calculateBondEnergy();

#  endif /* __STDC__ */


#endif /* ENERGY_MODULE */

/*
 *  LOCAL STRUCTURES
 */


#ifdef ENERGY_MODULE




#endif /* ENERGY_MODULE */


