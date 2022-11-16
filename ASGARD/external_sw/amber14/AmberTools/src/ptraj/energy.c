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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/energy.c,v 10.0 2008/04/15 23:24:11 case Exp $
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



#include <stdio.h>
#include <string.h>
#include <math.h>

#define ENERGY_MODULE
#include "ptraj.h"


   void
initializeAtomInfo(atomInfo *atom)
{
  atom->energy = 0.0;
  atom->force = 0.0;
  atom->bonds = 0;
  atom->angles = 0;
  atom->dihedrals = 0;
}



   rtfInfo *
loadEnergyInfoFromPrmtop(char *filename)
{
  FILE *fp;
  Parm *oldparm, *newparm;
  rtfInfo *rtf;

  int i, ai, aj;

  /*
   *  create a link to the current "global" parm
   */
  oldparm = parm;

  /*
   *  open up the file and set up a new parm
   */
  if ( openFile(&fp, filename, "r") == 0 )
    error("loadEnergyInfoFromPrmtop()", "Cannot open prmtop file %s", filename);
  newparm = (Parm *) safe_malloc(sizeof(Parm));
  newparm->fp = fp;
  newparm->filename = (char *) safe_malloc(sizeof(char) * (strlen(filename)+1));
  strcpy(newparm->filename, filename);
  parm = newparm;
  readParm();
  parm = oldparm;

  /*
   *  allocate rtfInfo *
   */
  rtf = (rtfInfo *) safe_malloc(sizeof(rtfInfo));
  INITIALIZE_rtfInfo(rtf);

  /*
   *  setup atom info
   */

  rtf->atoms = parm->NTOTAT;
  rtf->atom = (atomInfo *) safe_malloc(sizeof(atomInfo) * rtf->atoms);
  for (i=0; i<rtf->atoms; i++)
    initializeAtomInfo(&rtf->atom[i]);

  /*
   *  process bonds
   *    for each bond:
   *       copy over the parameters
   *       update the atom list which points back to the potential bonds
   *       zero out energy/force
   */

  rtf->bonds = newparm->NBONH + newparm->MBONA;
  rtf->bond = (Bond *) safe_malloc(sizeof(Bond) * rtf->bonds);
  for (i = 0; i < rtf->bonds; i++) {

    rtf->bond[i].rk = newparm->bond[i].rk;
    rtf->bond[i].req = newparm->bond[i].req;
    ai = newparm->bond[i].atom[0];
    rtf->bond[i].atom[0] = ai;
    aj = newparm->bond[i].atom[1];
    rtf->bond[i].atom[1] = aj;
    rtf->bond[i].scale = newparm->bond[i].scale;
    rtf->bond[i].active = 1;

    if ( rtf->atom[ai].bonds == MAX_BONDS_PER_ATOM ) {
      fprintf(stderr, "WARNING in loadEnergyInfoFromPrmtop().  Too many bonds per\n");
      fprintf(stderr, "atom were found for atom %i.  Max is %i\n.  Increase\n",
	      ai, MAX_BONDS_PER_ATOM);
      fprintf(stderr, "MAX_BONDS_PER_ATOM in energy.h and recompile\n");
    }
    if ( rtf->atom[aj].bonds == MAX_BONDS_PER_ATOM ) {
      fprintf(stderr, "WARNING in loadEnergyInfoFromPrmtop().  Too many bonds per\n");
      fprintf(stderr, "atom were found for atom %i.  Max is %i\n.  Increase\n",
	      aj, MAX_BONDS_PER_ATOM);
      fprintf(stderr, "MAX_BONDS_PER_ATOM in energy.h and recompile\n");
    }
	
    rtf->atom[ai].bondP[ rtf->atom[ai].bonds++ ] = i;
    rtf->atom[aj].bondP[ rtf->atom[aj].bonds++ ] = i;

    rtf->Ebond[i] = 0.0;

  }  

  return(rtf);


}

   rtfInfo *
loadEnergyInfoFromPSF()
{
  return(NULL);
}



   double
calculateBondEnergy(rtfInfo *rtf, double *x, double *y, double *z) 
{
  double rx, ry, rz, r2, r, delr;
  double totalE;

  int i, ai, aj;

  totalE = 0;
  for (i=0; i < rtf->bonds; i++) {

    if (rtf->bond[i].active) {

      ai = rtf->bond[i].atom[0];
      aj = rtf->bond[i].atom[1];

      rx = x[ai] - x[aj];
      ry = y[ai] - y[aj];
      rz = z[ai] - z[aj];

      r2 = rx*rx + ry*ry + rz*rz;
      r = sqrt(r2);
      delr = r - rtf->bond[i].req;

      rtf->Ebond[i] = rtf->bond[i].rk * delr * delr;
      totalE += rtf->Ebond[i];

      /*
      df = rtf->bond[i].rk * delr * 2.0 / r;
      fxi =  rx * df;
      fxj = -rx * df;
      fyi =  ry * df;
      fyj = -ry * df;
      fzi =  rz * df;
      fzj = -rz * df;
      */
    }
  }
  
  return (totalE);
}

