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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/experimental.c,v 10.0 2008/04/15 23:24:11 case Exp $
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
#include <ctype.h>
#include <math.h>

#define EXPERIMENTAL_MODULE
#include "ptraj.h"

/*
 *  this source file contains work in progress, hence the name "experimental"
 */

#define WATER_ATOMS 3
#define WATER_BONDS 3
#define WATER_EXCLUDED 4

     

/* count contiguous waters from the end of the molecule array...
 */

void countWaters(Parm *p, int *rwaters, int *rstart_atom, int *rstart_res)
{
  int i;
  int waters, start_atom, start_res;

  if ( ! parm->IFBOX ) {
    warning("countWaters()", "No box information in parm file!\n");
    return;
  }

  waters = 0;
  start_atom = p->NTOTAT;
  start_res = p->NTOTRS;
  for ( i = p->box->nspm; i >= p->box->nspsol; i-- ) {
    if ( p->box->nsp[i-1] == WATER_ATOMS ) {
      start_atom -= WATER_ATOMS;
      start_res -= 1;
      waters++;
      if ( strcmp(parm->residue[start_res].labres, "WAT ") != 0 ) {
	start_atom += WATER_ATOMS;
	start_res += 1;
	waters--;
	break;
      }
    } else
      break;
  }

  *rwaters = waters;
  *rstart_atom = start_atom;
  *rstart_res = start_res;
}


/* add or remove waters from the current parm file */

   Parm *
modifyTIP3P( int waters ) 
{
  Parm *newparm;
  int start_atom, start_res, natoms, nres;
  int old_waters;
  int i, j, len;
  double chargeO, chargeH1, chargeH2;
  double massO, massH1, massH2;
  Name igraphO, igraphH1, igraphH2;
  Name isymblO, isymblH1, isymblH2;
  Name itreeO,  itreeH1,  itreeH2;
  int iacO, iacH1, iacH2;
  int bondOH, bondHH;

  if ( ! parm->IFBOX ) {
    warning("modifyTIP3P()", 
	    "No box information in parm file; cannot strip waters!\n");
    return( (Parm *) NULL);
  }

  countWaters(parm, &old_waters, &start_atom, &start_res);

  fprintf(stdout, "Assuming the current parm file has %i waters\n", 
	  old_waters);
  fprintf(stdout, "starting at atom %i (residue %i)\n",
	  start_atom+1, start_res+1);

  if ( waters == old_waters ) {
    fprintf(stdout, "No waters need to be added or removed, returning...\n");
    return( (Parm *) NULL);
  } else if ( waters < old_waters ) {
    fprintf(stdout, "Attempting to remove %i waters\n", old_waters - waters);
  } else {
    fprintf(stdout, "Attempt to add %i waters\n", waters - old_waters);
  }

  natoms = start_atom + WATER_ATOMS * waters; 
  nres = start_res + waters;

  fprintf(stdout, "The first water starts with atom %i\n", start_atom+1);
  fprintf(stdout, "The first water residue is %i\n", start_res+1);

  fprintf(stdout, "The new parm file will have %i atoms and %i residues\n",
	  natoms, nres);

  fprintf(stdout, 
	  "Grabbing information to characterize the water parameters...\n");

  chargeO  = parm->atom[start_atom].chrg;
  chargeH1 = parm->atom[start_atom+1].chrg;
  chargeH2 = parm->atom[start_atom+2].chrg;

  massO  = parm->atom[start_atom].amass;
  massH1 = parm->atom[start_atom+1].amass;
  massH2 = parm->atom[start_atom+2].amass;

  strcpy(igraphO,  parm->atom[start_atom].igraph);
  strcpy(igraphH1, parm->atom[start_atom+1].igraph);
  strcpy(igraphH2, parm->atom[start_atom+2].igraph);

  strcpy(isymblO,  parm->atom[start_atom].isymbl);
  strcpy(isymblH1, parm->atom[start_atom+1].isymbl);
  strcpy(isymblH2, parm->atom[start_atom+2].isymbl);

  strcpy(itreeO,  parm->atom[start_atom].itree);
  strcpy(itreeH1, parm->atom[start_atom+1].itree);
  strcpy(itreeH2, parm->atom[start_atom+2].itree);

  iacO  = parm->atom[start_atom].iac;
  iacH1 = parm->atom[start_atom+1].iac;
  iacH2 = parm->atom[start_atom+2].iac;

  bondOH = parm->pbondH[parm->NBONH - old_waters * WATER_BONDS].icb;
  bondHH = parm->pbondH[parm->NBONH - old_waters * WATER_BONDS + 2].icb;

  fprintf(stdout, "Name:       %s      %s      %s\n", 
	  igraphO, igraphH1, igraphH2);
  fprintf(stdout, "Type:       %s      %s      %s\n", 
	  isymblO, isymblH1, isymblH2);
  fprintf(stdout, "Tree:       %s      %s      %s\n", 
	  itreeO,  itreeH1,  itreeH2);
  fprintf(stdout, "Charge: %9.5f %9.5f %9.5f\n", 
	  chargeO/18.222, chargeH1/18.222, chargeH2/18.222);
  fprintf(stdout, "Mass:   %9.5f %9.5f %9.5f\n", massO, massH1, massH2);
  fprintf(stdout, "IAC:        %4i      %4i     %4i\n", 
	  iacO, iacH1, iacH2);
  fprintf(stdout, "OW - HW bond pointer is %i (rk is %6.3f, tk is %6.3f)\n",
	  bondOH, parm->rk[bondOH], parm->req[bondOH]);
  fprintf(stdout, "HW - HW bond pointer is %i (rk is %6.3f, tk is %6.3f)\n",
	  bondHH, parm->rk[bondHH], parm->req[bondHH]);
  
  fprintf(stdout, "\nInitializing a new parm topology...\n");
  newparm = safe_malloc( sizeof( Parm ) );
  initializeParm( newparm );

  newparm->filename = "newparm";
  newparm->fp = NULL;
  if ( strncpy(newparm->title, parm->title, TITLE_LENGTH) == NULL )
    error("modifyTIP3P", "...scanning title from %s\n", parm->filename);
  
  newparm->NTOTAT = natoms;
  newparm->NTYPES = parm->NTYPES;
  newparm->NBONH  = parm->NBONH + ( waters - old_waters ) * WATER_BONDS;
  newparm->NBONA  = parm->NBONA;
  newparm->NTHETH = parm->NTHETH;
  newparm->NTHETA = parm->NTHETA;
  newparm->NPHIH  = parm->NPHIH;
  newparm->NPHIA  = parm->NPHIA;
  newparm->JHPARM = parm->JHPARM;
  newparm->JPARM  = parm->JPARM;
  newparm->NEXT   = parm->NEXT + ( waters - old_waters ) * WATER_EXCLUDED;
  newparm->NTOTRS = nres;
  newparm->MBONA  = parm->MBONA;
  newparm->MTHETS = parm->MTHETS;
  newparm->MPHIA  = parm->MPHIA;
  newparm->MUMBND = parm->MUMBND;
  newparm->MUMANG = parm->MUMANG;
  newparm->MPTRA  = parm->MPTRA;
  newparm->NATYP  = parm->NATYP;
  newparm->NHB    = parm->NHB;
  newparm->IFPERT = parm->IFPERT;
  if ( newparm->IFPERT ) {
    warning("modifyTIP3P()", 
	    "Code has not been verified with perturbation.\n");
    warning("modifyTIP3P()",
	    "Beware if there are any waters being perturbed!\n");
  }
  newparm->NBPER  = parm->NBPER;
  newparm->NGPER  = parm->NGPER;
  newparm->NDPER  = parm->NDPER;
  newparm->MBPER  = parm->MBPER;
  newparm->MGPER  = parm->MGPER;
  newparm->MDPER  = parm->MDPER;
  newparm->IFBOX  = parm->IFBOX;
  newparm->NMXRS  = parm->NMXRS;
  newparm->IFCAP  = parm->IFCAP;


  fprintf(stdout, "Updating old and new atom information...\n");

  /* allocate space for the atoms in rdparm format */
  newparm->atom = safe_malloc( sizeof( Atom ) * natoms );

  /* copy the "old" atom information to the new parm */
  memcpy(newparm->atom, parm->atom, (size_t) sizeof (Atom) * start_atom);

  /* add in the information for the waters */
  for (j = start_res, i = start_atom; 
       i < newparm->NTOTAT; 
       i += 3, j++) {
    newparm->atom[i].chrg   = chargeO;
    newparm->atom[i+1].chrg = chargeH1;
    newparm->atom[i+2].chrg = chargeH2;

    newparm->atom[i].amass   = massO;
    newparm->atom[i+1].amass = massH1;
    newparm->atom[i+2].amass = massH2;

    newparm->atom[i].chrg   = chargeO;
    newparm->atom[i+1].chrg = chargeH1;
    newparm->atom[i+2].chrg = chargeH2;

    newparm->atom[i].iac   = iacO;
    newparm->atom[i+1].iac = iacH1;
    newparm->atom[i+2].iac = iacH2;

    newparm->atom[i].numex   = 2;
    newparm->atom[i+1].numex = 1;
    newparm->atom[i+2].numex = 1;

    strcpy(newparm->atom[i].igraph,   igraphO);
    strcpy(newparm->atom[i+1].igraph, igraphH1);
    strcpy(newparm->atom[i+2].igraph, igraphH2);

    strcpy(newparm->atom[i].isymbl,   isymblO);
    strcpy(newparm->atom[i+1].isymbl, isymblH1);
    strcpy(newparm->atom[i+2].isymbl, isymblH2);

    strcpy(newparm->atom[i].itree,   itreeO);
    strcpy(newparm->atom[i+1].itree, itreeH1);
    strcpy(newparm->atom[i+2].itree, itreeH2);

    newparm->atom[i].res   = j;
    newparm->atom[i+1].res = j;
    newparm->atom[i+2].res = j;
  }

  newparm->nno = safe_malloc(sizeof(int)* newparm->NTYPES * newparm->NTYPES);
  memcpy(newparm->nno, parm->nno, (size_t)
	sizeof(int) * newparm->NTYPES * newparm->NTYPES);

  /* allocate residue information in rdparm format */
  newparm->residue = safe_malloc( sizeof( Residue ) * (nres+1) );

  /* copy old residue crap to new */
  memcpy(newparm->residue, parm->residue, (size_t) sizeof( Residue ) * start_res);


  for (i=start_atom, j=start_res; j < nres; i+= WATER_ATOMS, j++) {
    strcpy(newparm->residue[j].labres, "WAT ");
    newparm->residue[j].ipres = i+1;
    newparm->atom[i].res   = j;
    newparm->atom[i+1].res = j;
    newparm->atom[i+2].res = j;
  }
  newparm->residue[nres].ipres = 
    newparm->residue[nres-1].ipres + WATER_ATOMS;


  fprintf(stdout, "Copying parameter arrays...\n");
  fprintf(stdout, "RK...\n");
  newparm->rk  = (double *) safe_malloc(sizeof(double) * newparm->MUMBND);
  memcpy(newparm->rk,  parm->rk, (size_t) sizeof(double) * newparm->MUMBND);

  fprintf(stdout, "REQ...\n");
  newparm->req = (double *) safe_malloc(sizeof(double) * newparm->MUMBND);
  memcpy(newparm->req, parm->req, (size_t) sizeof(double) * newparm->MUMBND);

  fprintf(stdout, "TK...\n");
  newparm->tk  = (double *) safe_malloc(sizeof(double) * newparm->MUMANG);
  memcpy(newparm->tk,  parm->tk,  (size_t) sizeof(double) * newparm->MUMANG);

  fprintf(stdout, "TEQ...\n");
  newparm->teq = (double *) safe_malloc(sizeof(double) * newparm->MUMANG);
  memcpy(newparm->teq, parm->teq, (size_t) sizeof(double) * newparm->MUMANG);
  
  fprintf(stdout, "PK...\n");
  newparm->pk    = (double *) safe_malloc(sizeof(double) * newparm->MPTRA);
  memcpy(newparm->pk, parm->pk, (size_t) sizeof(double) * newparm->MPTRA);

  fprintf(stdout, "PN...\n");
  newparm->pn    = (double *) safe_malloc(sizeof(double) * newparm->MPTRA);
  memcpy(newparm->pn,  parm->pn, (size_t) sizeof(double) * newparm->MPTRA);

  fprintf(stdout, "PHASE...\n");
  newparm->phase = (double *) safe_malloc(sizeof(double) * newparm->MPTRA);
  memcpy(newparm->phase, parm->phase, (size_t) sizeof(double) * newparm->MPTRA);
  
  fprintf(stdout, "Copying solty array...\n");

  newparm->solty = (double *) safe_malloc(sizeof(double) * newparm->NATYP );
  memcpy(newparm->solty, parm->solty, (size_t) sizeof(double) * newparm->NATYP);
  
  fprintf(stdout, "Copying CN1 and CN2 arrays...\n");

  len = newparm->NTYPES;
  len = len * (len + 1) / 2;
  newparm->cn1 = (double *) safe_malloc( sizeof( double) * len );
  newparm->cn2 = (double *) safe_malloc( sizeof( double) * len );
  memcpy(newparm->cn1, parm->cn1, (size_t) sizeof(double) * len);
  memcpy(newparm->cn2, parm->cn2, (size_t) sizeof(double) * len);
  
  fprintf(stdout, "Updating bond arrays...\n");

  newparm->pbondH = safe_malloc( sizeof(ParmBond) * newparm->NBONH );
  memcpy(newparm->pbondH, parm->pbondH, (size_t)
	sizeof(ParmBond) * (parm->NBONH - old_waters * WATER_BONDS));
  for(i = parm->NBONH - old_waters * WATER_BONDS, j = start_atom;
      i < newparm->NBONH;
      i += 3, j+= 3) {
    
    newparm->pbondH[i].ib = j*3;
    newparm->pbondH[i].jb = (j+1)*3;
    newparm->pbondH[i].icb = bondOH;
    
    newparm->pbondH[i+1].ib = j*3;
    newparm->pbondH[i+1].jb = (j+2)*3;
    newparm->pbondH[i+1].icb = bondOH;

    newparm->pbondH[i+2].ib = (j+1)*3;
    newparm->pbondH[i+2].jb = (j+2)*3;
    newparm->pbondH[i+2].icb = bondHH;
  }

  newparm->pbond = safe_malloc( sizeof( ParmBond ) * newparm->MBONA );
  memcpy(newparm->pbond, parm->pbond, (size_t) sizeof(ParmBond) * newparm->MBONA);
  
  newparm->pangleH = safe_malloc( sizeof( ParmAngle ) * newparm->NTHETH );
  memcpy(newparm->pangleH, parm->pangleH, (size_t) sizeof(ParmAngle) * newparm->NTHETH);
  
  newparm->pangle = safe_malloc( sizeof( ParmAngle ) * newparm->MTHETS );
  memcpy(newparm->pangle, parm->pangle, (size_t) sizeof(ParmAngle) * newparm->MTHETS);
  
  newparm->pdihedralH = safe_malloc( sizeof( ParmDihedral ) * newparm->NPHIH );
  memcpy(newparm->pdihedralH, parm->pdihedralH, 
	(size_t) sizeof(ParmDihedral) * newparm->NPHIH);
  
  newparm->pdihedral = safe_malloc( sizeof( ParmDihedral ) * newparm->MPHIA );
  memcpy(newparm->pdihedral, parm->pdihedral,
	(size_t) sizeof(ParmDihedral) * newparm->MPHIA);
  
  /* The excluded atom list is hardcoded for WAT; LES will not
   * work if water is within an LES group
   */

  fprintf(stdout, "Working on excluded atom list...\n");
  
  newparm->natex = safe_malloc( sizeof( int ) * newparm->NEXT );
  for (j=1, i=0; i < start_atom; i++)
    j += newparm->atom[i].numex;
  j--;
  memcpy(newparm->natex, parm->natex, (size_t) sizeof(int) * j);

  /* for the excluded atom list, for the 1st atom we want
   * atoms 2 and 3, and for the 2nd atom 3 and for the 3rd 0,
   */
  for (i = start_atom; i < natoms; i += 3) {
    newparm->natex[j++] = i+2;
    newparm->natex[j++] = i+3;
    newparm->natex[j++] = i+3;
    newparm->natex[j++] = 0;
  }
  
  newparm->ag = safe_malloc( sizeof( double ) * newparm->NHB );
  newparm->bg = safe_malloc( sizeof( double ) * newparm->NHB );
  newparm->hbcut = safe_malloc( sizeof( double ) * newparm->NHB );
  memcpy(newparm->ag, parm->ag, (size_t) sizeof(double) * newparm->NHB);
  memcpy(newparm->bg, parm->bg, (size_t) sizeof(double) * newparm->NHB);
  memcpy(newparm->hbcut, parm->hbcut, (size_t) sizeof(double) * newparm->NHB);

  
  fprintf(stdout, "Join and rotate information...\n");

  for (i=start_atom; i < natoms; i += 3) {
    if ( i == start_atom ) 
      newparm->atom[i].join = i-2;
    else
      newparm->atom[i].join = newparm->atom[i-1].join;
    newparm->atom[i+1].join = i+1;
    newparm->atom[i+2].join = i+1;

    newparm->atom[i].irotat = i+3;
    newparm->atom[i+1].irotat = i+2;
    newparm->atom[i+2].irotat = i+3;
  }

  fprintf(stdout, "Box information...\n");

  newparm->box = safe_malloc( sizeof( Box ) );
  newparm->box->iptres = parm->box->iptres;
  newparm->box->nspm = parm->box->nspm - old_waters + waters;
  newparm->box->nspsol = parm->box->nspsol;
  newparm->box->nsp = safe_malloc( sizeof( int ) * newparm->box->nspm );
  memcpy(newparm->box->nsp, parm->box->nsp, 
	(size_t) sizeof(int) * (parm->box->nspm - old_waters));
  for (i = parm->box->nspm - old_waters; i < newparm->box->nspm; i++) {
    newparm->box->nsp[i] = WATER_ATOMS;
  }

  newparm->box->beta = parm->box->beta;
  newparm->box->box[0] = parm->box->box[0];
  newparm->box->box[1] = parm->box->box[1];
  newparm->box->box[2] = parm->box->box[2];

  if ( parm->IFPERT ) {

    fprintf(stdout, "Perturbation information...\n");

    newparm->pert = safe_malloc( sizeof( Pert ) );
    if ( parm->NBPER ) {
      newparm->pert->ibper = safe_malloc( sizeof( int ) * newparm->NBPER );
      newparm->pert->jbper = safe_malloc( sizeof( int ) * newparm->NBPER );
      newparm->pert->icbper = 
	safe_malloc( sizeof( int ) * 2 * newparm->NBPER );
      memcpy(newparm->pert->ibper, parm->pert->ibper, 
	    (size_t) sizeof(int)*newparm->NBPER);
      memcpy(newparm->pert->jbper, parm->pert->jbper, 
	    (size_t) sizeof(int)*newparm->NBPER);
      memcpy(newparm->pert->icbper, parm->pert->icbper, 
	    (size_t) sizeof(int)* 2 * newparm->NBPER);
    }
    if ( parm->NGPER ) {
      newparm->pert->itper = safe_malloc( sizeof( int ) * newparm->NGPER );
      newparm->pert->jtper = safe_malloc( sizeof( int ) * newparm->NGPER );
      newparm->pert->ktper = safe_malloc( sizeof( int ) * newparm->NGPER );
      parm->pert->ictper = safe_malloc( sizeof( int ) * 2 * newparm->NGPER );
      memcpy(newparm->pert->itper, parm->pert->itper, 
	    (size_t) sizeof(int)*newparm->NGPER);
      memcpy(newparm->pert->jtper, parm->pert->jtper, 
	    (size_t) sizeof(int)*newparm->NGPER);
      memcpy(newparm->pert->ktper, parm->pert->ktper, 
	    (size_t) sizeof(int)*newparm->NGPER);
      memcpy(newparm->pert->ictper, parm->pert->ictper, 
	    (size_t) sizeof(int) * 2 * newparm->NGPER);
    }
    if ( parm->NDPER ) {
      newparm->pert->ipper = safe_malloc( sizeof( int ) * newparm->NDPER ); 
      newparm->pert->jpper = safe_malloc( sizeof( int ) * newparm->NDPER );
      newparm->pert->kpper = safe_malloc( sizeof( int ) * newparm->NDPER );
      newparm->pert->lpper = safe_malloc( sizeof( int ) * newparm->NDPER );
      newparm->pert->icpper = 
	safe_malloc( sizeof( int ) * 2 * newparm->NDPER );
      memcpy(newparm->pert->ipper, parm->pert->ipper, 
	    (size_t) sizeof(int)*newparm->NDPER);
      memcpy(newparm->pert->jpper, parm->pert->jpper, 
	    (size_t) sizeof(int)*newparm->NDPER);
      memcpy(newparm->pert->kpper, parm->pert->kpper, 
	    (size_t) sizeof(int)*newparm->NDPER);
      memcpy(newparm->pert->lpper, parm->pert->lpper, 
	    (size_t) sizeof(int)*newparm->NDPER);
      memcpy(newparm->pert->icpper, parm->pert->icpper, 
	    (size_t) sizeof(int) * 2 * newparm->NDPER);
    }

    newparm->pert->labper = safe_malloc( sizeof( Name ) * newparm->NTOTRS );
    memcpy(newparm->pert->labper, parm->pert->labper, 
	  (size_t) sizeof(Name)* start_res);
    for (i=start_res; i < nres; i++)
      strcpy(newparm->pert->labper[i], "WAT ");

    newparm->pert->igrper = safe_malloc( sizeof( Name ) * newparm->NTOTAT );
    newparm->pert->ismper = safe_malloc( sizeof( Name ) * newparm->NTOTAT );
    memcpy(newparm->pert->igrper, parm->pert->igrper, (size_t) start_atom);
    memcpy(newparm->pert->ismper, parm->pert->ismper, (size_t) start_atom);
    for (i=start_atom; i < natoms; i+=3) {
      strcpy(parm->pert->igrper[i],   igraphO);
      strcpy(parm->pert->igrper[i+1], igraphH1);
      strcpy(parm->pert->igrper[i+2], igraphH2);

      strcpy(parm->pert->ismper[i],   isymblO);
      strcpy(parm->pert->ismper[i+1], isymblH1);
      strcpy(parm->pert->ismper[i+2], isymblH2);
    }

    newparm->pert->almper = safe_malloc( sizeof( double ) * newparm->NTOTAT );
    memcpy(newparm->pert->almper, parm->pert->almper, 
	  (size_t) sizeof(double) * parm->NTOTAT);

    newparm->pert->iaper  = safe_malloc( sizeof( int )    * newparm->NTOTAT );
    memcpy(newparm->pert->iaper, parm->pert->iaper, 
	  (size_t) sizeof(int)*parm->NTOTAT);
    for (i=parm->NTOTAT; i < newparm->NTOTAT; i++)
      newparm->pert->iaper[i] = 0;

    newparm->pert->iacper = safe_malloc( sizeof( int ) * newparm->NTOTAT );
    memcpy(newparm->pert->iacper, parm->pert->iacper, 
	  (size_t) sizeof(int)*parm->NTOTAT);
    for (i=parm->NTOTAT; i < newparm->NTOTAT; i++)
      newparm->pert->iacper[i] = newparm->atom[i].iac;

    newparm->pert->cgper  = safe_malloc( sizeof( double ) * newparm->NTOTAT );
    memcpy(newparm->pert->cgper, parm->pert->cgper, 
	  (size_t) sizeof(double) * parm->NTOTAT);
    for (i=parm->NTOTAT; i < newparm->NTOTAT; i++)
      parm->pert->cgper[i] = parm->atom[i].chrg;
  }

  /* LES crap */
  if (parm->JPARM) {

    newparm->nlestyp = parm->nlestyp;
    newparm->lestyp  = (int *) safe_malloc(sizeof(int) * parm->NTOTAT);
    newparm->lesfac  = (double *) safe_malloc(sizeof(double) * parm->nlestyp*parm->nlestyp);
    newparm->lescnum = (int *) safe_malloc(sizeof(int) * parm->NTOTAT);
    newparm->lessubsp= (int *) safe_malloc(sizeof(int) * parm->NTOTAT);

    for (i=0; i < parm->NTOTAT; i++) {
      newparm->lestyp[i] = parm->lestyp[i];
      newparm->lescnum[i] = parm->lescnum[i];
      newparm->lessubsp[i] = parm->lessubsp[i];
    }
    for (i=0; i < parm->nlestyp * parm->nlestyp; i++) 
      newparm->lesfac[i] = parm->lesfac[i];
  }


  fprintf(stdout, "Fixing up internal bond, angle and dihedral info...\n");

  /* fix up bond, angle, and dihedrals */
  newparm->bond = safe_malloc( sizeof( Bond ) * (newparm->NBONH + 
						 newparm->MBONA + 
						 newparm->NBPER) );
  for (i=0; i < newparm->NBONH; i++) {
    newparm->bond[i].atom[0] = unObfuscateAtom(newparm->pbondH[i].ib);
    newparm->bond[i].atom[1] = unObfuscateAtom(newparm->pbondH[i].jb);
    newparm->bond[i].rk = newparm->rk[newparm->pbondH[i].icb-1];
    newparm->bond[i].req = newparm->req[newparm->pbondH[i].icb-1];
    newparm->bond[i].scale = 1.0;
  }
  for (i=0; i < newparm->MBONA; i++) {
    newparm->bond[i+newparm->NBONH].atom[0] = 
      unObfuscateAtom(newparm->pbond[i].ib);
    newparm->bond[i+newparm->NBONH].atom[1] = 
      unObfuscateAtom(newparm->pbond[i].jb);
    newparm->bond[i+newparm->NBONH].rk = 
      newparm->rk[newparm->pbond[i].icb-1];
    newparm->bond[i+newparm->NBONH].req = 
      newparm->req[newparm->pbond[i].icb-1];
    newparm->bond[i+newparm->NBONH].scale = 1.0;
  }

  newparm->angle = safe_malloc( sizeof( Angle ) *
			       (newparm->NTHETH + newparm->MTHETS) );
  for (i=0; i < newparm->NTHETH; i++) {
    newparm->angle[i].atom[0] = unObfuscateAtom(newparm->pangleH[i].it);
    newparm->angle[i].atom[1] = unObfuscateAtom(newparm->pangleH[i].jt);
    newparm->angle[i].atom[2] = unObfuscateAtom(newparm->pangleH[i].kt);
    newparm->angle[i].tk = newparm->tk[newparm->pangleH[i].ict-1];
    newparm->angle[i].teq = newparm->teq[newparm->pangleH[i].ict-1];
    newparm->angle[i].scale = 1.0;                          
  }
  for (i=0; i < newparm->MTHETS; i++) {
    newparm->angle[i+newparm->NTHETH].atom[0] = 
      unObfuscateAtom(newparm->pangle[i].it);
    newparm->angle[i+newparm->NTHETH].atom[1] = 
      unObfuscateAtom(newparm->pangle[i].jt);
    newparm->angle[i+newparm->NTHETH].atom[2] = 
      unObfuscateAtom(newparm->pangle[i].kt);
    newparm->angle[i+newparm->NTHETH].tk = 
      newparm->tk[newparm->pangle[i].ict-1];   
    newparm->angle[i+newparm->NTHETH].teq = 
      newparm->teq[newparm->pangle[i].ict-1];
    newparm->angle[i+newparm->NTHETH].scale = 1.0;                          
  }

  newparm->dihedral = safe_malloc( sizeof( Dihedral ) *
				  (newparm->NPHIH + newparm->MPHIA) );
  for (i=0; i < newparm->NPHIH; i++) {
    newparm->dihedral[i].atom[0] = unObfuscateAtom(newparm->pdihedralH[i].ip);
    newparm->dihedral[i].atom[1] = unObfuscateAtom(newparm->pdihedralH[i].jp);
    newparm->dihedral[i].atom[2] = unObfuscateAtom(newparm->pdihedralH[i].kp);
    newparm->dihedral[i].atom[3] = unObfuscateAtom(newparm->pdihedralH[i].lp);
    newparm->dihedral[i].pk = newparm->pk[newparm->pdihedralH[i].icp-1];
    newparm->dihedral[i].pn = newparm->pn[newparm->pdihedralH[i].icp-1];
    newparm->dihedral[i].phase = newparm->phase[newparm->pdihedralH[i].icp-1];
    newparm->dihedral[i].scale = 1.0;                          
  }
  for (i=0; i < newparm->MPHIA; i++) {
    newparm->dihedral[i+newparm->NPHIH].atom[0] = 
      unObfuscateAtom(newparm->pdihedral[i].ip);
    newparm->dihedral[i+newparm->NPHIH].atom[1] = 
      unObfuscateAtom(newparm->pdihedral[i].jp);        
    newparm->dihedral[i+newparm->NPHIH].atom[2] = 
      unObfuscateAtom(newparm->pdihedral[i].kp);
    newparm->dihedral[i+newparm->NPHIH].atom[3] = 
      unObfuscateAtom(newparm->pdihedral[i].lp);
    newparm->dihedral[i+newparm->NPHIH].pk = 
      newparm->pk[newparm->pdihedral[i].icp-1];   
    newparm->dihedral[i+newparm->NPHIH].pn = 
      newparm->pn[newparm->pdihedral[i].icp-1];   
    newparm->dihedral[i+newparm->NPHIH].phase = 
      newparm->phase[newparm->pdihedral[i].icp-1];   
    newparm->dihedral[i+newparm->NPHIH].scale = 1.0;                          
  }


  return( newparm );

}

   void
testWater() 
{
  char buffer[BUFFER_SIZE];
  char *response;
  int waters, new_waters;
  int start_atom, start_res;
  Parm *newparm, *tmpparm;
  FILE *fpout;
  char *filenameReturn = NULL;

  if ( parm->IFBOX) 
    waters = parm->box->nspm - parm->box->nspsol + 1;
  else {
    waters = 0;
    fprintf(stdout, "No box information; returning...\n");
    return;
  }


  countWaters(parm, &waters, &start_atom, &start_res);

  fprintf(stdout, "This parm file appears to have %i contiguous waters\n", 
	  waters);
  fprintf(stdout, "starting at residue %i, atom %i\n", 
	  start_res+1, start_atom+1);
  if ( ! promptUserResponse(stdin, stdout,
			    "Is this correct? [yes] ", "yes", 1) ) {
    fprintf(stdout, 
	    "\nHmmmn.  This implies that either the waters are not\n");
    fprintf(stdout, 
	    "contiguous or the molecule array is somehow hosed.\n");
    fprintf(stdout, 
	    "If the molecule array is hosed, try using the\n");
    fprintf(stdout, 
	    "the \"modifymoleculeinfo\" command...  This program will\n");
    fprintf(stdout, "currently only handle contiguous waters...\n");
  }

  fprintf(stdout, "This routine will add or remove waters from the\n");
  fprintf(stdout, "current parm file.\n\n");
  sprintf(buffer, "How many waters for the modified parm file?  [%i] ", 
	  waters);
  response = promptUser(stdin, stdout, buffer);
		
  if (response == NULL || sscanf(response, "%i", &new_waters) != 1)
    new_waters = waters;
  safe_free( (void *) response );

  fprintf(stdout, "Attempting to modify parm file to have %i waters...\n",
	  new_waters);

  newparm = modifyTIP3P(new_waters);

  filenameReturn = promptToOpenFile(&fpout, "", "w",
		   "Input the name for the modified parm/topology file: ");

  tmpparm = parm;
  parm = newparm;
  writeParm( fpout, 1 ); 
  parm = tmpparm;

  safe_free( (void *) filenameReturn );

}



   int
findAtom(Name *atomName, int residue)
{
  int j;

  for (j = parm->residue[residue-1].ipres;
       j < parm->residue[residue].ipres; j++) {
    if ( isMatchResidue(parm->atom[j-1].igraph, atomName) == 1 ) {
      return j;
    }
  }
  return -1;
}
		


  
   void
checkGrid(char *gridFilename)
{
  FILE *gridFile;
  char *filename;
  char buffer[BUFFER_SIZE];
  int gridnx, gridlox, gridhix;
  int gridny, gridloy, gridhiy;
  int gridnz, gridloz, gridhiz;
  int *grid;
  int index;
  int i, j, k;
  double value;
  double box[6];
  int *histogram;
  int maxGrid, maxHistogram;
  int maxx, maxy, maxz;
  int total;

  filename = promptToOpenFile(&gridFile, gridFilename, "r",
			      "Input the name of XPLOR format grid file: ");
  
  /*
   * read first line of junk, then number of title lines
   * to read...
   */

  if (fgets(buffer, BUFFER_SIZE, gridFile) == NULL ||
      fgets(buffer, BUFFER_SIZE, gridFile) == NULL) {
    fprintf(stderr, "Unexpected EOF in %s\n", filename);
    return;
  }

  if (sscanf(buffer, "%i", &i) != 1) {
    fprintf(stderr, "Error scanning number of title lines%s\n", 
	    filename);
    return;
  }
  for (j = 0; j < i; j++) {
    if (fgets(buffer, BUFFER_SIZE, gridFile) == NULL) { 
      fprintf(stderr, "Unexpected EOF in %s\n", filename);
      return;
    }
  }
  fprintf(stdout, "Final title line is\n%s", buffer);

  /*
   * read grid size and box information
   */
  if (fgets(buffer, BUFFER_SIZE, gridFile) == NULL ||
      sscanf(buffer, "%i%i%i%i%i%i%i%i%i", 
	     &gridnx, &gridlox, &gridhix,
	     &gridny, &gridloy, &gridhiy,
	     &gridnz, &gridloz, &gridhiz) != 9 ) {
    fprintf(stdout, "Error scanning grid extents from %s\n", filename);
    return;
  }
  fprintf(stdout, "Grid size is %i by %i by %i\n", gridnx, gridny, gridnz);
  fprintf(stdout, "Grid runs from %i to %i by %i to %i by %i to %i\n",
	  gridlox, gridhix, gridloy, gridhiy, gridloz, gridhiz);

  grid = safe_malloc(sizeof(int) * gridnx * gridny * gridnz);
  for (i=0; i < gridnx * gridny * gridnz; i++) {
    grid[i] = 0;
  }
  
  if (fgets(buffer, BUFFER_SIZE, gridFile) == NULL ||
      sscanf(buffer, "%lf%lf%lf%lf%lf%lf",
	      &box[0], &box[1], &box[2], &box[3], &box[4], &box[5]) != 6 ) {
    fprintf(stdout, "Error scanning grid box information from %s\n",
	    filename);
    return;
  }

  /*
   * read the "ZYX"
   */
  if (fgets(buffer, BUFFER_SIZE, gridFile) == NULL ||
      strstr(buffer, "ZYX") == NULL) {
    fprintf(stdout, "%s", buffer);
    fprintf(stdout, "Error scanning ZYX keyword from %s\n", filename);
    return;
  }

  /*
   * read the slices
   */ 
  maxGrid = 0;
  total = 0;
  for (k = 0; k < gridnz; k++) {
    
    if (fscanf(gridFile, "%i", &i) != 1) {
      fprintf(stdout, "Error scanning slice marker for slice %i\n",
	      gridloz+k);
      return;
    }
    
    for (j = 0; j < gridny; j++) {
      for (i = 0; i < gridnx; i++) {
	if (fscanf(gridFile, "%lf", &value) != 1) {
	  fprintf(stdout, "Scanning %i %i %i, value is %i\n",
		  i, j, k, (int) value);
	  fprintf(stdout, "Error scanning grid value...\n");
	  return;
	}
	index = i * gridny * gridnz + j * gridnz + k;
	grid[index] = value;
	total += value;
	if ( grid[index] > maxGrid ) {
	  maxGrid = grid[index];
	  maxx = i;
	  maxy = j;
	  maxz = k;
	}
      }
    }
  }

  /*
   * histogram the grid
   */

  histogram = safe_malloc(sizeof(int) * maxGrid);
  for (i = 0; i < maxGrid; i++) {
    histogram[i] = 0;
  }

  for (k = 0; k < gridnz; k++) {
    for (j = 0; j < gridny; j++) {
      for (i = 0; i < gridnx; i++) {
	
	index = i * gridny * gridnz + j * gridnz + k;
	if (grid[index] > 0) {
	  histogram[ grid[index]-1 ]++;
	}
      }
    }
  }

  /*
   * dump the histogram
   */
  fprintf(stdout, "The maximum value of the grid is %i\n", maxGrid);

  maxHistogram = 0;
  for (i = 0; i < maxGrid; i++) {
    if (histogram[i] > maxHistogram) maxHistogram = histogram[i];
  }
  for (i = 0; i < maxGrid; i++) {
    fprintf(stdout, "%5i (%5i): ", i+1, histogram[i]);
    if (histogram[i] > 0) {
      for (j = 1; j < histogram[i]; j += (maxHistogram/60 + 1)) {
	fprintf(stdout, "*");
      }
    }
    fprintf(stdout, "\n");
  }
      
  for (i=0; i < maxGrid; i++) {
    for (j=i+1; j < maxGrid; j++) {
      histogram[i] += histogram[j];
    }
  }


  fprintf(stdout, "Dumping contour levels...\n");
  maxHistogram = 0;
  for (i = 0; i < maxGrid; i++) {
    if (histogram[i] > maxHistogram) maxHistogram = histogram[i];
  }
  maxHistogram /= 10;
  for (i = 0; i < maxGrid; i++) {
    fprintf(stdout, "%5i (%5i): ", i+1, histogram[i]);
    if (histogram[i] > 0) {
      if (histogram[i] > maxHistogram) {
	fprintf(stdout, "X");
      } else {
	for (j = 1; j < histogram[i]; j += (maxHistogram/60 + 1)) {
	  fprintf(stdout, "*");
	}
      }
    }
    fprintf(stdout, "\n");
  }

  fprintf(stdout, "Max grid value occurs at grid element (%i, %i, %i)\n",
	  maxx, maxy, maxz);
  fprintf(stdout, "Dumping nearby elements...\n");
  for (i = maxx-1; i <= maxx + 1; i++) {
    for (j = maxy-1; j <= maxy + 1; j++) {
      for (k = maxz-1; k <= maxz + 1; k++) {
	if ( i < gridnx && j < gridny && k < gridnz ) {

	  index = i * gridny * gridnz + j * gridnz + k;
	  fprintf(stdout, "(%4i, %4i, %4i) grid value is %i\n",
		  i, j, k, grid[index]);
	}
      }
    }
  }

  fprintf(stdout, "Occupancy estimate for 1 angstrom grid is %6.2f\n",
	  grid[maxx * gridny * gridnz + maxy * gridnz + maxz] +

	  grid[(maxx+1) * gridny * gridnz + maxy * gridnz + maxz] +
	  grid[(maxx-1) * gridny * gridnz + maxy * gridnz + maxz] +
	  grid[maxx * gridny * gridnz + (maxy+1) * gridnz + maxz] +
	  grid[maxx * gridny * gridnz + (maxy-1) * gridnz + maxz] +
	  grid[maxx * gridny * gridnz + maxy * gridnz + maxz+1] +
	  grid[maxx * gridny * gridnz + maxy * gridnz + maxz+1] +

	  0.3333333333 * (
          grid[(maxx+1) * gridny * gridnz + (maxy+1) * gridnz + maxz] +
	  grid[(maxx+1) * gridny * gridnz + (maxy-1) * gridnz + maxz] +
	  grid[(maxx-1) * gridny * gridnz + (maxy+1) * gridnz + maxz] +
	  grid[(maxx-1) * gridny * gridnz + (maxy-1) * gridnz + maxz] +
	  grid[(maxx+1) * gridny * gridnz + maxy * gridnz + maxz+1] +
	  grid[(maxx+1) * gridny * gridnz + maxy * gridnz + maxz-1] +
	  grid[(maxx-1) * gridny * gridnz + maxy * gridnz + maxz+1] +
	  grid[(maxx-1) * gridny * gridnz + maxy * gridnz + maxz-1] +
	  grid[maxx * gridny * gridnz + (maxy+1) * gridnz + maxz+1] +
	  grid[maxx * gridny * gridnz + (maxy+1) * gridnz + maxz-1] +
	  grid[maxx * gridny * gridnz + (maxy-1) * gridnz + maxz+1] +
	  grid[maxx * gridny * gridnz + (maxy-1) * gridnz + maxz-1] +

	  grid[(maxx+1) * gridny * gridnz + (maxy+1) * gridnz + maxz+1] +
	  grid[(maxx+1) * gridny * gridnz + (maxy+1) * gridnz + maxz-1] +
	  grid[(maxx+1) * gridny * gridnz + (maxy-1) * gridnz + maxz+1] +
	  grid[(maxx+1) * gridny * gridnz + (maxy-1) * gridnz + maxz-1] +
	  grid[(maxx-1) * gridny * gridnz + (maxy+1) * gridnz + maxz+1] +
	  grid[(maxx-1) * gridny * gridnz + (maxy+1) * gridnz + maxz-1] +
	  grid[(maxx-1) * gridny * gridnz + (maxy-1) * gridnz + maxz+1] +
	  grid[(maxx-1) * gridny * gridnz + (maxy-1) * gridnz + maxz-1]));

  fprintf(stdout, "The sum over all cells in the grid is %i\n", total);
  fprintf(stdout, "(divide this by #frames to get an estimate of the\n");
  fprintf(stdout, " number of atoms in the grid!)\n");

}










