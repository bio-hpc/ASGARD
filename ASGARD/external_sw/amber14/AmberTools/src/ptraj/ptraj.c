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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/ptraj.c,v 10.45 2010/03/17 20:02:31 case Exp $
 *
 *  Revision: $Revision: 10.45 $
 *  Date: $Date: 2010/03/17 20:02:31 $
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
#include <ctype.h>
#include <assert.h>

#include <sys/stat.h>

#define PTRAJ_MODULE
#include "ptraj.h"


/*
 *  This is the main routine for the ptraj functionality and in general this
 *  code is likely not to be modified by general users/developers, except:
 *
 *  (*) Adding new "actions" or routines to process/analyze coordinate sets.
 *      In this case, ptrajSetup() will likely need modification for any
 *      "action" on the coordinates or ptrajSetupAnalyze() for any routine that
 *      proprocesses data accumulated by the actions.
 *  (*) Adding new coordinate formats.  Modify checkCoordinates(), 
 *      ptrajSetupIO(), printCoordinateInfo(), ptrajPreProcessInputCoordinates(),
 *      and ptrajProcessInputCoordinates()
 *  (*) Modifying the code to parse "mask" strings: all the *mask* routines
 *
 *  For more information on adding new actions, see the details comments in
 *  actions.c.
 */


   int
atomToResidue(int atom, int residues, int *ipres)
{
  int i;

  if (ipres == NULL) return -1;

  for (i = 0; i < residues; i++)
    if (atom >= ipres[i] && atom < ipres[i+1])
      return (i+1);

  return -1;
}



   int
atomToSolventMolecule(int atom, int molecules, int *start, int *stop)
{
  int i;

  if (start == NULL || stop == NULL) return -1;

  for (i = 0; i < molecules; i++)
    if (atom <= start[i])
      return -1;
    else if (atom > start[i] && atom <= stop[i])
      return (i+1);

  return -1;
}



   int
atomToMolecule(int atom, int molecules, int *mol)
{
  int i, a;

  if (mol == NULL) return -1;

  a = 0;
  for (i = 0; i < molecules; i++) {
    a += mol[i];
    if (atom <= a) 
      return (i+1);
  }
  return -1;
}






   int
isActiveDetailed(int atom, int residue, int *mask, int atoms, int residues,
		 Name *atomName, Name *residueName, int *ipres)
{
  int i;
 
  if (residue >= residues || residue < 0) {
    printf("WARNING: residue out of range in isActiveDetailed, res %i (total %i)\n",
	   residue, residues);
    return 0;
  }
  for (i = ipres[residue]-1; i < ipres[residue+1]-1; i++)
    if ( mask[i] && strcmp(atomName[i], atomName[atom]) == 0 ) 
      return 1;

  return 0;

}


   int
isActive(int atom, int residue, int *mask, ptrajState *state)
{
  return( isActiveDetailed(atom, residue, mask,
			   state->atoms, 
			   state->residues,
			   state->atomName,
			   state->residueName,
			   state->ipres) );
}



   int
isActiveResidueDetailed(int residue, int *mask, int atoms, int residues,
			Name *atomName, Name *residueName, int *ipres)
{
  int i, total;

  total = 0;
  if (residue >= residues || residue < 0) {
    printf("WARNING: residue out of range in isActiveResidueDetailed(), res %i (total %i)\n",
	   residue, residues);
    return 0;
  }
  for (i = ipres[residue]-1; i < ipres[residue+1]-1; i++)
    if ( ! mask[i] )
      return 0;
    else
      total++;

  return total;

}


   int
isActiveResidue(int residue, int *mask, ptrajState *state)
{
  return( isActiveResidueDetailed(residue, mask,
				  state->atoms, 
				  state->residues,
				  state->atomName,
				  state->residueName,
				  state->ipres) );
}



   void
printAtomMaskDetailed(FILE *file, int *mask, int atoms, int residues,
		      Name *atomName, Name *residueName, int *ipres)
{
  int i, j, curres;
  char tmpatom[20];
  int *resactive, *ressimilar;
  int printed, numactive, numresactive, numressimilar;
  int incurres, innextres;

  printed = 0;
  numactive = 0;
  numresactive = 0;
  numressimilar = 0;

  /*
   *  This routine is kind of junky since in general I want to avoid printing
   *  as much detail as possible.  Therefore, we check to see is certain ranges 
   *  residues have all atoms active, etc. to avoid printing each atom in a residue.
   *  This makes it ugly and obtuse.
   */


  if (mask == NULL) {
    fprintf(file, "[No atoms are selected]");
    return;
  }

  j=0;
  for (i=0; i < atoms; i++)
    if (mask[i]) j++;

  if (j == 0) {
    fprintf(file, "[No atoms are selected]");
    return;
  }

     /*
      *  check if all atoms are active and if so print an asterisk
      */

  j = 0;
  for (i=0; i < atoms; i++) {
    j += mask[i];
  }
  if ( j == atoms ) {
    fprintf(file, "  * (All atoms are selected)");
    return;
  }
  numactive = j;

     /*
      *  determine which residues have all the atoms in that residue active
      */

  resactive = (int *) safe_malloc(sizeof(int) * residues);
  for (i=0; i < residues; i++) {
    resactive[i] = 0.0;
  }

  curres = 0;
  j = 0;
  for (i=0; i < atoms; i++) {
    if (i == ipres[curres+1]-1) {
      if (j == ipres[curres+1] - ipres[curres]) {
	resactive[curres] = 1.0;
	numresactive++;
      }
      j = 0;
      curres++;
    }
    if (mask[i])
      j++;
  }
  if (j == ipres[curres+1] - ipres[curres]) {
    resactive[curres] = 1.0;
  }

     /*
      *  determine the range over which the residues are fully active
      */

  for (curres = residues-2; curres >= 0; curres--) {
    if (resactive[curres]) {
      resactive[curres] += resactive[curres+1];
      numresactive--;
    }
  }


  /*
   *  determine ranges over which residues have the same atoms active
   *  as the next residue
   */
  ressimilar = (int *) safe_malloc(sizeof(int) * residues);
  for (i=0; i < residues; i++) {
    ressimilar[i] = 0.0;
  }

  for (curres = residues-2; curres >=0; curres--) {

    incurres = 0;
    innextres = 0;
    for (i = ipres[curres]-1; i < ipres[curres+2]-1; i++) { /* check current and next residue */
      if ( mask[i] ) {
	incurres++;
      }

      if (isActiveDetailed(i, i<ipres[curres+1]-1 ? curres+1 : curres, mask,
			   atoms, residues, atomName, residueName, ipres))
	innextres++;

      if (incurres != innextres) /* select counterparts in next residues too! */
	break;
    }
    if (incurres && innextres == incurres) {
      ressimilar[curres] = ressimilar[curres+1] + 1;
    } else {
      numressimilar++;
    }

  }

   
     /*
      *  do the actual printing
      */

  j = 0;
  for (curres = 0; curres < residues; curres++) {

    if (resactive[curres] ) {

      /*
       *  If all of the atoms are active in this residue print either the
       *  residue number or range as appropriate
       */

      if (resactive[curres] > 2) {
	if (j!=0 && j%10 != 0) fprintf(file, ",");
	fprintf(file, ":%i-%i", curres+1, curres+resactive[curres]);
	curres += resactive[curres]-1;
      } else {
	if (j!=0 && j%10 != 0) fprintf(file, ",");
	fprintf(file, ":%i", curres+1);
      }
      j++;
      if (j != 0 && j % 10 == 0) {
	fprintf(file, "\n    ");
	j = 0;
      }
    } else if (ressimilar[curres]) {

      /*
       *  If there is a set of residues with a similar atom selection...
       */

      printed = 0;
      if (ressimilar[curres] >= 1) {
	if (j!=0 && j%10 != 0) fprintf(file, ",");
	fprintf(file, ":%i-%i", curres+1, curres+ressimilar[curres]+1);
	curres += ressimilar[curres];
      } else {
	if (j!=0 && j%10 != 0) fprintf(file, ",");
	fprintf(file, ":%i", curres+1);
      }

      for (i = ipres[curres]-1; i < ipres[curres+1]-1; i++) {
	if ( mask[i] ) {
	  if (printed) 
	    fprintf(file, ",");
	  else {
	    fprintf(file, "@");
	    printed = 1;
	  }
	  strcpy(tmpatom, atomName[i]);
	  fprintf(file, "%s", strtok(tmpatom, " "));
	}
      }
      j++;
      if (j != 0 && j % 10 == 0) {
	fprintf(file, "\n    ");
	j = 0;
      }


    } else {

      /*
       *  Print individual atoms
       */

      if (numactive > 10 && numressimilar > 10 && numresactive > 10) {
	fprintf(file, "\n    ");
	numactive = 0;
      }
      for (i = ipres[curres]-1; i < ipres[curres+1]-1; i++) {
	if ( mask[i] ) {
	  if (j!=0 && j%10 != 0) fprintf(file, ",");
	  
	  strcpy(tmpatom, atomName[i]);
	  fprintf(file, ":%i@%s", curres+1, strtok(tmpatom, " "));
	  j++;
	}
	if (j != 0 && j % 10 == 0) {
	  fprintf(file, "\n    ");
	  j = 0;
	}
      }
    }
  }
  /*
    if (j!=0 && j%10 != 0) fprintf(file, "\n");
  */
  safe_free(resactive);
  safe_free(ressimilar);
}


   void
printAtom(FILE *fpout, int atom, ptrajState *state)
{
  int curres;
  char buffer[50];


  curres = atomToResidue(atom+1, state->residues, state->ipres)-1;
  sprintf(buffer, ":%i               ", curres+1);
  sprintf(buffer+6, "%s %s", state->atomName[atom], state->residueName[curres]);
  fprintf(fpout, "%s", buffer);
}


   void
printAtomCompact(FILE *fpout, int atom, ptrajState *state)
{
  int curres;
  char buffer[50], *bufferp;

  curres = atomToResidue(atom+1, state->residues, state->ipres)-1;

  sprintf(buffer, ":%i", curres+1);
  bufferp = buffer+strlen(buffer);
  sprintf(bufferp, "@%s", state->atomName[atom]);

  fprintf(fpout, "%10s", buffer);


}

   void
printAtomCompact2(FILE *fpout, int atom, ptrajState *state)
{
  int curres;
  char buffer[50], *bufferp;

  curres = atomToResidue(atom+1, state->residues, state->ipres)-1;

  sprintf(buffer, ":%i", curres+1);
  bufferp = buffer+strlen(buffer);
  sprintf(bufferp, "@%s", state->atomName[atom]);

  fprintf(fpout, "%s", buffer);


}






   parseEntry *
parseToken(char **textp, int operator)
{
  char *text;
  int i, j;
  parseEntry *p;

  text = *textp;

  for (i=0;; i++) {
    if (text[i] == (char) 0 ||
	text[i] == ',' ||
	text[i] == ':' ||
	text[i] == '@' ||
	(i>0 && text[i] == '-' && !isalpha(text[i-1])) ||
	isspace(text[i]))
      break;
  }

  p = (parseEntry *) safe_malloc(sizeof(parseEntry));
  p->isnumber = 0;
  
  if ( i > 0 ) 
    p->token = (char *) safe_malloc(sizeof(char) * (i+1));

  for (j=0; j < i; j++) {
    p->token[j] = text[j];
  }
  p->token[i] = (char) 0;


  p->operator = PARSE_NOOP;
  switch ( text[i] ) {

  case '-':
    p->operator = PARSE_CONTINUATION;

  case ',':

    text++;

  }

  if ( isdigit( p->token[0] ) ) {
    p->isnumber = 1;
    for (j=1; j < i; j++) {
      if ( isdigit(p->token[j]) == 0 ) {
	p->isnumber = 0;
      }
    }
  }

  text = text+i;
  skipWhitespace(text);
  /*
   *  extra special check to handle extra spacing between continuation
   *  operators
   */
  if (text[0] == '-') {
    p->operator = PARSE_CONTINUATION;
    text += 1;
    skipWhitespace(text);
  } else if (text[0] == ',') {
    text += 1;
    skipWhitespace(text);
  }

  *textp = text;
  return(p);


}


   int
isMatch(char *s1, char *s2) 
{
  int i;

  /*
   *  straight match
   */

  if ( strcmp(s1, s2) == 0 ) return 1;

  /*
   *  fuzzy match: NOTE this will break if s1 > s2
   */

  if ( strchr(s1, '?') != NULL || strchr(s1, '*') != NULL ) {

    /*
     *  look over the minimal map between the two strings
     */
    for (i=0; i < strlen(s1) && i < strlen(s2); i++) {

      if ( s1[i] != s2[i] ) {
	switch( s1[i] ) {

	case '*':           /* wild card, multiple characters */
	  return 1;
	case '?':           /* wild card, single character    */
	  if ( isspace(s2[i]) ) return 0;
	  break;
	default:            /* mismatch                       */
	  return 0;
	}
      }
    }
    return 1;
  }
  return 0;
}



   void
sortIndex(double *values, int *index, int size) 
{

#ifdef SHELL_SORT
  double *sortValue;
  double approx_ln_to_log2, ftmp;
  int logsplit, itmp, psort;
  int i, j, k;

  approx_ln_to_log2 = 1.0 / log( (double) 2.0 ) + 0.000001;

  sortValue = (double *) safe_malloc(sizeof(double) * size);
  for (i=0; i < size; i++)
    sortValue[i] = values[i];

  logsplit = ( log( (double) size ) * approx_ln_to_log2 );
  i = size;

  for (psort = 1; psort <= logsplit; psort++) {
    i >>= 1;
    for (j = i+1; j <= size; j++) {
      k = j - i;
      ftmp = sortValue[j-1];
      itmp = index[j-1];
      while (k >= 1 && sortValue[k-1] > ftmp ) {
	sortValue[k+i-1] = sortValue[k-1];
	index[k+i-1] = index[k-1];
	k -= i;
      }
      sortValue[k+i-1] = ftmp;
      index[k+i-1] = itmp;
    }
  }

  safe_free(sortValue);
#else
  double *sortValue;
  double ftmp;
  int itmp;
  int i,j;
    
  sortValue = (double *) safe_malloc(sizeof(double) * size);
  for (i=0; i < size; i++)
    sortValue[i] = values[i];

  for (i=1; i < size; i++) {
    ftmp = sortValue[i];
    itmp = index[i];
    j = i-1;
    while (j >= 0 && sortValue[j] > ftmp) {
      sortValue[j+1] = sortValue[j];
      index[j+1] = index[j];
      j--;
    }
    sortValue[j+1] = ftmp;
    index[j+1] = itmp;
  }
  safe_free(sortValue);
#endif
}


   void
sortIndexFloat(float *values, int *index, int size) 
{

#ifdef SHELL_SORT
  float *sortValue;
  float approx_ln_to_log2, ftmp;
  int logsplit, itmp, psort;
  int i, j, k;

  approx_ln_to_log2 = 1.0 / log( (double) 2.0 ) + 0.000001;

  sortValue = (float *) safe_malloc(sizeof(float) * size);
  for (i=0; i < size; i++)
    sortValue[i] = values[i];

  logsplit = ( log( (double) size ) * approx_ln_to_log2 );
  i = size;

  for (psort = 1; psort <= logsplit; psort++) {
    i >>= 1;
    for (j = i+1; j <= size; j++) {
      k = j - i;
      ftmp = sortValue[j-1];
      itmp = index[j-1];
      while (k >= 1 && sortValue[k-1] > ftmp ) {
	sortValue[k+i-1] = sortValue[k-1];
	index[k+i-1] = index[k-1];
	k -= i;
      }
      sortValue[k+i-1] = ftmp;
      index[k+i-1] = itmp;
    }
  }

  safe_free(sortValue);
#else
  float *sortValue;
  float ftmp;
  int itmp;
  int i,j;
    
  sortValue = (float *) safe_malloc(sizeof(float) * size);
  for (i=0; i < size; i++)
    sortValue[i] = values[i];

  for (i=1; i < size; i++) {
    ftmp = sortValue[i];
    itmp = index[i];
    j = i-1;
    while (j >= 0 && sortValue[j] > ftmp) {
      sortValue[j+1] = sortValue[j];
      index[j+1] = index[j];
      j--;
    }
    sortValue[j+1] = ftmp;
    index[j+1] = itmp;
  }
  safe_free(sortValue);
#endif
}





/*
 *  The goal of this routine is to create a new ptrajState (newstate)
 *  based on the old ptrajState (oldstate) deleting atoms that are
 *  set in the mask array.  This will modify the value of newstate.
 */

#undef  ROUTINE
#define ROUTINE "modifyStateByMask()"

   void
modifyStateByMask(ptrajState **newstatep, ptrajState **oldstatep,
		  int *mask, int strip)
{
  int i, ires, isol, imol;
  int j, jres, jsol, jmol;
  int curres, cursol, curmol; 
  int k;
  Name *atomName, *residueName;
  double *charges, *masses;
  int *startsol, *stopsol, *ipres, *moleculeInfo;
  ptrajState *newstate, *oldstate;

  oldstate = *oldstatep;

     /*
      *  allocate space for the new state
      */
  newstate = (ptrajState *) safe_malloc(sizeof(ptrajState));
  INITIALIZE_ptrajState(newstate);

     /*
      *  allocate space for temporary arrays and perform initialization
      */

  atomName     = (Name *)   safe_malloc(sizeof(Name)   * oldstate->atoms);
  charges      = (double *) safe_malloc(sizeof(double) * oldstate->atoms);
  masses       = (double *) safe_malloc(sizeof(double) * oldstate->atoms);
  residueName  = (Name *)   safe_malloc(sizeof(Name)   * oldstate->residues);
  ipres        = (int *)    safe_malloc(sizeof(int)    * (oldstate->residues+1));
  if (oldstate->molecules) 
    moleculeInfo = (int *)    safe_malloc(sizeof(int)    * oldstate->molecules);
  if (oldstate->solventMolecules) {
    startsol = (int *) safe_malloc(sizeof(int) * oldstate->solventMolecules);
    stopsol  = (int *) safe_malloc(sizeof(int) * oldstate->solventMolecules);

    for (i=0; i < oldstate->solventMolecules; i++)
      startsol[i] = -1;
  }
  
  j = 0; 
  jres = -1; jsol = -1; jmol = -1;
  ires = -1; isol = -1; imol = -1;

  /*
   *  loop over all atoms and set up information for the newstate if the atom is 
   *  not to be deleted...
   */
  for (i=0; i < oldstate->atoms; i++) {

    curres = atomToResidue(i+1, oldstate->residues, oldstate->ipres)-1;

    if ( mask[i] == (strip ? 0 : 1) ) {
      /*
       *  this atom is not to be deleted
       */

         /*
          *  copy over atom information
          */
      strcpy(atomName[j], oldstate->atomName[i]);
      charges[j] = oldstate->charges[i];
      masses[j] = oldstate->masses[i];

         /*
          *  check to see if we are in the same residue or not
          *  and copy relevant information
          */
      if (ires == -1 || ires != curres) {
	jres++;
	strcpy(residueName[jres], oldstate->residueName[curres]);
	ipres[jres] = j+1;
	ires = curres;
	
      }

         /*
          *  deal with the molecule information
          */
      if (oldstate->molecules) {
	curmol = atomToMolecule(i+1, oldstate->molecules, oldstate->moleculeInfo)-1;
	if (imol == -1 || imol != curmol) {
	  jmol++;
	  moleculeInfo[jmol] = oldstate->moleculeInfo[curmol];
	  for (k=1; k < oldstate->moleculeInfo[curmol]; k++)
	    if (mask[i+k] == (strip ? 1 : 0))
	      moleculeInfo[jmol]--;
	  imol = curmol;
	}
      }

         /*
          *  deal with the solvent information
          */
      if (oldstate->solventMolecules) {

	cursol = atomToSolventMolecule(i+1, oldstate->solventMolecules, 
				       oldstate->solventMoleculeStart,
				       oldstate->solventMoleculeStop)-1;

	if (cursol >= 0 && (isol == -1 || isol != cursol)) {
	  jsol++;
	  startsol[jsol] = j;
	  stopsol[jsol] = j; 
	  for (k=i; k < oldstate->solventMoleculeStop[cursol]; k++)
	    if (mask[k] == (strip ? 0 : 1)) stopsol[jsol]++;
	  isol = cursol;
	}
      }
         /*
          *  increment the new atom counter
          */
      j++;
    }
  }

     /*
      *  fix up IPRES
      */
  ipres[jres+1] = j+1;

     /*
      *  set up the newstate using the date placed into the temporary arrays;
      *  free up unneeded space as we go...
      */
  newstate->atoms = j;
  newstate->masses   = (double *) safe_malloc(sizeof(double) * newstate->atoms);
  newstate->charges  = (double *) safe_malloc(sizeof(double) * newstate->atoms);
  newstate->atomName = (Name *)   safe_malloc(sizeof(Name)   * newstate->atoms);

  for (i=0; i < newstate->atoms; i++) {
    newstate->masses[i]  = masses[i];
    newstate->charges[i] = charges[i];
    strcpy(newstate->atomName[i], atomName[i]);
  }
  safe_free(masses);
  safe_free(charges);
  safe_free(atomName);

  newstate->residues = jres+1;
  newstate->ipres = (int *) safe_malloc(sizeof(int) * (newstate->residues+1));
  newstate->residueName = (Name *) safe_malloc(sizeof(Name) * newstate->residues);
  for (i=0; i < newstate->residues; i++) {
    newstate->ipres[i] = ipres[i];
    strcpy(newstate->residueName[i], residueName[i]);
  }
  newstate->ipres[i] = ipres[i];

  safe_free(ipres);
  safe_free(residueName);

  newstate->IFBOX = oldstate->IFBOX;
  if (jsol < 0) {
    jsol = 0;
    newstate->solventMolecules = jsol;
    newstate->solventMoleculeStart = NULL;
    newstate->solventMoleculeStop = NULL;
    newstate->solventMask = NULL;
  } else {
    newstate->solventMolecules = jsol+1;
    newstate->solventMoleculeStart = (int *) safe_malloc(sizeof(int) * (jsol+1));
    newstate->solventMoleculeStop  = (int *) safe_malloc(sizeof(int) * (jsol+1));
    for (i=0; i < newstate->solventMolecules; i++) {
      newstate->solventMoleculeStart[i] = startsol[i];
      newstate->solventMoleculeStop[i]  = stopsol[i];
    }
    newstate->solventMask = (int *) safe_malloc(sizeof(int) * newstate->atoms);
    newstate->solventAtoms = 0;

    cursol = 0;
    for (i=0; i < newstate->atoms; i++) {
      newstate->solventMask[i] = 0;

      if (i == newstate->solventMoleculeStart[cursol]) {
	for (j=i; j < newstate->solventMoleculeStop[cursol]; j++) {
	  newstate->solventAtoms++;
	  newstate->solventMask[j] = 1;
	}
	i = newstate->solventMoleculeStop[cursol]-1;
	cursol++;
      }
    }
  }

  if (oldstate->solventMolecules) {
    safe_free(startsol);
    safe_free(stopsol);
  }

  newstate->maxFrames = oldstate->maxFrames;
  if (oldstate->molecules) {
    newstate->molecules = jmol+1;
    newstate->moleculeInfo = (int *) safe_malloc(sizeof(int) * newstate->molecules);
    for (i=0; i<newstate->molecules; i++)
      newstate->moleculeInfo[i] = moleculeInfo[i];
    safe_free(moleculeInfo);
  }

  for (i=0; i<6; i++)
    newstate->box[i] = oldstate->box[i];

  *newstatep = newstate;
}




   void
printAtomMask(FILE *file, int *mask, ptrajState *state)
{
  printAtomMaskDetailed(file,
			mask, 
			state->atoms, 
			state->residues,
			state->atomName,
			state->residueName,
			state->ipres);
}




   void
printHBondMask(int *mask, int *maskH1, int *maskH2, int *maskH3, ptrajState *state)
{
  int i;

  fprintf(stdout, "    Atom#  Residue#  Name   --   Atom#  Residue#  Name\n");
  for (i=0; i < state->atoms; i++) {

    if (mask[i] != 0) {
      
      fprintf(stdout, "  %5i    %4i       %s ",
	      i+1, atomToResidue(i+1, state->residues, state->ipres), state->atomName[i]);
      
      if (maskH1 != NULL) {
	if (maskH1[i] >= 0) {
	  
	  fprintf(stdout, " -- %5i    %4i       %s ",
		  maskH1[i]+1, atomToResidue(maskH1[i]+1, state->residues, state->ipres),
		  state->atomName[maskH1[i]]);
	}
	if (maskH2[i] >= 0) {
	  
	  fprintf(stdout, "\n  %5i    %4i       %s ",
		  i+1, atomToResidue(i+1, state->residues, state->ipres), state->atomName[i]);
	  fprintf(stdout, " -- %5i    %4i       %s ",
		  maskH2[i]+1, atomToResidue(maskH2[i]+1, state->residues, state->ipres),
		  state->atomName[maskH2[i]]);
	}
	
	if (maskH3[i] >= 0) {
	  fprintf(stdout, "\n  %5i    %4i       %s ",
		  i+1, atomToResidue(i+1, state->residues, state->ipres), state->atomName[i]);
	  fprintf(stdout, " -- %5i    %4i      %s ",
		  maskH3[i]+1, atomToResidue(maskH3[i]+1, state->residues, state->ipres),
		  state->atomName[maskH3[i]]);
	}
      }
      fprintf(stdout, "\n");
    }
  }
}


#undef  ROUTINE
#define ROUTINE "printHBondInfo()"

   void
printHBondInfo(int numdonor, int *donor, int numacceptor, 
	       int *acceptor, int *acceptorH, ptrajState *state)
{
  int i;

  if (numdonor > 0 && donor != NULL) {
    fprintf(stdout, "  HBOND DONOR LIST:\n");
    fprintf(stdout, "    Atom#  Residue#  Name\n");
    for (i=0; i < numdonor; i++) {

      fprintf(stdout, "  %5i    %4i       %s\n",
	      donor[i]+1, atomToResidue(donor[i]+1, state->residues, state->ipres), 
	      state->atomName[donor[i]]);
    }
  }

  if (numacceptor > 0 && acceptor != NULL) {
    fprintf(stdout, "  HBOND ACCEPTOR LIST:\n");
    fprintf(stdout, "    Atom#  Residue#  Name   --   Atom#  Residue#  Name\n");
    for (i=0; i < numacceptor; i++) {

      fprintf(stdout, "  %5i    %4i       %s ",
	      acceptor[i]+1, atomToResidue(acceptor[i]+1, state->residues, state->ipres), 
	      state->atomName[acceptor[i]]);

      fprintf(stdout, " -- %5i    %4i       %s\n",
	      acceptorH[i]+1, atomToResidue(acceptorH[i]+1, state->residues, state->ipres),
	      state->atomName[acceptorH[i]]);
    }
  }
}


   int *
processAtomMaskDetailed( char *maskString, int atoms, int residues,
			Name *atomName, Name *residueName, int *ipres)
{
  int actualAtoms;
  int not = 0;
  char *maskp;
  char *tmp;

  int *residueMask;
  int residueMaskActive = 0;
  int *atomMask;
  int atomMaskActive = 0;
  int res1, res2;
  int atom1, atom2;
  int continuation = 0;
  int i, j;
  stackType *residueStack = NULL;
  stackType *atomStack = NULL;
  parseEntry *pp;

  Name name;

  maskp = maskString;
  skipWhitespace(maskp);

  /*
   *  allocate mask strings
   */
  atomMask = (int *) safe_malloc(sizeof(int) * atoms);
  residueMask = (int *) safe_malloc(sizeof(int) * residues);
  memset(atomMask, 0, sizeof(int) * atoms);
  memset(residueMask, 0, sizeof(int) * residues);

  /*
   *  special case, choose ALL atoms
   */
  if ( maskp[0] == (char) 0 || maskp[0] == '*' ) {
    for (i=0; i < atoms; i++)
      atomMask[i] = 1;
    goto clean_return;
  }


  /*
   *  get rid of all NOT characters "~"; only one is significant
   *  and set NOT status if one is found...
   */
  while ( (tmp = strchr(maskString, '~' )) != NULL ) {
    not = 1;
    tmp[0] = ' ';
  }

  /*
   *  check for error
   */
  if (strchr(maskp, ':') == NULL &&
      strchr(maskp, '@') == NULL) {
    fprintf(stdout, "WARNING: Error in mask string, no \"@\" or \":\" present (%s)\n", 
	    maskString);
    safe_free( atomMask );
    safe_free( residueMask );
    return NULL;
  }

  /*
   *  do the main "parsing"
   */
  skipWhitespace(maskp);
  while ( maskp[0] != (char) 0 ) {
    
    if ( maskp[0] == ':' ) {
      maskp++;
      for (;;) {

	skipWhitespace(maskp);
	pp = parseToken(&maskp, PARSE_RESIDUE);
	pushBottomStack(&residueStack, (void *) pp);
	if (maskp[0] == (char) 0 || 
	    maskp[0] == '@' || 
	    maskp[0] == ':') break;
      }
    }

    if ( maskp[0] == '@' ) {
      maskp++;
      for (;;) {

	skipWhitespace(maskp);
	pp = parseToken(&maskp, PARSE_ATOM);
	pushBottomStack(&atomStack, (void *) pp);
	if ( maskp[0] == (char) 0 || maskp[0] == ':' ) break;
	if ( maskp[0] == '@' )
	  maskp++;
      }
    }

    /*
     *  now process the atomStack and residueStack
     */


    if ( not ) {
      for (i=0; i < atoms; i++)
	atomMask[i] = 1;
      for (i=0; i < residues; i++)
	residueMask[i] = 1;
    }

    while ( residueStack != NULL ) {

      if ( continuation ) {
	res1 = res2;
	res2 = -1;
      }

      pp = (parseEntry *) popStack( &residueStack );
      if ( pp->isnumber ) {
	if ( sscanf(pp->token, "%i", &res2) != 1 ) {
	  fprintf(stdout, "WARNING: error parsing atom mask\n");
	  safe_free( atomMask );
	  safe_free( residueMask );
	  return NULL;
	}
	res2--;

	if (continuation) {
	  continuation = 0;
	  if (res1 < 0) res1 = 0;
	  if (res2 >= residues) res2 = residues-1;
	  if (res2 < res1) res2 = res1;

	  for (i = res1; i <= res2; i++) {
	    residueMask[i] = (not ? 0 : 1);
	    residueMaskActive = 1;
	  }
	} else {
	  residueMask[res2] = (not ? 0 : 1);
	  residueMaskActive = 1;
	}

	if ( pp->operator == PARSE_CONTINUATION )
	  continuation = 1;

      }	else {

	continuation = 0;
	strcpy(name, "    ");
	for (i=0; i < strlen(pp->token) && i < NAME_SIZE; i++) {
	  name[i] = pp->token[i];
	}
	name[NAME_SIZE-1] = (char) 0;

	for (i=0; i < residues; i++)
	  if ( isMatch(name, residueName[i]) ) {
	    residueMask[i] = (not ? 0 : 1);
	    residueMaskActive = 1;
	  } 
      }
      safe_free(pp->token);
      pp->token=NULL;
      safe_free(pp);
    }

    while ( atomStack != NULL ) {

      if ( continuation ) {
	atom1 = atom2;
	atom2 = -1;
      }

      pp = (parseEntry *) popStack( &atomStack );
      if ( pp->isnumber ) {
	if ( sscanf(pp->token, "%i", &atom2) != 1 ) {
	  fprintf(stdout, "WARNING: error parsing atom mask\n");
	  safe_free( atomMask );
	  safe_free( residueMask );
	  return NULL;
	}
	atom2--;

	if (continuation) {
	  continuation = 0;
	  if (atom1 < 0) atom1 = 0;
	  if (atom2 > atoms) atom2 = atoms-1;
	  if (atom2 < atom1) atom2 = atom1;

	  for (i = atom1; i <= atom2; i++) {
	    atomMask[i] = (not ? 0 : 1);
	    atomMaskActive = 1;
	  }
	} else {
	  atomMask[atom2] = (not ? 0 : 1);
	  atomMaskActive = 1;
	}

	if ( pp->operator == PARSE_CONTINUATION )
	  continuation = 1;

      }	else {

	continuation = 0;
	strcpy(name, "    ");
	for (i=0; i < strlen(pp->token) && i < NAME_SIZE; i++) {
	  name[i] = pp->token[i];
	}
	name[NAME_SIZE-1] = (char) 0;

	for (i=0; i < atoms; i++)
	  if ( isMatch(name, atomName[i]) ) {
	    atomMask[i] = (not ? 0 : 1);
	    atomMaskActive = 1;
	  }
      }
      safe_free(pp->token);
      pp->token = NULL;
      safe_free(pp);
    }
    
    if ( atomMaskActive && residueMaskActive ) {
      for (i=0; i < residues; i++) {

	if ( residueMask[i] == 0 ) {
	  for (j = ipres[i]-1; 
	       j < ipres[i+1]-1;
	       j++) {
	    atomMask[j] = 0;
	  }
	}
      }
    } else if ( residueMaskActive ) {
      for (i=0; i < residues; i++) {

	for (j = ipres[i]-1;
	     j < ipres[i+1] - 1;
	     j++) {
	  if ( residueMask[i] )
	    atomMask[j] = 1;
	  else if (not)
	    atomMask[j] = 0;
	}
      }
    }
    
    atomMaskActive = 0;
    residueMaskActive = 0;
  }



 clean_return:

  actualAtoms = 0;
  for (i=0; i < atoms; i++)
    if ( atomMask[i] ) actualAtoms++;


  if ( (tmp = strchr(maskString, '\n')) )
    tmp[0] = (char) 0;

  if (actualAtoms > 0) {
    fprintf(stdout, "Mask %s%s] represents %i atoms\n", 
	    (not ? "[~" : "["), maskString, actualAtoms);
  } else  {
    fprintf(stdout, "Mask %s%s] represents %i atoms ",
            (not ? "[~" : "["), maskString, actualAtoms);
    fprintf(stdout, "!!!NO ATOMS DETECTED!!!\n");
    safe_free(atomMask);
    atomMask = NULL;
  }

  safe_free( residueMask );
  return ( atomMask );

}


   int *
processAtomMaskDetailedVH( char *maskString, int atoms, int residues,
			   Name *atomName, Name *residueName, int *ipres, void *x, void *y, void *z, char type)
{
  char *charmask;
  int *intmask;
  int i, actualAtoms;

  intmask = (int *) safe_malloc( atoms * sizeof(int) );

  if ( strcmp(maskString, "") == 0) {
    maskString = strcpy(maskString, "*");
  }

  charmask = parseMaskString(maskString, atoms, residues, atomName, residueName, ipres, x, y, z, type);
  
  /*
   *  the new mask parsing routine returns a character array
   *  rather than integer; this makes more sense, however in the meantime
   *  we need to convert between the two.
   */

  actualAtoms = 0;
  for (i = 0; i < atoms; i++) {
    if (charmask[i] == 'T') {
      intmask[i] = 1;
      actualAtoms++;
    } else
      intmask[i] = 0;
  }

  if (actualAtoms > 0) {
    fprintf(stdout, "Mask [%s] represents %i atoms\n", 
	    maskString, actualAtoms);
  } else  {
    fprintf(stdout, "Mask [%s] represents %i atoms ",
            maskString, actualAtoms);
    fprintf(stdout, "!!!NO ATOMS DETECTED!!!\n");
    safe_free(intmask);
    intmask = NULL;
  }

  if (prnlev > 2) {
    fprintf(stdout, "Parsed mask string matches:\n");
    printAtomMaskDetailed(stdout, intmask, atoms, residues, atomName, residueName, ipres);
  }

  safe_free(charmask);
  return(intmask);

}


   int *
processAtomMask( char *maskString, ptrajState *state )
{

  if ( maskString && (strchr(maskString, '<') != NULL ||
		      strchr(maskString, '>') != NULL)) {

    if (referenceInfo) /* referenceInfo is a global variable? */
      return (processAtomMaskDetailedVH(maskString, state->atoms, state->residues, 
					state->atomName, state->residueName, state->ipres, 
					referenceInfo->x, referenceInfo->y, referenceInfo->z, 'd'));
   } else  
        /*
	 *  using NULL array for coordinate in distance comparison will ensure error, may be fix later. 
	 */
      return (processAtomMaskDetailedVH(maskString, state->atoms, state->residues, 
					state->atomName, state->residueName, state->ipres, 
					NULL, NULL, NULL, 'd'));
  return NULL;
}


   int *
processAtomMaskWrapper(char *buffer, ptrajState *state, int printit, int returnit)
{
  int *mask;

  mask = processAtomMask(buffer, state);
  if (printit) {
    printAtomMaskDetailed(stdout,
			  mask,
			  state->atoms,
			  state->residues,
			  state->atomName,
			  state->residueName,
			  state->ipres);
    fprintf(stdout, "\n");
  }

  if (returnit) {
    return(mask);
  } else {
    safe_free(mask);
    return(NULL);
  }

}



   int *
processAtomMaskWrapperOLD(char *buffer, int printit, int returnit)
{
  int *mask;
  int i;
  ptrajState *state;

  state = (ptrajState *) safe_malloc( sizeof(ptrajState) );
  INITIALIZE_ptrajState(state);

  state->atoms = parm->NTOTAT;
  state->atomName = (Name *) safe_malloc( sizeof(Name) * state->atoms);
  state->masses = (double *) safe_malloc(sizeof(double) * state->atoms);
  state->charges = NULL;
  for (i=0; i < state->atoms; i++) {
    strcpy(state->atomName[i], parm->atom[i].igraph);
    state->masses[i] = parm->atom[i].amass;
  }

  state->residues = parm->NTOTRS;
  state->residueName = (Name *) safe_malloc (sizeof(Name) * state->residues);
  for (i=0; i < state->residues; i++) {
    strcpy(state->residueName[i], parm->residue[i].labres);
  }
  state->ipres = (int *) safe_malloc(sizeof(int) * (state->residues+1));
  for (i=0; i <= state->residues; i++) {
    state->ipres[i] = parm->residue[i].ipres;
  }

  mask = processAtomMaskDetailed(buffer,
				 state->atoms,
				 state->residues,
				 state->atomName,
				 state->residueName,
				 state->ipres);

  if (printit) {
    printAtomMaskDetailed(stdout,
			  mask,
			  state->atoms,
			  state->residues,
			  state->atomName,
			  state->residueName,
			  state->ipres);
    fprintf(stdout, "\n");
  }

  if ( state != NULL ) {
    safe_free( state->atomName );
    safe_free( state->residueName );
    safe_free( state->ipres );
    safe_free( state->masses );
    safe_free( state );
  }

  if (returnit) {
    return(mask);
  } else {
    safe_free(mask);
    return(NULL);
  }

}




  scalarInfo *
scalarStackGetName(stackType **scalarStackP, char *name)
{
  stackType *s;
  scalarInfo *info, *match;

  match = NULL;
  for (s = *scalarStackP; s != NULL; s = s->next) {
    info = (scalarInfo *) s->entry;
    if ( strcmp(info->name, name) == 0 )
      match = info;
  }
  return (match);
}

  transformHBondInfo *
hbondInfoStackGetName(stackType **scalarStackP, char *name)
{
  stackType *s;
  transformHBondInfo *info, *match;

  match = NULL;
  for (s = *scalarStackP; s != NULL; s = s->next) {
    info = (transformHBondInfo *) s->entry;
    if ( strcmp(info->name, name) == 0 )
      match = info;
  }
  return (match);
}


  transformMatrixInfo *
matrixInfoStackGetName(stackType **matrixStackP, char *name)
{
  stackType *s;
  transformMatrixInfo *info, *match;

  match = NULL;
  for (s = *matrixStackP; s != NULL; s = s->next) {
    info = (transformMatrixInfo *) s->entry;
    if ( strcmp(info->name, name) == 0 )
      match = info;
  }
  return (match);
}

  modesInfo *
modesInfoStackGetName(stackType **modesStackP, char *name)
{
  stackType *s;
  modesInfo *info, *match;

  match = NULL;
  for (s = *modesStackP; s != NULL; s = s->next) {
    info = (modesInfo *) s->entry;
    if ( strcmp(info->name, name) == 0 )
      match = info;
  }
  return (match);
}

   void
boxToRecip(double box[6], double ucell[9], double recip[9])
{
  double u12x,u12y,u12z;
  double u23x,u23y,u23z;
  double u31x,u31y,u31z;
  double volume;

  ucell[0] = box[0];
  ucell[1] = 0.0;
  ucell[2] = 0.0;
  ucell[3] = box[1]*cos(DEGRAD*box[5]);
  ucell[4] = box[1]*sin(DEGRAD*box[5]);
  ucell[5] = 0.0;
  ucell[6] = box[2]*cos(DEGRAD*box[4]);
  ucell[7] = (box[1]*box[2]*cos(DEGRAD*box[3]) - ucell[6]*ucell[3]) / ucell[4];
  ucell[8] = sqrt(box[2]*box[2] - ucell[6]*ucell[6] - ucell[7]*ucell[7]);

  u23x = ucell[4]*ucell[8] - ucell[5]*ucell[7];
  u23y = ucell[5]*ucell[6] - ucell[3]*ucell[8];
  u23z = ucell[3]*ucell[7] - ucell[4]*ucell[6];
  u31x = ucell[7]*ucell[2] - ucell[8]*ucell[1];
  u31y = ucell[8]*ucell[0] - ucell[6]*ucell[2];
  u31z = ucell[6]*ucell[1] - ucell[7]*ucell[0];
  u12x = ucell[1]*ucell[5] - ucell[2]*ucell[4];
  u12y = ucell[2]*ucell[3] - ucell[0]*ucell[5];
  u12z = ucell[0]*ucell[4] - ucell[1]*ucell[3];
  volume=ucell[0]*u23x + ucell[1]*u23y + ucell[2]*u23z;

  recip[0] = u23x/volume;
  recip[1] = u23y/volume;
  recip[2] = u23z/volume;
  recip[3] = u31x/volume;
  recip[4] = u31y/volume;
  recip[5] = u31z/volume;
  recip[6] = u12x/volume;
  recip[7] = u12y/volume;
  recip[8] = u12z/volume;

}

/*
 *  Calculate the minimum possible distance between periodic images.
 *  This routine assumes that the ucell and recip information have
 *  previously been set via a call to boxToRecip
 */

   double
calculateMinImagedDistance2(double *box, double *ucell, double *recip,
			    double x1, double y1, double z1,
			    double x2, double y2, double z2,
			    int *rix, int *riy, int *riz, int origin)
{
  double X, Y, Z;
  double min, dist;
  double fx, fy, fz;
  double f2x, f2y, f2z;
  int ix, iy, iz, i;

  if (box == NULL)
    return(0.0);

  min = 100.0 * (box[0]*box[0]+box[1]*box[1]+box[2]*box[2]);

  if (prnlev > 6) {
    fprintf(stdout, "ATOM      0  XXX A1      1     %7.3f %7.3f %7.3f\n",
	    x1, y1, z1);
    fprintf(stdout, "ATOM      1  XXX A2      1     %7.3f %7.3f %7.3f\n",
	    x2, y2, z2);

  }

  fx = x1*recip[0] + y1*recip[1] + z1*recip[2];
  fy = x1*recip[3] + y1*recip[4] + z1*recip[5];
  fz = x1*recip[6] + y1*recip[7] + z1*recip[8];

  f2x = x2*recip[0] + y2*recip[1] + z2*recip[2];
  f2y = x2*recip[3] + y2*recip[4] + z2*recip[5];
  f2z = x2*recip[6] + y2*recip[7] + z2*recip[8];

  if (origin) {
    fx += 0.5;
    fy += 0.5;
    fz += 0.5;
    f2x += 0.5;
    f2y += 0.5;
    f2z += 0.5;
  }

  fx = fx - floor(fx);
  fy = fy - floor(fy);
  fz = fz - floor(fz);

  f2x = f2x - floor(f2x);
  f2y = f2y - floor(f2y);
  f2z = f2z - floor(f2z);

  if (prnlev > 6) {
    i = 2;
    for (ix = -1; ix <= 1; ix++) {
      for (iy = -1; iy <= 1; iy++) {
	for (iz = -1; iz <= 1; iz++) {

	  X = (fx+ix)*ucell[0] + (fy+iy)*ucell[3] + (fz+iz)*ucell[6];
	  Y = (fx+ix)*ucell[1] + (fy+iy)*ucell[4] + (fz+iz)*ucell[7];
	  Z = (fx+ix)*ucell[2] + (fy+iy)*ucell[5] + (fz+iz)*ucell[8];

	  fprintf(stdout, "ATOM    %3i  XXX B%-2i     1     %7.3f %7.3f %7.3f\n",
		  i, i, X, Y, Z);
	  i++;
	}
      }
    }
  }

  i = 1;
  for (ix = -1; ix <= 1; ix++) {
    for (iy = -1; iy <= 1; iy++) {
      for (iz = -1; iz <= 1; iz++) {

	X = ((fx+ix)*ucell[0] + (fy+iy)*ucell[3] + (fz+iz)*ucell[6]) - (f2x*ucell[0] + f2y*ucell[3] + f2z*ucell[6]);
	Y = ((fx+ix)*ucell[1] + (fy+iy)*ucell[4] + (fz+iz)*ucell[7]) - (f2x*ucell[1] + f2y*ucell[4] + f2z*ucell[7]);
	Z = ((fx+ix)*ucell[2] + (fy+iy)*ucell[5] + (fz+iz)*ucell[8]) - (f2x*ucell[2] + f2y*ucell[5] + f2z*ucell[8]);

	dist = X*X + Y*Y + Z*Z;

	if (prnlev > 6) {
	  fprintf(stdout, "  IMAGE FAMILIAR  distance %3i: %6.3f  (%5i %5i %5i)\n",
		  i++, dist, ix, iy, iz);
	}

	if (dist < min) {
	  min = dist;
	  *rix = ix;
	  *riy = iy;
	  *riz = iz;
	}
      }
    }
  }
  if (prnlev > 4) {
    fprintf(stdout, "  IMAGE FAMILIAR, min distance is %6.3f (%5i %5i %5i)\n",
	    min, *rix, *riy, *riz);
  }

  return(min);
}



/*
 *  This routine will calculate a distance squared performing imaging,
 *  such that the minimum imaged distance (squared) is returned.  If
 *  the box is orthorhomic, simply subtract off multiple of the box lengthes.
 *  For non-orthorhomic, this is a little more tricky since this procedure
 *  is no longer applicable.  In this case, the ucell, recip and closest2
 *  information is used.  "closest2" represents a cutoff distance, such that if
 *  the calculate distance**2 is less than this return it (without calculating
 *  all possible image distances).  All possible images in each direction are 
 *  investigated currently.  This should be made "smarter" to only search in the
 *  appropriate octant based on the angle values.
 */
   double
calculateDistance2(int i, int j, double *x, double *y, double *z, 
		   double *box, double *ucell, double *recip, double closest2, int noimage)
     /*
       i -- index for first atom
       j -- index for second atom
       x -- coordinate arrays
       y
       z
       box -- box coordinates
       ucell
       recip
       closest2 -- minimum distance between a periodic image 
     */
{
  double X, Y, Z, W;
  double fx, fy, fz, f2x, f2y, f2z, min, dist;
  int ix, iy, iz;

  X = x[i] - x[j];
  Y = y[i] - y[j];
  Z = z[i] - z[j];

  if (box == NULL || box[0] == 0.0 || noimage > 0)
    return(X*X + Y*Y + Z*Z);
    
  if (box[3] == 90.0 && box[4] == 90.0 && box[5] == 90.0) {
    /*
     *  DO ORTHORHOMIBIC IMAGING (this is faster!)
     */

       /*
        *  rid sign information
        */
    if ( X < 0 ) X = -X;
    if ( Y < 0 ) Y = -Y;
    if ( Z < 0 ) Z = -Z;

       /*
        *  rid multiples of the box length and image
        */
    while ( X > box[0] ) X = X - box[0];
    while ( Y > box[1] ) Y = Y - box[1];
    while ( Z > box[2] ) Z = Z - box[2];

       /*
        *  find shortest distance in the periodic reference
        */
    W = box[0] - X;
    if ( W < X ) X = W;

    W = box[1] - Y;
    if ( W < Y ) Y = W;

    W = box[2] - Z;
    if ( W < Z ) Z = W;

    return ( X*X + Y*Y + Z*Z );

  } else {

    /*
     *  NON-ORTHORHOMBIC CASE: find shortest distance in periodic reference
     *  This is a brute force check requiring up to 26 distance evaluations.
     *  It has been adapted to be smarter by returning the first distance that
     *  is shorter than the minimum possible distance between images.
     */

    fx = x[j]*recip[0] + y[j]*recip[1] + z[j]*recip[2];
    fy = x[j]*recip[3] + y[j]*recip[4] + z[j]*recip[5];
    fz = x[j]*recip[6] + y[j]*recip[7] + z[j]*recip[8];

    f2x = x[i]*recip[0] + y[i]*recip[1] + z[i]*recip[2];
    f2y = x[i]*recip[3] + y[i]*recip[4] + z[i]*recip[5];
    f2z = x[i]*recip[6] + y[i]*recip[7] + z[i]*recip[8];

    fx = fx - floor(fx);
    fy = fy - floor(fy);
    fz = fz - floor(fz);

    f2x = f2x - floor(f2x);
    f2y = f2y - floor(f2y);
    f2z = f2z - floor(f2z);


    X = (fx*ucell[0] + fy*ucell[3] + fz*ucell[6]) - (f2x*ucell[0] + f2y*ucell[3] + f2z*ucell[6]); 
    Y = (fx*ucell[1] + fy*ucell[4] + fz*ucell[7]) - (f2x*ucell[1] + f2y*ucell[4] + f2z*ucell[7]); 
    Z = (fx*ucell[2] + fy*ucell[5] + fz*ucell[8]) - (f2x*ucell[2] + f2y*ucell[5] + f2z*ucell[8]);
    
    min = X*X + Y*Y + Z*Z;

    if (closest2 != 0.0 && min < closest2) return (min);


    for (ix = -1; ix <= 1; ix++) {
      for (iy = -1; iy <= 1; iy++) {
	for (iz = -1; iz <= 1; iz++) {

	  if (! (ix == 0 && iy == 0 && iz == 0) ) {
	    X = ((fx+ix)*ucell[0] + (fy+iy)*ucell[3] + (fz+iz)*ucell[6]) - (f2x*ucell[0] + f2y*ucell[3] + f2z*ucell[6]);
	    Y = ((fx+ix)*ucell[1] + (fy+iy)*ucell[4] + (fz+iz)*ucell[7]) - (f2x*ucell[1] + f2y*ucell[4] + f2z*ucell[7]);
	    Z = ((fx+ix)*ucell[2] + (fy+iy)*ucell[5] + (fz+iz)*ucell[8]) - (f2x*ucell[2] + f2y*ucell[5] + f2z*ucell[8]);
	    dist = X*X + Y*Y + Z*Z;

	    if (prnlev > 3) {
	      printf("DISTANCE + %2i*X %2i*Y %2i*Z unit cells is %8.3f\n", ix, iy, iz, sqrt(dist));
	    }
	    if (dist < min) {
	      min = dist;
	      if (closest2 != 0.0 && min < closest2) 
		return(min);
	    }
	  }
	}
      }
    }
    return(min);
  }
}



   ptrajState **
ptrajCurrentState()
{
  return(&ptrajCurrentStatePointer);
}

   ptrajState *
ptrajCopyState(ptrajState **stateinp)
{
  ptrajState *state, *statein;
  int i;

  /*
   *  Make a copy of the current state
   */
  statein = *stateinp;

  state = (ptrajState *) safe_malloc( sizeof(ptrajState) );
  INITIALIZE_ptrajState(state);

  state->atoms = statein->atoms;
  state->atomName = (Name *) safe_malloc( sizeof(Name) * state->atoms);
  state->masses = (double *) safe_malloc(sizeof(double) * state->atoms);
  state->charges = (double *) safe_malloc(sizeof(double) * state->atoms);
  for (i=0; i < state->atoms; i++) {
    strcpy(state->atomName[i], statein->atomName[i]);
    state->masses[i] = statein->masses[i];
    state->charges[i] = statein->charges[i];
  }

  state->residues = statein->residues;
  state->residueName = (Name *) safe_malloc (sizeof(Name) * state->residues);
  for (i=0; i < state->residues; i++) {
    strcpy(state->residueName[i], statein->residueName[i]);
  }
  state->ipres = (int *) safe_malloc(sizeof(int) * (state->residues+1));
  for (i=0; i <= state->residues; i++) {
    state->ipres[i] = statein->ipres[i];
  }

  state->solventMolecules = statein->solventMolecules;
  state->solventAtoms     = statein->solventAtoms;

  if (statein->solventMolecules) {

    state->solventMask = (int *) safe_malloc(sizeof(int) * statein->atoms);
    for (i=0; i < state->atoms; i++)
      state->solventMask[i] = statein->solventMask[i];

    state->solventMoleculeStart = (int *) 
      safe_malloc(sizeof(int) * statein->solventMolecules);
    state->solventMoleculeStop  = (int *) 
      safe_malloc(sizeof(int) * statein->solventMolecules);
    for (i=0; i < state->solventMolecules; i++) {
      state->solventMoleculeStart[i] = statein->solventMoleculeStart[i];
      state->solventMoleculeStop[i]  = statein->solventMoleculeStop[i];
    }
  } else {
    state->solventMoleculeStart = NULL;
    state->solventMoleculeStop  = NULL;
    state->solventMask = NULL;
  }

  state->IFBOX = statein->IFBOX;
  for (i=0; i < 6; i++)
    state->box[i] = statein->box[i];

  if ( statein->molecules > 0 ) {

    state->molecules = statein->molecules;
    state->moleculeInfo = (int *)
      safe_malloc( sizeof(int) * state->molecules );
    for (i=0; i < state->molecules; i++) {
      state->moleculeInfo[i] = statein->moleculeInfo[i];
    }
  }

  return(state);
}



   void
checkAtomMask(char *buffer)
{
  ptrajState **statep;

  statep = ptrajCurrentState();
  processAtomMaskWrapper(buffer, *statep, 1, 0);
}



   int *
returnAtomMask(char *buffer)
{
  ptrajState **statep;

  statep = ptrajCurrentState();
  return( processAtomMaskWrapper(buffer, *statep, 0, 1) );
}


   void
atomMaskIsActive(int *mask, ptrajState *state, int *activep, int *firstp)
{
  int i, active, first;

  if (mask == NULL) {
    *activep = 0;
    return;
  }

  active = 0;
  for (i=0; i < state->atoms; i++) {
    if (mask[i] != 0) {
      active++;
      if (active == 1) first = i;
    }
  }
  *activep = active;
  *firstp = first;
}



#undef  ROUTINE
#define ROUTINE "ptrajPrintState()"

   void
ptrajPrintState(ptrajState *state)
{
  int i,j;
  int curres;

  printf("Dumping state information...\n");
  printf("  atoms:      %i\n", state->atoms);
  printf("  residues:   %i\n", state->residues);
  printf("  box length: %8.3f %8.3f %8.3f\n", state->box[0], state->box[1], state->box[2]);
  printf("  box angles: %8.3f %8.3f %8.3f\n", state->box[3], state->box[4], state->box[5]);
  printf("  molecules:  %i\n", state->molecules);
  printf("  max frames: %i\n", state->maxFrames);
  if (prnlev > 1) {
    printf("  ATOM   NAME  RESIDUE NAME  CHARGE    MASS\n");
    curres = 0;
    for (i=0; i < state->atoms; i++) {
    
      if (i == state->ipres[curres+1]-1) curres++;

      printf("  %8i %s %8i %s %8.3f %8.3f\n", 
	     i+1, state->atomName[i], 
	     curres+1, state->residueName[curres],
	     state->charges[i], state->masses[i]);
    }
    printf("  Molecule information (atoms in each molecule):\n");
    for (i=0; i < state->molecules; i++) {
      printf("  %5i", state->moleculeInfo[i]);
      if (i!=0 && i%10==0) printf("\n");
    }
    printf("\n");
  }
  if (state->solventMolecules > 0) {
    printf("  solvent molecules: %i (%i atoms)\n", 
	   state->solventMolecules, state->solventAtoms);
    printf("  solvent mask is: ");
    printAtomMask(stdout, state->solventMask, state);
    fprintf(stdout, "\n");
    if (prnlev > 1) {
      for (i=0; i < state->solventMolecules; i++) {
	printf("  SOLVENT %4i: ", i+1);
	for (j=state->solventMoleculeStart[i]; j < state->solventMoleculeStop[i]; j++) {
	  printf(" (%s %5i)", state->atomName[j], j+1);
	}
	printf("\n");
      }
    }
  }
}


   void
ptrajClearState(ptrajState **stateinp)
{
  ptrajState *state;

  state = *stateinp;

  if (state != NULL) {
    safe_free( state->atomName );
    safe_free( state->residueName );
    safe_free( state->ipres );
    safe_free( state->masses );
    safe_free( state->charges );

    safe_free( state->solventMoleculeStart );
    safe_free( state->solventMoleculeStop );
    safe_free( state->solventMask );

    safe_free( state->moleculeInfo );
    INITIALIZE_ptrajState(state);
    safe_free( state );
  }
}


/*
 *  checkCoordinates(): a routine to determine what type of coordinate
 *  file specified by "filename" is and how many coordinate frames the 
 *  file represents.  The file is opened up, checked, then the file
 *  is closed.  It is the responsibility of later routines to reopen
 *  the file.  [This requirement is so the file limit will not be blown
 *  in the case of a user processing a boatload of files; in this way,
 *  sequential file access is maintained.]
 *  
 *  On success, a "coordinateInfo *" will be returned with the filename,
 *  format and start/stop represented.  
 *
 *  On failure, NULL is returned.
 * NOTE: Clean up memory if NULL is returned
 */

#undef  ROUTINE
#define ROUTINE "checkCoordinates()"

   coordinateInfo *
checkCoordinates(char *filename, int totalAtoms)
{
  struct stat frame_stat;
  int returnValue;
  double *x = NULL;
  float junk1, junk2, junk3, junk4, junk5, junk6, junk7, junk8, junk9;
  int i, actualAtoms, binposeof;
  int isBox = 0;
  int numBoxCoords = 0;
  int isVelocity = 0;
  int lines_per_set;
  int start = 1;
  int stop = 1;
  int frame_lines, seekable;
  int maxi, sizeFound; // For large gzip file calc, Amber Traj
  long long int file_size, frame_size, tmpfsize;
  float *binposScratch;
  FILE *fp;
  pdb_record r;
  coordType type = COORD_UNKNOWN;
  char buffer1[BUFFER_SIZE], buffer2[BUFFER_SIZE], buffer3[BUFFER_SIZE];
  Restart *restrt;
  coordinateInfo *trajInfo;
  charmmTrajectoryInfo **charmmTrajectoryp = NULL;
  charmmTrajectoryInfo  *charmmTrajectory = NULL;

  trajInfo = (coordinateInfo *) safe_malloc(sizeof(coordinateInfo));
  INITIALIZE_coordinateInfo(trajInfo);
  trajInfo->filename = copyString(filename);
  trajInfo->accessMode = 0;           // READ access
  trajInfo->isMPI=0;                  // Use standard file ops for reading
  if (worldsize>1) trajInfo->isMPI=1; // Use MPI file ops for reading
  

  /* Initialize Buffers  */
  memset(buffer1,' ',BUFFER_SIZE);
  memset(buffer2,' ',BUFFER_SIZE);
  memset(buffer3,' ',BUFFER_SIZE);

  /* Identify file based on the hex signature. Will currently detect 
   * compressed files and NETCDF files. Compression is
   * dealt with by the openFile routine in io.c.
   */
  returnValue = id_Filesig(filename,NULL);
  switch (returnValue) {
    case 0: // Unknown file
      trajInfo->isNetcdf=0;
      break;
    case 1: // Gzip
    case 2: // Bzip
    case 3: // Zip
      trajInfo->isNetcdf=0;
      // Reading compressed files wont work with MPI file ops
      trajInfo->isMPI=0;
      trajInfo->compressType=returnValue;
      break;
    case 4:
      printfone("NETCDF file:\n");
#ifdef BINTRAJ
      trajInfo->isNetcdf=1;
      type = COORD_AMBER_NETCDF;
#else
      printfone("WARNING: Ptraj was compiled without NETCDF support. Recompile with -DBINTRAJ.\n");
      safe_free(trajInfo->filename);
      safe_free(trajInfo);
      return NULL;
#endif
      break;
    case -1: // Could not open file - not always a catastrophic error.
      if (prnlev>0) printfone("Could not open %s for signature identification.\n",filename);
      return NULL;
      break;
    case -2: // Id_filesig internal error
      printfone("Error identifying %s file signature.\n",filename);
      return NULL;
      break;
    default: // Should not get anything else!
      printfone("Unknown error in Id_filesig().\n");
      return NULL;
  }

  /* Attempt to open the trajectory file.There are currently 3 basic
   * types of file: FILE, MPI_File, and netcdf files. Pointers for each are
   * trajInfo->file, trajInfo->mfp,  or in trajinfo->NCinfo respectively. 
   */
  if (openTraj(trajInfo)!=0) {
    printfone("trajin %s ignored; could not open file (%s)\n", filename, filename);
    safe_free(trajInfo->filename);
    safe_free(trajInfo);
    return NULL;
  }

  /*
   *  Try to determine what "type" of file the input file is if it is
   *  -not- a NetCDF file.
   *
   *  NOTE: if new coordType's are implemented, this section may
   *  need to be modified.
   *
   *    (a) read two lines of text
   *    (b) check for pdb  [this requires reading both lines of
   *        text since the comment from an AMBER restrt or 
   *        trajectory may have text which "fools" the pdb routines...]
   *    (c) look for the BINPOS magic string
   *    (d) look for CHARMM trajectory magic string
   *    (e) if not a pdb or BINPOS, check the second line to see if it 
   *        has three numbers which implies it is an AMBER trajectory file.
   */
  if (trajInfo->isNetcdf==0) {
    fp=trajInfo->file; // This FILE pointer is needed for file types other than AMBER_TRAJECTORY

    if (trajFile_fgets(buffer1, BUFFER_SIZE, trajInfo) == NULL ||
        trajFile_fgets(buffer2, BUFFER_SIZE, trajInfo) == NULL ||
        trajFile_fgets(buffer3, BUFFER_SIZE, trajInfo) == NULL) {
      printfone("%s: EOF encountered prematurely (%s)\n", ROUTINE, filename);
      safe_free(trajInfo->filename);
      safe_free(trajInfo);
      return NULL;
    }

      /*
       *  If we have reached this point and the file is a NETCDF file 
       *  this is a compressed NETCDF which is currently not supported.
       */
    if (buffer1[0] == 'C' && 
	buffer1[1] == 'D' && 
	buffer1[2] == 'F') {
      printfone("    *** Compressed NetCDF trajectories are currently not supported by\n");
      printfone("    *** PTRAJ. Decompress %s before processing.\n",filename);
      printfone("ERROR: trajin %s ignored: Compressed NetCDF not supported.\n",filename);
      safe_free(trajInfo->filename);
      safe_free(trajInfo);
      return NULL;
    } 
      /*
       *  Is this file a PDB file?
       */
    r = pdb_read_string(buffer1);

    if (r.record_type != PDB_UNKNOWN) {
      r = pdb_read_string(buffer2);
      if (prnlev>0) printf("Reading string (%s) R = %i\n", buffer2, r.record_type);
      if (r.record_type != PDB_UNKNOWN)
	type = COORD_PDB;
    } else if (strncmp(buffer1, "REMARK", 6) == 0)
      /*
       *  this is a hack to allow CHARMM pdb files to be read or other PDB's 
       *  which have malformed remarks
       */
      type = COORD_PDB;


      /*
       *  Does this file have the magic header for the Scripps binary format?
       */
    if (buffer1[0] == 'f' && 
	buffer1[1] == 'x' && 
	buffer1[2] == 'y' &&
	buffer1[3] == 'z') {
      type = COORD_BINPOS;
    }

      /*
       *  Does file file have the magic header for a CHARMM trajectory???
       */
    if (buffer1[4] == 'C' &&
	buffer1[5] == 'O' &&
	buffer1[6] == 'R' &&
	buffer1[7] == 'D') {
      type = COORD_CHARMM_TRAJECTORY;
    }

      /*
       *  Does the second line say "REMD" and therefore represent an AMBER REMD trajectory?
       *  DAN ROE: For reading HREMD trajectories, the only difference format-wise is that
       *           the header line begins HREMD instead of REMD. By treating the lambda
       *           value as if it were a temperature this code should work fine. May want
       *           to change this in the future.
       */
    if (((buffer2[0] == 'R')&&(buffer2[1] == 'E')&&
         (buffer2[2] == 'M')&&(buffer2[3] == 'D')) ||
        ((buffer2[0] == 'H')&&(buffer2[1] == 'R')&&
         (buffer2[2] == 'E')&&(buffer2[3] == 'M')&&
         (buffer2[4] == 'D')) ||
        ((buffer2[0] == 'R')&&(buffer2[1] == 'X')&&
         (buffer2[2] == 'S')&&(buffer2[3] == 'G')&&
         (buffer2[4] == 'L')&&(buffer2[5] == 'D')))
    {
      type = COORD_AMBER_REMD;
    }

      /*
       * DAN ROE - Modified so that ptraj will recognize restarts from
       * REMD, which have 3 numbers in the first line.
       *  Can we scan in four numbers from the first line?  If so this is 
       *  an AMBER trajectory otherwise assume it is an AMBER restrt file
       * tec3: mod to recognize a trajectory with only a single atom for mmpbsa
       */
    if ( type == COORD_UNKNOWN ) {
      if ( (i = sscanf(buffer2, "%f%f%f%f", &junk1, &junk2, &junk3, &junk4)) == 4 )
	type = COORD_AMBER_TRAJECTORY;
      else if ( (i = sscanf(buffer3, "%f%f%f%f", &junk1, &junk2, &junk3, &junk4)) == 3 )
	type = COORD_AMBER_TRAJECTORY;
      else
	type = COORD_AMBER_RESTART;
    }
  }

  /*
   *  EXIT if we have not figured out what kind of file this is
   */

  if ( type == COORD_UNKNOWN ) {
    printfone("  WARNING: trajin %s, unsupported file format...\n",
	    filename);
    safe_free(trajInfo->filename);
    safe_free(trajInfo);
    return NULL;
  }

  /*
   *  REOPEN THE FILE (if the file is not a NetCDF file)
   */

  if ( type != COORD_AMBER_NETCDF ) {
    if (trajFile_rewind(trajInfo)!=0 ) { 
      fprintf(stdout,"Error: checkCoordinates(): Unable to reset trajectory file!\n");
      safe_free(trajInfo->filename);
      safe_free(trajInfo);
      return NULL;
    }
  }

     /*
      *  Read through the file, depending on type, and make sure
      *  that they have the same number of atoms and casually check
      *  to make sure the file isn't corrupted.  In the case of a 
      *  trajectory, check to see how many frames are present in the
      *  file...
      */

#ifdef MPI
  if ( (worldsize > 1) && (type != COORD_AMBER_TRAJECTORY && type != COORD_AMBER_NETCDF) ) {
    if (worldrank == 0)
      fprintf(stderr, "ERROR: multiptraj currently only supports AMBER TRAJECTORY and AMBER NETCDF files.\n");
    return NULL;
  }
#endif

  /* Now that the trajectory has been identified, set up the trajectory
   * information.
   */
  actualAtoms = 0;

  switch ( type ) {

  case COORD_AMBER_NETCDF:
    if (NETCDF_setup(trajInfo,&actualAtoms)!=0) {
      printfone("Error setting up NETCDF trajectory for reading.\n");
      cleanTraj(trajInfo);
      trajInfo=NULL;
      return NULL;
    }
    // DEBUG info
    if (prnlev>0) {
      NETCDF_info_debug(trajInfo->NCInfo,trajInfo->filename); 
      dan_netcdf_debug(trajInfo->NCInfo->ncid);
    }
    // Needed for later MPI call
    // NOTE: shouldnt be needed anymore
    stop=trajInfo->stop;
    isBox=trajInfo->isBox;

    break;

  case COORD_PDB:

       /*
        *  read PDB records until the EOF is encountered counting
        *  the number of atoms encountered along the way...
        *  NOTE: this code will NOT properly handle pdbs that are
        *  strung together into a single file.
        */

    for (;;) {
      r = pdb_read_record(fp);

      if ( r.record_type == PDB_ATOM || r.record_type == PDB_HETATM )
	actualAtoms++;

      if ( r.record_type == PDB_END ) break;
    }

    break;

  case COORD_AMBER_RESTART:


    restrt = readAmberRestart(totalAtoms, filename, fp);
    if (restrt == NULL) return NULL;
    actualAtoms = restrt->natoms;
    isBox = restrt->isbox;
    isVelocity = restrt->restart;
    /* DAN ROE: Free restart structure, not used anymore  */
    FREE_Restart(restrt);

    break;

  case COORD_AMBER_TRAJECTORY:

    /*
     * Note: Recent improvements have been made to the processing of AMBER trajectory files.
             These include the use of frames during reading and the ability to seek in such
             files. This significantly improves the I/O speed when not all frames in a 
             trajectory are to be processed. 

             Additionally support has been added for reading of such trajectory files in
             parallel.

             Performance contributions by: Eric Absgarten (Stony Brook) and Paul Frybarger (SDSC) 
     */
    
      /* Read in title, get size in bytes */
      if (trajFile_fgets(buffer2,BUFFER_SIZE,trajInfo) == NULL) {
	fprintf(stderr, "WARNING in %s: EOF encountered during reading of\n", ROUTINE);
	fprintf(stderr, "   title from (%s)\n", filename);
	return NULL;
      }
      trajInfo->titleSize = strlen(buffer2);
      // Calculate the length of each coordinate frame in bytes
      frame_lines = (totalAtoms * 3) / 10;
      if (((totalAtoms * 3) % 10) > 0)
	frame_lines++;
      trajInfo->frameSize = ((totalAtoms * 3 * 8) + frame_lines);
      if (prnlev>0) fprintf(stdout,"[%i] Each frame has %i chars plus %i newlines. Total %i.\n",
                            worldrank,totalAtoms*3*8,frame_lines,trajInfo->frameSize);
      trajInfo->buffer=(char*) calloc((trajInfo->frameSize), sizeof(char));
      // Read the first frame of coordinates
      if (trajFile_fread(trajInfo)==0) {
        fprintf(stdout,"ERROR in read of Coords frame 1\n");
        return NULL;
      }
      // DEBUG - Print first line of coords
      if (prnlev>0) {
        memcpy(buffer1,trajInfo->buffer,81);
        buffer1[80]='\0';
        fprintf(stdout,"[%i] First line: %s\n",worldrank,buffer1);
      }
      /* Check for box coordinates. If present, update the frame size and 
       * reallocate the frame buffer.
       */
      if (trajFile_fgets(buffer1,BUFFER_SIZE,trajInfo)!=NULL) {
        if (prnlev>0) fprintf(stdout,"DEBUG: Box Line: (%s)\n",buffer1);
        numBoxCoords = sscanf(buffer1, "%f%f%f%f%f%f%f%f%f",
  			      &junk1, &junk2, &junk3, &junk4, &junk5, &junk6,
			      &junk7, &junk8, &junk9);
        if (totalAtoms > 2 && numBoxCoords < 9) {
  	  isBox = 1;
	  trajInfo->frameSize += (numBoxCoords * 8 + 1);
          trajInfo->buffer=(char*) realloc(trajInfo->buffer,(trajInfo->frameSize) * sizeof(char));
        }
      }

      /* Calculate number of frames. If not possible this could be a compressed
       * file and frames will be read until EOF. 
       * NOTE: It is necessary to use the stat command to get the file size 
       * instead of seeking in case the file has been popen'd.
       */
      if (stat(filename, &frame_stat) != 0) {
	fprintf(stderr, "WARNING in %s: Could not find file status for %s\n", ROUTINE, filename);
	return NULL;
      }

      // Determine Uncompressed File Size for Nframes calculation
      file_size=0;
      if (trajInfo->compressType==1)      // Gzip
        file_size=gzipFileSize(filename);
      else if (trajInfo->compressType==2) // Bzip2
        file_size=bzip2FileSize(filename);
      else if (trajInfo->compressType==3) // Zip
        file_size=zipFileSize(filename);
      if (file_size<0) {
        fprintf(stdout,
                "ERROR in %s: Could not calculate uncompressed file size for %s\n",
                ROUTINE, filename);
        return NULL;
      }
      if (file_size==0) file_size=frame_stat.st_size;

      if (prnlev>0)
        fprintf(stdout,"[%s] titleSize=%u frameSize=%u UncompfileSize=%lli FileSize=%lli\n",
                filename,trajInfo->titleSize,trajInfo->frameSize,file_size,frame_stat.st_size);
      frame_size = (long long int) trajInfo->titleSize;
      file_size = file_size - frame_size; // Subtract title size from file total size.
      frame_size = (long long int) trajInfo->frameSize;
      trajInfo->Nframes = (int) (file_size / frame_size);

      // Frame calculation for large gzip files
      // If uncompressed size is less than compressed size, uncompressed
      // size is likely > 4GB.
      if (trajInfo->compressType == 1 && file_size < (long long int)frame_stat.st_size) {
        // Since this is gzip compressed, if the file_size % frame size != 0, 
        // it could be that the uncompressed filesize > 4GB. Since 
        //   ISIZE = uncompressed % 2^32, 
        // try ((file_size + (2^32 * i)) % frame_size) and see if any are 0.
        if ( (file_size % frame_size) != 0) {
          // Determine the maximum number of iterations to try based on the
          // fact that Amber trajectories typically compress about 3x with
          // gzip. If the number of frames cannot accurately be calculated 
          // use the max estimated file size to estimate # frames so that
          // ptraj actions allocate enough memory.
          tmpfsize = (long long int) frame_stat.st_size;
          tmpfsize *= 4;
          tmpfsize = (tmpfsize - file_size) / 4294967296LL;
          maxi = (int) tmpfsize;
          maxi++;
          if (prnlev>0)
            printf("\tLooking for uncompressed gzip size > 4GB, %i iterations.\n",maxi);
          tmpfsize = 0;
          sizeFound=0;
          for (i = 0; i < maxi; i++ ) {
            tmpfsize = (4294967296LL * i) + file_size;
            if ( (tmpfsize % frame_size) == 0) {sizeFound=1; break;}
          }
          if (!sizeFound) {
            printf("Warning: Cannot accurately determine # of frames in gzipped trajectory %s.\n",
                   filename);
            printf("         This usually indicates the trajectory is corrupted.\n");
            printf("         Ptraj will attempt to estimate the correct number of frames.\n");
          }
          file_size = tmpfsize;
          trajInfo->Nframes = (int) (file_size / frame_size);
        }
      }

      if (prnlev>0) fprintf(stdout,"    File has %i frames.\n",trajInfo->Nframes);
      if ( (file_size % frame_size) == 0 ) {
	seekable = 1;
      } else {
	seekable = 0;
	fprintf(stderr, "%s: Could not predict number of frames for AMBER trajectory file: %s\n", ROUTINE, filename);
	fprintf(stderr, "\tIf this is not a compressed file then there is a problem\n");
      }
    stop = trajInfo->Nframes;

    actualAtoms = totalAtoms;
    //trajInfo->frameSize = frame_size;
    //trajInfo->titleSize = title_size;
    trajInfo->numBox = isBox ? numBoxCoords: 0;
    trajInfo->seekable = seekable;

    fprintf(stdout, "Rank: %i Atoms: %i FrameSize: %i TitleSize: %i NumBox: %i Seekable %i\n\n", 
            worldrank, actualAtoms, trajInfo->frameSize, trajInfo->titleSize, 
            trajInfo->numBox, trajInfo->seekable);

    break;

  case COORD_AMBER_REMD:
       /*
        * Just like reading a 'normal' Amber trajectory except there is one
        * extra line (at beginning of each coord set) to get through.
        */
    actualAtoms = totalAtoms;
    if (fgets(buffer1, BUFFER_SIZE, fp) == NULL) {
      printfone("WARNING in %s: EOF encountered during reading of\n", ROUTINE);
      printfone("   title from (%s)\n", filename);
      return NULL;
    }
    lines_per_set = (int) (actualAtoms * 3) / 10;
    if ( (actualAtoms * 3) % 10 ) lines_per_set += 1;
      /*
       *  Add one more for the extra REMD line 
       */
    lines_per_set++;
    for (i=0; fgets(buffer1, BUFFER_SIZE, fp) != NULL; i++) {
      if ( strchr(buffer1, (int) '*') != NULL ) {
	printfone("WARNING in %s, AMBER trajectory file is corrupted: '*' detected (%s)\n",
		  ROUTINE, filename);
        break;
      }
      if (i == lines_per_set) {
        /* DAN ROE: Modified to read HREMD header as well*/
        if (((buffer1[0] == 'R')&&(buffer1[1] == 'E')&&
             (buffer1[2] == 'M')&&(buffer1[3] == 'D')) ||
            ((buffer1[0] == 'H')&&(buffer1[1] == 'R')&&
             (buffer1[2] == 'E')&&(buffer1[3] == 'M')&&
             (buffer1[4] == 'D'))||
            ((buffer1[0] == 'R')&&(buffer1[1] == 'X')&&
             (buffer1[2] == 'S')&&(buffer1[3] == 'G')&&
             (buffer1[4] == 'L')&&(buffer1[5] == 'D')))
          isBox = 0;
        else
          isBox = 1;
      }
    }
    trajInfo->linesperset=lines_per_set+isBox;
    stop = i / (lines_per_set + isBox);
    // Allocate space for REMD TRAJ 0
    trajInfo->REMDtraj=(coordinateInfo**) malloc(sizeof(coordinateInfo*));
    trajInfo->REMDtraj[0]=trajInfo; 
    trajInfo->numREMDTRAJ=1;
    break;

  case COORD_BINPOS:

    if (openbinpos(fp) < 0) break;

    binposScratch = safe_malloc(sizeof(float) * totalAtoms * 3);
    binposeof = 0;
    for (i=0; binposeof == 0; i++) {

      if (readbinpos(fp, &actualAtoms, binposScratch, &binposeof) < 0)
        binposeof = 1;
    }
    stop = i-1;
    safe_free(binposScratch);
	break; 

  case COORD_CHARMM_TRAJECTORY:

    /*
     *  pre-process
     */
    stop = -1;
    x = NULL;
    charmmTrajectoryp = (charmmTrajectoryInfo **) safe_malloc(sizeof(charmmTrajectoryInfo *));
    *charmmTrajectoryp = NULL;
    readCharmmTrajectory(fp, charmmTrajectoryp, x, x, x, x, stop);
    charmmTrajectory = *charmmTrajectoryp;

    x = (double *) safe_malloc(sizeof(double) * charmmTrajectory->natrec);

    /*
     *  read in sets until there are no more
     */
    stop = 0;
    while ( readCharmmTrajectory(fp, charmmTrajectoryp, x, x, x, x, stop+1) )
      stop++;
    actualAtoms = charmmTrajectory->natrec;
    isVelocity = 0;
    isBox = charmmTrajectory->icntrl[10];

    if (charmmTrajectory->icntrl[0] != stop) {
      printfone("NOTE: this charmm trajectory contains %i sets, expecting %i\n",
		stop, charmmTrajectory->icntrl[0]);
    }

    safe_free(x);
    x = NULL;
    break;

  case COORD_UNKNOWN:
    safe_free(trajInfo->filename);
    safe_free(trajInfo);
    return NULL;
  }


  if (actualAtoms != totalAtoms) {
    printfone("  WARNING in %s: The actual number of atoms (%d)\n", ROUTINE, actualAtoms );
    printfone("  does not match the expected number of atoms (%d) in (%s)\n", 
	    totalAtoms, filename);
    printfone("  With this version of the code, this will likely lead to program failure!!!\n");
  }

  // Set number of frames. Amber trajectory sets Nframes above
  if (trajInfo->Nframes==0 && stop>0) trajInfo->Nframes=stop;
  trajInfo->start = start;
  trajInfo->stop = stop;
  trajInfo->offset = 1;
  trajInfo->isBox = isBox;
  trajInfo->isVelocity = isVelocity;
  trajInfo->type = type;
  // If compressed, no seek possible
  if (trajInfo->compressType > 0) trajInfo->seekable=0;

  // This should eventually get its own structure
  if (type==COORD_CHARMM_TRAJECTORY) 
    trajInfo->info = (void *) charmmTrajectory;
  
  if (closeTraj(trajInfo)!=0)
    printfone("WARNING: Error closing traj file %s!\n",trajInfo->filename);

  return ( trajInfo );
}

#undef  ROUTINE
#define ROUTINE "checkCoordinatesWrap()"

   void
checkCoordinatesWrap(char *filename)
{
  coordinateInfo *info;
  charmmTrajectoryInfo *charmmTraj;
  byte u;

  info = checkCoordinates(filename, parm->NTOTAT);

  if (info == NULL) {
    printfone("WARNING in %s.  Encountered a problem analyzing\n", ROUTINE);
    printfone("coordinates from file %s\n", filename);
    return;
  }

  switch ( info->type ) {

  case COORD_AMBER_NETCDF:
    printfone("File (%s) is a NetCDF AMBER trajectory%s",
	    info->filename, (info->isBox ? " with box coordinates" : ""));
    if (info->isVelocity > 0)
      printfone(" with velocities");
    if (info->stop > 0)
      printfone(" representing %i sets\n", info->stop);
    else
      printfone("\n");
    break;

  case COORD_PDB:
    printfone("File (%s) is a PDB file\n", info->filename);
    break;

  case COORD_AMBER_TRAJECTORY:
    
      printfone("File (%s) is an AMBER trajectory%s",
	      info->filename, (info->isBox ? " with box coordinates" : ""));
      if (info->stop > 0)
	printfone(" representing %i sets\n", info->stop);
      else
	printfone("\n");
      break;

  case COORD_AMBER_REMD:
      /*
       *  Not complete!
       */
    printfone("File (%s) is an AMBER REMD (new format) trajectory%s",
            info->filename, (info->isBox ? " with box coordinates" : ""));
    if (info->stop > 0)
      printfone(" representing %i sets\n", info->stop);
    else
      printfone("\n");
    break;

  case COORD_CHARMM_TRAJECTORY:

    charmmTraj = (charmmTrajectoryInfo *) info->info;
    fprintf(stdout,
	    "File (%s) is a CHARMM trajectory in %s endian binary format %s",
	    info->filename, (charmmTraj->byteorder ? "little" : "big"),
	    (info->isBox ? "with box coordinates" : ""));
    if (info->stop > 0)
      printfone("representing %i sets\n", info->stop);
    else
      printfone("\n");
    if (prnlev > 2) {
      printf("  NFILE = %i\n", charmmTraj->icntrl[0]); /* number of coordinate sets in file   */
      printf("  ISTEP = %i\n", charmmTraj->icntrl[1]); /* number of previous dynamics steps   */
      printf("  NINTV = %i\n", charmmTraj->icntrl[2]); /* frequency for saving coordinates    */
      printf("  NSAVC = %i\n", charmmTraj->icntrl[3]); /* number of steps for creation runs   */
      printf("  NSAVV = %i\n", charmmTraj->icntrl[4]);
      printf("  NDEGF = %i\n", charmmTraj->icntrl[7]);
      printf("  NFREA = %i\n", charmmTraj->icntrl[8]);
      u.i = charmmTraj->icntrl[9];
      printf("  DELTA = %f\n", u.f);
      printf("  QCRYS = %i\n", charmmTraj->icntrl[10]);
      printf("  QDIM4 = %i\n", charmmTraj->icntrl[11]);
      printf("  VERNU = %i\n", charmmTraj->icntrl[19]);
      printf("  Dumping the title:\n");
      printStack(&charmmTraj->titleStack, printString, NULL);
    }
    break;

  case COORD_AMBER_RESTART:
    printfone("File (%s) is an AMBER restart file ", info->filename);
    if ( info->isBox && info->isVelocity )
      printfone("with box and velocity information\n");
    else if (info->isBox)
      printfone("with box information\n");
    else if (info->isVelocity)
      printfone("with velocity information\n");
    else
      printfone("\n");
    break;
  }
}



#undef  ROUTINE
#define ROUTINE "loadCharmmPSF()"

   ptrajState *
loadCharmmPSF(FILE *fd, int skipheader)
{
  char *bufferp, *buffer, *pbuffer;
  int natom;

  double *charges, *masses;
  float charge, mass;
  int  *tmp_ipres, *tmp_molinfo;
  Name *atomName;
  Name residueNumber;
  Name *tmp_resName;
  Name segid, cursegid, vdwt;
  int atom;
  int residue;
  int fixed;
  int vdwp;
  int i,j,k,len, scannedValues;
  int curres, curresidx, curmol, numcurmol;
  int solventMolecules, *solventMoleculeStart, *solventMoleculeStop;
  int *solventMask, solventAtoms;
  int *tmp_solventMoleculeStart, *tmp_solventMoleculeStop;
  int ntitle;
  ptrajState *state;

  buffer = (char *) safe_malloc(sizeof(char) * BUFFER_SIZE);
  pbuffer = (char *) safe_malloc(sizeof(char) * BUFFER_SIZE);

     /*
      *  Open up the PSF file
      */
  if (fd == NULL) {
    error(ROUTINE, "Couldn't open file psf\n");
  }

     /*
      *  grab header and make sure it matches "PSF"
      */
  if (skipheader != 1) {
    bufferp = fgets(buffer, BUFFER_SIZE, fd);
    if (bufferp == NULL || strncmp(buffer, "PSF", 3) != 0) {
      error(ROUTINE, "File is not a CHARMM PSF file!\n");
    }
  }
  printfone("Reading in CHARMM PSF file\n");

     /*
      *  skip the next line, read the title and the blank line that follows
      */
  bufferp = fgets(buffer, BUFFER_SIZE, fd);
  bufferp = fgets(buffer, BUFFER_SIZE, fd);
  if (strstr(bufferp, "!NTITLE")) {
    if ( sscanf(bufferp, "%i", &ntitle) != 1 ) 
      ntitle = 0;
  }


  printfone("Reading in the title...\n\n");
  if (ntitle) {
    for (i=0; i<=ntitle; i++) {
      bufferp = fgets(buffer, BUFFER_SIZE, fd);
      printfone("%s", buffer);
    }

  } else {

    bufferp = fgets(buffer, BUFFER_SIZE, fd);
    while (buffer[0] == '*')
      bufferp = fgets(buffer, BUFFER_SIZE, fd);
      printfone("%s", buffer);
  }

     /*
      *  get the total number of atoms
      */
  bufferp = fgets(buffer, BUFFER_SIZE, fd);

  if (sscanf(buffer, "%i", &natom) != 1) error(ROUTINE, "Scanning natoms\n");
  printfone("Total number is atoms is %d\n", natom);


  /*
   *  read in the atom segid, residue, name and charge/mass info
   */
  charges     = (double *) safe_malloc(sizeof(double) * natom);
  masses      = (double *) safe_malloc(sizeof(double) * natom);
  tmp_ipres   = (int *)    safe_malloc(sizeof(double) * natom);
  tmp_molinfo = (int *)    safe_malloc(sizeof(double) * natom);

  atomName    = (Name *)   safe_malloc(sizeof(Name) * natom);
  tmp_resName = (Name *)   safe_malloc(sizeof(Name) * natom);

  tmp_solventMoleculeStart = (int *) safe_malloc(sizeof(int) * natom);
  tmp_solventMoleculeStop  = (int *) safe_malloc(sizeof(int) * natom);
  solventMask = (int *) safe_malloc(sizeof(int) * natom);
  memset(solventMask, 0, sizeof(int) * natom);
  memset(tmp_solventMoleculeStart, 0, sizeof(int) * natom);
  memset(tmp_solventMoleculeStop,  0, sizeof(int) * natom);

  curres = 1;
  curresidx = 1;
  tmp_ipres[0] = 1;
  solventMolecules = 0;
  solventAtoms = 0;

  printfone("Reading in the atom information...\n");

  for (i=0; i < natom; i++) {

    bufferp = fgets(buffer, BUFFER_SIZE, fd);

       /*
        * atom  seg res# res atom  vdwp  charge    mass    fixed?
        *  1    WAT  1   WAT  OH    2   -0.820000 15.9994    0
        */
    scannedValues = sscanf(buffer, "%8i %4s %4i %4s %4s %4i %f %f%8i",
			   &atom, (char *) &segid, &residue, (char *) &tmp_resName[i], 
			   (char *) &atomName[i], &vdwp, &charge, &mass, &fixed);

    if (scannedValues < 9) {

      /*
       *  assume this is an XPLOR psf file w/ text for the vdw type...
       *  NOTE: there was a funny problem with reading residues as numbers
       *  since when they get over 9999 they revert to 0001 and the sscanf choked.
       *  Replace the number with a string and then convert and it seems to work.
       */

      scannedValues = sscanf(buffer, "%8i %4s %4s %4s %4s %4s %f %f %i",
			     &atom, (char *) &segid, residueNumber, (char *) &tmp_resName[i], 
			     (char *) &atomName[i], (char *) &vdwt, &charge, &mass, &fixed);
      if (scannedValues == 9) {	
	residue = atoi(residueNumber);
      } else {

	/*
	 *  assume this is an XPLOR psf file w/o SEGID information...
	 */
	if (sscanf(buffer, "%8i %8i %4s %4s %4s %f %f%8i",
		   &atom, &residue, 
		   (char *) &tmp_resName[i], (char *) &atomName[i], (char *) &vdwt,
		   &charge, &mass, &fixed) < 8) {
	  
	  error(ROUTINE, "Scanning atom information\n");
	}
      }
    }

    atomName[i][NAME_SIZE-1] = (char) 0;
    tmp_resName[i][NAME_SIZE-1] = (char) 0;
    for (j=NAME_SIZE-2; atomName[i][j] == (char) 0; j--) {
      atomName[i][j] = ' ';
    }
    for (j=NAME_SIZE-2; tmp_resName[i][j] == (char) 0; j--) {
      tmp_resName[i][j] = ' ';
    }

    if (i==0) {
      strcpy(cursegid, segid);
      curmol = 0;
      numcurmol = 1;
    } else {

      if (strcmp(cursegid, segid) == 0) {
	numcurmol++;
      } else {
	tmp_molinfo[curmol] = numcurmol;
	curmol++;
	numcurmol = 1;
	strcpy(cursegid, segid);
      }
    }


    charges[i] = charge;
    masses[i] = mass;

    if (prnlev > 1) 
      printf("Atom %4i (%s) seg (%4s) res %4i (%s) charge: %7.3f mass: %7.3f fixed? %s %i\n",
	     atom, atomName[i], segid, residue, tmp_resName[i], charge, mass,
	     (fixed ? "yes" : "no"), vdwp);


    if (residue != curres) {
      tmp_ipres[curresidx] = i+1;
      curresidx++;
      curres = residue;
    }
  }

  tmp_ipres[curresidx] = i+1;
  tmp_molinfo[curmol] = numcurmol;

  for (i=0; i < curresidx; i++) {
    j = tmp_ipres[i]-1;
    if ( strcmp(tmp_resName[j], "WAT ") == 0 || 
	 strcmp(tmp_resName[j], "IP3 ") == 0 ||
	 strcmp(tmp_resName[j], "TIP3") == 0) {
      tmp_solventMoleculeStart[solventMolecules] = j;
      tmp_solventMoleculeStop[solventMolecules] = j+3;
      solventMolecules++;
    }
  }

  solventMoleculeStart = (int *) safe_malloc(sizeof(int) * (solventMolecules+1));
  solventMoleculeStop  = (int *) safe_malloc(sizeof(int) * (solventMolecules+1));
  for (i = 0; i < solventMolecules; i++) {

    solventMoleculeStart[i] = tmp_solventMoleculeStart[i];
    solventMoleculeStop[i]  = tmp_solventMoleculeStop[i];

    for (j=solventMoleculeStart[i]; j < solventMoleculeStop[i]; j++) {
      solventMask[j] = 1;
      solventAtoms++;
    }
  }

  if (prnlev > 1) {
    printf("The total number of molecules is %i\n", curmol+1);
    printf("The total number of residues is %i\n", curresidx);
  }

  printfone("Dumping out residue names:\n");
  j = 0;
  k = 0;
  len = 0;
  pbuffer[0] = (char) 0;
  for (i=0; i < curresidx; i++) {

    k = tmp_ipres[i]-1;
    sprintf(buffer+j, "%5s", tmp_resName[k]);
    j += 5;
    if ( i == curresidx - 1 || ( i != 0 && (i+1) % 10 == 0 ) ) {
      buffer[j] = (char) 0;
      if ( strcmp(pbuffer, buffer) == 0 ) {
	if ( ! len ) {
	  printfone(" ...\n");
	  len = 1;
	}
      } else {
	printfone("%s\n", buffer);
	len = 0;
	strcpy(pbuffer, buffer);
      }
      j = 0;
    }
  }

  state = (ptrajState *) safe_malloc(sizeof(ptrajState));
  INITIALIZE_ptrajState(state);

  state->box[0] = 0.0;
  state->box[1] = 0.0;
  state->box[2] = 0.0;
  state->box[3] = 90.0;
  state->box[4] = 90.0;
  state->box[5] = 90.0;
  state->atoms = natom;
  state->residues = curresidx;
  state->charges = charges;
  state->masses = masses;
  state->ipres = (int *) safe_malloc(sizeof(int) * (state->residues+1));
  state->residueName = (Name *) safe_malloc(sizeof(Name) * state->residues);
  for (i=0; i <= state->residues; i++)
    state->ipres[i] = tmp_ipres[i];
  for (i=0; i < state->atoms; i++) {
    curresidx = atomToResidue(i+1, state->residues, state->ipres)-1;
    strcpy(state->residueName[curresidx], tmp_resName[i]);
    i = state->ipres[curresidx+1]-2;
  }

  state->IFBOX = 0;
  state->molecules = curmol+1;
  state->moleculeInfo = (int *) safe_malloc(sizeof(int) * state->molecules);
  for (i=0; i < state->molecules; i++) {
    state->moleculeInfo[i] = tmp_molinfo[i];
  }

  state->solventMask = solventMask;
  state->solventAtoms = solventAtoms;
  state->solventMolecules = solventMolecules;
  state->solventMoleculeStart = solventMoleculeStart;
  state->solventMoleculeStop = solventMoleculeStop;

  state->atomName = atomName;

  safe_free(buffer);
  safe_free(tmp_ipres);
  safe_free(tmp_molinfo);
  safe_free(tmp_resName);
  safe_free(tmp_solventMoleculeStart);
  safe_free(tmp_solventMoleculeStop);

  if (prnlev >= 0) 
    ptrajPrintState(state);
  return state;
  
}



#undef  ROUTINE
#define ROUTINE "ptrajInitializeState()"

   void
ptrajInitializeState(ptrajState **statep, char *filename)
{
  FILE *fp;
  char *buffer;
  ptrajState *state;
  pdb_record pdbr, *pdb;
  double *X, *Y, *Z;
  int i, j, k, pdbRecords, *start, *stop;
  int ires, curres;

  /*
   *  open up the filename (if it exists, otherwise prompt the user)
   *  and determine whether it is a AMBER prmtop or CHARMM PSF file
   *  or a PDB file...
   */

  if (filename == NULL) {
    filename = promptToOpenFile(&fp, "", "r", 
             "Input the name of an AMBER prmtop, CHARMM PSF or PDB file: ");
  } else {
    if ( openFile(&fp, filename, "r") == 0 )
      error(ROUTINE, "Attempting to open parameter/topology file %s",
	    filename);
  }
  buffer = (char *) safe_malloc(sizeof(char)* BUFFER_SIZE);
  if (fgets(buffer, BUFFER_SIZE, fp) == NULL)
    error(ROUTINE, "Attempting to read parameter/topology file");

  /*
   *  check to see if the file has "PSF" or "PSF " as the header
   *  Examples of CHARMM file headers are:
   *  PSF
   *  PSF CMAP
   *  PSF CMAP CHEQ
   *  indicating a CHARMM PSF file, load it and return...
   */
  /* Look for "PSF\n" or "PSF " - we do this separately since
   * It can help with false positives if we just searched for
   * PSF for example.
   */
  if ( strcmp(buffer, "PSF\n") == 0 ||
       strcmp(buffer, "PSF \n") == 0 ||
       strncmp(buffer, "PSF ", 4) == 0 ) {
    state = loadCharmmPSF(fp, 1);
    *statep = state;
    /* DAN ROE: Make sure buffer is freed  */
    safe_free(buffer);
    return;
  }
  safe_free(buffer);

  /*
   *  check to see if the file is a valid PDB record
   */
  /*safe_freopen(fp);  */
  /* DAN ROE: Why is the file constantly reopened? Should probably just be 
   * rewound.  safe_freopen does a lot of allocating/deallocating which gets 
   * messy and leads to memory holes.
   */
  rewind(fp);

  pdbr = pdb_read_record(fp);
  if ( pdbr.record_type != PDB_UNKNOWN ) {

    /*safe_freopen(fp);   */
    rewind(fp);
    pdbRecords = loadPdb(fp, &pdb);
    state = (ptrajState *) safe_malloc( sizeof(ptrajState) );
    INITIALIZE_ptrajState(state);
    X = NULL;
    Y = NULL;
    Z = NULL;
    state->atoms = getCoordinatesFromPdb(pdb, X, Y, Z);
    printfone("Reading PDB file \"%s\": %i atoms found!!!\n", filename, state->atoms);

    state->atomName = (Name *) safe_malloc( sizeof(Name) * state->atoms);
    state->masses   = (double *) safe_malloc(sizeof(double) * state->atoms);
    state->charges  = (double *) safe_malloc(sizeof(double) * state->atoms);

    i = 0;
    ires = 0;
    curres = -1;
    for (j=0; j < pdbRecords; j++) {
      if ( pdb[j].record_type == PDB_ATOM ||
	   pdb[j].record_type == PDB_HETATM ) {

	state->atomName[i][0] = pdb[j].pdb.atom.name[1];
	if (pdb[j].pdb.atom.name[2] == (char) 0)
	  state->atomName[i][1] = ' ';
	else
	  state->atomName[i][1] = pdb[j].pdb.atom.name[2];

	if (pdb[j].pdb.atom.name[3] == (char) 0)
	  state->atomName[i][2] = ' ';
	else
	  state->atomName[i][2] = pdb[j].pdb.atom.name[3];

	state->atomName[i][3] = pdb[j].pdb.atom.name[0];
	state->atomName[i][4] = (char) 0;

	/*
	printf("Loaded up atom named -%s-\n", state->atomName[i]);
	*/
	state->masses[i] = 0.0;
	state->charges[i] = 0.0;
	if ( pdb[j].pdb.atom.residue.seq_num != curres ) {
	  ires++;
	  curres = pdb[j].pdb.atom.residue.seq_num;
	}
	i++;
      }

    }
    printf("Found %i residues\n", ires);

    state->residues = ires;
    state->residueName = (Name *) safe_malloc (sizeof(Name) * state->residues);
    state->ipres = (int *) safe_malloc(sizeof(int) * (state->residues+1));

    i = 0;
    ires = -1;
    curres = -1;
    for (j=0; j < pdbRecords; j++) {
      if ( pdb[j].record_type == PDB_ATOM ||
	   pdb[j].record_type == PDB_HETATM ) {
	i++;
	if ( pdb[j].pdb.atom.residue.seq_num != curres ) {
	  ires++;
	  curres = pdb[j].pdb.atom.residue.seq_num;
	  strcpy(state->residueName[ires], pdb[j].pdb.atom.residue.name);
	  if (ires >= 0) {
	    state->ipres[ires] = i;
	  }
	}
      }
    }
    state->ipres[state->residues] = state->atoms+1;

    state->molecules = 0;
    state->solventMolecules = 0;
    state->box[0] = 0.0;
    state->box[1] = 0.0;
    state->box[2] = 0.0;
    state->box[3] = 90.0;
    state->box[4] = 90.0;
    state->box[5] = 90.0;
    state->maxFrames = 0;

    /*

  for (i=0; i <= state->residues; i++) {

    printf("RESIDUE %i, ipres: %i\n", i+1, state->ipres[i]);

  }
    */


    *statep = state;
    return;
  }

  /*
   *  assume the file points to an AMBER prmtop and load it...
   */

  /*safe_freopen(fp); */
  rewind(fp);

  clearParm( parm );
  parm = safe_malloc( sizeof( Parm ) );
  initializeParm(parm);
  /* DAN ROE: Changed so that filename is copied to its own memory location */
  parm->filename=(char*) safe_malloc( (strlen(filename)+1)*sizeof(char));
  /*parm->filename = filename;  */
  strcpy(parm->filename,filename);
  parm->fp = fp;
  readParm();

  state = (ptrajState *) safe_malloc( sizeof(ptrajState) );
  INITIALIZE_ptrajState(state);

  state->atoms = parm->NTOTAT;
  state->atomName = (Name *) safe_malloc( sizeof(Name) * state->atoms);
  state->masses = (double *) safe_malloc(sizeof(double) * state->atoms);
  state->charges = (double *) safe_malloc(sizeof(double) * state->atoms);
  for (i=0; i < state->atoms; i++) {
    strcpy(state->atomName[i], parm->atom[i].igraph);
    state->masses[i] = parm->atom[i].amass;
    state->charges[i] = parm->atom[i].chrg / CHARGE_TO_KCALS;
  }

  state->residues = parm->NTOTRS;
  state->residueName = (Name *) safe_malloc (sizeof(Name) * state->residues);
  for (i=0; i < state->residues; i++) {
    strcpy(state->residueName[i], parm->residue[i].labres);
  }
  state->ipres = (int *) safe_malloc(sizeof(int) * (state->residues+1));
  for (i=0; i <= state->residues; i++) {
    state->ipres[i] = parm->residue[i].ipres;
  }

  state->box[0] = 0.0;
  state->box[1] = 0.0;
  state->box[2] = 0.0;
  state->box[3] = 90.0;
  state->box[4] = 90.0;
  state->box[5] = 90.0;
  state->maxFrames = 0;

  state->IFBOX = parm->IFBOX;
  if ( state->IFBOX ) {

    state->molecules = parm->box->nspm;
    state->moleculeInfo = (int *)
      safe_malloc( sizeof(int) * state->molecules );
    for (i=0; i < state->molecules; i++) {
      state->moleculeInfo[i] = parm->box->nsp[i];
    }

    state->box[0] = parm->box->box[0];
    state->box[1] = parm->box->box[1];
    state->box[2] = parm->box->box[2];
    if (parm->box->beta != 90.0) {

      if (parm->box->beta > 109.47 && parm->box->beta < 109.48) {

	state->box[3] = 2.0*acos(1.0/sqrt(3.0))*RADDEG;
	state->box[4] = state->box[3];
	state->box[5] = state->box[3];
	printfone(" Setting box to be an exact truncated octahedron, angle is %f\n",
		state->box[3]);

      } else if (parm->box->beta == 60.0) {

	printfone(" Setting box to be a rhombic dodecahedron, i.e. alpha=gamma=60.0, beta=90.0\n");
	state->box[3] = 60.0;
	state->box[4] = 90.0;
	state->box[5] = 60.0;
      }
      
    }
  }

  state->solventMolecules = 0;
  state->solventAtoms = 0;
  state->solventMask = NULL;
  state->solventMoleculeStart = NULL;
  state->solventMoleculeStop = NULL;


  start = (int *) safe_malloc(sizeof(int) * state->atoms);
  stop  = (int *) safe_malloc(sizeof(int) * state->atoms);
  state->solventMask = (int *) safe_malloc(sizeof(int) * state->atoms);
  for (i=0; i < state->atoms; i++)
    state->solventMask[i] = 0;

  if ( ! parm->IFBOX) {
    /*
     *  treat all the molecules starting with parm->box->nspsol
     *  as solvent molecules IF there is box information
     */

    j = 0;
    for (i=0; i < state->molecules; i++) {

      if (i+1 >= parm->box->nspsol) {
	/*
	 *  add this molecule to the solvent list
	 */

	state->solventAtoms += state->moleculeInfo[i];

	for (k = j; k < j+state->moleculeInfo[i]; k++)
	  state->solventMask[k] = 1;
	start[state->solventMolecules] = j;
	stop[ state->solventMolecules] = j+state->moleculeInfo[i];
	state->solventMolecules++;
      }

      j += state->moleculeInfo[i];
    }


  } else {
    /*
     *  treat all the residues named "WAT " as solvent molecules
     */

    for (i=0; i < state->residues; i++) {
      if (strcmp("WAT ", state->residueName[i]) == 0) {
	/*
	 * add this residue to the list of solvent
	 */
	j = state->ipres[i+1]-state->ipres[i];

	state->solventAtoms += j;
	start[state->solventMolecules] = state->ipres[i]-1;
	stop[ state->solventMolecules] = state->ipres[i+1]-1;
	state->solventMolecules++;

	for (k=state->ipres[i]-1; k < state->ipres[i+1]-1; k++)
	  state->solventMask[k] = 1;
      }
    }
  }

  state->solventMoleculeStart = start;
  state->solventMoleculeStop  = stop;


  /*
  for (i=0; i <= state->residues; i++) {

    printf("RESIDUE %i, ipres: %i\n", i+1, state->ipres[i]);

  }
  */

  *statep = state;
  return;
}




#undef  ROUTINE
#define ROUTINE printCoordinateInfo

   void
printCoordinateInfo( void *entry )
{
  coordinateInfo *p;
  charmmTrajectoryInfo *charmmTraj;
  //netcdfTrajectoryInfo *NCInfo;
  byte u;
  int i;

  p = (coordinateInfo *) entry;
  if (p == NULL) {
    printfone("  *** No coordinate file of this type has been defined\n");
    return;
  }

  switch ( p->type ) {

  case COORD_BINPOS:
    printfone("  File (%s) is a BINPOS file ", p->filename);
    if (p->stop > 0) {
      i = (p->stop - p->start)/p->offset+1;
      if (i != p->stop)
	printfone(" with %i sets (processing only %i)\n", p->stop, i);
      else
	printfone(" with %i sets\n", i);
    } else {
      printfone("\n");
    }
    break;

  case COORD_PDB:
    printfone("  File (%s) is a PDB file", p->filename);
    if (p->option2 == 1) 
      printfone(" with no atom wrapping");
    if (p->option1 == 1)
      printfone(": AMBER charges and radii in prmtop to occupancy and temp factor columns");
    else if (p->option1 == 2)
      printfone(": AMBER charges and PARSE radii to occupancy and temp factor columns");
    else if (p->option1 == 3)
      printfone(": AMBER charges and vdw radii (r*) to occupancy and temp factor columns");

    if (p->append) {
      printfone(" appended\n");
    } else
      printfone("\n");

    break;

  case COORD_AMBER_TRAJECTORY:

      printfone("  File (%s) is an AMBER trajectory%s",
	      p->filename, (p->isBox ? " (with box info)" : ""));
      if (p->stop > 0) {
	i = (p->stop - p->start)/p->offset+1;
	if (i != p->stop)
	  printfone(" with %i sets (processing only %i)", p->stop, i);
	else
	  printfone(" with %i sets", i);
      
      if (p->append) {
	printfone(" appended\n");
      } else
	printfone("\n");
    }

    break;

  case COORD_AMBER_REMD:
    printfone("  File (%s) is an AMBER REMD (new format) trajectory%s",
            p->filename, (p->isBox ? " (with box info)" : ""));
    if (p->stop > 0) {
      i = (p->stop - p->start)/p->offset+1;
      if (i != p->stop)
        printfone(" with %i sets (processing only %i)\n", p->stop, i);
      else
        printfone(" with %i sets\n", i);
    } else {
      printfone("\n");
    }
    break;

  case COORD_AMBER_NETCDF:
    
    printfone("  File (%s) is a NetCDF AMBER trajectory%s%s",
	    p->filename, 
	    (p->isBox ? " with box info" : ""),
	    (p->isVelocity ? " and velocities" : ""));

    if (p->NCInfo!=NULL) {
      //NCInfo=(netcdfTrajectoryInfo*) p->info;
      if (p->NCInfo->TempVarID!=-1) 
        printfone(" with replica temperatures");
    }

    if (p->stop > 0) {
      i = (p->stop - p->start)/p->offset+1;
      if (i != p->stop)
	printfone(" with %i sets (processing only %i)\n", p->stop, i);
      else
	printfone(" with %i sets\n", i);
    } else {
      printfone("\n");
    }
    if (prnlev > 2) {
      if (p->title != NULL)
	printfone("    title:        \"%s\"\n", p->title);
      if (p->application != NULL) 
	printfone("    application:  \"%s\"\n", p->application);
      if (p->program != NULL) 
	printfone("    program:      \"%s\"\n", p->program);
      if (p->version != NULL) 
	printfone("    version:      \"%s\"\n", p->version);
    }

    break;

  case COORD_CHARMM_TRAJECTORY:

    charmmTraj = (charmmTrajectoryInfo *) p->info;
    fprintf(stdout,
	    "  File (%s) is a CHARMM trajectory in %s endian binary format\n",
	    p->filename, (charmmTraj->byteorder ? "little" : "big"));
    printfone("%s", (p->isBox ? "    with box information " : "    "));
    if (p->stop > 0) {
      i = (p->stop - p->start)/p->offset+1;
      if (i != p->stop)
	printfone("representing %i sets (processing only %i)\n", p->stop, i);
      else
	printfone("representing %i sets\n", p->stop);
    } else
      printfone("\n");

    /*
     *  NFILE  -- number of coordinate sets in file
     *  ISTEP  -- number of previous dynamics steps
     *  NINTV  -- frequency for saving coordinates
     *  NSAVC  -- number of steps for creation run
     */  
    printfone("  NFILE = %8i ISTEP = %8i NINTV = %8i NSAVC = %8i NSAVV = %8i\n", 
	    charmmTraj->icntrl[0], charmmTraj->icntrl[1],
	    charmmTraj->icntrl[2], charmmTraj->icntrl[3],
	    charmmTraj->icntrl[4]);

    u.i = charmmTraj->icntrl[9];
    printf("  NDEGF = %8i NFREA = %8i DELTA = %8.4f QCRYS = %8i QDIM4 = %8i\n", 
	   charmmTraj->icntrl[7], charmmTraj->icntrl[8], u.f,
	   charmmTraj->icntrl[10], charmmTraj->icntrl[11]);
    printf("  VERNU = %8i NATREC= %8i NFREAT = %8i\n", 
	   charmmTraj->icntrl[19], charmmTraj->natrec, charmmTraj->nfreat);
    printf("  Dumping the title:\n");
    printStack(&charmmTraj->titleStack, printString, NULL);

    break;

  case COORD_AMBER_RESTART:
    printfone("  File (%s) is an AMBER restart file ", p->filename);
    if ( p->isBox && p->isVelocity )
      printfone("with box and velocity information\n");
    else if (p->isBox)
      printfone("with box information\n");
    else if (p->isVelocity)
      printfone("with velocity information\n");
    else
      printfone("\n");
    break;
  }

  /* REMDTRAJ information */
  if (p->isREMDTRAJ>0) {
    if (p->stop>0) {
      fprintf(stdout,"  Replica processing by temperature will occur.\n");
      fprintf(stdout,"  %i files total (First index is %0*i), ",p->numREMDTRAJ,
            p->EXTwidth, p->firstREMDTRAJ);
      fprintf(stdout,"frames at %lf K will be used.\n",p->remdtrajtemp);
    } 
  } /*else
      fprintf(stdout,"only frames from this file will be used (remdtraj not specified), ");*/

  if (prnlev > 4) {
    if ( p->file == NULL ) {
      printfone("    [the FILE is currently closed]\n");
    } else if (p->file == stdin) {
      printfone("    [the FILE is standard input]\n");
    } else if (p->file == stdout) {
      printfone("    [the FILE is standard output]\n");
    } else if (p->file == stderr) {
      printfone("    [the FILE is standard error]\n");
    }
  }

  if ( p->mask != NULL )
    printfone("    [The file has an active atom mask]\n");

}


#undef  ROUTINE
#define ROUTINE "ptrajCleanup()"

   void
ptrajCleanup()
{
  actionInformation *a;
  analyzeInformation *an;
  coordinateInfo *f;
  int *mask;

  /*
   *  Clean up the INPUT files, transformFileStack
   */

  while (transformFileStack != NULL  &&
	 (f = (coordinateInfo *) popStack(&transformFileStack)) != NULL ) {
   
    cleanTraj(f); 
/*    safe_free(f->filename);
    safe_free(f->info);
    safe_free(f->NCInfo); // DEBUG
    safe_free(f->mask);
    safe_free(f->x);
    safe_free(f->y);
    safe_free(f->z);
  
    safe_free(f->title);
    safe_free(f->program);
    safe_free(f->application);
    safe_free(f->version);

    // REMDTRAJ Cleanup 
    safe_free(f->remdtrajfiles);
    safe_free(f->remdncids);
    safe_free(f->baseFilename);
    safe_free(f->compressEXT);

    INITIALIZE_coordinateInfo(f);
    safe_free(f);*/
    
  }
  transformFileStack = NULL;

  /*
   *  Clean up the ACTIONS
   */

  while (transformActionStack != NULL &&
	 (a = (actionInformation *) popStack(&transformActionStack)) != NULL) {

       /*
        *  Free any associated mask
        */

    if (a->mask) 
      safe_free(a->mask);

       /*
        *  Free any associated state information
        */

    if (a->state != NULL) {
      ptrajClearState(&a->state);
      a->state = NULL;
    }

       /*
        *  Clean up any of the complex arguments as necessary; this
        *  is done by the associated action function in the PTRAJ_CLEANUP
        *  mode
        */

    if (a->type != TRANSFORM_TRANSFORM &&
	a->type != TRANSFORM_NOOP) {

      if (a->fxn != NULL) {
	a->fxn(a, NULL, NULL, NULL, NULL, PTRAJ_CLEANUP);
      }
      INITIALIZE_actionInformation(a);
    }

       /*
        *  Free up the action
        */

    safe_free(a);

  }

  transformActionStack = NULL;

  /*
   *  Clean up the ANALYZEs
   */

  while (transformAnalyzeStack != NULL &&
	 (an = (analyzeInformation *) popStack(&transformAnalyzeStack)) != NULL) {

       /*
        *  Clean up any of the complex arguments as necessary; this
        *  is done by the associated action function in the PTRAJ_CLEANUP
        *  mode
        */

    if (an->fxn != NULL) {
      an->fxn(an, NULL, PTRAJ_CLEANUP);
    }
    INITIALIZE_analyzeInformation(an);

       /*
        *  Free up the action
        */
    safe_free(an);

  }

  transformAnalyzeStack = NULL;

  /*
   *  Free up outInfo
   */

  if (globalOutInfo != NULL) {
    if (globalOutInfo->buffer!=NULL)
      free(globalOutInfo->buffer);
    safe_free(globalOutInfo->filename);
    globalOutInfo->filename = NULL;
    safe_free(globalOutInfo->mask);
    globalOutInfo->mask = NULL;
    safe_free(globalOutInfo->info);
    globalOutInfo->info = NULL;
    safe_free(globalOutInfo->NCInfo); // DEBUG
    globalOutInfo->NCInfo=NULL;       // DEBUG
    safe_free(globalOutInfo);
    globalOutInfo = NULL;
  }

  /*
   *  Free up referenceInfo
   */

  while (transformReferenceStack != NULL  &&
	 (f = (coordinateInfo *) popStack(&transformReferenceStack)) != NULL ) {
    
    safe_free(f->filename);
    safe_free(f->info);
    safe_free(f->NCInfo); // DEBUG
    safe_free(f->mask);
    safe_free(f->x);
    safe_free(f->y);
    safe_free(f->z);
    INITIALIZE_coordinateInfo(f);
    safe_free(f);
    
  }
  transformReferenceStack = NULL;
  referenceInfo = NULL;

  /*
   *  Clean up hbond donor/acceptor information
   */

  while (hbondDonorStack != NULL &&
	 (mask = (int *) popStack(&hbondDonorStack)) != NULL) {
    safe_free(mask);
  }
  while (hbondAcceptorStack != NULL &&
	 (mask = (int *) popStack(&hbondAcceptorStack)) != NULL) {
    safe_free(mask);
  }
  while (hbondAcceptorH1Stack != NULL &&
	 (mask = (int *) popStack(&hbondAcceptorH1Stack)) != NULL) {
    safe_free(mask);
  }
  while (hbondAcceptorH2Stack != NULL &&
	 (mask = (int *) popStack(&hbondAcceptorH2Stack)) != NULL) {
    safe_free(mask);
  }
  while (hbondAcceptorH3Stack != NULL &&
	 (mask = (int *) popStack(&hbondAcceptorH3Stack)) != NULL) {
    safe_free(mask);
  }

  /* DAN ROE: Close any files left on fileStack and free parm struct  */
  // NOTE: Is the fileStack obsolete now? Can we ditch it?
  fileStack_clear(); 
  clearParm( parm );

}


#undef  ROUTINE
#define ROUTINE "ptrajPreprocessInputCoordinates()"

   int
ptrajPreprocessInputCoordinates(coordinateInfo *currentCoordinateInfo)
{
  charmmTrajectoryInfo **charmmTrajp;
  char buffer[BUFFER_SIZE];
  int err;
  int i,j;

  /*
   *  open up the file. Exit if any errors occur.
   */
  err=openTraj(currentCoordinateInfo);
  if (err!=0) {
    printfone("WARNING in ptrajPreprocessInputCoordinates(): Error on opening\n");
    printfone("input coordinate file (%s)\n", currentCoordinateInfo->filename);
    fprintf(stdout,"ERROR: Rank %i failed to open file.\n",worldrank);
  }
#ifdef MPI
  // Check that all threads could open file - exit if any of them could not.
  i=0;
  MPI_Allreduce(&err,&i,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if (i>0) return 1;  
#endif 
  if (err!=0) return 1;


  /*
   *  preprocess the input coordinates as necessary 
   *  (i.e. to remove titles, etc.)
   */

  switch(currentCoordinateInfo->type) {

  case COORD_PDB:
    printfone("\nProcessing PDB file %s\n\n",
	    currentCoordinateInfo->filename);
    break;

  case COORD_BINPOS:
    printfone("\nProcessing BINPOS file %s\n\n",
	    currentCoordinateInfo->filename);
    break;

  case COORD_AMBER_RESTART:
    printfone("\nProcessing AMBER restart file %s\n\n",
	    currentCoordinateInfo->filename);
    break;

  case COORD_AMBER_TRAJECTORY:
    if (worldrank == 0)
      printfone("\nProcessing AMBER trajectory file %s\n\n",
	      currentCoordinateInfo->filename);
    break;

  case COORD_AMBER_REMD:
    printfone("\nProcessing AMBER REMD trajectory (new format) ");
    if (currentCoordinateInfo->isREMDTRAJ==0)
      fprintf(stdout,"file %s",currentCoordinateInfo->filename);
    fprintf(stdout,"\n\n");
    break;

  case COORD_CHARMM_TRAJECTORY:
    printfone("\nProcessing CHARMM trajectory file %s\n\n",
	    currentCoordinateInfo->filename);
    break;
  } 

  /* REMDTRAJ : If replica temperature processing was requested, 
   *            need to open the other replica files. 
   */
  if (currentCoordinateInfo->isREMDTRAJ>0) {
    printfone("  REMDTRAJ: Opening files %s -> %s%0*i\n",
            currentCoordinateInfo->filename,currentCoordinateInfo->baseFilename,
            currentCoordinateInfo->EXTwidth,
            currentCoordinateInfo->firstREMDTRAJ+(currentCoordinateInfo->numREMDTRAJ-1));

    /*NOTE: Currently spot 0 is not used for not NETCDF; inefficient */
    for (i=1; i<currentCoordinateInfo->numREMDTRAJ; i++) {
      fprintf (stdout,"    Opening %s\n",currentCoordinateInfo->REMDtraj[i]->filename);
      if (openTraj(currentCoordinateInfo->REMDtraj[i])!=0) {
        fprintf(stdout,"    REMDTRAJ: Error opening file.");
        return 1;
      }
        
     
    }
    fprintf(stdout,"  Done.\n");
  } 
  /* END REMDTRAJ */
      
  switch ( currentCoordinateInfo->type ) {

  case COORD_AMBER_TRAJECTORY:
    
    /*
     *  we need to read the title line to set up for processing
     *  by readAmberTrajectory()
     */
    /* DAN ROE - Right now all threads are reading title...may not be necessary
     * or matter since all threads now seek anyway.
     */
    //if ( fgets(buffer, BUFFER_SIZE, currentCoordinateInfo->file) == NULL ) {
    if ( trajFile_fgets(buffer, BUFFER_SIZE, currentCoordinateInfo) == NULL ) {
      fprintf(stdout, "WARNING: Error on processing the title from the AMBER\n");
      fprintf(stdout, "trajectory (%s)\n", currentCoordinateInfo->filename);
      return 1;
    }

    break;

  case COORD_AMBER_REMD:
      /*
       *  First, we need to read the title line from the original file
       */
    if ( fgets(buffer, BUFFER_SIZE, currentCoordinateInfo->file) == NULL ) {
      printfone("WARNING: Error on processing the title from the AMBER\n");
      printfone("REMD trajectory (%s)\n", currentCoordinateInfo->filename);
      return 1;
    }
    /* 
     * Now, if replica temperature trajectory processing was requested, read 
     * the title line from the other replica files.
     */
    if (currentCoordinateInfo->isREMDTRAJ>0) {
      for (i=1; i<currentCoordinateInfo->numREMDTRAJ; i++) {
        if ( fgets(buffer, BUFFER_SIZE, currentCoordinateInfo->REMDtraj[i]->file) == NULL ) {
          printfone("WARNING: Error on processing the title from the AMBER\n");
          printfone("REMD trajectory (%s)\n",currentCoordinateInfo->REMDtraj[i]->filename);
          return 1;
        }
      }
    }
    break;

  case COORD_BINPOS:

	if ( openbinpos( currentCoordinateInfo->file) < 0 ){
      printfone("Error on opening BINPOS file %s\n",
         currentCoordinateInfo->filename);
      return 1;
    }
    break;

  case COORD_CHARMM_TRAJECTORY:

    /*
     *  read in all the header information; this is done by calling with a
     *  negative value for the current set...
     */

    charmmTrajp = (charmmTrajectoryInfo **) safe_malloc(sizeof(charmmTrajectoryInfo *));
    *charmmTrajp = NULL;
    readCharmmTrajectory(currentCoordinateInfo->file, charmmTrajp,
			 NULL, NULL, NULL, NULL, -1);


    break;

  case COORD_UNKNOWN:

    printfone("WARNING: Attempting to process a coordinate file of unknown type\n");
    printfone("in %s, ignoring...\n", ROUTINE);
    return 1;
  }
  
  return 0;

}    

#undef  ROUTINE
#define ROUTINE "ptrajProcessInputCoordinates()"     

   void
ptrajProcessInputCoordinates(coordinateInfo *currentCoordinateInfo, 
			     ptrajState *state, 
			     double *X, double *Y, double *Z,
			     double *box, int set,
			     int *readCoordinates, int *processCoordinates)
{
  pdb_record *pdb;
  float *binpos;
  float fbox[6];
  int i,j,n_atoms,eof,err;
  charmmTrajectoryInfo **charmmTrajp;
  netcdfTrajectoryInfo *NCInfo;
  char xyz[3];
  float time;
  double repTemp;
  //FILE *currentRep;
  char buffer[BUFFER_SIZE];
  int k,currentNCID;
#ifdef MPI
  MPI_Offset start[3], count[3];
#else
  size_t start[3], count[3];
#endif

  /*
   *  NOTE: it is assumed that the box information is set to valid values on entry
   */



  /*
   *  READ IN THE CURRENT FILE, ONE SET AT A TIME
   */
  switch( currentCoordinateInfo->type ) {

  case COORD_PDB:

    i = loadPdb(currentCoordinateInfo->file, &pdb);
    j = getCoordinatesFromPdb(pdb, X, Y, Z);
    if (j != state->atoms) {

      /*
       *  on error, print warning but do not stop processing of this set.  This
       *  allows the specification of multiple reference sets
       */
      printfone("\nWARNING in %s: Unexpected number of atoms\n", ROUTINE);
      printfone(" encountered when reading PDB!\n");
      printfone("  coordinates read for %i atoms, expecting %i\n", j, state->atoms);
      
    }


    /*
     *  We assume that a PDB only contains ONE set of coordinates!!!
     */
    safe_fclose(currentCoordinateInfo->file);
    currentCoordinateInfo->file = NULL;
    *readCoordinates = 0;
    *processCoordinates = 1;
    break;
  
  case COORD_AMBER_RESTART:

    if ( getCoordinatesFromRestart(currentCoordinateInfo->file, 
				   X, Y, Z,
				   &box[0], &box[1], &box[2],
				   &box[3], &box[4], &box[5]) 
	 != state->atoms ) {
      printfone("WARNING in ptrajProcessInputCoordinates(): Unexpected\n");
      printfone("number of atoms encountered reading restrt file\n");
      ptrajCleanup();
      *readCoordinates = 0;
      *processCoordinates = 0;
      return;
    }

    /*
     *  get default box information if necessary
     */
    if ( state->IFBOX && (box[0] == 0.0 && box[1] == 0.0 && box[2] == 0.0) ) {
      box[0] = state->box[0];
      box[1] = state->box[1];
      box[2] = state->box[2];
    }

    if ( state->IFBOX && (box[3] == 0.0 && box[4] == 0.0 && box[5] == 0.0) ) {
      box[3] = state->box[3];
      box[4] = state->box[4];
      box[5] = state->box[5];
    }

    /*
     *  AMBER restrt files only contain a single set of coordinates
     */
    safe_fclose(currentCoordinateInfo->file);
    currentCoordinateInfo->file = NULL;
    *readCoordinates = 0;
    *processCoordinates = 1;
    break;

  case COORD_AMBER_TRAJECTORY:

    /*
     *  check to see if we've already loaded enough of the current
     *  trajectory file
     */

    if ( set > currentCoordinateInfo->stop && currentCoordinateInfo->stop != -1 ) {
      *readCoordinates = 0;
      
    } else {

      /*
       *  seek to next frame and read it in if file is not compressed
       *  NOTE: fseeko and long long are used for large file support.
       *  Might want to put in check to make sure this is supported by the system
       */
      if (currentCoordinateInfo->seekable)
        trajFile_fseek(currentCoordinateInfo,set-1);
      *readCoordinates=trajFile_fread(currentCoordinateInfo);
      bufferToXYZ(currentCoordinateInfo->buffer,currentCoordinateInfo->frameSize,X,Y,Z,
                  currentCoordinateInfo->numBox,box);

      /*if (currentCoordinateInfo->seekable)
	fseeko(currentCoordinateInfo->file, currentCoordinateInfo->titleSize + (long long) (set - 1) * currentCoordinateInfo->frameSize, SEEK_SET);

      *readCoordinates = readAmberTrajectory(currentCoordinateInfo->file,
					     state->atoms, X, Y, Z, box, set,
					     currentCoordinateInfo);*/
      
      if ( state->IFBOX && (box[0] == 0.0 && box[1] == 0.0 && box[2] == 0.0) ) {
	box[0] = state->box[0];
	box[1] = state->box[1];
	box[2] = state->box[2];
      }

      if ( state->IFBOX && (box[3] == 0.0 && box[4] == 0.0 && box[5] == 0.0) ) {
	box[3] = state->box[3];
	box[4] = state->box[4];
	box[5] = state->box[5];
      }

      *processCoordinates = 1;
    }

    if ( *readCoordinates == 0 ) {
      *processCoordinates = 0;
      //safe_fclose(currentCoordinateInfo->file);
      //currentCoordinateInfo->file = NULL;
      closeTraj(currentCoordinateInfo);
    }

    break;

  case COORD_AMBER_REMD:
      /*
       *  Check to see if we've already loaded enough of the current
       *  trajectory file.
       */
    if ( set > currentCoordinateInfo->stop ) {
      *readCoordinates = 0;        
    } else {
      /*
       *  scan through all replica files for target temp 
       *  if isREMDTRAJ==0 then we just want this replica 
       */
      err=1;
      //fprintf(stdout,"REMD DEBUG: remtrajtemp= %lf\n",currentCoordinateInfo->remdtrajtemp);
      if (currentCoordinateInfo->isREMDTRAJ==0) 
        j=1;
      else
        j=currentCoordinateInfo->numREMDTRAJ;
      for (i=0; i<j; i++) {
        /* REMDtraj 0 points to this coordinateInfo structure
	 *  Read in REMD header line 
         * <REMD> <Replica> <Exch#> <Step#> <Temp0>
	 */
        fscanf(currentCoordinateInfo->REMDtraj[i]->file,"%*s %*s %*s %*s %lf",&repTemp);
        if (prnlev>3) fprintf(stdout,"REMD DEBUG: Rep%03i temp= %lf\n",i,repTemp); 
     
	  /*
	   *  if this is our target temp or dont care about temp, read coords 
	   */
        if ( (repTemp==currentCoordinateInfo->remdtrajtemp)||
             (currentCoordinateInfo->isREMDTRAJ==0)       ) {
          if ((prnlev>3) && (currentCoordinateInfo->isREMDTRAJ>0))
            fprintf(stdout,
	            "  REMD target temperature (%lf) found (%lf) at replica index %i\n",
	     	    currentCoordinateInfo->remdtrajtemp,repTemp,i);
          err=0;
          *readCoordinates = readAmberTrajectory_nobuffer(currentCoordinateInfo->REMDtraj[i]->file,
                                                 state->atoms, X, Y, Z, box, set,
                                                 currentCoordinateInfo->isBox);
          // Store current replica Temp in ptraj state for possible write out
          state->temp0=repTemp;
        } else {
            // Empty read past coords and box info 
          for (k=0; k<currentCoordinateInfo->linesperset; k++)
            fgets(buffer,BUFFER_SIZE,currentCoordinateInfo->REMDtraj[i]->file);
        }
      } // End for loop over replica files
      /* 
       * If we found the target temperature, err will be 0. If not, generate
       * an error message and exit.
       */
      if (err==1) {
        *readCoordinates = 0;
        *processCoordinates = 0;
        fprintf(stdout,"\nREMDTRAJ: Final repTemp value read= %lf, set %i\n",repTemp,set);
        fprintf(stdout,"Could not find target %lf in any of the replica trajectories.\n",
                currentCoordinateInfo->remdtrajtemp);
        fprintf(stdout,"Check that all replica trajectory files were found and that\n");
        fprintf(stdout,
                "none of the trajectories are corrupted (e.g. missing a temperature).\n");
        fprintf(stderr,"%s: Target replica temperature not found in traj!",ROUTINE);
        return;
      }

      if ( state->IFBOX && (box[0] == 0.0 && box[1] == 0.0 && box[2] == 0.0) ) {
        box[0] = state->box[0];
        box[1] = state->box[1];
        box[2] = state->box[2];
      }
      if ( state->IFBOX && (box[3] == 0.0 && box[4] == 0.0 && box[5] == 0.0) ) {
        box[3] = state->box[3];
        box[4] = state->box[4];
        box[5] = state->box[5];
      }
      *processCoordinates = 1;
    }
    break;

  case COORD_AMBER_NETCDF:

    /*
     *  Prepare to load up a frame of NetCDF data and
     *  assume initially that we have successfully
     *  loaded up a set of coordinates
     */

#ifdef BINTRAJ
    NCInfo = (netcdfTrajectoryInfo *) currentCoordinateInfo->NCInfo;
    *readCoordinates = 1;
    *processCoordinates = 1;

      /*
       *  Check to see if we've loaded enough frames
       */

    if ( set > currentCoordinateInfo->stop ) {

      *readCoordinates = 0;
      *processCoordinates = 0;
  
    } else {

        /*
         *  unnecessary (but useful) sanity check: we should never trigger this since 
         *  this routine should not be called if coordinates are not present
         */

      if (NCInfo->coordinateVID == 0) {
	*readCoordinates = 0;
	printfone("Expecting coordinates in NetCDF file %s, none present...\n",
		currentCoordinateInfo->filename);
      }

      if ( *readCoordinates ) {
	if ( X == NULL || Y == NULL || Z == NULL ) 
	  error(ROUTINE, "coordinate arrays are NULL\n");

	/*
	 *  allocate space for coordinates from NetCDF file if necessary
         * NOTE: Should this be done somewhere else?
	 */
	if ( NCInfo->R == NULL ) {
	  NCInfo->R = (float *) safe_malloc( sizeof(float) * state->atoms * 3 );
	}

        /* DEFAULT NCID: If this is REMDTRAJ currentNCID will be set to the 
         * file containing the target temp.
         */
        currentNCID = NCInfo->ncid;
        /* DAN ROE: If this netcdf file has temperatures and we are not 
         * processing replica trajectories, get the temperature and store
         * it for possible write out.
         */
        if (currentCoordinateInfo->isREMDTRAJ==0 && NCInfo->TempVarID!=-1) {
          start[0]=set-1;
          count[0]=1;
#ifdef MPI
          err=ncmpi_get_vara_double(NCInfo->ncid, NCInfo->TempVarID,start,count,&repTemp);
#else
          err=nc_get_vara_double(NCInfo->ncid, NCInfo->TempVarID,start,count,&repTemp);
#endif
          if (err!=NC_NOERR)
            fprintf(stderr,"Could not get replica temperature from file.\n");
          /* DAN TEST: Store repTemp in ptraj state for possible writeout */
          state->temp0=repTemp;

        /* DAN REMDTRAJ - Which ncid has our target temperature? 
         *                Note: Temperatures have been checked
         *                for in checkCoordinates and setupREMDTRAJ.
         */
        } else if (currentCoordinateInfo->isREMDTRAJ>0) {
          /*fprintf(stdout,"\nTEMP SEARCH\n"); *DEBUG */
          start[0]=set-1;
          count[0]=1;

          j=-1;  // j==-1 after loop indicates error
          for (i=0; i<currentCoordinateInfo->numREMDTRAJ; i++) {
            k=currentCoordinateInfo->REMDtraj[i]->NCInfo->ncid;

#ifdef MPI
            err=ncmpi_get_vara_double(k, NCInfo->TempVarID,start,count,&repTemp);
#else
            err=nc_get_vara_double(k, NCInfo->TempVarID,start,count,&repTemp);
#endif
            if (err!=NC_NOERR) {
              fprintf(stderr,"Could not get replica temperature from file# %i\n",i);
              break;
            }
            if (prnlev>0) fprintf(stdout,"  Replica %i mytemp=%lf  ",i,repTemp);  //DEBUG 
            if (repTemp==currentCoordinateInfo->remdtrajtemp) {
              if (prnlev>0) fprintf(stdout,
                    "  REMD target temperature (%lf) found (%lf) at replica index %i\n",
                    currentCoordinateInfo->remdtrajtemp,repTemp,i);
              currentNCID=k; // This is the NCID that has our temperature 
              // Store repTemp in ptraj state for possible writeout
              state->temp0=repTemp;
              j=0;
              break;
            }
          } // End for loop over replica trajectories 
          /*fprintf(stdout,"\n"); * DEBUG */
          /* If we didn't find out target temperature, bad news. */
          if (j==-1) {
            *readCoordinates = 0;
            *processCoordinates = 0;
            fprintf(stdout,"\nREMDTRAJ: Final repTemp value read= %lf, set %i\n",repTemp,set);
            fprintf(stdout,"Could not find target %lf in any of the replica trajectories.\n",
                    currentCoordinateInfo->remdtrajtemp);
            fprintf(stdout,"Check that all replica trajectory files were found and that\n");
            fprintf(stdout,
                    "none of the trajectories are corrupted (e.g. missing a temperature).\n");
            fprintf(stderr,"%s: Target replica temperature not found in traj!",ROUTINE);
            return;
          } 
        } 
    
	/*
	 *  scan in coordinates for this frame
	 *
	 *  !! NOTE: C arrays are opposite the F90/Fortran so reverse
	 *     the array bounds on start and count.
	 */


	/*
	 *  PF - This only needs to be done once when first opening the file.
	 *  Seems like a waste to do it for every frame.
	 */

	start[0] = 0;
	count[0] = 3;
	xyz[0] = (char) 0;
#ifdef MPI
	err = ncmpi_get_vara_text(currentNCID, NCInfo->spatialVID, start, count, xyz);
#else
	err = nc_get_vara_text(currentNCID, NCInfo->spatialVID, start, count, xyz);
#endif
	if (err != NC_NOERR) {
	  printfone("Yikes!  We should see X, Y and Z chars in the spatial VID: %s\n",
		  nc_strerror(err));
	} else if ( xyz[0] != 'x' || xyz[1] != 'y' || xyz[2] != 'z' ) {
	  printfone("Whoa Nellie & Cripes!!!  We should see 'x', 'y' and 'y' characters\n");
	  printfone("in the spatial variables of the AMBER NetCDF trajectory file unless\n");
	  printfone("someone altered the spec.  We see '%c', '%c' and '%c'\n",
		  xyz[0], xyz[1], xyz[2]);
	}

	start[0] = set-1;
	start[1] = 0;
	start[2] = 0;
	count[0] = 1;
	count[1] = state->atoms;
	count[2] = 3;
      
#ifdef MPI
	err = ncmpi_get_vara_float(currentNCID, NCInfo->coordinateVID, start, count, NCInfo->R);
#else
	err = nc_get_vara_float(currentNCID, NCInfo->coordinateVID, start, count, NCInfo->R);
#endif
	if (err != NC_NOERR) {
	  printfone("Error detected upon reading NetCDF coordinates on set %i: %s", 
		  set, nc_strerror(err));
	  *readCoordinates = 0;
	  *processCoordinates = 0;
           return;
	} else {
	
	  /*
	   *  PF - There should be a more efficient way to read in the coords
	   *  directly to the X, Y, Z arrays rather than copying.
	   */

	  for (i = 0, j = 0; i < state->atoms; i++, j += 3) {
	    X[i] = NCInfo->R[ j   ];
	    Y[i] = NCInfo->R[ j+1 ];
	    Z[i] = NCInfo->R[ j+2 ];
	  }
	}
      }

      if ( *readCoordinates && state->IFBOX) {

	start[0] = set-1;
	start[1] = 0;
	start[2] = 0;
	count[0] = 1;
	count[1] = 3;
	count[2] = 0;

#ifdef MPI
	err = ncmpi_get_vara_double(currentNCID, NCInfo->cellLengthVID, start, count, box);
	if (err != NC_NOERR) 
	  error(ROUTINE, "Getting cell lengths: %s", nc_strerror(err));
	err = ncmpi_get_vara_double(currentNCID, NCInfo->cellAngleVID, start, count, &box[3]);
	if (err != NC_NOERR)
	  error(ROUTINE, "Getting cell angles: %s", nc_strerror(err));
#else
	err = nc_get_vara_double(currentNCID, NCInfo->cellLengthVID, start, count, box);
	if (err != NC_NOERR) 
	  error(ROUTINE, "Getting cell lengths: %s", nc_strerror(err));
	err = nc_get_vara_double(currentNCID, NCInfo->cellAngleVID, start, count, &box[3]);
	if (err != NC_NOERR)
	  error(ROUTINE, "Getting cell angles: %s", nc_strerror(err));
#endif
	/*
	printfone("\nGOT BOX INFO: %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
		box[0], box[1], box[2], box[3], box[4], box[5]);
	*/


      }

      /*
       *  get time
       */
      start[0] = set-1;
      count[0] = 1;
#ifdef MPI
      err = ncmpi_get_vara_float(currentNCID, NCInfo->timeVID, start, count, &time);
#else
      err = nc_get_vara_float(currentNCID, NCInfo->timeVID, start, count, &time);
#endif
      if (err != NC_NOERR)
	printfone("Warning in %s, getting the time from the NetCDF file: %s\n",
		ROUTINE, nc_strerror(err));

      /*
      printfone("\nREAD IN TIME %10.5f\n", time);
      */


    }

    if ( *readCoordinates == 0 ) {
      *processCoordinates = 0;
      //safe_fclose(currentCoordinateInfo->file);
      //currentCoordinateInfo->file = NULL;
      closeTraj(currentCoordinateInfo);
    }
    
#endif

    break;

  case COORD_CHARMM_TRAJECTORY:

    if ( set > currentCoordinateInfo->stop ) {
      *readCoordinates = 0;
      *processCoordinates = 0;
    } else {
      charmmTrajp = (charmmTrajectoryInfo ** ) safe_malloc(sizeof(charmmTrajectoryInfo *));
      *charmmTrajp = (charmmTrajectoryInfo *) currentCoordinateInfo->info;

      *readCoordinates = 
        readCharmmTrajectory(currentCoordinateInfo->file, charmmTrajp,
                             X, Y, Z, box, set);
    }
    if ( *readCoordinates == 0 ) {
      *processCoordinates = 0;
      closeTraj(currentCoordinateInfo);
    }
    break;

  case COORD_BINPOS:

    n_atoms = state->atoms;
    binpos = safe_malloc(sizeof(float) * n_atoms * 3);
    *readCoordinates = readbinpos( currentCoordinateInfo->file, 
				   &n_atoms, binpos, &eof);
    j=0;
    for (i=0; i<n_atoms; i++) {
		X[i] = binpos[j];
		Y[i] = binpos[j+1];
		Z[i] = binpos[j+2];
        j += 3;
    }
    safe_free(binpos);

    if ( (eof == 1) || (set > currentCoordinateInfo->stop) ) {
	*readCoordinates = 0;
    }
    if ( *readCoordinates == 0 ) {
      *processCoordinates = 0;
      safe_fclose(currentCoordinateInfo->file);
      currentCoordinateInfo->file = NULL;
    } else {
      *processCoordinates = 1;
    }
    break;

  default:
    *readCoordinates = 0;
    *processCoordinates = 0;
    printfone("WARNING in ptrajProcessInputCoordinates(): Attempting to process\n");
    printfone("a file of unknown type (%s)", currentCoordinateInfo->filename);
    printfone("Ignoring this file...\n");
    
  }
}

   int
lesSize( int atoms )
{
    int natomCL = 0, i=0;
    
    for( i=0; i < atoms; i++ )
    {
        if( parm->lescnum[i] == 0 || parm->lescnum[i] == 1 ) 
        {
            natomCL++;
        }
    }
    
    return natomCL;
}
   int*
lesMask( int atoms )
{
    int i;
    int* mask;

    mask = (int *) safe_malloc( sizeof(int) * atoms );

    for( i=0; i < atoms; i++ )
    {    
        if( parm->lescnum[i] == 0 || parm->lescnum[i] == 1 )
        {
            mask[i] = 1;
        }
        else
        { 
            mask[i] = 0;
        }
    }

    return mask;
}

   void
lesSplit( int atoms, double *X, double *Y, double *Z, 
          int icopy, double* xrep, double* yrep, double* zrep )
{
    int ia=0, i=0;
    
    for( i=0; i < atoms; i++ )
    {
        if( parm->lescnum[i] == 0 || parm->lescnum[i] == icopy + 1 )
        {
            xrep[ia] = X[i];
            yrep[ia] = Y[i];
            zrep[ia] = Z[i];
            ia++;
        }
    }
}

   void
lesAverage( int atoms, int nlescopy, double *X, double *Y, double *Z, 
	    double* xrep, double* yrep, double* zrep )
{
    int ia=0, i=0;
    double sum_x = 0.0;
    double sum_y = 0.0;
    double sum_z = 0.0;
    
    for( i=0; i < atoms; i++ )
    {
        if( parm->lescnum[i] == 0 )
        {
            xrep[ia] = X[i];
            yrep[ia] = Y[i];
            zrep[ia] = Z[i];
            ia++;
        }
        else
        {
            if( parm->lescnum[i] == 1 )
            {
                sum_x = 0.0;
                sum_y = 0.0;
                sum_z = 0.0;
            }
            
            sum_x += X[i];
            sum_y += Y[i];
            sum_z += Z[i];

            if( parm->lescnum[i] == nlescopy )
            {
                xrep[ia] = sum_x / nlescopy;
                yrep[ia] = sum_y / nlescopy;
                zrep[ia] = sum_z / nlescopy;
                ia++;
            }
        }
    }
}



#undef  ROUTINE
#define ROUTINE "ptrajOutputCoordinates()"

   void
ptrajOutputCoordinates(coordinateInfo *outInfo, ptrajState *state,
		       int set, int append, int first, int last, int atoms,
		       double *X, double *Y, double *Z, double *box)
{
  /*
   *  outInfo:  a pointer to the coordinateInfo file structure for the file to be output
   *  state:    the current ptrajState information
   *  set:      the current set or frame being processed
   *  append:   a flag (if non-zero) specifying that the outInfo file should be appended
   *  first:    a flag (if non-zero) specifying that this is the first call for this file
   *  last:     a flag (if non-zero) specifying that this is the last call for this file
   *  atoms:    the number of atoms
   *  X, Y, Z:  the coordinates
   *  box:      the box information (if not NULL)
   */
 
  char buffer[BUFFER_SIZE];
  pdb_record *pdb;
  netcdfTrajectoryInfo *NCInfo;
  int i, j, err, status;
  int bufferSize;
#ifdef BINTRAJ
  int dimensionID[NC_MAX_VAR_DIMS];
#ifdef MPI
  MPI_Offset start[3], count[3];
#else
  size_t start[3], count[3];
#endif
  char xyz[3];
  char abc[15] = { 'a', 'l', 'p', 'h', 'a', 
		   'b', 'e', 't', 'a', ' ', 
		   'g', 'a', 'm', 'm', 'a' };
#endif
  float time;
#ifdef MPI
  MPI_Offset offset;
#endif
  /*
   *  Four basic tasks are performed by this routine that enable output of the
   *  processed coordinates into different trajectory formats as specified by the
   *  (coordType) outInfo->type values which are defined in trajectory.h.
   *
   *  If you add a new coordType, implement code for each of these functions!
   *
   *    If (last) {
   *        (1) POST-PROCESS AND CLOSE FILE
   *        return;
   *    }
   *
   *    If (first) {
   *        (2) OPEN FILE AND PRE-PROCESS AND/OR SET UP TO APPEND [AND DO NOT return]
   *    }
   *    (3) DUMP CURRENT FRAME/SET
   */

    /*
     *  Perform final clean-up (closing) of trajectory file
     */
  if (last) {

    switch(outInfo->type) {

    case COORD_AMBER_TRAJECTORY:
#ifdef MPI
      if (outInfo->isMPI) {
	//safe_free(outInfo->buffer);
        /* DAN ROE: The get position and set size calls seem to cause problems, maybe
         * because they are collective routines and are being called by master. Why
         * are these needed at all?
         */
	/*if (worldrank == 0) {
	  MPI_File_get_position_shared(*(MPI_File *) outInfo->file, &offset);
	  MPI_File_set_size(*(MPI_File *) outInfo->file, offset);
	}*/
	MPI_File_close((MPI_File *) outInfo->file);
        // DAN ROE: Free memory!
        safe_free(outInfo->file);
	outInfo->file = NULL;
	return;
      }
#endif
    case COORD_CHARMM_TRAJECTORY:
    case COORD_BINPOS:

      safe_fclose(outInfo->file);
      outInfo->file = NULL;
      return;

    case COORD_AMBER_NETCDF:
      if (closeTraj(outInfo)==1)
        printfone("Error closing NetCDF file.\n");
      return;

    default:
      return;
    }

  }

    /*
     *  Do recursive calls as necessary for outputing LES trajectories
     */
  if ( outInfo->les_action != LES_NONE && outInfo->les_status == LES_READY ) {
    int natomCL = lesSize( atoms );
  
    double* xrep = safe_malloc( sizeof( double ) * natomCL );
    double* yrep = safe_malloc( sizeof( double ) * natomCL );
    double* zrep = safe_malloc( sizeof( double ) * natomCL );

    outInfo->les_status = LES_DONE;
    if ( outInfo->les_action == LES_SPLIT ) {
      int icopy=0;
      for ( icopy=0; icopy < outInfo->nlescopy; icopy++ ) {
	if ( icopy > 0 && first != 0 ) first = 0;
	lesSplit( atoms, X, Y, Z, icopy, xrep, yrep, zrep );

	ptrajOutputCoordinates( outInfo, state, (set-1) * outInfo->nlescopy + icopy + 1, 
				append, first, last, natomCL, xrep, yrep, zrep, box );
      }
    } else {
      assert( outInfo->les_action == LES_AVERAGE );
          
      lesAverage( atoms, outInfo->nlescopy, X, Y, Z, xrep, yrep, zrep );

      ptrajOutputCoordinates( outInfo, state, set, append, first, last, natomCL, xrep, yrep, zrep, box );
    }

    safe_free( xrep );
    safe_free( yrep );
    safe_free( zrep );

    outInfo->les_status = LES_READY;

    return;
  }

  if (first) {

    switch( outInfo->type ) {
	    
    case COORD_PDB:
      /*
       *  file handling is done per frame
       */
      break;

    case COORD_AMBER_NETCDF:

       /*
        *  On first call, preprocess to create the NetCDF file and to dump header information,
        */
       if (openTraj(outInfo)==1) {
         printfone("Error opening NETCDF output trajectory.\n");
         outInfo->type = COORD_UNKNOWN;
         break;
       }
       if (NETCDF_setupOutput(outInfo,atoms)==1) {
         printfone("Error setting up NETCDF output trajectory.\n");
         outInfo->type = COORD_UNKNOWN;
         break;
       }

    break;

    case COORD_AMBER_TRAJECTORY:
      
      /* set up buffer  */
      bufferSize = (atoms * 3 * 8) + (atoms * 3 / 10 + 2);
      if (box != NULL)
	bufferSize += 6 * 8 + 1;
      outInfo->buffer = safe_malloc( bufferSize * sizeof *outInfo->buffer );
      
      /*
       *  preprocessing to open file, if necessary
       */
      if (append) {
	
#ifdef MPI
        // DAN ROE: No cast to malloc call, potentially dangerous.
	if (outInfo->isMPI) {
	  outInfo->file = safe_malloc(sizeof(MPI_File));
	  if ( MPI_File_open( MPI_COMM_WORLD, outInfo->filename, MPI_MODE_RDWR | MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, (MPI_File *) outInfo->file) ) {
	    printfone("WARNING in ptrajOutputCoordinates(): Error opening\n");
	    printfone("output coordinate file in append mode (%s)\n", outInfo->filename);
	    outInfo->type = COORD_UNKNOWN;
             // DAN ROE: Just setting to NULL doesnt free mem!
            safe_free(outInfo->file);
	    outInfo->file = NULL;
	  }
	} else {
#endif
	  if ( ! openFile(&outInfo->file, outInfo->filename, "a") ) {
	    printfone("WARNING in ptrajOutputCoordinates(): Error opening\n");
	    printfone("output coordinate file in append mode (%s)\n", outInfo->filename);
	    /* DAN ROE: Safe fclose not needed up file open failed!
	     *safe_fclose(outInfo->file);
	     * DAN ROE: If file couldnt open set type to unknown to avoid further processing
         */
	    outInfo->type = COORD_UNKNOWN;
	    outInfo->file = NULL;
	  }
#ifdef MPI
	}
#endif
      } else {
#ifdef MPI
	if (outInfo->isMPI) {
	  outInfo->file = safe_malloc(sizeof(MPI_File));
	  if ( MPI_File_open( MPI_COMM_WORLD, outInfo->filename, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, (MPI_File *) outInfo->file) ) {
	    printfone("WARNING in ptrajOutputCoordinates(): Error opening\n");
	    printfone("output coordinate file (%s)\n", outInfo->filename);
	    outInfo->type = COORD_UNKNOWN;
            // DAN ROE: FREE MEM!
            safe_free(outInfo->file);
	    outInfo->file = NULL;
	  } else {
	    /* Write trajectory title if file opened correctly, only master rank */
	    if (worldrank == 0) {
	      if (outInfo->title == NULL)
		sprintf(buffer, "ptraj generated trajectory\n");
	      else
		sprintf(buffer, "%s\n", outInfo->title);
	      MPI_File_write_shared(*(MPI_File *) outInfo->file, buffer, strlen(buffer), MPI_CHAR, MPI_STATUS_IGNORE);
	    }
	  }
	} else {
#endif
	  if ( ! openFile(&outInfo->file, outInfo->filename, "w") ) {
	    printfone("WARNING in ptrajOutputCoordinates(): Error opening\n");
	    printfone("output coordinate file (%s)\n", outInfo->filename);
	    /*safe_fclose(outInfo->file); */
	    outInfo->type = COORD_UNKNOWN;
	    outInfo->file = NULL;
	  } else {
	    /* Write trajectory title if file opened correctly */
	    if (outInfo->title == NULL) 
	      fprintf(outInfo->file, "ptraj generated trajectory\n");
	    else
	      fprintf(outInfo->file, "%s\n", outInfo->title);
	  }
#ifdef MPI
	}
#endif
      }

      /* DAN TEST: REMD trajectory output  */
      if (outInfo->isREMDTRAJ) 
        fprintf(stdout,"TRAJOUT: Replica temperature will be written to output trajectory.\n");
    
      break;

    case COORD_CHARMM_TRAJECTORY:

       /*
        *  preprocessing to open file and write header
        */
      if (append)
	error(ROUTINE, "Appending CHARMM trajectories is not implemented\n");

      if ( ! openFile(&outInfo->file, outInfo->filename, "w") ) {
	printfone("WARNING in ptrajOutputCoordinates(): Error opening\n");
	printfone("output coordinate file (%s)\n", outInfo->filename);
        outInfo->type = COORD_UNKNOWN;
        outInfo->file=NULL;
      }
      dumpCharmmTrajectory(outInfo->file, (charmmTrajectoryInfo *) outInfo->info,
			   atoms, NULL, NULL, NULL, NULL, -1);
      break;

    case COORD_AMBER_RESTART:

      if (append) {
	sprintf(buffer, "%s", outInfo->filename);
	status = openFile(&outInfo->file, buffer, "a");
      } else if (set < 0) {
	sprintf(buffer, "%s", outInfo->filename);
	status = openFile(&outInfo->file, buffer, "w");
      } else {
	sprintf(buffer, "%s.%i", outInfo->filename, set);
	status = openFile(&outInfo->file, buffer, "w");
      }

      if (status == 0) {
	printfone("WARNING in ptrajOutputCoordinates(): Could not open\n");
	printfone("output coordinate file %s for %s.  Not dumping to output file.\n",
		buffer, (append ? "appending" : "writing")); 
	safe_fclose(outInfo->file);
	outInfo->file = NULL;
      }
  
      break;
    
    case COORD_BINPOS:

       /*
        *  preprocessing; open input file and write magic header
        */

      if ( ! openFile(&outInfo->file, outInfo->filename, (append ? "a" : "w")) ) {
	printfone("WARNING in ptrajOutputCoordinates(): Error opening\n");
	printfone("output coordinate file (%s) for %s\n", outInfo->filename,
		(append ? "appending" : "writing"));
        outInfo->type = COORD_UNKNOWN;
      }
      fwrite( "fxyz", 4, 1, outInfo->file );

      break;
    
    default:

      break;
    }
    /*
     *  end of pre-processing on first call
     */
  }
  
  switch( outInfo->type ) {
	    
  case COORD_PDB:

    if (append) {
      sprintf(buffer, "%s", outInfo->filename);
      status = openFile(&outInfo->file, buffer, "a");
    } else if (set < 0) {
      sprintf(buffer, "%s", outInfo->filename);
      status = openFile(&outInfo->file, buffer, "w");
    } else {
      sprintf(buffer, "%s.%i", outInfo->filename, set);
      status = openFile(&outInfo->file, buffer, "w");
    }
    if (status == 0) {
      printfone("ptrajOutputCoordinates(): Could not %s PDB file %s\n", 
	      (append ? "append" : "write" ), buffer);
    } else {

      /*
       *  if append, add MODEL records
       */
      if (append) {
	fprintf(outInfo->file, "MODEL %8d\n", set);
      }

      if (outInfo->title != NULL)
	fprintf(outInfo->file, "REMARK %3d %47s (set %5d)\n", 1, outInfo->title, set);
      /*
	fprintf(outInfo->file, "TITLE     %s (set %i)\n", outInfo->title, set);
      */

      if( outInfo->les_action != LES_NONE )
      {
          int* mask = lesMask( state->atoms );
          pdb = ptrajStateToPdb(state, mask, outInfo->option2);
          free( mask );
      }
      else
      {   
          pdb = ptrajStateToPdb(state, NULL, outInfo->option2);
      }

      putCoordinatesInPdb(pdb, atoms, X, Y, Z);

      if (outInfo->option1 > 0) {
	   /*
            *  include charges and radii and dump out in higher precision
            */
	putQandRInPdb(pdb, outInfo->option1, outInfo->option2, 0);
	savePdbHigh(outInfo->file, pdb);
      } else {
	savePdb(outInfo->file, pdb);
      }
      if (append) fprintf(outInfo->file, "TER\nENDMDL\n");
      safe_fclose(outInfo->file);
      outInfo->file = NULL;
      safe_free( (void *) pdb );

    }
    break;

  case COORD_AMBER_NETCDF:

#ifdef BINTRAJ
    NCInfo = (netcdfTrajectoryInfo *) outInfo->NCInfo;
    if (NCInfo == NULL) return;

    if (NCInfo->R == NULL)
      NCInfo->R = safe_malloc(sizeof(double) * atoms * 3);

    for (i = 0, j = 0; i < atoms; i++, j += 3) {
      NCInfo->R[ j   ] = X[i];
      NCInfo->R[ j+1 ] = Y[i];
      NCInfo->R[ j+2 ] = Z[i];
    }

    start[0] = NCInfo->currentFrame;
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = atoms;
    count[2] = 3;
#ifdef MPI
    err = ncmpi_put_vara_float(NCInfo->ncid, NCInfo->coordinateVID, start, count, NCInfo->R);
#else
    err = nc_put_vara_float(NCInfo->ncid, NCInfo->coordinateVID, start, count, NCInfo->R);
#endif
    if (err != NC_NOERR)
      printfone("NetCDF error on output of coordinates in set %i: %s\n", set, nc_strerror(err));

    if (outInfo->isBox) {

      start[0] = NCInfo->currentFrame;
      start[1] = 0;
      start[2] = 0;
      count[0] = 1;
      count[1] = 3;
      count[2] = 0;

      if (prnlev > 4) {
	printfone("\nPUT BOX INFO: %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
		box[0], box[1], box[2], box[3], box[4], box[5]);
      }

#ifdef MPI
      err = ncmpi_put_vara_double(NCInfo->ncid, NCInfo->cellLengthVID, start, count, box);
#else
      err = nc_put_vara_double(NCInfo->ncid, NCInfo->cellLengthVID, start, count, box);
#endif
      if (err != NC_NOERR)
	printfone("NetCDF error on output of cell lengths in set %i: %s\n", set, nc_strerror(err));

#ifdef MPI
      err = ncmpi_put_vara_double(NCInfo->ncid, NCInfo->cellAngleVID, start, count, &box[3]);
#else
      err = nc_put_vara_double(NCInfo->ncid, NCInfo->cellAngleVID, start, count, &box[3]);
#endif
      if (err != NC_NOERR)
	printfone("NetCDF error on output of cell angles in set %i: %s\n", set, nc_strerror(err));

    }

    start[0] = NCInfo->currentFrame;
    count[0] = 1;
    /*
     *  TODO: Figure out time
     */
    time = 0.0;
#ifdef MPI
    err = ncmpi_put_vara_float(NCInfo->ncid, NCInfo->timeVID, start, count, &time);
#else
    err = nc_put_vara_float(NCInfo->ncid, NCInfo->timeVID, start, count, &time);
#endif
    if (err != NC_NOERR)
        printfone("NetCDF error on output of time in set %i: %s\n", set, nc_strerror(err));

    /* DAN TEST: Put Temp0, dont be lazy with pointers eventually  */
    if (outInfo->isREMDTRAJ==1) { 
#ifdef MPI
      err = ncmpi_put_vara_double(NCInfo->ncid, NCInfo->TempVarID, start, count, &(state->temp0));
#else
      err = nc_put_vara_double(NCInfo->ncid, NCInfo->TempVarID, start, count, &(state->temp0));
#endif
      if (err != NC_NOERR)
        printfone("NetCDF error on output of temperature in set %i: %s\n", set, nc_strerror(err));
    }

    NCInfo->currentFrame += worldsize;
#ifdef MPI
    ncmpi_sync(NCInfo->ncid);
#else
    nc_sync(NCInfo->ncid);
#endif

#endif
    break;


  case COORD_AMBER_TRAJECTORY:

    if (outInfo->isREMDTRAJ==1) {
      /* Format: REMD  current replica#, exchange#, step#, and mytargettemp */
      fprintf(outInfo->file,"REMD  %8i %8i %8i %8.3lf\n",0,set,set,state->temp0);
    }
    dumpAmberTrajectory(outInfo, atoms, X, Y, Z, outInfo->isBox ? box : NULL);
    break;

  case COORD_CHARMM_TRAJECTORY:

    dumpCharmmTrajectory(outInfo->file, (charmmTrajectoryInfo *) outInfo->info,
			 atoms, X, Y, Z, box, set);
    break;

  case COORD_AMBER_RESTART:

    if (append) {
      sprintf(buffer, "%s", outInfo->filename);
      status = openFile(&outInfo->file, buffer, "a");
    } else {
      sprintf(buffer, "%s.%i", outInfo->filename, set); 
      status = openFile(&outInfo->file, buffer, "w");
    }
    
    if (status == 0) {
      printfone("ptrajOutputCoordinates(): Could not %s coordinate file %s, not dumping...\n",
	      (append ? "append" : "write"), outInfo->filename);
      break;
    }

    dumpAmberRestart(outInfo->file, atoms, X, Y, Z, NULL, NULL, NULL, 
		     (outInfo->isBox ? box : NULL), outInfo->title);
    safe_fclose(outInfo->file);
    outInfo->file = NULL;
    break;
    
  case COORD_BINPOS:

       /*
        *  preprocessing; open input file and write magic header
        */

    writebinpos(outInfo->file, atoms, X, Y, Z);

    break;
    
  default:
	      
    return;
  }
}



#undef  ROUTINE
#define ROUTINE "parseHBondDonor()"

   void
parseHBondDonor(stackType **argumentStackPointer, ptrajState *state, int *hbondDonor)
{
  char *buffer, *buffer1;
  int *mask, i, j;
  Name atom, res;

  buffer = argumentStackKeyToString( argumentStackPointer, "mask", NULL );
     /*
      *  mask
      */
  if (buffer != NULL) {
    mask = processAtomMask(buffer, state);
    atomMaskIsActive(mask, state, &i, &j);
    if (i==0) {
      printfone("WARNING in ptraj, donor: No atoms selected (%s), ignoring...\n",
	      buffer);
      safe_free(buffer);
      safe_free(mask);
      return;
    }

    for (i=0; i < state->atoms; i++) {
      if (mask[i] != 0) {
	hbondDonor[i] = 1;
	if (prnlev > 0)
	  printfone("  DONOR: adding atom %5i, residue %4i, atom name %s\n", 
		  i+1, atomToResidue(i+1, state->residues, state->ipres), state->atomName[i]);
      }
    }
    safe_free(buffer);
    safe_free(mask);
    return;
  }
    
     /*
      *  resname atomname
      */
  buffer  = getArgumentString( argumentStackPointer, NULL );
  buffer1 = getArgumentString( argumentStackPointer, NULL );
  if (buffer != NULL && buffer1 != NULL) {

       /*
        *  copy in atom and residue names and pad with spaces as appropriate
        */
    for (i=0; i < 4; i++) {
      res[i] = (char) 0;
      atom[i] = (char) 0;
    }
    strncpy(res,  buffer,  5);
    strncpy(atom, buffer1, 5);
    res[4]  = (char) 0;
    atom[4] = (char) 0;
    i = 3;
    while (res[i] == (char) 0) {
      res[i] = ' ';
      i--;
    }
    i = 3;
    while (atom[i] == (char) 0) {
      atom[i] = ' ';
      i--;
    }
    
    for (i=0; i < state->atoms; i++) {
      j = atomToResidue(i+1, state->residues, state->ipres)-1;
      if (strcmp(atom, state->atomName[i]) == 0 &&
	  strcmp(res, state->residueName[j]) == 0) {
	hbondDonor[i] = 1;
	if (prnlev > 0)
	  printfone("  DONOR: adding atom %5i, residue %4i, atom name %s\n", 
		  i+1, j+1, state->atomName[i]);
      }
    }
    safe_free(buffer);
    safe_free(buffer1);
  } else {
    printfone("WARNING in ptraj, donor: either the residue or atom name specification is blank!\n");
  }
  return;
}


#undef  ROUTINE
#define ROUTINE "parseHBondAcceptor()"

   void
parseHBondAcceptor(stackType **argumentStackPointer, ptrajState *state, int *hbondAcceptor,
		   int *hbondAcceptorH1, int *hbondAcceptorH2, int *hbondAcceptorH3)
{
  char *buffer, *buffer1;
  int *mask, *mask1, i, j, i1, j1, k, stop;
  Name atom, atom1, res;

  buffer1 = NULL;
  buffer  = argumentStackKeyToString( argumentStackPointer, "mask", NULL );
  if (buffer != NULL) {
    buffer1 = getArgumentString( argumentStackPointer, NULL );
  }
     /*
      *  mask
      */
  if (buffer != NULL && buffer1 != NULL) {
    mask = processAtomMask(buffer, state);
    atomMaskIsActive(mask, state, &i, &i1);
    mask1 = processAtomMask(buffer1, state);
    atomMaskIsActive(mask1, state, &j, &j1);

    stop = 0;
    if (i==0) {
      printfone("WARNING in ptraj, acceptor: No heavy atom was selected (%s), ignoring...\n",
		buffer);
      stop = 1;
    }
    if (j==0) {
      fprintf(stdout,
	      "WARNING in ptraj, acceptor: No hydrogen atom was selected (%s), ignoring...\n",
	      buffer1);
      stop = 1;
    }

    if (i != j) {
      fprintf(stdout,
	      "WARNING in ptraj, acceptor: There is not a 1-1 correspondence between the\n");
      printfone("atom selection in the two masks %s and %s which contain %i and %i\n",
	      buffer, buffer1, i, j);
      printfone("atoms respectively.  Ignoring...\n");
      stop = 1;
    }

    if (stop) {
      safe_free(buffer);
      safe_free(buffer1);
      safe_free(mask);
      safe_free(mask1);
      return;
    }


    stop = i;
    for (i=0; i < stop; i++) {
    
      hbondAcceptor[i1] = 1;
      if (hbondAcceptorH1[i1] >= 0) {
	if (hbondAcceptorH2[i1] >= 0) {
	  if (hbondAcceptorH3[i1] >= 0) {
	    fprintf(stdout,
		    "WARNING in ptraj, acceptor: More than three hydrogens have been selected\n");
	    printfone("as acceptors for the atom %5i.  Ignoring the latest...\n",
		    i1+1);
	  } else
	    hbondAcceptorH3[i1] = j1;
	} else
	  hbondAcceptorH2[i1] = j1;
      } else
	hbondAcceptorH1[i1] = j1;
    
      if (prnlev > 0) {
	printfone("  ACCEPTOR: adding atom %5i, res %4i, atom name %s -- ",
		i1+1, atomToResidue(i1+1, state->residues, state->ipres), state->atomName[i1]);
	printfone("atom %5i, res %4i, atom name %s\n",
		j1+1, atomToResidue(j1+1, state->residues, state->ipres), state->atomName[j1]);
      }
      mask[i1] = 0;
      mask1[j1] = 0;
      while (mask[i1] == 0 && i1 < state->atoms) i1++;
      while (mask1[j1] == 0 && j1 < state->atoms) j1++;
      if (i1 == state->atoms || j1 == state->atoms) return;
    }
    safe_free(mask);
    safe_free(mask1);
    safe_free(buffer);
    safe_free(buffer1);
    return;
  } else {
    if (buffer1 != NULL || buffer != NULL) {
      
      printfone(
	      "WARNING in ptraj, acceptor: Error in mask specification.  Ignoring...\n");
      safe_free(buffer);
      safe_free(buffer1);
      return;
    }
  }
    
     /*
      *  resname atomname atomname
      */
  stop = 0;
  buffer  = getArgumentString( argumentStackPointer, NULL );
  if (stop || buffer == NULL)
    stop = 1;
  else {
    for (i=0; i < 4; i++) {
      res[i] = (char) 0;
    }
    strncpy(res,  buffer,  5);
    res[4]  = (char) 0;
    i = 3;
    while (res[i] == (char) 0) {
      res[i] = ' ';
      i--;
    }
    safe_free(buffer);
  }
  
  buffer  = getArgumentString( argumentStackPointer, NULL );
  if (stop || buffer == NULL)
    stop = 1;
  else {
    for (i=0; i < 4; i++) {
      atom[i] = (char) 0;
    }
    strncpy(atom,  buffer,  5);
    atom[4]  = (char) 0;
    i = 3;
    while (atom[i] == (char) 0) {
      atom[i] = ' ';
      i--;
    }
    safe_free(buffer);
  }

  buffer  = getArgumentString( argumentStackPointer, NULL );
  if (stop || buffer == NULL)
    stop = 1;
  else {
    for (i=0; i < 4; i++) {
      atom1[i] = (char) 0;
    }
    strncpy(atom1,  buffer,  5);
    atom1[4]  = (char) 0;
    i = 3;
    while (atom1[i] == (char) 0) {
      atom1[i] = ' ';
      i--;
    }
    safe_free(buffer);
  }

  if (stop == 0) {
    
    for (i=0; i < state->atoms; i++) {
      j = atomToResidue(i+1, state->residues, state->ipres)-1;
      if (strcmp(atom, state->atomName[i]) == 0 &&
	  strcmp(res, state->residueName[j]) == 0) {
	
	for (k=state->ipres[j]-1; k < state->ipres[j+1]-1; k++) {
	  if (strcmp(atom1, state->atomName[k]) == 0) {
	    hbondAcceptor[i] = 1;
	    if (hbondAcceptorH1[i] >= 0)
	      if (hbondAcceptorH2[i] >= 0)
		if (hbondAcceptorH3[i] >= 0) {
		  printfone("WARNING in ptraj, acceptor: More than three hydrogens ");
		  printfone("have been selected\nas acceptors for atom %i. ", i+1);
		  printfone("Ignoring the latest...\n");
		} else
		  hbondAcceptorH3[i] = k;
	      else
		hbondAcceptorH2[i] = k;
	    else
	      hbondAcceptorH1[i] = k;

	    if (prnlev > 0) {
	      printfone("  ACCEPTOR: adding atom %5i, res %4i, atom name %s -- ",
		      i+1, atomToResidue(i+1, state->residues, state->ipres), 
		      state->atomName[i]);
	      printfone("atom %5i, res %4i, atom name %s\n",
		      k+1, atomToResidue(k+1, state->residues, state->ipres), 
		      state->atomName[k]);
	    }
	  }
	}
      }
    }
  } else 
    printfone(
	    "WARNING in ptraj, acceptor: error in specification of res/atom selection\n");
}


/* DAN ROE 
 * setupREMDTRAJ(): This is called during a TRAJIN transform in ptrajSetupIO
 * to set up this trajectory for processing as an REMD trajectory. All other
 * replica trajectory files are searched for, and the base filename, starting
 * replica, and number of replicas is recorded so the files may be opened
 * later during processing.
 *
 * NOTE: took out filename arg, should be in info structure.
 * UPDATE: coordinateInfo now contains a REMDtraj field which is itself
 * a coordinateInfo**. This routine should now allocate an array of 
 * coordinateInfo pointers.
 * NOTE: Make REMDtraj[0] just point to info? Dangerous? 
 */
   int
setupREMDTRAJ(coordinateInfo *info, int totalAtoms) {
  coordinateInfo *trajInfo;
  char* EXT;             /* Store Replica traj filename extension   */
  int loop;              /* For looping over replica files          */
  char* repFilename;     /* Store replica base filename             */
  int i,j,k;

  /* DAN ROE: Check that if remdtraj was specified in trajin this format
   *           is suitable for REMD processing.
   */
  if (info->type==COORD_AMBER_NETCDF) {
    if (info->NCInfo->TempVarID==-1) {
      fprintf(stdout,"WARNING: setupREMDTRAJ(): remdtraj requested in trajin but this NETCDF trajectory does not contain temperatures!\n");
      return 1;
    }
  } else if (info->type!=COORD_AMBER_REMD) {
    fprintf(stdout,"WARNING: setupREMDTRAJ(): remdtraj requested in trajin but this is not an amber REMD trajectory!\n"); 
    return 1;
  }

  fprintf(stdout,"  REMDTRAJ: Using specified file as lowest replica: %s\n",info->filename);
  fprintf(stdout,"  REMDTRAJ: Frames at %lf K will be processed.\n",info->remdtrajtemp);
  /*
   *  Scan for additional REMD traj files.
   *  Assume the extension of given trajectory is the number
   *  of the lowest replica, and that the other files are in
   *  sequence (e.g. rem.000, rem.001, rem.002 etc).
   */
  loop=strlen(info->filename);
  /* Was the file zipped? */
  info->compressEXT=(char*) malloc(5*sizeof(char));
  strcpy(info->compressEXT,"");
  if ( (info->filename[loop-3]=='.')&&(info->filename[loop-2]=='g')&&(info->filename[loop-1]=='z') ) {
    fprintf(stdout,"  REMDTRAJ: File is gzipped.\n");
    strcpy(info->compressEXT,".gz");
    loop=loop-3;
  } else if (   (info->filename[loop-4]=='.')&&(info->filename[loop-3]=='b')
              &&(info->filename[loop-2]=='z')&&(info->filename[loop-1]=='2') ) {
    fprintf(stdout,"  REMDTRAJ: File is bzipped.\n");
    strcpy(info->compressEXT,".bz2");
    loop=loop-4;
  }

  /* First, find location of last '.' and store it in i*/
  for (j=0; j<loop; j++)
    if (info->filename[j]=='.') i=j;

  if (prnlev>0) {
    fprintf(stdout,"  REMDDEBUG: Last . in %s located at %i\n",info->filename,i);
    fprintf(stdout,"  REMDDEBUG: Allocating %i for extension\n",loop-i);
    fprintf(stdout,"  REMDDEBUG: EXTwidth=%i\n",loop-i-1);
  }

  /* Get filename extension */
  EXT=(char*) safe_malloc((loop-i)*sizeof(char));
  info->EXTwidth=loop-i-1;
  k=0;
  for (j=i+1; j<loop; j++)
    EXT[k++]=info->filename[j];
  EXT[k]='\0';

  if (prnlev>0) printfone("  REMDDEBUG: Replica extension is %s\n",EXT);

  /* Check that all digits in extension are numbers */
  for (j=0; j<k; j++) {
    if (isdigit(EXT[j])==0) {
      fprintf(stdout,
              "REMDTRAJ: WARNING: Character #%i (%c) in extension %s is not a number!\n",
              j,EXT[j],EXT);
      safe_free(EXT);
      return 1;  
    }
  }

  /* Look for the other replica files, assuming the name is basefilename.num or
   *   basefilename.num.gz
   */
  j=atoi(EXT);
  if (prnlev>0) fprintf(stderr,"  REMDDEBUG: index of first replica = %i\n",j);
  safe_free(EXT);
  info->firstREMDTRAJ=j;

  /*
   *  Allocate memory for the replica filenames.
   */
  k=strlen(info->filename);
  info->baseFilename=(char*) safe_malloc((k+1)*sizeof(char));
  repFilename=(char*) safe_malloc((k+1)*sizeof(char));

  /*
   *  Store base filename. 
   *  Variable 'i' still contains location of '.' before # extension
   */
  strncpy(info->baseFilename,info->filename,i+1);
  if (prnlev>0) fprintf(stderr,"  REMDDEBUG: base filename = %s\n",info->baseFilename);

  /* 
   * Search for a replica number lower than this. Correct functioning
   * of the replica code requires the file specified by trajin be the 
   * lowest # replica.
   */
  sprintf(repFilename,"%s%0*i%s",info->baseFilename,info->EXTwidth,j-1,info->compressEXT);
  trajInfo=checkCoordinates(repFilename,totalAtoms);
  if (trajInfo!=NULL) {
    fprintf(stdout,
            "  REMDTRAJ: WARNING: Replica# found lower than file specified with trajin!\n");
    printfone("            (Found %s)\n",repFilename);
    printfone("            trajin <file> remdtraj requires lowest # replica!\n");
    cleanTraj(trajInfo);
  }

  /*
   *  Now count REMD files 
   */
  // Allocate array for REMD trajs - set 0 to this info
  if (info->REMDtraj==NULL) {
    info->REMDtraj=(coordinateInfo**) malloc(sizeof(coordinateInfo*));
    info->REMDtraj[0]=info;
  }

  i=1;    /* # of replica files counted.                    */
  loop=1; /* Continue loop                                  */
  j++;    /* Replica file index, start at next replica file */
  fprintf(stdout,"  REMDTRAJ: Scanning for other REMD files.\n");
  while  ( loop == 1 ) {
    /*
     *  filename to scan for 
     */
    sprintf(repFilename,"%s%0*i%s",info->baseFilename,info->EXTwidth,j,info->compressEXT);
    trajInfo=checkCoordinates(repFilename,totalAtoms);
    if (trajInfo==NULL) {
      loop=0;
      if (prnlev>0) fprintf(stderr,"  REMDDEBUG: %s not found, exiting loop.\n",repFilename);
    } else {
      info->REMDtraj=(coordinateInfo**) realloc(info->REMDtraj, (i+1) * sizeof(coordinateInfo*));
      info->REMDtraj[i]=trajInfo;
      if (prnlev>0) fprintf(stderr,"    REMDDEBUG: Found %s\n",info->REMDtraj[i]->filename);
      i++; j++;
    }
    
  } /* End loop to find replica trajs */
  fprintf(stdout,"  REMDTRAJ: Found %i replica traj files.\n",i);
  info->numREMDTRAJ=i;
  safe_free(repFilename);
  // NOTE: should check here that all trajs have temperatures

  return 0;
}

/*
 *  ptrajSetupIO(): This is called upon receipt of the trigger strings
 *  that specify input and output files, reference structures and
 *  global variables that may be used by various actions (such as solvent
 *  information or hydrogen bonding donor/acceptor atoms).
 *  For I/O, the transformFileStack list of input files is set up as well
 *  as the output trajectory file.
 */

#undef  ROUTINE
#define ROUTINE "ptrajSetupIO()"

   void
ptrajSetupIO(stackType *argumentStack, actionType type)
{
  char *filename;
  char *buffer;
  char *title, *titlenew;
  int i, j, k, start, stop, offset;
  coordinateInfo *info;
  ptrajState *state, **statep;
  int *mask, byres, bytype, curres, *startsol, *stopsol;
  Name res;
  charmmTrajectoryInfo *charmmTraj, *ctrj;
  stackType *sp;
  coordinateInfo *infiles;
  double boxtmp;

  statep = ptrajCurrentState();
  state = *statep;

  /*
   *  Input/Output setup for ptraj.  The following command/arguments are parsed:
   *
   *  TRANSFORM_BOX
   *
   *    box [x value] [y value] [z value] [alpha value] [beta value] [gamma value] 
   *        [fixx] [fixy] [fixz] [fixalpha] [fixbeta] [fixgamma]
   *
   *  TRANSFORM_BENCHMARK
   *
   *    benchmark [out <filename>] [short]
   *
   *  TRANSFORM_PRNLEV
   *
   *    prnlev <value>
   *
   *  TRANSFORM_REFERENCE
   *
   *    reference filename
   *
   *  TRANSFORM_SOLVENT
   *
   *    solvent mask1 [mask2] ... [maskN] [byres | bytype | byname]
   *
   *  TRANSFORM_TRAJIN
   *
   *    trajin filename [start] [stop] [delta]
   *
   *  TRANSFORM_TRAJOUT
   *
   *    trajout filename [nobox] [ PDB | RESTART | BINPOS | NETCDF ] [append]
   *      [title <title>] [application <application>] [program <program>]
   *
   *  TRANSFORM_DONOR
   *
   *    donor {print | clear | <resname> <atomname> | mask <mask>}
   *
   *  TRANSFORM_ACCEPTOR
   *
   *    acceptor print | clear | 
   *       { <resname> <atomname1> <atomname2> } | 
   *       { mask <mask1> <mask2> }
   *
   */


  switch ( type ) {

  case TRANSFORM_PRNLEV:

    prnlev = getArgumentInteger(&argumentStack, 0.0);
    printfone("  PRNLEV: value is now %i\n", prnlev);
    return;

  case TRANSFORM_BENCHMARK:

    bench = 1;
    benchshort = argumentStackContains( &argumentStack, "short" );
    benchfile = argumentStackKeyToString( &argumentStack, "out", NULL );
    return;

  case TRANSFORM_BOX:

    boxtmp = argumentStackKeyToDouble(&argumentStack, "x", 0.0);
    if (boxtmp != 0.0) {
      state->IFBOX = 1;
      state->box[0] = boxtmp;
      printfone("  BOX X = %f\n", boxtmp);
    }
    boxtmp = argumentStackKeyToDouble(&argumentStack, "y", 0.0);
    if (boxtmp != 0.0) {
      state->IFBOX = 1;
      state->box[1] = boxtmp;
      printfone("  BOX Y = %f\n", boxtmp);
    }
    boxtmp = argumentStackKeyToDouble(&argumentStack, "z", 0.0);
    if (boxtmp != 0.0) {
      state->IFBOX = 1;
      state->box[2] = boxtmp;
      printfone("  BOX Z = %f\n", boxtmp);
    }
    boxtmp = argumentStackKeyToDouble(&argumentStack, "alpha", 0.0);
    if (boxtmp != 0.0) {
      state->IFBOX = 1;
      state->box[3] = boxtmp;
      printfone("  BOX ALPHA = %f\n", boxtmp);
    }
    boxtmp = argumentStackKeyToDouble(&argumentStack, "beta", 0.0);
    if (boxtmp != 0.0) {
      state->IFBOX = 1;
      state->box[4] = boxtmp;
      printfone("  BOX BETA  = %f\n", boxtmp);
    }
    boxtmp = argumentStackKeyToDouble(&argumentStack, "gamma", 0.0);
    if (boxtmp != 0.0) {
      state->IFBOX = 1;
      state->box[5] = boxtmp;
      printfone("  BOX GAMMA = %f\n", boxtmp);
    }

    if (argumentStackContains( &argumentStack, "fixx" ))
      state->boxfixed[0] = 1;
    if (argumentStackContains( &argumentStack, "fixy" ))
      state->boxfixed[1] = 1;
    if (argumentStackContains( &argumentStack, "fixz" ))
      state->boxfixed[2] = 1;
    if (argumentStackContains( &argumentStack, "fixalpha" ))
      state->boxfixed[3] = 1;
    if (argumentStackContains( &argumentStack, "fixbeta" ))
      state->boxfixed[4] = 1;
    if (argumentStackContains( &argumentStack, "fixgamma" ))
      state->boxfixed[5] = 1;

    return;

  case TRANSFORM_DONOR:

       /*
        *  print
        */
    if (argumentStackContains( &argumentStack, "print" )) {
      if (hbondDonor != NULL) {
	printHBondMask(hbondDonor, NULL, NULL, NULL, state);
      }
      return;
    }

       /*
        *  clear
        */
    if (argumentStackContains( &argumentStack, "clear" )) {
      if (hbondDonor != NULL) {
	pushBottomStack(&hbondDonorStack, (void *) hbondDonor);
	hbondDonor = NULL;
      }
      return;
    }

       /*
        *  allocate space for the hbondDonor if necessary...
        */
    if (hbondDonor == NULL) {
      hbondDonor = (int *) safe_malloc(sizeof(int) * state->atoms);
      for (k=0; k < state->atoms; k++)
	hbondDonor[k] = 0;
    }

    parseHBondDonor(&argumentStack, state, hbondDonor);

       /*
        *  make sure at least one donor has been chosen or free the memory
        */
    atomMaskIsActive(hbondDonor, state, &i, &j);
    if (i == 0) {
      safe_free(hbondDonor);
      hbondDonor = NULL;
    }
    return;


  case TRANSFORM_ACCEPTOR:

       /*
        *  print
        */
    if (argumentStackContains( &argumentStack, "print" )) {
      if (hbondAcceptor != NULL) {
	printHBondMask(hbondAcceptor, hbondAcceptorH1, hbondAcceptorH2, hbondAcceptorH3, state);
      }
      return;
    }

       /*
        *  clear
        */
    if (argumentStackContains( &argumentStack, "clear" )) {
      if (hbondAcceptor != NULL) {
	pushBottomStack(&hbondAcceptorStack, (void *) hbondAcceptor);
	pushBottomStack(&hbondAcceptorH1Stack, (void *) hbondAcceptorH1);
	pushBottomStack(&hbondAcceptorH2Stack, (void *) hbondAcceptorH2);
	pushBottomStack(&hbondAcceptorH3Stack, (void *) hbondAcceptorH3);
	hbondAcceptor = NULL;
	hbondAcceptorH1 = NULL;
	hbondAcceptorH2 = NULL;
	hbondAcceptorH3 = NULL;
      }
      return;
    }

       /*
        *  allocate space for the hbondAcceptor lists if necessary...
        */
    if (hbondAcceptor == NULL) {
      hbondAcceptor = (int *) safe_malloc(sizeof(int) * state->atoms);
      hbondAcceptorH1 = (int *) safe_malloc(sizeof(int) * state->atoms);
      hbondAcceptorH2 = (int *) safe_malloc(sizeof(int) * state->atoms);
      hbondAcceptorH3 = (int *) safe_malloc(sizeof(int) * state->atoms);
      for (k=0; k < state->atoms; k++) {
	hbondAcceptor[k] = 0;
	hbondAcceptorH1[k] = -1;
	hbondAcceptorH2[k] = -1;
	hbondAcceptorH3[k] = -1;
      }
    }

    parseHBondAcceptor(&argumentStack, state, hbondAcceptor, 
		       hbondAcceptorH1, hbondAcceptorH2, hbondAcceptorH3);

       /*
        *  make sure at least one acceptor has been chosen or free the memory
        */
    atomMaskIsActive(hbondAcceptor, state, &i, &j);
    if (i == 0) {
      safe_free(hbondAcceptor);
      safe_free(hbondAcceptorH1);
      safe_free(hbondAcceptorH2);
      safe_free(hbondAcceptorH3);
      hbondAcceptor = NULL;
      hbondAcceptorH1 = NULL;
      hbondAcceptorH2 = NULL;
      hbondAcceptorH3 = NULL;
    }
    return;

  case TRANSFORM_REFERENCE:

    filename = getArgumentString( &argumentStack, NULL );
    if (filename == NULL) {
      printfone("WARNING in ptrajSetupIO(): reference command lacks a filename!\n");
      return;
    }

    info = checkCoordinates(filename, state->atoms);
    safe_free(filename);

    if ( info == NULL )
      return;

    info->mask = NULL;
    info->offset = 1;
    info->start = 1;
    info->stop = -1;
    info->option1 = 0;
    info->option2 = 0;

    if ( ptrajPreprocessInputCoordinates(info) ) {
      safe_free(info->filename);
      safe_free(info->info);
      safe_free(info->NCInfo); // DEBUG
      return;
    }

    info->x = (double *) safe_malloc(sizeof(double) * state->atoms);
    info->y = (double *) safe_malloc(sizeof(double) * state->atoms);
    info->z = (double *) safe_malloc(sizeof(double) * state->atoms);

    /*
     *  note, we initialize the memory here since the checks to disable processing
     *  of truncated sets have been removed, i.e. it is possible to read in a truncated
     *  PDB reference structure noting that the truncated part will be replaced by zeros.
     *  If the set is truncated, the user will be warned...
     */
    for (i=0; i < state->atoms; i++) {
      info->x[i] = 0.0;
      info->y[i] = 0.0;
      info->z[i] = 0.0;
    }
    ptrajProcessInputCoordinates(info, state, info->x, info->y, info->z, state->box,
				 1, &i, &j);
    if (j == 0) {
      safe_free(info->x);
      safe_free(info->y);
      safe_free(info->z);
      cleanTraj(info);
      /*safe_free(info->filename);
      safe_free(info->info);
      safe_free(info->NCInfo); // DEBUG*/
      return;
    }

    referenceInfo = info;
    pushBottomStack( &transformReferenceStack, (void *) info );
    
    return;

  case TRANSFORM_SOLVENT:

    /*
     *  solvent [byres | bytype | byname] mask1 [mask2] [mask3] ...
     */

    byres  = argumentStackContains( &argumentStack, "byres");
    bytype = argumentStackContains( &argumentStack, "bytype");
    if ( argumentStackContains( &argumentStack, "byname") == 1) {
      byres  = 0;
      bytype = 0;
    }

    if (bytype == 1) {
      printfone("WARNING in ptrajSetupIO(): solvent \"bytype\" option has not yet\n");
      printfone("been implemented.  Defaulting to byname...\n");
      bytype = 0;
    }

    buffer = (char *) argumentStack->entry;

       /*
        *  if we've exhausted the argument stack, assume this is a call to
	*  solvent with no arguments; in this case we do not alter the solvent
	*  information, we just print a summary of it...
        */
    if (buffer[0] == (char) 0) {
      ptrajPrintState(state);
      return;
    }

    /*
     *  re-initialize the solvent information
     */
    state->solventMolecules = 0;
    state->solventAtoms = 0;
    if (state->solventMask != NULL) safe_free(state->solventMask);
    state->solventMask = NULL;
    if (state->solventMoleculeStart != NULL) safe_free(state->solventMoleculeStart);
    state->solventMoleculeStart = NULL;
    if (state->solventMoleculeStop != NULL) safe_free(state->solventMoleculeStop);
    state->solventMoleculeStop = NULL;

    state->solventMask = (int *) safe_malloc(sizeof(int) * state->atoms);
    for (i=0; i < state->atoms; i++)
      state->solventMask[i] = 0;
    startsol = (int *) safe_malloc(sizeof(int) * state->atoms);
    stopsol  = (int *) safe_malloc(sizeof(int) * state->atoms);

    buffer = NULL;
    mask = NULL;
    while ( (buffer = getArgumentString(&argumentStack, NULL)) != NULL ) {

      if (byres) {
	/*
	 *  load up solvent by residue using the residues specified in
	 *  the input masks...
	 */

	mask = processAtomMask(buffer, state);
	printfone("       Searching for solvent by mask ");
	printAtomMask(stdout, mask, state);
	printfone("\n");
	
	for (i=0; i < state->atoms; i++) {

	  if (mask[i]) {
	    curres = atomToResidue(i+1, state->residues, state->ipres)-1;
	    j = isActiveResidue(curres, mask, state);

	    if (j > 0) {
	      /*
	       *  all the atoms in "curres" are active, hence add this
	       *  solvent molecule to the list of solvent.  Note that
	       *  isActiveResidue returns the number of active atoms
	       *  found (if not zero)...
	       */

	      state->solventAtoms += j;
	      startsol[state->solventMolecules] = i;
	      stopsol[ state->solventMolecules] = i+j;
	      state->solventMolecules++;
	      for (k=i; k < i+j; k++)
		state->solventMask[k] = 1;
	      i += j-1;
	    }
	  }
	}
    
      } else if (bytype) {

	/*
	 *  search for solvent using the mask to define
	 *  what a representative solvent molecule is...
	 *
	 *  CURRENTLY NOT IMPLEMENTED
	 */
	mask = processAtomMask(buffer, state);

      } else {
	/*
	 *  search for solvent by residue name
	 */

           /*
            *  copy the residue name from the buffer into "res" and pad
            *  with spaces to conform to standard atom/residue naming
            */
	for (i=0; i < 4; i++)
	  res[i] = (char) 0;
	strncpy(res, buffer, 5);
	res[4] = (char) 0;
	i = 3;
	while (res[i] == (char) 0) {
	  res[i] = ' ';
	  i--;
	}

	printfone("       Searching for solvent by residue name %s\n", res);

           /*
            *  loop over all residues and check for matches
            */
	for (i=0; i < state->residues; i++) {
	  if (strcmp(res, state->residueName[i]) == 0) {
	    /*
	     * add this residue to the list of solvent
	     */
	    j = state->ipres[i+1]-state->ipres[i];

	    state->solventAtoms += j;
	    startsol[state->solventMolecules] = state->ipres[i]-1;
	    stopsol[ state->solventMolecules] = state->ipres[i+1]-1;
	    state->solventMolecules++;
	    for (k=state->ipres[i]-1; k < state->ipres[i+1]-1; k++)
	      state->solventMask[k] = 1;
	  }
	}

      }

      if (mask!=NULL) free(mask);
      safe_free(buffer);
    }

       /*
        *  update the solvent information and free any unnecessary memory
        */
    state->solventMoleculeStart = (int *) 
      safe_malloc(sizeof(int) * state->solventMolecules);
    state->solventMoleculeStop = (int *) 
      safe_malloc(sizeof(int) * state->solventMolecules);
    for (i=0; i < state->solventMolecules; i++) {
      state->solventMoleculeStart[i] = startsol[i];
      state->solventMoleculeStop[i]  = stopsol[i];
    }

    safe_free(buffer);
    safe_free(stopsol);
    safe_free(startsol);
    if (prnlev > 0)
      ptrajPrintState(state);
    return;

  case TRANSFORM_TRAJIN:

    filename = getArgumentString( &argumentStack, NULL );
    if (filename == NULL) {
      printfone("WARNING in ptrajSetupIO: trajin command lacks a filename!\n");
      return;
    }

    /*
     *  check the coordinates and push the information to the bottom
     *  of the transformFileStack
     */
    printfone("  Checking coordinates: %s\n", filename);
    fflush(stdout);

#ifdef MPI
    t1 = MPI_Wtime();
#else
    t1 = clock();
#endif

    info = checkCoordinates(filename, state->atoms);

#ifdef MPI
    t2 = MPI_Wtime();
    checkInputTime += (t2 - t1);
#else
    t2 = clock();
    checkInputTime += (t2 - t1) / CLOCKS_PER_SEC;
#endif    

    if (info == (coordinateInfo *) NULL) {
      printfone("\tCould not process trajectory %s\n",filename);
      safe_free(filename);
      return;
    }

    /* DAN ROE: Set isREMDTRAJ and remdtrajtemp here since the previous call to 
     *            checkCoordinates() gives a new info structure.
     *          isREMDTRAJ specifies that we want to look for frames at the temperature
     *            given by remdtrajtemp (i.e. it is processed as a temperature trajectory). 
     *            Otherwise the trajectory is treated normally, i.e. as a replica trajectory.
     */
    info->isREMDTRAJ=argumentStackContains(&argumentStack,"remdtraj");
    info->remdtrajtemp=argumentStackKeyToDouble(&argumentStack, "remdtrajtemp", 0.0);
    /* wxw: check for a RXSGLD trajectory if it is not a REMD trajectory */
    if(info->isREMDTRAJ==0)info->isREMDTRAJ=argumentStackContains(&argumentStack,"rxsgldtraj");
    if(info->remdtrajtemp==0.0)info->remdtrajtemp=argumentStackKeyToDouble(&argumentStack, "rxsgldid", 0.0);

    /* DAN ROE: Moved REMD setup routine to its own place */
    if (info->isREMDTRAJ==1) { 
      if ( setupREMDTRAJ(info,state->atoms)!=0 ) {
        info->isREMDTRAJ=0;
        fprintf(stdout,"                   Trajectory will be processed as replica traj.\n");
      }
    }

    /*
     *  set start, stop and offset if they were specified and relevent
     */

    start  = getArgumentInteger( &argumentStack,  1 );
    stop   = getArgumentInteger( &argumentStack, -1 );
    offset = getArgumentInteger( &argumentStack,  1 );

    if (stop > 0) {
      if (info->stop == -1)
	info->stop = stop;
      else if (stop > info->stop) {
	printfone(
		"FYI %s: trajin stop value (%i) is greater than the number of sets read in.\n", 
		ROUTINE, stop);
	printfone("Setting stop to the maximum value (%i)\n",
		info->stop);
      } else {
	info->stop = stop;
      }
    }

    if (start > 1) {
      if (info->stop > 0 && start > info->stop) {
	printfone("WARNING in %s: trajectory start is > stop; no\n", ROUTINE);
	printfone("configurations will be processed\n");
      }
      info->start = start;
    }

    info->offset = offset;
    if (info->offset != 1 && info->offset > info->stop - info->start) {
      printfone("WARNING in %s: Offset is so large that only 1 set\n", ROUTINE);
      printfone("  will be processed...\n");
    }

#ifdef MPI
    /*
     *  Current limitation of parallel ptraj is that the number of input frames must
     *  be evenly divisible by offset * number of ranks. The reason is that file i/o
     *  takes place using MPI I/O and lines are printed with MPI_File_write_ordered
     *  where each thread must make a call. If there is an odd number of frames then
     *  at least one of the threads will not write (call) and it will hang waiting
     */
    if ( ! checkDivisibility(&info->stop, info->start, info->offset) ) {
      if (worldrank == 0) {
	printf("WARNING in ptraj(): Number of sets to be processed is not a multiple of\n");
	printf("offset * number of ranks for input file %s.\n", filename);
	printf("This is a current limitation of parallel ptraj. Stop frame being set to: %d\n", info->stop);
      }
    }
#endif    

    pushBottomStack( &transformFileStack, (void *) info );
    return;

  case TRANSFORM_TRAJOUT:

    /*
     *  Set up the coordinateInfo structure for the file to be
     *  output.  Note that we do not actually open the filename
     *  yet as this is done by ptrajOutputCoordinates()
     */

    filename = getArgumentString( &argumentStack, NULL );
    if (filename == NULL) {
      printfone("WARNING in %s: trajout command lacks a filename!\n", ROUTINE);
      safe_free(filename);
      return;
    }

    info = (coordinateInfo *) safe_malloc(sizeof(coordinateInfo));
    INITIALIZE_coordinateInfo(info);
    info->filename = copyString(filename);
    if (argumentStackContains( &argumentStack, "nobox" ))
      info->isBox = 0;
    else
      info->isBox = state->IFBOX;
    info->append = argumentStackContains( &argumentStack, "append" );
    // Set access mode
    if (info->append==1)
      info->accessMode=2;
    else
      info->accessMode=1;

    /* DAN ROE: Specify write replica temperature in output trajectory  */
    info->isREMDTRAJ=argumentStackContains(&argumentStack,"remdtraj");
    // Disable replica temperature out for multiptraj
    if (info->isREMDTRAJ==1 && worldsize>1) {
      printfone("TRAJOUT: remdtraj does not work correctly for multiprocessor writes.\n"); 
      printfone("         Temperatures will NOT be written to output trajectory.\n");
      info->isREMDTRAJ=0;
    }

    /*
     *  TODO: Integrate LES options into processing aka strip so that actions can act on
     *  LES average or subset trajectories
     */


    /* 
     *  check if there are les options
     */

    info->les_action = LES_NONE;
    if (argumentStackContains( &argumentStack, "les" ) )
    {
        info->nlescopy = ( parm->nlestyp == 1 ) ? (int)( parm->lesfac[0] + 0.1 ) : (int)( parm->lesfac[3] + 0.1 );

        if( argumentStackContains( &argumentStack, "split" ) )
        {
            info->les_action = LES_SPLIT;
	    info->les_status = LES_READY;
        }
        else if( argumentStackContains( &argumentStack, "average" ) )
        {
            info->les_action = LES_AVERAGE;
	    info->les_status = LES_READY;
        }
        else
        {
            error( "setup_les_output", "unknown les action" );
        } 

    }
    
       /*
        *  check to see if a format other than amber trajectory is wanted
        */
    if (argumentStackContains( &argumentStack, "pdb" ))
      info->type = COORD_PDB;
    else if (argumentStackContains( &argumentStack, "restart" ))
      info->type = COORD_AMBER_RESTART;
    else if (argumentStackContains( &argumentStack, "restrt" ))
      info->type = COORD_AMBER_RESTART;
    else if (argumentStackContains( &argumentStack, "rest" ))
      info->type = COORD_AMBER_RESTART;
    else if (argumentStackContains( &argumentStack, "binpos" ))
      info->type = COORD_BINPOS;
    else if (argumentStackContains( &argumentStack, "charmm" ))
      info->type = COORD_CHARMM_TRAJECTORY;
#ifdef BINTRAJ
    else if (argumentStackContains( &argumentStack, "netcdf" ))
      info->type = COORD_AMBER_NETCDF;
    else if (argumentStackContains( &argumentStack, "cdf" ))
      info->type = COORD_AMBER_NETCDF;
#else
    else if (argumentStackContains( &argumentStack, "netcdf" )) {
      info->type = COORD_AMBER_TRAJECTORY;
      printfone("trajout: NetCDF support is not compiled into this version (Add -DBINTRAJ)\n");
      printfone("         defaulting to AMBER trajectory format...\n");
    } else if (argumentStackContains( &argumentStack, "cdf" )) {
      info->type = COORD_AMBER_TRAJECTORY;
      printfone("trajout: NetCDF support is not compiled into this version (Add -DBINTRAJ)\n");
      printfone("         defaulting to AMBER trajectory format...\n");
    }
#endif
    else
      info->type = COORD_AMBER_TRAJECTORY;

      /*
       *  set up program/version information
       */
    switch (info->type) {
    case COORD_PDB:
      info->title = argumentStackKeyToString( &argumentStack, "title", "PDB file generated by ptraj" );
      break;
    case COORD_AMBER_NETCDF:
      /*
      info->title = argumentStackKeyToString( &argumentStack, "title", "NetCDF trajectory generated by ptraj" );
      */
      info->title = argumentStackKeyToString( &argumentStack, "title", "" );
      info->isNetcdf=1;
      break;

    default:
      info->title = argumentStackKeyToString( &argumentStack, "title", "trajectory generated by ptraj" );
    }
    info->application = argumentStackKeyToString( &argumentStack, "application", "AMBER" );
    /*
    info->program =     argumentStackKeyToString( &argumentStack, "program", "ptraj" );
    */
    info->program =     argumentStackKeyToString( &argumentStack, "program", "sander" );
    /*
    info->version =     argumentStackKeyToString( &argumentStack, "version", PTRAJ_VERSION_STRING );
    */
    info->version =     argumentStackKeyToString( &argumentStack, "version", "9.0" );

#ifdef MPI
    /* If output is AMBER trajectory and using MPI, set isMPI to 1 */
    // Only use MPI file ops with more than 1 thread.
    if (info->type == COORD_AMBER_TRAJECTORY && worldsize>1)
      info->isMPI = 1;
#endif      

       /*
        *  check to see if we want charges/radii dumped to pdb
        */
    if (info->type == COORD_PDB) {

      info->option1 = 0;
      if (argumentStackContains( &argumentStack, "dumpq" ))
	info->option1 = 1;
      else if (argumentStackContains( &argumentStack, "parse" ))
	info->option1 = 2;
      else if (argumentStackContains( &argumentStack, "dumpr*" ))
	info->option1 = 3;

      if (argumentStackContains( &argumentStack, "nowrap" ))
	info->option2 = 1;
    }

       /*
        *  check to see if other CHARMM related information is present
        */
    if (info->type == COORD_CHARMM_TRAJECTORY) {

      /*
       *  search through the list of input files to setup the CHARMM information
       *  structure defaults; if none is present, make it up and/or modify according
       *  to what the user specifies...
       */
      charmmTraj = (charmmTrajectoryInfo *) safe_malloc(sizeof(charmmTrajectoryInfo));
      INITIALIZE_charmmTrajectoryInfo(charmmTraj);

      ctrj = NULL;
      for (sp = transformFileStack; ctrj == NULL && sp != NULL; sp = sp->next) {
	infiles = (coordinateInfo *) sp->entry;
	if (infiles->type == COORD_CHARMM_TRAJECTORY) 
	  ctrj = (charmmTrajectoryInfo *) infiles->info;
      }

      if (ctrj != NULL) {
	charmmTraj->byteorder = ctrj->byteorder;
	charmmTraj->magic = ctrj->magic;
	for (i=0; i < 20; i++)
	  charmmTraj->icntrl[i] = ctrj->icntrl[i];
	charmmTraj->ntitle = ctrj->ntitle;
	for (sp = ctrj->titleStack; sp != NULL; sp = sp->next) {
	  title = (char *) sp->entry;
	  titlenew = (char *) safe_malloc(sizeof(char) * (strlen(title) + 1));
	  strcpy(titlenew, title);
	  pushBottomStack(&charmmTraj->titleStack, (void *) title);
	}
	charmmTraj->natrec = ctrj->natrec;
	charmmTraj->nfreat = ctrj->nfreat;
	if (ctrj->nfreat != ctrj->natrec) {
	  charmmTraj->freeat = (int *) safe_malloc(sizeof(int) * ctrj->nfreat);
	  for (i=0; i < ctrj->nfreat; i++) 
	    charmmTraj->freeat[i] = ctrj->freeat[i];
	}
	for (i=0; i<6; i++)
	  charmmTraj->xtlabc[i] = ctrj->xtlabc[i];

      } else {
	charmmTraj->byteorder = 0;
	charmmTraj->magic.c[0] = 'C';
	charmmTraj->magic.c[1] = 'O';
	charmmTraj->magic.c[2] = 'R';
	charmmTraj->magic.c[3] = 'D';
	for (i=0; i<20; i++)
	  charmmTraj->icntrl[i] = 0;
	charmmTraj->icntrl[19] = 26;

	if (state->IFBOX) {
	  charmmTraj->icntrl[10] = 1;  /* QCRYS */
	}
      }

      if (info->isBox == 0) charmmTraj->icntrl[10] = 0;


      if (argumentStackContains( &argumentStack, "big" ))
	charmmTraj->byteorder = 0;
      else if (argumentStackContains( &argumentStack, "little" ))
	charmmTraj->byteorder = 1;
	
      info->info = (void *) charmmTraj;
    }

    safe_free(filename);
    globalOutInfo = info;
    return;
  }
}


/*
 *  ptrajSetup(): This routine is called for every trigger that is related
 *  to coordinate processing (i.e. not those commands that are I/O 
 *  related, such as trajin, trajout or reference and not those that
 *  are involved with postprocessing any acculated data).  This creates
 *  the transformActionStack stack of "actions" to be performed and 
 *  is called by dispatchToken() upon receipt of the appropriate trigger.
 *  Most of the actual setup of the action function is performed by the
 *  actual action function itself in the PTRAJ_SETUP mode.  See the
 *  detailed comments in actions.c for more information.
 */

#undef  ROUTINE
#define ROUTINE "ptrajSetup()"

   void
ptrajSetup(stackType *argumentStack, actionType type)
{
  actionInformation *action;
  int ierr;

     /*
      *  Allocate and initialize a actionInformation structure
      */

  action = (actionInformation *)
    safe_malloc(sizeof(actionInformation));
  INITIALIZE_actionInformation(action);

     /*
      *  Place a copy of the current state into the action
      */

  action->state = ptrajCopyState(ptrajCurrentState());

     /*
      *  Set the action type
      */

  action->type = type;

     /*
      *  Set the action function
      */

  switch ( type ) {

  case TRANSFORM_ANGLE:

    action->type = TRANSFORM_ANGLE;
    action->fxn  = (actionFunction) transformAngle;
    break;

  case TRANSFORM_ATOMICFLUCT:

    action->type = TRANSFORM_ATOMICFLUCT;
    action->fxn  = (actionFunction) transformAtomicFluct;
    break;

  case TRANSFORM_ATOMICFLUCT3D:

    action->type = TRANSFORM_ATOMICFLUCT3D;
    action->fxn  = (actionFunction) transformAtomicFluct3D;
    break;

  case TRANSFORM_AVERAGE:

    action->type = TRANSFORM_AVERAGE;
    action->fxn  = (actionFunction) transformAverage;
    break;

  case TRANSFORM_CENTER:

    action->type = TRANSFORM_CENTER;
    action->fxn  = (actionFunction) transformCenter;
    break;

  case TRANSFORM_CHECKOVERLAP:

    action->type = TRANSFORM_CHECKOVERLAP;
    action->fxn  = (actionFunction) transformCheckOverlap;
    break;

  case TRANSFORM_CLOSESTWATERS:

    action->type = TRANSFORM_CLOSESTWATERS;
    action->fxn  = (actionFunction) transformClosestWaters;
    break;

  case TRANSFORM_CLUSTER:

    action->type = TRANSFORM_CLUSTER;
    action->fxn  = (actionFunction) transformCluster;
    break;

  case TRANSFORM_CLUSTERATTRIBUTE:

    action->type = TRANSFORM_CLUSTERATTRIBUTE;
    action->fxn  = (actionFunction) transformClusterAttribute;
    break;

  case TRANSFORM_CORRELATION:

    action->type = TRANSFORM_CORRELATION;
    action->fxn  = (actionFunction) transformCorr;
    break;

  case TRANSFORM_CONTACTS:

    action->type = TRANSFORM_CONTACTS;
    action->fxn  = (actionFunction) transformContacts;
    break;

  case TRANSFORM_DIHEDRAL:

    action->type = TRANSFORM_DIHEDRAL;
    action->fxn  = (actionFunction) transformDihedral; 
    break;

  case TRANSFORM_DIHEDRALCLUSTER:

    action->type = TRANSFORM_DIHEDRALCLUSTER;
    action->fxn  = (actionFunction) transformDihedralCluster; 
    break;

  case TRANSFORM_DIFFUSION:

    action->type = TRANSFORM_DIFFUSION;
    action->fxn  = (actionFunction) transformDiffusion; 
    break;

  case TRANSFORM_DIPOLE:

    action->type = TRANSFORM_DIPOLE;
    action->fxn  = (actionFunction) transformDipole;
    break;

  case TRANSFORM_DISTANCE:

    action->type = TRANSFORM_DISTANCE;
    action->fxn  = (actionFunction) transformDistance;
    break;

  case TRANSFORM_DNAIONTRACKER:

    action->type = TRANSFORM_DNAIONTRACKER;
    action->fxn  = (actionFunction) transformDNAiontracker;
    break;

  case TRANSFORM_ECHO:

    action->type = TRANSFORM_ECHO;
    action->fxn = (actionFunction) transformEcho;
    break;

  case TRANSFORM_ENERGY:

    action->type = TRANSFORM_ENERGY;
    action->fxn  = (actionFunction) transformEnergy;
    break;

  case TRANSFORM_GRID:

    action->type = TRANSFORM_GRID;
    action->fxn  = (actionFunction) transformGrid;
    break;

  case TRANSFORM_HBOND:

    action->type = TRANSFORM_HBOND;
    action->fxn  = (actionFunction) transformHBond;
    break;

  case TRANSFORM_IMAGE:

    action->type = TRANSFORM_IMAGE;
    action->fxn  = (actionFunction) transformImage;
    break;

  case TRANSFORM_MATRIX:

    action->type = TRANSFORM_MATRIX;
    action->fxn  = (actionFunction) transformMatrix;
    break;

  case TRANSFORM_PRINCIPAL:

    action->type = TRANSFORM_PRINCIPAL;
    action->fxn  = (actionFunction) transformPrincipal;
    break;

  case TRANSFORM_PROJECTION:

    action->type = TRANSFORM_PROJECTION;
    action->fxn  = (actionFunction) transformProjection;
    break;

  case TRANSFORM_PUCKER:

    action->type = TRANSFORM_PUCKER;
    action->fxn  = (actionFunction) transformPucker;
    break;

  case TRANSFORM_RADIAL:

    action->type = TRANSFORM_RADIAL;
    action->fxn  = (actionFunction) transformRadial;
    break;

  case TRANSFORM_RADIUSOFGYRATION:

    action->type = TRANSFORM_RADIUSOFGYRATION;
    action->fxn  = (actionFunction) transformRadiusOfGyration;
    break;

  case TRANSFORM_RANDOMIZEIONS:

    action->type = TRANSFORM_RANDOMIZEIONS;
    action->fxn  = (actionFunction) transformRandomizeIons;
    break;

  case TRANSFORM_RMS:
      
    action->type = TRANSFORM_RMS;
    action->fxn  = (actionFunction) transformRMS;
    break;

  case TRANSFORM_RUNNINGAVERAGE:
      
    action->type = TRANSFORM_RUNNINGAVERAGE;
    action->fxn  = (actionFunction) transformRunningAverage;
    break;

  case TRANSFORM_SCALE:
      
    action->type = TRANSFORM_SCALE;
    action->fxn  = (actionFunction) transformScale;
    break;

  case TRANSFORM_SECONDARYSTRUCT:

    action->type = TRANSFORM_SECONDARYSTRUCT;
    action->fxn  = (actionFunction) transformSecondaryStruct;
    break;

  case TRANSFORM_STRIP:
      
    action->type = TRANSFORM_STRIP;
    action->fxn  = (actionFunction) transformStrip;
    break;

  case TRANSFORM_TRANSLATE:

    action->type = TRANSFORM_TRANSLATE;
    action->fxn  = (actionFunction) transformTranslate;
    break;

  case TRANSFORM_TRUNCOCT:

    action->type = TRANSFORM_TRUNCOCT;
    action->fxn  = (actionFunction) transformTruncOct;
    break;

  case TRANSFORM_UNWRAP:

    action->type = TRANSFORM_UNWRAP;
    action->fxn  = (actionFunction) transformUnwrap;
    break;

  case TRANSFORM_VECTOR:

    action->type = TRANSFORM_VECTOR;
    action->fxn  = (actionFunction) transformVector;
    break;

  case TRANSFORM_WATERSHELL:

    action->type = TRANSFORM_WATERSHELL;
    action->fxn  = (actionFunction) transformWatershell;
    break;

  case TRANSFORM_2DRMS:

    action->type = TRANSFORM_2DRMS;
    action->fxn  = (actionFunction) transform2dRMS;
    break;

  case TRANSFORM_TRANSFORM:
  case TRANSFORM_NOOP:

    action->type = type;
    action->fxn = NULL;
    break;

  default:

    printfone("%s: Attempting to setup an unknown action type %i\n", ROUTINE, type);
    error(ROUTINE, "There is no way you should be here!  Terminating...\n");

  }

     /*
      *  Parse the arguments.  This is done by the action->fxn in the
      *  PTRAJ_SETUP mode.  The argumentStack is placed into the
      *  complex argument 1 slot for this.
      */

  ierr = 0;
  if (action->type != TRANSFORM_TRANSFORM && 
      action->type != TRANSFORM_NOOP) {

    action->carg1 = (void *) &argumentStack;
    ierr = action->fxn(action, NULL, NULL, NULL, NULL, PTRAJ_SETUP);

  }

     /*
      *  If the setup fails, -1 is returned and therefore this action
      *  should not be placed on the action stack and the associated
      *  memory should be freed
      */

  if (ierr < 0) {
    safe_free(action->state);
    action->state= NULL;
    safe_free(action);
    action = NULL;
  } else {

     /*
      *  Place the now setup action structure onto the transformActionStack
      */

    pushBottomStack( &transformActionStack, (void *) action );
  
  }


}

int
checkDivisibility(int *stop, int start, int offset) {
  /*
   *  Check to see if evenly divisible, if not change stop so it is.
   *  Return 0 if stop has been changed, else return 1
   */
#ifdef MPI
  if ( (*stop - start + 1) % (offset * worldsize) != 0) {
    *stop = ((*stop - start + 1) / (offset * worldsize)) * (offset * worldsize) + (start - 1);
    return 0;
  }
#endif    
    
  return 1;
}

#undef  ROUTINE
#define ROUTINE "ptrajSetupAnalyze()"

   void
ptrajSetupAnalyze(stackType *argumentStack, actionType type)
{
  analyzeInformation *analyze;
  char *buffer;
  int ierr;

     /*
      *  Make sure that this is indeed an "analyze" action (which really isn't
      *  necessary)
      */
  if (type != TRANSFORM_ANALYZE) {
    printfone("Error in ptrajSetupAnalyze(): Called with the wrong type!\n");
    printfone("Ignoring this command...\n");
    return;
  }

     /*
      *  Grab the first argument off the argument stack.  This is the "trigger" for
      *  the analyze function
      */
  buffer = getArgumentStringLower(&argumentStack, NULL);
  if (buffer == NULL) {
    printfone("ptrajSetupAnalyze(): No command passed to analyze, ignoring...\n");
    return;
  }

     /*
      *  Allocate and initialize a analyzeInformation structure and set the type
      */
  analyze = (analyzeInformation *)
    safe_malloc(sizeof(analyzeInformation));
  INITIALIZE_analyzeInformation(analyze);

     /*
      *  search for a match to the trigger (stored in buffer)
      */

  if (strncmp(buffer, "correlationcoe", 14) == 0) {
    analyze->type = ANALYZE_CORRELATIONCOEFFICIENT;
    analyze->fxn  = (analyzeFunction) analyzeCorrelationCoefficient;

  } else if (strncmp(buffer, "crank", 5) == 0) {
    analyze->type = ANALYZE_CRANKSHAFT;
    analyze->fxn  = (analyzeFunction) analyzeCrankshaft;

  } else if (strncmp(buffer, "hbond", 5) == 0) {
    analyze->type = ANALYZE_HBOND;
    analyze->fxn  = (analyzeFunction) analyzeHBond;

  } else if (strncmp(buffer, "matrix", 6) == 0) {
    analyze->type = ANALYZE_MATRIX;
    analyze->fxn  = (analyzeFunction) analyzeMatrix;

  } else if (strncmp(buffer, "modes", 5) == 0) {
    analyze->type = ANALYZE_MODES;
    analyze->fxn  = (analyzeFunction) analyzeModes;

  } else if (strncmp(buffer, "set", 3) == 0) {
    analyze->type = ANALYZE_SET;
    analyze->fxn  = (analyzeFunction) analyzeSet;

  } else if (strncmp(buffer, "stat", 4) == 0) {
    analyze->type = ANALYZE_STATISTICS;
    analyze->fxn  = (analyzeFunction) analyzeStatistics;

  } else if (strncmp(buffer, "timecorr", 8) == 0) {
    analyze->type = ANALYZE_TIMECORR;
    analyze->fxn  = (analyzeFunction) analyzeTimecorr;

  } else if (strncmp(buffer, "test", 4) == 0) {
    analyze->type = ANALYZE_TEST;
    analyze->fxn  = (analyzeFunction) analyzeTest;

  } else {

    printfone("WARNING in ptrajSetupAnalyze(): unknown analyze type %i\n", type);
    safe_free(analyze);
    return;

  }

     /*
      *  Parse the arguments.  This is done by the analyze->fxn itself in the
      *  PTRAJ_SETUP mode.  The argumentStack is placed into the
      *  complex argument 1 slot for this.
      */

  ierr = 0;
  if (analyze->type != ANALYZE_NOOP) {

    analyze->carg1 = (void *) &argumentStack;
    ierr = analyze->fxn(analyze, scalarStack, PTRAJ_SETUP);

  }

     /*
      *  If the setup fails, -1 is returned and therefore this action
      *  should not be placed on the action stack and the associated
      *  memory should be freed
      */

  if (ierr < 0) {
    safe_free(analyze);
    analyze = NULL;
  } else {

     /*
      *  Place the now setup analyze structure onto the transformAnalyzeStack
      */

    pushBottomStack( &transformAnalyzeStack, (void *) analyze );
  
  }


}




   void
ptrajProcessTrajectoryFiles(double* X, double* Y, double* Z, 
			    ptrajMode ActionMode, ptrajState *startingState,
			    double box[6], int boxfixed[6], 
			    int outputTrajectory, int* set, int* processed)
{
  /*
   *  --------- MAIN LOOP FOR COORDINATE PROCESSING -------------
   *
   *  loop over each of the files representing the coordinates.
   *
   *  "set"        -- the global counter (for this thread) over all sets, all files
   *  "local_set"  -- the counter over sets in each individual file
   *  "output_set" -- the output frame counter
   */

  int firstOutput;
  int readCoordinates;
  int processCoordinates;
  int local_set;
  int output_set;
  stackType *sp;
  coordinateInfo *currentCoordinateInfo;
  int suppressProcessing;
  int i,per;
  double boxnew[6];
  actionInformation *action;
  ptrajState *currentState;
  stackType* actionStackTemp;
  int startSet, stopSet, offset, frameskip;

  output_set = 0;
  firstOutput = 1;
  for (sp = transformFileStack; sp != NULL; sp = sp->next) {

    currentCoordinateInfo = (coordinateInfo *) sp->entry;
    per = 0;

    /*
     *  ------------- PREPROCESS COORDINATE FILES-----------------
     */

    /*
     *  open up the file and preprocess the coordinates
     */
    
#ifdef MPI
    t1 = MPI_Wtime();
#else
    t1 = clock();
#endif    
    
    if (ptrajPreprocessInputCoordinates(currentCoordinateInfo)) 
      {
	/* Failure!  Bail out now. */
	ptrajCleanup();
	return;
      }
    
#ifdef MPI
    t2 = MPI_Wtime();
    inputTime += (t2 - t1);
#else
    t2 = clock();
    inputTime += (t2 - t1) / CLOCKS_PER_SEC;
#endif    

    startSet = currentCoordinateInfo->start;
    stopSet = currentCoordinateInfo->stop;
    offset = currentCoordinateInfo->offset;
    /*
     *  -------- READ IN, PROCESS, AND OUTPUT COORDINATES --------
     *
     *  The following two variables control whether to continue
     *  reading configurations from this file or not and also
     *  whether the coordinates should be processed after this
     *  read...
     *
     *  "readCoordinates"    > 0 if there are more coordinate sets 
     *                           to read in this file
     *
     *  "processCoordinates" > 0 if this coordinate set should be
     *                           processed
     */

    readCoordinates = 1;
    processCoordinates = 1;
    currentState = startingState; /* init */

    /*
     *  If the file is seekable, set local_set to the starting frame
     *  and skip by offset, otherwise start from frame 1 and read sequentially
     */

    if (currentCoordinateInfo->seekable || currentCoordinateInfo->isNetcdf) {
      local_set = (startSet + worldrank * offset);
      frameskip = (offset * worldsize);
    } else {
      local_set = 1;
      frameskip = 1;
    }
    if (output_set==0) output_set=local_set;

    while (readCoordinates) {
      /*      if (local_set > stopSet)
       *	break;
       *
       *  read in the current file of coordinates, a single set at a time.
       */
      for (i=0; i<6; i++)
	boxnew[i] = box[i];
      
#ifdef MPI
      t1 = MPI_Wtime();
#else
      t1 = clock();
#endif    
      
      ptrajProcessInputCoordinates(currentCoordinateInfo, startingState, X, Y, Z, boxnew,
				   local_set, &readCoordinates, &processCoordinates);
#ifdef MPI
      t2 = MPI_Wtime();
      inputTime += (t2 - t1);
#else
      t2 = clock();
      inputTime += (t2 - t1) / CLOCKS_PER_SEC;
#endif    
      
      for (i=0; i<6; i++)
	if (boxfixed[i] == 0)
	  box[i] = boxnew[i];

      
      /*
       *  process coordinates if necessary
       */

      if (processCoordinates) {

	/*
	 *  check to see if this snapshot is within bounds/offset
	 */

	if ((local_set >= startSet &&
	     (local_set <= stopSet ||
	      stopSet == -1)) &&
	    ((offset == 1) ||
	     ((local_set - startSet) %
	      offset == 0))) {

	  /*
	   *  TRAVERSE THE ACTION STACK to perform each
	   *  action on each coordinate set.  Note that a 
	   *  particular action can suppress processing of
	   *  further actions on the stack and prevent output
	   *  by setting the suppressProcessing flag in the
	   *  actionInformation structure.  This is useful 
	   *  when various sets are to be accumulated (for
	   *  example when calculating running averages) prior
	   *  to output and further processing...
	   */

	  suppressProcessing = 0;

#ifdef MPI
	  t1 = MPI_Wtime();
#else
	  t1 = clock();
#endif    

	  if (transformActionStack) {

	    for (actionStackTemp = transformActionStack;
		 actionStackTemp != NULL;
		 actionStackTemp = actionStackTemp->next) {
        
	      action = (actionInformation *) actionStackTemp->entry;
	      if (action->type != TRANSFORM_NOOP &&
		  action->type != TRANSFORM_TRANSFORM &&
		  suppressProcessing == 0) {

		for (i=0; i<6; i++)
		  boxnew[i] = box[i]; /* protect box coordinates */

		action->fxn(action, X, Y, Z, boxnew, ActionMode);

		for (i=0; i<6; i++)
		  if (boxfixed[i] == 0)
		    box[i] = boxnew[i];


		/*
		 *  update the current state
		 */
		currentState = action->state;
            
		/*
		 *  check if any of the actions have suppressed output
		 */

		if (suppressProcessing == 0)
		  suppressProcessing = action->suppressProcessing;
	      }
	    }
	  }

#ifdef MPI
	  t2 = MPI_Wtime();
	  actionTime += (t2 - t1);
#else
	  t2 = clock();
	  actionTime += (t2 - t1) / CLOCKS_PER_SEC;
#endif    

	  /*
	   *  perform output as necessary
	   */

#ifdef MPI
	  t1 = MPI_Wtime();
#else
	  t1 = clock();
#endif    

	  if (outputTrajectory && suppressProcessing == 0) {

	    ptrajOutputCoordinates(globalOutInfo, currentState, output_set, globalOutInfo->append, firstOutput, 0,
				   currentState->atoms, X, Y, Z, box);
	    firstOutput = 0;

	  }

#ifdef MPI
	  t2 = MPI_Wtime();
	  outputTime += (t2 - t1);
#else
	  t2 = clock();
	  outputTime += (t2 - t1) / CLOCKS_PER_SEC;
#endif    

	  (*processed)++;
    
	  if (stopSet == -1) {
	    if (local_set % 50 == 0 || local_set == 1)
	      fprintf(stdout, "\nSet %6i ", local_set);
	    if (local_set % 50 != 0)
	      fprintf(stdout, ".");
	  }
	} else if (stopSet == -1) {
	  if (local_set % 50 == 0 || local_set == 1)
	    fprintf(stdout, "\nSet %6i ", local_set);
	  if (local_set % 50 != 0)
	    fprintf(stdout, " ");
	}
	if (stopSet != -1) {
	  if (currentCoordinateInfo->seekable) {
	    if (startSet != stopSet && ((local_set - startSet) * 100) / (stopSet - startSet) != per) {
	      per = ((local_set - startSet) * 100) / (stopSet - startSet);
	      if (per % 25 == 0 || per == 1) {
		printfone(" %d%% ", per);
	      } else if (per % 2 == 0) {
		printfone(".");
	      }
	    }
	  } else {
	    if ((local_set * 100) / stopSet != per) {
	      per = (local_set * 100) / stopSet;
	      if (per % 25 == 0 || per == 1) {
		printfone(" %d%% ", per);
	      } else if (per % 2 == 0) {
		printfone(".");
	      }
	    }
	  }
	}
    
	fflush(stdout);
	*set += 1;
        output_set += frameskip;
      }   /* IF (processCoordinates) */
      local_set += frameskip;
    }     /* WHILE (readCoordinates) */
    if (per != 100)
      printfone(" 100%%");
    printfone("\n");
  }     /* FOR over transformFileStack */
}

#undef  ROUTINE
#define ROUTINE "ptraj()"

   void
ptraj(char *filenamep)
{
  FILE *infile, *bfile;
  char buffer[BUFFER_SIZE];
  char *bufferp;
  coordinateInfo *currentCoordinateInfo;
  double *X = NULL, *Y = NULL, *Z = NULL;
  double box[6], boxnew[6];
  int boxfixed[6];
  actionInformation *action;
  analyzeInformation *analyze;
  int set;
  int local_set;
  int processed;
  int i,j;
  int suppressProcessing;
  stackType *actionStackTemp = NULL;
  int outputTrajectory = 0;
  stackType *sp, *argumentStack;
  ptrajState *startingState, *currentState, **statep;
  char *continuation;
  int totalFrames;

  int readCoordinates;
  int processCoordinates;
  int firstOutput;
  int SecondPassRequired;

  checkInputTime = inputTime = outputTime = actionTime = t1 = t2 = 0.0;

    /*
     *  --------------- INPUT FILE PROCESSING --------------------
     *
     *  if an input file was specified, open it up, 
     *  else use standard input...
     */
  if ( filenamep == NULL || 
       strcmp(filenamep, "") == 0 || 
       strcmp(filenamep, "stdin") == 0 ) {

    printfone("\nPTRAJ: Processing input from \"STDIN\" ...\n");
    infile = stdin;

  } else if ( openFile(&infile, filenamep, "r") == 0 ) {
    printfone("WARNING in %s: Could not open input file (%s), exiting\n", ROUTINE,
	      filenamep);
    ptrajCleanup();
    return;

  } else {
    printfone("\nPTRAJ: Processing input from file %s\n", filenamep);
  }


    /*
     *  ----------------- SETUP INITIAL STATE --------------------
     *
     *  this gives a snapshot of the current "state" based on the 
     *  appropriate parameter/topology information (i.e. AMBER prmtop)
     *  upon entry.  Note this assumes that the GLOBAL ptrajState
     *  information was previously set in the main routine (main.c)
     *  via a call to ptrajInitializeState().
     */

  statep = ptrajCurrentState();
  startingState = *statep;
  currentState = startingState;

    /*
     *  --------------- PROCESS THE INPUT FILE -------------------
     *
     *  Read the input file, line by line, using the "dispatchToken()"
     *  routine (dispatch.c) to find a match for each command, ignoring
     *  comments in the input file.  Each "command" typed by the user 
     *  has an associated routine that is called for setup.  Currently,
     *  this is ptrajSetup() for "actions" and ptrajSetupIO for input/
     *  output functions.  Note that ptrajSetup() will actually call the
     *  "action" routine to perform the setup and parse arguments.  Note
     *  also that the global state pointer may be altered!!!
     *
     *  Lines of text are processed until EOF is encountered.
     */

  argumentStack = NULL;
  while ( (bufferp = fgets(buffer, BUFFER_SIZE, infile)) != NULL ) {

    continuation = bufferp;
    while ( continuation != NULL ) {
      if (strlen(buffer) >= BUFFER_SIZE)
	continuation = NULL;
      else {
	continuation = strrchr(buffer, '\\');
   
	if (continuation)
	  continuation = fgets(continuation, (BUFFER_SIZE - strlen(buffer) - 1), infile);
      }
    }

    skipWhitespace(bufferp);

    // Exit on keyword 'go'
    if (strncmp(bufferp,"go\n",3) == 0) break;

      /*
       *  skip blank lines and/or comments denoted by "#" or "%"
       */
    if (bufferp[0] != (char) 0 && 
	bufferp[0] != '#' &&
	bufferp[0] != '%') {  

      printfone("\nPTRAJ: %s", buffer);
      dispatchToken((Token *) &ptrajTokenlist, argumentStack, (char *) buffer);
    }
  }

  if (infile != stdin)
    safe_fclose(infile);

    /*
     *  ----------------- ERROR CHECKING ---------------------
     */

  if (transformFileStack == NULL) {
    printfone("WARNING in %s: No input trajectories specified (trajin), aborting...\n",ROUTINE);
    ptrajCleanup();
    return;
  }

  if (globalOutInfo == NULL) {
    printfone("[No output trajectory specified (trajout)]\n");
  } else {
    outputTrajectory = 1;
  }

    /*
     *  set default box information.  This will allow, at setup time, specification
     *  of the current box sizes (if none were specified such as can happen when
     *  reading CHARMM PSF files) and also the option to FIX box sizes.
     */
  if (currentState->IFBOX) {
    for (i=0; i<6; i++) {
      box[i] = currentState->box[i];
      boxfixed[i] = currentState->boxfixed[i];
    }
  } else {
    box[0] = 0.0;
    box[1] = 0.0;
    box[2] = 0.0;
    box[3] = 90.0;
    box[4] = 90.0;
    box[5] = 90.0;
    boxfixed[0] = 0; boxfixed[1] = 0; boxfixed[2] = 0;
    boxfixed[3] = 0; boxfixed[4] = 0; boxfixed[5] = 0;
  }

    /*
     * --------- CHECK HOW MANY FRAMES WILL BE PROCESSED -------
     */

  /*
   *  TO DO: FIX THIS AS IT IS BROKEN
   */

  startingState->maxFrames = 0;
  totalFrames=0;
  for (sp = transformFileStack; sp != NULL; sp = sp->next) {
    currentCoordinateInfo = (coordinateInfo *) sp->entry;

    fprintf(stdout,"  %s: %i frames.\n",currentCoordinateInfo->filename,currentCoordinateInfo->Nframes);

    startingState->maxFrames += (currentCoordinateInfo->stop - 
				 currentCoordinateInfo->start) /
                                 currentCoordinateInfo->offset + 1;
    // In case we dont know how many frames will be processed, use total number of frames
    totalFrames += currentCoordinateInfo->Nframes;
  }

  printfone("\nPTRAJ: Successfully read the input file.\n");
  if (startingState->maxFrames == -1) {
    printfone("       Coordinate processing will occur until EOF (unknown number of frames).\n");
    // Use total frames for action memory alloc.
    startingState->maxFrames = totalFrames;
  } else
    printfone("       Coordinate processing will occur on %i frames.\n", 
	      startingState->maxFrames);
  printfone("       Summary of I/O and actions follows:\n\n");
  
  
  /*
   * ----- PRINT A SUMMARY OF THE FILES IN/OUT AND ACTIONS ------
   */
  
  printfone("INPUT COORDINATE FILES\n");
  
  printStack( &transformFileStack, printCoordinateInfo, NULL );

  if ( globalOutInfo == NULL ) {
    printfone("\nNO OUTPUT COORDINATE FILE WAS SPECIFIED\n");
  } else {
    printfone("\nOUTPUT COORDINATE FILE\n");
    printCoordinateInfo( (void *) globalOutInfo );
  }

  if (transformReferenceStack != NULL) {
    printfone("\nREFERENCE FILE\n");
    printStack( &transformReferenceStack, printCoordinateInfo, NULL );
  }

  if ( transformActionStack == NULL ) {

    printfone("\nNO ACTIONS WERE SPECIFIED\n");

  } else {

    printfone("\nACTIONS\n");
    i = 1;
    for (actionStackTemp = transformActionStack;
	 actionStackTemp != NULL;
	 actionStackTemp = actionStackTemp->next) {

      action = (actionInformation *) actionStackTemp->entry;
         /*
          *  With each action's state variable local, it is necessary
          *  to update each with the current states maxFrames value!!!
          */
      action->state->maxFrames = startingState->maxFrames;
         /*
          *  call the action function with the mode PTRAJ_STATUS
          */
      if (action->type != TRANSFORM_NOOP  &&
	  action->type != TRANSFORM_TRANSFORM) {
	printfone("%3i>", i);
	for (j=0; j < 6; j++)
	  boxnew[j] = box[j]; /* protect the current box info */
	action->fxn(action, X, Y, Z, boxnew, PTRAJ_STATUS);
	i++;
      }
    }
    printfone("\n");
  }

  if (transformAnalyzeStack != NULL) {
    printfone("\nANALYZE\n");
    i = 1;
    for (actionStackTemp = transformAnalyzeStack;
	 actionStackTemp != NULL;
	 actionStackTemp = actionStackTemp->next) {

      analyze = (analyzeInformation *) actionStackTemp->entry;
      if (analyze->type != ANALYZE_NOOP) {
	printfone("%3i>", i);
	analyze->fxn(analyze, scalarStack, PTRAJ_STATUS);
	i++;
      }
    }
    printfone("\n");
  }

    /*
     *  ---------- ALLOCATE SPACE FOR COORDINATES ------------
     *
     *  perform and initial setup necessary prior to reading in
     *  the coordinates
     */

  X = (double *) safe_malloc(sizeof(double) * startingState->atoms);
  Y = (double *) safe_malloc(sizeof(double) * startingState->atoms);
  Z = (double *) safe_malloc(sizeof(double) * startingState->atoms);

    /*
     *  Perform one pass over each trajectory file: 
     */

  processed = 0;
  set = 1;
  ptrajProcessTrajectoryFiles(X,Y,Z,PTRAJ_ACTION,startingState,box,boxfixed,outputTrajectory,&set,&processed);

#ifdef MPI
  /*
   * Set up a barrier to sync threads
   */

  MPI_Barrier(MPI_COMM_WORLD);
#endif

    /*
     *  Tell each action that the first pass is complete 
     */

  for (actionStackTemp = transformActionStack; actionStackTemp != NULL; actionStackTemp = actionStackTemp->next) 
  {
      action = (actionInformation *) actionStackTemp->entry;
      action->fxn(action, NULL, NULL, NULL, NULL, PTRAJ_FIRSTPASS);
  }

    /*
     *  If necessary, perform a second pass over each trajectory file: 
     */

  SecondPassRequired = 0;
  for (actionStackTemp = transformActionStack; actionStackTemp != NULL; actionStackTemp = actionStackTemp->next) 
  {
    action = (actionInformation *) actionStackTemp->entry;
    if (action->performSecondPass)
    {
        SecondPassRequired = 1;
        break;
    }
  }
  if (SecondPassRequired)
  {
    set = 1;
    processed = 0;
    outputTrajectory = 0;
    ptrajProcessTrajectoryFiles(X,Y,Z,PTRAJ_SECONDPASS,
        startingState,box,boxfixed,outputTrajectory,&set,&processed);
  }


    /*
     *  -------------- FINAL POSTPROCESSING -----------------
     */

  if (outputTrajectory) {
    ptrajOutputCoordinates(globalOutInfo, currentState, set, globalOutInfo->append, 0, 1, 
			   currentState->atoms, NULL, NULL, NULL, NULL);
  }

    /*
     *  -------------- DUMP ACCUMULATED DATA -------------------
     */

#ifdef MPI  
  MPI_Barrier(MPI_COMM_WORLD);
  printf("\nRank %d successful. Read in %i sets and processed %i sets.\n",
	  worldrank, set-1, processed);
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
#else
  printf("\n\nPTRAJ: Successfully read in %i sets and processed %i sets.\n",
	  set-1, processed);
#endif
  printfone("\nDumping accumulated results (if any)\n\n");
  for (actionStackTemp = transformActionStack;
       actionStackTemp != NULL;
       actionStackTemp = actionStackTemp->next) {

    action = (actionInformation *) actionStackTemp->entry;
    if (action->type != TRANSFORM_NOOP &&
	action->type != TRANSFORM_TRANSFORM) {

      /* Update maxFrames of each action with # sets actually processed
       * Only do this if unable to calculate the total number of frames,
       * by definition this means this is not an MPI run.
       */
      if (action->state->maxFrames<0)
        action->state->maxFrames = processed;
      for (i=0; i<6; i++)
	boxnew[i] = box[i];
      action->fxn(action, X, Y, Z, boxnew, PTRAJ_PRINT);
    }
  }

    /*
     *  ----------- PERFORM ANY REQUESTED ANALYSIS -------------
     */

  if (transformAnalyzeStack != NULL) {

    printfone("\nPTRAJ: Analyzing accumulated data\n");
    for (actionStackTemp = transformAnalyzeStack;
	 actionStackTemp != NULL;
	 actionStackTemp = actionStackTemp->next) {

      analyze = (analyzeInformation *) actionStackTemp->entry;
      if (analyze->type) {
	analyze->fxn(analyze, scalarStack, PTRAJ_ACTION);
	analyze->fxn(analyze, scalarStack, PTRAJ_PRINT);
      }
    }
  }

  /*
   *  Print benchmark stats to stdout and file bench.out
   *  If using MPI, send timings to rank 0, which
   *  will print out information for each process
   */

  if (bench) {
    printBench(stdout, 0);

    if (benchfile != NULL) {
      if (worldrank == 0)
	openFile(&bfile, benchfile, "w");
      printBench(bfile, benchshort);
      if (worldrank == 0)
	safe_fclose(bfile);
    }
  }


  ptrajCleanup();
  ptrajClearState(&startingState);
  safe_free(X);
  safe_free(Y);
  safe_free(Z);
  X = NULL;
  Y = NULL;
  Z = NULL;

}

void
printBench(FILE *fp, int bshort) {

  if (worldsize == 1) {

    if (fp == stdout || bshort == 0) {
      fprintf(fp, "Timings...\n");
      fprintf(fp, "\n-------------------------------\n");
      fprintf(fp, "| Check Input Time | %8.3f |\n", checkInputTime);
      fprintf(fp, "| Input Time       | %8.3f |\n", inputTime);
      fprintf(fp, "| Output Time      | %8.3f |\n", outputTime);
      fprintf(fp, "| Action Time      | %8.3f |\n", actionTime);
      fprintf(fp, "|------------------|----------|\n");
      fprintf(fp, "| Total Time       | %8.3f |\n", checkInputTime + inputTime + outputTime + actionTime);
      fprintf(fp, "-------------------------------\n");
      fprintf(fp, "\n");
    } else {
      fprintf(fp, "%8.3f\n", checkInputTime);
      fprintf(fp, "%8.3f\n", inputTime);
      fprintf(fp, "%8.3f\n", outputTime);
      fprintf(fp, "%8.3f\n", actionTime);
      fprintf(fp, "%8.3f\n", checkInputTime + inputTime + outputTime + actionTime);
    }

  } else {

#ifdef MPI

    int i;
    double *bufferCIT, *bufferIT, *bufferOT, *bufferAT;
    double *longest, *sum;
    
    /*
     *  Only rank 0 needs to initialize arrays to hold all timings
     */
    if (worldrank == 0) {
      bufferCIT = safe_malloc(worldsize * sizeof *bufferCIT);
      bufferIT  = safe_malloc(worldsize * sizeof *bufferIT);
      bufferOT  = safe_malloc(worldsize * sizeof *bufferOT);
      bufferAT  = safe_malloc(worldsize * sizeof *bufferAT);
      
      longest = safe_malloc(4 * sizeof *longest);
      sum = safe_malloc(4 * sizeof *sum);
      
      memset( longest, '\0', 4 );
      memset( sum, '\0', 4 );
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather( &checkInputTime, 1, MPI_DOUBLE, bufferCIT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather( &inputTime,      1, MPI_DOUBLE, bufferIT,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather( &outputTime,     1, MPI_DOUBLE, bufferOT,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather( &actionTime,     1, MPI_DOUBLE, bufferAT,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if (worldrank == 0) {
      for (i = 0; i < worldsize; i++) {
	if (bufferCIT[i] > longest[0])
	  longest[0] = bufferCIT[i];
	if (bufferIT[i] > longest[1])
	  longest[1] = bufferIT[i];
	if (bufferOT[i] > longest[2])
	  longest[2] = bufferOT[i];
	if (bufferAT[i] > longest[3])
	  longest[3] = bufferAT[i];
	
	sum[0] += bufferCIT[i];
	sum[1] += bufferIT[i];
	sum[2] += bufferOT[i];
	sum[3] += bufferAT[i];
      }
      
      if (fp == stdout || bshort == 0) {
	fprintf(fp, "Timings...\n");
	fprintf(fp, "\n");
	fprintf(fp, "--------------------");
	for (i = 0; i < worldsize; i++)
	  fprintf(fp, "-----------");
	fprintf(fp, "\n");
	fprintf(fp, "| Rank             |");
	for (i = 0; i < worldsize; i++)
	  fprintf(fp, "    %2d    |", i);
	fprintf(fp, "\n");
	fprintf(fp, "|------------------|");
	for (i = 0; i < worldsize; i++)
	  fprintf(fp, "----------|");
	fprintf(fp, "\n");
	fprintf(fp, "| Check Input Time |");
	for (i = 0; i < worldsize; i++)
	  fprintf(fp, " %8.3f |", bufferCIT[i]);
	fprintf(fp, "\n");
	fprintf(fp, "| Input Time       |");
	for (i = 0; i < worldsize; i++)
	  fprintf(fp, " %8.3f |", bufferIT[i]);
	fprintf(fp, "\n");
	fprintf(fp, "| Output Time      |");
	for (i = 0; i < worldsize; i++)
	  fprintf(fp, " %8.3f |", bufferOT[i]);
	fprintf(fp, "\n");
	fprintf(fp, "| Action Time      |");
	for (i = 0; i < worldsize; i++)
	  fprintf(fp, " %8.3f |", bufferAT[i]);
	fprintf(fp, "\n");
	fprintf(fp, "|------------------|");
	for (i = 0; i < worldsize; i++)
	  fprintf(fp, "----------|");
	fprintf(fp, "\n");
	fprintf(fp, "| Total Time       |");
	for (i = 0; i < worldsize; i++)
	  fprintf(fp, " %8.3f |", bufferCIT[i] + bufferIT[i] + bufferOT[i] + bufferAT[i]);
	fprintf(fp, "\n");
	fprintf(fp, "--------------------");
	for (i = 0; i < worldsize; i++)
	  fprintf(fp, "-----------");
	fprintf(fp, "\n");
      
	fprintf(fp, "\n");
	fprintf(fp, "-----------------------------------------------------");
	fprintf(fp, "\n");
	fprintf(fp, "|                  | Average  | Longest  |  Total   |");
	fprintf(fp, "\n");
	fprintf(fp, "|------------------|----------|----------|----------|");
	fprintf(fp, "\n");
	fprintf(fp, "| Check Input Time | %8.3f | %8.3f | %8.3f |", sum[0]/worldsize, longest[0], sum[0]);
	fprintf(fp, "\n");
	fprintf(fp, "| Input Time       | %8.3f | %8.3f | %8.3f |", sum[1]/worldsize, longest[1], sum[1]);
	fprintf(fp, "\n");
	fprintf(fp, "| Output Time      | %8.3f | %8.3f | %8.3f |", sum[2]/worldsize, longest[2], sum[2]);
	fprintf(fp, "\n");
	fprintf(fp, "| Action Time      | %8.3f | %8.3f | %8.3f |", sum[3]/worldsize, longest[3], sum[3]);
	fprintf(fp, "\n");
	fprintf(fp, "|------------------|----------|----------|----------|");
	fprintf(fp, "\n");
	fprintf(fp, "| Total Time       | %8.3f | %8.3f | %8.3f |", (sum[0] + sum[1] + sum[2] + sum[3]) / worldsize, longest[0] + longest[1] + longest[2] + longest[3], sum[0] + sum[1] + sum[2] + sum[3]);
	fprintf(fp, "\n");
	fprintf(fp, "-----------------------------------------------------");
	fprintf(fp, "\n");
      } else {
	fprintf(fp, "%8.3f", sum[0]/worldsize);
	fprintf(fp, "%8.3f", sum[1]/worldsize);
	fprintf(fp, "%8.3f", sum[2]/worldsize);
	fprintf(fp, "%8.3f", sum[3]/worldsize);
	fprintf(fp, "%8.3f", (sum[0] + sum[1] + sum[2] + sum[3]) / worldsize);
      }
    }
#endif    
  }
}
