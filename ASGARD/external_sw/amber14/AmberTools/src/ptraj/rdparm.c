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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/rdparm.c,v 10.5 2010/02/22 17:22:07 sbrozell Exp $
 *
 *  Revision: $Revision: 10.5 $
 *  Date: $Date: 2010/02/22 17:22:07 $
 *  Last checked in by $Author: sbrozell $
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
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#define RDPARM_MODULE
#include "ptraj.h"

/*
 *  This source file contains the bulk of the subroutines for the 
 *  rdparm functionality of the rdparm/ptraj programs.
 */


/*----------------------------------------------------------------

The following subroutines are defined...

   void        [global]
verbosity()
   ...toggles the verbosity status [verboseStatus] which is a 
   global, externally referenced variable.  This routine is *not* 
   a macro since it may be called externally and referenced by
   a function pointer.

   void        [global]
parmInfo()
   ...prints out information about the current parm file to
   stdout.

   int         [local]
unObfuscateAtom( int )
   ...inside the topology file, PARM for some unknown reason 
   obfuscates each atom number by (when positive), at = ( at + 3 ) / 3
   This routine, when passed an atom number as read from the parm 
   file converts it to the true number.

   char *      [local]
rdparmPrintAtom( int, char *, int )
   ...prints out the atom number, atom residue number, and 
   atom name.

   void        [local]
parmScanDouble( double *, int, char * )
   ...scan a double from the file specified by parm->fp

   void        [local]
parmScanString( char *, int, char * )
   ...scan a string from the file specified by parm->fp

   void        [local]
verbose( char *, char * )
   ...print out information, subject to the verbosity
   status [verboseStatus]

   void        [global]
readParm()
   ...reads in the parm file specified by parm->fp into the
   structure Parm *parm.  parm must be pre-allocated and
   the parm->fp must also be initialized.

   void        [global]
writeParm( FILE *, new )
   ...writes out the current Parm *parm structure (which may 
   have been modified) to the opened file specified by the 
   argument.  If new != 0 write a newer style prmtop.

   void        [global]
printBonds()
   ...pretty prints out all of the current bonds from the 
   Parm *parm structure to stdout.

   void        [global]
printPerturbedBonds()
   ...pretty prints out all of the current perturbed bonds
   in the current Parm *parm structure to stdout.

   void        [global]
printAngles()
   ...pretty prints out all of the current angles from the 
   Parm *parm structure to stdout.

   void        [global]
printPerturbedAngles()
   ...pretty prints out all of the current perturbed angles
   in the current Parm *parm structure to stdout.

   void        [global]
printDihedrals()
   ...pretty prints out all of the current perturbed
   dihedral angles from the Parm *parm structure to stdout.

   void        [global]
printPerturbedDihedrals()
   ...pretty prints out all of the current perturbed dihedrals
   in the current Parm *parm structure to stdout.

   void        [global]
printAtomInfo()
   ...pretty prints out atom information from the Parm *parm 
   structure to stdout.

   void        [global]
deleteBAD( char *, int )
   ...deletes the bond | angle | dihedral (which one depends on
   the first letter of the first argument) represented by the
   second argument.  NOTE: this is the number printed by
   pretty printing (i.e. printBonds) and the numbers will
   change after each deletion!

   void        [global]
deletePerturbedBAD( char *, int )
   ...deletes the perturbed bond | angle | dihedral (which one 
   depends on the first letter of the first argument) represented 
   by the second argument.  NOTE: this is the number printed by
   pretty printing (i.e. printPerturbedBonds) and the numbers will
   change after each deletion!

   int         [local]
processAtomInput(FILE *, FILE *, char *)
   ...prompts the user to input an atom specification. 
   NOTE: currently this only handles atom numbers, not
   midas style atom specification...

   double      [local]  <-- could be made global
processDoubleInput(FILE *, FILE *, char *)
   ...prompts the user (via the 1st argument, with output to the
   second) for a double value.  The third argument is a string
   which follows the "Please input " string.

   void        [global]
restrainBAD(char *)
   ...similar to deleteBAD, however this routine adds restraints
   for bonds | angles | dihedrals.  The user is prompted for the
   necessary information.

initializeParm( Parm *parm )
   ...initialize the recently allocated parm structure

   void        [global]
getParm( char * )
   ...prompts for filename and calls readParm to read in a AMBER
   parameter/topology file.

   void        [global]
putParm( char *)
   ...prompts for filename and calls writeParm to write out the
   current Parm *parm into an AMBER parameter/topology file.

   void        [global]
quit()
   ...quits the program by calling exit().

   void
doMardi2Sander [global]
   ...converts mardigras style constraints to sander style constraints.


   void
modifyBoxInfo()      [global]
   ...allows the user to change the box information for the current
   parm file

   void
modifyMolInfo  [global]
   ...analyzes the current molecule information searching for
   inconsistencies, such as the EDIT ADD water option screw up
   which places all the add waters into a single molecule.

   void
modifyTIP3P( int ) [global]
   ...strips or adds water to topology file specified by the global
   variable parm


-----------------------------------------------------------------*/


/*----------- GLOBAL VARIABLES, EXTERNALLY REFERENCED -----------*/
int verboseStatus = 0;
int isModifiedParm = FALSE;
int isWrittenParm = TRUE;
int isErrorParm = 0;
int newparm;
Parm *parm;


/*------------------------ LOCAL MACROS -------------------------*/

/* 
The following three macro's write out data specified by xxx starting
with i = srt and ending with i = stp-1.  xxx should represent an array
indexed by i, such as xxx = foo[i].bar.fum These are placed into
macros since they are used frequently to dump the FORTRAN style
formatted output.  They could potentially be replaced by a general
FORTRAN style output routines, however since the xxx can represent an
arbitrary structure indexed by i perhaps this means is better... 
*/

#define write20A4(fp, xxx,srt,stp) count=0; for (i=(srt); i<(stp); i++) { \
   fprintf(fp, "%4s", (xxx)); if (++count == 20) { count = 0; \
   fprintf(fp, "\n");  } } if ( (count) || !(stp)) fprintf(fp, "\n");

#define write5E16(fp, xxx,srt,stp) count=0; for (i=(srt); i<(stp); i++) { \
   fprintf(fp, "%16.8E", (xxx)); if (++count == 5) { count = 0; \
   fprintf(fp, "\n"); }} if ( (count) || !(stp) ) fprintf(fp, "\n");

#define write12I6(fp, xxx,srt,stp) count=0; for (i=(srt); i<(stp); i++) { \
   fprintf(fp, "%6d", (xxx)); if (++count == 12) { count = 0; \
   fprintf(fp, "\n"); }} if ( (count) || !(stp) ) fprintf(fp, "\n");

#define write10I8(fp, xxx,srt,stp) count=0; for (i=(srt); i<(stp); i++) { \
   fprintf(fp, "%8d", (xxx)); if (++count == 10) { count = 0; \
   fprintf(fp, "\n"); }} if ( (count) || !(stp) ) fprintf(fp, "\n");

/*------------------------ SUBROUTINES ---------------------------*/
/*             see comments above for each subroutine             */


void
prnlevSet(char *level)
{
  sscanf(level, "%i", &prnlev);
}
  

   void
NAME_COPY( FILE *fp, char *stringp, char *errorp)
{
  int i, j;
  int ierror;

  for (i = 0; i < NAME_SIZE - 1; i++) {

    ierror = fscanf(fp, "%c", &stringp[i]);
    if ( ierror == EOF ) 
      error("NAME_COPY", "... EOF encounters on scanning %s", errorp);

    while (i == 0 && stringp[i] == '\r' )
      ierror = fscanf(fp, "%c", &stringp[i]);

    while (i == 0 && stringp[i] == '\n' )
      ierror = fscanf(fp, "%c", &stringp[i]);

    if (i > 0 && i < NAME_SIZE-2 && stringp[i] == '\r') {
      fprintf(stdout, "Warning: Scanned \\r for character %i; scanning again\n", i+1);
      ierror = fscanf(fp, "%c", &stringp[i]);
    }
    if (i > 0 && i < NAME_SIZE-2 && stringp[i] == '\n') {
      fprintf(stdout, "Warning: Scanned \\n for character %i; replacing with space\n", i+1);
      for (j = i; j < NAME_SIZE-1; j++)
	stringp[j] = ' ';
      i = NAME_SIZE - 2;
    }
  }

  if ( stringp[NAME_SIZE-2] == '\r' ) {
    fprintf(stdout, "Warning: Scanned \\r for character %i; scanning again\n", NAME_SIZE-1);
    ierror = fscanf(fp, "%c", &stringp[i]);
  }
  if ( stringp[NAME_SIZE-2] == '\n' ) {
    fprintf(stdout, "Warning: Scanned \\n for character %i; replacing with space\n", NAME_SIZE-1);
    stringp[NAME_SIZE-2] = ' ';
  }

  stringp[NAME_SIZE-1] = (char) 0;


}


   void            /* global, referenced by dispatch() */
verbosity()
{
  if ( verboseStatus ) 
    verboseStatus = 0;
  else 
    verboseStatus = 1;
  fprintf(stdout, "Verbosity is now turned %s\n", 
	  (verboseStatus ? "ON" : "OFF"));
}

                   
   void            /* global, referenced by dispatch() */
parmInfo()   
{
  int i, count;

  fprintf(stdout, "\n  The parm file was loaded originally from a file\n");
  fprintf(stdout, "  named %s.\n", parm->filename);
  if ( isErrorParm ) {
    fprintf(stdout, "  There was an error upon reading in the parm file\n");
    fprintf(stdout, "  hence the following information may not be relevant\n");
  }
  if ( isModifiedParm ) {
    fprintf(stdout, "  The parm file has been modified.\n");
    if ( !isWrittenParm ) {
      fprintf(stdout, "  A modified parm file has not been written since\n");
      fprintf(stdout, "  last modification\n");
    } else
      fprintf(stdout, "  No modifications since last writeparm\n");
  }
  fprintf(stdout, "\n  The title is...\n\n%s\n\n", parm->title);
  fprintf(stdout, "  The control variables:\n\n");
  fprintf(stdout, "NTOTAT is %-9i (total number of atoms)\n", parm->NTOTAT);
  fprintf(stdout, "NTYPES is %-9i (number of atom types)\n", parm->NTYPES);
  fprintf(stdout, "NBONH  is %-9i (number of bonds w/ hydrogen)\n", 
	  parm->NBONH);
  fprintf(stdout, "NBONA  is %-9i (number of bonds w/out hydrogen)\n", 
	  parm->NBONA);
  fprintf(stdout, "NTHETH is %-9i (number of angles w/ hydrogen)\n", 
	  parm->NTHETH);
  fprintf(stdout, "NTHETA is %-9i (number of angles w/out hydrogen)\n", 
	  parm->NTHETA);
  fprintf(stdout, "NPHIH  is %-9i (number of dihedrals w/ hydrogen)\n", 
	  parm->NPHIH);
  fprintf(stdout, "NPHIA  is %-9i (number of dihedrals w/out hydrogen)\n", 
	  parm->NPHIA);
  if (parm->JPARM == 1) {
    fprintf(stdout, "NPARM  is %-9i (LES is ACTIVE)\n", 
	    parm->JPARM);
  }
  fprintf(stdout, "NEXT   is %-9i (number of excluded atoms)\n", 
	  parm->NEXT);
  fprintf(stdout, "NTOTRS is %-9i (number of residues)\n", parm->NTOTRS);
  fprintf(stdout, "MBONA  is %-9i (NBONA + number of constriant bonds)\n", 
	  parm->MBONA );
  fprintf(stdout, "MTHETS is %-9i (NTHETA + number of constraint angles)\n", 
	  parm->MTHETS);
  fprintf(stdout, "MPHIA  is %-9i (NPHIA + number of constraint dihedrals)\n",
	  parm->MPHIA);
  fprintf(stdout, "MUMBND is %-9i (number of unique bond types)\n", 
	  parm->MUMBND);
  fprintf(stdout, "MUMANG is %-9i (number of unique angle types)\n", 
	  parm->MUMANG);
  fprintf(stdout, "MPTRA  is %-9i (number of unique dihedral types)\n", 
	  parm->MPTRA);
  fprintf(stdout, "NATYP  is %-9i (# of ``atoms'' defined in parm topo)\n", 
	  parm->NATYP);
  fprintf(stdout, "NHB    is %-9i (# of type of h-bond pair interactions)\n", 
	  parm->NHB);

  fprintf(stdout, "IFBOX  is %-9i (== 1 IF periodic box info)\n", 
	 parm->IFBOX);
  if ( parm->IFBOX == 1 ) {
    fprintf(stdout, " BOX SIZE is %g by %g by %g angstroms\n", 
	    parm->box->box[0], parm->box->box[1], parm->box->box[2]);
    fprintf(stdout, " IPTRES is %-8i (final residue that is part of solute)\n",
	   parm->box->iptres);
    fprintf(stdout, " NSPM   is %-8i (the total number of molecules)\n",
	   parm->box->nspm);
    fprintf(stdout, " NSPSOL is %-8i (first solvent molecule)\n",
	   parm->box->nspsol);
    fprintf(stdout, " NUMBER of atoms in each molecule, molecule %i to %i\n",
	   1, parm->box->nspm);
    count=1;
    for (i= 0; i < parm->box->nspm; i++) {
      fprintf(stdout, "%4i ", parm->box->nsp[i]);
      if ( count != 0 && count % 20 == 0 ) fprintf(stdout, "\n"); 
      count++;
    }
    if ( count % 10 ) fprintf(stdout, "\n");
    fprintf(stdout, " BETA   is %-8.2f (angle of box)\n", parm->box->beta);
  }
  fprintf(stdout, "NMXRS  is %-9i (number of atoms in largest residue)\n", 
	 parm->NMXRS);
  fprintf(stdout, "IFCAP  is %-9i (== 1 IF cap option was used in edit)\n", 
	 parm->IFCAP);

  /* write perturbation info if IFPERT > 0 */
  fprintf(stdout, 
	  "IFPERT is %-9i (== 1 IF perturbation info, else 0)\n", 
	  parm->IFPERT);
  if ( parm->IFPERT == 1 ) {
    fprintf(stdout, 
	    " NBPER  is %-8i (number of bonds to be perturbed)\n", 
	    parm->NBPER);
    fprintf(stdout, 
	    " NGPER  is %-8i (number of angles to be perturbed)\n", 
	    parm->NGPER);
    fprintf(stdout, 
	    " NDPER  is %-8i (number of dihedrals to be perturbed)\n", 
	    parm->NDPER);
    fprintf(stdout, 
	    " MBPER  is %-8i (# of pert bonds across bounds to non-pert)\n", 
	    parm->MBPER);
    fprintf(stdout, 
	    " MGPER  is %-8i (# of pert angles across bounds to non-pert)\n", 
	    parm->MGPER);
    fprintf(stdout, 
	    " MDPER  is %-8i (# of pert diheds across bounds to non-pert)\n", 
	    parm->MDPER);
  }
  fprintf(stdout, "\n\n");
}


   int             /* local, referenced by rdparm() */
unObfuscateAtom( int at ) 
{
  if ( at < 0 )
    at = - ( -at + 3 ) / 3;
  else
    at = ( at + 3 ) / 3;
  return( at );
}


   char *          /* local, referenced by multiple routines */
rdparmPrintAtom( int at, char *str ) 
{
  int rs;
  
  rs = parm->atom[at-1].res+1;
  sprintf(str, ":%i@%s", rs, parm->atom[at-1].igraph);
  return(str);
}


   void            /* local, referenced by rdparm() */
parmScanDouble(double *value, int isError, char *err_string) 
{
  if ( scanDouble(parm->fp, value, isError, err_string ) )
    if ( isErrorParm == 0 ) 
      isErrorParm = 1;
}


   void            /* local, referenced by rdparm() */
parmScanString(char *value, int isError, char *err_string) 
{
  if ( scanString(parm->fp, value, isError, err_string ) )
    if ( isErrorParm == 0 ) 
      isErrorParm = 1;
}


char *buffer1_12I6 = NULL;
char *buffer2_12I6 = NULL;
int   current_12I6 = -1;

   int            /* local, referenced by rdparm() */
loadAndReturn12I6(FILE *fp, int status, char *errorString)
{
  int i, j, currentValue = 0;

  if( newparm ){

    if( status ) fscanf(parm->fp, "%d", &currentValue);

  } else {

    /*
     * initialize buffers if necessary.  buffer1_12I6 will hold the initial 
     * line of text and buffer2_12I6 will hold the ' ' padded text to properly
     * define 6 character integer boundaries
     */
    if (buffer1_12I6 == NULL) {
      buffer1_12I6 = (char *) safe_malloc(sizeof(char) * BUFFER_SIZE);
      buffer1_12I6[0] = (char) 0;
    }
    if (buffer2_12I6 == NULL) {
      buffer2_12I6 = (char *) safe_malloc(sizeof(char) * BUFFER_SIZE);
      buffer2_12I6[0] = (char) 0;
    }
  
    /*
     * this is to check for the initialization status and this routine should
     * be called with this parameter after every series of reads
     */
    if (status == 0) {
      current_12I6 = -1;
      return -1;
    }
    /*
     * if the status < 0, then a line of text needs to be read
     */
    if (status < 0) {
      current_12I6 = -1;
      if ( fgets(buffer1_12I6, 120, fp ) == NULL )
        fprintf(stderr,"Error scanning a line of 12I6: %s\n", errorString);
      return -1;
    }
  
    /*
     * if this is the first number to be read from the line, grab a line of 
     * text and pad it appropriately...
     */
    while ( current_12I6    == -1 || current_12I6    == 12 ) {
      if ( fgets(buffer1_12I6, 120, fp ) == NULL )
        fprintf(stderr,"Error scanning a line of 12I6: %s\n", errorString);
      /*
       * remove spurious newlines
       */
      if (buffer1_12I6[0] == '\n' || buffer1_12I6[0] == '\r')
        if ( fgets(buffer1_12I6, 120, fp ) == NULL )
  	fprintf(stderr,"Error scanning a line of 12I6: %s\n", errorString);
  
      j = 0;
      for (i=0; i < 12; i++) {
        strncpy(buffer2_12I6+j, buffer1_12I6+(i*6), 6);
        j = j + 6;
        buffer2_12I6[j++] = ' ';
      }
      buffer2_12I6[j] = (char) 0;
      current_12I6 = 0;
    }
    
  
    if (sscanf(buffer2_12I6+(current_12I6*7), "%6d", &currentValue) < 1)
      fprintf(stderr, "Error scanning a value: %s\n", errorString);
    
    current_12I6++;
  }
  return(currentValue);
}
  
  



   void            /* local */
verbose(char *info_string, char *err_string) 
{
  if (verboseStatus)
    fprintf(stdout, "%s\n", info_string);
  if (isErrorParm == 1) {
    fprintf(stdout,"WARNING: an error upon reading the parm file has\n");
    fprintf(stdout,"been reported.  Proceed with caution since some of\n");
    fprintf(stdout,"the data for the current parm file may be undefined\n");
    if ( strcmp( err_string, "" ) != 0 ) 
      fprintf(stdout, "ERROR: %s\n", err_string);
    isErrorParm = 2;
  }
}


   void
clearParm(Parm *parm)
{ /* routine to clear out an old parm structure... */
  isErrorParm = 0;
  isModifiedParm = 0;
  if ( parm == NULL )
    return;

  safe_free((void *) parm->filename);
  safe_free((void *) parm->atom);
  safe_free((void *) parm->nno);
  safe_free((void *) parm->residue);          
  safe_free((void *) parm->rk);
  safe_free((void *) parm->req);
  safe_free((void *) parm->tk);
  safe_free((void *) parm->teq);
  safe_free((void *) parm->pk);
  safe_free((void *) parm->pn);
  safe_free((void *) parm->phase);   
  safe_free((void *) parm->solty);             
  safe_free((void *) parm->cn1);
  safe_free((void *) parm->cn2);                 
  safe_free((void *) parm->pbondH);
  safe_free((void *) parm->pbond);
  safe_free((void *) parm->pangleH);
  safe_free((void *) parm->pangle);         
  safe_free((void *) parm->pdihedralH);
  safe_free((void *) parm->pdihedral);
  safe_free((void *) parm->natex);                
  safe_free((void *) parm->ag);                
  safe_free((void *) parm->bg);                
  safe_free((void *) parm->hbcut);
  if (parm->IFBOX) {
    safe_free((void *) parm->box->nsp);
    safe_free((void *) parm->box);
  }
  if (parm->IFCAP) 
    safe_free((void *) parm->cap);

  if (parm->IFPERT) {
    safe_free((void *) parm->pert->ibper);
    safe_free((void *) parm->pert->jbper);        
    safe_free((void *) parm->pert->icbper);               
    safe_free((void *) parm->pert->itper);
    safe_free((void *) parm->pert->jtper);
    safe_free((void *) parm->pert->ktper);
    safe_free((void *) parm->pert->ictper);                       
    safe_free((void *) parm->pert->ipper);
    safe_free((void *) parm->pert->jpper);
    safe_free((void *) parm->pert->kpper);
    safe_free((void *) parm->pert->lpper);
    safe_free((void *) parm->pert->icpper);                       
    safe_free((void *) parm->pert->labper);                      
    safe_free((void *) parm->pert->igrper);                      
    safe_free((void *) parm->pert->ismper);                      
    safe_free((void *) parm->pert->almper);                    
    safe_free((void *) parm->pert->iaper);                        
    safe_free((void *) parm->pert->iacper);                       
    safe_free((void *) parm->pert->cgper);                         
    safe_free((void *) parm->pert);
  }

  safe_free((void *) parm->bond);
  safe_free((void *) parm->angle);
  safe_free((void *) parm->dihedral);
  if ( parm->coords != NULL ) {
    safe_free((void *) parm->coords->title);
    safe_free((void *) parm->coords->x);
    safe_free((void *) parm->coords->y);
    safe_free((void *) parm->coords->z);
    safe_free((void *) parm->coords->vx);
    safe_free((void *) parm->coords->vy);
    safe_free((void *) parm->coords->vz);
    safe_free((void *) parm->coords);
  }
  safe_free((void *) parm);
  parm = NULL;
}


   static int    /* local; used by readParm() and find_flag()   */
s_getline( char *line, int max, FILE *fp )
{
  /*
   *   "safe" version of getline: exits with error in EOF is found
   */

  if( fgets( line, max, fp ) == NULL )
    return( 0 );
  else
    return strlen( line );

}

  static int    /* local, used by readParm()  */
find_flag( FILE *fp, char *label )
{
   /*
    *  Do not assume flags are in the proper order...
    */

  char buffer[ BUFFER_SIZE ];
  char *bufferp;
  int status;

  status = 2;
  while (status > 0) {

    if ( s_getline( buffer, BUFFER_SIZE, parm->fp ) != 0 ) {

      if (strncmp(buffer, "%FLAG", 5) == 0) {

	bufferp = buffer+5;
	skipWhitespace(bufferp);
	if ( strncmp( bufferp, label, strlen(label) ) == 0 ) {
	  /*
	   *  found match so read another line of text and exit
	   */
	  s_getline( buffer, BUFFER_SIZE, parm->fp );
	  return 1;
	}
      }
    } else {
      /*
       *  reopen the file and search from the beginning if status > 0 else exit
       */

      if (status > 1) {
        /*parm->fp = safe_freopen(parm->fp);  */
        /* DAN ROE: Just rewind the file pointer instead  */
        rewind(parm->fp);
	if (prnlev > 3)
	  printf("RE-OPENED PRMTOP FILE\n");
      }
      status--;
    }
  }

  if (worldrank == 0)
    fprintf(stdout, " PRMTOP does not contain %%FLAG %s\n", label);
  return 0;
}

   void            /* global, referenced by dispatch() */
readParm() 
{
  int i, j, len;
  char *buffer;
  char outputbuffer[121];
  char previousbuffer[121];

  isErrorParm = 0;

  /*
   *  This routine will read in an AMBER parm file of format
   *  at least from AMBER 3A to 4.1.  For more information on 
   *  the parameter file format, see
   *
   *  http://www.amber.ucsf.edu/amber/formats.html
   *
   */


  /* 
   *  Read in the title...
   *
   *  NOTE: although the parameter file supposedly has titles
   *  of up to 20A4 or 80 characters only, in practice this is not
   *  always the case (such as with programs such as INTERFACE).
   *  therefore we read in the first complete line and grab out the
   *  title from this...
   *
   */

  buffer = (char *) safe_malloc(sizeof(char) * BUFFER_SIZE);
  s_getline( buffer, BUFFER_SIZE, parm->fp );

  /*  see if this is a "new" format parm file:    */

  if( !strncmp( buffer, "%VERSION", 8 ) ){
     newparm = 1;

     /* TITLE = regular prmtop, CTITLE = Chamber prmtop. */
     if ( (find_flag( parm->fp, "TITLE" ) == 0) && (find_flag( parm->fp, "CTITLE" ) == 0) ) {
       error("readParm", "...failed to find TITLE\n");
     }
     len = s_getline( buffer, BUFFER_SIZE, parm->fp );
     strncpy(parm->title, buffer, len-1);
  } else {
     newparm = 0;
     strncpy(parm->title, buffer, TITLE_LENGTH);
  }
  verbose("Scanned in title...", "");

  /*
   *  Read in the integer control variables (3 lines) 
   */
  if ( newparm && find_flag( parm->fp, "POINTERS" ) == 0 ) {
      error("readParm", "...failed to find POINTERS\n");
  }

  parm->NTOTAT = loadAndReturn12I6(parm->fp, 1, "parm->NTOTAT");
  parm->NTYPES = loadAndReturn12I6(parm->fp, 1, "parm->NTYPES");
  parm->NBONH  = loadAndReturn12I6(parm->fp, 1, "parm->NBONH");
  parm->NBONA  = loadAndReturn12I6(parm->fp, 1, "parm->NBONA");
  parm->NTHETH = loadAndReturn12I6(parm->fp, 1, "parm->NTHETH");
  parm->NTHETA = loadAndReturn12I6(parm->fp, 1, "parm->NTHETA");
  parm->NPHIH  = loadAndReturn12I6(parm->fp, 1, "parm->NPHIH");
  parm->NPHIA  = loadAndReturn12I6(parm->fp, 1, "parm->NPHIA");
  parm->JHPARM = loadAndReturn12I6(parm->fp, 1, "parm->JHPARM");
  parm->JPARM  = loadAndReturn12I6(parm->fp, 1, "parm->JPARM");
  parm->NEXT   = loadAndReturn12I6(parm->fp, 1, "parm->NEXT");
  parm->NTOTRS = loadAndReturn12I6(parm->fp, 1, "parm->NTOTRS");
  loadAndReturn12I6(parm->fp, 0, "");

  parm->MBONA  = loadAndReturn12I6(parm->fp, 1, "parm->MBONA");
  parm->MTHETS = loadAndReturn12I6(parm->fp, 1, "parm->MTHETS");
  parm->MPHIA  = loadAndReturn12I6(parm->fp, 1, "parm->MPHIA");
  parm->MUMBND = loadAndReturn12I6(parm->fp, 1, "parm->MUMBND");
  parm->MUMANG = loadAndReturn12I6(parm->fp, 1, "parm->MUMANG");
  parm->MPTRA  = loadAndReturn12I6(parm->fp, 1, "parm->MPTRA");
  parm->NATYP  = loadAndReturn12I6(parm->fp, 1, "parm->NATYP");
  parm->NHB    = loadAndReturn12I6(parm->fp, 1, "parm->NHB");
  parm->IFPERT = loadAndReturn12I6(parm->fp, 1, "parm->IFPERT");
  parm->NBPER  = loadAndReturn12I6(parm->fp, 1, "parm->NBPER");
  parm->NGPER  = loadAndReturn12I6(parm->fp, 1, "parm->NGPER");
  parm->NDPER  = loadAndReturn12I6(parm->fp, 1, "parm->NDPER");
  loadAndReturn12I6(parm->fp, 0, "");

  parm->MBPER  = loadAndReturn12I6(parm->fp, 1, "parm->MBPER");
  parm->MGPER  = loadAndReturn12I6(parm->fp, 1, "parm->MGPER");
  parm->MDPER  = loadAndReturn12I6(parm->fp, 1, "parm->MDPER");
  parm->IFBOX  = loadAndReturn12I6(parm->fp, 1, "parm->IFBOX");
  parm->NMXRS  = loadAndReturn12I6(parm->fp, 1, "parm->NMXRS");
  parm->IFCAP  = loadAndReturn12I6(parm->fp, 1, "parm->IFCAP");
  if( newparm )
    parm->NUMEXTRA  = loadAndReturn12I6(parm->fp, 1, "parm->NUMEXTRA");
  else
    parm->NUMEXTRA = 0;
  loadAndReturn12I6(parm->fp, 0, "");

  if (prnlev > 2)
    fprintf(stdout, "Read in control variables\n");
  if (prnlev > 4)
    fprintf(stdout, "%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i\n", 
	    parm->NTOTAT, parm->NTYPES, parm->NBONH, 
	    parm->NBONA,  parm->NTHETH, parm->NTHETA, 
	    parm->NPHIH,  parm->NPHIA,  parm->JHPARM, 
	    parm->JPARM,  parm->NEXT,   parm->NTOTRS);
  if (prnlev > 4)
    fprintf(stdout, "%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i\n", 
	  parm->MBONA,  parm->MTHETS, parm->MPHIA,  
	  parm->MUMBND, parm->MUMANG, parm->MPTRA,
	  parm->NATYP,  parm->NHB,    parm->IFPERT,
	  parm->NBPER,  parm->NGPER,  parm->NDPER);
                   
  if (prnlev > 4)
    fprintf(stdout,"%6i%6i%6i%6i%6i%6i\n",
	  parm->MBPER,  parm->MGPER,  parm->MDPER,  
	  parm->IFBOX,  parm->NMXRS,  parm->IFCAP);


  /*
   *  Allocate Atom space and read in atom names. 
   */

  if ( newparm == 0 || find_flag( parm->fp, "ATOM_NAME" ) ) {
    parm->atom = safe_malloc( sizeof( Atom ) * parm->NTOTAT );
    for (i=0; i < parm->NTOTAT; i++) {
      NAME_COPY(parm->fp, parm->atom[i].igraph, "igraph");
    }

    if (prnlev > 2) 
      fprintf(stdout, "Read in atom names...\n");
    if (prnlev > 4) {
      for (i=0; i < parm->NTOTAT; i++) {
	fprintf(stdout, "%4s", parm->atom[i].igraph);
	if ( (i+1)%20 == 0 ) {
	  fprintf(stdout, "\n");
	}
      }
      fprintf(stdout, "\n");
    }
  }

  /*
   *  Read in charges 
   */
  if ( newparm == 0 || find_flag( parm->fp, "CHARGE" ) ) {
    for (i=0; i < parm->NTOTAT; i++) 
      parmScanDouble(&parm->atom[i].chrg, FALSE, "chrg");

    if (prnlev > 2)
      fprintf(stdout, "Read in charges...\n");
    if (prnlev > 4) {
      for (i=0; i < parm->NTOTAT; i++) {
	fprintf(stdout, "%16.8E", parm->atom[i].chrg);
	if ( (i+1)%5 == 0 ) fprintf(stdout, "\n");
      }
      fprintf(stdout, "\n");
    }
  }

  /*
   *  Read in masses 
   */
  if ( newparm == 0 || find_flag( parm->fp, "MASS" ) ) {
    for (i=0; i < parm->NTOTAT; i++) 
      parmScanDouble(&parm->atom[i].amass, FALSE, "amass");

    if (prnlev > 2)
      fprintf(stdout, "Read in masses...\n");
    if (prnlev > 4) {
      for (i=0; i < parm->NTOTAT; i++) {
	fprintf(stdout, "%16.8E", parm->atom[i].amass);
	if ( (i+1)%5 == 0 ) fprintf(stdout, "\n");
      }
      fprintf(stdout, "\n");
    }
  }


  /*
   *  Read in iac array
   */
  if ( newparm == 0 || find_flag( parm->fp, "ATOM_TYPE_INDEX" ) ) {
    for (i=0; i < parm->NTOTAT; i++) 
      parm->atom[i].iac = loadAndReturn12I6(parm->fp, 1, "iac");
    loadAndReturn12I6(parm->fp, 0, "");

    if (prnlev > 2)
      fprintf(stdout, "Read in IAC (atoms involved in L-J)...\n");
    if (prnlev > 4) {
      for (i=0; i < parm->NTOTAT; i++) {
	fprintf(stdout, "%6i", parm->atom[i].iac);
	if ( (i+1)%12 == 0 ) fprintf(stdout, "\n");
      }
      fprintf(stdout, "\n");
    }
  }

  /*
   *  Read in numex 
   */
  if ( newparm == 0 || find_flag( parm->fp, "NUMBER_EXCLUDED_ATOMS" ) ) {
    for (i=0; i < parm->NTOTAT; i++) 
      parm->atom[i].numex = loadAndReturn12I6(parm->fp, 1, "numex");
    loadAndReturn12I6(parm->fp, 0, "");

    if (prnlev > 2)
      fprintf(stdout, "Read in NUMEX (index to excl atom list)...\n");
    if (prnlev > 4) {
      for (i=0; i < parm->NTOTAT; i++) {
	fprintf(stdout, "%6i", parm->atom[i].numex);
	if ( (i+1)%12 == 0 ) fprintf(stdout, "\n");
      }
      fprintf(stdout, "\n");
    }
  }


  /*
   *  Read in nno (ico) after allocating memory 
   */
  if ( newparm == 0 || find_flag( parm->fp, "NONBONDED_PARM_INDEX" ) ) {
    len = parm->NTYPES; 
    len = len * len;
    parm->nno = safe_malloc( sizeof( int ) * len );
    for (i=0; i < len; i++) 
      parm->nno[i] = loadAndReturn12I6(parm->fp, 1, "nno");
    loadAndReturn12I6(parm->fp, 0, "");

    if (prnlev > 2)
      fprintf(stdout, "Read in NNO (index for nonbond of @type)...\n");
    if (prnlev > 4) {
      for (i=0; i < len; i++) {
	fprintf(stdout, "%6i", parm->nno[i]);
	if ( (i+1)%12 == 0 ) fprintf(stdout, "\n");
      }
      fprintf(stdout, "\n");
    }
  }


  /*
   *  Read in residue labels after allocating memory 
   */
  if ( newparm == 0 || find_flag( parm->fp, "RESIDUE_LABEL" ) ) {
    parm->residue = safe_malloc( sizeof( Residue ) * (parm->NTOTRS + 1) );
    for (i=0; i < parm->NTOTRS; i++) {
      NAME_COPY(parm->fp, parm->residue[i].labres, "labres");
    }

    if (prnlev < 2) 
      if (worldrank == 0)
	fprintf(stdout, "Residue labels:\n\n");
    if (prnlev > 3) {

      j = 0;
      len = 0;
      previousbuffer[0] = (char) 0;
      for (i=0; i < parm->NTOTRS; i++) {

	sprintf(outputbuffer+j, "%5s", parm->residue[i].labres);
	j += 5;
	if ( i == parm->NTOTRS - 1 || ( i != 0 && (i+1) % 10 == 0 ) ) {
	  outputbuffer[j] = (char) 0;
	  if ( strcmp(previousbuffer, outputbuffer) == 0 ) {
	    if ( ! len ) {
	      fprintf(stdout, " ...\n");
	      len = 1;
	    }
	  } else {
	    fprintf(stdout, "%s\n", outputbuffer);
	    len = 0;
	    strcpy(previousbuffer, outputbuffer);
	  }
	  j = 0;
	}
      }
    }
  
    if (prnlev < 1) {

      j = 0;
      len = 0;
      previousbuffer[0] = (char) 0;
      for (i=0; i < parm->NTOTRS; i++) {
      
	sprintf(outputbuffer+j, "%5s", parm->residue[i].labres);
	j += 5;
	if ( i == parm->NTOTRS - 1 || ( i != 0 && (i+1) % 10 == 0 ) ) {
	  outputbuffer[j] = (char) 0;
	  if ( strcmp(previousbuffer, outputbuffer) == 0 ) {
	    if ( ! len ) {
	      fprintf(stdout, " ...\n");
	      len = 1;
	    }
	  } else {
	    if (worldrank == 0)
	      fprintf(stdout, "%s\n", outputbuffer);
	    len = 0;
	    strcpy(previousbuffer, outputbuffer);
	  }
	  j = 0;
	}
      }
      fprintf(stdout, "\n");
    }
  }

  /*
   *  Read in ipres 
   */
  if ( newparm == 0 || find_flag( parm->fp, "RESIDUE_POINTER" ) ) {
    for (i=0; i < parm->NTOTRS; i++) 
      parm->residue[i].ipres = loadAndReturn12I6(parm->fp, 1, "ipres");
    loadAndReturn12I6(parm->fp, 0, "");

    parm->residue[parm->NTOTRS].ipres = parm->NTOTAT + 1;
    /*
     *  fill in residues for each atom 
     */
    for (i=0; i < parm->NTOTRS; i++) 
      for (j=parm->residue[i].ipres-1; j < parm->residue[i+1].ipres - 1; j++)
	parm->atom[j].res = i;
    parm->atom[parm->NTOTAT-1].res = parm->NTOTRS-1;

    if (prnlev > 2)
      fprintf(stdout, "Read in the residue to atom pointer list...\n");
    if (prnlev > 4) {
      for (i=0; i < parm->NTOTRS; i++) {
	fprintf(stdout, "%6i", parm->residue[i].ipres);
	if ( (i+1)%12 == 0 ) fprintf(stdout, "\n");
      }
      fprintf(stdout, "\n");
    }
  }

  /*
   *  Read in bond parameters after allocating memory 
   */
  parm->rk = safe_malloc( sizeof( double ) * parm->MUMBND );
  parm->req = safe_malloc( sizeof( double ) * parm->MUMBND );
  if ( newparm == 0 || (newparm && find_flag( parm->fp, "BOND_FORCE_CONSTANT" )) ) {
    for (i=0; i < parm->MUMBND; i++)
      parmScanDouble(&parm->rk[i], FALSE, "rk");
  }

  if ( newparm == 0 || find_flag( parm->fp, "BOND_EQUIL_VALUE" ) ) {
    for (i=0; i < parm->MUMBND; i++)
      parmScanDouble(&parm->req[i], FALSE, "req");

    if (prnlev > 2)
      fprintf(stdout, "Read in bond parameters RK and REQ...\n");
    if (prnlev > 4) {
      for (i=0; i < parm->MUMBND; i++) {
	fprintf(stdout, "%16.8E", parm->rk[i]);
	if ( (i+1)%5 == 0 ) fprintf(stdout, "\n");
      }
      if ( i%5 ) fprintf(stdout, "\n");
      
      for (i=0; i < parm->MUMBND; i++) {
	fprintf(stdout, "%16.8E", parm->req[i]);
	if ( (i+1)%5 == 0 ) fprintf(stdout, "\n");
      }
      if ( i%5 ) fprintf(stdout, "\n");
    }
  }


  /*
   *  Read in angle parameters after allocating memory 
   */
  parm->tk = safe_malloc( sizeof( double ) * parm->MUMANG );
  parm->teq = safe_malloc( sizeof( double ) * parm->MUMANG );
  if ( newparm == 0 || find_flag( parm->fp, "ANGLE_FORCE_CONSTANT" ) ) {
    for (i=0; i < parm->MUMANG; i++)
      parmScanDouble(&parm->tk[i], FALSE, "tk");
  }

  if ( newparm == 0 || find_flag( parm->fp, "ANGLE_EQUIL_VALUE" ) ) {
    for (i=0; i < parm->MUMANG; i++)
      parmScanDouble(&parm->teq[i], FALSE, "teq");

    if (prnlev > 2)
      fprintf(stdout, "Read in angle parameters TK and TEQ...\n");
    if (prnlev > 4) {
      for (i=0; i < parm->MUMANG; i++) {
	fprintf(stdout, "%16.8E", parm->tk[i]);
	if ( (i+1)%5 == 0 ) fprintf(stdout, "\n");
      }
      if ( i%5 ) fprintf(stdout, "\n");
      for (i=0; i < parm->MUMANG; i++) {
	fprintf(stdout, "%16.8E", parm->teq[i]);
	if ( (i+1)%5 == 0 ) fprintf(stdout, "\n");
      }
      if ( i%5 ) fprintf(stdout, "\n");
    }
  }

  /*
   *  Read in dihedral parameters after allocating memory 
   */
  parm->pk = safe_malloc( sizeof( double ) * parm->MPTRA );
  parm->pn = safe_malloc( sizeof( double ) * parm->MPTRA );
  parm->phase = safe_malloc( sizeof( double ) * parm->MPTRA );

  if ( newparm == 0 || find_flag( parm->fp, "DIHEDRAL_FORCE_CONSTANT" ) ) {
    for (i=0; i < parm->MPTRA; i++) 
      parmScanDouble(&parm->pk[i], FALSE, "pk");
  }

  if ( newparm == 0 || find_flag( parm->fp, "DIHEDRAL_PERIODICITY" ) ) {
    for (i=0; i < parm->MPTRA; i++) 
      parmScanDouble(&parm->pn[i], FALSE, "pn");
  }

  if ( newparm == 0 || find_flag( parm->fp, "DIHEDRAL_PHASE" ) ) {
    for (i=0; i < parm->MPTRA; i++) 
      parmScanDouble(&parm->phase[i], FALSE, "phase");

    if (prnlev > 2)
      fprintf(stdout, "Read in dihedral parameters PK, PN and PHASE...\n");
    if (prnlev > 4) {
      for (i=0; i < parm->MPTRA; i++) {
	fprintf(stdout, "%16.8E", parm->pk[i]);
	if ( (i+1)%5 == 0 ) fprintf(stdout, "\n");
      }
      if ( i%5 ) fprintf(stdout, "\n");
      for (i=0; i < parm->MPTRA; i++) {
	fprintf(stdout, "%16.8E", parm->pn[i]);
	if ( (i+1)%5 == 0 ) fprintf(stdout, "\n");
      }
      if ( i%5 ) fprintf(stdout, "\n");
      for (i=0; i < parm->MPTRA; i++) {
	fprintf(stdout, "%16.8E", parm->phase[i]);
	if ( (i+1)%5 == 0 ) fprintf(stdout, "\n");
      }
      if ( i%5 ) fprintf(stdout, "\n");
    }
  }

  /*
   *  Read in solty 
   */
  if ( newparm == 0 || find_flag( parm->fp, "SOLTY" ) ) {
    parm->solty = safe_malloc( sizeof ( double ) * parm->NATYP );
    for (i=0; i < parm->NATYP; i++)
      parmScanDouble(&parm->solty[i], FALSE, "solty");

    if (prnlev > 2)
      fprintf(stdout, "Read in SOLTY...\n");
    if (prnlev > 4) {
      for (i=0; i < parm->NATYP; i++) {
	fprintf(stdout, "%16.8E", parm->solty[i]);
	if ( (i+1)%5 == 0 ) fprintf(stdout, "\n");
      }
      fprintf(stdout, "\n");
    }
  }


  /*
   *  Read in L-J parameters CN1 and CN2 
   */
  len = parm->NTYPES;
  len = len * (len + 1) / 2;

  parm->cn1 = safe_malloc( sizeof( double) * len );
  if ( newparm == 0 || find_flag( parm->fp, "LENNARD_JONES_ACOEF" ) ) {
    for (i=0; i < len; i++)
      parmScanDouble(&parm->cn1[i], FALSE, "cn1");
  }


  parm->cn2 = safe_malloc( sizeof( double) * len );
  if ( newparm == 0 || find_flag( parm->fp, "LENNARD_JONES_BCOEF" ) ) {
    for (i=0; i < len; i++)
      parmScanDouble(&parm->cn2[i], FALSE, "cn2");

    if (prnlev > 2)
      fprintf(stdout, "Read in L-J parameters CN1 and CN2...\n");
    if (prnlev > 4) {
      for (i=0; i < len; i++) {
	fprintf(stdout, "%16.8E", parm->cn1[i]);
	if ( (i+1)%5 == 0 ) fprintf(stdout, "\n");
      }
      if (i%5) fprintf(stdout,"\n");
      for (i=0; i < len; i++) {
	fprintf(stdout, "%16.8E", parm->cn2[i]);
	if ( (i+1)%5 == 0 ) fprintf(stdout, "\n");
      }
      fprintf(stdout, "\n");
    }
  }


  /*
   *  Read in info for bonds with hydrogen, post allocating memory 
   */
  parm->pbondH = safe_malloc( sizeof( ParmBond ) * parm->NBONH );

  if ( newparm == 0 || find_flag( parm->fp, "BONDS_INC_HYDROGEN" ) ) {
    for (i=0; i < parm->NBONH; i++) {

      parm->pbondH[i].ib  = loadAndReturn12I6(parm->fp, 1, "bonds with hydrogen, IB");
      parm->pbondH[i].jb  = loadAndReturn12I6(parm->fp, 1, "bonds with hydrogen, JB");
      parm->pbondH[i].icb = loadAndReturn12I6(parm->fp, 1, "bonds with hydrogen, ICB");

    }
    loadAndReturn12I6(parm->fp, (parm->NBONH ? 0 : -1), "");

    if (prnlev > 2) 
      fprintf(stdout, "Read in info for bonds w/ hydrogen...\n");
    if (prnlev > 4) {
      j = 0;
      for (i=0; i < parm->NBONH; i++) {
	fprintf(stdout, "%6i", parm->pbondH[i].ib);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
	fprintf(stdout, "%6i", parm->pbondH[i].jb);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
	fprintf(stdout, "%6i", parm->pbondH[i].icb);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
      }
      fprintf(stdout, "\n");
    }
  }

  /*
   *  Read in info for bonds withOUT hydrogen, post allocating memory 
   */
  parm->pbond = safe_malloc( sizeof( ParmBond ) * parm->MBONA );

  if ( newparm == 0 || find_flag( parm->fp, "BONDS_WITHOUT_HYDROGEN" ) ) {
    for (i=0; i < parm->MBONA; i++) {

      parm->pbond[i].ib  = loadAndReturn12I6(parm->fp, 1, "bonds without hydrogen, IB");
      parm->pbond[i].jb  = loadAndReturn12I6(parm->fp, 1, "bonds without hydrogen, JB");
      parm->pbond[i].icb = loadAndReturn12I6(parm->fp, 1, "bonds without hydrogen, ICB");

    }
    loadAndReturn12I6(parm->fp, (parm->MBONA ? 0 : -1), "");

    if (prnlev > 2) 
      fprintf(stdout, "Read in info for bonds w/out hydrogen...\n");
    if (prnlev > 4) {
      j = 0;
      for (i=0; i < parm->MBONA; i++) {
	fprintf(stdout, "%6i", parm->pbond[i].ib);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
	fprintf(stdout, "%6i", parm->pbond[i].jb);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
	fprintf(stdout, "%6i", parm->pbond[i].icb);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
      }
      fprintf(stdout, "\n");
    }
  }

  /*
   *  Read in info for angles with hydrogen, post allocating memory 
   */
  parm->pangleH = safe_malloc( sizeof( ParmAngle ) * parm->NTHETH );

  if ( newparm == 0 || find_flag( parm->fp, "ANGLES_INC_HYDROGEN" ) ) {
    for (i=0; i < parm->NTHETH; i++) {
      parm->pangleH[i].it  = loadAndReturn12I6(parm->fp, 1, "angles with hydrogen, IT");
      parm->pangleH[i].jt  = loadAndReturn12I6(parm->fp, 1, "angles with hydrogen, JT");
      parm->pangleH[i].kt  = loadAndReturn12I6(parm->fp, 1, "angles with hydrogen, KT");
      parm->pangleH[i].ict = loadAndReturn12I6(parm->fp, 1, "angles with hydrogen, ICT");
    }
    loadAndReturn12I6(parm->fp, (parm->NTHETH ? 0 : -1), "");

    if (prnlev > 2) 
      fprintf(stdout, "Read in info for angles w/ hydrogen...\n");
    if (prnlev > 4) {
      j = 0;
      for (i=0; i < parm->NTHETH; i++) {
	fprintf(stdout, "%6i", parm->pangleH[i].it);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
	fprintf(stdout, "%6i", parm->pangleH[i].jt);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
	fprintf(stdout, "%6i", parm->pangleH[i].kt);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
	fprintf(stdout, "%6i", parm->pangleH[i].ict);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
      }
      fprintf(stdout, "\n");
    }
  }

  /*
   *  Read in info for angles without hydrogen, post allocating memory 
   */
  parm->pangle = safe_malloc( sizeof( ParmAngle ) * parm->MTHETS );
  if ( newparm == 0 || find_flag( parm->fp, "ANGLES_WITHOUT_HYDROGEN" ) ) {
    for (i=0; i < parm->MTHETS; i++) {
      parm->pangle[i].it  = loadAndReturn12I6(parm->fp, 1, "angles without hydrogen, IT");
      parm->pangle[i].jt  = loadAndReturn12I6(parm->fp, 1, "angles without hydrogen, JT");
      parm->pangle[i].kt  = loadAndReturn12I6(parm->fp, 1, "angles without hydrogen, KT");
      parm->pangle[i].ict = loadAndReturn12I6(parm->fp, 1, "angles without hydrogen, ICT");
    }
    loadAndReturn12I6(parm->fp, (parm->MTHETS ? 0 : -1), "");

    if (prnlev > 2) 
      fprintf(stdout, "Read in info for angles w/out hydrogen...\n");
    if (prnlev > 4) {
      j = 0;
      for (i=0; i < parm->MTHETS; i++) {
	fprintf(stdout, "%6i", parm->pangle[i].it);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
	fprintf(stdout, "%6i", parm->pangle[i].jt);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
	fprintf(stdout, "%6i", parm->pangle[i].kt);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
	fprintf(stdout, "%6i", parm->pangle[i].ict);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
      }
      fprintf(stdout, "\n");
    }
  }

  /*
   *  Read in info for dihedrals with hydrogen, post allocating memory 
   */
  if ( newparm == 0 || find_flag( parm->fp, "DIHEDRALS_INC_HYDROGEN" ) ) {
    parm->pdihedralH = safe_malloc( sizeof( ParmDihedral ) * parm->NPHIH );
    for (i=0; i < parm->NPHIH; i++) {
      parm->pdihedralH[i].ip  = loadAndReturn12I6(parm->fp, 1, "dihedral w/ hydrogen, IP");
      parm->pdihedralH[i].jp  = loadAndReturn12I6(parm->fp, 1, "dihedral w/ hydrogen, JP");
      parm->pdihedralH[i].kp  = loadAndReturn12I6(parm->fp, 1, "dihedral w/ hydrogen, KP");
      parm->pdihedralH[i].lp  = loadAndReturn12I6(parm->fp, 1, "dihedral w/ hydrogen, LP");
      parm->pdihedralH[i].icp = loadAndReturn12I6(parm->fp, 1, "dihedral w/ hydrogen, ICP");
    }
    loadAndReturn12I6(parm->fp, (parm->NPHIH ? 0 : -1), "");

    if (prnlev > 2) 
      fprintf(stdout, "Read in info for dihedrals w/ hydrogen...\n");
    if (prnlev > 4) {
      j = 0;
      for (i=0; i < parm->NPHIH; i++) {
	fprintf(stdout, "%6i", parm->pdihedralH[i].ip);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
	fprintf(stdout, "%6i", parm->pdihedralH[i].jp);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
	fprintf(stdout, "%6i", parm->pdihedralH[i].kp);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
	fprintf(stdout, "%6i", parm->pdihedralH[i].lp);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
	fprintf(stdout, "%6i", parm->pdihedralH[i].icp);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
      }
      fprintf(stdout, "\n");
    }
  }

  /*
   *  Read in info for dihedrals with hydrogen, post allocating memory 
   */
  parm->pdihedral = safe_malloc( sizeof( ParmDihedral ) * parm->MPHIA );
  if ( newparm == 0 || find_flag( parm->fp, "DIHEDRALS_WITHOUT_HYDROGEN" ) ) {
    for (i=0; i < parm->MPHIA; i++) {
      parm->pdihedral[i].ip  = loadAndReturn12I6(parm->fp, 1, "dihedral w/out hydrogen, IP");
      parm->pdihedral[i].jp  = loadAndReturn12I6(parm->fp, 1, "dihedral w/out hydrogen, JP");
      parm->pdihedral[i].kp  = loadAndReturn12I6(parm->fp, 1, "dihedral w/out hydrogen, KP");
      parm->pdihedral[i].lp  = loadAndReturn12I6(parm->fp, 1, "dihedral w/out hydrogen, LP");
      parm->pdihedral[i].icp = loadAndReturn12I6(parm->fp, 1, "dihedral w/out hydrogen, ICP");
    }
    loadAndReturn12I6(parm->fp, (parm->MPHIA ? 0 : -1), "");

    if (prnlev > 2) 
      fprintf(stdout, "Read in info for dihedrals w/out hydrogen...\n");
    if (prnlev > 4) {
      j = 0;
      for (i=0; i < parm->MPHIA; i++) {
	fprintf(stdout, "%6i", parm->pdihedral[i].ip);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
	fprintf(stdout, "%6i", parm->pdihedral[i].jp);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
	fprintf(stdout, "%6i", parm->pdihedral[i].kp);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
	fprintf(stdout, "%6i", parm->pdihedral[i].lp);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
	fprintf(stdout, "%6i", parm->pdihedral[i].icp);
	j++; if ( j%12 == 0 ) fprintf(stdout, "\n");
      }
      fprintf(stdout, "\n");
    }
  }


  /*
   *  Read in excluded atom list, post allocating memory 
   */
  parm->natex = safe_malloc( sizeof( int ) * parm->NEXT );
  if ( newparm == 0 || find_flag( parm->fp, "EXCLUDED_ATOMS_LIST" ) ) {
    for (i=0; i < parm->NEXT; i++) 
      parm->natex[i] = loadAndReturn12I6(parm->fp, 1, "natex");

    loadAndReturn12I6(parm->fp, 0, "");

    if (prnlev > 2) 
      fprintf(stdout, "Read in excluded atom list...\n");
    if (prnlev > 4) {
      for (i=0; i < parm->NEXT; i++) {
	fprintf(stdout, "%6i", parm->natex[i]);
	if ( (i+1)%12 == 0 ) fprintf(stdout, "\n");
      }
      fprintf(stdout, "\n");
    }
  }

  /*
   *  Read in h-bond r**12 and r**10 terms 
   *  (and hbcut) after alloc. memory 
   */
  parm->ag = safe_malloc( sizeof( double ) * parm->NHB );
  parm->bg = safe_malloc( sizeof( double ) * parm->NHB );
  parm->hbcut = safe_malloc( sizeof( double ) * parm->NHB );
  if ( newparm == 0 || find_flag( parm->fp, "HBOND_ACOEF" ) ) {
    for (i=0; i < parm->NHB; i++) 
      parmScanDouble(&parm->ag[i],    FALSE, "ag");
  }
  if ( newparm == 0 || find_flag( parm->fp, "HBOND_BCOEF" ) ) {
    for (i=0; i < parm->NHB; i++) 
      parmScanDouble(&parm->bg[i],    FALSE, "bg");
  }
  if( newparm == 0 || find_flag( parm->fp, "HBCUT" ) ) {
    for (i=0; i < parm->NHB; i++) 
      parmScanDouble(&parm->hbcut[i], FALSE, "hbcut");

    if (prnlev > 2) 
      fprintf(stdout, "Read in h-bond parameters: AG, BG, and HBCUT...\n");
    if (prnlev > 4) {
      for (i=0; i < parm->NHB; i++) {
	fprintf(stdout, "%16.8E", parm->ag[i]);
	if ( (i+1)%5 == 0 ) fprintf(stdout, "\n");
      }
      if (i%5) fprintf(stdout,"\n");
      for (i=0; i < parm->NHB; i++) {
	fprintf(stdout, "%16.8E", parm->bg[i]);
	if ( (i+1)%5 == 0 ) fprintf(stdout, "\n");
      }
      if (i%5) fprintf(stdout,"\n");
      for (i=0; i < parm->NHB; i++) {
	fprintf(stdout, "%16.8E", parm->hbcut[i]);
	if ( (i+1)%5 == 0 ) fprintf(stdout, "\n");
      }
      if (i%5) fprintf(stdout,"\n");
    }
  }


  /*
   * ISYMBL (TYPES)
   */
  if ( newparm == 0 || find_flag( parm->fp, "AMBER_ATOM_TYPE" ) ) {
    for (i=0; i < parm->NTOTAT; i++) {
      NAME_COPY(parm->fp, parm->atom[i].isymbl, "isymbl");
    }

    if (prnlev > 2) 
      fprintf(stdout, "Read in atomic symbols (types)...\n");
    if (prnlev > 4) {
      for (i=0; i < parm->NTOTAT; i++) {
	fprintf(stdout, "%4s", parm->atom[i].isymbl);
	if ( (i+1)%20 == 0 ) {
	  fprintf(stdout, "\n");
	}
      }
      fprintf(stdout, "\n");
    }
  }

  /*
   * TREE
   */
  if ( newparm == 0 || find_flag( parm->fp, "TREE_CHAIN_CLASSIFICATION" ) ) {
    for (i=0; i < parm->NTOTAT; i++) {
      NAME_COPY(parm->fp, parm->atom[i].itree, "itree");
    }


    if (prnlev > 2) 
      fprintf(stdout, "Read in tree information...\n");
    if (prnlev > 4) {
      for (i=0; i < parm->NTOTAT; i++) {
	fprintf(stdout, "%4s", parm->atom[i].itree);
	if ( (i+1)%20 == 0 ) {
	  fprintf(stdout, "\n");
	}
      }
      fprintf(stdout, "\n");
    }
  }

  /*
   * JOIN
   */
  if ( newparm == 0 || find_flag( parm->fp, "JOIN_ARRAY" ) ) {
    for (i=0; i < parm->NTOTAT; i++) 
      parm->atom[i].join = loadAndReturn12I6(parm->fp, 1, "join info");
    loadAndReturn12I6(parm->fp, 0, "");

    if (prnlev > 2)
      fprintf(stdout, "Read in the JOIN info...\n");
    if (prnlev > 4) {
      for (i=0; i < parm->NTOTAT; i++) {
	fprintf(stdout, "%6i", parm->atom[i].join);
	if ( (i+1)%12 == 0 ) fprintf(stdout, "\n");
      }
      fprintf(stdout, "\n");
    }
  }

  /*
   * IROTAT
   */
  if ( newparm == 0 || find_flag( parm->fp, "IROTAT" ) ) {
    for (i=0; i < parm->NTOTAT; i++) 
      parm->atom[i].irotat = loadAndReturn12I6(parm->fp, 1, "irotat info");
    loadAndReturn12I6(parm->fp, 0, "");

    if (prnlev > 2)
      fprintf(stdout, "Read in the IROTAT info...\n");
    if (prnlev > 4) {
      for (i=0; i < parm->NTOTAT; i++) {
	fprintf(stdout, "%6i", parm->atom[i].irotat);
	if ( (i+1)%12 == 0 ) fprintf(stdout, "\n");
      }
      fprintf(stdout, "\n");
    }
  }


  /*
   *  Get BOX information if IFBOX > 0 
   */
  if ( parm->IFBOX ) {
    if (prnlev > 2)
      fprintf(stdout, "Scanning Box\n");
    parm->box = safe_malloc( sizeof( Box ) );
    if ( newparm == 0 || find_flag( parm->fp, "SOLVENT_POINTERS" ) ) {

      parm->box->iptres = loadAndReturn12I6(parm->fp, 1, "box information, iptres");
      parm->box->nspm   = loadAndReturn12I6(parm->fp, 1, "box information, nspm");
      parm->box->nspsol = loadAndReturn12I6(parm->fp, 1, "box information, nspsol");
      loadAndReturn12I6(parm->fp, 0, "");
    }

    if ( newparm == 0 || find_flag( parm->fp, "ATOMS_PER_MOLECULE" ) ) {
      len = parm->box->nspm;
      parm->box->nsp = safe_malloc( sizeof( int ) * len );
      for (i=0; i < len; i++)
	parm->box->nsp[i] = loadAndReturn12I6(parm->fp, 1, "box information, nsp");
      loadAndReturn12I6(parm->fp, 0, "");
    }

    if ( newparm == 0 || find_flag( parm->fp, "BOX_DIMENSIONS" ) ) {
      parmScanDouble(&parm->box->beta,   FALSE, "box, beta");
      parmScanDouble(&parm->box->box[0], FALSE, "box, box[0]");
      parmScanDouble(&parm->box->box[1], FALSE, "box, box[1]");
      parmScanDouble(&parm->box->box[2], FALSE, "box, box[2]");
    }

    if (prnlev > 2) 
      fprintf(stdout, "Read in box information...\n");
    if (prnlev > 4) {
      fprintf(stdout, "%6d%6d%6d\n", 
	      parm->box->iptres, parm->box->nspm, parm->box->nspsol);

      for (i=0; i < len; i++) {
	fprintf(stdout, "%6i", parm->box->nsp[i]);
	if ( (i+1)%12 == 0 ) fprintf(stdout, "\n");
      }
      if (i%12) fprintf(stdout, "\n");

      fprintf(stdout, "%16.8E%16.8E%16.8E%16.8E\n", 
	      parm->box->beta, parm->box->box[0], parm->box->box[1], parm->box->box[2]);
    }
  }


  /*
   *  Get RADII and SCREEN terms if present
   */
  parm->RADII = 0;
  parm->RADIUS_SET[0] = (char) 0;
  if ( newparm && find_flag( parm->fp, "RADIUS_SET" ) ) {
    s_getline( parm->RADIUS_SET, 81, parm->fp );
  }
  if ( newparm && find_flag( parm->fp, "RADII" ) ) {
    parm->RADII = 1;
    for (i=0; i < parm->NTOTAT; i++) 
      parmScanDouble(&parm->atom[i].radii, FALSE, "radii");

    if (prnlev > 2)
      fprintf(stdout, "Read in radii...\n");
    if (prnlev > 4) {
      for (i=0; i < parm->NTOTAT; i++) {
	fprintf(stdout, "%16.8E", parm->atom[i].radii);
	if ( (i+1)%5 == 0 ) fprintf(stdout, "\n");
      }
      fprintf(stdout, "\n");
    }
  }

  parm->SCREEN = 0;
  if ( newparm && find_flag( parm->fp, "SCREEN" ) ) {
    parm->SCREEN = 1;
    for (i=0; i < parm->NTOTAT; i++) 
      parmScanDouble(&parm->atom[i].screen, FALSE, "screen");

    if (prnlev > 2)
      fprintf(stdout, "Read in screening parameters...\n");
    if (prnlev > 4) {
      for (i=0; i < parm->NTOTAT; i++) {
	fprintf(stdout, "%16.8E", parm->atom[i].screen);
	if ( (i+1)%5 == 0 ) fprintf(stdout, "\n");
      }
      fprintf(stdout, "\n");
    }
  }


  /*
   *  Scan cap info if IFCAP > 0 
   */
  if ( parm->IFCAP ) {
    fprintf(stdout, "Scanning Cap; this has not been debugged!\n");

    if ( newparm == 0 || find_flag( parm->fp, "CAP_INFO" ) ) {
      parm->cap = safe_malloc( sizeof( Cap ) );
      parm->cap->natcap = loadAndReturn12I6(parm->fp, 1, "cap information, cap");
      loadAndReturn12I6(parm->fp, 0, "");
    }

    if ( newparm == 0 || find_flag( parm->fp, "CAP_INFO2" ) ) {
      parmScanDouble(&parm->cap->cutcap, FALSE, "cap, cutcap");
      parmScanDouble(&parm->cap->xcap,   FALSE, "cap, xcap");
      parmScanDouble(&parm->cap->ycap,   FALSE, "cap, ycap");
      parmScanDouble(&parm->cap->zcap,   FALSE, "cap, zcap");
    }
    verbose("Read in CAP information...", "");
  }

  /*
   *  Scan perturbation info if IFPERT > 0 
   */
  if ( parm->IFPERT ) {
    fprintf(stdout, "Scanning Perturbation info...\n");
    parm->pert = safe_malloc( sizeof( Pert ) );
    if ( newparm == 0 || find_flag( parm->fp, "PERT_BOND_ATOMS" ) ) {
      if ( parm->NBPER ) {
	fprintf(stdout, "  Bonds...\n");
	parm->pert->ibper = safe_malloc( sizeof( int ) * parm->NBPER );
	parm->pert->jbper = safe_malloc( sizeof( int ) * parm->NBPER );
	for (i=0; i < parm->NBPER; i++) {
	  parm->pert->ibper[i] = loadAndReturn12I6(parm->fp, 1, "pert, ibper");
	  parm->pert->jbper[i] = loadAndReturn12I6(parm->fp, 1, "pert, jbper");
	}
	loadAndReturn12I6(parm->fp, 0, "");

	parm->pert->icbper = safe_malloc( sizeof( int ) * 2 * parm->NBPER );
	if ( newparm == 0 || find_flag( parm->fp, "PERT_BOND_PARAMS" ) ) {
	  for (i=0; i < 2 * parm->NBPER; i++) 
	    parm->pert->icbper[i] = loadAndReturn12I6(parm->fp, 1, "pert, icbper");
	  loadAndReturn12I6(parm->fp, 0, "");
	}
      }
    }

    if ( parm->NGPER ) {
      fprintf(stdout, "  Angles...\n");
      parm->pert->itper = safe_malloc( sizeof( int ) * parm->NGPER );
      parm->pert->jtper = safe_malloc( sizeof( int ) * parm->NGPER );
      parm->pert->ktper = safe_malloc( sizeof( int ) * parm->NGPER );
      if ( newparm == 0 || find_flag( parm->fp, "PERT_ANGLE_ATOMS" ) ) {
	for (i=0; i < parm->NGPER; i++) {
	  parm->pert->itper[i] = loadAndReturn12I6(parm->fp, 1, "pert, itper");
	  parm->pert->jtper[i] = loadAndReturn12I6(parm->fp, 1, "pert, jtper");
	  parm->pert->ktper[i] = loadAndReturn12I6(parm->fp, 1, "pert, ktper");
	}
	loadAndReturn12I6(parm->fp, 0, "");

	parm->pert->ictper = safe_malloc( sizeof( int ) * 2 * parm->NGPER );
	if ( newparm == 0 || find_flag( parm->fp, "PERT_ANGLE_PARAMS" ) ) {
	  for (i=0; i < 2 * parm->NGPER; i++)
	    parm->pert->ictper[i] = loadAndReturn12I6(parm->fp, 1, "pert, ictper");
	  loadAndReturn12I6(parm->fp, 0, "");
	}
      }
    }

    if ( parm->NDPER ) {
      fprintf(stdout, "  Dihedrals...\n");
      parm->pert->ipper = safe_malloc( sizeof( int ) * parm->NDPER ); 
      parm->pert->jpper = safe_malloc( sizeof( int ) * parm->NDPER );
      parm->pert->kpper = safe_malloc( sizeof( int ) * parm->NDPER );
      parm->pert->lpper = safe_malloc( sizeof( int ) * parm->NDPER );
      if ( newparm == 0 || find_flag( parm->fp, "PERT_DIHEDRAL_ATOMS" ) ) {
	for (i=0; i < parm->NDPER; i++) {
	  parm->pert->ipper[i] = loadAndReturn12I6(parm->fp, 1, "pert, ipper");
	  parm->pert->jpper[i] = loadAndReturn12I6(parm->fp, 1, "pert, jpper");
	  parm->pert->kpper[i] = loadAndReturn12I6(parm->fp, 1, "pert, kpper");
	  parm->pert->lpper[i] = loadAndReturn12I6(parm->fp, 1, "pert, lpper");
	}
	loadAndReturn12I6(parm->fp, 0, "");

	parm->pert->icpper = safe_malloc( sizeof( int ) * 2 * parm->NDPER );
	if ( newparm == 0 || find_flag( parm->fp, "PERT_DIHEDRAL_PARAMS" ) ) {
	  for (i=0; i < 2 * parm->NDPER; i++)
	    parm->pert->icpper[i] = loadAndReturn12I6(parm->fp, 1, "pert, icpper");
	  loadAndReturn12I6(parm->fp, 0, "");
	}
      }
    }

    fprintf(stdout, "  Extra perturbation information...\n");
    parm->pert->labper = safe_malloc( sizeof( Name ) * parm->NTOTRS );
    if ( newparm == 0 || find_flag( parm->fp, "PERT_RESIDUE_NAME" ) ) {
      for (i=0; i < parm->NTOTRS; i++) {
	NAME_COPY(parm->fp, parm->pert->labper[i], "pert, labper");
      }
    }

    parm->pert->igrper = safe_malloc( sizeof( Name ) * parm->NTOTAT );
    if ( newparm == 0 || find_flag( parm->fp, "PERT_ATOM_NAME" ) ) {
      for (i=0; i < parm->NTOTAT; i++) {
	NAME_COPY(parm->fp, parm->pert->igrper[i], "pert, igrper");
      }
    }

    parm->pert->ismper = safe_malloc( sizeof( Name ) * parm->NTOTAT );
    if ( newparm == 0 || find_flag( parm->fp, "PERT_ATOM_SYMBOL" ) ) {
      for (i=0; i < parm->NTOTAT; i++) {
	NAME_COPY(parm->fp, parm->pert->ismper[i], "pert, ismper");
      }
    }

    parm->pert->almper = safe_malloc( sizeof( double ) * parm->NTOTAT );
    if ( newparm == 0 || find_flag( parm->fp, "ALMPER" ) ) {
      for (i=0; i < parm->NTOTAT; i++)
	parmScanDouble(&parm->pert->almper[i], FALSE, "pert, almper");
    }

    parm->pert->iaper  = safe_malloc( sizeof( int )    * parm->NTOTAT );
    if ( newparm == 0 || find_flag( parm->fp, "IAPER" ) ) {
      for (i=0; i < parm->NTOTAT; i++)
	parm->pert->iaper[i] = loadAndReturn12I6(parm->fp, 1, "pert, iaper");
      loadAndReturn12I6(parm->fp, 0, "");
    }

    parm->pert->iacper = safe_malloc( sizeof( int )    * parm->NTOTAT );
    if ( newparm == 0 || find_flag( parm->fp, "PERT_ATOM_TYPE_INDEX" ) ) {
      for (i=0; i < parm->NTOTAT; i++)
	parm->pert->iacper[i] = loadAndReturn12I6(parm->fp, 1, "pert, iacper");
      loadAndReturn12I6(parm->fp, 0, "");
    }

    parm->pert->cgper  = safe_malloc( sizeof( double ) * parm->NTOTAT );
    if ( newparm == 0 || find_flag( parm->fp, "PERT_CHARGE" ) ) {
      for (i=0; i < parm->NTOTAT; i++)
	parmScanDouble(&parm->pert->cgper[i],  FALSE, "pert, cgper");
    }

    verbose("Read in perturbation information...", "");
  } /* END of: if ( parm->IFPERT ) */

  if (parm->JPARM == 1) {

    fprintf(stdout, "Reading in LES information: ");

    if ( newparm == 0 || find_flag( parm->fp, "LES_NTYP" ) ) {
      parm->nlestyp = loadAndReturn12I6(parm->fp, 1, "LES, nlestyp");
      loadAndReturn12I6(parm->fp, 0, "");
      fprintf(stdout, " %i copies found.\n", parm->nlestyp);

      parm->lestyp  = (int *) safe_malloc(sizeof(int) * parm->NTOTAT);
      parm->lesfac  = (double *) safe_malloc(sizeof(double) * parm->nlestyp*parm->nlestyp);
      parm->lescnum = (int *) safe_malloc(sizeof(int) * parm->NTOTAT);
      parm->lessubsp= (int *) safe_malloc(sizeof(int) * parm->NTOTAT);
    }

    if ( newparm == 0 || find_flag( parm->fp, "LES_TYPE" ) ) {
      for (i=0; i < parm->NTOTAT; i++)
	parm->lestyp[i] = loadAndReturn12I6(parm->fp, 1, "LES, lestyp");
      loadAndReturn12I6(parm->fp, 0, "");
    }

    if ( newparm == 0 || find_flag( parm->fp, "LES_FAC" ) ) {
      for (i=0; i < parm->nlestyp * parm->nlestyp; i++) 
	parmScanDouble(&parm->lesfac[i],  FALSE, "LES, lesfac");
    }

    if ( newparm == 0 || find_flag( parm->fp, "LES_CNUM" ) ) {
      for (i=0; i < parm->NTOTAT; i++)
	parm->lescnum[i] = loadAndReturn12I6(parm->fp, 1, "LES, lescnum");
      loadAndReturn12I6(parm->fp, 0, "");
    }

    if ( newparm == 0 || find_flag( parm->fp, "LES_ID" ) ) {
      for (i=0; i < parm->NTOTAT; i++)
	parm->lessubsp[i] = loadAndReturn12I6(parm->fp, 1, "LES, lessubsp");
      loadAndReturn12I6(parm->fp, 0, "");
    }
  }

  /* END of reading in parameter file */

  /* fix up bond, angle, and dihedrals */
  parm->bond = safe_malloc( sizeof( Bond ) * 
			   (parm->NBONH + parm->MBONA + parm->NBPER) );
  for (i=0; i < parm->NBONH; i++) {
    parm->bond[i].atom[0] = unObfuscateAtom(parm->pbondH[i].ib);
    parm->bond[i].atom[1] = unObfuscateAtom(parm->pbondH[i].jb);
    parm->bond[i].rk = parm->rk[parm->pbondH[i].icb-1];
    parm->bond[i].req = parm->req[parm->pbondH[i].icb-1];
    parm->bond[i].scale = 1.0;
  }
  for (i=0; i < parm->MBONA; i++) {
    parm->bond[i+parm->NBONH].atom[0] = 
      unObfuscateAtom(parm->pbond[i].ib);
    parm->bond[i+parm->NBONH].atom[1] = 
      unObfuscateAtom(parm->pbond[i].jb);
    parm->bond[i+parm->NBONH].rk = parm->rk[parm->pbond[i].icb-1];
    parm->bond[i+parm->NBONH].req = parm->req[parm->pbond[i].icb-1];
    parm->bond[i+parm->NBONH].scale = 1.0;
  }
/*
  for (i=0; i < parm->NBPER; i++) {
    parm->bond[i+parm->NBONH+parm->MBONA].atom[0] =
      unObfuscateAtom(parm->pert->ibper[i]);
    parm->bond[i+parm->NBONH+parm->MBONA].atom[1] =
      unObfuscateAtom(parm->pert->jbper[i]);
    parm->bond[i+parm->NBONH+parm->MBONA].rk =
      parm->rk[parm->pert->icbper[

*/
  parm->angle = safe_malloc( sizeof( Angle ) *
			    (parm->NTHETH + parm->MTHETS) );
  for (i=0; i < parm->NTHETH; i++) {
    parm->angle[i].atom[0] = unObfuscateAtom(parm->pangleH[i].it);
    parm->angle[i].atom[1] = unObfuscateAtom(parm->pangleH[i].jt);
    parm->angle[i].atom[2] = unObfuscateAtom(parm->pangleH[i].kt);
    parm->angle[i].tk = parm->tk[parm->pangleH[i].ict-1];
    parm->angle[i].teq = parm->teq[parm->pangleH[i].ict-1];
    parm->angle[i].scale = 1.0;                          
  }
  for (i=0; i < parm->MTHETS; i++) {
    parm->angle[i+parm->NTHETH].atom[0] = 
      unObfuscateAtom(parm->pangle[i].it);
    parm->angle[i+parm->NTHETH].atom[1] = 
      unObfuscateAtom(parm->pangle[i].jt);
    parm->angle[i+parm->NTHETH].atom[2] = 
      unObfuscateAtom(parm->pangle[i].kt);
    parm->angle[i+parm->NTHETH].tk = parm->tk[parm->pangle[i].ict-1];   
    parm->angle[i+parm->NTHETH].teq = parm->teq[parm->pangle[i].ict-1];
    parm->angle[i+parm->NTHETH].scale = 1.0;                          
  }

  parm->dihedral = safe_malloc( sizeof( Dihedral ) *
			    (parm->NPHIH + parm->MPHIA) );
  for (i=0; i < parm->NPHIH; i++) {
    parm->dihedral[i].atom[0] = unObfuscateAtom(parm->pdihedralH[i].ip);
    parm->dihedral[i].atom[1] = unObfuscateAtom(parm->pdihedralH[i].jp);
    parm->dihedral[i].atom[2] = unObfuscateAtom(parm->pdihedralH[i].kp);
    parm->dihedral[i].atom[3] = unObfuscateAtom(parm->pdihedralH[i].lp);
    parm->dihedral[i].pk = parm->pk[parm->pdihedralH[i].icp-1];
    parm->dihedral[i].pn = parm->pn[parm->pdihedralH[i].icp-1];
    parm->dihedral[i].phase = parm->phase[parm->pdihedralH[i].icp-1];
    parm->dihedral[i].scale = 1.0;                          
  }
  for (i=0; i < parm->MPHIA; i++) {
    parm->dihedral[i+parm->NPHIH].atom[0] = 
      unObfuscateAtom(parm->pdihedral[i].ip);
    parm->dihedral[i+parm->NPHIH].atom[1] = 
      unObfuscateAtom(parm->pdihedral[i].jp);        
    parm->dihedral[i+parm->NPHIH].atom[2] = 
      unObfuscateAtom(parm->pdihedral[i].kp);
    parm->dihedral[i+parm->NPHIH].atom[3] = 
      unObfuscateAtom(parm->pdihedral[i].lp);
    parm->dihedral[i+parm->NPHIH].pk = parm->pk[parm->pdihedral[i].icp-1];   
    parm->dihedral[i+parm->NPHIH].pn = parm->pn[parm->pdihedral[i].icp-1];   
    parm->dihedral[i+parm->NPHIH].phase = 
      parm->phase[parm->pdihedral[i].icp-1];   
    parm->dihedral[i+parm->NPHIH].scale = 1.0;                          
  }

  if ( isErrorParm == 0 ) {
    if (prnlev > 2) 
      fprintf(stdout, "Successfully completed readParm.\n");
  } else
    fprintf(stdout, "There were errors upon reading parm file!!!\n");
  
  safe_free(buffer);
}

   void            /* global, referenced by dispatch() */
writeParmOLD( FILE *fp ) 
{
  int count, i, len;

  /* write title */
  parm->title[80] = (char) 0;
  fprintf(fp, "%s%s", parm->title,
	  (strchr(parm->title, '\n') == NULL ? "\n" : "" ));
  fflush(fp);

  /* write integer control variables */
  fprintf(fp, "%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n",
	  parm->NTOTAT, parm->NTYPES, parm->NBONH, 
	  parm->NBONA,  parm->NTHETH, parm->NTHETA, 
	  parm->NPHIH,  parm->NPHIA,  parm->JHPARM, 
	  parm->JPARM,  parm->NEXT,   parm->NTOTRS);

  fprintf(fp, "%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d%6d\n",
	  parm->MBONA,  parm->MTHETS, parm->MPHIA,  
	  parm->MUMBND, parm->MUMANG, parm->MPTRA,
	  parm->NATYP,  parm->NHB,    parm->IFPERT,
	  parm->NBPER,  parm->NGPER,  parm->NDPER);

  fprintf(fp, "%6d%6d%6d%6d%6d%6d\n",
	  parm->MBPER,  parm->MGPER,  parm->MDPER,  
	  parm->IFBOX,  parm->NMXRS,  parm->IFCAP);

  /* write atom names */
  write20A4(fp, parm->atom[i].igraph, 0, parm->NTOTAT);

  /* write charges */
  write5E16(fp, parm->atom[i].chrg, 0, parm->NTOTAT);

  /* write masses */
  write5E16(fp, parm->atom[i].amass, 0, parm->NTOTAT);

  /* write out iac */
  write12I6(fp, parm->atom[i].iac, 0, parm->NTOTAT);

  /* write out numex */
  write12I6(fp, parm->atom[i].numex, 0, parm->NTOTAT);
  
  fflush(fp);

  /* write out nno */
  len = parm->NTYPES; 
  len = len * len;
  write12I6(fp, parm->nno[i], 0, len);

  /* write out residue labels */
  write20A4(fp, parm->residue[i].labres, 0, parm->NTOTRS);

  /* write out ipres */
  write12I6(fp, parm->residue[i].ipres, 0, parm->NTOTRS);

  /* write out bond, angle, dihedral parameters */
  write5E16(fp, parm->rk[i], 0, parm->MUMBND);
  write5E16(fp, parm->req[i], 0, parm->MUMBND);

  write5E16(fp, parm->tk[i], 0, parm->MUMANG);
  write5E16(fp, parm->teq[i], 0, parm->MUMANG);

  write5E16(fp, parm->pk[i], 0, parm->MPTRA);
  write5E16(fp, parm->pn[i], 0, parm->MPTRA);
  write5E16(fp, parm->phase[i], 0, parm->MPTRA);


  fflush(fp);

  /* write out solty */
  write5E16(fp, parm->solty[i], 0, parm->NATYP);

  /* write out cn1 and cn2 */
  len = parm->NTYPES;
  len = len * (len + 1) / 2;
  write5E16(fp, parm->cn1[i], 0, len);
  write5E16(fp, parm->cn2[i], 0, len);

  fflush(fp);

  /* write bonds */
  for (i=0; i < parm->NBONH; i++) {
    fprintf(fp, "%6d%6d%6d", parm->pbondH[i].ib, parm->pbondH[i].jb,
	    parm->pbondH[i].icb);
    if ( (i+1)%4 == 0) fprintf(fp, "\n");
  } 
  if ( i%4 || parm->NBONH == 0 ) fprintf(fp, "\n");

  for (i=0; i < parm->MBONA; i++) {
    fprintf(fp, "%6d%6d%6d", parm->pbond[i].ib, parm->pbond[i].jb,
	    parm->pbond[i].icb);
    if ((i+1)%4 == 0) fprintf(fp, "\n");
  }
  if ( i%4 || parm->MBONA == 0 ) fprintf(fp, "\n");

  fflush(fp);

  /* write angles */
  for (i=0; i < parm->NTHETH; i++) {
    fprintf(fp, "%6d%6d%6d%6d", parm->pangleH[i].it, parm->pangleH[i].jt, 
	    parm->pangleH[i].kt, parm->pangleH[i].ict);
    if ( (i+1) % 3 == 0 ) fprintf(fp, "\n");
  }
  if ( i%3 || parm->NTHETH == 0 ) fprintf(fp, "\n");

  for (i=0; i < parm->MTHETS; i++) {
    fprintf(fp, "%6d%6d%6d%6d", parm->pangle[i].it, parm->pangle[i].jt,
	    parm->pangle[i].kt, parm->pangle[i].ict);
    if ( (i+1) % 3 == 0 ) fprintf(fp, "\n");
  }
  if ( i%3 || parm->MTHETS == 0 ) fprintf(fp, "\n");

  fflush(fp);

  /* write dihedrals */
  count = 0;
  for (i=0; i < parm->NPHIH; i++) {
    fprintf(fp, "%6d", parm->pdihedralH[i].ip);
    count++; if ( count == 12 ) { count = 0; fprintf(fp, "\n"); }
    fprintf(fp, "%6d", parm->pdihedralH[i].jp);
    count++; if ( count == 12 ) { count = 0; fprintf(fp, "\n"); }
    fprintf(fp, "%6d", parm->pdihedralH[i].kp);
    count++; if ( count == 12 ) { count = 0; fprintf(fp, "\n"); }
    fprintf(fp, "%6d", parm->pdihedralH[i].lp);
    count++; if ( count == 12 ) { count = 0; fprintf(fp, "\n"); }
    fprintf(fp, "%6d", parm->pdihedralH[i].icp);
    count++; if ( count == 12 ) { count = 0; fprintf(fp, "\n"); }
  } 
  if ( (count % 12) || (parm->NPHIH == 0) ) fprintf(fp, "\n");
  fflush(fp);

  count = 0;
  for (i=0; i < parm->MPHIA; i++) {
    fprintf(fp, "%6d", parm->pdihedral[i].ip);
    count++; if ( count == 12 ) { count = 0; fprintf(fp, "\n"); }
    fprintf(fp, "%6d", parm->pdihedral[i].jp);
    count++; if ( count == 12 ) { count = 0; fprintf(fp, "\n"); }
    fprintf(fp, "%6d", parm->pdihedral[i].kp);
    count++; if ( count == 12 ) { count = 0; fprintf(fp, "\n"); }
    fprintf(fp, "%6d", parm->pdihedral[i].lp);
    count++; if ( count == 12 ) { count = 0; fprintf(fp, "\n"); }
    fprintf(fp, "%6d", parm->pdihedral[i].icp);
    count++; if ( count == 12 ) { count = 0; fprintf(fp, "\n"); }
  }
  if ( (count % 12) || (parm->MPHIA == 0) ) fprintf(fp, "\n");

  fflush(fp);

  /* write out exclued atom list */
  write12I6(fp, parm->natex[i], 0, parm->NEXT);

  /* write out h-bond r**12 and r**10 terms (and hbcut) */
  write5E16(fp, parm->ag[i], 0, parm->NHB);
  write5E16(fp, parm->bg[i], 0, parm->NHB);
  write5E16(fp, parm->hbcut[i], 0, parm->NHB);

  /* write out atomic and tree symbols, tree join, and rotat info */
  write20A4(fp, parm->atom[i].isymbl, 0, parm->NTOTAT);
  write20A4(fp, parm->atom[i].itree, 0, parm->NTOTAT);
  write12I6(fp, parm->atom[i].join, 0, parm->NTOTAT);
  write12I6(fp, parm->atom[i].irotat, 0, parm->NTOTAT);

  /* write box info if IFBOX > 0 */
  if ( parm->IFBOX ) {
    fprintf(fp, "%6d%6d%6d\n", parm->box->iptres, parm->box->nspm,
	    parm->box->nspsol);
    write12I6(fp, parm->box->nsp[i], 0, parm->box->nspm);
    fprintf(fp, "%16.8E%16.8E%16.8E%16.8E\n",
	    parm->box->beta, parm->box->box[0],
	    parm->box->box[1], parm->box->box[2]);
  }

  /* write cap info if IFCAP > 0 */
  if ( parm->IFCAP ) {
    fprintf(fp, "%6d\n%16.8E%16.8E%16.8E%16.8E\n",
	    parm->cap->natcap, parm->cap->cutcap, parm->cap->xcap,
	    parm->cap->ycap,   parm->cap->zcap);
  }

  /* write perturbation info if IFPERT > 0 */
  if ( parm->IFPERT ) {
    /* perturbed bonds and parameters */
    for (i=0; i < parm->NBPER; i++) {
      fprintf(fp, "%6d%6d", parm->pert->ibper[i], parm->pert->jbper[i]);
      if ( (i > 0) && ((i+1)%6 == 0) ) fprintf(fp, "\n");
    }
    if ( i % 6 != 0 ) fprintf(fp, "\n");
    write12I6(fp, parm->pert->icbper[i], 0, 2*parm->NBPER);

    /* perturbed angles and parameters */
    for (i=0; i < parm->NGPER; i++) {
      fprintf(fp, "%6d%6d%6d", parm->pert->itper[i],
	      parm->pert->jtper[i], parm->pert->ktper[i]);
      if ( (i > 0) && ((i+1)%4 == 0) ) fprintf(fp, "\n");
    }
    if ( i % 4 != 0 ) fprintf(fp, "\n");
    write12I6(fp, parm->pert->ictper[i], 0, 2*parm->NGPER);

    /* perturbed dihedrals and parameters */
    for (i=0; i < parm->NDPER; i++) {
      fprintf(fp, "%6d%6d%6d%6d", parm->pert->ipper[i],
	      parm->pert->jpper[i], parm->pert->kpper[i], 
	      parm->pert->lpper[i]);
      if ( (i > 0) && ( (i+1)%3 == 0) ) fprintf(fp, "\n");
    }
    if ( i+1 % 3 != 0 ) fprintf(fp, "\n");
    write12I6(fp, parm->pert->icpper[i], 0, 2*parm->NDPER);

    write20A4(fp, parm->pert->labper[i], 0, parm->NTOTRS);

    write20A4(fp, parm->pert->igrper[i], 0, parm->NTOTAT);

    write20A4(fp, parm->pert->ismper[i], 0, parm->NTOTAT);

    write5E16(fp, parm->pert->almper[i], 0, parm->NTOTAT);

    write12I6(fp, parm->pert->iaper[i], 0, parm->NTOTAT);

    write12I6(fp, parm->pert->iacper[i], 0, parm->NTOTAT);

    write5E16(fp, parm->pert->cgper[i], 0, parm->NTOTAT);
  }

  /* write out LES stuff */
  if (parm->JPARM == 1) {

    fprintf(fp, "%6i\n", parm->nlestyp);
    write12I6(fp, parm->lestyp[i], 0, parm->NTOTAT);
    write5E16(fp, parm->lesfac[i], 0, parm->nlestyp*parm->nlestyp);
    write12I6(fp, parm->lescnum[i], 0, parm->NTOTAT);
    write12I6(fp, parm->lessubsp[i], 0, parm->NTOTAT);
  }
  fprintf(stdout, "Done...\n");
}

   void            /* global, referenced by dispatch() */
writeParm( FILE *fp, int new ) 
{
  int count, i, len;
  char buffer[81];
  char timebuffer[20];

  time_t curtime;
  struct tm *loctime;

  if (new == 0) {
    writeParmOLD( fp );
    return;
  }


  curtime = time(NULL);
  loctime = localtime(&curtime);
  strftime(timebuffer, 20, "%m/%d/%y  %H:%M:%S", loctime);

  /*
   *  Dump version info and write title 
   */

  sprintf(buffer, "%%VERSION  VERSION_STAMP = V0001.000  DATE = %s CREATED BY PTRAJ", timebuffer);
  fprintf(fp, "%-80s\n", buffer);

  fprintf(fp, "%-80s\n", "%FLAG TITLE");
  fprintf(fp, "%-80s\n", "%FORMAT(20a4)");
  parm->title[80] = (char) 0;
  fprintf(fp, "%s%s", parm->title,
	  (strchr(parm->title, '\n') == NULL ? "\n" : "" ));
  fflush(fp);
                                                                                
  fprintf(fp, "%-80s\n", "%FLAG POINTERS");
  fprintf(fp, "%-80s\n", "%FORMAT(10I8)");

  /*
   *  write integer control variables 
   */
  fprintf(fp, "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
	  parm->NTOTAT, parm->NTYPES, parm->NBONH, parm->NBONA,  parm->NTHETH, 
	  parm->NTHETA, parm->NPHIH,  parm->NPHIA, parm->JHPARM, parm->JPARM);

  fprintf(fp, "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
	  parm->NEXT,   parm->NTOTRS, parm->MBONA,  parm->MTHETS, parm->MPHIA,  
	  parm->MUMBND, parm->MUMANG, parm->MPTRA,  parm->NATYP,  parm->NHB);

  fprintf(fp, "%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d\n",
	  parm->IFPERT, parm->NBPER,  parm->NGPER,  parm->NDPER,  parm->MBPER,  
	  parm->MGPER,  parm->MDPER,  parm->IFBOX,  parm->NMXRS,  parm->IFCAP);

  fprintf(fp, "%8d\n",
	  parm->NUMEXTRA);

  /*
   *  write atom names 
   */
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG ATOM_NAME", "%FORMAT(20a4)");                                                                   
  write20A4(fp, parm->atom[i].igraph, 0, parm->NTOTAT);

  /* 
   *  write charges 
   */
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG CHARGE", "%FORMAT(5E16.8)");
  write5E16(fp, parm->atom[i].chrg, 0, parm->NTOTAT);

  /*
   *  write masses 
   */
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG MASS", "%FORMAT(5E16.8)");
  write5E16(fp, parm->atom[i].amass, 0, parm->NTOTAT);

  /*
   *  write out iac 
   */
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG ATOM_TYPE_INDEX", "%FORMAT(10I8)");
  write10I8(fp, parm->atom[i].iac, 0, parm->NTOTAT);

  /* 
   *  write out numex 
   */
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG NUMBER_EXCLUDED_ATOMS", "%FORMAT(10I8)");
  write10I8(fp, parm->atom[i].numex, 0, parm->NTOTAT);
  
  /*
   *  write out nno 
   */
  len = parm->NTYPES; 
  len = len * len;
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG NONBONDED_PARM_INDEX", "%FORMAT(10I8)");
  write10I8(fp, parm->nno[i], 0, len);

  /*
   *  write out residue labels 
   */
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG RESIDUE_LABEL", "%FORMAT(20a4)");
  write20A4(fp, parm->residue[i].labres, 0, parm->NTOTRS);

  /*
   *  write out ipres 
   */
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG RESIDUE_POINTER", "%FORMAT(10I8)");
  write10I8(fp, parm->residue[i].ipres, 0, parm->NTOTRS);

  /*
   *  write out bond, angle, dihedral parameters 
   */
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG BOND_FORCE_CONSTANT", "%FORMAT(5E16.8)");
  write5E16(fp, parm->rk[i], 0, parm->MUMBND);
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG BOND_EQUIL_VALUE", "%FORMAT(5E16.8)");
  write5E16(fp, parm->req[i], 0, parm->MUMBND);

  fprintf(fp, "%-80s\n%-80s\n", "%FLAG ANGLE_FORCE_CONSTANT", "%FORMAT(5E16.8)");
  write5E16(fp, parm->tk[i], 0, parm->MUMANG);
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG ANGLE_EQUIL_VALUE", "%FORMAT(5E16.8)");
  write5E16(fp, parm->teq[i], 0, parm->MUMANG);

  fprintf(fp, "%-80s\n%-80s\n", "%FLAG DIHEDRAL_FORCE_CONSTANT", "%FORMAT(5E16.8)");
  write5E16(fp, parm->pk[i], 0, parm->MPTRA);
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG DIHEDRAL_PERIODICITY", "%FORMAT(5E16.8)");
  write5E16(fp, parm->pn[i], 0, parm->MPTRA);
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG DIHEDRAL_PHASE", "%FORMAT(5E16.8)");
  write5E16(fp, parm->phase[i], 0, parm->MPTRA);

  /*
   *  write out solty 
   */
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG SOLTY", "%FORMAT(5E16.8)");
  write5E16(fp, parm->solty[i], 0, parm->NATYP);

  /*
   *  write out cn1 and cn2 
   */
  len = parm->NTYPES;
  len = len * (len + 1) / 2;
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG LENNARD_JONES_ACOEF", "%FORMAT(5E16.8)");
  write5E16(fp, parm->cn1[i], 0, len);
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG LENNARD_JONES_BCOEF", "%FORMAT(5E16.8)");
  write5E16(fp, parm->cn2[i], 0, len);

  /*
   *  write bonds
   */
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG BONDS_INC_HYDROGEN", "%FORMAT(10I8)");
  count = 0;
  for (i=0; i < parm->NBONH; i++) {
    fprintf(fp, "%8d", parm->pbondH[i].ib);
    count++; if (count%10 == 0) fprintf(fp, "\n");
    fprintf(fp, "%8d", parm->pbondH[i].jb);
    count++; if (count%10 == 0) fprintf(fp, "\n");
    fprintf(fp, "%8d", parm->pbondH[i].icb);
    count++; if (count%10 == 0) fprintf(fp, "\n");
  } 
  if ( count%10 || parm->NBONH == 0 ) fprintf(fp, "\n");

  fprintf(fp, "%-80s\n%-80s\n", "%FLAG BONDS_WITHOUT_HYDROGEN", "%FORMAT(10I8)");
  count = 0;
  for (i=0; i < parm->MBONA; i++) {
    fprintf(fp, "%8d", parm->pbond[i].ib);
    count++; if (count%10 == 0) fprintf(fp, "\n");
    fprintf(fp, "%8d", parm->pbond[i].jb);
    count++; if (count%10 == 0) fprintf(fp, "\n");
    fprintf(fp, "%8d", parm->pbond[i].icb);
    count++; if (count%10 == 0) fprintf(fp, "\n");
  }
  if ( count%10 || parm->MBONA == 0 ) fprintf(fp, "\n");

  /*
   *  write angles 
   */
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG ANGLES_INC_HYDROGEN", "%FORMAT(10I8)");
  count = 0;
  for (i=0; i < parm->NTHETH; i++) {
    fprintf(fp, "%8d", parm->pangleH[i].it);
    count++; if (count%10 == 0) fprintf(fp, "\n");
    fprintf(fp, "%8d", parm->pangleH[i].jt);
    count++; if (count%10 == 0) fprintf(fp, "\n");
    fprintf(fp, "%8d", parm->pangleH[i].kt);
    count++; if (count%10 == 0) fprintf(fp, "\n");
    fprintf(fp, "%8d", parm->pangleH[i].ict);
    count++; if (count%10 == 0) fprintf(fp, "\n");
  }
  if ( count%10 || parm->NTHETH == 0 ) fprintf(fp, "\n");

  fprintf(fp, "%-80s\n%-80s\n", "%FLAG ANGLES_WITHOUT_HYDROGEN", "%FORMAT(10I8)");
  count = 0;
  for (i=0; i < parm->MTHETS; i++) {
    fprintf(fp, "%8d", parm->pangle[i].it);
    count++; if (count%10 == 0) fprintf(fp, "\n");
    fprintf(fp, "%8d", parm->pangle[i].jt);
    count++; if (count%10 == 0) fprintf(fp, "\n");
    fprintf(fp, "%8d", parm->pangle[i].kt);
    count++; if (count%10 == 0) fprintf(fp, "\n");
    fprintf(fp, "%8d", parm->pangle[i].ict);
    count++; if (count%10 == 0) fprintf(fp, "\n");
  }
  if ( count%10 || parm->MTHETS == 0 ) fprintf(fp, "\n");

  /*
   *  write dihedrals 
   */
  count = 0;
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG DIHEDRALS_INC_HYDROGEN", "%FORMAT(10I8)");
  for (i=0; i < parm->NPHIH; i++) {
    fprintf(fp, "%8d", parm->pdihedralH[i].ip);
    count++; if ( count == 10 ) { count = 0; fprintf(fp, "\n"); }
    fprintf(fp, "%8d", parm->pdihedralH[i].jp);
    count++; if ( count == 10 ) { count = 0; fprintf(fp, "\n"); }
    fprintf(fp, "%8d", parm->pdihedralH[i].kp);
    count++; if ( count == 10 ) { count = 0; fprintf(fp, "\n"); }
    fprintf(fp, "%8d", parm->pdihedralH[i].lp);
    count++; if ( count == 10 ) { count = 0; fprintf(fp, "\n"); }
    fprintf(fp, "%8d", parm->pdihedralH[i].icp);
    count++; if ( count == 10 ) { count = 0; fprintf(fp, "\n"); }
  } 
  if ( (count % 10) || (parm->NPHIH == 0) ) fprintf(fp, "\n");
  fflush(fp);

  count = 0;
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG DIHEDRALS_WITHOUT_HYDROGEN", "%FORMAT(10I8)");
  for (i=0; i < parm->MPHIA; i++) {
    fprintf(fp, "%8d", parm->pdihedral[i].ip);
    count++; if ( count == 10 ) { count = 0; fprintf(fp, "\n"); }
    fprintf(fp, "%8d", parm->pdihedral[i].jp);
    count++; if ( count == 10 ) { count = 0; fprintf(fp, "\n"); }
    fprintf(fp, "%8d", parm->pdihedral[i].kp);
    count++; if ( count == 10 ) { count = 0; fprintf(fp, "\n"); }
    fprintf(fp, "%8d", parm->pdihedral[i].lp);
    count++; if ( count == 10 ) { count = 0; fprintf(fp, "\n"); }
    fprintf(fp, "%8d", parm->pdihedral[i].icp);
    count++; if ( count == 10 ) { count = 0; fprintf(fp, "\n"); }
  }
  if ( (count % 10) || (parm->MPHIA == 0) ) fprintf(fp, "\n");

  /*
   *  write out excluded atom list 
   */
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG EXCLUDED_ATOMS_LIST", "%FORMAT(10I8)");
  write10I8(fp, parm->natex[i], 0, parm->NEXT);

  /*
   *  write out h-bond r**12 and r**10 terms (and hbcut) 
   */
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG HBOND_ACOEF", "%FORMAT(5E16.8)");
  write5E16(fp, parm->ag[i], 0, parm->NHB);
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG HBOND_BCOEF", "%FORMAT(5E16.8)");                                                                 
  write5E16(fp, parm->bg[i], 0, parm->NHB);
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG HBCUT", "%FORMAT(5E16.8)");
  write5E16(fp, parm->hbcut[i], 0, parm->NHB);


  /*
   *  write out atomic and tree symbols, tree join, and rotat info 
   */
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG AMBER_ATOM_TYPE", "%FORMAT(20a4)");
  write20A4(fp, parm->atom[i].isymbl, 0, parm->NTOTAT);
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG TREE_CHAIN_CLASSIFICATION", "%FORMAT(20a4)");
  write20A4(fp, parm->atom[i].itree, 0, parm->NTOTAT);
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG JOIN_ARRAY", "%FORMAT(10I8)");
  write10I8(fp, parm->atom[i].join, 0, parm->NTOTAT);
  fprintf(fp, "%-80s\n%-80s\n", "%FLAG IROTAT", "%FORMAT(10I8)");
  write10I8(fp, parm->atom[i].irotat, 0, parm->NTOTAT);

  /*
   *  write box info if IFBOX > 0 
   */
  if ( parm->IFBOX ) {


    fprintf(fp, "%-80s\n%-80s\n", "%FLAG SOLVENT_POINTERS", "%FORMAT(3I8)");
    fprintf(fp, "%8d%8d%8d\n", parm->box->iptres, parm->box->nspm,
	    parm->box->nspsol);
    fprintf(fp, "%-80s\n%-80s\n", "%FLAG ATOMS_PER_MOLECULE", "%FORMAT(10I8)");
    write10I8(fp, parm->box->nsp[i], 0, parm->box->nspm);
    fprintf(fp, "%-80s\n%-80s\n", "%FLAG BOX_DIMENSIONS", "%FORMAT(5E16.8)");
    fprintf(fp, "%16.8E%16.8E%16.8E%16.8E\n",
	    parm->box->beta, parm->box->box[0],
	    parm->box->box[1], parm->box->box[2]);
  }



  /* write cap info if IFCAP > 0 */
  if ( parm->IFCAP ) {
    fprintf(fp, "%-80s\n%-80s\n", "%FLAG CAP_INFO", "%FORMAT(10I8)");                                                                 
    fprintf(fp, "%8d\n", parm->cap->natcap);
    fprintf(fp, "%-80s\n%-80s\n", "%FLAG CAP_INFO2", "%FORMAT(5E16.8)");                                                                 
    fprintf(fp, "%16.8E%16.8E%16.8E%16.8E\n",
	    parm->cap->cutcap, parm->cap->xcap,
	    parm->cap->ycap,   parm->cap->zcap);
  }

  /*
   *  write radii and screen
   */
  if (parm->RADIUS_SET[0] != (char) 0)
    fprintf(fp, "%-80s\n%-80s\n%-80s\n", "%FLAG RADIUS_SET", "%FORMAT(1a80)", parm->RADIUS_SET);
  if (parm->RADII) {
    fprintf(fp, "%-80s\n%-80s\n", "%FLAG RADII", "%FORMAT(5E16.8)");                                                                 
    write5E16(fp, parm->atom[i].radii, 0, parm->NTOTAT);
  }
  if (parm->SCREEN) {
    fprintf(fp, "%-80s\n%-80s\n", "%FLAG SCREEN", "%FORMAT(5E16.8)");                                                                 
    write5E16(fp, parm->atom[i].screen, 0, parm->NTOTAT);
  }


  /*
   *  write perturbation info if IFPERT > 0 
   */
  if ( parm->IFPERT ) {
    /*
     *  perturbed bonds and parameters 
     */
    for (i=0; i < parm->NBPER; i++) {
      fprintf(fp, "%8d%8d", parm->pert->ibper[i], parm->pert->jbper[i]);
      if ( (i > 0) && ((i+1)%5 == 0) ) fprintf(fp, "\n");
    }
    if ( i % 5 != 0 ) fprintf(fp, "\n");
    write10I8(fp, parm->pert->icbper[i], 0, 2*parm->NBPER);

    /* perturbed angles and parameters */
    count = 0;
    for (i=0; i < parm->NGPER; i++) {
      fprintf(fp, "%8d", parm->pert->itper[i]);
      count++; if ( count == 10 ) { count = 0; fprintf(fp, "\n"); }
      fprintf(fp, "%8d", parm->pert->jtper[i]);
      count++; if ( count == 10 ) { count = 0; fprintf(fp, "\n"); }
      fprintf(fp, "%8d", parm->pert->ktper[i]);
      count++; if ( count == 10 ) { count = 0; fprintf(fp, "\n"); }
    }
    if ( count % 10 != 0 ) fprintf(fp, "\n");
    write10I8(fp, parm->pert->ictper[i], 0, 2*parm->NGPER);

    /* perturbed dihedrals and parameters */
    count = 0;
    for (i=0; i < parm->NDPER; i++) {
      fprintf(fp, "%8d", parm->pert->ipper[i]);
      count++; if ( count == 10 ) { count = 0; fprintf(fp, "\n"); }
      fprintf(fp, "%8d", parm->pert->jpper[i]);
      count++; if ( count == 10 ) { count = 0; fprintf(fp, "\n"); }
      fprintf(fp, "%8d", parm->pert->kpper[i]);
      count++; if ( count == 10 ) { count = 0; fprintf(fp, "\n"); }
      fprintf(fp, "%8d", parm->pert->lpper[i]);
      count++; if ( count == 10 ) { count = 0; fprintf(fp, "\n"); }
    }
    if ( count % 10 != 0 ) fprintf(fp, "\n");
    write10I8(fp, parm->pert->icpper[i], 0, 2*parm->NDPER);

    write20A4(fp, parm->pert->labper[i], 0, parm->NTOTRS);

    write20A4(fp, parm->pert->igrper[i], 0, parm->NTOTAT);

    write20A4(fp, parm->pert->ismper[i], 0, parm->NTOTAT);

    write5E16(fp, parm->pert->almper[i], 0, parm->NTOTAT);

    write10I8(fp, parm->pert->iaper[i], 0, parm->NTOTAT);

    write10I8(fp, parm->pert->iacper[i], 0, parm->NTOTAT);

    write5E16(fp, parm->pert->cgper[i], 0, parm->NTOTAT);
  }

  /* write out LES stuff */
  if (parm->JPARM == 1) {
    fprintf(fp, "%-80s\n%-80s\n", "%FLAG LES_NTYP", "%FORMAT(10I8)");                                                                 
    fprintf(fp, "%8i\n", parm->nlestyp);
    fprintf(fp, "%-80s\n%-80s\n", "%FLAG LES_TYPE", "%FORMAT(10I8)");                                                                 
    write10I8(fp, parm->lestyp[i], 0, parm->NTOTAT);
    fprintf(fp, "%-80s\n%-80s\n", "%FLAG LES_FAC", "%FORMAT(5E16)");                                                                 
    write5E16(fp, parm->lesfac[i], 0, parm->nlestyp*parm->nlestyp);
    fprintf(fp, "%-80s\n%-80s\n", "%FLAG LES_CNUM", "%FORMAT(10I8)");                                                                 
    write10I8(fp, parm->lescnum[i], 0, parm->NTOTAT);
    fprintf(fp, "%-80s\n%-80s\n", "%FLAG LES_ID", "%FORMAT(10I8)");                                                                 
    write10I8(fp, parm->lessubsp[i], 0, parm->NTOTAT);
  }
  fprintf(stdout, "Done...\n");
}

   void
dumpBond(char type, int bondnum, int atom1, int atom2, 
	 double rk, double req, double scale)
{
  char bondNames[50];
  
  sprintf(bondNames, "  :%i@%s :%i@%s", 
	  parm->atom[atom1-1].res+1, parm->atom[atom1-1].igraph,
	  parm->atom[atom2-1].res+1, parm->atom[atom2-1].igraph);
  
  fprintf(stdout, "%c %6i:  %6.2f %6.3f %s (%i,%i)\n", 
	 type, bondnum, rk, req, bondNames, atom1, atom2);
}

#define BOND_PRINT_STRING \
"\n    Bond     Kb     Req       atom names   (numbers)\n"

   void            /* global, referenced by dispatch() */
printBonds(char *atomSpec )
{
  int i, total, nbonh;
  int *mask = NULL;

  if (atomSpec) 
    mask = returnAtomMask(atomSpec);

  total = parm->NBONH + parm->MBONA;
  nbonh = parm->NBONH;
  if ( total ) {
    fprintf(stdout, BOND_PRINT_STRING);
    for (i=0; i < total; i++) {
      if (i == nbonh) fprintf(stdout, "\n");
      if ( mask == NULL || mask[ parm->bond[i].atom[0]-1 ] ||
                           mask[ parm->bond[i].atom[1]-1 ] )
      dumpBond(' ', i+1, parm->bond[i].atom[0], parm->bond[i].atom[1], 
	       parm->bond[i].rk,
	       parm->bond[i].req,
	       parm->bond[i].scale);
    }
  } else
    fprintf(stdout, " No bonds present in this parameter file\n");
  if (atomSpec) {
    safe_free(mask);
  }
}
#undef BOND_PRINT_STRING


#define BOND_PRINT_STRING \
"\n   Bond\nlambda       Kb     Req       atom names   (numbers)\n"


   void            /* global, referenced by dispatch() */
printPerturbedBonds(char *atomSpec )
{
  int i;
  int *mask = NULL;

  if ( parm->NBPER ) {

    if (atomSpec) {
      mask = returnAtomMask(atomSpec);
    }

    fprintf(stdout, BOND_PRINT_STRING);
    for (i=0; i < parm->NBPER; i++) {
      if (mask == NULL ||
	  mask[ parm->pert->ibper[i]/3 ] ||
	  mask[ parm->pert->jbper[i]/3 ]) {
	dumpBond('0', i+1, parm->pert->ibper[i]/3+1, parm->pert->jbper[i]/3+1, 
		 parm->rk[parm->pert->icbper[i]-1],
		 parm->req[parm->pert->icbper[i]-1],
		 1.0);
	dumpBond('1', i+1, parm->pert->ibper[i]/3+1, parm->pert->jbper[i]/3+1, 
		 parm->rk[parm->pert->icbper[i+parm->NBPER]-1],
		 parm->req[parm->pert->icbper[i+parm->NBPER]-1],
		 1.0);
	fprintf(stdout, "\n");
      }
    }
    if (atomSpec)
      safe_free(mask);
    
  } else
    fprintf(stdout, " No perturbed bonds\n");
}
#undef BOND_PRINT_STRING

   void
dumpAngle(char type, int anglenum,
	  int atom1, int atom2, int atom3, 
	  double tk, double teq, double scale)
{
  char angleNames[50];
  double conversion = 180 / (4.0 * (double) atan( (double) 1.0 ));

  sprintf(angleNames, "  :%i@%s :%i@%s :%i@%s", 
	  parm->atom[atom1-1].res+1, parm->atom[atom1-1].igraph,
	  parm->atom[atom2-1].res+1, parm->atom[atom2-1].igraph,
	  parm->atom[atom3-1].res+1, parm->atom[atom3-1].igraph);

  fprintf(stdout, "%c %5i:  %6.3f %6.2f %s (%i,%i,%i)\n", 
	  type, anglenum, tk, teq * conversion, angleNames, atom1, atom2, atom3);
}

#define ANGLE_PRINT_STRING \
"\n  Angle   Kthet  degrees        atom names        (numbers)\n"



   void            /* global, referenced by dispatch() */
printAngles(char *atomSpec)
{
  int i, total,ntheth;
  int *mask = NULL;

  total = parm->NTHETH + parm->MTHETS;
  ntheth = parm->NTHETH;

  if ( total ) {

    if (atomSpec)
      mask = returnAtomMask(atomSpec);

    fprintf(stdout, ANGLE_PRINT_STRING);
    for (i=0; i < total; i++) {
      if (i == ntheth) fprintf(stdout, "\n");
      if (mask == NULL ||
	  mask[ parm->angle[i].atom[0]-1 ] ||
	  mask[ parm->angle[i].atom[1]-1 ] ||
	  mask[ parm->angle[i].atom[2]-1 ]) {
	dumpAngle(' ', i+1,
		  parm->angle[i].atom[0],
		  parm->angle[i].atom[1],
		  parm->angle[i].atom[2],
		  parm->angle[i].tk,
		  parm->angle[i].teq,
		  parm->angle[i].scale);
      }
    }
    if (atomSpec) safe_free(mask);
  } else
    fprintf(stdout, " No angles defined in this parameter file\n");
}
#undef ANGLE_PRINT_STRING

#define ANGLE_PRINT_STRING \
"\n   Angle\nlambda     Kthet  Req             atom names      (numbers)\n"

   void            /* global, referenced by dispatch() */
printPerturbedAngles(char *atomSpec)
{
  int i;
  int *mask = NULL;

  if ( parm->NGPER ) {
    if (atomSpec)
      mask = returnAtomMask(atomSpec);

    fprintf(stdout, ANGLE_PRINT_STRING);
    for (i=0; i < parm->NGPER; i++) {
      if (mask == NULL ||
	  mask[ parm->pert->itper[i]/3 ] ||
	  mask[ parm->pert->jtper[i]/3 ] ||
	  mask[ parm->pert->ktper[i]/3 ] ) {
	dumpAngle('0', i+1,
		  parm->pert->itper[i]/3+1, 
		  parm->pert->jtper[i]/3+1, 
		  parm->pert->ktper[i]/3+1, 
		  parm->tk[parm->pert->ictper[i]-1],
		  parm->teq[parm->pert->ictper[i]-1],
		  1.0);
	dumpAngle('1', i+1,
		  parm->pert->itper[i]/3+1, 
		  parm->pert->jtper[i]/3+1, 
		  parm->pert->ktper[i]/3+1, 
		  parm->tk[parm->pert->ictper[i+parm->NGPER]-1],
		  parm->teq[parm->pert->ictper[i+parm->NGPER]-1],
		  1.0);
	fprintf(stdout, "\n");
      }
    }
    if (atomSpec)
      safe_free(mask);
  } else
    fprintf(stdout, " No perturbed angles\n");
}
#undef ANGLE_PRINT_STRING


   void
dumpDihedral(FILE *outFile, char type, int anglenum,
	     int atom1, int atom2, int atom3, int atom4,
	     double pk, double pn, double phase, double scale)
{
  char dihedralNames[80];

  if ( anglenum < 0 )
    fprintf(outFile, "%c          %6.3f  %4.2f %4.1f\n",
	   type, pk, phase, pn);
  else {
    sprintf(dihedralNames, "  :%i@%s :%i@%s :%i@%s :%i@%s", 
	    parm->atom[atom1-1].res+1, parm->atom[atom1-1].igraph,
	    parm->atom[atom2-1].res+1, parm->atom[atom2-1].igraph,
	    parm->atom[atom3-1].res+1, parm->atom[atom3-1].igraph,
	    parm->atom[atom4-1].res+1, parm->atom[atom4-1].igraph);
    
    fprintf(outFile, "%c %6i:  %6.3f  %4.2f %4.1f ",
	    type, anglenum, pk, phase, pn);
    fprintf(outFile, " %s (%i,%i,%i,%i)\n",
	    dihedralNames, atom1, atom2, atom3, atom4);
  }
}

#define DIHEDRAL_PRINT_STRING \
"\nDihedral    pk     phase pn                atoms\n"


   void            /* global */
PrintDihedral(int dihedralNumber, FILE *outFile, int *mask)
{
  int i, j, total, nphih;
  int atom1, atom2, atom3, atom4;
  int params;
  char type;
  total = parm->NPHIH + parm->MPHIA;
  nphih = parm->NPHIH;

  if ( total ) {
    if ( dihedralNumber == -1 )
      fprintf(outFile, DIHEDRAL_PRINT_STRING);
    for (i= (dihedralNumber == -1 ? 0 : dihedralNumber-1); 
	 i < (dihedralNumber == -1 ? total : dihedralNumber); i++) {
      
      atom1 = parm->dihedral[i].atom[0];      
      atom2 = parm->dihedral[i].atom[1];     
      atom3 = parm->dihedral[i].atom[2];
      atom4 = parm->dihedral[i].atom[3];     

      if (mask == NULL || 
	  mask[atom1-1] ||
	  mask[atom2-1] ||
	  mask[iabs(atom3)-1] ||
	  mask[iabs(atom4)-1]) {
	  

	/* NOTE: what is all this type crap?  
	 *   if atom3 < 0 --> end atoms shouldn't have nonbonds between them
	 *   if atom4 < 0 --> improper dihedral
	 */
	type = ' ';
	if ( atom3 < 0 && atom4 < 0 ) 
	  type = 'B';
	else if ( atom3 < 0 )
	  type = 'E';
	else if ( atom4 < 0 )
	  type = 'I';
	j = (i<nphih ? i+1: i-nphih+1);
	if ( i == nphih ) fprintf(outFile, "\n");
	dumpDihedral(outFile, type, i+1, 
		     atom1, atom2, iabs(atom3), iabs(atom4),
		     parm->dihedral[i].pk, 
		     parm->dihedral[i].pn,
		     parm->dihedral[i].phase, 
		     parm->dihedral[i].scale);
	if ( parm->dihedral[i].pn < 0 ) {
	  if ( i < nphih ) {
	    params = parm->pdihedralH[i].icp - 1;
	  } else
	    params = parm->pdihedral[i-nphih].icp - 1;
	  j = 0;
	  do {
	    j++;
	    dumpDihedral(outFile, type, -1, -1, -1, -1, -1, 
			 parm->pk[params+j],
			 parm->pn[params+j],
			 parm->phase[params+j],
			 parm->dihedral[i].scale);
	  } while ( parm->pn[params+j] < 0 );
	}
      }
    }
  } else
    fprintf(outFile, " No dihedrals defined in this parameter file\n");
}




   void            /* global */
dumpDihedrals()
{
  int i;

  fprintf(stdout, "Dihedral Type     PEAK    PN   PHASE\n");
  for (i=0; i < parm->MPTRA; i++) {
    fprintf(stdout, "   %5i    %10.4f  %10.4f  %10.4f\n",
	   i,
	   parm->pk[i],
	   parm->pn[i],
	   parm->phase[i]);
  }
}




   void            /* global, referenced by dispatch() */
printDihedrals(char *atomSpec)
{
  int *mask = NULL;
  if (atomSpec) {
    mask = returnAtomMask(atomSpec);
  }
  PrintDihedral(-1, stdout, mask);
  if (atomSpec) {
    safe_free(mask);
  }
}



   void            /* global, referenced by dispatch() */
printPerturbedDihedrals(char *atomSpec)
{
  int i, atom1, atom2, atom3, atom4;
  char type;
  int j, params;
  int *mask = NULL;

  if ( parm->NDPER ) {

    if (atomSpec)
      mask = returnAtomMask(atomSpec);

    fprintf(stdout, DIHEDRAL_PRINT_STRING);
    for (i=0; i < parm->NDPER; i++) {

      atom1 = parm->pert->ipper[i]/3+1;
      atom2 = parm->pert->jpper[i]/3+1;
      if ( parm->pert->kpper[i] < 0 ) 
	atom3 = (-parm->pert->kpper[i])/3+1;
      else
	atom3 = parm->pert->kpper[i]/3+1;

      if ( parm->pert->lpper[i] < 0 ) 
	atom4 = (-parm->pert->lpper[i])/3+1;
      else
	atom4 = parm->pert->lpper[i]/3+1;

      if (mask == NULL ||
	  mask[atom1-1] ||
	  mask[atom2-1] ||
	  mask[atom3-1] ||
	  mask[atom4-1]) {
	type = ' ';
	if ( parm->pert->kpper[i] < 0 && parm->pert->lpper[i] < 0 ) 
	  type = 'B';
	else if ( parm->pert->kpper[i] < 0 )
	  type = 'E';
	else if ( parm->pert->lpper[i] < 0 )
	  type = 'I';
	if ( type == ' ' ) type = '0';
	dumpDihedral(stdout, type, i+1, atom1, atom2, atom3, atom4,
		     parm->pk[parm->pert->icpper[i]-1],
		     parm->pn[parm->pert->icpper[i]-1],
		     parm->phase[parm->pert->icpper[i]-1],
		     1.0);
	if ( parm->pn[parm->pert->icpper[i]-1] < 0 ) {
	  params = parm->pert->icpper[i] - 1;
	  j = 1;
	  do {
	    dumpDihedral(stdout, type, -1, -1, -1, -1, -1, 
			 parm->pk[params+j],
			 parm->pn[params+j],
			 parm->phase[params+j],
			 1.0);
	  } while ( parm->pn[j++] < 0 );
	}

	if ( type == '0' ) type = '1';
	dumpDihedral(stdout, type, i+1, atom1, atom2, iabs(atom3), iabs(atom4),
		     parm->pk[parm->pert->icpper[i+parm->NDPER]-1],
		     parm->pn[parm->pert->icpper[i+parm->NDPER]-1],
		     parm->phase[parm->pert->icpper[i+parm->NDPER]-1],
		     1.0);

	if ( parm->pn[parm->pert->icpper[i+parm->NDPER]-1] < 0 ) {
	  params = parm->pert->icpper[i+parm->NDPER] - 1;
	  j = 1;
	  do {
	    dumpDihedral(stdout, type, -1, -1, -1, -1, -1, 
			 parm->pk[params+j],
			 parm->pn[params+j],
			 parm->phase[params+j],
			 1.0);
	  } while ( parm->pn[j++] < 0 );
	}
	fprintf(stdout, "\n");
      }
      if (atomSpec)
	safe_free(mask);
    }
    
  } else
    fprintf(stdout, " No perturbed dihedrals\n");
}
#undef DIHEDRAL_PRINT_STRING

   void            /* global, referenced by dispatch() */
printExcluded()
{
  int i, j, exclude, counter;

  fprintf(stdout, "The excluded atom list, NEXT is %i\n", parm->NEXT);

  exclude = 0;
  for (i=0; i < parm->NTOTAT; i++) {
    counter = 0;
    fprintf(stdout, "%4s : ", parm->atom[i].igraph);
    for (j = exclude; 
	 j < exclude + parm->atom[i].numex; 
	 j++) {
      counter++;
      if (parm->natex[j] != 0)
	fprintf(stdout, " %4s", parm->atom[parm->natex[j]-1].igraph);
      if ( counter % 10 == 0 && i != parm->NTOTAT-1 ) fprintf(stdout,"\n       ");
    }
    fprintf(stdout, "\n");
    exclude += parm->atom[i].numex;
  }
/*
	 j < (i == parm->NTOTAT-1 : parm->NEXT ? parm->atom[i+1].numex-1);
*/
}

   void            /* global, referenced by dispatch() */
printLJ()
{
  int i, j, k;
  int ntypes;
  int index;
  int atom1, atom2;
  ntypes = parm->NTYPES;

  fprintf(stdout, "  Lennard-Jones\n");
  fprintf(stdout, "  Types             A             C\n");
  for (i=0; i < ntypes; i++) {
    for (j=i; j < ntypes; j++) {

      index = parm->nno[ntypes*i+j]-1;
      atom1 = parm->NTOTAT;
      atom2 = parm->NTOTAT;
      for (k=0; k < parm->NTOTAT; k++) {
	if (i+1 == parm->atom[k].iac)
	  atom1 = k;
	if (j+1 == parm->atom[k].iac)
	  atom2 = k;
      }
      if (atom1 == parm->NTOTAT || atom2 == parm->NTOTAT) {
	if ( parm->IFPERT ) {
	  for (k=0; k < parm->NTOTAT; k++) {
	    if (i+1 == parm->pert->iacper[k] && atom1 == parm->NTOTAT)
	      atom1 = -k;
	    if (j+1 == parm->pert->iacper[k] && atom2 == parm->NTOTAT)
	      atom2 = -k;
	  }
	  if ( atom1 == 0 ) atom1 = -parm->NTOTAT;
	  if ( atom2 == 0 ) atom2 = -parm->NTOTAT;
	}
      }
      if (atom1 != parm->NTOTAT && atom2 != parm->NTOTAT) {
	fprintf(stdout, "%c  %s %s  %12.2lf  %12.2lf\n",
               (index < 0 ? 'H' : ' '),
	       (atom1 < 0 ? 
		   (atom1 == -parm->NTOTAT ?
                      parm->pert->ismper[0] :
		      parm->pert->ismper[-atom1]) :
		   parm->atom[atom1].isymbl),
	       (atom2 < 0 ? 
		   (atom2 == -parm->NTOTAT ?
		      parm->pert->ismper[0] :
		      parm->pert->ismper[-atom2]) :
		   parm->atom[atom2].isymbl),
               (index < 0 ? parm->ag[-index] : parm->cn1[index]),
               (index < 0 ? parm->bg[-index] : parm->cn2[index]));
      } else {
	warning("printLJ", "Atom type not found for current index %i",
		index);
      }
    }
  }
}


/*
 *  Dump out each of the atom types and its symbol (including the
 *  perturbation types
 */

   void
printAtomTypes()
{
  int i, j, index;
  int atom;
  char *type;
  double rstar, eps, a, c;
  

  if (parm->NTYPES > 0) {
    fprintf(stdout, "NOTE: if either A or C is zero (*), we cannot infer\n");
    fprintf(stdout, "the value of r* or epsilon and assume zero...\n\n");
    fprintf(stdout, "  Type     r*      eps\n");
  }
  for (i=0; i < parm->NTYPES; i++) {

    index = parm->nno[parm->NTYPES*i+i]-1;
    atom = parm->NTOTAT;
    for (j = 0; j < parm->NTOTAT; j++) {
      if (i+1 == parm->atom[j].iac)
	atom = j;
    }
    
    if ( atom == parm->NTOTAT && parm->IFPERT ) {
      /*
       *  this implies that the atom is in the perturbed list
       */

      for (j = 0; j < parm->NTOTAT; j++) {
	if (i+1 == parm->pert->iacper[j])
	  atom = -j;
      }
    }

    if ( atom != parm->NTOTAT ) {

      if ( atom < 0 )
	
	type = (char *) parm->pert->ismper[-atom];

      else

	type = (char *) parm->atom[atom].isymbl;

      if ( atom < 0 ) {
	index = parm->nno[ parm->NTYPES * (parm->pert->iacper[-atom]-1) +
			  parm->pert->iacper[-atom]-1 ] - 1;
      } else {
	index = parm->nno[ parm->NTYPES * (parm->atom[atom].iac-1) + 
			  parm->atom[atom].iac-1]-1;
      }

      if (index < 0) {
	a = parm->ag[-index];
	c = parm->bg[-index];
      } else {
	a = parm->cn1[index];
	c = parm->cn2[index];
      }
      if ( a == 0 || c == 0 ) {
	rstar = 0.0;
	eps = 0.0;
      } else {
	rstar = pow( (double) 2.0 * a / c, (double) 1.0/6.0) / 2.0;
	eps =(c * c) / (4.0 * a);
      }


      fprintf(stdout, "%c  %s %7.4f %7.4f\n", 
	      ( (a == 0 || c == 0) ? '*' : ' ' ), type, rstar, eps);
    }
  }
}


   int
getAtomTypes(Name atomName, 
	     double *rstar, 
	     double *eps,
	     int mode)
{
  int i, index;
  double a, c;
  int atom;

  /*
   *  mode = 0   --> match by name (igraph)
   *  otherwise  --> match by type (isymbl)
   */


  atom = -1;
  for (i=0; i < parm->NTOTAT; i++) {
  
    if ( strcmp((mode == 0 ? parm->atom[i].igraph : parm->atom[i].isymbl),
		atomName) == 0 ) {
      atom = i;
      index = parm->nno[parm->NTYPES * (parm->atom[atom].iac-1) + 
			parm->atom[atom].iac-1] - 1;

      break;
    }
  }
  if ( atom == -1 && parm->IFPERT ) {

    for (i = 0; i < parm->NTOTAT; i++) {
	
      if ( strcmp((mode == 0 ? parm->pert->igrper[i] : 
		   parm->pert->ismper[i]),
		  atomName) == 0 ) {
	atom = i;
	index = parm->nno[parm->NTYPES * (parm->pert->iacper[atom]-1) + 
			  parm->pert->iacper[atom]-1] - 1;
	break;
      }
    }
  }

  if ( atom < 0 ) {
    warning("getAtomType()", "Atom %s (%s) not found\n", 
	    (mode == 0 ? "name" : "type"), atomName);
    return 0;
  }

  if (index < 0) {
    a = parm->ag[-index];
    c = parm->bg[-index];
  } else {
    a = parm->cn1[index];
    c = parm->cn2[index];
  }
  if ( a == 0 || c == 0 ) {
    *rstar = 0.0;
    *eps = 0.0;
  } else {
    *rstar = pow( (double) 2.0 * a / c, (double) 1.0/6.0) / 2.0;
    *eps =(c * c) / (4.0 * a);
  }

  return 1;
}


   int
getAtomCharge(Name atomName,
	      int residue,
	      double *charge, 
	      int mode)
{
  int i;
  int atom;
  int res;

  /*
   *  mode = 0   --> match by name (igraph)
   *  otherwise  --> match by type (isymbl)
   */


  atom = -1;
  res = residue - 1;

  for (i=parm->residue[res].ipres-1; i < parm->residue[res+1].ipres; i++) {
    if ( strcmp((mode == 0 ? parm->atom[i].igraph : parm->atom[i].isymbl),
		atomName) == 0 && parm->atom[i].res == res ) {
      *charge = parm->atom[i].chrg;
      atom = i;
      break;
    }
  }

  if ( atom == -1 && parm->IFPERT ) {
    for (i=parm->residue[res].ipres-1; i < parm->residue[res+1].ipres; i++) {
      if ( strcmp((mode == 0 ? parm->pert->igrper[i] : 
		   parm->pert->ismper[i]),
		  atomName) == 0 && parm->atom[i].res == res ) {
	*charge = parm->pert->cgper[i];
	atom = i;
	break;
      }
    }
  }

  if ( atom < 0 ) {
    warning("getAtomType()", "Atom %s (%s), residue %i not found\n", 
	    (mode == 0 ? "name" : "type"), atomName, residue);
    return 0;
  } else
    return 1;
}

   int
getAtomRadii(Name atomName,
	     int residue,
	     double *radii, 
	     int mode)
{
  int i;
  int atom;
  int res;

  /*
   *  mode = 0   --> match by name (igraph)
   *  otherwise  --> match by type (isymbl)
   */

  atom = -1;
  res = residue - 1;

  for (i=parm->residue[res].ipres-1; i < parm->residue[res+1].ipres; i++) {
    if ( strcmp((mode == 0 ? parm->atom[i].igraph : parm->atom[i].isymbl),
		atomName) == 0 && parm->atom[i].res == res ) {
      if (parm->RADII) 
	*radii = parm->atom[i].radii;
      else
	*radii = 0.0;
      atom = i;
      break;
    }
  }

  if ( atom < 0 ) {
    warning("getAtomRadii()", "Atom %s (%s), residue %i not found\n", 
	    (mode == 0 ? "name" : "type"), atomName, residue);
    return 0;
  } else
    return 1;
}


   void
checkVelocity()
{
  char buffer[BUFFER_SIZE];
  FILE *fpin;
  char *filename;
  int maxi;
  double velocity, avgvel, maxvel;
  double *x,*y,*z;
  double box[6];
  int i, set, status;
  char atom1[20];
  coordinateInfo *trajInfo;

  printf("This routine will read in a velocity file and check\n");
  printf("for the largest entries...\n\n");
  filename = promptToOpenFile(&fpin, "", "r",
                  "Input the name of an AMBER trajectory file: ");
  safe_fclose(fpin);

  x = (double *) safe_malloc( sizeof(double) * parm->NTOTAT );
  y = (double *) safe_malloc( sizeof(double) * parm->NTOTAT );
  z = (double *) safe_malloc( sizeof(double) * parm->NTOTAT );

  /*
   *  Preprocess AMBER trajectory file by calling checkCoordinates
   */
  trajInfo = checkCoordinates(filename, parm->NTOTAT);

  if ( trajInfo->type != COORD_AMBER_TRAJECTORY )
    error("checkVelocity()", "File: %s is not an AMBER trajectory\n", filename);

  openFile(&fpin, filename, "r");

  if ( fgets(buffer, BUFFER_SIZE, fpin) == NULL )
    error("checkVelocity()", "fgets returned NULL\n");

  /* LOOP over coordinate sets until status = 0
   * status should be set by the program responsible for reading
   * in a particular coordinate set
   */
  status = 1;
  for (set = 1; status ; set++) {
    
    status = readAmberTrajectory(fpin, parm->NTOTAT, x, y, z, box, set, trajInfo);

    if (status == 0) break;

    avgvel = 0.0;
    maxvel = 0.0;
    for (i=0; i < parm->NTOTAT; i++) {

      velocity = x[i]*x[i]+y[i]*y[i]+z[i]*z[i];
      velocity = sqrt ( velocity );
      avgvel += velocity;

      if ( maxvel < velocity ) {
	maxi = i;
	maxvel = velocity;
      }
    }

    avgvel = avgvel / parm->NTOTAT;

    printf("------------------------------------------------\n");
    printf("Set %i, average velocity is %8.3f\n", set, avgvel);
    rdparmPrintAtom(maxi, atom1);
    printf("Maximum occured for atom %i [%s] (%8.3f)\n", maxi+1, atom1, maxvel);
    for (i=0; i < parm->NTOTAT; i++) {

      velocity = x[i]*x[i]+y[i]*y[i]+z[i]*z[i];
      velocity = sqrt( velocity );
      if (velocity > (maxvel * 0.95)) {
	rdparmPrintAtom(i, atom1);
	printf("   Atom %i [%s] has velocity (%8.3f) > 0.95 times the average\n",
	       i+1, atom1, velocity);
      }
    }
    
  }
  set -= 1;
  printf("Successfully read in a total of %5i coordinate sets\n", set);

}


   void
testit()
{
  pdb_record *pdb;

  pdb = parmToPdb(parm, NULL);	

  savePdb(stdout, pdb);

}




   void            /* global, referenced by dispatch() */
printDelphiCharge(char *atomSpec)
{
  int i;
  int *mask = NULL;

  if (atomSpec) {
    mask = returnAtomMask(atomSpec);
  }

  fprintf(stdout, "!Amber rdparm charges for %s\n", parm->filename);
  fprintf(stdout, "atom__resnumbc_charge_\n");

  for (i=0; i < parm->NTOTAT; i++) {
    if (mask == NULL || mask[i]) {

      fprintf(stdout, "%s  %s%3i  %7.4f\n",
                parm->atom[i].igraph,
		parm->residue[parm->atom[i].res].labres,
		parm->atom[i].res+1,
                parm->atom[i].chrg / CHARGE_TO_KCALS);
    }
  }
  if (atomSpec) 
    safe_free(mask);
}




   void            /* global, referenced by dispatch() */
printAtomInfo(char *atomSpec)
{
  int i;
  int *mask = NULL;

  if (atomSpec) {
    mask = returnAtomMask(atomSpec);
  }

  if (parm == NULL) return;

  if (parm->IFPERT) {
    fprintf(stdout, "Perturbation is ON                                  ");
    fprintf(stdout, "Perturbed atom information\n");
  }
  fprintf(stdout, "Number:  Atom   Charge Mass ( Residue ) Type Tree    %s\n",
	 ( parm->IFPERT ? " Atom   Charge ( Residue ) Type" : ""));
  for (i=0; i < parm->NTOTAT; i++) {
    if (mask == NULL || mask[i]) {
      if ( parm->IFPERT && parm->pert->iaper[i] == 1 ) {
	/* atom is being perturbed */
	fprintf(stdout, "%6i:  %s %8.5lf %4.1lf  (%4i:%s)  %2s  %1s  %s %7.4lf  (%4i:%s)  %2s\n",
		i+1, parm->atom[i].igraph,
		parm->atom[i].chrg / CHARGE_TO_KCALS,
		parm->atom[i].amass,
		parm->atom[i].res+1,
		parm->residue[parm->atom[i].res].labres,
		parm->atom[i].isymbl, parm->atom[i].itree,
		parm->pert->igrper[i],
		parm->pert->cgper[i] / CHARGE_TO_KCALS,
		parm->atom[i].res+1,
		parm->pert->labper[parm->atom[i].res],
		parm->pert->ismper[i]);
      
      } else {
	fprintf(stdout, "%6i:  %s %8.5lf %4.1lf (%4i:%s)  %2s  %1s\n",
		i+1, parm->atom[i].igraph,
		parm->atom[i].chrg / CHARGE_TO_KCALS,
		parm->atom[i].amass,
		parm->atom[i].res+1,
		parm->residue[parm->atom[i].res].labres,
		parm->atom[i].isymbl, parm->atom[i].itree);
      }
    }
  }
  if (atomSpec) 
    safe_free(mask);
}


   void            /* global, referenced by dispatch() */
countAtoms(char *atomSpec)
{
  int i, count;
  int *mask = NULL;

  if (atomSpec) mask = returnAtomMask(atomSpec);
  count = 0;

  for (i=0; i < parm->NTOTAT; i++)
    if (mask[i] == 1) count++;

  printf("Mask %s represents %i atoms: ", atomSpec, count);

  if (atomSpec) safe_free(mask);
}


   void            /* global, referenced by dispatch() */
chargeOnAtoms(char *atomSpec)
{
  int i, count;
  int *mask = NULL;
  double charge;

  charge = 0.0;
  count = 0;

  if (atomSpec) mask = returnAtomMask(atomSpec);

  for (i=0; i < parm->NTOTAT; i++)
    if (mask[i] == 1) {
      count++;
      charge += parm->atom[i].chrg;
    }

  printf("Mask %s represents %i atoms and as charge %10.3f: \n",
	 atomSpec, count, charge / CHARGE_TO_KCALS);

  if (atomSpec) safe_free(mask);
}




   void
deleteBAD(char *type, int i)
{
  int j, shift;

  switch ( type[0] ) {

  case 'B':
  case 'b':

    fprintf(stdout, "\nAttempting to delete bond %i\n\n", i);

    if ( i > 0 && i <= parm->NBONH ) {  
      /* case: 1 -> NBONH; bond with hydrogen */

      if ( i != parm->NBONH ) {
	shift = sizeof(ParmBond) * (parm->NBONH - i);
	shiftArray( parm->pbondH+i, parm->pbondH+i-1, shift );
      }

      shift = sizeof(Bond) * (parm->NBONH + parm->MBONA - i);
      shiftArray( parm->bond+i, parm->bond+i-1, shift );

      parm->NBONH -= 1;

    } else if ( i > parm->NBONH && i < parm->NBONH + parm->MBONA ) {
      /* case: NBONH+1 -> NBONH + MBONA; bonds without hydrogen,
       * except not final bond */

      j = i - parm->NBONH;
      shift = sizeof(ParmBond) * (parm->MBONA - j);
      shiftArray( parm->pbond+j, parm->pbond+j-1, shift );

      shift = sizeof(Bond) * (parm->NBONH + parm->MBONA - i);
      shiftArray( parm->bond+i, parm->bond+i-1, shift );

      if ( i <= parm->NBONA + parm->NBONH) parm->NBONA -= 1;
      parm->MBONA -= 1;

    } else if ( i == (parm->NBONH + parm->MBONA) ) {
      /* case: NBONH + MBONA + 1; final bond */

      if (parm->MBONA == parm->NBONA) parm->NBONA -= 1;
      parm->MBONA -= 1;

    } else {
      warning("deletebond", "No bond of number %i exists\n", i);
      return;
    }
    break;


  case 'A':
  case 'a':

    fprintf(stdout, "\nAttempting to delete angle %i\n\n", i);

    if ( i > 0 && i <= parm->NTHETH ) {
      /* case: 1 -> NTHETH; angle with hydrogen */

      if ( i != parm->NTHETH ) {
	shift = sizeof(ParmAngle) * (parm->NTHETH - i);
	shiftArray( parm->pangleH+i, parm->pangleH+i-1, shift );
      }

      shift = sizeof(Angle) * (parm->NTHETH + parm->MTHETS - i);
      shiftArray( parm->angle+i, parm->angle+i-1, shift );
	    
      parm->NTHETH -= 1;

    } else if ( i > parm->NTHETH && i < parm->NTHETH + parm->MTHETS ) {
      /* case: NTHETH+1 -> NTHETH + MTHETS; angle without hydrogen,
       * except not final bond */

      j = i - parm->NTHETH;
      shift = sizeof(ParmAngle) * (parm->MTHETS - j);
      shiftArray( parm->pangle+j, parm->pangle+j-1, shift );

      shift = sizeof(Angle) * (parm->NTHETH + parm->MTHETS - i);
      shiftArray( parm->angle+i, parm->angle+i-1, shift );
	    
      if ( i <= parm->NTHETA + parm->NTHETH) parm->NTHETA -= 1;
      parm->MTHETS -= 1;

    } else if ( i == (parm->NTHETH + parm->MTHETS) ) {
      /* case: NTHETH + MTHETS + 1; final angle */

      if ( parm->NTHETA == parm->MTHETS) parm->NTHETA -= 1;
      parm->MTHETS -= 1;

    } else {
      warning("deleteangle", "No angle of number %i exists\n", i);
      return;
    }
    break;

  case 'D':
  case 'd':

    fprintf(stdout, "\nAttempting to delete dihedral %i\n\n", i);

    if ( i > 0 && i <= parm->NPHIH ) {
      /* case: 1 -> NPHIH; dihedral with hydrogen */

      if ( i != parm->NPHIH ) {
	shift = sizeof(ParmDihedral) * (parm->NPHIH - i);
	shiftArray( parm->pdihedralH+i, parm->pdihedralH+i-1, shift );
      }

      shift = sizeof(Dihedral) * (parm->NPHIH + parm->MPHIA - i);
      shiftArray( parm->dihedral+i, parm->dihedral+i-1, shift );
	    
      parm->NPHIH -= 1;

    } else if ( i > parm->NPHIH && i < parm->NPHIH + parm->MPHIA ) {
      /* case: NPHIH+1 -> NPHIH + MPHIA; dihedral without hydrogen,
       * except not final bond */

      j = i - parm->NPHIH;
      shift = sizeof(ParmDihedral) * (parm->MPHIA - j);
      shiftArray( parm->pdihedral+j, parm->pdihedral+j-1, shift );

      shift = sizeof(Dihedral) * (parm->NPHIH + parm->MPHIA - i);
      shiftArray( parm->dihedral+i, parm->dihedral+i-1, shift );
	    
      if ( i <= parm->NPHIA + parm->NPHIH) parm->NPHIA -= 1;
      parm->MPHIA -= 1;

    } else if ( i == (parm->NPHIH + parm->MPHIA) ) {
      /* case: NPHIH + MPHIA + 1; final dihedral */

      if ( parm->NPHIA == parm->MPHIA ) parm->NPHIA -= 1;
      parm->MPHIA -= 1;

    } else {
      warning("deletedihedral", "No dihedral of number %i exists\n", i);
      return;
    }

    break;

  }
  isModifiedParm = TRUE;
  isWrittenParm = FALSE;
}



   void
deletePerturbedBAD(char *type, int i)
{
  int shift;

  if ( ! parm->IFPERT ) {
    warning("deletePerturbedBAD()", 
	    "This is not a perturbation file!\n");
    return;
  }

  switch ( type[0] ) {

  case 'B':
  case 'b':

    fprintf(stdout, "\nAttempting to delete perturbing bond %i\n\n", i);

    if ( i > 0 && i <= parm->NBPER ) {  

      if ( i != parm->NGPER ) {
	shift = sizeof(double) * (parm->NBPER - i);
	shiftArray(parm->pert->ibper+i,
		   parm->pert->ibper+i-1, shift);
	shiftArray(parm->pert->jbper+i,
		   parm->pert->jbper+i-1, shift);
	shift = sizeof(double) * (parm->NBPER*2 - i);
	shiftArray(parm->pert->icbper+i,
		   parm->pert->icbper+i-1, shift);
	shift = sizeof(double) * (parm->NBPER*2 - i + parm->NBPER);
	shiftArray(parm->pert->icbper+i+parm->NBPER,
		   parm->pert->icbper+i+parm->NBPER-1, shift);
      }

      parm->NBPER -= 1;

    } else if ( i == parm->NBPER ) {

      parm->NBPER -= 1;
      
    } else {
      warning("deletePerturbedBAD()", 
	      "No perturbed bond of number %i exists\n", i);
      return;
    }
    break;


  case 'A':
  case 'a':

    fprintf(stdout, "\nAttempting to delete perturbed angle %i\n\n", i);

    if ( i > 0 && i <= parm->NGPER ) {  

      if ( i != parm->NGPER ) {
	shift = sizeof(double) * (parm->NGPER - i);
	shiftArray(parm->pert->itper+i,
		   parm->pert->itper+i-1, shift);
	shiftArray(parm->pert->jtper+i,
		   parm->pert->jtper+i-1, shift);
	shiftArray(parm->pert->ktper+i,
		   parm->pert->ktper+i-1, shift);
	shift = sizeof(double) * (parm->NGPER*2 - i);
	shiftArray(parm->pert->ictper+i,
		   parm->pert->ictper+i-1, shift);
	shift = sizeof(double) * (parm->NGPER*2 - i + parm->NGPER);
	shiftArray(parm->pert->ictper+i+parm->NGPER,
		   parm->pert->ictper+i+parm->NGPER-1, shift);
      }

      parm->NGPER -= 1;

    } else if ( i == parm->NGPER ) {

      parm->NGPER -= 1;
      
    } else {
      warning("deletePerturbedBAD()", 
	      "No perturbed angle of number %i exists\n", i);
      return;
    }
    break;

  case 'D':
  case 'd':

    fprintf(stdout, "\nAttempting to delete perturbed dihedral %i\n\n", i);

    if ( i > 0 && i <= parm->NDPER ) {  

      if ( i != parm->NDPER ) {
	shift = sizeof(double) * (parm->NDPER - i);
	shiftArray(parm->pert->ipper+i,
		   parm->pert->ipper+i-1, shift);
	shiftArray(parm->pert->jpper+i,
		   parm->pert->jpper+i-1, shift);
	shiftArray(parm->pert->kpper+i,
		   parm->pert->kpper+i-1, shift);
	shiftArray(parm->pert->lpper+i,
		   parm->pert->lpper+i-1, shift);
	shift = sizeof(double) * (parm->NDPER*2 - i);
	shiftArray(parm->pert->icpper+i,
		   parm->pert->icpper+i-1, shift);
	shift = sizeof(double) * (parm->NDPER*2 - i + parm->NDPER);
	shiftArray(parm->pert->icpper+i+parm->NDPER,
		   parm->pert->icpper+i+parm->NDPER-1, shift);
      }

      parm->NDPER -= 1;

    } else if ( i == parm->NDPER ) {

      parm->NDPER -= 1;
      
    } else {
      warning("deletePerturbedBAD()", 
	      "No perturbed angle of number %i exists\n", i);
      return;
    }
    break;

  }

  if ( parm->NDPER == 0 && parm->NGPER == 0 && parm->NBPER == 0 ) {
    /* this implies there are NO perturbed angles, dihedrals or bonds */
    warning("deletePerturbedBAD()", 
	    "There are no perturbed bonds, angles or dihedrals!\n");
  }

  isModifiedParm = TRUE;
  isWrittenParm = FALSE;
}



   int
processAtomInput(FILE *fpin, FILE *fpout, char *info) 
{
  char line[LINE_SIZE], *linep;
  int atomNumber;

  for ( ; ; ) {
    fprintf(fpout, "Please enter an atom (%s): ", info);
    fflush(fpout);
    if ( fgets(line, LINE_SIZE, fpin) != NULL ) {
      linep = line;
      skipWhitespace(linep);
      if (isalpha(linep[0]) ||
	  strchr(linep, ':' ) != NULL ||
	  strchr(linep, '@' ) != NULL ) {

#ifdef DEBUG
	fprintf(stdout, "Searching for atom spec %s\n", line);
#endif

	/* NOTE: this part is currently not implemented, i.e.
	 * searching for midas style atom spec, hence the following:
         */
	atomNumber = -1;

      } else {
	if (sscanf(linep, "%i", &atomNumber) != 1) 
	  error("processAtomInput()", 
		"error in scanning atom number from %s\n",
		linep);
      }
      if ( (atomNumber < 1) || (atomNumber > parm->NTOTAT) ) 
	fprintf(fpout, "Atom number %i is out of range\n", atomNumber);
      else
	return atomNumber;
    }
  }
}


   double
processDoubleInput(FILE *fpin, FILE *fpout, char *info) 
{
  char line[LINE_SIZE], *linep;
  double value;

  for ( ; ; ) {
    fprintf(fpout, "Please enter %s: ", info);
    fflush(fpout);
    if ( fgets(line, LINE_SIZE, fpin) != NULL ) {
      linep = line;
      skipWhitespace(linep);
      if ( sscanf(linep, "%lf", &value ) != 1) {
	fprintf(fpout, "Error scanning double value from %s\n", linep);
      } else
	return value;
    }
  }
}

   void  
restrainBAD(char *type) 
{
  int at1, at2, at3, at4, tmp1, tmp2;
  size_t size;
  int parameter_index;
  double p1, p2;
  char atom1[20], atom2[20], atom3[20], atom4[20];
  double conversion = 4.0 * (double) atan( (double) 1.0 ) / 180.0;
  double tmp = atan( (double) 1.0 );

  switch ( type[0] ) {

  case 'b':
  case 'B':
    /* A bond constraint */
    at1 = processAtomInput(stdin, stdout, "bond atom 1");
    at2 = processAtomInput(stdin, stdout, "bond atom 2");
    p1 = processDoubleInput(stdin, stdout, "the force constant");
    p2 = processDoubleInput(stdin, stdout, "the equilibrium length");

    parameter_index = parm->MUMBND;
#ifdef OVERWRITE_PARAMETERS
    for (i=0; i < parm->MUMBND; i++) {
      if ( (p1 == parm->rk[i]) && (p2 == parm->req[i]) ) {
	parameter_index = i;
      }
    }
#endif    
    /* re-allocate rdparm bond storage */
    size = sizeof( Bond ) * (parm->NBONH + parm->MBONA);
    parm->bond = safe_realloc(parm->bond, size, sizeof( Bond ));
    /* re-allocate rdparm PARM type bond storage */
    size = sizeof( ParmBond ) * parm->MBONA;
    parm->pbond = safe_realloc(parm->pbond, size, sizeof( ParmBond ));
    /* re-allocate space for PARM parameters if necessary */
    if ( parameter_index == parm->MUMBND ) {
      size = sizeof( double ) * parm->MUMBND;
      parm->rk  = safe_realloc( parm->rk,  size, sizeof(double) );
      parm->req = safe_realloc( parm->req, size, sizeof(double) );
      parm->MUMBND += 1;
    }
    /* load up rdparm bond */
    parm->bond[parm->NBONH + parm->MBONA].atom[0] = at1;
    parm->bond[parm->NBONH + parm->MBONA].atom[1] = at2;
    parm->bond[parm->NBONH + parm->MBONA].rk      = p1;
    parm->bond[parm->NBONH + parm->MBONA].req     = p2;
    parm->bond[parm->NBONH + parm->MBONA].scale   = 1.0;

    /* load up PARM type bond storage */
    parm->rk [parameter_index] = p1;
    parm->req[parameter_index] = p2;
    parm->pbond[parm->MBONA].ib = 3 * at1 - 3;
    parm->pbond[parm->MBONA].jb = 3 * at2 - 3;
    parm->pbond[parm->MBONA].icb = parameter_index+1;

    /* update control variables */
    parm->MBONA += 1;
    
    fprintf(stdout, "Added BOND restraint with force constant %g and length %g\n",
	    p1, p2);
    fprintf(stdout, "Atoms: %s --- %s\n",
	    rdparmPrintAtom(at1, atom1), 
	    rdparmPrintAtom(at2, atom2));
    
    break;
  case 'a':
  case 'A':
    /* An angle constraint */
    at1 = processAtomInput(stdin, stdout, "angle atom 1");
    at2 = processAtomInput(stdin, stdout, "angle atom 2");
    at3 = processAtomInput(stdin, stdout, "angle atom 3");
    p1 = processDoubleInput(stdin, stdout, "the force constant");
    p2 = processDoubleInput(stdin, stdout, "the equilibrium angle");

    /* check to see if parameters represented are already in parms */
    parameter_index = parm->MUMANG;
#ifdef OVERWRITE_PARAMETERS
    for (i=0; i < parm->MUMANG; i++) {
      if ( (p1 == parm->tk[i]) && (p2 * conversion == parm->teq[i]) ) {
	parameter_index = i;
      }
    }
#endif

    /* re-allocate rdparm angle storage */
    size = sizeof( Angle ) * (parm->NTHETH + parm->MTHETS);
    parm->angle = safe_realloc(parm->angle, size, sizeof( Angle ));

    /* re-allocate rdparm PARM type angle storage */
    size = sizeof( ParmAngle ) * parm->MTHETS;
    parm->pangle = 
      safe_realloc(parm->pangle, size, sizeof( ParmAngle ));

    /* re-allocate space for PARM parameters, if necessary */
    if ( parameter_index == parm->MUMANG ) {
      size = sizeof( double ) * parm->MUMANG;
      parm->tk  = safe_realloc( parm->tk,  size, sizeof(double) );
      parm->teq = safe_realloc( parm->teq, size, sizeof(double) );
      parm->MUMANG += 1;
    }

    /* load up rdparm angle */
    parm->angle[parm->NTHETH + parm->MTHETS].atom[0] = at1;
    parm->angle[parm->NTHETH + parm->MTHETS].atom[1] = at2;
    parm->angle[parm->NTHETH + parm->MTHETS].atom[2] = at3;
    parm->angle[parm->NTHETH + parm->MTHETS].tk      = p1;
    parm->angle[parm->NTHETH + parm->MTHETS].teq     = p2 * conversion;
    parm->angle[parm->NTHETH + parm->MTHETS].scale   = 1.0;

    /* load up PARM type angle storage */
    parm->tk [parameter_index] = p1;
    parm->teq[parameter_index] = p2 * conversion;
    parm->pangle[parm->MTHETS].it  = 3 * at1 - 3;
    parm->pangle[parm->MTHETS].jt  = 3 * at2 - 3;
    parm->pangle[parm->MTHETS].kt  = 3 * at3 - 3;
    parm->pangle[parm->MTHETS].ict = parameter_index+1;

    /* update control variables */
    parm->MTHETS += 1;
    
    fprintf(stdout, "Added ANGLE restraint with force constant %g and angle %g\n",
	    p1, p2);
    fprintf(stdout, "Atoms: %s --- %s --- %s\n",
	    rdparmPrintAtom(at1, atom1), 
	    rdparmPrintAtom(at2, atom2),
	    rdparmPrintAtom(at3, atom3));

    break;
  case 'd':
  case 'D':
  case 't':
  case 'T':
    /* A dihedral constraint */
    at1 = processAtomInput(stdin, stdout, "dihedral atom 1");
    at2 = processAtomInput(stdin, stdout, "dihedral atom 2");
    at3 = processAtomInput(stdin, stdout, "dihedral atom 3");
    at4 = processAtomInput(stdin, stdout, "dihedral atom 4");
    p1 = processDoubleInput(stdin, stdout, "the force constant");
    p2 = processDoubleInput(stdin, stdout, "the equilibrium angle");
    if (p2 < 0.0) p2 += 360.0;

    /* check to see if parameters represented are already in parms */
    parameter_index = parm->MPTRA;
#ifdef OVERWRITE_PARAMETERS
    for (i=0; i < parm->MPTRA; i++) {
      if ((p1 == parm->pk[i]) && 
	  ((p2 + 180) * conversion == parm->phase[i]) &&
	  (parm->pn[i] == 1.0 ) ) {
	parameter_index = i;
      }
    }
#endif

    /* re-allocate rdparm dihedral storage */
    size = sizeof( Dihedral ) * (parm->NPHIH + parm->MPHIA);
    parm->dihedral = 
      safe_realloc(parm->dihedral, size, sizeof( Dihedral ));

    /* re-allocate rdparm PARM type dihedral storage */
    size = sizeof( ParmDihedral ) * parm->MPHIA;
    parm->pdihedral = 
      safe_realloc(parm->pdihedral, size, sizeof( ParmDihedral ));

    /* re-allocate space for PARM parameters, if necessary */
    if ( parameter_index == parm->MPTRA ) {
      size = sizeof( double ) * parm->MPTRA;
      parm->pk  =   safe_realloc( parm->pk,    size, sizeof(double) );
      parm->pn  =   safe_realloc( parm->pn,    size, sizeof(double) );
      parm->phase = safe_realloc( parm->phase, size, sizeof(double) );
      parm->MPTRA += 1;
    }
    /* reverse dihedral if necessary */
    if ( at3 == 1 ) {
      tmp1 = at3;
      tmp2 = at4;
      at4 = at1;
      if ( tmp2 < 0 ) at4 = -at4;
      at3 = at2;
      if ( tmp1 < 0 ) at3 = -at3;
      at1 = iabs(tmp2);
      at2 = iabs(tmp1);
    }
	
    /* load up rdparm dihedral */
    parm->dihedral[parm->NPHIH + parm->MPHIA].atom[0] = at1;
    parm->dihedral[parm->NPHIH + parm->MPHIA].atom[1] = at2;
    parm->dihedral[parm->NPHIH + parm->MPHIA].atom[2] = at3;
    parm->dihedral[parm->NPHIH + parm->MPHIA].atom[3] = at4;
    parm->dihedral[parm->NPHIH + parm->MPHIA].pk      = p1;
    parm->dihedral[parm->NPHIH + parm->MPHIA].pn      = 1.0;
    parm->dihedral[parm->NPHIH + parm->MPHIA].phase   = 
      (p2 + 180.0) * conversion;
    parm->dihedral[parm->NPHIH + parm->MPHIA].scale   = 1.0;

    /* load up PARM type dihedral storage */
    parm->pk[parameter_index]    = p1;
    parm->pn[parameter_index]    = 1.0;
    parm->phase[parameter_index] = (p2 + 180.0) * conversion;
    parm->pdihedral[parm->MPHIA].ip  = 3 * at1 - 3;
    parm->pdihedral[parm->MPHIA].jp  = 3 * at2 - 3;
    parm->pdihedral[parm->MPHIA].kp  = - (3 * at3 - 3);
    parm->pdihedral[parm->MPHIA].lp  = 3 * at4 - 3;
    parm->pdihedral[parm->MPHIA].icp = parameter_index+1;

    /* update control variables */
    parm->MPHIA += 1;
    
    fprintf(stdout, "Added DIHEDRAL restraint with force constant %g and dihedral %g\n",
	    p1, p2);
    fprintf(stdout, "Atoms: %s --- %s --- %s --- %s\n",
	   rdparmPrintAtom(at1, atom1), 
	   rdparmPrintAtom(at2, atom2),
	   rdparmPrintAtom(at3, atom3),
	   rdparmPrintAtom(at4, atom4));

    break;
  default:
    warning("restrainBAD", "Unknown thing (%s) to constrain, expecting %s\n",
	    type, "bond, angle, or dihedral");
  }

}

   void
translateRestart(char *filenamep)
{
  double xc, yc, zc;
  int i;

  FILE *fpin, *fpout;
  char *coordin = NULL;
  char *coordout = NULL;

  if ( parm->IFBOX ) {
    fprintf(stdout, "   This routine is not for periodic boxes!\n");
    return;
  }
  
  /* open AMBER coordinate/restart file */
  coordin = promptToOpenFile(&fpin, filenamep, "r",
	    "Input the name of coordinate file to fix: ");
  
  /* read in the data */
  fprintf(stdout, "\n   Reading in the coordinate/restart file...\n");
  safe_free( (void *)  parm->coords );
  parm->coords = readAmberRestart(parm->NTOTAT, coordin, fpin);

  /* accumulate the center of geometry */
  xc = 0.0;
  yc = 0.0;
  zc = 0.0;
  for (i=0; i < parm->coords->natoms; i++) {
    xc += parm->coords->x[i];
    yc += parm->coords->y[i];
    zc += parm->coords->z[i];
  }
  xc /= parm->coords->natoms;
  yc /= parm->coords->natoms;
  zc /= parm->coords->natoms;

  fprintf(stdout, "\n   Transforming coordinates...\n");
  for (i=0; i < parm->coords->natoms; i++) {
    parm->coords->x[i] -= xc;
    parm->coords->y[i] -= yc;
    parm->coords->z[i] -= zc;
  }    
  
  coordout = promptToOpenFile(&fpout, "", "w",
            "Output transformed coordinates to file: ");
  
  writeAmberRestart(parm->coords, coordout, fpout);

  safe_fclose( fpin );
  safe_fclose( fpout );
  
  safe_free((void *) parm->coords->title);
  safe_free((void *) parm->coords->x);
  safe_free((void *) parm->coords->y);
  safe_free((void *) parm->coords->z);
  safe_free((void *) parm->coords->vx);
  safe_free((void *) parm->coords->vy);
  safe_free((void *) parm->coords->vz);
  safe_free((void *) parm->coords);

  fprintf(stdout, "   done...\n\n");

  safe_free( (void *) coordin );
  safe_free( (void *) coordout );
}



   void
translateBox(char *filenamep)
{
  int i;
  double boxx, boxy, boxz;
  double xc, yc, zc;
  int isAmber;
  int isSpasms;
  int outsideBounds = 0; 
  FILE *fpin, *fpout;
  char *coordin = NULL;
  char *coordout = NULL;

  if ( ! parm->IFBOX ) {
    fprintf(stdout, "   This parm file does not contain any box information!\n");
    return;
  }
  
  isAmber = 0;
  isSpasms = 0;
  
  /* open AMBER coordinate/restart file */
  coordin = promptToOpenFile(&fpin, filenamep, "r",
	    "Input the name of coordinate box to fix: ");
  
  /* read in the data */
  fprintf(stdout, "\n   Reading in the coordinate/restart file...\n");
  safe_free( (void *)  parm->coords );
  parm->coords = readAmberRestart(parm->NTOTAT, coordin, fpin);

  /* check out the box */
  fprintf(stdout, "   Checking box...\n");
  if ( parm->coords->box ) {
    /* box coordinates are in restart file, make sure parm has them
     * and check to see if they match...
     */
    if ((parm->box->box[0] != parm->coords->box[0]) ||
	(parm->box->box[1] != parm->coords->box[1]) ||
	(parm->box->box[2] != parm->coords->box[2])) {
      fprintf(stdout, "   NOTE: Box coordinates from restart do not match topology!\n");
      fprintf(stdout, "   The parm box coordinates are:\n%10.3f  %10.3f  %10.3f\n",
	     parm->box->box[0], parm->box->box[1], parm->box->box[2]);
    }
    fprintf(stdout,"   Box coordinates were specified in the restart file.\n");
    fprintf(stdout,"   Most likely, this implies an AMBER constant pressure run.\n");
    fprintf(stdout,"   These restart file coordinates will be used in contrast\n");
    fprintf(stdout,"   to the parm topology coordinates...\n");
    boxx = parm->coords->box[0];
    boxy = parm->coords->box[1];
    boxz = parm->coords->box[2];
  } else {
    /* box coordinates were not in the restart file
     */
    fprintf(stdout,"\n   NOTE: Box coordinates were not specified in the restart\n");
    fprintf(stdout,"   file.  The values from the parm topology will be used...\n");
    boxx = parm->box->box[0];
    boxy = parm->box->box[1];
    boxz = parm->box->box[2];
  }
  fprintf(stdout,"\n   The box coordinates are:\n%10.3f  %10.3f  %10.3f\n",
	  boxx, boxy, boxz);

  /* accumulate the center of geometry */
  xc = 0.0;
  yc = 0.0;
  zc = 0.0;
  for (i=0; i < parm->coords->natoms; i++) {
    xc += parm->coords->x[i];
    yc += parm->coords->y[i];
    zc += parm->coords->z[i];
  }
  xc /= parm->coords->natoms;
  yc /= parm->coords->natoms;
  zc /= parm->coords->natoms;

  /* check to see if center of geometry suggests whether this is an
   * AMBER or SPASMS style box */
  if ((xc < boxx/4.0 && xc > -boxx/4.0) &&
      (yc < boxy/4.0 && yc > -boxy/4.0) &&
      (zc < boxz/4.0 && zc > -boxz/4.0))
    isSpasms = 1;

  if ((xc < 3 * boxx/4.0 && xc > boxx/4.0) &&
      (yc < 3 * boxy/4.0 && yc > boxy/4.0) &&
      (zc < 3 * boxz/4.0 && zc > boxz/4.0))
    isAmber = 1;

  /* check for anomolous case where c-of-g falls outside bounds... */
  if ( isSpasms == 0 && isAmber == 0 ) {
    fprintf(stdout,"\n   This is quite a strange coordinate box since the center\n");
    fprintf(stdout,"   of geometry doesn't fall into a range which suggests\n");
    fprintf(stdout,"   it is either an AMBER or SPASMS style box...\n\n");
    isAmber = promptUserResponse(stdin, stdout,
				 "   Which style of box is it?  [AMBER] ",
				 "amber", 1);
    isSpasms = ( isAmber ? 0 : 1 );
  }
  fprintf(stdout,"\n   Based on the location of the center of geometry of this\n");
  fprintf(stdout,"   box, it appears to be a %s style box...\n",
	 (isSpasms ? "SPASMS" : "AMBER"));
  /* check to see how many atoms have coordinates outside of the box */
  for (i=0; i < parm->coords->natoms; i++) {
    if (parm->coords->x[i] < (isSpasms ? -boxx/2.0 : 0.0) ||
	parm->coords->y[i] < (isSpasms ? -boxy/2.0 : 0.0) ||
	parm->coords->z[i] < (isSpasms ? -boxz/2.0 : 0.0))
      outsideBounds += 1;

    if (parm->coords->x[i] > (isSpasms ? boxx/2.0 : boxx) ||
	parm->coords->y[i] > (isSpasms ? boxy/2.0 : boxy) ||
	parm->coords->z[i] > (isSpasms ? boxz/2.0 : boxz))
      outsideBounds += 1;
  }
  if ( outsideBounds ) {
    fprintf(stdout,"   Approximately %i atoms fall outside the boundaries\n",
	    outsideBounds);
    fprintf(stdout,"   of the %s style box.\n", ( isSpasms ? "SPASMS" : "AMBER"));
  }

  if (! promptUserResponse(stdin, stdout,
          (isSpasms ?
	  "\n   Do you want to convert your box from SPASMS to AMBER? [yes] " :
	  "\n   Do you want to convert your box from AMBER to SPASMS? [yes] "),
          "yes", 1))
    return;
  
  fprintf(stdout,"\n   Transforming coordinates...\n");
  /* NOTE: xc, yc, zc are being reused for scaling... */
  xc = (isSpasms ? boxx/2.0 : -boxx/2.0);
  yc = (isSpasms ? boxy/2.0 : -boxy/2.0);
  zc = (isSpasms ? boxz/2.0 : -boxz/2.0);

  for (i=0; i < parm->coords->natoms; i++) {
    parm->coords->x[i] += xc;
    parm->coords->y[i] += yc;
    parm->coords->z[i] += zc;
  }    
  
  coordout = promptToOpenFile(&fpout, "", "w",
            "Output transformed coordinates to file: ");
  
  writeAmberRestart(parm->coords, coordout, fpout);

  safe_fclose( fpin );
  safe_fclose( fpout );
  
  safe_free((void *) parm->coords->title);
  safe_free((void *) parm->coords->x);
  safe_free((void *) parm->coords->y);
  safe_free((void *) parm->coords->z);
  safe_free((void *) parm->coords->vx);
  safe_free((void *) parm->coords->vy);
  safe_free((void *) parm->coords->vz);
  safe_free((void *) parm->coords);

  fprintf(stdout, "   done...\n\n");

  safe_free( (void *) coordin );
  safe_free( (void *) coordout );

}





   void
initializeParm(Parm *p)
{

  if (p == NULL) return;

  p->filename = NULL;
  p->fp = NULL;
  p->atom = NULL;
  p->nno = NULL;
  p->residue = NULL;
  p->rk = NULL;
  p->req = NULL;
  p->tk = NULL;
  p->teq = NULL;
  p->pk = NULL;
  p->pn = NULL;
  p->phase = NULL;
  p->solty = NULL;
  p->cn1 = NULL;
  p->cn2 = NULL;
  p->pbondH = NULL;
  p->pbond = NULL;
  p->pangleH = NULL;
  p->pangle = NULL;
  p->pdihedral = NULL;
  p->pdihedralH = NULL;
  p->natex = NULL;
  p->ag = NULL;
  p->bg = NULL;
  p->hbcut = NULL;
  p->box = NULL;
  p->cap = NULL;
  p->pert = NULL;
  p->bond = NULL;
  p->angle = NULL;
  p->dihedral = NULL;
  p->coords = NULL;
}


   void
getParm( char *filenamep )
{

  clearParm( parm );
  parm = safe_malloc( sizeof( Parm ) ); 
  initializeParm(parm);

  parm->filename = promptToOpenFile(&parm->fp, filenamep, "r", 
                   "Input the name of an AMBER parm file: ");
/*
  fprintf(stdout, "Opened parm file %s\n", parm->filename);
*/
  readParm();
  
}

   void
putParm(char *filenamep )
{
  FILE *fp;
  char *filename = NULL;

  filename = promptToOpenFile(&fp, filenamep, "w",
             "Input the name of an AMBER parm file to write: ");	      
  fprintf(stdout, "Attemping to write (modified?) parm file (%s) to file: %s\n",
	 parm->filename, filename);
  writeParm( fp, 1 );
  isWrittenParm = TRUE;

  safe_free( (void *) filename );
}

   void
quit()
{
  fprintf(stdout, "EXITING program...\n");
  exit(1);
}


/*

NOTES:

how to print 16.8 doubles:   printf("%16.8E\n", <double_value>);

*/



   int
isMatchResidue(Name residue, char *rs) 
{
  /* SPECIAL CASES: These come first, even though it is less efficient
   * to catch some strange AMBER special cases.  NOTE THAT THIS ROUTINE
   * has not been updated for AMBER4.1 residue naming conventions as of yet...
   */

  /* HIS in AMBER can have many different names... */
  if ( strcmp(rs, "HIS") == 0 ) {
    if (strncmp(residue, "HI", 2) == 0) {
      return 1;
    }
  }

  /* STANDARD search */
  if (strncmp(residue, rs, strlen(rs)) == 0) {
    if (strlen(rs) > 1)
      return 1;
    else {
/*
      warning("isMatchResidue()", 
	"Returning a match (%s to %s) even though this is questionable...\n",
	      residue, rs);
*/
      return 1;
    }
  }
  return 0;
}

/* This routine will take a list of atom Names (not necessary padded
 * correctly) and coerce it AMBER to names where simple naming
 * convention differences are known.  It is intended for mardigras to
 * AMBER names...
 *
 * NOTE: this matching is not yet complete!  FURTHERMORE, it has
 * not been updated to parm94.dat naming conventions...
 */
   void
fixAtomNames(Name *atom1, Name *atom2, int *res1, int *res2,
	     double *upper, double *lower, int total, int offset) 
{
  Name *atoms;
  int *res;
  int i, j;
  int isHA3;
  int isHB3;
  int isHD1, isHD2, isHD3, isHD4, isHD5, isHD6;
  int isHE1, isHE2;
  int isHG1, isHG2;
  int isHZ1, isHZ2;
  int isHH,  isHG;
  int isMK;

  isHA3 = 0;
  isHB3 = 0;
  isHD1 = 0; isHD2 = 0; isHD3 = 0; isHD4 = 0; isHD5 = 0; isHD6 = 0;
  isHE1 = 0; isHE2 = 0; 
  isHG1 = 0; isHG2 = 0; 
  isHZ1 = 0; isHZ2 = 0;
  isHH = 0;  isHG = 0;
  isMK = 0;

  for (i=0; i < total; i++) {
    if ( strcmp(atom1[i], "HA3") == 0 ||
         strcmp(atom2[i], "HA3") == 0 )
      isHA3 = 1;
    if ( strcmp(atom1[i], "HB3") == 0 ||
         strcmp(atom2[i], "HB3") == 0 )
      isHB3 = 1;
    if ( strcmp(atom1[i], "HD1") == 0 ||
	 strcmp(atom2[i], "HD1") == 0 )
      isHD1 = 1;		
    if ( strcmp(atom1[i], "HD2") == 0 ||
	 strcmp(atom2[i], "HD2") == 0 )
      isHD2 = 1;		
    if ( strcmp(atom1[i], "HD3") == 0 ||
	 strcmp(atom2[i], "HD3") == 0 )
      isHD3 = 1;		
    if ( strcmp(atom1[i], "HD4") == 0 ||
	 strcmp(atom2[i], "HD4") == 0 )
      isHD4 = 1;		
    if ( strcmp(atom1[i], "HD5") == 0 ||
	 strcmp(atom2[i], "HD5") == 0 )
      isHD5 = 1;		
    if ( strcmp(atom1[i], "HD6") == 0 ||
	 strcmp(atom2[i], "HD6") == 0 )
      isHD6 = 1;		
    if ( strcmp(atom1[i], "HE1") == 0 ||
         strcmp(atom2[i], "HE1") == 0 )
      isHE1 = 1;		
    if ( strcmp(atom1[i], "HE2") == 0 ||
         strcmp(atom2[i], "HE2") == 0 )
      isHE2 = 1;		
    if ( strcmp(atom1[i], "HG1") == 0 ||
	 strcmp(atom2[i], "HG1") == 0 )
      isHG1 = 1;		
    if ( strcmp(atom1[i], "HG2") == 0 ||
         strcmp(atom2[i], "HG2") == 0 )
      isHG2 = 1;		
    if ( strcmp(atom1[i], "HZ1") == 0 ||
	 strcmp(atom2[i], "HZ1") == 0 )
      isHZ1 = 1;		
    if ( strcmp(atom1[i], "HZ2") == 0 ||
	 strcmp(atom2[i], "HZ2") == 0 )
      isHZ2 = 1;		
    if ( strcmp(atom1[i], "HH") == 0 ||
         strcmp(atom2[i], "HH") == 0 )
      isHH = 1;			
    if ( strcmp(atom1[i], "HG") == 0 ||
	 strcmp(atom2[i], "HG") == 0 )
      isHG = 1;
    if ( strcmp(atom1[i], "MK") == 0 ||
	 strcmp(atom2[i], "MK") == 0 )
      isMK = 1;
  }


  for (j=0; j < 2 ; j++) {
    if (j == 0) {
      atoms = atom1;
      res = res1;
    }
    if (j == 1) {
      atoms = atom2;
      res = res2;
    }

  for (i=0; i < total; i++) {

    /* ACE:
     *   convert MK to CH3 and add/subtract 1.0 angstroms
     */
    if (strcmp(parm->residue[res[i]+offset-1].labres, "ACE ") == 0) {
      if ( isMK ) {
	if ( strcmp(atoms[i], "MK") == 0) {
	  fprintf(stdout, 
		  "Detected a MK atom in a ACE residue.  Assuming it is\n");
	  fprintf(stdout, "a psuedo atom and hence we want to change it to\n");
	  fprintf(stdout, "a CH3 and modify the upper/lower bounds...\n");
	  if ( promptUserResponse(stdin, stdout,
				  "Does you want to do this?  [yes] ", 
				  "yes", 1)) {
	    warning("fixAtomNames", 
		    "Modified MK to CH3 in constraint %i (%s)",
		    i+1, atoms[i]);
	    strcpy(atoms[i], "CH3");
	    upper[i] += 1.0;
	    lower[i] -= 1.0;
	    if (lower[i] < 0.0) lower[i] = 0.0;
	  }
	}
      }
    }

    /* ASP */
    if (strcmp(parm->residue[res[i]+offset-1].labres, "ASP ") == 0) {
      if ( isHB3 == 0 ) {
	if ( strcmp(atoms[i], "HB2") == 0) {
	  warning("fixAtomNames", 
		  "ASP: Modified HB2 to HB3 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HB3");
	} else if ( strcmp(atoms[i], "HB1") == 0) {
	  warning("fixAtomNames", 
		  "ASP: Modified HB1 to HB2 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HB2");
	}
      }
    }

    /* GLN */
    if (strcmp(parm->residue[res[i]+offset-1].labres, "GLN ") == 0) {
      if ( isHB3 == 0 ) {
	if ( strcmp(atoms[i], "HB2") == 0) {
	  warning("fixAtomNames", 
		  "GLN: Modified HB2 to HB3 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HB3");
	} else if ( strcmp(atoms[i], "HB1") == 0) {
	  warning("fixAtomNames", 
		  "GLN: Modified HB1 to HB2 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HB2");
	}
      }
/*
      if ( strcmp(atoms[i], "HE21") == 0) {
	warning("fixAtomNames", 
		"GLN: Modified HE21 to HNE1 in constraint %i (%s)",
		i+1, atoms[i]);
	strcpy(atoms[i], "HNE1");
      } else if ( strcmp(atoms[i], "HE22") == 0) {
	warning("fixAtomNames", 
		"GLN: Modified HE22 to HNE2 in constraint %i (%s)",
		i+1, atoms[i]);
	strcpy(atoms[i], "HNE2");
      }
*/
      if ( isHG1 || isHG2 ) {
	if ( strcmp(atoms[i], "HG2") == 0) {
	  warning("fixAtomNames", 
		  "GLN: Modified HG2 to HG3 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HG3");
	} else if ( strcmp(atoms[i], "HG1") == 0) {
	  warning("fixAtomNames", 
		  "GLN: Modified HG1 to HG2 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HG2");
	}
      }
    } 
  
    /* GLY */
    if (strcmp(parm->residue[res[i]+offset-1].labres, "GLY ") == 0) {
      if ( isHA3 == 0 ) {
	if ( strcmp(atoms[i], "HA2") == 0) {
	  warning("fixAtomNames", 
		  "GLY: Modified HA2 to HA3 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HA3");
	} else if ( strcmp(atoms[i], "HA1") == 0) {
	  warning("fixAtomNames", 
		  "GLY: Modified HA1 to HA2 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HA2");
	}
      }
    }

    /* HIS, HIP, HID, HIE */
    if (strcmp(parm->residue[res[i]+offset-1].labres, "HIS ") == 0 ||
	strcmp(parm->residue[res[i]+offset-1].labres, "HIP ") == 0 ||
	strcmp(parm->residue[res[i]+offset-1].labres, "HIE ") == 0 ||
	strcmp(parm->residue[res[i]+offset-1].labres, "HID ") == 0) {
      if ( isHB3 == 0 ) {
	if ( strcmp(atoms[i], "HB2") == 0) {
	  warning("fixAtomNames", 
		  "HIS: Modified HB2 to HB3 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HB3");
	} else if ( strcmp(atoms[i], "HB1") == 0) {
	  warning("fixAtomNames", 
		  "HIS: Modified HB1 to HB2 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HB2");
	}
      }
      if ( isHD1 || isHD2 ) {
	if ( strcmp(atoms[i], "HD1") == 0) {
	  warning("fixAtomNames", 
		  "HIS: Modified HD1 to HND in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HND");
	} else if ( strcmp(atoms[i], "HD2") == 0) {
	  warning("fixAtomNames", 
		  "HIS: Modified HD2 to HD in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HD");
	}
      }
      if ( isHE1 ) {
	if ( strcmp(atoms[i], "HE1") == 0) {
	  warning("fixAtomNames", 
		  "HIS: Modified HE1 to HE in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HE");
	}
      }
    }
     
    /* LEU */
    if (strcmp(parm->residue[res[i]+offset-1].labres, "LEU ") == 0) {
      if ( isHB3 == 0 ) {
	if ( strcmp(atoms[i], "HB2") == 0) {
	  warning("fixAtomNames", 
		  "LEU: Modified HB2 to HB3 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HB3");
	} else if ( strcmp(atoms[i], "HB1") == 0) {
	  warning("fixAtomNames", 
		  "LEU: Modified HB1 to HB2 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HB2");
	}
      }
      if ( isHD1 || isHD2 || isHD3 || isHD4 || isHD5 || isHD6 ) {
	if ( strcmp(atoms[i], "HD1") == 0) {
	  warning("fixAtomNames", 
		  "LEU: Modified HD1 to HD11 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HD11");
	} else if ( strcmp(atoms[i], "HD2") == 0) {
	  warning("fixAtomNames", 
		  "LEU: Modified HD2 to HD12 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HD12");
	} else if ( strcmp(atoms[i], "HD3") == 0) {
	  warning("fixAtomNames", 
		  "LEU: Modified HD3 to HD13 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HD13");
	} else if ( strcmp(atoms[i], "HD4") == 0) {
	  warning("fixAtomNames", 
		  "LEU: Modified HD4 to HD21 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HD21");
	} else if ( strcmp(atoms[i], "HD5") == 0) {
	  warning("LEU: fixAtomNames", 
		  "Modified HD5 to HD22 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HD22");
	}
	else if ( strcmp(atoms[i], "HD6") == 0) {
	  warning("fixAtomNames", 
		  "LEU: Modified HD6 to HD23 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HD23");
	}
      }
    }

    /* LYS */
    if (strcmp(parm->residue[res[i]+offset-1].labres, "LYS ") == 0) {
      if ( isHB3 == 0 ) {
	if ( strcmp(atoms[i], "HB2") == 0) {
	  warning("fixAtomNames", 
		  "LYS: Modified HB2 to HB3 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HB3");
	} else if ( strcmp(atoms[i], "HB1") == 0) {
	  warning("fixAtomNames", 
		  "LYS: Modified HB1 to HB2 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HB2");
	}
      }
      if ( isHD3 == 0 ) {
	if ( strcmp(atoms[i], "HD2") == 0) {
	  warning("fixAtomNames", 
		  "LYS: Modified HD2 to HD3 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HD3");
	} else if ( strcmp(atoms[i], "HD1") == 0) {
	  warning("fixAtomNames", 
		  "LYS: Modified HD1 to HD2 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HD2");
	}
      }
      if ( isHE1 || isHE2 ) {
	if ( strcmp(atoms[i], "HE2") == 0) {
	  warning("fixAtomNames", 
		  "LYS: Modified HE2 to HE3 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HE3");
	} else if ( strcmp(atoms[i], "HE1") == 0) {
	  warning("fixAtomNames", 
		  "LYS: Modified HE1 to HE2 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HE2");
	}
      }
      if ( isHG1 || isHG2 ) {
	if ( strcmp(atoms[i], "HG2") == 0) {
	  warning("fixAtomNames", 
		  "LYS: Modified HG2 to HG3 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HG3");
	} else if ( strcmp(atoms[i], "HG1") == 0) {
	  warning("fixAtomNames", 
		  "LYS: Modified HG1 to HG2 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HG2");
	}
      }
      if ( isHZ1 || isHZ2 ) {
	if ( strcmp(atoms[i], "HZ1") == 0) {
	  warning("fixAtomNames", 
		  "LYS: Modified HZ1 to HNZ1 in constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HNZ1");
	} else if ( strcmp(atoms[i], "HZ2") == 0) {
	  warning("fixAtomNames", 
		  "LYS: Modified HZ2 to HNZ2 constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HNZ2");
	}
      }
    } 

    /* MET */
    if (strcmp(parm->residue[res[i]+offset-1].labres, "MET ") == 0) {
      if ( isHB3 == 0 ) {
	if ( strcmp(atoms[i], "HB2") == 0) {
	  warning("fixAtomNames", 
		  "MET: Modified HB2 to HB3 constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HB3");
	} else if ( strcmp(atoms[i], "HB1") == 0) {
	  warning("fixAtomNames", 
		  "MET: Modified HB1 to HB2 constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HB2");
	}
      }
      if ( isHG1 || isHG2 ) {
	if ( strcmp(atoms[i], "HG2") == 0) {
	  warning("fixAtomNames", 
		  "MET: Modified HG2 to HG3 constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HG3");
	} else if ( strcmp(atoms[i], "HG1") == 0) {
	  warning("fixAtomNames", 
		  "MET: Modified HG1 to HG2 constraint %i (%s)",
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HG2");
	}
      }
    }

    /* SER */
    if (strcmp(parm->residue[res[i]+offset-1].labres, "SER ") == 0) {
      if ( isHB3 == 0 ) {
	if ( strcmp(atoms[i], "HB2") == 0) {
	  warning("fixAtomNames", 
		  "SER: Modified HB2 to HB3 constraint %i (%s)", 
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HB3");
	} else if ( strcmp(atoms[i], "HB1") == 0) {
	  warning("fixAtomNames", 
		  "SER: Modified HB1 to HB2 constraint %i (%s)", 
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HB2");
	}
      }
      if ( isHG == 0 ) {
	if ( strcmp(atoms[i], "HG") == 0) {
	  warning("fixAtomNames", 
		  "SER: Modified HG to HOG constraint %i (%s)", i+1, atoms[i]);
	  strcpy(atoms[i], "HOG");
	}
      }
    }

    /* TYR */
    if (strcmp(parm->residue[res[i]+offset-1].labres, "TYR ") == 0) {
      if ( isHB3 == 0 ) {
	if ( strcmp(atoms[i], "HB2") == 0) {
	  warning("fixAtomNames", 
		  "TYR: Modified HB2 to HB3 constraint %i (%s)", 
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HB3");
	} else if ( strcmp(atoms[i], "HB1") == 0) {
	  warning("fixAtomNames", 
		  "TYR: Modified HB1 to HB2 constraint %i (%s)", 
		  i+1, atoms[i]);
	  strcpy(atoms[i], "HB2");
	}
      }
      if ( isHH == 0 ) {
	if ( strcmp(atoms[i], "HH") == 0) {
	  warning("fixAtomNames", 
		  "TYR: Modified HH to HOH constraint %i (%s)", i+1, atoms[i]);
	  strcpy(atoms[i], "HOH");
	}
      }
    }
  }
  }
}


   void
doMardi2Sander(char *filenameIn)
{
  FILE *fpin, *fpout;
  char buffer[BUFFER_SIZE], junk[50], junk1[50];
  int allocated, space, constraints;
  double *lower, *upper, *average;
  Name *resName1, *resName2, *atomName1, *atomName2;
  int *res1, *res2;
  int offset = 0;
  int ret;
  int i, j, k, badMatch;
  float rk2, rk3;
  char *filenameInReturn = NULL;
  char *filenameOutReturn = NULL;

  fprintf(stdout, "\nThis routine will attempt to read in a mardigras\n");
  fprintf(stdout, "style distance constraint file, matching the atom names\n");
  fprintf(stdout, "to the atom names of the currently loaded parameter\n");
  fprintf(stdout, "file and then will output sander style constraints to\n");
  fprintf(stdout, "file which can later be appended to your default sander\n");
  fprintf(stdout, "input file...\n\n");

  filenameInReturn = promptToOpenFile(&fpin, filenameIn, "r",
		"Input the name of a mardigras type constraint file: ");

  space = 1000;
  lower   = (double *) safe_malloc(sizeof(double) * space);
  upper   = (double *) safe_malloc(sizeof(double) * space);
  average = (double *) safe_malloc(sizeof(double) * space);
  atomName1 = (Name *) safe_malloc(sizeof(Name) * space);
  atomName2 = (Name *) safe_malloc(sizeof(Name) * space);
  resName1  = (Name *) safe_malloc(sizeof(Name) * space);
  resName2  = (Name *) safe_malloc(sizeof(Name) * space);
  res1 = (int *) safe_malloc(sizeof(int) * space);
  res2 = (int *) safe_malloc(sizeof(int) * space);
  allocated = space;

  /* READ all lines of text in the file searching for constraint lines */

  fprintf(stdout, "Reading in constraints...\n");
  constraints = 0; 

  while ( fgets(buffer, BUFFER_SIZE, fpin) != NULL ) {
    if ( strlen(buffer) > 1 && strstr(buffer, "CONSTRAINT") != NULL) {
      /* allocate more space if necessary... */
      if (constraints == allocated) {
	lower   = (double *) 
	  safe_realloc(lower,   
		       sizeof(double) * allocated, sizeof(double) * space);
	upper   = (double *)
	  safe_realloc(upper,   
		       sizeof(double) * allocated, sizeof(double) * space);
	average = (double *)
	  safe_realloc(average, 
		       sizeof(double) * allocated, sizeof(double) * space);
	atomName1 = (Name *) 
	  safe_realloc(atomName1, 
		       sizeof(Name) * allocated, sizeof(Name) * space);
	atomName2 = (Name *) 
	  safe_realloc(atomName2, 
		       sizeof(Name) * allocated, sizeof(Name) * space);
	resName1 = (Name *) 
	  safe_realloc(resName1, 
		       sizeof(Name) * allocated, sizeof(Name) * space);
	resName2 = (Name *) 
	  safe_realloc(resName2, 
		       sizeof(Name) * allocated, sizeof(Name) * space);
	res1 = (int *) 
	  safe_realloc(res1, sizeof(int) * allocated, sizeof(int) * space);
	res2 = (int *) 
	  safe_realloc(res2, sizeof(int) * allocated, sizeof(int) * space);
	allocated += space;
      }
      /* scan if the information */
      ret = sscanf(buffer, "%s%d%s%s%d%s%s%lf%lf%lf",
		   junk, &res1[constraints], resName1[constraints], 
		   atomName1[constraints], &res2[constraints], 
		   resName2[constraints], atomName2[constraints],
		   &upper[constraints], &lower[constraints], 
		   &average[constraints]);
      if (ret < 9) {
	warning("mardi2sander()", 
           "Error on scanning constraint, skipping it, line follows\n%s\n",
	   buffer);
      } else {
	/* scan successful */
	if ( ret == 9 ) average[constraints] = 0.0;
	constraints++;

      }
    }
  }
  safe_fclose(fpin);
  fprintf(stdout, "Successfully read in %i constraints\n", constraints);
  

  fprintf(stdout, 
     "Attempting to match residue names...  If the first residue in\n");
  fprintf(stdout, 
     "the parm file doesn't match, the second is tried and so on, in\n");
  fprintf(stdout, 
     "case the numbering is shifted (by the addition of extra residues\n");
  fprintf(stdout, 
     "at the beginning of the sequence...\nThis is currently *not*");
  fprintf(stdout, " robust!\n\n");
  badMatch = 0;
  for (i=0; i < constraints; i++) {
    if (isMatchResidue(parm->residue[res1[i]-1].labres, resName1[i]) == 0 ||
	isMatchResidue(parm->residue[res2[i]-1].labres, resName2[i]) == 0) {

      fprintf(stdout, 
	"\nConstraint %i: failed to match mardigras residues %i:%s or %i:%s\n",
	      i, res1[i], resName1[i], res2[i], resName2[i]);
      fprintf(stdout, "to parm residues %i:%s or %i:%s...\n",
	      res1[i], parm->residue[res1[i]-1].labres, 
	      res2[i], parm->residue[res2[i]-1].labres);
      badMatch = 1;
      break;
    }      
  }
  offset = 0;
  if (badMatch) {
    fprintf(stdout, 
       "\nHmmmn, let's try to see if end residues were added by shifting\n");
    fprintf(stdout, "the numbering...\n");
    i = 0; 
    offset = 0;
    while (res1[i]-1+offset < parm->NTOTRS && 
	   res2[i]-1+offset < parm->NTOTRS ) {
      j = isMatchResidue(parm->residue[res1[i]-1+offset].labres, resName1[i]);
      k = isMatchResidue(parm->residue[res2[i]-1+offset].labres, resName2[i]);
      if (j && k )
	break;
      offset++;
    }
    if (res1[i]-1+offset == parm->NTOTRS ||
	res2[i]-1+offset == parm->NTOTRS) {
      fprintf(stdout, 
	      "Nope, shifting the offset didn't work.  This suggests\n");
      fprintf(stdout, 
	      "that the naming between the mardigras constraint file or\n");
      fprintf(stdout, 
	      "residue numbering is significantly off...\n");
      return;
    }

    badMatch = 0;
    for (i=0; i < constraints; i++) {

      if (isMatchResidue(parm->residue[res1[i]-1+offset].labres, 
			 resName1[i]) == 0 ||
	  isMatchResidue(parm->residue[res2[i]-1+offset].labres, 
			 resName2[i]) == 0) {

	fprintf(stdout, 
           "\nConstraint %i: failed to match mardigras %i:%s or %i:%s\n",
		i, res1[i], resName1[i], res2[i], resName2[i]);
	fprintf(stdout, "to parm residues %i:%s or %i:%s...\n",
		res1[i]+offset, parm->residue[res1[i]-1+offset].labres, 
		res2[i]+offset, parm->residue[res2[i]-1+offset].labres);
	badMatch = 1;
	break;
      }      
    }
  }
  if (badMatch) {
    fprintf(stdout, 
       "\nDarn!  The program wasn't able to fix the problem by simply\n");
    fprintf(stdout, 
       "shifting the numbering of the parm residues.  Perhaps\n");
    fprintf(stdout, 
       "non-standard residue naming is apparent or the residue numbering\n");
    fprintf(stdout, 
       "is significantly different.  Please modify the constraint file...\n");
    return;
  } else {
    fprintf(stdout, 
	    "Yes!  We were able to match the residue names with an offset\n");
    fprintf(stdout, 
	    "from the parm residue to the mardigras numbering of %i...\n", 
	    offset);
  }


  fprintf(stdout, "\nAttemping now to match atom names...\n");
  fixAtomNames(atomName1, atomName2, res1, res2,
	       upper, lower, constraints, offset);
  
  fprintf(stdout, "WARNING: This code is still under development and may\n");
  fprintf(stdout, "not be fully functional.  Following is a list of\n");
  fprintf(stdout, "constraints found with the amber naming conventions\n");
  fprintf(stdout, "Please check out this list...  If you find any\n");
  fprintf(stdout, "inconsistencies or get any warning's, please send them\n");
  fprintf(stdout, "on by!  NOTE: not all residues have been set up yet!\n");
  fprintf(stdout, "(cheatham@cgl.ucsf.edu)\n\n");

  for (i=0; i < constraints; i++) {
    for (j = parm->residue[res1[i]-1+offset].ipres; 
	 j < parm->residue[res1[i]+offset].ipres; j++) {
      if ( isMatchResidue(parm->atom[j-1].igraph, atomName1[i]) == 1 )
      break;
    }
    for (k = parm->residue[res2[i]-1+offset].ipres; 
	 k < parm->residue[res2[i]+offset].ipres; k++) {
      if ( isMatchResidue(parm->atom[k-1].igraph, atomName2[i]) == 1 )
	break;
    }

    if (j == parm->residue[res1[i]+offset].ipres)
      fprintf(stdout, 
	 "Constraint %i: no match, residue %i (%s) for atom1 %s\n",
	      i, res1[i], parm->residue[res1[i]].labres, atomName1[i]);
    else if (k == parm->residue[res2[i]+offset].ipres)
      fprintf(stdout, 
         "Constraint %i: no match, residue %i (%s) for atom2 %s\n",
	      i, res2[i], parm->residue[res2[i]].labres, atomName2[i]);
    else {
      fprintf(stdout, 
	 "Constraint %i, %8s -- %8s %7.3f and %7.3f\n",
	      i, rdparmPrintAtom(j, junk), rdparmPrintAtom(k, junk1), lower[i], upper[i]);
    }
  }

  if ( !promptUserResponse(stdin, stdout, 
			  "Does this list look O.K.?  [yes] ", "yes", 1))
    return;

  fprintf(stdout, "Currently we are assuming that the lower bound\n");
  fprintf(stdout, "for the parabolic (to linear) constraint is 0.0\n");
  fprintf(stdout, "The upper bound (after which the constraint becomes\n");
  fprintf(stdout, "linear is the upper bound + 2.0\n");
  fprintf(stdout, "What do you want as the force constant for the constraint\n");
  fprintf(stdout, "lower bound? (rk2)  ");
  fflush(stdout);
  if (fscanf(stdin, "%f", &rk2) != 1) {
    warning("doMardi2Sander()", "Error on scanning a floating point value\n");
    return;
  }
  fprintf(stdout, "What do you want as the force constant for the constraint\n");
  fprintf(stdout, "upper bound? (rk3)  ");
  fflush(stdout);
  if (fscanf(stdin, "%f", &rk3) != 1) {
    warning("doMardi2Sander()", "Error on scanning a floating point value\n");
    return;
  }

  filenameOutReturn = promptToOpenFile(&fpout, "", "w",
		"Input the name of a file to output constraint info: ");
  
  for (i=0; i < constraints; i++) {
    for (j = parm->residue[res1[i]-1+offset].ipres; 
	 j < parm->residue[res1[i]+offset].ipres; j++) {
      if ( isMatchResidue(parm->atom[j-1].igraph, atomName1[i]) == 1 )
	break;
    }
    for (k = parm->residue[res2[i]-1+offset].ipres; 
	 k < parm->residue[res2[i]+offset].ipres; k++) {
      if ( isMatchResidue(parm->atom[k-1].igraph, atomName2[i]) == 1 )
	break;
    }

    if (j != parm->residue[res1[i]+offset].ipres &&
        k != parm->residue[res2[i]+offset].ipres) {


      if ( lower[i] > upper[i] ) {
	fprintf(stdout, "WARNING, constraint %i: lower > upper\n",
		i+1);
      }

      fprintf(fpout, " &rst iat= %i, %i, 0, 0,\n",
	      res1[i]+offset, res2[i]+offset);
      fprintf(fpout, 
	      "    atnam(1)='%s',atnam(2)='%s',\n",
	      parm->atom[j-1].igraph, parm->atom[k-1].igraph);
      fprintf(fpout, 
	      "    iresid=1, r1=0.0,r2=%7.3f,r3=%7.3f,r4=%7.3f,rk2=%6.2f,rk3=%6.2f &end\n",
	      lower[i], upper[i], upper[i]+2.0,rk2,rk3);
    }
  }
  fprintf(fpout, " &rst iat=0 &end\n");

  /* clean-up */
  safe_free( (void *) lower );
  safe_free( (void *) upper );
  safe_free( (void *) average );
  safe_free( (void *) atomName1 );
  safe_free( (void *) atomName2 );
  safe_free( (void *) resName1 );
  safe_free( (void *) resName2 );

  safe_free( (void *) filenameInReturn );
  safe_free( (void *) filenameOutReturn );
}




/* routine to dump out the current box information and prompt the
 * user for any changes...
 */
   void
modifyBoxInfo()
{
  float boxx, boxy, boxz;
  float beta;
  int iptres;
  char buffer[BUFFER_SIZE];

  if (parm->IFBOX) {
    fprintf(stdout, "\n NOTE: This routine is for displaying and changing the\n");
    fprintf(stdout, " basic box information.  The arrays holding the information\n");
    fprintf(stdout, " on the molecules will not be printed here.  This information\n");
    fprintf(stdout, " is available with the parminfo command or may be analyzed\n");
    fprintf(stdout, " with the modifyMoleculeInformation command.\n\n");
    fprintf(stdout, " BOX SIZE is %g by %g by %g angstroms\n",
           parm->box->box[0], parm->box->box[1], parm->box->box[2]);
    fprintf(stdout, " BETA   is %-8.2f (angle of box)\n", parm->box->beta);
    fprintf(stdout, " IPTRES is %-8i (final residue that is part of solute)\n",
           parm->box->iptres);
    fprintf(stdout, " NSPM   is %-8i (the total number of molecules)\n",
           parm->box->nspm);
    fprintf(stdout, " NSPSOL is %-8i (first solvent molecule)\n",
           parm->box->nspsol);

    if ( ! promptUserResponse(stdin, stdout,
        " Do you want to modify the box information?  [no] ", "no", 1) ) {
      fprintf(stdout, "NOTE: changes will only be saved if you write the parm file!\n");
      fprintf(stdout, "New box in X direction: [%.4f]  ",
             parm->box->box[0]);
      fflush(stdout);
      fgets(buffer, BUFFER_SIZE, stdin);
      if (sscanf(buffer, "%f", &boxx) == 1)
        parm->box->box[0] = boxx;
      fprintf(stdout, "New box in Y direction: [%.4f]  ",
             parm->box->box[1]);
      fflush(stdout);
      fgets(buffer, BUFFER_SIZE, stdin);
      if (sscanf(buffer, "%f", &boxy) == 1)
        parm->box->box[1] = boxy;
      fprintf(stdout, "New box in Z direction: [%.4f]  ",
             parm->box->box[2]);
      fflush(stdout);
      fgets(buffer, BUFFER_SIZE, stdin);
      if (sscanf(buffer, "%f", &boxz) == 1)
        parm->box->box[2] = boxz;

      fprintf(stdout, "Box angle: [%.4f]  ",
             parm->box->beta);
      fflush(stdout);
      fgets(buffer, BUFFER_SIZE, stdin);
      if (sscanf(buffer, "%f", &beta) == 1)
        parm->box->beta = beta;

      fprintf(stdout, "Last residue that is part of solute: [%i]  ",
             parm->box->iptres);
      fflush(stdout);
      fgets(buffer, BUFFER_SIZE, stdin);
      if (sscanf(buffer, "%i", &iptres) == 1)
        parm->box->iptres = iptres;

      fprintf(stdout, "\n NEW Box information is...\n\n");
      fprintf(stdout, " BOX SIZE is %g by %g by %g angstroms\n",
             parm->box->box[0], parm->box->box[1], parm->box->box[2]);
      fprintf(stdout, " BETA   is %-8.2f (angle of box)\n", parm->box->beta);
      fprintf(stdout, " IPTRES is %-8i (final residue that is part of solute)\n",
             parm->box->iptres);
      fprintf(stdout, " NSPM   is %-8i (the total number of molecules)\n",
             parm->box->nspm);
      fprintf(stdout, " NSPSOL is %-8i (first solvent molecule)\n",
             parm->box->nspsol);
    }

  } else
    fprintf(stdout, "\n No box information is present in the topology file!\n");
}


   void
modifyMolInfo()
{
  int i, j, imol, count, iat, numwat, addmol, icont;
  int *residues;
  Name prevResidue;

  if (! parm->IFBOX) {
    fprintf(stdout, " There is no box information in this parm file!\n");
    return;
  }

  fprintf(stdout, "\nHere is the current molecule information:\n");
  fprintf(stdout, " NUMBER of atoms in each molecule, molecule %i to %i\n",
         1, parm->box->nspm);
  count=1;
  for (i= 0; i < parm->box->nspm; i++) {
    fprintf(stdout, "%5i ", parm->box->nsp[i]);
    if ( count % 10 == 0 ) fprintf(stdout, "\n");
    count++;
  }
  if ( ! count % 10 == 0 ) fprintf(stdout, "\n");

  fprintf(stdout, "Analyzing the molecules...\n");

  /* create a list of every residue each atom is in...
   */
  residues = safe_malloc(sizeof(int) * (parm->NTOTAT+1));
  for (i=0; i < parm->NTOTRS; i++) {
    for (j=parm->residue[i].ipres-1; j < parm->residue[i+1].ipres-1; j++) {
      residues[j] = i;
    }
  }
  residues[parm->NTOTAT] = parm->NTOTRS;

  /* loop over molecules, with "iat" as running index and print out what
   * residues are in each molecule...
   */
  fprintf(stdout, "Printing out the residues in each molecule.  NOTE: if\n");
  fprintf(stdout, "the molecule only has one residue and repeats, \"...\" is\n");
  fprintf(stdout, "printed to save space...\n\n");
  prevResidue[0] = (char) 0;
  addmol = 0;
  for (iat = 0, imol = 0;
       iat < parm->NTOTAT;
       iat += parm->box->nsp[imol], imol++) {

    /*
     *  add a check on imol to prevent blowing the array in case of
     *  a corrupted parm file...
     */
    if (imol >= parm->box->nspm) {
      warning("modifyMolInfo()", 
	      "Blew the molecule array (nsp); molecule %i >= total molecules (%i)\n",
	      imol, parm->box->nspm);
      return;
    }
    
    if (addmol == 0) numwat = 0;

    if (prevResidue[0] != (char) 0 &&
        residues[iat+parm->box->nsp[imol]]-residues[iat] == 1 &&
        strcmp(prevResidue, parm->residue[ residues[iat] ].labres) == 0) {
      if (icont == 0) {
        fprintf(stdout, "                   ...\n");
        icont = 1;
      }
    } else
      icont = 0;


    if (imol == parm->box->nspm-1 || icont == 0)
      fprintf(stdout, "Molecule %5i: ", imol+1);
    for (j = residues[iat], count = 1;
         j < residues[iat+parm->box->nsp[imol]]; j++) {
      if (imol == parm->box->nspm-1 || icont == 0) {
        fprintf(stdout, "%s ", parm->residue[j].labres);
        if (count % 10 == 0) fprintf(stdout, "\n                ");
        count++;
      }

      if ( strcmp(parm->residue[j].labres, "WAT ") == 0 && addmol == 0 )
        numwat += 1;

    }
    if (imol == parm->box->nspm-1 || (icont == 0 && (count-1) % 10 != 0))
      fprintf(stdout, "\n");

    if (residues[iat+parm->box->nsp[imol]]-residues[iat] == 1)
      strcpy(prevResidue, parm->residue[ residues[iat] ].labres);
    else {
      prevResidue[0] = (char) 0;
      icont = 0;
    }
    if (numwat > 1 &&
        numwat == residues[iat+parm->box->nsp[imol]]-residues[iat] &&
        addmol == 0)
      addmol = imol+1;
  }
  fprintf(stdout, " NSPM   is %-8i (the total number of molecules)\n",
	  parm->box->nspm);
  fprintf(stdout, " IPTRES is %-8i (final residue that is part of solute)\n",
	  parm->box->iptres);
  fprintf(stdout, " NSPSOL is %-8i (first solvent molecule)\n",
	  parm->box->nspsol);

  if (addmol > 0) {
    fprintf(stdout, 
	    "\nDetected %i water residues in a single molecule, molecule %i\n",
	    numwat, addmol);
    fprintf(stdout, 
	    "This is an error and will lead to incorrect periodic imaging\n");
    fprintf(stdout, "of these waters, improper pressure scaling and incorrect\n");
    fprintf(stdout, "densities.  Perhaps you added waters using the ADD water\n");
    fprintf(stdout, 
	    "option to edit which generates this incorrect behavior...\n\n");

    if ( promptUserResponse(stdin, stdout,
      "Do you want to try and fix this?  [no] ", "no", 1) ) {
      return;
    }

    parm->box->nsp = (int *)
      safe_realloc(parm->box->nsp, sizeof(int) * parm->box->nspm,
                   sizeof(int) * numwat-1);
    for (i=parm->box->nspm; i < parm->box->nspm + numwat-1; i++) {
      parm->box->nsp[i] = 3;
    }
    parm->box->nsp[addmol-1] = 3;

    for (imol=0, iat = 0; imol < addmol-1; imol++)
      iat += parm->box->nsp[imol];
    parm->box->iptres = residues[iat-1]+1;
    parm->box->nspm += numwat-1;
    if (parm->box->nspsol > addmol)
      parm->box->nspsol = addmol;


    fprintf(stdout, "Done!  Current information follows:\n");
    fprintf(stdout, " NUMBER of atoms in each molecule, molecule %i to %i\n",
           1, parm->box->nspm);
    count=1;
    for (i= 0; i < parm->box->nspm; i++) {
      fprintf(stdout, "%5i ", parm->box->nsp[i]);
      if ( count % 10 == 0 ) fprintf(stdout, "\n");
      count++;
    }
    if ( ! count % 10 == 0 ) fprintf(stdout, "\n");
    prevResidue[0] = (char) 0;
    fprintf(stdout, "Printing out the residues in each molecule.  NOTE: if\n");
    fprintf(stdout, "the molecule only has one residue and repeats, \"...\"\n");
    fprintf(stdout, "is printed to save space...\n\n");
    for (iat = 0, imol = 0;
         iat < parm->NTOTAT;
         iat += parm->box->nsp[imol], imol++) {

      if (prevResidue[0] != (char) 0 &&
          residues[iat+parm->box->nsp[imol]]-residues[iat] == 1 &&
          strcmp(prevResidue, parm->residue[ residues[iat] ].labres) == 0) {
        if (icont == 0) {
          fprintf(stdout, "                   ...\n");
          icont = 1;
        }
      } else
        icont = 0;

      if (imol == parm->box->nspm-1 || icont == 0) {
        fprintf(stdout, "Molecule %5i: ", imol+1);
        for (j = residues[iat], count = 1;
             j < residues[iat+parm->box->nsp[imol]]; j++) {
          fprintf(stdout, "%s ", parm->residue[j].labres);
          if (count % 10 == 0) fprintf(stdout, "\n                ");
          count++;
        }
        if ( ! (count % 10 == 0) ) fprintf(stdout, "\n");
      }

      if (residues[iat+parm->box->nsp[imol]]-residues[iat] == 1)
        strcpy(prevResidue, parm->residue[ residues[iat] ].labres);
      else {
        prevResidue[0] = (char) 0;
        icont = 0;
      }
    }
    fprintf(stdout, "\n NSPM   is %-8i (the total number of molecules)\n",
           parm->box->nspm);
    fprintf(stdout, " IPTRES is %-8i (final residue that is part of solute)\n",
           parm->box->iptres);
    fprintf(stdout, " NSPSOL is %-8i (first solvent molecule)\n", 
	    parm->box->nspsol);
    fprintf(stdout, "\n\nRemember that to save these changes, you must write\n");
    fprintf(stdout, "a new parm file (see writeparm).\n");
  }
}
