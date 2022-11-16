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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/help.c,v 10.0 2008/04/15 23:24:11 case Exp $
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

#define HELP_MODULE
#include "ptraj.h"

/*
 *  This source file contains the interactive help information for the
 *  rdparm functionality
 */

   void
help(char *help_trigger)
{
  help_trigger = toLowerCase( help_trigger );

  if ( strcmp( help_trigger, "" ) != 0 ) {

    /* HELP for PRINT{bonds,angles,dihedrals} */
    if (strncmp(help_trigger, "printb",  (size_t) 6) == 0 ||
	strncmp(help_trigger, "bonds",   (size_t) 5) == 0 ||
	strncmp(help_trigger, "printan", (size_t) 7) == 0 ||
	strncmp(help_trigger, "angle",   (size_t) 5) == 0 ||
	strncmp(help_trigger, "printd",  (size_t) 6) == 0 ||
	strncmp(help_trigger, "dihed",   (size_t) 5) == 0) {
      printf("\n  bonds, printbonds\n");
      printf("  angles, printangles\n");
      printf("  dihedrals, printdihedrals\n\n");
      printf("  Displays a list of information about the bonds or angles\n");
      printf("  or dihedrals listed in the current AMBER topology file.\n");
      printf("  The format is the number (used for deleting), the parm\n");
      printf("  reference number, the atoms involved (index into parm\n");
      printf("  arrays and Midas style format) and the parameters...\n\n");
      return;

      /* HELP for PRINTPERTURBED{bonds,angles,dihedrals} */
    } else if (strncmp(help_trigger, "printb",  (size_t) 6) == 0 ||
	       strncmp(help_trigger, "perturbedbonds",   (size_t) 10) == 0 ||
	       strncmp(help_trigger, "perturbedangles",  (size_t) 10) == 0 ||
	       strncmp(help_trigger, "perturbeddihedral",(size_t) 10) == 0 ||
	       strncmp(help_trigger, "pertbond",         (size_t)  5) == 0 ||
	       strncmp(help_trigger, "pertangle",        (size_t)  5) == 0 ||
	       strncmp(help_trigger, "pertdihedral",     (size_t)  5) == 0) {
      printf("\n  pertbonds, perturbedbonds\n");
      printf("  pertangles, perturbedangles\n");
      printf("  pertdihedrals, perturbeddihedrals\n\n");
      printf("  Displays a list of information about the bonds or angles\n");
      printf("  or dihedrals perturbed in the current AMBER topology file.\n");
      printf("  The format is the number (used for deleting), the\n");
      printf("  parameters and the atom numbers and names (Midas style\n");
      printf("  format)...\n\n");
      return;

      /* HELP for PRINTATOMS, ATOMS */
    } else if (strncmp(help_trigger, "atoms",  (size_t) 5) == 0 ||
	       strncmp(help_trigger, "printa", (size_t) 6) == 0) {
      printf("\n  atoms, printatoms\n\n");
      printf("  Displays all of the atoms defined in the current AMBER\n");
      printf("  topology file.  This include the atom number (index),\n");
      printf("  name, charge, residue (number and name), atom type, and\n");
      printf("  tree symbol and additionally perturbation information\n");
      printf("  if present...\n\n");
      return;

      /* HELP for PRINTLENNARDJONES, PRINTLJ */
    } else if (strncmp(help_trigger, "printl",  (size_t) 6) == 0) {
      printf("\n  printLennardJones\n\n");
      printf("  Displays all the Lennard-Jones parameters for each\n");
      printf("  possible atom type interaction...\n\n");
      return;

      /* HELP for TRANSFORM */
    } else if (strncmp(help_trigger, "transf", (size_t) 6) == 0) {
      printf("\n  transform <trajectory>\n\n");
      printf("  Transforms an AMBER trajectory file by translating\n");
      printf("  each coordinate set to place its center of geometry\n");
      printf("  or mass at the origin.  The user is prompted for an\n");
      printf("  output filename...\n\n");
      return;

      /* HELP for RMS */
    } else if (strncmp(help_trigger, "rms", (size_t) 3) == 0) {
      printf("\n  rms <trajectory>\n\n");
      printf("  This routine allows one to rms structures within\n");
      printf("  an AMBER trajectory file.  Currently the help file is\n");
      printf("  out of date.  If a trajectory file is not specified\n");
      printf("  you will be prompted, etc...\n\n");
      return;

      /* HELP for DELETE */
    } else if (strncmp(help_trigger,  "delete", (size_t) 6) == 0) {
      printf("\n  delete bond <number>\n");
      printf("  delete angle <number>\n");
      printf("  delete dihedral <number>\n\n");
      printf("  This routine (perhaps slightly dangerously) allows one\n");
      printf("  to remove bonds, angles, and dihedrals from the parm file\n");
      printf("  The number specified (<number>) is the value printed in\n");
      printf("  in the first column when the bonds, angles, or dihedrals\n");
      printf("  are printed.  NOTE that after each deletion, the entries\n");
      printf("  are re-numbered.  The original parm file read in is not\n");
      printf("  changed.  In order to write out the changes, use the\n");
      printf("  writeparm\n\n");
      return;

      /* HELP for DELETE PERTURBED */
    } else if (strncmp(help_trigger,  "delperturbed", (size_t) 6) == 0) {
      printf("\n  delperturbed bond <number>\n");
      printf("  delperturbed angle <number>\n");
      printf("  delperturbed dihedral <number>\n\n");
      printf("  This routine (perhaps slightly dangerously) allows one\n");
      printf("  to remove perturbed bonds, angles, and dihedrals from\n");
      printf("  the parm file.  The number specified (<number>) is the\n");
      printf("  value printed in the first column when the perturbed\n");
      printf("  bonds, angles, or dihedrals are printed (see pertbond,\n");
      printf("  pertangle and pertdihedral).  NOTE the numbering changes\n");
      printf("  after each deletion.  The original parm file read in is\n");
      printf("  not changed.  In order to write out the changes, use the\n");
      printf("  writeparm\n\n");
      return;

      /* HELP for RESTRAIN */
    } else if (strncmp(help_trigger,  "restrain", (size_t) 8) == 0) {
      printf("\n  restrain bond\n");
      printf("  restrain angle\n");
      printf("  restrain dihedral\n\n");
      printf("  This routine allows one to add restraints to the\n");
      printf("  topology file, as normally performed by PARM.  The user\n");
      printf("  is prompted for the relevant information.\n\n");
      return;

      /* HELP for VERBOSE */
    } else if (strncmp(help_trigger, "verbose", (size_t) 7) == 0) {
      printf("\n  verbose\n\n");
      printf("  This command toggles the verbosity mode, by default on\n");
      printf("  Currently verbose is %s\n\n", 
	     ( verboseStatus ? "on" : "off" ));
      return;

      /* HELP FOR SYSTEM */
    } else if (strncmp(help_trigger, "system", (size_t) 6) == 0) {
      printf("\n  system <string>\n\n");
      printf("  This command simply executes the UNIX system command which\n");
      printf("  will execute the <string> command in the shell.\n");
      printf("  NOTE: in order to process a multiple word string as a\n");
      printf("  single argument, the string must be quoted with either\n");
      printf("  \"\"'s or \\\\'s or ``'s\n\n");
      return;

      /* HELP for OPENPARM */
    } else if (strncmp(help_trigger, "openparm", (size_t) 8) == 0) {
      printf("\n  openparm <filename>\n\n");
      printf("  This command reads in a AMBER topology/parameter file.\n");
      printf("  Note that the most currently read/modified parm file is\n");
      printf("  is used by the various routines herein.\n\n");
      return;

      /* HELP for WRITEPARM */
    } else if (strncmp(help_trigger, "writeparm", (size_t) 9) == 0) {
      printf("\n  writeparm <filename>\n\n");
      printf("  This command writes out the most recently read/modified\n");
      printf("  AMBER topology/parameter file to the filename specified.\n\n");
      return;
      
      /* HELP for PARMINFO */
    } else if (strncmp(help_trigger, "parminfo", (size_t) 8) == 0) {
      printf("\n  status\n\n");
      printf("  This command prints information about the current\n");
      printf("  AMBER parameter/topology file\n\n");
      return;

      /* HELP for DDRIVE */
    } else if (strncmp(help_trigger, "ddrive", (size_t) 6) == 0) {
      printf("\n  ddrive <filename>\n\n");
      printf("  This command facilitates the creation of SPASMS dihedral\n");
      printf("  driver input files.  The user is prompted for necessary\n");
      printf("  information, such as the particular dihedral angle to\n");
      printf("  drive (use the numbers reported by print dihedrals)\n\n");
      return;

      /* HELP for TRANSLATEBOX */
    } else if (strncmp(help_trigger, "transl", (size_t) 6) == 0) {
      printf("\n  translatebox <filename>\n\n");
      printf("  This command analyzes a AMBER or SPASMS coordinate set\n");
      printf("  with box coordinates.  It allows conversion of the box\n");
      printf("  from SPASMS to AMBER or AMBER to SPASMS reference frames\n\n");
      return;

      /* HELP for ANALYZE */
    } else if (strncmp(help_trigger, "analyze", (size_t) 7) == 0) {
      printf("\n  analyze <coordinate filename>\n\n");
      printf("  This command analyzes a AMBER or SPASMS coordinate set\n");
      printf("  or trajectory file.  It is currently under development.\n\n");
      return;

      /* HELP for MARDI2SANDER */
    } else if (strncmp(help_trigger, "mardi2sander", (size_t) 6) == 0) {
      printf("\n  mardi2sander <constraint filename>\n\n");
      printf("  This command analyzes a set of mardigras constraints\n");
      printf("  for the current parm file, outputting a file that has the\n");
      printf("  constraints placed in sander format.  The user is prompted for\n");
      printf("  necessary information...  (Currently under development).\n\n");
      return;

      /* HELP for PRINTEXCLUDED */
    } else if (strncmp(help_trigger, "printexcluded", (size_t) 7) == 0) {
      printf("\n  printExcluded\n\n");
      printf("  This command dumps out the excluded atom list\n\n");
      return;

      /* HELP for PRINTLENNARDJONES */
    } else if (strncmp(help_trigger, "printlennard", (size_t) 7) == 0) {
      printf("\n  printLennardJones\n\n");
      printf("  This command dumps out the various non bond interaction\n");
      printf("  and parameters...\n\n");
      return;

      /* HELP for PRINTTYPES */
    } else if (strncmp(help_trigger, "printtypes", (size_t) 6) == 0) {
      printf("\n  printTypes\n\n");
      printf("  This command dumps out the various atom types from the\n");
      printf("  convoluted LJ parameters in the parm file in terms of\n");
      printf("  r* and epsilon...\n\n");
      return;

      /* HELP for MODIFYBOXINFO */
    } else if (strncmp(help_trigger, "modifyboxinfo", (size_t) 7) == 0) {
      printf("\n  modifyboxinfo\n\n");
      printf("  This command dumps out the current box information and\n");
      printf("  prompts the user for changes...\n\n");
      return;

      /* HELP for MODIFYMOLECULEINFO */
    } else if (strncmp(help_trigger, "modifymoleculeinfo", (size_t) 7) == 0) {
      printf("\n  modifyMoleculeInfo\n\n");
      printf("  This command dumps out the current molecule information\n");
      printf("  and looks for inconsistencies...\n\n");
      return;

      /* HELP for CHECKCOORDS */
    } else if (strncmp(help_trigger, "checkcoords", (size_t) 6) == 0) {
      printf("\n  checkcoords mdcrd\n\n");
      printf("  This command reads in the trajectory file and determines\n");
      printf("  if it is corrupted and tells how many sets are active...\n\n");
      return;

      /* HELP for STRIPWATER */
    } else if (strncmp(help_trigger, "stripwater", (size_t) 6) == 0) {
      printf("\n  stripwater\n\n");
      printf("  This command will remove waters from the trajectory file.  It\n");
      printf("  assumes that all the water is located sequentially at the end of\n");
      printf("  the prmtop file and that the water is a three point water.\n\n");
      return;

      /* HELP for PTRAJ (or newtransform) */
    } else if (strncmp(help_trigger, "ptraj", (size_t) 5) == 0 ||
	       strncmp(help_trigger, "newtransform", (size_t) 6) == 0) {
      printf("\n  ptraj\n\n");
      printf("  This command parses a script to perform a variety of analysis and\n");
      printf("  modifications to AMBER trajectory files.  See either the HTML\n");
      printf("  manual ptraj.html or the printed manual for more information\n\n");
      return;
    }


    /* ADD MORE HELP here */
    else {
      printf("\n\n   UNKNOWN help string %s\n", help_trigger);
    }
  }
  
  /* HELP for HELP */
  printf("\n   The following commands are currently available:\n\n");
  printf("   help, ?\n");
  printf("   atoms, printAtoms\n");
  printf("   bonds, printBonds\n");
  printf("   angles, printAngles\n");
  printf("   dihedrals, printDihedrals\n");
  printf("   pertbonds, perturbedBonds\n");
  printf("   pertangles, perturbedAngles\n");
  printf("   pertdihedrals, perturbedDihedrals\n");
  printf("   printExluded\n");
  printf("   printLennardJones\n");
  printf("   printTypes\n");
  printf("   parmInfo\n");
  printf("   checkcoords\n");
  printf("   ddrive       <filename>\n");
  printf("   delete       <bond || angle || dihedral> <number>\n");
  printf("   delperturbed <bond || angle || dihedral> <number>\n");
  printf("   restrain     <bond || angle || dihedral>\n");
  printf("   openparm     <filename>\n");
  printf("   writeparm    <filename>\n");
  printf("   system       <string>\n");
  printf("   mardi2sander <constraint file>\n");
  printf("   analyze      <AMBER trajectory || AMBER coordinates>\n");
  printf("   rms          <AMBER trajectory>\n");
  printf("   stripwater\n");
  printf("   transform    <AMBER trajectory>\n");
  printf("   translateBox <AMBER coords>\n");
  printf("   modifyBoxInfo\n");
  printf("   modifyMolInfo\n");
  printf("   quit, exit\n");
  printf("\n");
}


