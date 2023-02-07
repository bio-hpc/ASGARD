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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/main.c,v 10.1 2009/08/24 23:58:48 rcw Exp $
 *
 *  Revision: $Revision: 10.1 $
 *  Date: $Date: 2009/08/24 23:58:48 $
 *  Last checked in by $Author: rcw $
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


/*  _______________________________________________________________________
 *
 *
 *              THIS IS THE MAIN PROGRAM OF PTRAJ AND RDPARM
 *
 *         This program is used for reading and processing AMBER
 *         (and potentially other) topology and trajectory files.  
 *
 *  The program is command driven via a simple textual interface converting
 *  text typed by the user into a commands and arguments run to process
 *  trajectory and topology/parameter files.  For more information, read the
 *  documents provided with AMBER.  Information is also maintained in a HTML
 *  file "rdparm.html" in this directory.  Developers may also want to read
 *  the detailed comments provided in the source code and in particular in
 *  files ptraj.c, actions.c and dispatch.c and the associated header files.
 *
 *  In order to create an executable, various ancillary code is required and 
 *  not all of it is described in detail here.
 *
 *   o  interface.c    --- code defining the basic user interface.  Note that
 *                         much of the underlying functionality is implemented
 *                         in dispatch.c
 *
 *   o  utility.c      --- defines error routines, safe memory allocation, etc.
 *
 *   o  rdparm.c       --- defines much of the functionality for reading and 
 *                         processing AMBER (and other) topologies/parameters.
 *                         User access to many of these routines is accessed
 *                         indirectly through the ptraj interface with the 
 *                         command rdparm or directly if the executable is 
 *                         named rdparm.
 *
 *   o  ptraj.c        --- the main code for trajectory processing.
 *
 *   o  actions.c      --- the code implementing the ptraj "actions"
 *
 *   o  dispatch.c     --- defines the token lookup, and string stack processing,
 *                         and handles the conversion from strings typed by the 
 *                         user to subroutines run.
 *
 *   o  io.c           --- input/output routines, file handing.  Note:
 *                         coordinate/trajectory IO is in trajectory.c
 *
 *   o  trajectory.c   --- special input/output routines for 
 *                         trajectories/coordinates
 *
 *   o  rms.c          --- code to calculate root-mean-squared deviations between
 *                         conformations.
 *
 *   o  display.c      --- code for printing 2D RMS plots
 *
 *   o  torsion.c      --- code for calculating torsions
 *
 *   o  experimental.c --- new stuff
 *
 *   o  pdb/<files>    --- code to read/write PDB files from the Computer Graphics
 *                         lab at UCSF
 *
 *   o  main.  c       --- this routine
 *
 ***********************************************************************/


#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAIN_MODULE
#include "ptraj.h"

/* PF - multiptraj - include mpi header */
#ifdef MPI
#include <mpi.h>
#endif

   void
calculateConstants()
{
  PI = 4.0 * (double) atan( (double) 1.0 );
  RADDEG = (double) 180.0 / PI;
  DEGRAD = PI / (double) 180.0;
  SMALL = 0.00000000000001;
}


   int
main(int argCount, char **argPointer)
{
  ptrajState **statep;
  char *name;

  if (argCount > 1) {
  /* --version: Print version and exit */
    if (strcmp(argPointer[1],"--version")==0 ||
        strcmp(argPointer[1],"-V"       )==0   ) 
    {
      printf("PTRAJ: Version %s\n",PTRAJ_VERSION_STRING);
      return 0; 
    }
  /* --defines: Print compiler defines and exit */
    if (strcmp(argPointer[1],"--defines")==0) {
      printf("Compiled with:");
#ifdef MPI
      printf(" -DMPI");
#endif
#ifdef BINTRAJ
      printf(" -DBINTRAJ");
#endif
      printf("\n");
      return 0;
    }
  }

#ifdef MPI
  /* MULTIPTRAJ */
  MPI_Init(&argCount, &argPointer);
  MPI_Comm_size(MPI_COMM_WORLD, &worldsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &worldrank);
  /* MULTIPTRAJ */
#else
  worldrank=0;
  worldsize=1;
#endif

  calculateConstants();

  name = strrchr(argPointer[0], (char) '/');
  if (name == NULL) name = argPointer[0];

    /*
     *  If the executable name contains the string rdparm, start up the
     *  rdparm parser
     */

  if (strstr(name, "rdparm") != NULL) {

#ifdef MPI
    /* MULTIPTRAJ */
    if (worldsize > 1) {
      if (worldsize == 0)
	fprintf(stderr, "RDPARM cannot be run with > 1 processor.\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(0);
    }
    /* MULTIPTRAJ */
#endif

    fprintf(stdout, "\n");
    fprintf(stdout, "  \\-/ \n");
    fprintf(stdout, "  -/-   Welcome to RDPARM (an interactive utility for reading AMBER topology files)\n"); 
    fprintf(stdout, "  /-\\                  Type \"help\" or \"?\" for more information\n");
    fprintf(stdout, "  \\-/ \n");
    fprintf(stdout, "  -/-   Version: %s\n", PTRAJ_VERSION_STRING);
    fprintf(stdout, "  /-\\ \n\n");

    statep = ptrajCurrentState();
    if (argCount > 1) {
      ptrajInitializeState( statep, argPointer[1] );
    } else {
      printf("        Usage: %s <parm-file>\n\n", argPointer[0]);
      ptrajInitializeState( statep, NULL );
    }

    interface(INTERFACE_RDPARM, NULL);

  } else {

    /*
     *  otherwise, start up the PTRAJ parser
     */

    statep = ptrajCurrentState();

    if (worldrank==0) {
      fprintf(stdout, "\n");
      fprintf(stdout, "  \\-/  \n");
      fprintf(stdout, "  -/-   PTRAJ: a utility for processing trajectory files\n");
      fprintf(stdout, "  /-\\  \n");
      fprintf(stdout, "  \\-/   Version: %s\n", PTRAJ_VERSION_STRING);
      fprintf(stdout, "  -/-   Executable is: \"%s\"\n", argPointer[0]);
      fprintf(stdout, "  /-\\   Running on %i processor(s)\n",worldsize);
      fprintf(stdout, "  \\-/   ");
      fflush(stdout);
    }

    /* Single processor */
    if (worldsize==1) {
      if (argCount > 1)
        ptrajInitializeState( statep, argPointer[1] ); /* At least ptraj <top> */
      else
        ptrajInitializeState( statep, NULL ); /* Just ptraj */

      if (argCount > 2)
        interface(INTERFACE_PTRAJ, argPointer[2]); /* ptraj <top> <in> */
      else
        interface(INTERFACE_PTRAJ, NULL); /* ptraj <top> */
    }

#ifdef MPI
    else {
      /* Multiprocessor - must be run with top and input file */
      /*
	if(worldsize != 2) {
        if (worldrank==0) fprintf(stderr,"Currently multiptraj only works for 2 nodes.\n");
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(0);
      }
      */
      if (argCount < 3 ) {
        if (worldrank==0) fprintf(stderr,"MULTIPTRAJ must be run with <TOP> and <INPUT>.\n");
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(0);
      }
      ptrajInitializeState( statep, argPointer[1] );
      interface(INTERFACE_PTRAJ, argPointer[2] );
    } /* worldsize */
#endif

  } /* rdparm/ptraj */

#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif

  exit(0);
}
