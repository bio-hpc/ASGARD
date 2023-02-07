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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/interface.c,v 10.0 2008/04/15 23:24:11 case Exp $
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
 *  interface.c
 *
 *  This is the basic user interface for rdparm and ptraj; the major 
 *  functionality is implemented in the file dispatch.h
 *
 *  Currently only two user interface modes are supported, that for 
 *  "rdparm" and that for "ptraj".  The "ptraj" input file processing is
 *  performed by ptraj() in ptraj.c and that for "rdparm" is performed here.
 *  The "rdparm" interface is slightly more interactive in that the user is
 *  prompted, text is processed and a command immediately performed.
 *  The "ptraj" interface works more in the mode of processing an input
 *  file until a particular command is reached (or EOF) and then a whole
 *  series of commands are performed.  In both cases, conversion of text
 *  typed by the user to commands run by this program are performed by 
 *  the dispatchToken() routine defined in dispatch.c.
 */

#include <stdio.h>
#include <string.h>

#include "ptraj.h"


char *rdparm_header = "\n   RDPARM MENU.  Please enter commands.  Use \"?\" or \"help\"\n   for more info.  \"exit\" or \"quit\" to leave program...\n\n";
char *rdparm_prompt = "\nRDPARM MENU:  ";


   void
interface(interfaceMode mode, char *filename)
{
  stackType *argumentStack;
  char buffer[BUFFER_SIZE];
  Token *tokenlist;

  argumentStack = NULL;

  switch(mode) {

  case INTERFACE_PTRAJ:

    /*
     *  ptraj input file processing is performed in ptraj(), ptraj.c
     */

    tokenlist = (Token *) &ptrajTokenlist;
    ptraj(filename);
    break;

  case INTERFACE_RDPARM:
    

    tokenlist = (Token *) &rdparmTokenlist;

    fprintf(stdout, rdparm_header);
    fprintf(stdout, rdparm_prompt);
    while (1) {

      fflush(stdout);
      if ( fgets(buffer, BUFFER_SIZE, stdin) == NULL ) {

	warning("interface()", "NULL input\n");

      } else {

	dispatchToken( tokenlist, argumentStack, (char *) buffer);

      }
    
      fprintf(stdout, rdparm_prompt);
    }
  }
}


