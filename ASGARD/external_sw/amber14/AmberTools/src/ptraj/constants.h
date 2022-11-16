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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/constants.h,v 10.0 2008/04/15 23:24:11 case Exp $
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




#ifndef LINE_SIZE
#define LINE_SIZE 512
#endif

#ifndef BUFFER_SIZE
#define BUFFER_SIZE 2048
#endif

#ifndef FILENAME_SIZE
#define FILENAME_SIZE 512
#endif


/*
 * Conversion factor for a point charge to kcals/mol.
 * This is the value defined in AMBER Parm and LEaP.
 */
#ifndef CHARGE_TO_KCALS
#define CHARGE_TO_KCALS 18.2223
#endif

#ifndef DEFAULT_GRID_SIZE 
#define DEFAULT_GRID_SIZE 100
#endif

#ifndef DEFAULT_GRID_SPACING 
#define DEFAULT_GRID_SPACING 0.5
#endif


/*
 *  USEFUL MACROS
 */

#define ABS(a) ((a)>=0 ? (a) : -(a))



#  ifdef MAIN_MODULE
double PI, RADDEG, DEGRAD, SMALL;
#  else
extern double PI, RADDEG, DEGRAD, SMALL;
#  endif

