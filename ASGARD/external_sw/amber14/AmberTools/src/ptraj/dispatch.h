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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/dispatch.h,v 10.0 2008/04/15 23:24:11 case Exp $
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


typedef void (*fxnPntr)();


typedef struct _Token {
   char *str;      /* The trigger string                                  */
   int minMatch;   /* minimum characters required to match the trigger    */
   int numArg;     /* number of arguments following the trigger           */
   int dispatch;   /* a number representing the expected argument types   */
   fxnPntr fxn;    /* a pntr to a function to be executed when triggered  */
} Token;


#  ifdef __STDC__
typedef void (*DispatchFxnType)( int, int, fxnPntr );
#  else
typedef void (*DispatchFxnType)();
#  endif


#  ifdef __STDC__
typedef void (*UserFxnType)( stackType, int, void *);
#  else
typedef void (*UserFxnType)();
#  endif


#ifndef DISPATCH_MODULE 


/*
 * Prototypes for externally visible functions
 */

#  ifdef __STDC__

extern void dispatchToken( Token *, stackType *, char *);
extern void dispatch(int, int, fxnPntr);
extern Token *searchTokenList(Token *, char *);
extern Token rdparmTokenlist[];
extern Token ptrajTokenlist[];

extern char  *getArgumentString(stackType **, char *);
extern char  *getArgumentStringLower(stackType **, char *);
extern int    getArgumentInteger(stackType **, int);
extern float  getArgumentFloat(stackType **, float);
extern double getArgumentDouble(stackType **, double);

extern int    argumentStringContains(stackType **, char *);
extern int    argumentStackContains(stackType **, char *);
extern char  *argumentStackKeyToString(stackType **, char *, char *);
extern int    argumentStackKeyToInteger(stackType **, char *, int);
extern float  argumentStackKeyToFloat(stackType **, char *, float);
extern double argumentStackKeyToDouble(stackType **, char *, double);


#  else

extern void dispatchToken();
extern void dispatch();
extern Token *searchTokenList();
extern Token rdparmTokenlist[];
extern Token ptrajTokenlist[];

extern char  *getArgumentString();
extern char  *getArgumentStringLower();
extern int    getArgumentInteger();
extern float  getArgumentFloat();
extern double getArgumentDouble();

extern int    argumentStringContains();
extern int    argumentStackContains();
extern char  *argumentStackKeyToString();
extern int    argumentStackKeyToInteger();
extern float  argumentStackKeyToFloat();
extern double argumentStackKeyToDouble();

#  endif

#endif /* DISPATCH_MODULE */

