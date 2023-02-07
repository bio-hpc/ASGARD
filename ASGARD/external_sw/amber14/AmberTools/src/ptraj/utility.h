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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/utility.h,v 10.3 2010/01/28 17:03:53 droe Exp $
 *
 *  Revision: $Revision: 10.3 $
 *  Date: $Date: 2010/01/28 17:03:53 $
 *  Last checked in by $Author: droe $
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



typedef struct _arrayType {
  int length;
  void *entry;
} arrayType;

typedef struct _stackType {
  void *entry;
  struct _stackType *next;
} stackType;

#  ifdef __STDC__
typedef void (*fxnPrintStackEntry)(void *);
#  else
typedef void (*fxnPrintStackEntry)();
#  endif

#ifndef UTILITY_MODULE

extern void warning();
extern void error();

#  ifdef __STDC__

extern void * safe_malloc(int);
extern void * SafeMalloc(char *, int, size_t);
extern void * safe_realloc(void *, int, int);
extern void   safe_free(void *);

extern int scanString( FILE *,   char *, int, char *);
extern int scanInt(    FILE *,    int *, int, char *);
extern int scanDouble( FILE *, double *, int, char *);

extern void shiftArray( void *, void *, int );
extern int  stringContains(char *, char *);
extern int  stringMatch(char *, char *);
extern char ** stringSplit( char *, char *);
extern char * toLowerCase( char * );

extern void pushBottomStack( stackType **, void * );
extern void pushStack( stackType **, void * );
extern void *popStack( stackType ** );
extern void clearStack( stackType ** );
extern void printStack( stackType **, fxnPrintStackEntry, char * );
extern char *copyString( char *);
extern void printfone( char *, ... );

#  else

extern void * safe_malloc();
extern void * SafeMalloc();
extern void * safe_realloc();
extern void   safe_free();
extern int scanString();
extern int scanInt();
extern int scanDouble();

extern void shiftArray();
extern int  stringContains();
extern int  stringMatch();
extern char ** stringSplit();
extern char * toLowerCase();

extern void pushBottomStack();
extern void pushStack();
extern void *popStack();
extern void clearStack();
extern void printStack();
extern char *copyString();
extern void printfone();

#  endif
#endif

#define iabs(a) ( a > 0 ? a : -a )

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef skipWhitespace   
#define skipWhitespace( xxx ) \
  while( ((xxx[0] == '\n') || isspace(xxx[0])) && strlen(xxx) ) xxx++;
#endif













