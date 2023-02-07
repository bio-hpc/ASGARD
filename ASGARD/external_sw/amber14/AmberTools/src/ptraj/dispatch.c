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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/dispatch.c,v 10.3 2009/12/28 20:24:04 isjoung Exp $
 *
 *  Revision: $Revision: 10.3 $
 *  Date: $Date: 2009/12/28 20:24:04 $
 *  Last checked in by $Author: isjoung $
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

#define DISPATCH_MODULE
#include "ptraj.h"


/*  The dispatch and token lookup routines defined herein provide
 *  a consistent programmer interface between text typed by
 *  a ``user'' and the underlying code/function calls (and their
 *  arguments).  This is handled via the definition of Token
 *  structures that map specific ``trigger'' strings to specific 
 *  functions to be called.  In addition, this structure contains 
 *  information about expected arguments (more on this below).
 *
 *  Currently two Token list structures are defined which specify
 *  the interface for the "rdparm" and "ptraj" routines; both are defined
 *  in the code below.
 *
 *    "rdparm" --- rdparmTokenlist
 *    "ptraj"  --- ptrajTokenlist
 *  
 *  Each of these Token list structures contains:
 *
 *  "str"      -->  A ``trigger string'' (which matches to a specific function).
 *
 *  "minMatch" -->  The minimum number of characters typed by the user that
 *                  are necessary to match the trigger.
 *
 *  "numArg"   -->  The number of arguments (following the trigger string)
 *                  expected for that function.  These arguments are simply
 *                  additional white space delineated strings (that are not
 *                  treated as trigger strings).  Quotes can be used to
 *                  create strings containing whitespace characters.
 *
 *  "dispatch" -->  An enumerated type representing the ``types'' of the 
 *                  arguments (i.e. how to process the strings) or simply
 *                  codes for the function that will be executed.
 *
 *  "fxn"      -->  A function pointer which points to the function which is
 *                  to be called when that trigger is encountered.
 *
 *  To search the Token list, process and convert the arguments
 *  (if necessary) and execute the specified function, the routine
 *  dispatchToken() is called:
 *
 *     dispatchToken(tokenlist, arguments, text)
 *
 *  "tokenlist" -- the Token structure to be searched
 *  "arguments" -- a stack for processing/converting the string arguments
 *  "text"      -- represents the line of text type by the user
 *
 *  Argument processing can be done in two ways:
 *
 *  (1) If the number of arguments is negative, then the remaining set of
 *      whitespace delineated strings on that line of text is assumed to be
 *      arguments to that function and will be placed onto a stack.  This stack
 *      must then be processed by the called (dispatched) function.  It is this
 *      functions responsibility to cleanup the stack (i.e. free any associated
 *      strings) and perform necessary conversion.  In this code, the "ptraj"
 *      interface (ptrajTokenlist) uses this functionality.

 *  (2) If the number of arguments > 0 then only this many arguments are
 *      grabbed off of the line of text.  Based on the value of the enumerated
 *      type "dispatch", argument processing and conversion are performed prior
 *      to calling the user function.  Processing in this manner allows multiple
 *      triggers to be referenced on the same line of text.  In this code,
 *      the "rdparm" interface (rdparmTokenlist) is implemented in this
 *      manner.
 *
 *  To add a new function to a particular Token structure, simply modify the 
 *  particular list.  Make sure that the list remains NULL terminated and
 *  furthermore that the ``minimum characters to match'' do not overlap with 
 *  any other functions (i.e. you don't want to map a single trigger string to
 *  two different functions since only the first will match).  If you are using
 *  the dispatchToken() routine to convert string arguments and you would like
 *  to implement a new argument types, modify the DispatchArgType enumerated 
 *  type to represent these arguments and modify the dispatchToken() routine to 
 *  process these arguments.  Easier is probably having the function called 
 *  process the arguments off of an argumentStack (as is done with the "ptraj"
 *  functionality)...
 *
 *  To process information on the argument stack, the following routines are 
 *  defined:
 *
 *  getArgumentString(): Return the next element on the stack as a string.
 *
 *     char  *getArgumentString(stackType **, char *)
 *
 *  getArgumentStringLower(): Return the next element on the stack as a
 *     lower case string.
 *
 *     char  *getArgumentStringLower(stackType **, char *);
 *
 *  getArgumentInteger(): Return the next element on the stack as an integer
 *  getArgumentFloat():                                           a float
 *  getArgumentDouble():                                          a double
 *
 *     int    getArgumentInteger(stackType **, int);
 *     float  getArgumentFloat(stackType **, float);
 *     double getArgumentDouble(stackType **, double);
 *
 *  argumentStringContains(): Remove the next element on the stack if it 
 *     contains a specific string, returning 1 otherwise do not perturb the
 *     stack and return 0.
 *
 *     int    argumentStringContains(stackType **, char *);
 *
 *  argumentStackContains(): Like the previous routine, except search the 
 *     whole stack.
 *
 *     int    argumentStackContains(stackType **, char *);
 *
 *  argumentStackKeyToXXX(): If the stack contains a particular string,
 *     grab the following argument on the stack and convert it to a string,
 *     integer, float or double as appropriate.  Both elements are removed
 *     from the stack.  If the specific string is not on the stack, do
 *     not perturb the stack and return the specified default value.
 *  
 *     char  *argumentStackKeyToString(stackType **, char *, char *);
 *     int    argumentStackKeyToInteger(stackType **, char *, int);
 *     float  argumentStackKeyToFloat(stackType **, char *, float);
 *     double argumentStackKeyToDouble(stackType **, char *, double);
 */



/*
 *  Enumerated type representing the expected arguments 
 *  for the functions in the Token structure for use when
 *  the expected number of arguments is > 0.
 *
 *  DP_NOARGS             --- no arguments
 *  DP_STRING             --- a string argument
 *  DP_STR_STRINT         --- string and an integer within a string
 *  DP_STR_STRING         --- two strings
 *  DP_OPTSTR             --- optionally a string argument
 */

enum DispatchArgType { DP_NOARGS, DP_STRING, DP_STR_STRINT,
		       DP_STR_STRING, DP_OPTSTR };


/*
 *  Define the set of Token structures for your application.
 *  (these should perhaps be done on the fly to allow aliasing)
 *
 */

/*
 *  The following is the Token structure for the "main" menu in
 *  the rdparm program.  Make sure that if you add new triggers that
 *  you check the other triggers and min to match to make sure you
 *  can't also match another trigger...
 *
 * "trigger", min to match, # args, DispatchArgType, function pointer
 */

Token rdparmTokenlist[] = {
  {"?",                    1,  1,  DP_OPTSTR,       help},
  {"angles",               4,  1,  DP_OPTSTR,       printAngles},
  {"atoms",                4,  1,  DP_OPTSTR,       printAtomInfo},
  {"bonds",                4,  1,  DP_OPTSTR,       printBonds},
  {"checkcoordinates",     6,  1,  DP_STRING,       checkCoordinatesWrap},
  {"checkgrid",            6,  1,  DP_STRING,       checkGrid},
  {"checkmask",            6,  1,  DP_STRING,       checkAtomMask},
  {"checkvelocity",        6,  1,  DP_NOARGS,       checkVelocity},
  {"charge",               4,  1,  DP_OPTSTR,       chargeOnAtoms},
  {"count",                4,  1,  DP_OPTSTR,       countAtoms},
  {"delete",               6,  2,  DP_STR_STRINT,   deleteBAD},
  {"delperturbed",         6,  2,  DP_STR_STRINT,   deletePerturbedBAD},
  {"dihedrals",            4,  1,  DP_OPTSTR,       printDihedrals},
  {"dumpdihedrals",        5,  1,  DP_OPTSTR,       dumpDihedrals},
  {"exit",                 4,  0,  DP_NOARGS,       quit},
  {"help",                 4,  1,  DP_OPTSTR,       help},
  {"info",                 4,  0,  DP_NOARGS,       parmInfo},
  {"mardi2sander",         6,  1,  DP_STRING,       doMardi2Sander},
  {"modifyboxinfo",        7,  0,  DP_NOARGS,       modifyBoxInfo},
  {"modifymoleculeinfo",   7,  0,  DP_NOARGS,       modifyMolInfo},
  {"newtransform",         5,  1,  DP_STRING,       ptraj},
  {"openparm",             4,  1,  DP_STRING,       getParm},
  {"parminfo",             5,  0,  DP_NOARGS,       parmInfo},
  {"pertbonds",            5,  1,  DP_OPTSTR,       printPerturbedBonds},
  {"perturbedbonds",      10,  1,  DP_OPTSTR,       printPerturbedBonds},
  {"printbonds",           6,  1,  DP_OPTSTR,       printBonds},
  {"printangles",          7,  1,  DP_OPTSTR,       printAngles},
  {"pertangles",           5,  1,  DP_OPTSTR,       printPerturbedAngles},
  {"pertdihedrals",        5,  1,  DP_OPTSTR,       printPerturbedDihedrals},
  {"perturbedangles",     10,  1,  DP_OPTSTR,       printPerturbedAngles},
  {"perturbeddihedrals",  10,  1,  DP_OPTSTR,       printPerturbedDihedrals},
  {"printatoms",           7,  1,  DP_OPTSTR,       printAtomInfo},
  {"printdelphi",          7,  1,  DP_OPTSTR,       printDelphiCharge},
  {"printdihedrals",       7,  1,  DP_OPTSTR,       printDihedrals},
  {"printexcluded",        6,  0,  DP_STRING,       printExcluded},
  {"printlennardjones",    6,  0,  DP_NOARGS,       printLJ},
  {"printtypes",           6,  0,  DP_NOARGS,       printAtomTypes},
  {"prnlev",               6,  1,  DP_OPTSTR,       prnlevSet},
  {"ptraj",                5,  1,  DP_STRING,       ptraj},
  {"quit",                 1,  0,  DP_NOARGS,       quit},
  {"rdparm",               6,  1,  DP_STRING,       getParm},
  {"readparm",             4,  1,  DP_STRING,       getParm},
  {"restrain",             4,  1,  DP_STRING,       restrainBAD},
  {"scalebox",             6,  1,  DP_STRING,       translateBox},
  {"stripwater",           6,  0,  DP_NOARGS,       testWater},
  {"system",               6,  1,  DP_STRING,       doSystem},
  {"testit",               6,  0,  DP_NOARGS,       testit},
  {"testwater",            5,  0,  DP_NOARGS,       testWater},
  {"translatebox",        10,  1,  DP_STRING,       translateBox},
  {"translaterestart",    10,  1,  DP_STRING,       translateRestart},
  {"verbose",              3,  0,  DP_NOARGS,       verbosity},
  {"writeparm",            4,  1,  DP_STRING,       putParm},
  {NULL,                   0,  0,          0,       NULL}
};



/*
 *  The "ptraj" command list.  Currently these commands are either 
 *  processed by ptrajSetupIO (for input/output operations) or
 *  by ptrajSetup() (when an action on the trajectory is
 *  performed).  In the latter case, the argumentStack is simply 
 *  passed to the appropriate action routine (based on the type, i.e.
 *  TRANSFORM_CENTER ==> transformCenter() to be parsed.  See
 *  ptrajSetup() in ptraj.c and actions.c for more information.
 *
 * "trigger", min to match, # args, DispatchArgType, function pointer
 */

Token ptrajTokenlist[] = {
  {"acceptor",       4, -1, TRANSFORM_ACCEPTOR,        ptrajSetupIO},
  {"analyze",        4, -1, TRANSFORM_ANALYZE,         ptrajSetupAnalyze},
  {"angle",          4, -1, TRANSFORM_ANGLE,           ptrajSetup},
  {"atomicfluct3D", 12, -1, TRANSFORM_ATOMICFLUCT3D,   ptrajSetup},
  {"atomicfluct",    5, -1, TRANSFORM_ATOMICFLUCT,     ptrajSetup},
  {"average",        5, -1, TRANSFORM_AVERAGE,         ptrajSetup},
  {"benchmark",      5, -1, TRANSFORM_BENCHMARK,       ptrajSetupIO},
  {"box",            3, -1, TRANSFORM_BOX,             ptrajSetupIO},
  {"center",         5, -1, TRANSFORM_CENTER,          ptrajSetup},
  {"checkoverlap",   6, -1, TRANSFORM_CHECKOVERLAP,    ptrajSetup},
  {"closestwaters",  5, -1, TRANSFORM_CLOSESTWATERS,   ptrajSetup},
  {"clusterdihedr",  8, -1, TRANSFORM_DIHEDRALCLUSTER, ptrajSetup},
  {"clusterattribu", 8, -1, TRANSFORM_CLUSTERATTRIBUTE,ptrajSetup},
  {"cluster",        5, -1, TRANSFORM_CLUSTER,         ptrajSetup},
  {"correlation",    6, -1, TRANSFORM_CORRELATION,     ptrajSetup},
  {"contacts",       6, -1, TRANSFORM_CONTACTS,        ptrajSetup},
  {"diffusion",      4, -1, TRANSFORM_DIFFUSION,       ptrajSetup},
  {"dihedral",       4, -1, TRANSFORM_DIHEDRAL,        ptrajSetup},
  {"dipole",         4, -1, TRANSFORM_DIPOLE,          ptrajSetup},
  {"distance",       4, -1, TRANSFORM_DISTANCE,        ptrajSetup},
  {"dnaiontracker",  7, -1, TRANSFORM_DNAIONTRACKER,   ptrajSetup},
  {"donor",          4, -1, TRANSFORM_DONOR,           ptrajSetupIO},
  {"echo",           4, -1, TRANSFORM_ECHO,            ptrajSetup},
  {"energy",         4, -1, TRANSFORM_ENERGY,          ptrajSetup},
  {"grid",           4, -1, TRANSFORM_GRID,            ptrajSetup},
  {"hbond",          4, -1, TRANSFORM_HBOND,           ptrajSetup},
  {"image",          5, -1, TRANSFORM_IMAGE,           ptrajSetup},
  {"matrix",         4, -1, TRANSFORM_MATRIX,          ptrajSetup},
  {"principal",      5, -1, TRANSFORM_PRINCIPAL,       ptrajSetup},
  {"prnlev",         6, -1, TRANSFORM_PRNLEV,          ptrajSetupIO},
  {"projection",    10, -1, TRANSFORM_PROJECTION,      ptrajSetup},
  {"pucker",         4, -1, TRANSFORM_PUCKER,          ptrajSetup},
  {"radial",         4, -1, TRANSFORM_RADIAL,          ptrajSetup},
  {"radgyr",         4, -1, TRANSFORM_RADIUSOFGYRATION,ptrajSetup},
  {"randomizeions",  10,-1, TRANSFORM_RANDOMIZEIONS,   ptrajSetup},
  {"rms2d",          5, -1, TRANSFORM_2DRMS,           ptrajSetup},
  {"rms",            3, -1, TRANSFORM_RMS,             ptrajSetup},
  {"runningaver",    3, -1, TRANSFORM_RUNNINGAVERAGE,  ptrajSetup},
  {"scale",          5, -1, TRANSFORM_SCALE,           ptrajSetup},
  {"secstruct",      3, -1, TRANSFORM_SECONDARYSTRUCT, ptrajSetup},
  {"strip",          5, -1, TRANSFORM_STRIP,           ptrajSetup},
  {"transform",      6, -1, TRANSFORM_TRANSFORM,       ptrajSetup},
  {"translate",      6, -1, TRANSFORM_TRANSLATE,       ptrajSetup},
  {"truncoct",       5, -1, TRANSFORM_TRUNCOCT,        ptrajSetup},
  {"test",           4, -1, TRANSFORM_TEST,            ptrajSetup},
  {"torsion",        4, -1, TRANSFORM_DIHEDRAL,        ptrajSetup},
  {"unwrap",         4, -1, TRANSFORM_UNWRAP,          ptrajSetup},
  {"vector",         4, -1, TRANSFORM_VECTOR,          ptrajSetup},
  {"watershell",     6, -1, TRANSFORM_WATERSHELL,      ptrajSetup},
  {"2drms",          5, -1, TRANSFORM_2DRMS,           ptrajSetup},
  {"go",             2, -1, TRANSFORM_TRANSFORM,       ptrajSetup},
  {"reference",      5, -1, TRANSFORM_REFERENCE,       ptrajSetupIO},
  {"solvent",        4, -1, TRANSFORM_SOLVENT,         ptrajSetupIO},
  {"trajin",         5, -1, TRANSFORM_TRAJIN,          ptrajSetupIO},
  {"trajout",        5, -1, TRANSFORM_TRAJOUT,         ptrajSetupIO},
  {NULL,             0,  0, TRANSFORM_NOOP,            NULL}
};




/*
 *     Token *
 *  searchTokenList( Token *token, char *text )
 *
 *  ...searches the Token *token structure for a trigger matching
 *  the text.  When found, a pointer to this Token entry is 
 *  returned.  If not found, this implied you've hit the NULL
 *  terminator at the end of the Token list (which is critical!)
 */


   Token *
searchTokenList( Token *token, char *text ) 
{
  Token *t;
  int minMatch;

  text = toLowerCase(text);
  for ( t=token; t->str != NULL; t++ ) {
    minMatch = strlen(t->str);
    if ( (t->minMatch > 0) && (t->minMatch < minMatch) )
      minMatch = t->minMatch;
    if ( strncmp(text, t->str, minMatch) == 0 ) 
      return(t);
  }
  warning("dispatchToken", "Token string \"%s\" not found in tokenlist\n", 
	  text);
  return( (Token *) NULL );
}


/*
 *     char *
 *  processArgument(char *, stackType **)
 *
 *  This routine parses the string typed by the user and places the
 *  string on the stackType *argumentStack.  Argument strings may be 
 *  quoted.  This routine is kind of ugly due to the checking for
 *  quoted strings...
 */

   char *
processArgument(char *strp, stackType **argumentStack) 
{
  int length;
  int isQuoted = 0;
  char *argument;
  char *null_string = "";
  char quoted;

  /*
   * check if the leading character is indicative of a quoted string
   */

  quoted = strp[0];
  switch ( quoted ) {
  case '"':
  case '\'':
  case '`':
    isQuoted = 1;
    break;
  }

  /*
   * if it is quoted and missing an end quote, ignore the leading
   * quote...
   */

  if ( isQuoted ) 
    if ( strchr( strp+1, quoted ) == NULL ) {
      warning("dispatchToken", "Missing end quote on quoted string %s\n",
	      strp);
      isQuoted = 0;
      strp = strp+1;
    }     

  /* 
   * process the string.  Note that a new copy of the string is
   * allocated for further use on the stack rather than a pointer
   * to the crap passed in since this may be freed at any time...
   * Note that the string is push'ed to the bottom of the stack
   * since we want the arguments in the proper order and since
   * we know there is nothing on "top" of the stack that is bad since
   * the stack was cleared prior to entry of this routine...
   */
  
  if ( isQuoted ) {
    for (length=1; strp[length] != quoted; length++);
    if (length == 1) {
      argument = safe_malloc( sizeof(char) * (strlen(null_string) + 1) );
      argument = (char *) strcpy(argument, (char *) null_string);
      warning("dispatchToken", 
	      "Expecting argument (in quotes), none found! (using \"\")\n");
      strp = strp + 2;
    } else {
      argument = safe_malloc( sizeof(char) * (length + 2) );
      strncpy(argument, strp+1, length-1);
      argument[length] = '\n';
      argument[length+1] = (char) 0;
      strp = strp+length+1;
    }
  } else {
    for (length=0; (strlen( strp+length ) > 0) && 
	         ! ( isspace( strp[length] ) || (strp[length] == '\n') ); 
	 length++);
    if ( length == 0 ) {
      argument = safe_malloc( sizeof(char) * (strlen(null_string) + 1) );
      argument = (char *) strcpy(argument, (char *) null_string);
    } else {
      argument = safe_malloc( sizeof(char) * (length+1) );
      sscanf(strp, "%s", argument);
      strp = strp+length;
    }
  }
  pushBottomStack( argumentStack, (void *) argument);
  return( strp );
}  


/*
 *  argumentStringContains(argumentStack, match)
 *
 *  grab the next argument off of the argument stack and check if
 *  it matches the string "match" up to the length of the match
 *  string.  If there is a match, free the argument and return 1.
 *  If it does not match, push the argument back on the stack and 
 *  return 0.  
 *
 *  THE MATCHING IS CASE INSENSITIVE!!!
 */

   int
argumentStringContains(stackType **argumentStack, char *match)
{
  char *argument;

  argument = popStack(argumentStack);
  if (argument == NULL) {
    return 0;
  } else if ( stringMatch(argument, match) ) {
    safe_free(argument);
    return 1;
  } else {
    pushStack(argumentStack, argument);
    return 0;
  }
}


/*
 *  Like argumentStringContains above except that the entire stack is 
 *  searched.
 */

   int
argumentStackContains(stackType **argumentStack, char *match)
{
  char *argument;
  int isMatch;
  stackType *a;

  a = NULL;
  isMatch = 0;

  while ( (argument = popStack(argumentStack)) != NULL ) {

    if ( isMatch == 0 && stringMatch(argument, match) ) {
      isMatch = 1;
      safe_free(argument);
    } else {
      pushBottomStack(&a, argument);
    }
  }
  *argumentStack = a;
  return(isMatch);
}

/*
 *  Grab a string off the argumentStack and return it
 */

   char *
getArgumentString(stackType **argumentStack, char *defvalue)
{
  char *returnValue;

  returnValue = (char *) popStack(argumentStack);
  if (returnValue == NULL || returnValue[0] == (char) 0) {

    if (defvalue != NULL) {
      returnValue = (char *) safe_malloc(sizeof(char) * (strlen(defvalue)+1));
      strcpy(returnValue, defvalue);
    } else
      returnValue = NULL;
  }
  return(returnValue);
}


/*
 *  Grab a lower case string off the argumentStack and return it
 */

   char *
getArgumentStringLower(stackType **argumentStack, char *defvalue)
{
  char *returnValue;

  returnValue = (char *) popStack(argumentStack);
  if (returnValue == NULL || returnValue[0] == (char) 0) {

    if (defvalue != NULL) {
      returnValue = (char *) safe_malloc(sizeof(char) * (strlen(defvalue)+1));
      strcpy(returnValue, defvalue);
    } else
      returnValue = NULL;
  }
  if (returnValue != NULL) {
    returnValue = toLowerCase(returnValue);
  }
  return(returnValue);
}


/*
 *  Grab a string off the argument stack, convert it to an integer,
 *  free up the entry and return the integer
 */

   int
getArgumentInteger(stackType **argumentStack, int defvalue)
{
  char *argument;
  int returnValue;

  argument = popStack(argumentStack);
  if ( argument == NULL || argument[0] == (char) 0 ) {
    returnValue = defvalue;
  } else {
    if (sscanf(argument, "%i", &returnValue) < 1) {
      warning("getArgumentInteger()", "Scanning integer argument failed");
    }
  }
  safe_free(argument);
  return( returnValue );
}


/*
 *  Grab a string off the argument stack, convert it to an float,
 *  free up the entry and return the float
 */

   float
getArgumentFloat(stackType **argumentStack, float defvalue)
{
  char *argument;
  float returnValue;

  argument = popStack(argumentStack);
  if (argument == NULL || argument[0] == (char) 0) {
    returnValue = defvalue;
  } else {

    if (sscanf(argument, "%f", &returnValue) < 1) {
      warning("getArgumentFloat()", "Scanning floating point argument failed");
    }
  }
  safe_free(argument);
  return( returnValue );
}

/*
 *  Grab a string off the argument stack, convert it to an double,
 *  free up the entry and return the double
 */


   double
getArgumentDouble(stackType **argumentStack, double defvalue)
{
  char *argument;
  double returnValue;

  argument = popStack(argumentStack);
  if (argument == NULL || argument[0] == (char) 0) {
    returnValue = defvalue;
  } else {

    if (sscanf(argument, "%lf", &returnValue) < 1) {
      warning("getArgumentDouble()", "Scanning double precision argument failed");
    }
  }
  safe_free(argument);
  return( returnValue );
}


/*
 *  If the string "match" is on the stack, grab the argument that directly follows
 *  it, return it as a string, and free up the two entries on the stack, otherwise
 *  do not perturb the stack.  
 *  Matching is case insensitive and based on the length of the "match" string.
 */

   char *
argumentStackKeyToString(stackType **argumentStack, char *match, char *defvalue)
{
  char *argument;
  stackType *a;
  char *returnValue;
  int notset;

  a = NULL;

  notset = 1;
  returnValue = defvalue;

  while ( (argument = popStack(argumentStack)) != NULL ) {

    if ( notset && stringMatch(argument, match) ) {
      safe_free(argument);
      returnValue = getArgumentString(argumentStack, defvalue);
      notset = 0;
    } else {
      pushBottomStack(&a, argument);
    }
  }
  *argumentStack = a;
  return(returnValue);
}

/*
 *  Like argumentStackKeyToString above except that the next argument on the
 *  stack is assumed to be an integer
 */
   int
argumentStackKeyToInteger(stackType **argumentStack, char *match, int defvalue)
{
  char *argument;
  stackType *a;
  int returnValue;
  int notset;

  a = NULL;
  notset = 1;
  returnValue = defvalue;
  while ( (argument = popStack(argumentStack)) != NULL ) {

    if ( notset && stringMatch(argument, match) ) {
      safe_free(argument);
      returnValue = getArgumentInteger(argumentStack, defvalue);
      notset = 0;
    } else {
      pushBottomStack(&a, argument);
    }
  }
  *argumentStack = a;
  return(returnValue);
}


/*
 *  Like argumentStackKeyToString above except that the next argument on the
 *  stack is assumed to be a float
 */

   float
argumentStackKeyToFloat(stackType **argumentStack, char *match, float defvalue)
{
  char *argument;
  stackType *a;
  float returnValue;
  int notset;

  a = NULL;
  notset = 1;
  returnValue = defvalue;
  while ( (argument = popStack(argumentStack)) != NULL ) {

    if ( notset && stringMatch(argument, match) ) {
      safe_free(argument);
      returnValue = getArgumentFloat(argumentStack, defvalue);
      notset = 0;
    } else {
      pushBottomStack(&a, argument);
    }
  }
  *argumentStack = a;
  return(returnValue);
}


/*
 *  Like argumentStackKeyToString above except that the next argument on the
 *  stack is assumed to be a double
 */

   double
argumentStackKeyToDouble(stackType **argumentStack, char *match, double defvalue)
{
  char *argument;
  stackType *a;
  double returnValue;
  int notset;

  a = NULL;
  notset = 1;
  returnValue = defvalue;
  while ( (argument = popStack(argumentStack)) != NULL ) {

    if ( notset && stringMatch(argument, match) ) {
      safe_free(argument);
      returnValue = getArgumentDouble(argumentStack, defvalue);
      notset = 0;
    } else {
      pushBottomStack(&a, argument);
    }
  }
  *argumentStack = a;
  return(returnValue);
}







/* 
 *  dispatchToken( Token *list, stackType *argumentStack, char *text)
 *
 *  This routine is essentially the main routine in the dispatch
 *  abstraction.
 *
 *    (1) Grab a while space delineated trigger string out of the 
 *        "text" and search the Token "list" using 
 *        searchTokenList().
 *    (2) If a match is found, use the information returned from the
 *        token to dump arguments onto a stack using
 *        processArguments().
 *    (3) Convert the arguments if necessary
 *    (4) Call the dispatch routine "fxn".
 */

   void
dispatchToken( Token *tokenlist, stackType *argumentStack, char *text)
{
  int skip;
  char *strp;
  Token *t;
  char trigger[LINE_SIZE];
  char *str1, *str2;
  int int1;

  t = NULL;
  if ( text == NULL ) {
    warning("dispatchToken", 
	    "Passed NULL string as argument, returning...\n");
    return;
  }

  /*
   *  skip all the leading whitespace
   */
  skipWhitespace(text);

  /*
   *  move over each whitespace delineated word in the input string
   *  grabbing the trigger strings (to search for) or arguments to
   *  process until the text is exhausted
   */
  for (strp = text; strcmp( strp, "" ) != 0; ) {

    /*
     *  scan in a "trigger"
     */
    if ( sscanf(strp, "%s", trigger) != 1 )
      error("dispatchToken()", "Error scanning trigger from string (%s)\n",
	    strp);

#ifdef DEBUG_DISPATCH
    fprintf(stderr, "dispatchToken(): trigger is (%s)\n", trigger);
#endif

    /*
     *  search the tokenlist for the trigger.  If the trigger is
     *  found, increment the string pointer strp to jump over the
     *  current trigger
     */

    t = searchTokenList( tokenlist, trigger );
    strp = strp + strlen(trigger);

    /*
     *  if a valid match was found process the arguments as necessary
     *  onto the argument stack.
     */

    if (t != NULL) {
      
      if (argumentStack != NULL)         /* clean up stack if necessary */
	clearStack( &argumentStack ); 

      if (t->numArg < 0) {

	/*
	 *  If the expected argument count is < 0 treat the rest of the 
	 *  text in strp as arguments, place on argument stack, call the
	 *  associated function and return
	 */

	do {
	  skipWhitespace( strp );
	  strp = processArgument(strp, &argumentStack);
	} while ( strp[0] != (char) 0 );

	t->fxn(argumentStack, t->dispatch);
	return;
      } else {

	/*
	 *  Process only as much text is expected based on the number
	 *  of expected arguments
	 */

	for (skip = t->numArg; skip > 0; skip--) {
	  skipWhitespace( strp );
	  strp = processArgument(strp, &argumentStack);
	}

	switch ( t->dispatch ) {

	case DP_NOARGS:
	  t->fxn(); 
	  break;

	case DP_STRING:
	case DP_OPTSTR:
	  str1 = popStack( &argumentStack );
	  t->fxn(str1);
	  safe_free( (void *) str1 );
	  break;

	case DP_STR_STRINT:
	  str1 = popStack( &argumentStack );
	  str2 = popStack( &argumentStack ); 
	  if ( sscanf(str2, "%d", &int1) != 1)
	    error("dispatch()", "Failure on scan of integer from %s\n", str2);
	  t->fxn(str1, int1);

	  safe_free( (void *) str1 ); 
	  safe_free( (void *) str2 );
	  break;

	case DP_STR_STRING:
	  str1 = popStack( &argumentStack );
	  str2 = popStack( &argumentStack ); 

	  t->fxn(str1, str2);

	  safe_free( (void *) str1 ); 
	  safe_free( (void *) str2 );
	  break;

	default:
	  warning("dispatchToken()", "DispatchArgType (%i) not in list, cannot execute\n",
		  t->dispatch);
	}
      }
    }

    skipWhitespace ( strp );

  }
}


