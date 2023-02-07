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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/utility.c,v 10.4 2010/01/28 17:03:53 droe Exp $
 *
 *  Revision: $Revision: 10.4 $
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <stdarg.h>

#define UTILITY_MODULE
#include "utility.h"

/*
 *  The following functions are defined (see also the supporting
 *  header file, utility.h):
 *
 *
 *  error(char *function, char *format, ... )
 *     ...a varargs function that prints out the first argument 
 *     (which typically represents the calling function name for
 *     debugging) followed by the error string which may contain 
 *     printf format characters, followed by an arbitrary number of
 *     arguments representing the printf style formats.  
 *     Execution of the program is terminated when this function 
 *     is called.
 *
 *  warning(char *function, char *format, ... )
 *     ...same as above, but execution is not terminated.
 *
 *     void *
 *  safe_malloc(size)
 *     ... allocates size bytes of memory, initializes the values
 *     to zero, and performs error checking.
 *
 *     void *
 *  safe_realloc(memory, cur_size, increase)
 *     ...a safe realloc command, increasing the space pointed to 
 *     by memory by increase size.
 *
 *  safe_free(memory)
 *     ...frees the memory pointed to by memory (with checking 
 *     to make sure the pointer isn't NULL
 *
 *  scanString(), scanInt(), scanDouble()
 *     ...these functions scan values from a file, printing a 
 *     string if an error occurs.  See the functions themselves for 
 *     more information.
 *
 *     char *
 *  toLowerCase(string)
 *     ...converts the string to lower case
 *
 *  shiftArray(src, dest, length)
 *     ...shifts the src array by length elements putting the result
 *     in the dest array.
 *
 *  pushStack( stackType **stackp, void *entry )
 *     ...push a pointer to "entry" on the stackType ** stack or
 *     &stack
 *
 *  pushBottomStack( stackType **stackp, void *entry )
 *     ...push a pointer to "entry" on to the *bottom* of the
 *     stackType ** structure &stack
 *
 *     void *
 *  popStack( stackType **stackp )
 *     ...pop the top element off the &stack and delete the
 *     associated stack element
 *
 *  clearStack( stackType **stackp )
 *     ...clear out the whole &stack and free up the associated
 *     memory
 *
 *  printStack( stackType **stackp, (*fxnPrintStackEntry)(void *), char * )
 *     ...a debugging routine to print out a stack; the function
 *     passed in (fxnPrintStackEntry)() is set up to print the
 *     entry of interest.  The (char *) is optional text to print before
 *     the entry (along with the stack entry number).
 *
 *     char *
 *  copyString( char *input )
 *     ...returns a copy of the string "input" (after allocating new space
 *     and copying the contents)...
 *
 */


void error(char *function, char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  fprintf(stderr, "\nERROR in %s: ", function);
  vfprintf(stderr, fmt, args);
  fprintf(stderr, "\n" );
  va_end(args);
  exit(1);
}

void warning(char *function, char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  fprintf(stderr, "\nWARNING in %s: ", function);
  vfprintf(stderr, fmt, args);
  fprintf(stderr, "\n" );
  va_end(args);
}

   void *
SafeMalloc(char* filename, int line, size_t size)
{
  void *mem;
  if ( size ) {
    if ( (mem = (void *) malloc( size )) == NULL ) {
      fprintf(stderr, "Error in file %s line %i when calling SafeMalloc() to alloc of %i bytes.\n ", filename, line, (int) size);
    }
    mem = memset(mem, (char) 0, size);
  } else 
    mem = NULL;
  return( mem );
}

   void *
safe_malloc(size_t size)
{
  void *mem;
  if ( size ) {
    if ( (mem = (void *) malloc( size )) == NULL )
      error( "safe_malloc", "Error in alloc of %i bytes", size);
    mem = memset(mem, (char) 0, size);
  } else 
    mem = NULL;
  return( mem );
}


/* 
 *  Performs error checking on realloc() and zero initializes new memeory.
 *  NOTE: the semantics of the call are slightly different than realloc()!
 */
   void *
safe_realloc(void *mem, size_t cur_size, size_t increase)
{
  void *temp;

  if ( cur_size == 0 ) 
    return ( safe_malloc( increase ) );
  else if ( (temp = (void *) 
	     realloc((void *) mem, (size_t) (cur_size + increase))) == NULL )
    error( "safe_realloc", "Couldn't allocate %i bytes more", increase);

  /* cast temp to avoid pointer arithmetic on void* which has unknown size */
  /* use char* as a byte level pointer per traditional malloc/realloc usage */
  memset( (char *) temp + cur_size, (char) 0, increase);

  mem = temp;
  return(mem);
}



   void
safe_free(void *pointer)
{
  if ( pointer != NULL ) free(pointer);
}


/*
 *  Scan for a string from the file (fp).  If failure occurs
 *  print a warning message and exit if isError is TRUE.
 */
   int
scanString(FILE *fp, char *ip, int isError, char *err_string) 
{
  if ( ! fscanf( fp, "%s", ip ) ) {
    if ( isError == TRUE ) 
      error("scanString", "...scanning %s\n", err_string);
    else {
      warning("scanString", "...scanning %s\n", err_string);
      return TRUE;
    }
  }
  return FALSE;
}

/*
 *  Scan for a integer from the file (fp).  If failure occurs
 *  print a warning message and exit if Error is TRUE.
 */
   int
scanInt(FILE *fp, int *ip, int isError, char *err_string) 
{
  if ( ! fscanf( fp, "%d", ip ) ) {
    if ( isError == TRUE ) 
      error("scanInt", "...scanning %s\n", err_string);
    else {
      warning("scanInt", "...scanning %s\n", err_string);
      return TRUE;
    }
  }
  return FALSE;
}

/*
 *  Scan for a double from the file (fp).  If failure occurs
 *  print a warning message and exit if isError is TRUE.
 */
   int
scanDouble(FILE *fp, double *ip, int isError, char *err_string) 
{
  if ( ! fscanf( fp, "%lf", ip ) ) {
    if ( isError == TRUE ) 
      error("scanDouble", "...scanning %s\n", err_string);
    else {
      warning("scanDouble", "...scanning %s\n", err_string);
      return TRUE;
    }
  }
  return FALSE;
}


/* Convert the string to lower case.  This has not been tested on
 * all character sets.  For this to work, it is necessary that 
 * there be a constant difference in values between a lower case
 * letter and it's corresponding upper case value.
 */
   char *
toLowerCase( char *string_in ) 
{
  int i, length;

  length = strlen( string_in );
  for (i=0; i < length; i++ ) {
    string_in[i] = tolower( string_in[i] );
  }
  return string_in;
}


/*
 *  stringContains(string, match) -- this is like strstr except that instead
 *  of returning a pointer to the first occurrence of "match" in "string"
 *  1 is returned if the "match" is within "string" otherwise 0 is returned.
 *  AN ADDITIONAL DIFFERENCE is that this matching is not case sensitive!!!
 */

   int
stringContains(char *string, char *match)
{
  char *lowerString, *lowerMatch;

  if (string == NULL || match == NULL) return 0;

  lowerString = safe_malloc( (size_t) sizeof(char) * (strlen(string)+1));
  strcpy(lowerString,string);
  lowerMatch  = safe_malloc( (size_t) sizeof(char) * (strlen(match)+1));
  strcpy(lowerMatch,match);
  lowerString = toLowerCase(lowerString);
  lowerMatch  = toLowerCase(lowerMatch);

  if (strstr(lowerString, lowerMatch) != NULL) {
    safe_free(lowerString);
    safe_free(lowerMatch);
    return 1;
  } else {
    safe_free(lowerString);
    safe_free(lowerMatch);
    return 0;
  }
}

/*
 *  this is like strcmp, however the match is only check on the length
 *  of the string "match" and the comparison is case insenstive
 */

   int
stringMatch(char *string, char *match)
{
  char *lowerString, *lowerMatch;

  if (string == NULL || match == NULL) return 0;

  lowerString = safe_malloc( (size_t) sizeof(char) * (strlen(string)+1));
  strcpy(lowerString,string);
  lowerMatch  = safe_malloc( (size_t) sizeof(char) * (strlen(match)+1));
  strcpy(lowerMatch,match);
  lowerString = toLowerCase(lowerString);
  lowerMatch  = toLowerCase(lowerMatch);

  /*
   *  this should be an exact, not partial match
   *     if (strncmp(lowerMatch, lowerString, strlen(lowerMatch)) == 0) {
   */
  if (strcmp(lowerMatch, lowerString) == 0) {
    safe_free(lowerString);
    safe_free(lowerMatch);
    return 1;
  } else {
    safe_free(lowerString);
    safe_free(lowerMatch);
    return 0;
  }
}


/*
 *  This function matches wild character '*' and '?'. 
 *  '*' can match 0 or more characters, '?' can match any single character.
 *  Wild characters can NOT be in *string, only in *match.
 *  The comparison is case senstive
 */

   int
stringMatchWild(char *string, char *match)
{
  char *s, *m;
  
  if (string == NULL || match == NULL) return 0;
  
  s = string;
  m = match;
  
  while (*s || *m) {
    switch (*m) {
      case '?':
        if (!(*s)) return 0;
        s++; m++;
        break;
      case '*':
        m++;
        while (*s) {
          if (stringMatchWild(s, m) == 1) return 1;
          s++;
        }
        return(stringMatchWild(s, m)); /* Now *s == '\0' */
        break;
      default:
        if (*s != *m) return 0;
        s++; m++;
    }
  }
  return 1; /* when *string == *match == '\0' */
  
}


  /*
   *  stringSplit()
   *
   *  This routine will split str1 by any characters in str2 into a token array.
   *  The returned array, and its contents are newly assigned, so you need to take care of the releasing of memory.
   *  The characters in str2 will not appear in the newly created token array.
   *  Example: stringSplit ("ABC, ,EFG  HI", ", ") will return 
   *  char **string = {"ABC","EFG","HI"}
   */
   char**
stringSplit(char *str1, char *str2)
{
  char **string;
  int i = 0;
  char *c;
  char *oldc;
  int tokens; 
  
  if (*str1) tokens = 1;
  
  for (c = str1, oldc = c; *c; c++) {
    if (strchr(str2, *c)) {
      if (oldc < c) tokens++;
      oldc = c + 1;
    }	
  }
  if (oldc < c) tokens++;

  string = safe_malloc(sizeof(char*) * tokens);
  c = oldc = str1;
  for (c = str1, oldc = c; *c; c++) {
    if (strchr(str2, *c)) {
      if (oldc < c) {
        string[i] = safe_malloc(sizeof(char) * (c - oldc + 1));
        strncpy(string[i], oldc, c - oldc);
        string[i][c-oldc] = '\0';
        i++;
      }
      oldc = c + 1;
    }	
  }
  if (oldc < c) {
    string[i] = safe_malloc(sizeof(char) * (c - oldc + 1));
    strncpy(string[i], oldc, c - oldc);
    string[i][c-oldc] = '\0';
  }
  string[tokens - 1] = NULL;
  
  return string;
}	  

/*
 *  Performs a safety check on a memcpy() call
 */
   void
shiftArray( void *src, void *dest, int length )
{
  if ( length < 0 ) 
    warning("shiftArray", "Cannot currently shift array < 0 bytes (%i)",
	  length);
  else if ( length > 0 ) {
    memcpy( dest, src, (size_t) length );
  } 
}



   void
pushStack( stackType **stackp, void *entry )
{
  stackType *sp;

  sp = safe_malloc( (size_t) sizeof(stackType) );
  sp->entry = entry;
  sp->next = *stackp;
  *stackp = sp;
}


   void
pushBottomStack( stackType **stackp, void *entry )
{
  stackType *sp;

  if ( *stackp == NULL )
    pushStack( stackp, entry );
  else {
    sp = *stackp;
    while ( sp->next != NULL ) sp = sp->next;
    sp->next = safe_malloc( (size_t) sizeof(stackType) );
    sp->next->entry = entry;
    sp->next->next = NULL;
  }
}


   void *
popStack( stackType **stackp )
{
  void *entry;
  stackType *sp;

  if ( *stackp == NULL ) {
    return( (char *) NULL );
  }

  sp = *stackp;
  entry = sp->entry;

  *stackp = sp->next;
  sp->next=NULL;
  free( sp );
  return( entry );
}



   void
clearStack( stackType **stackp ) 
{
  stackType *sp, *tmp;

  if (stackp == NULL) {
    return;
  }

  tmp = NULL;
  for (sp = *stackp; sp != NULL; ) {
    if ( tmp != NULL ) safe_free( (void *) tmp);
    safe_free( (void *) sp->entry);
    sp->entry = NULL;
    tmp = sp;
    sp = sp->next;
  }
  *stackp = NULL;
}


   void
printStack( stackType **stackp, fxnPrintStackEntry fxn, char *babble) 
{
  stackType *p;
  int i;
  for (p = *stackp, i = 1; p != NULL; p = p->next, i++) {
    if ( babble != NULL ) 
      fprintf(stdout, "%s (%i)\n", babble, i);
    fxn(p->entry);
  }

}

   char *
copyString( char *input )
{
  char *output;

  output = (char *) safe_malloc( sizeof(char) * (strlen(input) + 1) );
  if (output != NULL)
    strcpy(output, input);

  return output;
}


