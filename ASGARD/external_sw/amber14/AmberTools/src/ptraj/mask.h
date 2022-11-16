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
 *  AMBER CVS tracking:
 *
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/mask.h,v 10.0 2008/04/15 23:24:11 case Exp $ 
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
 *  Enhanced atom selection parser: Viktor Hornak, Stony Brook University
 */
#define  MAXSELE 1000
#define  ALL      0
#define  NUMLIST  1
#define  NAMELIST 2
#define  TYPELIST 3
#define  ELEMLIST 4

typedef struct StackNodeTag {
  char *entry;        
  struct StackNodeTag  *next;         
} StackNode;

extern char * parseMaskString(char *, int, int, Name *, Name *, int *, void *, void *, void *, char);

#  ifdef __STDC__
int isEmpty (void);
StackNode * MakeNode(char *);
void Push(char *);
char * Pop(void);
char * Retrieve(void);
void  PrintStack(void);
int isOperator(char);
int isOperand(char);
int priority(char);
void tokenize(char *, char *);
void torpn(char *, char *);

char * eval(char *, int, int, Name *, Name *, int *, void *, void *, void *, char);
char * selectDist(char *, char *, int, int, int *, void *, void *, void *, char);
char * binop(char, char *, char *, int);
char * neg(char *, int);
int isElemMatch(char *, char *);
int isNameMatch(char *, char *); 
void resnum_select(int, int, char *, int, int *);
void resname_select(char *, char *, int, Name *, int *);
void all_select(char *, int);
void atnum_select(int, int, char *, int);
void atname_select(char *, char *, int, Name *);
void attype_select(char *, char *, int, Name *);
void atelem_select(char *, char *, int, Name *);
void residue_numlist(char *, char *, int, int *);
void residue_namelist(char *, char *, int, Name *, int *);
void atom_numlist(char *, char *, int);
void atom_namelist(char *, char *, int, Name *);
void atom_typelist(char *, char *, int, Name *);
void atom_elemlist(char *, char *, int, Name *);
char * selectElemMask(char *, int, int, Name *, Name *, int *);
char * parseMaskString(char *, int, int, Name *, Name *, int *, void *, void *, void *, char);
#  else
int isEmpty ();
StackNode * MakeNode();
void Push();
char * Pop();
char * Retrieve();
void  PrintStack();
int isOperator();
int isOperand();
int priority();
void tokenize();
void torpn();

char * eval();
char * binop();
char * neg();
int isElemMatch();
int isNameMatch(); 
void resnum_select();
void resname_select();
void all_select();
void atnum_select();
void atname_select();
void attype_select();
void atelem_select();
void residue_numlist();
void residue_namelist();
void atom_numlist();
void atom_namelist();
void atom_typelist();
void atom_elemlist();
char * selectElemMask();
char * parseMaskString();
#  endif
#define MASK_MODULE
#ifndef MASK_MODULE 

/*
 * Prototypes for externally visible functions
 */
#  ifdef __STDC__
extern char * parseMaskString(char *, int, int, Name *, Name *, int *, void *, void *, void *, char);
#  else
extern char * parseMaskString();
#  endif

#endif /* ifndef MASK_MODULE */

