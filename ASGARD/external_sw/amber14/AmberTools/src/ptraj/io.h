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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/io.h,v 10.3 2009/08/24 23:58:48 rcw Exp $
 *
 *  Revision: $Revision: 10.3 $
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



typedef struct _fileType {
  char *buffer;
  FILE *file;
  char *mode;
  int popen;
} fileType;

#ifdef IO_MODULE

/*
 * GLOBAL: fileStack -- contains a circular list fileType structures in
 * order to maintain a list of opened files and respective filenames, 
 * modes, and popen status.
 *
 */

stackType *fileStack = (stackType *) NULL;

#  ifdef __STDC__
extern int safe_fclose_buffer(char *);
#  else
extern int safe_fclose_buffer();
#  endif

#else /* IO_MODULE */

extern stackType *fileStack;

#  ifdef __STDC__

extern void doSystem(char *);
extern char * promptToOpenFile( FILE **, char *, char *, char *);
extern long long int gzipFileSize(char *);
extern long long int bzip2FileSize(char *);
extern long long int zipFileSize(char *);
extern int id_Filesig(char *,FILE *);
extern int openFile( FILE **, char *, char *);
extern int  promptUserResponse(FILE *, FILE *, char *, char *, int);
extern char * promptUser(FILE *, FILE *, char *);
extern int safe_fclose(FILE *);
extern int safe_fclose_buffer(char *);
extern int safe_close(FILE *, char *);
extern void fileStack_clear();
extern FILE *safe_open(FILE *, char *, int);
extern FILE *safe_fopen(char *, char *);
extern FILE *safe_popen(FILE *, char *);
extern FILE *safe_freopen(FILE *);
extern void printfone(char *, ...);
     
#  else

extern void doSystem();
extern char * promptToOpenFile();
extern long long int gzipFileSize();
extern long long int bzip2FileSize();
extern long long int zipFileSize();
extern int id_Filesig();
extern int openFile();
extern int  promptUserResponse();
extern char * promptUser();
extern int safe_fclose();
extern int safe_fclose_buffer();
extern int safe_close();
extern void fileStack_clear();
extern FILE *safe_open();
extern FILE *safe_fopen();
extern FILE *safe_popen();
extern FILE *safe_freopen();
extern void printfone();
     
#  endif

#endif /* IO_MODULE */

