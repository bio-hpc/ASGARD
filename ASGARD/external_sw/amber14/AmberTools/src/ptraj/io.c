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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/io.c,v 10.7 2010/01/27 22:07:26 droe Exp $
 *
 *  Revision: $Revision: 10.7 $
 *  Date: $Date: 2010/01/27 22:07:26 $
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




#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define IO_MODULE
#include "ptraj.h"

/*
 * This routine defines the subroutines used in rdparm/ptraj for
 * various program input and output.
 */


/* -----------------------------------------------------------------

                   SUBROUTINE DEFINITIONS:


   void
doSystem(char *command)
...a wrapper for the system(command) call with no error checking

   FILE *
safe_open(char *buffer, char *mode, int is_popen)
...a routine to open a file of name "buffer" with mode "mode".  If
"is_popen" is not zero, then the popen() command is called instead
of fopen() with the "buffer" and "mode".  This routine maintains a 
fileType list of associated "buffer", "mode", FILE and popen() status
which is used to closing and/or reopening files.  If the file has
previously been opened and is maintained on the fileStack fileType list, 
then the file will be closed and reopened.  Files opened with safe_open
should be closed using safe_close(), safe_fclose() or safe_fclose_buffer()
in order to clean up the fileStack.

   FILE *
safe_freopen(FILE *file)
...reopen a previously "safe" opened file

   FILE *
safe_fopen(char *buffer, char *mode)
...a wrapper to safe_open above that calls it with is_popen == 0

   FILE *
safe_popen(char *buffer, char *mode)
...a wrapper to safe_open above that calls it with is_popen == 1

   int
safe_close(FILE *closed, char *buffer)
...will close FILE "closed" and update the fileStack if necessary.
If only the "buffer" is passed in (representing the file name or
popen command), then the fileStack will be searched for the appropriate
FILE to close.  If the FILE "closed" or "buffer" is present in the 
fileStack, the proper pclose() or fclose() will be called, else
fclose() is called.

   int
safe_fclose(FILE *closed)
...a wrapper to safe_close for the FILE "closed"

   int
safe_fclose_buffer(char *buffer)
n...a wrapper to safe_close() for the buffer "buffer"

promptToOpenFile(fpp, filename, mode, prompt)
... attempts to open filename with mode.  If a failure occurs, the
user is prompted (with prompt) and asked to input a new filename.
This continues until a sucessful open occurs, where upon the file is
returned in fpp.

   int
openFile(fpp, filename, mode)
...attempts to open the "filename" with "mode" returning 1 is successful

   int
promptUserResponse(fpin, fpout, response, min_to_match)
...writes ``prompt'' to fpout, fget()'ing a response from fpin.  The 
default response is in ``response'', where ``min_to_match'' is the
minimum number of characters to match.  1 is returned upon successful 
match, 0 otherwise.

   int
safe_fclose(FILE *closed)
...uses pclose or fclose to close a file in order to smoothly terminate
popen operationa...


-------------------------------------------------------------------------*/



   void
doSystem(char *command)
{
  system(command);
}



#ifdef DEBUG_IO

   void
printFileStack()
{
  stackType *stack;
  fileType *file;
  int i;

  for (stack = fileStack, i = 0; stack != NULL; i++) { 
    file = (fileType *) stack->entry;
    fprintf(stdout, "  printFileStack: file list (%x) entry %i is %x (%s)\n", 
	   stack, i+1, file->file, file->buffer);
    stack = stack->next;
    if ( stack == fileStack ) break;
  }
}

#endif


   fileType *
is_open_file_buffer(char *buffer)
{
  stackType *sp;
  fileType *file;

  if ( buffer == NULL ) return (NULL);
  if ( fileStack == NULL ) return (NULL);

  sp = fileStack;
  do {

    file = (fileType *) sp->entry;
    if ( strcmp(file->buffer, buffer) == 0 ) {

#ifdef DEBUG_IO
      fprintf(stdout, "Buffer search (%x:%s) has previously been %s\n",
	      file->file, file->buffer, (file->popen ? "popened" : "opened"));
#endif
      return (file);
    }
    sp = sp->next;
  } while ( sp != fileStack );

  return (NULL);

}


   stackType *
is_open_file_buffer_stack(char *buffer)
{
  stackType *sp;
  fileType *file;

  if ( buffer == NULL ) return (NULL);
  if ( fileStack == NULL ) return (NULL);

  sp = fileStack;
  do {

    file = (fileType *) sp->entry;
    if ( strcmp(file->buffer, buffer) == 0 ) {
      return (sp);
    }
    sp = sp->next;
  } while ( sp != fileStack );

  return (NULL);

}

   fileType *
is_open_file(FILE *f)
{
  stackType *sp;
  fileType *file;

  if ( f == NULL ) return NULL;
  if ( fileStack == NULL ) return NULL;

  sp = fileStack;
  do {

    file = (fileType *) sp->entry;
    if ( file->file == f ) {

#ifdef DEBUG_IO
      fprintf(stdout, "File search (%x:%s) has previously been %s\n",
	      file->file, file->buffer, (file->popen ? "popened" : "opened"));
#endif
      return (file);
    }
    sp = sp->next;
  } while ( sp != fileStack );

  return (NULL);

}


   stackType *
is_open_file_stack(FILE *f)
{
  stackType *sp;
  fileType *file;

  if ( f == NULL ) return NULL;
  if ( fileStack == NULL ) return NULL;

  sp = fileStack;
  do {

    file = (fileType *) sp->entry;
    if ( file->file == f ) {
      return (sp);
    }
    sp = sp->next;
  } while ( sp != fileStack );

  return (NULL);

}



   FILE *
safe_open(char *buffer, char *mode, int is_popen)
{
  stackType *tmpStack, *sp;
  fileType *file;
  FILE *fp;
  int length;


  /* check to make sure that the file is not already opened and
   * close it if it is...
   */
  if ( is_open_file_buffer(buffer) != NULL ) 
    safe_fclose_buffer(buffer);

  fp = NULL;
  if (is_popen) {
    fp = popen(buffer, mode);
  } else {
    fp = fopen(buffer, mode);
  }

  if ( fp == NULL) 
    return NULL;

  /*
   * file opened successfully, therefore save file information...
   */

  file = (fileType *) safe_malloc(sizeof(fileType));

  file->popen = is_popen;

  length = sizeof(char) * (strlen(buffer) + 1);
  file->buffer = (char *) safe_malloc( length );
  memcpy(file->buffer, buffer, (size_t) length);

  length = sizeof(char) * (strlen(mode) + 1);
  file->mode = (char *) safe_malloc( length );
  memcpy(file->mode, mode, (size_t) length);

  file->file = fp;

  /*
   * create an entry for the fileStack and add it
   */
  tmpStack = (stackType *) safe_malloc(sizeof(stackType));
  tmpStack->entry = (void *) file;
  tmpStack->next = NULL;

  if (fileStack == NULL) {
    /*
     * fileStack hasn't been initialized yet, therefore this is
     * the first entry and we need to make the list circular!
     */

    fileStack = tmpStack;
    fileStack->next = fileStack;

  } else {
    /*
     * fileStack exists, hence insert the current entry as the next
     * entry and make sure to reset the next entry to retain the
     * list integrity
     */

    if (fileStack->next == fileStack) {
      fileStack->next = tmpStack;
      tmpStack->next = fileStack;
    } else {
      sp = fileStack->next;
      fileStack->next = tmpStack;
      tmpStack->next = sp;
    }
  }

#ifdef DEBUG_IO
  fprintf(stdout, 
	  "Returning (%x) from safe_open(%s, %s, %i), fileStack follows:\n",
	  fp, buffer, mode, is_popen);
  printFileStack();
#endif

  return( fp );
}


/* wrappers to the safe_open call */

   FILE *
safe_freopen(FILE *f)
{
  fileType *file;
  FILE *returnFile;
  char *buffer;
  char *mode;

  returnFile = NULL;
  file = is_open_file(f);
  if ( file == NULL )
    warning("safe_freopen()", 
	    "Attempting to reopen a file not in the fileStack!\n");

  /*
   *  the safe_open will close the already open file and
   *  deallocate the pointers file->buffer and file->mode
   *  therefore copies of these need to be passed in which
   *  can be deallocated AFTER the safe_open call...
   */
  buffer = safe_malloc(sizeof(char) * (strlen(file->buffer)+1));
  strcpy(buffer, file->buffer);
  mode = safe_malloc(sizeof(char) * (strlen(file->mode)+1));
  strcpy(mode, file->mode);

  returnFile = safe_open(buffer, mode, file->popen);

  safe_free(buffer);
  safe_free(mode);
  return( returnFile );
}


   FILE *
safe_fopen(char *buffer, char *mode)
{
  return( safe_open(buffer, mode, 0) );
}

   FILE *
safe_popen(char *buffer, char *mode)
{
  return( safe_open(buffer, mode, 1) );
}


   int
safe_close(FILE *closed, char *buffer) 
{
  stackType *sp, *destroy;
  fileType *file;
  int return_value;
  destroy = NULL;

#ifdef DEBUG_IO
  fprintf(stdout, "safe_close(%x, %s) ", closed, buffer);
#endif

  if ( closed == (FILE *) NULL && buffer == (char *) NULL ) {
    warning("safe_close()", "Attempting to close a NULL FILE and buffer!\n");
    return -1;
  }

  if ( fileStack != NULL ) {

    if ( buffer == NULL ) {
      /* search based on FILE entry */

      destroy = is_open_file_stack(closed);
      
    } else {
      /* search based on buffer entry */

      destroy = is_open_file_buffer_stack(buffer);

    }
  }

  /*
   * if the FILE or buffer is not in the fileStack, close it anyways
   * if it is a file, otherwise simply return...
   */

  if ( destroy == NULL ) {

#ifdef DEBUG_IO
    fprintf(stdout, "entry not found in fileStack\n");
#endif

    if ( closed )
      return ( fclose( closed ) );
    else
      return -1;
  }
    
  /*
   * destroy the current entry...
   */

  file = (fileType *) destroy->entry;
  if ( file->popen ) {
    return_value = pclose(file->file);
  } else 
    return_value = fclose(file->file);

  file->file = NULL;
  safe_free( (void *) file->buffer );
  safe_free( (void *) file->mode );
  safe_free(file);

  if ( destroy == fileStack && destroy == fileStack->next ) {
    /*
     * if there is only one entry, set the fileStack to NULL
     */

    fileStack = NULL;

  } else {
    /*
     * if there are multiple entries, find the one previous
     * to destroy in order to properly set the next entries
     */

    for (sp = fileStack;;) {
      if (sp->next == destroy) break;
      sp = sp->next;
    }
    sp->next = destroy->next;
    if ( destroy == fileStack ) fileStack = destroy->next;

  }  
    
#ifdef DEBUG_IO
  fprintf(stdout,
	  "destroyed (%x), fileStack:\n", destroy);
  printFileStack();
#endif

  safe_free( (void *) destroy );

  return( return_value );
}

/* wrappers to safe_close */

   int
safe_fclose(FILE *closed)
{
  return ( safe_close(closed, (char *) NULL) );
}

   int
safe_fclose_buffer(char *buffer)
{
  return ( safe_close( (FILE *) NULL, buffer) );
}

/*
 * DAN ROE:
 * fileStack_clear()
 * To be called after ptraj finishes. Clean up any entries in the fileStack
 * which have not already been freed.
 */
void fileStack_clear() {
  fileType *file;

  while (fileStack!=NULL) {
    file=(fileType*) fileStack->entry;   
    safe_fclose_buffer(file->buffer);
  }

  return;
}


#define COMMAND_SIZE 100
/* Attempts to open filename with mode, returning the result in fpp.
 * If failure occurs, the prompt is printed to stdout and the user
 * is prompted on stdin to input another filename...
 *
 * WARNING: the return value contains the opened file name and
 * is allocated; it is the callers responsibility to free up the
 * memory...
 *
 */
   char *
promptToOpenFile(FILE **fpp, char *filename, char *mode, char *prompt) 
{
  char buffer[LINE_SIZE];
  char bufferp[LINE_SIZE+COMMAND_SIZE];
  char *filename_ret;
  char *is_compressed;
  int reqprompt = FALSE;

  /* copy the input filename into a buffer
   */
  if ( strcpy(buffer, filename) == NULL ) 
    error("openFile()", "strcpy failure\n");

  /* enter a "return" terminated loop which loops until
   * the file is successfully opened. 
   */
  for (;;) {

    /* print out the prompt if necessary...
     */
    is_compressed = NULL;
    while (strcmp(buffer, "") == 0 || buffer[0] == '\n' || reqprompt ) {
      printf("%s", prompt); fflush(stdout);
      if ( fgets(buffer, LINE_SIZE, stdin) == NULL )
	error("openFile()", "fgets returned NULL\n");
      *strchr(buffer, '\n') = (char) 0;
      reqprompt = FALSE;
    }
    
    /*
     *  Check to see if the files is compressed or supposed to be
     *  compressed by looking for the ".Z" or ".gz" or ".bz2" extension.
     *  If it is present, use popen to open up the file and
     *  the appropriate utility (depending on reading or writing)
     *  which should be in your path or you are hosed...
     */

    if ( (is_compressed = strstr(buffer, ".Z")) != NULL &&
	(is_compressed[2] == (char) 0 || isspace(is_compressed[2])) ) {
 
      /*
       *  .Z found: implied compressed with "compress", read with "zcat"
       */

      if (mode[0] == 'r') 
	sprintf(bufferp, "zcat %s", buffer);
      else if (mode[0] == 'w')
	sprintf(bufferp, "compress > %s", buffer);
      else if (mode[0] == 'a')
	sprintf(bufferp, "compress >> %s", buffer);
      else
	is_compressed = NULL;

      if (is_compressed) {
      
	fprintf(stdout, "Appended .Z detected, using popen(\"%s\")\n", bufferp);
	if ( ( *fpp = safe_popen(bufferp, mode) ) == NULL ) {
	  fprintf(stdout, "\n  Could not open compressed file (%s) with mode (%s), %s\n",
		  buffer, mode, "  try again...\n");
	  reqprompt = TRUE;
	}
      }
    } else if ( (is_compressed = strstr(buffer, ".gz")) != NULL &&
		(is_compressed[3] == (char) 0 || isspace(is_compressed[3])) ) {
 
      if (mode[0] == 'r') 
	sprintf(bufferp, "gunzip -c %s", buffer);
      else if (mode[0] == 'w')
	sprintf(bufferp, "gzip > %s", buffer);
      else if (mode[0] == 'a')
	sprintf(bufferp, "gzip >> %s", buffer);
      else
	is_compressed = NULL;
      fprintf(stdout, "Detected .gz extension, using popen(\"%s\")\n", bufferp);
      if ( ( *fpp = safe_popen(bufferp, mode) ) == NULL ) {
	fprintf(stdout, "\n  Could not open compressed file (%s) with mode (%s), %s\n",
		buffer, mode, "  try again...\n");
	reqprompt = TRUE;
      }
    } else if ( (is_compressed = strstr(buffer, ".bz2")) != NULL &&
                (is_compressed[4] == (char) 0 || isspace(is_compressed[4])) ) {
                                                                                                                                                      
      if (mode[0] == 'r')
        sprintf(bufferp, "bunzip2 -c %s", buffer);
      else if (mode[0] == 'w')
        sprintf(bufferp, "bzip2 > %s", buffer);
      else if (mode[0] == 'a')
        sprintf(bufferp, "bzip2 >> %s", buffer);
      else
        is_compressed = NULL;
      fprintf(stdout, "Detected .bz2 extension, using popen(\"%s\")\n", bufferp);
      if ( ( *fpp = safe_popen(bufferp, mode) ) == NULL ) {
	fprintf(stdout, "\n  Could not open compressed file (%s) with mode (%s), %s\n",
		buffer, mode, "  try again...\n");
	reqprompt = TRUE;
	}
    }


    if (is_compressed == NULL) {

      /*
       * otherwise, try simply opening up the file
       */

      *fpp = NULL;
      if (  ( *fpp = safe_fopen(buffer, mode) ) == NULL ) {
	fprintf(stdout, "\n  Could not open file (%s) with mode (%s), try again...\n",
		buffer, mode);
	reqprompt = TRUE;
      } 
    }

    /*
     * if the file was successfully opened up, then print out
     * summary information otherwise we go back into the loop...
     */
    if ( ! reqprompt ) {
      printf("Opened file \"%s\" with mode (%s)\n", buffer, mode);
      filename_ret = safe_malloc( sizeof( char ) * (strlen(buffer)+1) );
      if ( strcpy(filename_ret, buffer) == NULL )
	error("openFile()", "strcpy failure with %s\n", buffer);
      return( filename_ret );
    }
  }
}

/* gzipFileSize()
 * DRR: Return the uncompressed size in bytes of gzipped file by peeking 
 * at the last 4 bytes.
 * NOTE: long long int should be equivalent to off_t. 
 */
long long int gzipFileSize(char *filename) {
  FILE *infile;
  unsigned char b1,b2,b3,b4;
  long long int val,temp;

  if (filename==NULL) return -1;
  if ( (infile = fopen(filename,"rb"))==NULL ) {
    fprintf(stdout,"Error: gzipFileSize: Could not open %s for reading.\n",filename);
    return -1;
  }

  // Place 4 bytes from the end
  fseek(infile, -4, SEEK_END);

  b1=0; b2=0; b3=0; b4=0;
  fread(&b4,1,1,infile);
  fread(&b3,1,1,infile);
  fread(&b2,1,1,infile);
  fread(&b1,1,1,infile);

  val = 0;
  temp = (long long int) b1;
  temp <<= 24;
  val = val | temp;
  temp = (long long int) b2;
  temp <<= 16;
  val = val | temp;
  temp = (long long int) b3;
  temp <<= 8;
  val = val | temp;
  temp = (long long int) b4;
  val = val | temp;

  //val = (b1 << 24) | (b2 << 16) + (b3 << 8) + b4;

  fclose(infile);

  if (prnlev>0) fprintf(stdout,"gzipFileSize: Uncompressed size of %s: %lli\n",filename,val);

  return val;
}

/* 
 * bzip2FileSize()
 * DRR: Return the uncompressed size of bzip2 file in bytes by counting 
 * all characters using bzcat and wc.
 */
long long int bzip2FileSize(char *filename) {
  long long int val;
  char *command;
  FILE *pipe;

  if (filename==NULL) return -1;
  // Use bzcat <file> | wc -c to calc file size
  command=(char*) malloc( (15 + strlen(filename)) * sizeof(char));
  sprintf(command,"bzcat %s | wc -c",filename);
  if ((pipe=popen(command,"r"))==NULL) {
    fprintf(stdout,"Error: bzip2FileSize: Could not open %s for reading.\n",filename);
    fprintf(stdout,"       Check that bzcat and wc are present on your system.\n");
    return -1;
  }
  fscanf(pipe,"%lli",&val);
  pclose(pipe);

  if (prnlev>0) fprintf(stdout,"bzip2FileSize: Uncompressed size of %s: %lli\n",filename,val);

  return val;
}

/* 
 * zipFileSize()
 * DRR: Return the uncompressed size of zip file in bytes by counting 
 * all characters using unzip and wc.
 */
long long int zipFileSize(char *filename) {
  long long int val;
  char *command;
  FILE *pipe;

  if (filename==NULL) return -1;
  // Use unzip -p <file> | wc -c to calc file size
  command=(char*) malloc( (18 + strlen(filename)) * sizeof(char));
  sprintf(command,"unzip -p %s | wc -c",filename);
  if ((pipe=popen(command,"r"))==NULL) {
    fprintf(stdout,"Error: zipFileSize: Could not open %s for reading.\n",filename);
    fprintf(stdout,"       Check that unzip and wc are present on your system.\n");
    return -1;
  }
  fscanf(pipe,"%lli",&val);
  pclose(pipe);

  if (prnlev>0) fprintf(stdout,"zipFileSize: Uncompressed size of %s: %lli\n",filename,val);

  return val;
}

/* DAN ROE:
 * id_Filesig(): Attempt to identify the file type by first 3 hex vals.
 * A filename or an open stream should be supplied.
 * If a filename is supplied, the file is opened, checked, then closed.
 * If a stream is supplied, the stream is rewound, read, and rewound.
 * Returns a number signifying file type
 *  0: Unknown
 *  1: gzip 
 *  2: bzip2
 *  3: zip
 *  4: netcdf
 * A return value < 0 means an error occured.
 *  -1: File not found
 *  -2: Internal error
 */
int id_Filesig(char *filename, FILE *infile) {
  unsigned char *h;
  int i, type;

  /* Check that either filename or infile is specified, but not both */
  if ((filename==NULL)&&(infile==NULL)) {
    fprintf(stderr,"idFilesig: Both filename and infile are NULL.\n");
    return -2;
  }
  if ((filename!=NULL)&&(infile!=NULL)) {
    fprintf(stderr,"idFilesig: Both filename and infile are not NULL.\n");
    return -2;
  }

  /* Open filename if neccessary, otherwise rewind stream */
  if (infile==NULL) {
    if ( (infile=fopen(filename,"rb"))==NULL) {
      if (prnlev>0) fprintf(stdout,"idFilesig: Cannot open file for signature identification.\n");
      return -1;
    }
  } else
    rewind(infile);

  /* Read first 3 bytes from file */
  h=(unsigned char*) calloc(3,sizeof(unsigned char));
  fread(h,1,1,infile);
  fread(h+1,1,1,infile);
  fread(h+2,1,1,infile);

  if (prnlev>0) fprintf(stderr,"  idFilesig: Hex sig: %x %x %x\n",h[0],h[1],h[2]);

  /* Attempt to identify file  */
  type=0;
  if ((h[0]==0x1f) && (h[1]==0x8b) && (h[2]==0x8)) {
    if (prnlev>0) fprintf(stdout,"  idFilesig: Gzipped file\n"); 
    type=1;
  } else if ((h[0]==0x42) && (h[1]==0x5a) && (h[2]==0x68)) {
    if (prnlev>0) fprintf(stdout,"  idFilesig: Bzip2d file\n"); 
    type=2;
  } else if ((h[0]==0x50) && (h[1]==0x4b) && (h[2]==0x3)) {
    if (prnlev>0) fprintf(stdout,"  idFilesig: Zipped file\n"); 
    type=3;
  } else if ((h[0]==0x43) && (h[1]==0x44) && (h[2]==0x46)) {
    if (prnlev>0) fprintf(stdout,"  idFilesig: NETCDF file\n"); 
    type=4;
  }

  free(h);
  /* If no filename was specified, rewind the stream, otherwise close the file */
  if (filename==NULL)
    rewind(infile);
  else
    fclose(infile);

  return type;
}

/* DAN ROE:
 * id_compression_from_filename():
 * Attempt to identify if user wants compression based on filename
 * extension:
 *  .gz, gzip = 1 
 *  .bz2, bzip = 2
 *  .z, zip(pkzip) = 3
 *  .Z, compress(LZW) = 4
 * This routine is how openFile() used to identify compression.
 */
int id_compression_from_filename(char *filename) {
  char *is_compressed;
  int type;

  is_compressed = NULL;
  type=0;

  if (prnlev>0) fprintf(stdout,"  Attempting to identify compression from filename: %s,",filename); 

  if ( (is_compressed = strstr(filename, ".z")) != NULL &&
       (is_compressed[2] == (char) 0 || isspace(is_compressed[2])) ) {
    type=3;
  } else if ( (is_compressed = strstr(filename, ".gz")) != NULL &&
              (is_compressed[3] == (char) 0 || isspace(is_compressed[3])) ) {
    type=1;
  } else if ( (is_compressed = strstr(filename, ".bz2")) != NULL &&
              (is_compressed[4] == (char) 0 || isspace(is_compressed[4])) ) {
    type=2;
  } else if ( (is_compressed = strstr(filename, ".Z")) != NULL &&
              (is_compressed[2] == (char) 0 || isspace(is_compressed[2])) ) {
    type=4;
  }
 
  if (prnlev>0) fprintf(stdout," detected %i\n",type); 

  return type;
}

/*
 *  Open file specified by filename with "mode" to the file pointer
 *  fpp
 *  Returns 0 if an error occured, otherwise returns 1
 */

   int
openFile(FILE **fpp, char *filename, char *mode) 
{
  char buffer[BUFFER_SIZE];
  int filesig,hasCompression,filename_compresstype;
  char openMode[2];

  /* DAN ROE: openMode is a local copy of the requested file open mode
   * for use with popen. This is done in case popen needs to be called 
   * in append mode. Since popen has no append mode, openMode will be 
   * changed to w. Otherwise openMode will just be mode.  
   */
  strcpy(openMode,mode);

  /*
   *  Check for NULL filename
   */

  if ( strcmp(filename, "") == 0 || filename[0] == '\n' ) {
    *fpp = NULL;
    return 0;
  }

  /*
   *  Check to see if the file exists if we are opening it for reading
   */

  /*
  if ( mode[0] == 'r' ) {
    *fpp = fopen(filename, mode);
    if (*fpp == NULL) {
      return 0;
    } else {
      fclose(*fpp);
    }
  }
  */

  /*
   * DAN ROE:
   * Check to see if the file is compressed (or to be compressed)
   * If so use popen to open up the file via the appropriate command...
   * Currently assumes these commands are in the path.
   * Now using hex signature to ID compression from existing files.
   * When doing output (w|a), if file does not already exist try to ID
   * compression type from filename like before.
   * For backwards compatibility or if the file does not exist look at the
   * filename extension as well (moved to a separate routine.
   */
  hasCompression=0;
  filesig=id_Filesig(filename,NULL);
  /* For backwards compatibility with compress, 
     also look at filename extension  */
  filename_compresstype=id_compression_from_filename(filename);

  if (prnlev>0) fprintf(stdout,"Filesig: %i, Filename_compresstype: %i\n",
                        filesig,filename_compresstype);

  /* Read file - first detect compression via file signature only */
  if (mode[0] == 'r') {
    hasCompression=1;
    switch (filesig) {
      case 1: sprintf(buffer, "gunzip -c %s", filename); break;
      case 2: sprintf(buffer, "bunzip2 -c %s", filename); break;
      case 3: sprintf(buffer, "unzip -p %s", filename); break;
      default: hasCompression=0;
    } 
    /* DAN ROE: If no compression, check for .Z in filename, indicating 
     * UNIX compress. This is for backwards compatibility only since I
     * dont know the signature for this file type.
     */
    if ((hasCompression==0)&&(filename_compresstype==4)) {
      sprintf(buffer, "uncompress -c %s", filename);
      hasCompression=1;
    }

  /* Write/Overwrite file - If file exists, use detected compression (if any).
   * If file does not exist, use compression detected by filename.
   */
  } else if (mode[0] == 'w') {
    /* File not found or detection failed - 
       Use detect compression from filename  */
    if (filesig<=0) 
      filesig=filename_compresstype;
    
    hasCompression=1;
    switch (filesig) {
      case 1: sprintf(buffer, "gzip > %s", filename); break;
      case 2: sprintf(buffer, "bzip2 > %s", filename); break;
      case 3: sprintf(buffer, "zip > %s", filename); break;
      case 4: sprintf(buffer, "compress > %s", filename); break;
      default: hasCompression=0;
    }

  /* Append file - if file exists, use detected compression (if any).
   * If file does not exist, use compression detected by filename.
   */  
  } else if (mode[0] == 'a') {
    /* File not found or detection failed - 
       Attempt to detect compression from filename  */
    if (filesig<=0) 
      filesig=filename_compresstype;

    /* Currently appending to zip files does not work correctly.  */
    if (filesig==3) {
      fprintf(stderr,"openFile: Error: Unable to append to zip (.z) files.\n");
      return 0;
    }

    hasCompression=1;
    switch (filesig) {
      case 1: sprintf(buffer, "gzip >> %s", filename); break;
      case 2: sprintf(buffer, "bzip2 >> %s", filename); break;
      /*case 3: sprintf(buffer, "zip %s -", filename); break;  */
      case 4: sprintf(buffer, "compress >> %s", filename); break;
      default: hasCompression=0;
    }
    /* DAN ROE: If this file is compressed ptraj will use popen. However, popen does
     * not have an "a" mode, just r/w - change mode to w.
     */
    if (hasCompression==1) strcpy(openMode,"w"); 
    
  }
    
  if (hasCompression==1) {

    /* PF - multiptraj - compression not currently supported for parallel ptraj */
#ifdef MPI
    fprintf(stdout, "Compression is not currently supported for parallel ptraj.\n");
    return 0;
#endif

    /*
     *  trajectory is compressed, use popen()
     */

    if (prnlev > 2)
      fprintf(stdout, "Opening compressed file: %s with mode: %s\n", buffer,openMode);

    if ( ( *fpp = safe_popen(buffer, openMode) ) == NULL ) {
      if (prnlev > 4)
	fprintf(stdout, "Could not open compressed file (%s) with mode (%s)\n",
		buffer, openMode);
      return 0;
    }

  } else {

    /*
     *  trajectory is NOT compressed
     */

    if (prnlev > 2)
      fprintf(stdout, "Opening uncompressed file: %s with mode: %s\n",filename,mode);

    *fpp = NULL;
    if (  ( *fpp = safe_fopen(filename, mode) ) == NULL ) {

      if (*fpp == NULL && filename[0] == '~') {

	if (mode[0] == 'r')
	  sprintf(buffer, "cat %s", filename);
	else if (mode[0] == 'w')
	  sprintf(buffer, "cat > %s", filename);
	else if (mode[0] == 'a')
	  sprintf(buffer, "cat >> %s", filename);


	*fpp = safe_popen(buffer, mode);
	if (*fpp == NULL) {
	  if ( prnlev > 4) 
	    fprintf(stdout, "Could not open file (%s) with mode (%s)\n", filename, mode);
	  return 0;
	}
      } else
	return 0;
    }
  }  

  return 1;

}


   int
promptUserResponse(FILE *fpin, FILE *fpout, char *prompt,
		   char *response, int min_to_match)
{
  char buffer[BUFFER_SIZE];
  char *bufferp;

  fprintf(fpout, "%s", prompt);
  fflush(fpout);
  
  if ( fgets(buffer, BUFFER_SIZE, fpin) == NULL )
    error("promptUserResponse()", "fgets returned NULL\n");

  bufferp = buffer;
  skipWhitespace(bufferp);
  bufferp = toLowerCase(bufferp);

  if (strcmp(bufferp, "") == 0 || 
      strncmp(bufferp, response, min_to_match) == 0)
    return 1;
  else
    return 0;
}



   char *
promptUser(FILE *fpin, FILE *fpout, char *prompt)
{
  char buffer[BUFFER_SIZE];
  char *bufferp;

  fprintf(fpout, "%s", prompt);
  fflush(fpout);
  
  if ( fgets(buffer, BUFFER_SIZE, fpin) == NULL )
    error("promptUser()", "fgets returned NULL\n");
  if (strcmp(buffer, "\n") == 0)
    return( (char *) NULL);
  bufferp = (char *) safe_malloc(sizeof(char) * (strlen(buffer)+1) );
  strcpy(bufferp, buffer);
  return(bufferp);
}

void
printfone( char *fmt, ... ) {
  va_list argp;
  va_start(argp, fmt);
#ifdef MPI
  if (worldrank == 0)
    vprintf( fmt, argp );
#else
  vprintf( fmt, argp );
#endif
  va_end(argp);
}
