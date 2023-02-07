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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/trajectory.c,v 10.16 2010/03/17 15:53:13 case Exp $
 *
 *  Revision: $Revision: 10.16 $
 *  Date: $Date: 2010/03/17 15:53:13 $
 *  Last checked in by $Author: case $
 *  ________________________________________________________________________
 *
 *
 *  CONTACT INFO: To learn who the code developers are, who to contact for
 *  more information, and to know what version of the code this is, refer
 *  to the CVS information and the include files (contributors.h && version.h)
 *
 */


 /**** Recent Changes 2009 ****/
 /* Rewritten I/O to support frames - this allows for random seeking in 
    ASCII trajectory files giving exceptional I/O performance improvements
    when seeking in ASCII trajectories - Eric Absgarten (SUNY Stoney Brook)

    1) Significantly faster than old way of doing things
    2) Took Coordinate Checking of *s out of the function
       moving that into checkCoordinates or somewhere else entirely
       since it only really ever needs to be done once on a trajectory.

    Parallel I/O support - multiple MPI threads read different frames. Additionally
    add support for parallel netcdf. - Paul Frybarger (SDSC)
 */

#include "contributors.h"
#include "version.h"

/*  ________________________________________________________________________
 */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#define TRAJECTORY_MODULE
#include "ptraj.h"

/*
 *  This source file contains routines primarily for reading and writing
 *  various coordinate formats for use by rdparm/ptraj.  The routines
 *  (at least for reading/writing trajectories) should be easily
 *  extendable to different formats.
 */

/*----------------------------------------------------------------------

   The following subroutines are defined herein:

   Restart *coords
readAmberRestart(natoms, filename, fp)

   ...reads in an AMBER coordinate or restart file from the open
   file called "filename" (pointed to by the FILE *fp) into a
   Restart *coords structure.  This is returned...

   void
writeAmberRestart(coords, filename, fp)

   ...writes out the AMBER coordinate or restart file specified
   in the coords structure to the open file (FILE *fp) of
   "filename"

   int
getCoordinatesFromRestart(filep, x, y, z, box)

   ...grab out the coordinates from an AMBER restart file from the open file
   (FILE *filep).

   int
loadPdb(fpin, pdb_record **recordP)

   ...load up the PDB into a series of pdb_records

   void
savePdb(fpout, pdb_record *record)

   ...write out a PDB based on a series of pdb_records


   void
savePdbHigh(FILE *fpout, pdb_record *record)

   ...break the pdb format slightly by writing out a high precision version 
   of the temperature and occupancy records.

   int
getCoordinatesFromPdb(pdb_record *pdb, double *x, double *y, double *z)

   ...get the coordinates from a series of pdb_records (loaded up by loadPDB).

   pdb_record *
parmToPdb(Parm *p, int *mask) 

   ...using the atom mask passed in, convert a pdb_record to AMBER names, etc

   void
putCoordinatesInPdb(pdb_record *r, int atoms, double *x, double *y, double *z)

   ...put coordinates into a pdb_record structure

   void
putQandRInPdb(pdb_record *record, int option1, int option2, int debug)

   ...put not only the coordinates but also radii and charges into 
   occupancy and temp factor columns.  See the source below for more info.

   int         [global]
readAmberTrajectory(FILE *fpin, int natoms, 
                    double *x, double *y, double *z, double *box,
		    int set, int hasBoxCoordinates)
   ...reads in a single set from an AMBER-style trajectory dumps.
   The file "fpin" must have previously been opened.  If "set" < 0
   the file is preprocessed to remove the title otherwise the
   trajectory frame is loaded into "x", "y" and "z".
   FIX COMMENTS

   void        [global]
dumpAmberTrajectory(coordinateInfo *outInfo, FILE *fpout, int natoms,
                    double *x, double *y, double *z, double *box);
   ...dumps the coordinates from the arrays x, y and z to the 
   preopened and processed (title previously printed) file (fpout).


   void        [global]
writebinpos( FILE *fpout, int n_atom, double x[], double y[], double z[] )
   ...dumps the coordinates from the arrays x, y and z to the
   preopened and processed (magic numbers already written) file (fpout)
   in the Scripps/Case binary file format.


   int
readbinpos( FILE *fpin, int *n_atom, float apos[], int *eoflag )

   ...read in a BINPOS (Scripps/Case binary file format) file.


------------------------------------------------------------------------*/


/*
 *  routine to read in 4 bytes from the FILE *fd and order in both
 *  big and little endian format...  If the value is "expected",
 *  the byteorder is set (1: little, 0: big).  If the value for both
 *  is not expected, the smaller of the two orders is returned (but a
 *  warning is given).
 */

   int
binaryByteOrder(FILE *fd, int expected, int *byteorder)
{
  int i;
  byte a, b;
  
  i = fgetc(fd);
  ungetc((int) i, fd);

  a.c[0] = fgetc(fd);    /* LOAD UP AS BIG ENDIAN */
  a.c[1] = fgetc(fd);
  a.c[2] = fgetc(fd);
  a.c[3] = fgetc(fd);

  b.c[3] = a.c[0];       /* FLIP TO LITTLE ENDIAN */
  b.c[2] = a.c[1];
  b.c[1] = a.c[2];
  b.c[0] = a.c[3];

  if (prnlev > 2) {
    printf("SGI      (big endian) order value is %d\n", a.i);
    printf("Linux (little endian) order value is %d\n", b.i);
  }

  if (a.i == expected) {
    *byteorder = 0;
    return a.i;
  } else if (b.i == expected) {
    *byteorder = 1;
    return b.i;
  }

  /*
   *  If an unexpected value is found (i.e. not == expected), then don't die,
   *  just return the smaller of the two values and set the endian-ness as 
   *  appropriate...
   */
  fprintf(stderr, "binaryByteOrder() [binary I/O]: integer value read was not the\n");
  fprintf(stderr, "expected value of %i.  Big endian: %i  Little endian: %i\n",
	  expected, a.i, b.i);
  fprintf(stderr, "Returning the smaller value, endian-ness\n");
  
  if (a.i < b.i) {
    *byteorder = 0;
    return a.i;
  } else {
    *byteorder = 1;
    return b.i;
  }

}



   int
readBinaryInteger(FILE *fd, int byteorder)
{
  int i;
  byte u;

  if (byteorder)
    for (i=3; i>=0; i--)
      u.c[i] = fgetc(fd);
  else 
    for (i=0; i<4; i++) 
      u.c[i] = fgetc(fd);

  return(u.i);
}


   void
writeBinaryInteger(FILE *fd, int byteorder, int value)
{
  int i;
  byte u;

  u.i = value;

  if (byteorder)
    for (i=3; i>=0; i--)
      fputc( (int) u.c[i], fd );
  else 
    for (i=0; i<4; i++) 
      fputc( (int) u.c[i], fd );

}


   float
readBinaryFloat(FILE *fd, int byteorder)
{
  int i;
  byte u;

  if (byteorder)
    for (i=3; i>=0; i--)
      u.c[i] = fgetc(fd);
  else
    for (i=0; i<4; i++) 
      u.c[i] = fgetc(fd);

  return(u.f);
}


   void
writeBinaryFloat(FILE *fd, int byteorder, float value)
{
  int i;
  byte u;

  u.f = value;

  if (byteorder)
    for (i=3; i>=0; i--)
      fputc( (int) u.c[i], fd );
  else 
    for (i=0; i<4; i++) 
      fputc( (int) u.c[i], fd );

}


   double
readBinaryDouble(FILE *fd, int byteorder)
{
  int i;
  bytebyte u;

  if (byteorder) 
    for (i=7; i>=0; i--)
      u.c[i] = fgetc(fd);
  else
    for (i=0; i<8; i++)
      u.c[i] = fgetc(fd);

  return(u.d);
}

   void
writeBinaryDouble(FILE *fd, int byteorder, double value)
{
  int i;
  bytebyte u;

  u.d = value;

  if (byteorder) 
    for (i=7; i>=0; i--)
      fputc( (int) u.c[i], fd);
  else
    for (i=0; i<8; i++)
      fputc( (int) u.c[i], fd);
}


   void
printString(void *entry)
{
  char *title;

  title = (char *) entry;
  printf("%s", title);

}



   Restart *
readAmberRestart(int natoms, char *filename, FILE *fp)
{
  char buffer[BUFFER_SIZE];
  char number[13];
  int i, returned;
  Restart *coords;

  coords = (Restart *) safe_malloc( sizeof(Restart) );

  /*
   *  initialize restart structure values             
   */
  coords->isbox = 0;       /* assume no box          */
  coords->restart = 0;     /* assume no velocities   */

  coords->x = (double *) safe_malloc( sizeof(double) * natoms );
  coords->y = (double *) safe_malloc( sizeof(double) * natoms );
  coords->z = (double *) safe_malloc( sizeof(double) * natoms );

     /*
      *  read in the title line 
      */
  if ( fgets(buffer, BUFFER_SIZE, fp) == NULL )
    error("readAmberRestart()()", "fgets returned NULL on title read\n");
  coords->title = (char *) 
    safe_malloc(sizeof(char) * (strlen(buffer) + 1));
  coords->title = strcpy(coords->title, buffer);

     /* 
      *  read in NATOMS and time (into buffer) 
      */
  if ( fgets(buffer, BUFFER_SIZE, fp) == NULL )
    error("translateBox()", "fgets returned NULL\n");

  returned = sscanf(buffer, "%i%lf", 
		    &coords->natoms, 
		    &coords->time);
  
  if ( returned == 1 )
    coords->time = -1.0;   /* no time was present */
  else if ( returned != 2 ) {
    fprintf(stderr,
	    "readAmberRestart(): sscanf on atoms and time in restart file %s failed\n", filename);
    return NULL;
  }
  if ( coords->natoms != natoms ) {
    fprintf(stderr, "readAmberRestart(): topology/coordinates file are inconsist with\n");
    fprintf(stderr, "NATOMS = %i (%i)\n", coords->natoms, natoms);
    return NULL;
  }

     /*
      *  scan in coordinates 
      */
  for (i=0; i < natoms; i += 2) {
    if ( fgets(buffer, BUFFER_SIZE, fp) == NULL )
      error("readAmberRestart()", "fgets returned NULL on reading coords\n");

    if ( strchr(buffer, (int) '*') != NULL ) {
      fprintf(stderr, "readAmberRestart(): WARNING!  File is corrupted...\n");
      fprintf(stderr, "'*' detected in AMBER restrt (%s)\n", filename);
      return NULL;
    }

    number[12] = '\0';
    if ( i == natoms - 1 ) {
       
      strncpy( number, &buffer[0], 12 );
      coords->x[i] = strtod( number, (char **) NULL );
      strncpy( number, &buffer[12], 12 );
      coords->y[i] = strtod( number, (char **) NULL );
      strncpy( number, &buffer[24], 12 );
      coords->z[i] = strtod( number, (char **) NULL );

    } else {

      strncpy( number, &buffer[0], 12 );
      coords->x[i] = strtod( number, (char **) NULL );
      strncpy( number, &buffer[12], 12 );
      coords->y[i] = strtod( number, (char **) NULL );
      strncpy( number, &buffer[24], 12 );
      coords->z[i] = strtod( number, (char **) NULL );
      strncpy( number, &buffer[36], 12 );
      coords->x[i+1] = strtod( number, (char **) NULL );
      strncpy( number, &buffer[48], 12 );
      coords->y[i+1] = strtod( number, (char **) NULL );
      strncpy( number, &buffer[60], 12 );
      coords->z[i+1] = strtod( number, (char **) NULL );
  }
  }

  /*
   *  scan in velocities
   */
  if ( fgets(buffer, BUFFER_SIZE, fp) == NULL ) 
    return(coords);
  else {

    /*
     *  check to see if only the box information is present
     */
    returned = sscanf(buffer, "%lf%lf%lf%lf%lf%lf",
		      &coords->box[0],
		      &coords->box[1],
		      &coords->box[2],
		      &coords->box[3],
		      &coords->box[4],
		      &coords->box[5]);
    if (returned != 3 && returned != 6) {
      fprintf(stderr, "readAmberRestart(): error scanning box from %s\n", buffer);
      return NULL;
    }
    coords->isbox = 1;

    /*
     *  check to see if there is more text, if there is that means that there is
     *  velocity information and therefore get the first velocities from the box
     *  and reset the box information until we find out if it is there...
     */
    if ( fgets(buffer, BUFFER_SIZE, fp) == NULL ) 
      return(coords);

    coords->isbox = 0;

    /*
     *  allocate memory in structures as necessary... 
     */
    coords->vx = (double *) safe_malloc( sizeof(double) * natoms );
    coords->vy = (double *) safe_malloc( sizeof(double) * natoms );
    coords->vz = (double *) safe_malloc( sizeof(double) * natoms );

    coords->vx[0] = coords->box[0];
    coords->vy[0] = coords->box[1];
    coords->vz[0] = coords->box[2];
    coords->vx[1] = coords->box[3];
    coords->vy[1] = coords->box[4];
    coords->vz[1] = coords->box[5];

    for (i=2; i < natoms; i += 2) {
      if ( i > 2 )
	if ( fgets(buffer, BUFFER_SIZE, fp) == NULL )
	  error("readAmberRestart()", 
		"fgets returned NULL on reading velos\n");
      if ( i == natoms - 1 ) {
	if ( sscanf(buffer, "%lf%lf%lf", 
		    &coords->vx[i], 
		    &coords->vy[i], 
		    &coords->vz[i]) != 3 ) 
	  error("readAmberRestart()", "scanning coords from %s\n", buffer);
      } else 
	if ( sscanf(buffer, "%lf%lf%lf%lf%lf%lf", 
		    &coords->vx[i], 
		    &coords->vy[i], 
		    &coords->vz[i],
		    &coords->vx[i+1], 
		    &coords->vy[i+1], 
		    &coords->vz[i+1]) != 6 ) 
	  error("readAmberRestart()", "scanning coords from %s\n", buffer);
    }
    coords->restart = 1;
  }
  
  /* scan in box coordinates */
  if ( fgets(buffer, BUFFER_SIZE, fp) == NULL ) 
    return(coords);
  else {
    returned = sscanf(buffer, "%lf%lf%lf%lf%lf%lf",
		      &coords->box[0],
		      &coords->box[1],
		      &coords->box[2],
		      &coords->box[3],
		      &coords->box[4],
		      &coords->box[5]);
    if (returned != 3 && returned != 6) {
      fprintf(stderr, "readAmberRestart(): error scanning box from %s\n", buffer);
      return NULL;
    }
    coords->isbox = 1;
  }

  return(coords);

}


   void
writeAmberRestart(Restart *coords, char *filename, FILE *fp)
{
  int i;

  /* write out title */
  fprintf(fp, "%s\n", coords->title);
  
  /* write out natoms and time */
  fprintf(fp, "%5i", coords->natoms);
  if ( coords->time > 0 )
    fprintf(fp, "%15.7f\n", coords->time);
  else
    fprintf(fp, "\n");

  /* dump out transformed coordinates */
  for ( i=0; i < coords->natoms; i += 2 ) {
    if ( i == coords->natoms - 1 ) 
      fprintf(fp, "%12.7f%12.7f%12.7f\n", 
	      coords->x[i], 
	      coords->y[i], 
	      coords->z[i]);
    else
      fprintf(fp, "%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f\n", 
	      coords->x[i], 
	      coords->y[i], 
	      coords->z[i], 
	      coords->x[i+1], 
	      coords->y[i+1], 
	      coords->z[i+1]);
  }

  if ( coords->restart ) {
    for ( i=0; i < coords->natoms; i += 2 ) {
      if ( i == coords->natoms - 1 ) 
	fprintf(fp, "%12.7f%12.7f%12.7f\n", 
		coords->vx[i], 
		coords->vy[i], 
		coords->vz[i]);
      else
	fprintf(fp, "%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f\n", 
		coords->vx[i], 
		coords->vy[i], 
		coords->vz[i], 
		coords->vx[i+1], 
		coords->vy[i+1], 
		coords->vz[i+1]);
    }
  }
  
  if ( coords->isbox ) {
    fprintf(fp, "%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f\n", 
	    coords->box[0],
	    coords->box[1],
	    coords->box[2],
	    coords->box[3],
	    coords->box[4],
	    coords->box[5]);
  }
  fflush(fp);
}



   void
dumpAmberRestart(FILE *fpout, int atoms,
		 double *x, double *y, double *z, 
		 double *vx, double *vy, double *vz,
		 double *box, char *title)
{
  int i, vxa, vya, vza;
  Restart *restrt;  

  vxa = 0;
  vya = 0;
  vza = 0;

  if (vx == NULL) {
    vx = safe_malloc(sizeof(double) * atoms);
    for (i=0; i < atoms; i++) {
      vx[i] = 0.0;
    }
    vxa = 1;
  }

  if (vy == NULL) {
    vy = safe_malloc(sizeof(double) * atoms);
    for (i=0; i < atoms; i++) {
      vy[i] = 0.0;
    }
    vya = 1;
  }

  if (vz == NULL) {
    vz = safe_malloc(sizeof(double) * atoms);
    for (i=0; i < atoms; i++) {
      vz[i] = 0.0;
    }
    vza = 1;
  }

  restrt = (Restart *) safe_malloc( sizeof(Restart) );
  if (title == NULL) 
    restrt->title    = "restrt file generated by ptraj";
  else 
    restrt->title = title;
  restrt->natoms   = atoms;
  restrt->time     = 0.0;
  restrt->x        = x;
  restrt->y        = y;
  restrt->z        = z;
  restrt->vx       = vx;
  restrt->vy       = vy;
  restrt->vz       = vz;
  restrt->restart  = 1;
  if (box == NULL) {
    restrt->isbox    = 0;
  } else {
    restrt->isbox    = 1;
    restrt->box[0]   = box[0];
    restrt->box[1]   = box[1];
    restrt->box[2]   = box[2];
    restrt->box[3]   = box[3];
    restrt->box[4]   = box[4];
    restrt->box[5]   = box[5];
  }

  writeAmberRestart(restrt, NULL, fpout);

  safe_free(restrt);
  if (vxa) safe_free(vx);
  if (vya) safe_free(vy);
  if (vza) safe_free(vz);

}


   int
getCoordinatesFromRestart(FILE *fp, double *x, double *y, double *z, double *box)
{
  char buffer[BUFFER_SIZE];
  char number[13];
  int i, returned;
  int natoms;
  double time;

  box[0] = 0.0;
  box[1] = 0.0;
  box[2] = 0.0;
  box[3] = 90.0;
  box[4] = 90.0;
  box[5] = 90.0;

     /*
      *  read in the title line 
      */
  if ( fgets(buffer, BUFFER_SIZE, fp) == NULL )
    error("getCoordinatesFromRestart()()", 
	  "fgets returned NULL on title read\n");

     /*
      *  read in NATOMS and time (into buffer) 
      */
  if ( fgets(buffer, BUFFER_SIZE, fp) == NULL )
    error("getCoordinatesFromRestart()", "fgets returned NULL\n");

  returned = sscanf(buffer, "%i%lf", &natoms, &time);
  
  if ( returned == 1 )
    time = -1.0;   /* no time was present */
  else if ( returned != 2 ) 
    error("getCoordinatesFromRestart()",
	  "sscanf on atoms and time in restart file failed\n");

  if ( natoms != parm->NTOTAT )
    return natoms;

     /*
      *  allocate space for coordinates if necessary 
      */
  if ( x == NULL && y == NULL && z == NULL ) {
    x = (double *) safe_malloc(sizeof(double) * parm->NTOTAT);
    y = (double *) safe_malloc(sizeof(double) * parm->NTOTAT);
    z = (double *) safe_malloc(sizeof(double) * parm->NTOTAT);
  }

     /*
      *  scan in coordinates 
      */
  number[12] = '\0';
  for (i=0; i < natoms; i += 2) {

    if ( fgets(buffer, BUFFER_SIZE, fp) == NULL )
      error("getCoordinatesFromRestart()",
          "fgets returned NULL on reading coords\n");

    if ( strchr(buffer, (int) '*') != NULL ) {
      fprintf(stderr, "getCoordinatesFromRestart(): WARNING!  File is corrupted (contains '*' meaning overflow...)\n");
      return 1;
    }


    if ( i == natoms - 1 ) {
       
      strncpy( number, &buffer[0], 12 );
      x[i] = strtod( number, (char **) NULL );
      strncpy( number, &buffer[12], 12 );
      y[i] = strtod( number, (char **) NULL );
      strncpy( number, &buffer[24], 12 );
      z[i] = strtod( number, (char **) NULL );

    } else {

      strncpy( number, &buffer[0], 12 );
      x[i] = strtod( number, (char **) NULL );
      strncpy( number, &buffer[12], 12 );
      y[i] = strtod( number, (char **) NULL );
      strncpy( number, &buffer[24], 12 );
      z[i] = strtod( number, (char **) NULL );
      strncpy( number, &buffer[36], 12 );
      x[i+1] = strtod( number, (char **) NULL );
      strncpy( number, &buffer[48], 12 );
      y[i+1] = strtod( number, (char **) NULL );
      strncpy( number, &buffer[60], 12 );
      z[i+1] = strtod( number, (char **) NULL );
    }
  }

     /*
      *  scan in velocities 
      */
  if ( fgets(buffer, BUFFER_SIZE, fp) == NULL || buffer[0] == '\n' ) 
    return(natoms);
  else {
    for (i=0; i < natoms; i += 2) {

      if (i == 0) {
	returned = sscanf(buffer, "%lf%lf%lf%lf%lf%lf",
			  &box[0],&box[1],&box[2],&box[3],&box[4],&box[5]);
      } else {

	if ( fgets(buffer, BUFFER_SIZE, fp) == NULL ) {

	  if (i == 2) {
	    /*
	     *  we are bombing because there were no velocities, only box info...
	     */
	    return(natoms);
	  } else {
	    error("getCoordinatesFromRestart()", 
		  "fgets returned NULL on reading velos\n");
	  }
	}
      }
    }
  }

  box[0] = 0.0;
  box[1] = 0.0;
  box[2] = 0.0;
  box[3] = 90.0;
  box[4] = 90.0;
  box[5] = 90.0;
  
     /*
      *  scan in box coordinates
      */
  if ( fgets(buffer, BUFFER_SIZE, fp) == NULL ) 
    return(natoms);
  else {

    returned = sscanf(buffer, "%lf%lf%lf%lf%lf%lf",
		      &box[0],&box[1],&box[2],&box[3],&box[4],&box[5]);
    if (returned != 3 && returned != 6) {
      error("getCoordinatesFromRestart()", 
	    "scanning box info from %s\n", buffer);
    }
  }

  return(natoms);

}



   int
loadPdb(FILE *fpin, pdb_record **recordP)
{
  pdb_record *record;
  int status, line;
  int allocated;

  /* allocate space for the pdb_record's as needed and continuously 
   * call pdb_read_record() to read in the pdb file until EOF is
   * detected, implied by a return value of PDB_END from
   * pdb_read_record()...
   */
  record = (pdb_record *) safe_malloc(sizeof(pdb_record) * BUFFER_SIZE);
  allocated = BUFFER_SIZE;

  line = 0;
  status = 1;
  while ( status ) {
    if (line == allocated) {
      record = (pdb_record *) 
	safe_realloc(record, sizeof(pdb_record)*allocated,
		     sizeof(pdb_record)*BUFFER_SIZE);
      allocated += BUFFER_SIZE;
    }
    record[line] = pdb_read_record(fpin);
    if ( record[line].record_type == PDB_END ) 
      status = 0;
    else
      line++;
  }
  *recordP = record;
  return(line);
}

   void
savePdb(FILE *fpout, pdb_record *record)
{
  int i;
  for (i=0; record[i].record_type != PDB_END; i++)
    if ( record[i].record_type != PDB_UNKNOWN )
      pdb_write_record(fpout, &record[i], NULL, 0);
}


   void
savePdbHigh(FILE *fpout, pdb_record *record)
{
  char *atom_fmt = 
    "ATOM  %5d %-4s%c%-4s%c%4d%c   %8.3f%8.3f%8.3f%8.4f%8.4f\n";
  int i;
  struct pdb_atom *atom;

  for (i=0; record[i].record_type != PDB_END; i++)
    if ( record[i].record_type == PDB_ATOM ||
         record[i].record_type == PDB_HETATM ) {
      
      atom = (struct pdb_atom *) &record[i].pdb;
      fprintf(fpout,
	      atom_fmt,
	      atom->serial_num, 
	      atom->name, 
	      atom->alt_loc,
	      atom->residue.name, 
	      atom->residue.chain_id,
	      atom->residue.seq_num,
	      atom->residue.insert_code,
	      atom->x,
	      atom->y,
	      atom->z,
	      atom->occupancy,
	      atom->temp_factor);

    } else if ( record[i].record_type != PDB_UNKNOWN )
      pdb_write_record(fpout, &record[i], NULL, 0);
}



   int
getCoordinatesFromPdb(pdb_record *pdb, double *x, double *y, double *z)
{
  int entries, i;
  struct pdb_atom *atom;

  /*
   *  initial check to see how many atoms are in the pdb
   */
  for (i = 0, entries = 0; pdb[i].record_type != PDB_END; i++) {
    if (pdb[i].record_type == PDB_ATOM ||
	pdb[i].record_type == PDB_HETATM) {
      entries++;
    }
  }

  /*
   *  allocated space for the coordinates IF NECESSARY
   */
  if ( x == NULL && y == NULL && z == NULL ) {
    x = (double *) safe_malloc(sizeof(double) * entries);
    y = (double *) safe_malloc(sizeof(double) * entries);
    z = (double *) safe_malloc(sizeof(double) * entries);
  }

  /*
   *  load up the coordinates
   */
  for (i = 0, entries = 0; pdb[i].record_type != PDB_END; i++) {
    if (pdb[i].record_type == PDB_ATOM ||
	pdb[i].record_type == PDB_HETATM) {
      atom = (struct pdb_atom *) &pdb[i].pdb;
      x[entries] = (double) atom->x;
      y[entries] = (double) atom->y;
      z[entries] = (double) atom->z;
      entries++;
    }
  }
  return entries;
}

   pdb_record *
parmToPdb(Parm *p, int *mask) 
{
  pdb_record *r;
  int i;
  int atoms;

  if (mask == NULL)
    atoms = p->NTOTAT;
  else {
    atoms = 0;
    for (i=0; i < p->NTOTAT; i++) {
      if (mask[i]) atoms++;
    }
  }

  r = (pdb_record *) safe_malloc(sizeof(pdb_record) * (atoms + 1));

  atoms = 0;
  for (i=0; i < p->NTOTAT; i++) {
    if ( (mask == NULL) || mask[i] ) {
      r[atoms].record_type = PDB_ATOM;
      r[atoms].pdb.atom.serial_num = atoms+1;
      strcpy(r[atoms].pdb.atom.name, p->atom[i].igraph);
      r[atoms].pdb.atom.name[4] = (char) 0;
      r[atoms].pdb.atom.alt_loc = ' ';
      strcpy(r[atoms].pdb.atom.residue.name, 
	     p->residue[p->atom[i].res].labres);
      r[atoms].pdb.atom.residue.name[4] = (char) 0;
      r[atoms].pdb.atom.residue.chain_id = ' ';
      r[atoms].pdb.atom.residue.seq_num = p->atom[i].res+1;
      r[atoms].pdb.atom.residue.insert_code = ' ';
      r[atoms].pdb.atom.x = 0.0;
      r[atoms].pdb.atom.y = 0.0;
      r[atoms].pdb.atom.z = 0.0;
      r[atoms].pdb.atom.occupancy = 0.0;
      r[atoms].pdb.atom.temp_factor = 0.0;
      r[atoms].pdb.atom.ftnote_num = 0.0;
      atoms++;
    }
  }
  r[atoms].record_type = PDB_END;

  return r;
} 


   pdb_record *
ptrajStateToPdb(ptrajState *state, int *mask, int option) 
{
  pdb_record *r;
  int i, mol, molc;
  int atoms, curres;

  if (mask == NULL)
    atoms = state->atoms;
  else {
    atoms = 0;
    for (i=0; i < state->atoms; i++) {
      if (mask[i]) atoms++;
    }
  }

  if (state->IFBOX) {
    atoms += state->molecules;
  }

  r = (pdb_record *) safe_malloc(sizeof(pdb_record) * (atoms + 1));

  atoms = 0;
  for (i=0, mol=0, molc = 0; i < state->atoms; i++, molc++) {

    if (state->IFBOX && molc >= state->moleculeInfo[mol]) {
      mol++;
      r[atoms].record_type = PDB_TER;
      atoms++;
      molc = 0;
    }

    if ( (mask == NULL) || mask[i] ) {
      r[atoms].record_type = PDB_ATOM;

      /*
       *  Handle PDB overflow by restarting atom numbers > 100,000.
       */
      r[atoms].pdb.atom.serial_num = (i+1)%100000;

      if (option) {
          /*
           *  don't wrap the atom names
           */
	strcpy(r[atoms].pdb.atom.name, state->atomName[i]);
      } else {
	r[atoms].pdb.atom.name[1] = state->atomName[i][0];
	r[atoms].pdb.atom.name[2] = state->atomName[i][1];
	r[atoms].pdb.atom.name[3] = state->atomName[i][2];
	r[atoms].pdb.atom.name[0] = state->atomName[i][3];
      }
      r[atoms].pdb.atom.name[4] = (char) 0;
      r[atoms].pdb.atom.alt_loc = ' ';
      curres = atomToResidue(i+1, state->residues, state->ipres)-1;
      strcpy(r[atoms].pdb.atom.residue.name, state->residueName[curres]);
      r[atoms].pdb.atom.residue.name[4] = (char) 0;
      r[atoms].pdb.atom.residue.chain_id = ' ';
      r[atoms].pdb.atom.residue.seq_num = curres+1;
        /*
         *  Handle PDB residue number overflow by adding in chain_id and restarting count
         */
      if ( r[atoms].pdb.atom.residue.seq_num > 9999 ) {
	r[atoms].pdb.atom.residue.seq_num = r[atoms].pdb.atom.residue.seq_num % 10000;
	r[atoms].pdb.atom.residue.chain_id = (char) 
	  ( curres / 10000 + 64); /* 65 is ascii "A" */
      }
      r[atoms].pdb.atom.residue.insert_code = ' ';
      r[atoms].pdb.atom.x = 0.0;
      r[atoms].pdb.atom.y = 0.0;
      r[atoms].pdb.atom.z = 0.0;
      r[atoms].pdb.atom.occupancy = 0.0;
      r[atoms].pdb.atom.temp_factor = 0.0;
      r[atoms].pdb.atom.ftnote_num = 0.0;
      atoms++;
    }
  }
  r[atoms].record_type = PDB_END;

  return r;
} 


   void
putCoordinatesInPdb(pdb_record *r, int atoms, double *x, double *y, double *z)
{
  int i, j;

  j = 0;
  for (i=0; r[i].record_type != PDB_END; i++) {
    if (r[i].record_type == PDB_ATOM || r[i].record_type == PDB_HETATM) {
      r[i].pdb.atom.x = x[j];
      r[i].pdb.atom.y = y[j];
      r[i].pdb.atom.z = z[j];
      j++;
    }
    if (j > atoms) break;
  }
}




   void
putQandRInPdb(pdb_record *record, int option1, int option2, int debug)
{
  int i, j, k;
  struct pdb_atom *atom;
  double rstar, eps, chrg;
  Name atomName;
  double conversion;

  conversion = CHARGE_TO_KCALS;

  for (i=0; record[i].record_type != PDB_END; i++) {

    if (record[i].record_type == PDB_ATOM ||
	record[i].record_type == PDB_HETATM) {

      atom = (struct pdb_atom *) &record[i].pdb;
      strcpy(atomName, atom->name);

      if (option2) {
	/*
	 *  don't wrap atom names
	 */
	while ( atomName[0] == ' ' ) {
	  atomName[0] = atomName[1];
	  atomName[1] = atomName[2];
	  atomName[2] = atomName[3];
	  atomName[3] = atomName[4];
	  atomName[4] = (char) 0;
	}
      } else {

	atomName[4] = atomName[0];
	atomName[0] = atomName[1];
	atomName[1] = atomName[2];
	atomName[2] = atomName[3];
	atomName[3] = atomName[4];
	atomName[4] = (char) 0;
      }

      for ( k=0; atomName[k] != (char) 0; k++ ) ;
      if ( k < 5 )
	for (j=k; j<4; j++)
	  atomName[j] = ' ';

      rstar = 0.0;
      if (option1 == 1) {

	/*
	 * use AMBER radii from prmtop file
	 */

	if ( ! getAtomRadii(atomName, atom->residue.seq_num, &rstar, 0) && debug )
	  warning("", "Atom radii not found for atom %s\n", atomName);

      } else if (option1 == 3) {

	/*
	 * use AMBER vdw radii (r*)
	 */

	if ( ! getAtomTypes(atomName, &rstar, &eps, 0) && debug )
	  warning("", "Atom type (r*) not found for atom %s\n", 
		  atomName);

      } else if (option1 == 2) {

	/*
	 * use PARSE radii
	 */

	switch (atomName[0] ) {

	case 'H':
	  rstar = 1.0;
	  break;
	case 'C':
	  rstar = 1.7;
	  break;
	case 'N':
	  rstar = 1.5;
	  break;
	case 'O':
	  rstar = 1.4;
	  break;
	case 'P':
	  rstar = 2.0;
	  break;
	case 'S':
	  rstar = 1.85;
	  break;
	default:
	  warning("putQandRInPdb()", "Parse radii not present for atom name %s\n",
		  atomName);
	}
      }


      if ( ! getAtomCharge(atomName, atom->residue.seq_num, &chrg, 0) && debug )
	warning("", "Atom charge not found for atom %s\n", 
		atomName);

      if (debug) {
	printf("Atom %4s has charge %8.4f and radii %8.4f\n", 
	       atomName, chrg/conversion, rstar);
      }

      atom->occupancy = chrg / conversion;
      atom->temp_factor = rstar; 

    }
  }
}

readAmberTrajectory(FILE *fpin, int natoms,
		    double *x, double *y, double *z, double *box,
		    int set, coordinateInfo *trajInfo)
{
  char *buffer, *bufferptr, lastchar;
  int atom, coordnum, bi;

  if ( set < 0 ) {
    return 1;
  }

  if ( x == NULL || y == NULL || z == NULL )
    error("readAmberTrajectory()", "coordinate arrays are NULL\n");

  buffer = bufferptr = trajInfo->buffer;
  coordnum = 0;

  /* Read whole set into buffer */
  if ( fread (buffer, trajInfo->frameSize, 1, fpin) == 0)
    return 0;

  for (atom = 0; atom < natoms; atom++) {
    
    /*
     *  Save last character and write null character
     *  Important if the numbers run continuously in
     *  trajectory file, like: 1111.1112222.2223333.333
     */

    lastchar = buffer[8];
    buffer[8] = '\0';
    x[atom] = atof(buffer);
    if (x[atom] == 0.0)   /* atof returns 0.0 if it finds ***'s */
      if (memchr (buffer, '*', 8) != NULL) {   /* check for ***'s in buffer */
	warning("readAmberTrajectory()",
		"Set #%i has coordinate out of bounds (i.e. ******'s)\n", set);
	return 0;
      }
    buffer += 8;
    buffer[0] = lastchar;   /* reset saved character */
    coordnum++;
    if (coordnum % 10 == 0)   /* newline every 10th coord */
      buffer++;

    lastchar = buffer[8];
    buffer[8] = '\0';
    y[atom] = atof(buffer);
    if (y[atom] == 0.0)
      if (memchr (buffer, '*', 8) != NULL) {
	warning("readAmberTrajectory()",
		"Set #%i has coordinate out of bounds (i.e. ******'s)\n", set);
	return 0;
      }
    buffer += 8;
    buffer[0] = lastchar;
    coordnum++;
    if (coordnum % 10 == 0)
      buffer++;

    lastchar = buffer[8];
    buffer[8] = '\0';
    z[atom] = atof(buffer);
    if (z[atom] == 0.0)
      if (memchr (buffer, '*', 8) != NULL) {
	warning("readAmberTrajectory()",
		"Set #%i has coordinate out of bounds (i.e. ******'s)\n", set);
	return 0;
      }
    buffer += 8;
    buffer[0] = lastchar;
    coordnum++;
    if (coordnum % 10 == 0)
      buffer++;

  }

  buffer++;

  /*
   *  Read in box if box coords exist. numbox is the number
   *  of box coordinates, determined by checkCoordinates
   */

  if (trajInfo->numBox) {
    for (bi = 0; bi < trajInfo->numBox; bi++) {
      lastchar = buffer[8];
      buffer[8] = '\0';
      box[bi] = atof(buffer);
      buffer += 8;
      buffer[0] = lastchar;
    }
  }

  /* Reset buffer to beginning of memory */
  buffer = bufferptr;
  return 1;
}

/*  NOTE: this routine is slightly obfuscated in order to handle cases
 *  where trajectory files are corrupted.  Specifically, when running
 *  long time scale (in vacuo) calculations, the molecule can move
 *  far away from the origin, leading to rather large coordinate
 *  values.  Since the trajectory file contains numbers of the format
 *  f8.3 (or %8.3f), this means numbers larger than 9999.999 or smaller
 *  than -999.999 cause a ******** to appear in the output file.  So,
 *  to read in coordinates, a block of 24 characters is read in for
 *  each atom, this block is search for *'s, and then the coordinate
 *  values are extracted into the coordinate arrays.
 */
   int
readAmberTrajectory_nobuffer(FILE *fpin, int natoms, 
                    double *x, double *y, double *z, double *box,
                    int set, int hasBoxCoordinates)
{
  fpos_t fileMarker;         /* marker for current file postition */
  char *junk;
  char buffer[BUFFER_SIZE];
  char coords[26];    /* room for 24 characters (3f8.3) + \n + (char) 0  */
  char c;     
  int j,ret;
  double xx, yy, zz;          /* placeholders for the scanned values */
  int returnValue;
  double j1, j2, j3, j4, j5, j6, j7, j8, j9;

  for (j=0; j<6; j++) {
    box[j] = 0.0;
  }

  if ( set < 0 ) {
    /*
     *  perform initial setup, if not already done... 
     *  (read title)
     */
    if ( fgets(buffer, BUFFER_SIZE, fpin) == NULL )
      error("readAmberTrajectory()", "fgets returned NULL\n");
    return 1;
  }

  if ( x == NULL || y == NULL || z == NULL )
    error("readAmberTrajectory()",
          "coordinate arrays are NULL\n");

  for (j=0; j < natoms; j++) {
    /*
     *  scan in 24 characters representing the 3f8.3 values 
     */
    if ( returnValue = fscanf(fpin, "%24c", coords) ) {
      if ( returnValue != 1 ) {
        if (j != 0) {
          junk = (char *) coords;
          fprintf(stderr, "\nSet #%i appears corrupted (%s)\n", set, junk);
        }
        return 0;
      }
    }
    coords[24] = (char) 0;
    /* 
     *  if a newline is found, this suggests that one character too
     *  few has been read in; read in this extra character...
     */
    if ( strchr(coords, '\n') != NULL ) {
      fscanf(fpin, "%1c", &c);
      coords[24] = c;
      coords[25] = (char) 0;
    }


    if (prnlev > 7) {
      fprintf(stderr, "Set #%5i -- coordinates: %s\n", set, coords);
    }
    
    /* 
     *  check for corrupted trajectory 
     */
    if ( strchr(coords, '*') != NULL ) {
      warning("readAmberTrajectory()", 
              "Set #%i has coordinate out of bounds (i.e. *****'s)\n", set);
      return 0;
    }
    
    /* 
     *  place into coordinate arrays... 
     */
    ret = 0;
    if (sscanf(coords, "%lf", &xx) == 1) {
      ret = 1;
      junk = strchr(coords, '.');
      if (junk != NULL) junk += 4;
      if (junk != NULL && sscanf(junk, "%lf", &yy) == 1) {
        ret = 2;
        junk = strchr(junk, '.'); 
        if (junk != NULL) junk += 4;
        if (junk != NULL && sscanf(junk, "%lf", &zz) == 1) {
          ret = 3;
        }
      }
    }


    if (prnlev > 7) {
      fprintf(stderr, "                            %7.3f %7.3f %7.3f\n", xx, yy, zz);
    }


    /*
    ret = sscanf(coords, "%f%f%f", &xx, &yy, &zz);
    */
    if (ret != 3 ) {
      junk = (char *) coords;
      warning("readAmberTrajectory()",
              "Set #%i is corrupted (%s)...\n", set, junk);
      if (prnlev > 2) {
        fprintf(stdout, "Only %i values returned on scanning %s, expecting 3\n",
                ret, coords);
      }
      return 0;
    }
    x[j] = xx;
    y[j] = yy;
    z[j] = zz;
  }

  /* search for box coordinates,
   * note: we need > 2 atoms for this check to work, then we mark the current
   * file position, get a line of text.  If the line of text is not NULL
   * this implies that we are not at the end of the file and we can then
   * try to scan for three atoms worth of coordinates.  If this fails, this
   * implies we have box coordinates (up to six with ewald!) and we
   * do not need to reset the file position.  Otherwise we reset the
   * file position.  What a nightmare!  Doh!!!! This doesn't work with
   * popen()'ed files since fseek() doesn't work on these.  Therefore, we
   * add logic to allow the caller to specify that yes, this file DOES
   * have box coordinates to avoid the check...
   */

  if ( hasBoxCoordinates != 0 && natoms > 2 ) {
    if (hasBoxCoordinates == 1 || fgetpos(fpin, &fileMarker) == 0) {
      junk = fgets(buffer, 80, fpin);
      if (junk != NULL && (strcmp(buffer, "\n") == 0)) {
        junk = fgets(buffer, 80, fpin);
      }

      if (junk != NULL) {

        j = sscanf(buffer, "%lf%lf%lf%lf%lf%lf%lf%lf%lf", 
                   &j1, &j2, &j3, &j4, &j5, &j6, &j7, &j8, &j9);
        if (j == 3) {
          /* have box coordinates */
          box[0] = j1; 
          box[1] = j2; 
          box[2] = j3;
        } else if (j == 6) {
          /* ewald? */
          box[0] = j1; 
          box[1] = j2; 
          box[2] = j3;
          box[3] = j4; 
          box[4] = j5; 
          box[5] = j6;
        } else if ( hasBoxCoordinates == 1 || fsetpos(fpin, &fileMarker) != 0 ) {
            warning("readAmberTrajectory()", 
                    "fsetpos() failed on when looking for box coordinates!\n");
            return 0;
        }
      }
    } else {
      perror("fgetpos() failed, error string follows:\n");
    }
  }
  return 1;
}


/* the following are used by dump trajectory */

void write_one( FILE *fpout, double x, double y, double z, int* pcount )
{
    fprintf(fpout, "%8.3lf", x);
    (*pcount)++;
    if ( *pcount % 10 == 0 ) 
        fprintf(fpout, "\n");

    fprintf(fpout, "%8.3lf", y);
    (*pcount)++;
    if ( *pcount % 10 == 0 )
        fprintf(fpout, "\n");
    
    fprintf(fpout, "%8.3lf", z);
    (*pcount)++;
    if ( *pcount % 10 == 0 )
        fprintf(fpout, "\n");
}


void
dumpAmberTrajectory(coordinateInfo *outInfo, int natoms,
		    double *x, double *y, double *z, double *box)
{
  int count, j;

  static char formatStandard[] = "%8.3lf";
  static char formatTruncated[] = "%8.2lf";
  static char formatWayTruncated[] = "%8.1lf";
  char *format;
  char *buffer, *bufferptr;

  buffer = bufferptr = outInfo->buffer;

  format = (char *) formatStandard;
  for (j=0; j < natoms; j++) {
    if (x[j] > 9999.999 || x[j] < -999.999 ||
	y[j] > 9999.999 || y[j] < -999.999 ||
	z[j] > 9999.999 || z[j] < -999.999) {
      format = (char *) formatTruncated;
      for (count = j; count < natoms; count++) {
	if (x[count] > 99999.999 || x[count] < -9999.999 ||
	    y[count] > 99999.999 || y[count] < -9999.999 ||
	    z[count] > 99999.999 || z[count] < -9999.999) {
	format = (char *) formatWayTruncated;
	count = natoms;
	}
      }
      j = natoms;
    }
  }

  for (j=0, count = 0; j < natoms; j++) {
    sprintf(buffer, format, x[j]);
    buffer += 8;
    count++; 
    if ( (count % 10) == 0 ) {
      sprintf(buffer, "\n");
      buffer++;
    }

    sprintf(buffer, format, y[j]);
    buffer += 8;
    count++;
    if ( (count % 10) == 0 ) {
      sprintf(buffer, "\n");
      buffer++;
    }
    
    sprintf(buffer, format, z[j]);
    buffer += 8;
    count++;
    if ( (count % 10) == 0 ) {
      sprintf(buffer, "\n");
      buffer++;
    }

  }
  if ( count % 10 ) {
      sprintf(buffer, "\n");
      buffer++;
  }
  if (box != NULL) {
    sprintf(buffer, "%8.3lf%8.3lf%8.3lf\n", 
	    box[0], box[1], box[2]);
    buffer += 25;
    buffer[0] = '\0';
  }
  
#ifdef MPI
  if (outInfo->isMPI) {
    /* now print to file, in an ordered fashion so sets are ordered  */
    MPI_File_write_ordered(*(MPI_File *) outInfo->file, bufferptr, strlen(bufferptr), MPI_CHAR, MPI_STATUS_IGNORE);
  } else {
#endif
    fwrite(bufferptr, 1, strlen(bufferptr), outInfo->file);
#ifdef MPI
  }
#endif

  /* reset the buffer to point to beginning of memory  */
  buffer = bufferptr;
}

void dumpAmberBox(FILE* fpout, double* box)
{
  if (box != NULL) 
  {
    fprintf(fpout, "%8.3lf%8.3lf%8.3lf\n", box[0], box[1], box[2]);
  }
  fflush(fpout);
}

   int
openbinpos(FILE *fpin)
{
	char	magic[ 10 ];

	if( fread( magic, 1, 4, fpin ) != 4 ){
		fprintf( stderr, "Couldn't read magic number from BINPOS\n" );
		return( -1 );
	}

	magic[ 4 ] = '\0';
	if( strcmp( magic, "fxyz" ) != 0 ){
		fprintf( stderr, "bad magic number \"%s\"\n", magic );
		return( -1 );
	}
	return 0;
}

   int
readbinpos( FILE *fpin, int *n_atom, float apos[], int *eoflag )
{
	int	count, n_atom_in;
	*eoflag = 0 ;

	if( fread( &n_atom_in, sizeof( int ), 1, fpin ) != 1 ) {
		*eoflag = 1;
		return 0; 
	}
	*n_atom = n_atom_in;
/*	fprintf( stderr, "looking for %d atoms\n", n_atom_in);  */
	if( ( count = fread( apos, sizeof( float ), 3 * *n_atom, fpin ) )
		!= 3 * *n_atom ){
		fprintf( stderr, "Could only read %d of %d atoms requested\n",
			count / 3, *n_atom );
		return( -1 );
	}
	return(1);
}



   void
writebinpos( FILE *fpout, int n_atom, double x[], double y[], double z[] )
{
	int i,j;
	float *apos;

	fwrite( &n_atom, sizeof( int ), 1, fpout ) ;

	apos = (float *) safe_malloc(sizeof(float) * 3 * n_atom);
	j = 0;
	for( i=0; i<n_atom; i++ ){
		apos[j] = x[i];
		apos[j+1] = y[i];
		apos[j+2] = z[i];
		j += 3;
	}

	fwrite( apos, sizeof( float ), 3 * n_atom, fpout );

	safe_free(apos);
}


/*
 *  Convert between the symmetric shape matrix (xtl) and the unit cell 
 *  parameters (box).  If forward > 0 then convert xtl->box otherwise
 *  convert box->xtl
 *
 *  NOTE: this leads to values for the symmetric shape matrix that are
 *  *slightly* different than what CHARMM produces (in the 14th decimal
 *  place).  Therefore, any trajectories converted with the XTL information
 *  will be slightly different!!!
 */

   void
xtlabcToBox(double *xtl, double *box, int forward)
{
  double a, b, c, ab, bc, ca;
  double hth[3][3];
  double evalue[3], ev[3][3];

  if (forward > 0) {

    /*
      A  = SQRT(XTLABC(1)**2 + XTLABC(2)**2 + XTLABC(4)**2)
      B  = SQRT(XTLABC(2)**2 + XTLABC(3)**2 + XTLABC(5)**2)
      C  = SQRT(XTLABC(4)**2 + XTLABC(5)**2 + XTLABC(6)**2)
      AB = XTLABC(2)*(XTLABC(1) + XTLABC(3)) + XTLABC(4)*XTLABC(5)
      BC = XTLABC(5)*(XTLABC(3) + XTLABC(6)) + XTLABC(2)*XTLABC(4)
      CA = XTLABC(4)*(XTLABC(1) + XTLABC(6)) + XTLABC(2)*XTLABC(5)
      XUCELL(1) = A
      XUCELL(2) = B
      XUCELL(3) = C
      XUCELL(4) = ACOS(BC/(B*C))*RADDEG
      XUCELL(5) = ACOS(CA/(C*A))*RADDEG
      XUCELL(6) = ACOS(AB/(A*B))*RADDEG
     */

    box[0] = sqrt( xtl[0]*xtl[0] + xtl[1]*xtl[1] + xtl[3]*xtl[3] );
    box[1] = sqrt( xtl[1]*xtl[1] + xtl[2]*xtl[2] + xtl[4]*xtl[4] );
    box[2] = sqrt( xtl[3]*xtl[3] + xtl[4]*xtl[4] + xtl[5]*xtl[5] );

    ab = xtl[1]*(xtl[0] + xtl[2]) + xtl[3]*xtl[4];
    bc = xtl[4]*(xtl[2] + xtl[5]) + xtl[1]*xtl[3];
    ca = xtl[3]*(xtl[0] + xtl[5]) + xtl[1]*xtl[4];

    box[3] = acos( bc / (box[1]*box[2]) ) * RADDEG;
    box[4] = acos( ca / (box[2]*box[0]) ) * RADDEG;
    box[5] = acos( ab / (box[0]*box[1]) ) * RADDEG;

#ifdef DEBUG_XTLABC
      printf("\n\nBOX   : %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", 
	     box[0], box[1], box[2], box[3], box[4], box[5]);
      printf("XTLABC: %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", 
	     xtl[0], xtl[1], xtl[2], xtl[3], xtl[4], xtl[5]);
#endif


  } else {

    /*
    hth[0][0] = box[0]*box[0];
    
    if (box[5] - 90.0 > 0.0000000001)
      hth[1][0] = box[0]*box[1]*cos(DEGRAD*box[5]);
    else
      hth[1][0] = 0.0;
    hth[0][1] = hth[1][0];

    hth[1][1] = box[1]*box[1];

    if (box[4] - 90.0 > 0.0000000001)
      hth[2][0] = box[0]*box[2]*cos(DEGRAD*box[4]);
    else
      hth[2][0] = 0.0;
    hth[0][2] = hth[2][0];

    if (box[3] - 90.0 > 0.0000000001)
      hth[1][2] = box[1]*box[2]*cos(DEGRAD*box[3]);
    else
      hth[1][2] = 0.0;
    hth[2][1] = 0.0;

    hth[2][2] = box[2]*box[2];
    */


    hth[0][0] = box[0]*box[0];
    
    a = box[5]-90.0;
    if ( ABS(a) > 0.000000001)
      hth[1][0] = box[0]*box[1]*cos(DEGRAD*box[5]);
    else 
      hth[1][0] = 0.0;
    hth[0][1] = hth[1][0];

    hth[1][1] = box[1]*box[1];

    a = box[4]-90.0;
    if ( ABS(a) > 0.000000001)
      hth[2][0] = box[0]*box[2]*cos(DEGRAD*box[4]);
    else
      hth[2][0] = 0.0;
    hth[0][2] = hth[2][0];
 
    a = box[3]-90.0;
    if ( ABS(a) > 0.000000001)
      hth[1][2] = box[1]*box[2]*cos(DEGRAD*box[3]);
    else
      hth[1][2] = 0.0;
    hth[2][1] = 0.0;

    hth[2][2] = box[2]*box[2];

    jacobi(hth, 3, evalue, ev);

    if (evalue[0] < 0.000000001 ||
	evalue[1] < 0.000000001 ||
	evalue[2] < 0.000000001) {
      fprintf(stderr, "xtlabcToBox, error on conversion from box->xtlabc!\n");
    }

    a = sqrt(evalue[0]);
    b = sqrt(evalue[1]);
    c = sqrt(evalue[2]);

    xtl[0] = a*ev[0][0]*ev[0][0] + b*ev[0][1]*ev[0][1] + c*ev[0][2]*ev[0][2];
    xtl[2] = a*ev[1][0]*ev[1][0] + b*ev[1][1]*ev[1][1] + c*ev[1][2]*ev[1][2];
    xtl[5] = a*ev[2][0]*ev[2][1] + b*ev[2][1]*ev[2][1] + c*ev[2][2]*ev[2][2];
    xtl[1] = a*ev[0][0]*ev[1][0] + b*ev[0][1]*ev[1][1] + c*ev[0][2]*ev[1][2];
    xtl[3] = a*ev[0][0]*ev[2][0] + b*ev[0][1]*ev[2][1] + c*ev[0][2]*ev[2][2];
    xtl[4] = a*ev[1][0]*ev[2][0] + b*ev[1][1]*ev[2][1] + c*ev[1][2]*ev[2][2];

#ifdef DEBUG_XTLABC
    printf("\nXTLABC: %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", 
	   xtl[0], xtl[1], xtl[2], xtl[3], xtl[4], xtl[5]);
    printf("BOX   : %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n\n", 
	   box[0], box[1], box[2], box[3], box[4], box[5]);
#endif

  }

}



   int
readCharmmTrajectory(FILE *fpin, charmmTrajectoryInfo **trajp, 
		     double *x, double *y, double *z,
		     double *box, int set)
{
  int i, j, start, stop, loadfree, toload;
  int hdr;
  charmmTrajectoryInfo *traj;
  char *title;

  if (set >= 0 && *trajp == NULL) {
     /*
      *  make sure that the charmmTrajectoryInfo has been previously 
      *  setup; this is done by calling this routine with set == -1 
      */
    error("readCharmmTrajectory()", "charmmTrajectoryInfo is NULL");
  }

  traj = *trajp;

  if (set <= 0) {
       /*
        *  read in the header and icntrl information; setup the 
        *  charmmTrajectoryInfo structure
        */ 
    if (traj != NULL) {
      if (traj->titleStack != NULL) {
	clearStack(&traj->titleStack);
	traj->titleStack = NULL;
      }
      if (traj->nfreat != traj->natrec && traj->freeat != NULL) {
	safe_free(traj->freeat);
	traj->freeat = NULL;
	safe_free(traj->fixedx);
	traj->fixedx = NULL;
	safe_free(traj->fixedy);
	traj->fixedy = NULL;
	safe_free(traj->fixedz);
	traj->fixedz = NULL;
      }
      safe_free(traj);
      traj=NULL;
    }
    traj = (charmmTrajectoryInfo *) safe_malloc(sizeof(charmmTrajectoryInfo));
    INITIALIZE_charmmTrajectoryInfo(traj);
    *trajp = traj;

    for (i=0; i < 20; i++) traj->icntrl[i] = 0;

       /*
        * skip initial 4-bytes (FORTRAN record); this represents how many bytes 
	* are to be read for this record.  The value is written both before and
        * after the record that is to be read.  If it is -1, then we are at the end
        */
    hdr = binaryByteOrder(fpin, 84, &traj->byteorder);

       /*
        * load the magic string and check it out
        */
    traj->magic.i = readBinaryInteger(fpin, 0);
    if (traj->magic.c[0] != 'C' &&
	traj->magic.c[1] != 'O' &&
	traj->magic.c[2] != 'R' &&
	traj->magic.c[3] != 'D') {
      fprintf(stderr, 
	      "Attempting to read a CHARMM trajectory file: incorrect header! %c%c%c%c\n",
	      traj->magic.c[0], traj->magic.c[1], traj->magic.c[2], traj->magic.c[3]);
      return -1;
    }

       /*
        * load ICNTRL variables and print for debugging
        */
    for (i=0; i < 20; i++) traj->icntrl[i] = readBinaryInteger(fpin, traj->byteorder);

    if (traj->icntrl[11] == 1) {
      warning("readCharmmTrajectory", "DIM4 (fourth dimension) not supported yet!\n");
      return -1;
    }

    hdr =  readBinaryInteger(fpin, traj->byteorder);

       /*
        * read in TITLE
        */

    hdr = readBinaryInteger(fpin, traj->byteorder);
    traj->ntitle = readBinaryInteger(fpin, traj->byteorder);

    if (prnlev > 2) 
      printf("NTITLE = %6i\n", traj->ntitle);

    for (j=0; j < traj->ntitle; j++) {

      title = (char *) safe_malloc(sizeof(char) * 82);
      for (i=0; i<80; i++) {
	title[i] = fgetc(fpin);
      }
      title[80] = '\n';
      title[81] = (char) 0;
      pushBottomStack(&traj->titleStack, (void *) title);

    }

    if (prnlev > 2) {
      fprintf(stdout, "Dumping the titles loaded...\n");
      printStack(&traj->titleStack, printString, NULL);
    }

    hdr = readBinaryInteger(fpin, traj->byteorder);


       /*
        *  read in NATREC (the number of atom records written), set NFREAT
        */
    hdr = readBinaryInteger(fpin, traj->byteorder);
    traj->natrec = readBinaryInteger(fpin, traj->byteorder);
    if (prnlev > 2) 
      printf("NATREC = %6i\n", traj->natrec);

    hdr = readBinaryInteger(fpin, traj->byteorder);

    traj->nfreat = traj->natrec - traj->icntrl[8];
    if (prnlev > 2)
      printf("NFREAT = %6i\n", traj->nfreat);

       /*
        *  if NATREC != NFREAT read the list of free atoms
        */
    if (traj->nfreat != traj->natrec) {
      hdr = readBinaryInteger(fpin, traj->byteorder);
      traj->freeat = (int *) safe_malloc(sizeof(int) * traj->nfreat);
      for (i=0; i<traj->nfreat; i++) 
	traj->freeat[i] = readBinaryInteger(fpin, traj->byteorder);

      if (prnlev > 4) {
	printf("Dumping free atoms\n");
	for (i=0; i<traj->nfreat; i++) {
	  printf("%6i ", traj->freeat[i]);
	  if ( (i+1)%10 == 0) printf("\n");
	}
	if ( (i+1)%10 ) printf("\n");
      }

      hdr = readBinaryInteger(fpin, traj->byteorder);

    }

    /*
     *  Return from set < 0 (preprocessing)
     */
    return 1;
  }

  /*
   *  LOAD UP THE COORDINATES
   */

     /*
      * load crystal/box information if necessary
      */

  if (traj->icntrl[10] == 1) {

    hdr = readBinaryInteger(fpin, traj->byteorder);
    if (hdr == -1) return 0;
    else {

      for (i=0; i<6; i++) traj->xtlabc[i] = readBinaryDouble(fpin, traj->byteorder);

      if (traj->icntrl[19] >= 22)
	xtlabcToBox(traj->xtlabc, box, 1);

      hdr = readBinaryInteger(fpin, traj->byteorder);
      if (hdr == -1) return 0;
    }
  }

     /*
      *  load up the coordinates
      */

  loadfree = (traj->nfreat != traj->natrec && set != 0);
  toload = traj->natrec;
  if (loadfree) toload = traj->nfreat;

  hdr = readBinaryInteger(fpin, traj->byteorder);
  if (hdr == -1) return 0;
  for (i=0; i < toload; i++)
    x[ (loadfree ? traj->freeat[i]-1 : i) ] = (double) readBinaryFloat(fpin, traj->byteorder);
  hdr = readBinaryInteger(fpin, traj->byteorder);
  if (hdr == -1) return 0;

  hdr = readBinaryInteger(fpin, traj->byteorder);
  if (hdr == -1) return 0;
  for (i=0; i < toload; i++)
    y[ (loadfree ? traj->freeat[i]-1 : i) ] = (double) readBinaryFloat(fpin, traj->byteorder);
  hdr = readBinaryInteger(fpin, traj->byteorder);
  if (hdr == -1) return 0;

  hdr = readBinaryInteger(fpin, traj->byteorder);
  if (hdr == -1) return 0;
  for (i=0; i < toload; i++)
    z[ (loadfree ? traj->freeat[i]-1 : i) ] = (double) readBinaryFloat(fpin, traj->byteorder);
  hdr = readBinaryInteger(fpin, traj->byteorder);
  if (hdr == -1) return 0;

     /*
      *  if this is the first visit, all of the coordinates will be read;
      *  the charmmTrajectoryInfo fixed? arrays need to be updated
      */
  if (set == 1 && traj->nfreat != traj->natrec) {
    traj->fixedx = (double *) safe_malloc(sizeof(double) * (traj->natrec - traj->nfreat));
    traj->fixedy = (double *) safe_malloc(sizeof(double) * (traj->natrec - traj->nfreat));
    traj->fixedz = (double *) safe_malloc(sizeof(double) * (traj->natrec - traj->nfreat));

    start = 0;
    toload = 0;

    for (i=0; i < traj->nfreat+1; i++) {
      if (i == traj->nfreat)
	stop = traj->natrec+1;
      else
	stop = traj->freeat[i];

      for (j=start; j < stop-1; j++) {
	traj->fixedx[toload] = x[j];
	traj->fixedy[toload] = y[j];
	traj->fixedz[toload] = z[j];
	toload++;
      }
      start = stop;
    }
  }

     /*
      *  set coordinates for the fixed atoms if necessary
      */

  if (loadfree) {
    start = 0;
    toload = 0;
    for (i=0; i < traj->nfreat+1; i++) {
      if (i == traj->nfreat)
	stop = traj->natrec+1;
      else
	stop = traj->freeat[i];

      for (j=start; j < stop-1; j++) {
	x[j] = traj->fixedx[toload];
	y[j] = traj->fixedy[toload];
	z[j] = traj->fixedz[toload];
	toload++;
      }
      start = stop;
    }
  }

#ifdef DEBUG_CHARMM
  for (i=0; i < 5; i++)
    printf("ATOM %6i  OH  WAT %5i     %7.3f %7.3f %7.3f\n", i+1, i+1, 
	   x[i], y[i], z[i]);
  printf("=====\n");
#endif
  return 1;
  
}


   void
dumpCharmmTrajectory(FILE *fpout, charmmTrajectoryInfo *traj, int atoms,
		     double *x, double *y, double *z,
		     double *box, int set)
{

  int i, dumpfree, todump;
  char *title;
  stackType *t;
  /*
  char *ptrajHeader = "*  MODIFIED/WRITTEN BY PTRAJ";
  */
  if (traj == NULL) return;


  /*
   *  NOTE: all written "records" must be preceeded and followed by writing the
   *  number of bytes to be written...
   */


  if ( set < 0 ) {
    /*
     *  DUMP OUT THE HEADERS AND INFO WHICH PRECEEDS THE COORDINATES
     */

       /*
        *  If there is a mismatch between the expected number of
        *  atom records and what we plan on writting (likely due to
        *  a strip or closestwaters command), change the number of atom
        *  records to write and also do *not* treat any atoms as fixed 
        *  (since it will be a pain to trim the list of fixed atoms)
        */
    if (traj->natrec != atoms) {
      traj->natrec = atoms;
      traj->nfreat = atoms;
      traj->icntrl[8] = 0;
    }


       /*
        *  Record 1: the magic header and icntrl variables
        */
    writeBinaryInteger(fpout, traj->byteorder, 84);
    writeBinaryInteger(fpout, 0, traj->magic.i);
    for (i=0; i < 20; i++) 
      writeBinaryInteger(fpout, traj->byteorder, traj->icntrl[i]);
    writeBinaryInteger(fpout, traj->byteorder, 84);

       /*
        *  Record 2: the titles!
        */

    i = 0;
    for (t = traj->titleStack; t != NULL; t = t->next) {
      title = (char *) t->entry;
      i += strlen(title)-1;
      title[strlen(title)-1] = (char) 0;
    }
    /*
    i += strlen(ptrajHeader)-1;
    */
    writeBinaryInteger(fpout, traj->byteorder, i+4);
    writeBinaryInteger(fpout, traj->byteorder, traj->ntitle);
    for (t = traj->titleStack; t != NULL; t = t->next) {
      title = (char *) t->entry;
      fputs(title, fpout);
    }
    /*
    title = (char *) safe_malloc(sizeof(char) * (strlen(ptrajHeader)+1));
    strcpy(title, ptrajHeader);
    pushBottomStack(&traj->titleStack, (void *) title);
    fputs(title, fpout);
    */
    writeBinaryInteger(fpout, traj->byteorder, i+4);

       /*
        *  Record 3: NATREC
        */
    writeBinaryInteger(fpout, traj->byteorder, 4);
    writeBinaryInteger(fpout, traj->byteorder, traj->natrec);
    writeBinaryInteger(fpout, traj->byteorder, 4);

       /*
        *  Record 4: IF FIXED ATOMS
        */
    if (traj->natrec != traj->nfreat) {

      writeBinaryInteger(fpout, traj->byteorder, traj->nfreat*4);
      for (i=0; i < traj->nfreat; i++) 
	writeBinaryInteger(fpout, traj->byteorder, traj->freeat[i]);
      writeBinaryInteger(fpout, traj->byteorder, traj->nfreat*4);
    }

    /*
     *  END OF PRE-PROCESSING/HEADER INFORMATION
     */ 

  } else {

    if (traj->icntrl[10]) {
      /*
       *  write XTLABC information
       */

      writeBinaryInteger(fpout, traj->byteorder, 48);
      if (traj->icntrl[19] >= 22)
	xtlabcToBox(traj->xtlabc, box, 0);
      for (i=0; i<6; i++) 
	writeBinaryDouble(fpout, traj->byteorder, traj->xtlabc[i]);
      writeBinaryInteger(fpout, traj->byteorder, 48);
    }

    dumpfree = (traj->nfreat != traj->natrec && set != 0);
    todump = traj->natrec;
    if (dumpfree) todump = traj->nfreat;

    /* X */
    writeBinaryInteger(fpout, traj->byteorder, todump*4);
    for (i=0; i < todump; i++)
      writeBinaryFloat(fpout, traj->byteorder, (float) x[(dumpfree ? traj->freeat[i]-1 : i)]);
    writeBinaryInteger(fpout, traj->byteorder, todump*4);

    /* Y */
    writeBinaryInteger(fpout, traj->byteorder, todump*4);
    for (i=0; i < todump; i++)
      writeBinaryFloat(fpout, traj->byteorder, (float) y[(dumpfree ? traj->freeat[i]-1 : i)]);
    writeBinaryInteger(fpout, traj->byteorder, todump*4);

    /* Z */
    writeBinaryInteger(fpout, traj->byteorder, todump*4);
    for (i=0; i < todump; i++)
      writeBinaryFloat(fpout, traj->byteorder, (float) z[(dumpfree ? traj->freeat[i]-1 : i)]);
    writeBinaryInteger(fpout, traj->byteorder, todump*4);

  }
}

/* ============================== TRAJ ROUTINES =============================== 
 * DAN ROE: Idea is to abstract the handling of trajectories in oder to
 * simplify ptraj.c and make future modifications (adding traj types)
 * easier. These routines operate on the coordinateInfo structure,
 * which should contain enough information to do what is necessary.
 */

/* 
 * openTraj()
 * Given a coordinateInfo structure, open the trajectory.
 * Return 0 on ok, 1 on error.
 */
int openTraj(coordinateInfo *trajInfo) {
  int err;

  if (prnlev>0) fprintf(stdout,"openTraj(): Opening %s\n",trajInfo->filename);

  if (trajInfo->isNetcdf==1)
    return NETCDF_open(trajInfo);
  else if (trajInfo->isMPI==1) {
    switch (trajInfo->accessMode) {
      case 0: err=parallel_open_file_read(trajInfo, trajInfo->filename);  break;
      case 1: err=parallel_open_file_write(trajInfo, trajInfo->filename); break;
      case 2: 
        printfone("WARNING: openTraj(): ptraj.MPI does not support opening files in append mode.\n");
        return 1;
      break;
    }
  } else {
    switch (trajInfo->accessMode) {
      case 0: err=openFile(&(trajInfo->file), trajInfo->filename, "r"); break;
      case 1: err=openFile(&(trajInfo->file), trajInfo->filename, "w"); break;
      case 2: err=openFile(&(trajInfo->file), trajInfo->filename, "a"); break;
    }
    // openFile return of 0 means failure
    if (err==0) return 1; 
  }
  return 0;
}

/*
 * closeTraj()
 * Given a coordinateInfo structure, close the trajectory.
 * All error checking will be done by the various subroutines.
 * Return 0 on OK, 1 on error
 */
int closeTraj(coordinateInfo *trajInfo) {
  int err,i;

  // Close any REMD trajectories
  for (i=1; i<trajInfo->numREMDTRAJ; i++)
    closeTraj(trajInfo->REMDtraj[i]);

  if (prnlev>0) fprintf(stdout,"closeTraj(): Closing %s\n",trajInfo->filename);

  if (trajInfo->isNetcdf==1)
    return NETCDF_close(trajInfo);
  else if (trajInfo->isMPI==1)
    return parallel_close_file(trajInfo);
  else { 
    // On error, safe_fclose returns EOF for normal files, -1 for pipes.
    if (trajInfo->file==NULL) return 0;
    err=safe_fclose(trajInfo->file);
    if (err==EOF || err==-1) 
      return 1;
    else
      return 0;
  }
  
  // Should never get here
  return 1;
} 

/*
 * cleanTraj()
 * Clean up a coordinateInfo structure. Close the associated file if 
 * necessary and free all memory.
 */
void cleanTraj(coordinateInfo *f) {
  int i;

  if (prnlev>0) fprintf(stdout,"cleanTraj(): Cleaning %s\n",f->filename);
  // Close file just in case 
  //closeTraj(f);

  safe_free(f->filename);
  safe_free(f->buffer);
  safe_free(f->info);
  if (f->NCInfo != NULL) {
      safe_free(f->NCInfo->timeUnits);
      safe_free(f->NCInfo->coordinateUnits);
      safe_free(f->NCInfo->cellLengthUnits);
      safe_free(f->NCInfo->cellAngleUnits);
      safe_free(f->NCInfo->velocityUnits);
      safe_free(f->NCInfo->Conventions);
      safe_free(f->NCInfo->ConventionVersion);
      safe_free(f->NCInfo->R);
      safe_free(f->NCInfo); 
  }
  safe_free(f->mask);
  safe_free(f->x);
  safe_free(f->y);
  safe_free(f->z);

  safe_free(f->title);
  safe_free(f->program);
  safe_free(f->application);
  safe_free(f->version);

  /* REMDTRAJ Cleanup */
  safe_free(f->baseFilename);
  safe_free(f->compressEXT);
  if (f->REMDtraj!=NULL) { 
    // Start at 1; 0 points to this info structure
    for (i=1; i < f->numREMDTRAJ; i++)
      cleanTraj(f->REMDtraj[i]);
    free(f->REMDtraj);
  }

  INITIALIZE_coordinateInfo(f);
  safe_free(f);
  return;
}

/*
 * trajFile_fgets()
 * fgets wrapper. Makes no sense and should not be called for netcdf 
 */
char *trajFile_fgets(char *buffer, int num, coordinateInfo *C) {

  if (C==NULL) return NULL;

  if (C->isMPI==0) 
    return fgets(buffer,num,C->file);
  else 
    return parallel_fgets(buffer,num,C);

  return 0;
}

/* trajFile_fseek()
 * Wrapper for fseek. Seek to a given frame. 
 * fseeko is used for better compatibility with large files. To avoid
 * losing bits each variable is explicitly converted to off_t in offset 
 * calculation.
 */
int trajFile_fseek(coordinateInfo *C, int frame) {
  int err;
  off_t offset;

  if (C==NULL) return -1;

  if (C->isMPI==0) {
    offset = (off_t) frame;
    offset *= (off_t) C->frameSize;
    offset += (off_t) C->titleSize;
    err=fseeko(C->file, offset, SEEK_SET);
  } else
    err=parallel_fseek(C,frame);

  return err;
}

/* trajFile_rewind()
 */
int trajFile_rewind(coordinateInfo *C) {
  int err;

  if (C==NULL) return 1;

  if (C->compressType!=0) {
    // seek doesnt work on compressed files. Reopen
    err=closeTraj(C);
    if (err!=0) {
      fprintf(stdout,"Error: trajFile_rewind(): Could not properly close file.\n");
      return 1;
    }
    err=openTraj(C);
  } else if (C->isMPI==0)
    err=fseek(C->file, 0L, SEEK_SET);
  else
    err=parallel_rewind(C);

  return err;
}

/* trajFile_fseek_end()
 * Wrapper for fseek, Seek to end of file.
 */
int trajFile_fseek_end(coordinateInfo *C) {
  int err;

  if (C==NULL) return 1;

  if (C->isMPI==0) 
    err=fseek(C->file,0,SEEK_END);
  else
    err=parallel_fseek_end(C);

  if (err!=0) {
    perror("ERROR: trajFile_fseek_end:");
    return 1;
  }
  return 0;
}

/* trajFile_get_position()
 * Wrapper for get position/ftell 
 * NOTE: Careful about long conversions. MPI_Offset is long long I think 
 */
int trajFile_get_position(coordinateInfo *C, long int *offset) {

  if (C==NULL) return -1;

  if (C->isMPI==0)
    *offset=ftell(C->file);
  else
    parallel_get_position(C,offset);

  if (*offset<0) {
    fprintf(stdout,"ERROR: trajFile_get_position\n");
    return 1;
  }   

  return 0;
}

/* trajFile_fread()    
 * Wrapper for fread. Put a frame into the frame buffer.
 * To be consistent with older read routines, Return 0 on error, 1 on success. 
 */
int trajFile_fread(coordinateInfo *C) {

  if (C->isMPI==0) {
    return fread(C->buffer,C->frameSize,1,C->file);
    /*if (ferror(C->file)) {
      perror("ERROR: trajFile_fread:");
      return 0;
    }*/
  } else
    return parallel_fread(C);

  return 1;
}
    
/*
 * bufferToXYZ()
 * Take a char buffer of frameSize bytes and convert it to 3 arrays of
 * doubles X, Y and Z. Assumes numbers are 8 bytes long and arranged
 * X0 Y0 Z0 X1 Y1 Z1 etc.
 * If there are box coordinates read those into box.
 */
int bufferToXYZ(char *buffer, int frameSize, double *X, double *Y, double *Z, int numBox, double box[6]) {
  char *ptr,*end;
  char number[9];
  int i,j,atom,boxsize;
  
  number[8]='\0';
  boxsize=(numBox*8)+1; // may want to precalc
  end=buffer+frameSize-boxsize;
  i=0;    // XYZ index
  j=0;    // number index
  atom=0; // atom index
  for (ptr=buffer; ptr<end; ptr++) {
    if (*ptr=='\n') continue;
    if (*ptr=='*') return -1;
    number[j++]=*ptr;
    i++;
    if (i==8) {
      X[atom]=atof(number);
      j=0;
    } else if (i==16) {
      Y[atom]=atof(number);
      j=0;
    } else if (i==24) {
      Z[atom++]=atof(number);
      j=0;
      i=0;
    }
  }
  // Box stuff goes here
  if (numBox>0) {
    number[0]=ptr[8];
    ptr[8]='\0';
    box[0]=atof(ptr);
    ptr[8]=number[0];

    number[0]=ptr[16];
    ptr[16]='\0';
    box[1]=atof(ptr+8);
    ptr[16]=number[0];

    number[0]=ptr[24];
    ptr[24]='\0';
    box[2]=atof(ptr+16);
    ptr[24]=number[0];

    if (numBox>3) {
      ptr+=24;
      number[0]=ptr[8];
      ptr[8]='\0';
      box[0]=atof(ptr);
      ptr[8]=number[0];

      number[0]=ptr[16];
      ptr[16]='\0';
      box[1]=atof(ptr+8);
      ptr[16]=number[0];

      number[0]=ptr[24];
      ptr[24]='\0';
      box[2]=atof(ptr+16);
      ptr[24]=number[0];
    }
  }

  /*box[0]=0.0; 
  box[1]=0.0;
  box[2]=0.0;
  box[3]=0.0;
  box[4]=0.0;
  box[5]=0.0;*/

  return 0;
}



/*
 * getTrajTemp()
 * If processing as REMD files, search for the target temperature.
 * Otherwise get the temperature of the default file.
 */
/*int trajFile_getTemp(coordinateInfo *trajInfo, double *repTemp) {
  int err;

  err=0;
  switch (trajInfo->type) {
    case COORD_AMBER_REMD:
     if (fscanf(currentRep,"%*s %*s %*s %*s %lf",repTemp)==EOF)
       err=1;
     break;

    case COORD_AMBER_NETCDF:
      err=NETCDF_get_temperature(trajInfo,repTemp);
      break; 
 
    default:
      printfone("ERROR: Temperature read requested, but this trajectory type does not have temperature!\n");
      err=1;
  }
  return err;
}*/

/*
 * getReplicaTrajTempID()
 * If this is a replica trajInfo (i.e. the trajInfo structure contains
 * multiple trajectories from a single REMD run), scan through each
 * to find the target temperature (remdtrajtemp) and return the trajectory ID. 
 * Return -1 on error.
 */
/*int getReplicaTrajTempID(coordinateInfo *trajInfo) {


}*/
