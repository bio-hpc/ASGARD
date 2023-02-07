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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/display.c,v 10.0 2008/04/15 23:24:11 case Exp $
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


#include <stdio.h>
#include <math.h>
#include <string.h>

#define DISPLAY_MODULE
#include "ptraj.h"

/*
 *  Code to dump out graphs/plots, etc should be placed here
 */


/* the box size and margins in inches */
#define BOX_SIZE 6.0
#define X_MARGIN 1.25
#define Y_MARGIN 1.25

/*
 * dumpTopo(): this routine will dump a postscript file which represents a 2D
 * grid with the third dimension represented by various gray scale values over
 * the range.
 *
 * input parameters:
 *
 *  fileOut    -- the FILE * of the already opened file
 *  maxValues  -- the size of 1 dimension of the square array in values
 *  values     -- values[i+maxValues*j] is the height at point i,j
 *  entries    -- the actual number of set values (i.e. max i)
 *  shades     -- the number of shades to use (sort of...)
 *     = 0     -- implies continuous shading over entire range
 *     = n     -- implies "n" shading levels
 *  shadeCut   -- IF NULL, then "n" snading levels are used over the range.
 *                IF not NULL, then shadeCut(i)->shadeCut(i+1) is the
 *                interval for the shade.  As i gets larger, the color gets
 *                darker...
 *
 */

   void 
dumpTopo(FILE *fileOut, int max_values, float *values, 
	 int entries, int shades, float *shadeCut)
{
  int i, j, k;
  float box_size;
  float x, y, z;
  float max_entry, min_entry;
  float range;
  float value;

  max_entry = values[0];
  min_entry = values[0];
  for (i=0; i < entries; i++) {
    for (j=0; j < entries; j++ ) {
      if (max_entry < values[i+max_values*j])
	max_entry = values[i+max_values*j];
      if (min_entry > values[i+max_values*j])
	min_entry = values[i+max_values*j];
    }
  }
  range = max_entry - min_entry;

  box_size = BOX_SIZE / (float) entries;

  /* write out the head are start defining procedures */
  fprintf(fileOut, "%%!PS\n %% -- define procedures --\n");
  fprintf(fileOut, " /inch {72 mul} def\n");

  /* define a box */
  fprintf(fileOut, " /box\n");
  fprintf(fileOut, " { newpath moveto\n");
  fprintf(fileOut, " %8.4f inch %8.4f inch rlineto\n", box_size, 0.0);
  fprintf(fileOut, " %8.4f inch %8.4f inch rlineto\n", 0.0, box_size);
  fprintf(fileOut, " %8.4f inch %8.4f inch rlineto\n", -box_size, 0.0);
  fprintf(fileOut, " closepath } def\n");

  /* define how to fill a box with gray */
  fprintf(fileOut, " /fillbox { setgray fill } def\n");

  for (i=0; i < entries; i++) {
    for (j=0; j < entries; j++ ) {

      x = (float) i / (float) entries * BOX_SIZE + X_MARGIN;
      y = (float) j / (float) entries * BOX_SIZE + Y_MARGIN;

      if (shades == 0 ) {

	z = (max_entry - values[i+max_values*j]) / range;

      } else if ( shadeCut == NULL) {

	z = (max_entry - values[i+max_values*j]) / range;
	z = z * shades;
	z = (int) z / (float) shades;
	  
      } else {

	z = -1.0;
	value = values[i+max_values*j];
	for (k = 0; k < 2*shades; k += 2) {
	  if ( value >= shadeCut[k] && value <= shadeCut[k+1] )
	    z = (float) k / (float) (2.0 * shades);
	}
	if (z == -1.0) z = 1.0;
      }
	
      fprintf(fileOut, " %8.4f inch %8.4f inch box\n", x, y);
      fprintf(fileOut, " %8.4f fillbox\n", z);

    }
  }

  fprintf(fileOut, " showpage\n");

# ifdef DISPLAY_MODULE
  safe_fclose(fileOut);
# else
  fclose(fileOut);
#endif
}


