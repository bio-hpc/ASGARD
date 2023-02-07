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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/vector.h,v 10.0 2008/04/15 23:24:11 case Exp $
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




       /* T = U cross V 
        */
#define VOP_3D_COORDS_CROSS_PRODUCT(TX,TY,TZ,UX,UY,UZ,VX,VY,VZ) \
  TX = (UY * VZ) - (UZ * VY); \
  TY = (UZ * VX) - (UX * VZ); \
  TZ = (UX * VY) - (UY * VX)

       /* Multiply the 3x3 matrix times the coordinates specified
        * in x, y and z.  xx, yy and zz are temporary variables.
        * The x, y and z arrays are modified!
        */
#define VOP_3x3_TIMES_COORDS(matrix, x, y, z, xx, yy, zz) \
  xx = matrix[0][0] * x + \
       matrix[0][1] * y + \
       matrix[0][2] * z;  \
  yy = matrix[1][0] * x + \
       matrix[1][1] * y + \
       matrix[1][2] * z;  \
  zz = matrix[2][0] * x + \
       matrix[2][1] * y + \
       matrix[2][2] * z;  \
  x = xx; y = yy; z = zz


       /* Multiply the transpose of the 3x3 matrix times the 
        * coordinates specified in x, y and z.  xx, yy and zz 
        * are temporary variables. 
        * The x, y and z arrays are modified!
        */
#define VOP_3x3_TRANSPOSE_TIMES_COORDS(matrix, x, y, z, xx, yy, zz) \
  xx = matrix[0][0] * x + matrix[1][0] * y + matrix[2][0] * z;  \
  yy = matrix[0][1] * x + matrix[1][1] * y + matrix[2][1] * z;  \
  zz = matrix[0][2] * x + matrix[1][2] * y + matrix[2][2] * z;  \
  x = xx; y = yy; z = zz


#define VOP_3x3_TIMES_3x3(t, u, v) \
  t[0][0] = u[0][0] * v[0][0] + u[0][1] * v[1][0] + u[0][2] * v[2][0]; \
  t[0][1] = u[0][0] * v[0][1] + u[0][1] * v[1][1] + u[0][2] * v[2][1]; \
  t[0][2] = u[0][0] * v[0][2] + u[0][1] * v[1][2] + u[0][2] * v[2][2]; \
  t[1][0] = u[1][0] * v[0][0] + u[1][1] * v[1][0] + u[1][2] * v[2][0]; \
  t[1][1] = u[1][0] * v[0][1] + u[1][1] * v[1][1] + u[1][2] * v[2][1]; \
  t[1][2] = u[1][0] * v[0][2] + u[1][1] * v[1][2] + u[1][2] * v[2][2]; \
  t[2][0] = u[2][0] * v[0][0] + u[2][1] * v[1][0] + u[2][2] * v[2][0]; \
  t[2][1] = u[2][0] * v[0][1] + u[2][1] * v[1][1] + u[2][2] * v[2][1]; \
  t[2][2] = u[2][0] * v[0][2] + u[2][1] * v[1][2] + u[2][2] * v[2][2]
 
