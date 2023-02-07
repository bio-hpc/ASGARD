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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/torsion.c,v 10.0 2008/04/15 23:24:11 case Exp $
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
#include <string.h>
#include <math.h>

#define TORSION_MODULE
#include "ptraj.h"



/*  Given four connected points, A--B--C--D, a dihedral for
 *  the B--C bond can be constructed from...
 *
 *             -1
 *  angle = cos   ABxBC * CDx(-BC)
 *
 *  NOTE: the cross products should be normalized?
 */

   double
torsion(double x1, double y1, double z1,
	double x2, double y2, double z2,
	double x3, double y3, double z3,
	double x4, double y4, double z4)
{
  double Lx, Ly, Lz, Lnorm;
  double Rx, Ry, Rz, Rnorm;
  double Sx, Sy, Sz;
  double angle;

  VOP_3D_COORDS_CROSS_PRODUCT(     Lx,      Ly,      Lz,
			      (x2-x1), (y2-y1), (z2-z1),
			      (x3-x2), (y3-y2), (z3-z2));

  VOP_3D_COORDS_CROSS_PRODUCT(     Rx,      Ry,      Rz,
			      (x4-x3), (y4-y3), (z4-z3),
			      (x2-x3), (y2-y3), (z2-z3));

  Lnorm = sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
  Rnorm = sqrt(Rx*Rx + Ry*Ry + Rz*Rz);

  VOP_3D_COORDS_CROSS_PRODUCT(Sx, Sy, Sz,
			      Lx, Ly, Lz,
			      Rx, Ry, Rz);

  angle = (Lx*Rx + Ly*Ry + Lz*Rz) / (Lnorm * Rnorm);

  if ( angle > 1.0 ) angle = 1.0;
  if ( angle < -1.0 ) angle = -1.0;

  angle = acos( angle );
  angle = angle * RADDEG;

  if ( (Sx * (x3-x2) + Sy * (y3-y2) + Sz * (z3-z2)) < 0 )
    angle = -angle;
  return(angle);
}




   double
angle(double x1, double y1, double z1,
      double x2, double y2, double z2,
      double x3, double y3, double z3)
{
  double angle, xij, yij, zij, xkj, ykj, zkj, rij, rkj;

  xij = x1 - x2;
  yij = y1 - y2;
  zij = z1 - z2;

  xkj = x3 - x2;
  ykj = y3 - y2;
  zkj = z3 - z2;

  rij = xij*xij + yij*yij + zij*zij;
  rkj = xkj*xkj + ykj*ykj + zkj*zkj;

  if (rij > SMALL && rkj > SMALL) {
    angle = (xij*xkj + yij*ykj + zij*zkj) / sqrt(rij*rkj);
    if (angle > 1.0) 
      angle = 1.0;
    else if (angle < -1.0)
      angle = -1.0;
    angle = acos(angle) * RADDEG;
  } else
    angle = 0.0;

  return(angle);
}




   double
puckeras(double x1, double y1, double z1,
       double x2, double y2, double z2,
       double x3, double y3, double z3,
       double x4, double y4, double z4,
       double x5, double y5, double z5, 
       double *amplitude)
{
  double pucker;
  double v1, v2, v3, v4, v5, a, b;

  pucker = 0.0;
  *amplitude = 0.0;


  v4 = torsion(x4,y4,z4,x5,y5,z5,x1,y1,z1,x2,y2,z2);
  v5 = torsion(x5,y5,z5,x1,y1,z1,x2,y2,z2,x3,y3,z3);
  v1 = torsion(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4);
  v2 = torsion(x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5);
  v3 = torsion(x3,y3,z3,x4,y4,z4,x5,y5,z5,x1,y1,z1);

  a = (v1*cos(0.0) + 
       v2*cos( 4.0*PI/5.0) +
       v3*cos( 8.0*PI/5.0) +
       v4*cos(12.0*PI/5.0) +
       v5*cos(16.0*PI/5.0))*0.4;

  b = (v1*sin(0.0) + 
       v2*sin( 4.0*PI/5.0) +
       v3*sin( 8.0*PI/5.0) +
       v4*sin(12.0*PI/5.0) +
       v5*sin(16.0*PI/5.0))*-0.4;

  *amplitude = sqrt(a*a + b*b);

  if (*amplitude != 0.0)
    pucker = atan2(b,a)*RADDEG;
  if (pucker < 0) pucker += 360.0;

  return(pucker);

}

   double
puckercp(double x1i, double y1i, double z1i,
	 double x2i, double y2i, double z2i,
	 double x3i, double y3i, double z3i,
	 double x4i, double y4i, double z4i,
	 double x5i, double y5i, double z5i, 
       double *amplitude)
{
  double pucker, norm;
  double x1, y1, z1;
  double x2, y2, z2;
  double x3, y3, z3;
  double x4, y4, z4;
  double x5, y5, z5;
  double nx, ny, nz;
  double rcx, rcy, rcz;
  double r1x, r1y, r1z;
  double r2x, r2y, r2z;
  double sum1, sum2;

  pucker = 0.0;
  *amplitude = 0.0;

  x2 = x1i; y2 = y1i; z2 = z1i;
  x3 = x2i; y3 = y2i; z3 = z2i;
  x4 = x3i; y4 = y3i; z4 = z3i;
  x5 = x4i; y5 = y4i; z5 = z4i;
  x1 = x5i; y1 = y5i; z1 = z5i;
 
     /*
      *  calculate geometric center
      */
  rcx = (x1 + x2 + x3 + x4 + x5)/5.0;
  rcy = (y1 + y2 + y3 + y4 + y5)/5.0;
  rcz = (z1 + z2 + z3 + z4 + z5)/5.0;

  x1 -= rcx; y1 -= rcy; z1 -=rcz;
  x2 -= rcx; y2 -= rcy; z2 -=rcz;
  x3 -= rcx; y3 -= rcy; z3 -=rcz;
  x4 -= rcx; y4 -= rcy; z4 -=rcz;
  x5 -= rcx; y5 -= rcy; z5 -=rcz;

     /*
      *  calculate normal vectors
      */
  r1x = x1 * sin(0.0) +
        x2 * sin(2.0*PI/5.0) +
        x3 * sin(4.0*PI/5.0) +
        x4 * sin(6.0*PI/5.0) +
        x5 * sin(8.0*PI/5.0);
  r1y = y1 * sin(0.0) +
        y2 * sin(2.0*PI/5.0) +
        y3 * sin(4.0*PI/5.0) +
        y4 * sin(6.0*PI/5.0) +
        y5 * sin(8.0*PI/5.0);
  r1z = z1 * sin(0.0) +
        z2 * sin(2.0*PI/5.0) +
        z3 * sin(4.0*PI/5.0) +
        z4 * sin(6.0*PI/5.0) +
        z5 * sin(8.0*PI/5.0);

  r2x = x1 * cos(0.0) +
        x2 * cos(2.0*PI/5.0) +
        x3 * cos(4.0*PI/5.0) +
        x4 * cos(6.0*PI/5.0) +
        x5 * cos(8.0*PI/5.0);
  r2y = y1 * cos(0.0) +
        y2 * cos(2.0*PI/5.0) +
        y3 * cos(4.0*PI/5.0) +
        y4 * cos(6.0*PI/5.0) +
        y5 * cos(8.0*PI/5.0);
  r2z = z1 * cos(0.0) +
        z2 * cos(2.0*PI/5.0) +
        z3 * cos(4.0*PI/5.0) +
        z4 * cos(6.0*PI/5.0) +
        z5 * cos(8.0*PI/5.0);

  /*
   *  calculate vector normal to plane
   */



  VOP_3D_COORDS_CROSS_PRODUCT( nx,  ny,  nz,
			       r1x, r1y, r1z,
			       r2x, r2y, r2z );

  norm = sqrt(nx*nx + ny*ny + nz*nz);
  nx /= norm;
  ny /= norm;
  nz /= norm;

  /*
   *  rotate around Z axis
   */
  z1 = x1*nx + y1*ny + z1*nz;
  z2 = x2*nx + y2*ny + z2*nz;
  z3 = x3*nx + y3*ny + z3*nz;
  z4 = x4*nx + y4*ny + z4*nz;
  z5 = x5*nx + y5*ny + z5*nz;


  sum1 = z1 * cos(0.0) +
         z2 * cos( 4.0*PI/5.0) +
         z3 * cos( 8.0*PI/5.0) +
         z4 * cos(12.0*PI/5.0) +
         z5 * cos(16.0*PI/5.0);
  sum2 = -(z1 * sin(0.0) +
           z2 * sin( 4.0*PI/5.0) +
           z3 * sin( 8.0*PI/5.0) +
           z4 * sin(12.0*PI/5.0) +
           z5 * sin(16.0*PI/5.0));

  norm = sqrt(sum1*sum1 + sum2*sum2);
  *amplitude = norm * sqrt(2.0/5.0);
  pucker = asin( sum2 / norm );
  if (sum1 < 0.0)
    pucker = PI - pucker;
  else if (pucker < 0.0)
    pucker += 2.0*PI;

  pucker *= RADDEG;

  return(pucker);

}

