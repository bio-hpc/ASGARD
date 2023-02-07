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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/rms.c,v 10.1 2008/07/08 23:20:51 sbrozell Exp $
 *
 *  Revision: $Revision: 10.1 $
 *  Date: $Date: 2008/07/08 23:20:51 $
 *  Last checked in by $Author: sbrozell $
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
#include <math.h>
#include "utility.h"
#include "vector.h"

/*
 *  Define the RMS routines; other ptraj functionality isn't necessary!
 */


   void
normalize( double a[3] )
{
  double b;

  b = 1.0/sqrt((double)(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]));
  
  a[0] *= b;
  a[1] *= b;
  a[2] *= b;
}

   int
diagEsort(double *mat, double *Emat, double *Evec[], double *Eigenvalue)
{
  int njrot;
  int i, j, k, i3;
  double eigenvector[9], *eA, v;

  if (!jacobi3(mat, Eigenvalue, eigenvector, &njrot)) {
    printf("convergence failed\n");
    return(0);
  }

  for (i=i3=0; i<3; i++, i3+=3)
    for (j=0; j<3; j++)
      Emat[i3+j] = eigenvector[j*3+i];

  for (i=0; i<3; i++)
    Evec[i] = (double *) &Emat[i*3];

  for (i=0; i<2; i++) {
    v = Eigenvalue[k=i];
    for (j=i+1; j<3; j++)
      if (Eigenvalue[j] > v)  
        v = Eigenvalue[k=j];
    if (k != i) {

      Eigenvalue[k] = Eigenvalue[i];
      Eigenvalue[i] = v;
      eA = Evec[i];
      Evec[i] = Evec[k];
      Evec[k] = eA;
    }
  }
  return(1);
}

/*
 *  jacobi3() - get jacobian of 3x3 matrix
 */

#define ROTATE(ARR,MAJ1,MIN1,MAJ2,MIN2) { \
  g = ARR[MAJ1 + MIN1]; \
  h = ARR[MAJ2 + MIN2]; \
  ARR[MAJ1 + MIN1] = g - s*(h+g*tau); \
  ARR[MAJ2 + MIN2] = h + s*(g-h*tau); }

   int
jacobi3(double *a, double *d, double *v, int *nrot)
{
  int  i, j, ip, iq, p3, j3;
  double  tresh, theta, tau, t, sm, s, h, g, c, b[3], z[3];

  for (ip=p3=0; ip<3; ip++,p3+=3) {
    /*
     *  initialize the identity matrix 
     */
    for (iq=0; iq<3; iq++) 
      v[p3 + iq] = 0.0;
    v[p3 + ip] = 1.0;
    /* 
     *  initialize b and d to diagonal of a
     */
    b[ip] = d[ip] = a[p3 + ip];
    z[ip] = 0.0;
  }
  *nrot = 0;
  for (i=0; i<50; i++) {

    sm = 0.0;
    for (ip=p3=0; ip<2; ip++,p3+=3) {
      for (iq=ip+1; iq<3; iq++)
        sm += fabs(a[p3 + iq]);
    }

    if (sm == 0.0) {
      return(1);
    }
    if (i < 3) 
      tresh = sm * 0.2 / 9.0;   /* on 1st three sweeps... */
    else       
      tresh = 0.0;      /* thereafter... */
    for (ip=p3=0; ip<2; ip++,p3+=3) {
      for (iq=ip+1; iq<3; iq++) {
        g = 100.0 * fabs(a[p3 + iq]);

        if ( i > 3  &&  fabs(d[ip])+g == fabs(d[ip])&& fabs(d[iq])+g == fabs(d[iq])) {
          a[p3 + iq] = 0.0;
        } else if (fabs(a[p3 + iq]) > tresh) {
          h = d[iq]-d[ip];
          if (fabs(h)+g==fabs(h))
            t = a[p3 + iq] / h;
          else {
            theta = 0.5 * h / a[p3 + iq];
            t = 1.0 / (fabs(theta)+(double)sqrt(1.0+theta*theta));
            if (theta < 0.0) 
              t = -t;
          }
          c = 1.0 / (double)sqrt(1+t*t);
          s = t * c;
          tau = s / (1.0+c);
          h = t * a[p3 + iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[p3 + iq] = 0.0;
          for (j=j3=0; j<=ip-1; j++,j3+=3) 
            ROTATE(a,j3,ip,j3,iq)
	      for (j=ip+1; j<=iq-1; j++) 
		ROTATE(a,p3,j,j*3,iq)
		  for (j=iq+1; j<3; j++) 
		    ROTATE(a,p3,j,iq*3,j)

		      for (j3=0; j3<9; j3+=3) 
			ROTATE(v,j3,ip,j3,iq)

			  ++(*nrot);
        }
      }
    }
    for (ip=0; ip<3; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }
  printf("Too many iterations in routine JACOBI\n");
  return(0);
}


   double
rms(int n, int mode, 
double *mass, int *mask,
double *toFitX, double *toFitY, double *toFitZ,
double *X, double *Y, double *Z, 
double rotation[3][3], double translation[3], int fit)
{

  int ierr=0;
  int i, modifiedCount;
  char *err;
  double rms_return;
  double *weights;
  double rot[9], rtr[9];
  int i3, k3, j, k;
  double mwss;
  double b[9], U[9];
  double *Evector[3], Eigenvalue[3], Emat[9];
  double x, y, z, xx, yy, zz;
  double total_mass;
  double sig3;
  double cp[3];
  double cofmX, cofmY, cofmZ;
  double cofmX1, cofmY1, cofmZ1;

  weights = safe_malloc( sizeof(double) * n );
  total_mass = 0.0;

  if (fit == 0) {
    /*
     *  Don't do the fit, just calculate rmsd: don't calculate 
     *  any translation/rotation 
     */
    rms_return = 0.0;
    for (i=0; i < n; i++) {
      if (mask != NULL && mask[i] == 1) {
        if (mass != NULL)
          weights[i] = mass[i];
        else
          weights[i] = 1.0;
        total_mass += weights[i];
        xx = X[i] - toFitX[i];
        yy = Y[i] - toFitY[i];
        zz = Z[i] - toFitZ[i];
        rms_return += weights[i]*(xx*xx + yy*yy + zz*zz);
      }
    }
    rms_return = sqrt(rms_return / total_mass);
    safe_free(weights);
    return (double) rms_return;
  }


  /*
   *  the rest below is for fit=1, i.e. calculate translation and
   *  rotation matrix as well as rmsd value of the fitted region 
   */

  for (i=0, modifiedCount=n; i < n; i++) {
    if ( mask != NULL && mask[i] == 0 ) {
      weights[i] = 0.0;
      modifiedCount--;
    } else
    if (mass != NULL)
      weights[i] = mass[i];
    else
      weights[i] = 1.0;
    total_mass += weights[i];
  }

  if ( mode ) 
    if ( rotation == NULL || translation == NULL )
      error("rms", "rotation matrix and translation vector are NULL?");

  if ( modifiedCount > 2 ) {

    memset(rot,  0,   sizeof(double) * 9);
    memset(rtr,  0,   sizeof(double) * 9);
    memset(U,    0,   sizeof(double) * 9);

    cofmX =  0.0;
    cofmY =  0.0;
    cofmZ =  0.0;
    cofmX1 = 0.0;
    cofmY1 = 0.0;
    cofmZ1 = 0.0;

    /*
     *  First shift the center of mass of all the atoms to be fit to
     *  the origin for both trajectory and reference coordinates.
     */
    for (k=0; k < n; k++) {
      cofmX += weights[k] * toFitX[k];
      cofmY += weights[k] * toFitY[k];
      cofmZ += weights[k] * toFitZ[k];
      cofmX1 += weights[k] * X[k];
      cofmY1 += weights[k] * Y[k];
      cofmZ1 += weights[k] * Z[k];
    }

    cofmX /= total_mass;
    cofmY /= total_mass;
    cofmZ /= total_mass;
    cofmX1 /= total_mass;
    cofmY1 /= total_mass;
    cofmZ1 /= total_mass;

    for (k=0; k < n; k++) {
      toFitX[k] -= cofmX;
      toFitY[k] -= cofmY;
      toFitZ[k] -= cofmZ;

      X[k] -= cofmX1;
      Y[k] -= cofmY1;
      Z[k] -= cofmZ1;
    }

    mwss = 0.0;
    for (k=0; k < n; k++) {

      x  = toFitX[k];
      y  = toFitY[k];
      z  = toFitZ[k];
      xx = X[k];
      yy = Y[k];
      zz = Z[k];

      mwss += weights[k] * ( x*x + y*y + z*z + xx*xx + yy*yy + zz*zz );

      /*
       *  calculate the Kabsch matrix: R = (rij) = Sum(wn*yni*xnj) 
       */
      rot[0] += weights[k] * x * xx;
      rot[1] += weights[k] * x * yy;
      rot[2] += weights[k] * x * zz;

      rot[3] += weights[k] * y * xx;
      rot[4] += weights[k] * y * yy;
      rot[5] += weights[k] * y * zz;

      rot[6] += weights[k] * z * xx;
      rot[7] += weights[k] * z * yy;
      rot[8] += weights[k] * z * zz;
    }

    mwss *= 0.5;   /* E0 = 0.5*Sum(wn*(xn^2+yn^2)) */

    /*
     *  calculate Kabsch multiplied by its transpose: RtR 
     */
    rtr[0] = rot[0]*rot[0] + rot[1]*rot[1] + rot[2]*rot[2];
    rtr[1] = rot[0]*rot[3] + rot[1]*rot[4] + rot[2]*rot[5];
    rtr[2] = rot[0]*rot[6] + rot[1]*rot[7] + rot[2]*rot[8];
    rtr[3] = rot[3]*rot[0] + rot[4]*rot[1] + rot[5]*rot[2];
    rtr[4] = rot[3]*rot[3] + rot[4]*rot[4] + rot[5]*rot[5];
    rtr[5] = rot[3]*rot[6] + rot[4]*rot[7] + rot[5]*rot[8];
    rtr[6] = rot[6]*rot[0] + rot[7]*rot[1] + rot[8]*rot[2];
    rtr[7] = rot[6]*rot[3] + rot[7]*rot[4] + rot[8]*rot[5];
    rtr[8] = rot[6]*rot[6] + rot[7]*rot[7] + rot[8]*rot[8];


    if (!diagEsort(rtr, Emat, Evector, Eigenvalue))
      return(0);

    /*
     *  a3 = a1 x a2 
     */
    VOP_3D_COORDS_CROSS_PRODUCT(Evector[2][0], Evector[2][1], Evector[2][2], 
				Evector[0][0], Evector[0][1], Evector[0][2],
				Evector[1][0], Evector[1][1], Evector[1][2]);

    /*
     *  Evector dot transpose rot:  b = R.ak 
     */
    b[0] = Evector[0][0] * rot[0] + 
    Evector[0][1] * rot[3] + 
    Evector[0][2] * rot[6];
    b[1] = Evector[0][0] * rot[1] + 
    Evector[0][1] * rot[4] + 
    Evector[0][2] * rot[7];
    b[2] = Evector[0][0] * rot[2] + 
    Evector[0][1] * rot[5] + 
    Evector[0][2] * rot[8];
    normalize(&b[0]);
    b[3] = Evector[1][0] * rot[0] + 
    Evector[1][1] * rot[3] + 
    Evector[1][2] * rot[6];
    b[4] = Evector[1][0] * rot[1] + 
    Evector[1][1] * rot[4] + 
    Evector[1][2] * rot[7];
    b[5] = Evector[1][0] * rot[2] + 
    Evector[1][1] * rot[5] + 
    Evector[1][2] * rot[8];
    normalize(&b[3]);
    b[6] = Evector[2][0] * rot[0] + 
    Evector[2][1] * rot[3] + 
    Evector[2][2] * rot[6];
    b[7] = Evector[2][0] * rot[1] + 
    Evector[2][1] * rot[4] + 
    Evector[2][2] * rot[7];
    b[8] = Evector[2][0] * rot[2] + 
    Evector[2][1] * rot[5] + 
    Evector[2][2] * rot[8];
    normalize(&b[6]);

    /*
     *  b3 = b1 x b2 
     */
    VOP_3D_COORDS_CROSS_PRODUCT(cp[0], cp[1], cp[2],
    b[0],   b[1],  b[2],
    b[3],   b[4],  b[5]);

    if ( (cp[0] * b[6] + cp[1] * b[7] + cp[2] * b[8]) < 0.0 )
      sig3 = -1.0;
    else
      sig3 = 1.0;

    b[6] = cp[0]; b[7] = cp[1]; b[8] = cp[2];

    /*
     *  U has the best rotation 
     */
    for (k=k3=0; k<3; k++,k3+=3)
      for (i=i3=0;i<3; i++,i3+=3)
    for (j=0; j<3; j++) {
      U[i3 + j] += Evector[k][j] * b[k3 + i];
    }

    /*
     *  E = E0 - sqrt(mu1) - sqrt(mu2) - sig3*sqrt(mu3) 
     */
    rms_return = mwss 
      - sqrt(fabs(Eigenvalue[0]))
      - sqrt(fabs(Eigenvalue[1]))
      - sig3 * sqrt(fabs(Eigenvalue[2]));
    if ( rms_return < 0 ) {
      rms_return = 0.0;
    } else {
      rms_return = sqrt( (2.0 * rms_return) / total_mass);
    }

    if (mode) {
      /*
       *  Save rotation matrix which does the best overlap of trajectory
       *  coordinates to reference coordinates when they are both centered
       *  on their CMs. The actual modification (=rotation) of trajectory
       *  coords happens in the calling routine (actions.c::transformRMS())
       */
      rotation[0][0] = U[0];
      rotation[0][1] = U[1];
      rotation[0][2] = U[2];
      rotation[1][0] = U[3];
      rotation[1][1] = U[4];
      rotation[1][2] = U[5];
      rotation[2][0] = U[6];
      rotation[2][1] = U[7];
      rotation[2][2] = U[8];

      /*
       *  Move the reference back so that it stays unchanged. This is
       *  necessary to preserve the meaning of CM shift on next frame
       * iteration. 
       */
      for (k=0; k < n; k++) {
        X[k] += cofmX1;
        Y[k] += cofmY1;
        Z[k] += cofmZ1;
      }
      /*
       *  Once the reference coords are shifted back to its original
       *  position (the for-cycle above), we need to shift trajectory
       *  coordinates by the same amount (i.e. CM of the reference) 
       *  to get them overlapped with the reference. The actual
       *  translation of trajectory coordinates happens in the calling
       *  routine (actions.c::transformRMS() )
       */
      translation[0] = cofmX1;
      translation[1] = cofmY1;
      translation[2] = cofmZ1;

    }
  } else
  ierr = -1;

  if (ierr != 0) {
    switch (ierr) {
      case -1:
      err = "Number of atoms less than 2";
      break;
      case -2:  /* ierr is never set to -2 previously ?? */
      err = "Illegal weights";
      break;
      default:
      err = "Unknown error";
      break;
    }
    error("rms", "KRMS_ reported %s\n", err);
  }
  safe_free(weights);
  return (double) rms_return;
}


/* This is the new rms function. f stands for FLOAT and FIT, if fit is specified.
 * if mode == 0, dispersion is calculated, no fit will be done even if fit == 1, which should not occur.
 * if mode == 1 and fit == 0, rms deviation is calculated but no structure will move.
 * if mode == 1 and fit == 1, rms deviation is calculated, XYZ moves back, but toFixXYZ's centroid moved to (0,0,0), as original functionality. Alignment will be done in the calling function.
 * if mode == 1 and fit == 2, rms deviation is calculated and toFixXYZ will align to XYZ. 
 */
   double
rmsf(int n, int mode, 
double *mass, int *mask,
float *toFitX, float *toFitY, float *toFitZ,
float *X, float *Y, float *Z, 
float rotation[3][3], float translation[3], int fit)
{

  int ierr=0;
  int i, modifiedCount;
  char *err;
  double rms_return;
  double *weights;
  double rot[9], rtr[9];
  int i3, k3, j, k;
  double mwss;
  double b[9], U[9];
  double *Evector[3], Eigenvalue[3], Emat[9];
  double x, y, z, xx, yy, zz;
  double total_mass;
  double sig3;
  double cp[3];
  double cofmX, cofmY, cofmZ;
  double cofmX1, cofmY1, cofmZ1;
  double xtemp, ytemp, ztemp;

  weights = safe_malloc( sizeof(double) * n );
  total_mass = 0.0;

  if (mode == 0) {
    /*
     *  Don't do the fit, just calculate rmsd: don't calculate 
     *  any translation/rotation 
     */
    rms_return = 0.0;
    for (i=0; i < n; i++) {
      if (mask != NULL && mask[i] == 1) {
        if (mass != NULL)
          weights[i] = mass[i];
        else
          weights[i] = 1.0;
        total_mass += weights[i];
        xx = X[i] - toFitX[i];
        yy = Y[i] - toFitY[i];
        zz = Z[i] - toFitZ[i];
        rms_return += weights[i]*(xx*xx + yy*yy + zz*zz);
      }
    }
    rms_return = sqrt(rms_return / total_mass);
    safe_free(weights);
    return (double) rms_return;
  }


  /*
   *  the rest below is for fit=1, i.e. calculate translation and
   *  rotation matrix as well as rmsd value of the fitted region 
   */

  for (i=0, modifiedCount=n; i < n; i++) {
    if ( mask != NULL && mask[i] == 0 ) {
      weights[i] = 0.0;
      modifiedCount--;
    } else
    if (mass != NULL)
      weights[i] = mass[i];
    else
      weights[i] = 1.0;
    total_mass += weights[i];
  }

  if ( mode ) 
    if ( rotation == NULL || translation == NULL )
      error("rms", "rotation matrix and translation vector are NULL?");

  if ( modifiedCount > 2 ) {

    memset(rot,  0,   sizeof(double) * 9);
    memset(rtr,  0,   sizeof(double) * 9);
    memset(U,    0,   sizeof(double) * 9);

    cofmX =  0.0;
    cofmY =  0.0;
    cofmZ =  0.0;
    cofmX1 = 0.0;
    cofmY1 = 0.0;
    cofmZ1 = 0.0;

    /*
     *  First shift the center of mass of all the atoms to be fit to
     *  the origin for both trajectory and reference coordinates.
     */
    for (k=0; k < n; k++) {
      cofmX += weights[k] * toFitX[k];
      cofmY += weights[k] * toFitY[k];
      cofmZ += weights[k] * toFitZ[k];
      cofmX1 += weights[k] * X[k];
      cofmY1 += weights[k] * Y[k];
      cofmZ1 += weights[k] * Z[k];
    }

    cofmX /= total_mass;
    cofmY /= total_mass;
    cofmZ /= total_mass;
    cofmX1 /= total_mass;
    cofmY1 /= total_mass;
    cofmZ1 /= total_mass;

    for (k=0; k < n; k++) {
      toFitX[k] -= cofmX;
      toFitY[k] -= cofmY;
      toFitZ[k] -= cofmZ;

      X[k] -= cofmX1;
      Y[k] -= cofmY1;
      Z[k] -= cofmZ1;
    }

    mwss = 0.0;
    for (k=0; k < n; k++) {

      x  = toFitX[k];
      y  = toFitY[k];
      z  = toFitZ[k];
      xx = X[k];
      yy = Y[k];
      zz = Z[k];

      mwss += weights[k] * ( x*x + y*y + z*z + xx*xx + yy*yy + zz*zz );

      /*
       *  calculate the Kabsch matrix: R = (rij) = Sum(wn*yni*xnj) 
       */
      rot[0] += weights[k] * x * xx;
      rot[1] += weights[k] * x * yy;
      rot[2] += weights[k] * x * zz;

      rot[3] += weights[k] * y * xx;
      rot[4] += weights[k] * y * yy;
      rot[5] += weights[k] * y * zz;

      rot[6] += weights[k] * z * xx;
      rot[7] += weights[k] * z * yy;
      rot[8] += weights[k] * z * zz;
    }

    mwss *= 0.5;   /* E0 = 0.5*Sum(wn*(xn^2+yn^2)) */

    /*
     *  calculate Kabsch multiplied by its transpose: RtR 
     */
    rtr[0] = rot[0]*rot[0] + rot[1]*rot[1] + rot[2]*rot[2];
    rtr[1] = rot[0]*rot[3] + rot[1]*rot[4] + rot[2]*rot[5];
    rtr[2] = rot[0]*rot[6] + rot[1]*rot[7] + rot[2]*rot[8];
    rtr[3] = rot[3]*rot[0] + rot[4]*rot[1] + rot[5]*rot[2];
    rtr[4] = rot[3]*rot[3] + rot[4]*rot[4] + rot[5]*rot[5];
    rtr[5] = rot[3]*rot[6] + rot[4]*rot[7] + rot[5]*rot[8];
    rtr[6] = rot[6]*rot[0] + rot[7]*rot[1] + rot[8]*rot[2];
    rtr[7] = rot[6]*rot[3] + rot[7]*rot[4] + rot[8]*rot[5];
    rtr[8] = rot[6]*rot[6] + rot[7]*rot[7] + rot[8]*rot[8];


    if (!diagEsort(rtr, Emat, Evector, Eigenvalue))
      return(0);

    /*
     *  a3 = a1 x a2 
     */
    VOP_3D_COORDS_CROSS_PRODUCT(Evector[2][0], Evector[2][1], Evector[2][2], 
				Evector[0][0], Evector[0][1], Evector[0][2],
				Evector[1][0], Evector[1][1], Evector[1][2]);

    /*
     *  Evector dot transpose rot:  b = R.ak 
     */
    b[0] = Evector[0][0] * rot[0] + 
    Evector[0][1] * rot[3] + 
    Evector[0][2] * rot[6];
    b[1] = Evector[0][0] * rot[1] + 
    Evector[0][1] * rot[4] + 
    Evector[0][2] * rot[7];
    b[2] = Evector[0][0] * rot[2] + 
    Evector[0][1] * rot[5] + 
    Evector[0][2] * rot[8];
    normalize(&b[0]);
    b[3] = Evector[1][0] * rot[0] + 
    Evector[1][1] * rot[3] + 
    Evector[1][2] * rot[6];
    b[4] = Evector[1][0] * rot[1] + 
    Evector[1][1] * rot[4] + 
    Evector[1][2] * rot[7];
    b[5] = Evector[1][0] * rot[2] + 
    Evector[1][1] * rot[5] + 
    Evector[1][2] * rot[8];
    normalize(&b[3]);
    b[6] = Evector[2][0] * rot[0] + 
    Evector[2][1] * rot[3] + 
    Evector[2][2] * rot[6];
    b[7] = Evector[2][0] * rot[1] + 
    Evector[2][1] * rot[4] + 
    Evector[2][2] * rot[7];
    b[8] = Evector[2][0] * rot[2] + 
    Evector[2][1] * rot[5] + 
    Evector[2][2] * rot[8];
    normalize(&b[6]);

    /*
     *  b3 = b1 x b2 
     */
    VOP_3D_COORDS_CROSS_PRODUCT(cp[0], cp[1], cp[2],
    b[0],   b[1],  b[2],
    b[3],   b[4],  b[5]);

    if ( (cp[0] * b[6] + cp[1] * b[7] + cp[2] * b[8]) < 0.0 )
      sig3 = -1.0;
    else
      sig3 = 1.0;

    b[6] = cp[0]; b[7] = cp[1]; b[8] = cp[2];

    /*
     *  U has the best rotation 
     */
    for (k=k3=0; k<3; k++,k3+=3)
      for (i=i3=0;i<3; i++,i3+=3)
        for (j=0; j<3; j++) {
          U[i3 + j] += Evector[k][j] * b[k3 + i];
        }

    /*
     *  E = E0 - sqrt(mu1) - sqrt(mu2) - sig3*sqrt(mu3) 
     */
    rms_return = mwss 
      - sqrt(fabs(Eigenvalue[0]))
      - sqrt(fabs(Eigenvalue[1]))
      - sig3 * sqrt(fabs(Eigenvalue[2]));
    if ( rms_return < 0 ) {
      rms_return = 0.0;
    } else {
      rms_return = sqrt( (2.0 * rms_return) / total_mass);
    }

    if (mode) {
      /*
       *  Save rotation matrix which does the best overlap of trajectory
       *  coordinates to reference coordinates when they are both centered
       *  on their CMs. The actual modification (=rotation) of trajectory
       *  coords happens in the calling routine (actions.c::transformRMS())
       */
      rotation[0][0] = U[0];
      rotation[0][1] = U[1];
      rotation[0][2] = U[2];
      rotation[1][0] = U[3];
      rotation[1][1] = U[4];
      rotation[1][2] = U[5];
      rotation[2][0] = U[6];
      rotation[2][1] = U[7];
      rotation[2][2] = U[8];

      /*
       *  Move the reference back so that it stays unchanged. This is
       *  necessary to preserve the meaning of CM shift on next frame
       * iteration. 
       */
      for (k=0; k < n; k++) {
        X[k] += cofmX1;
        Y[k] += cofmY1;
        Z[k] += cofmZ1;
      }
      /*
       *  Once the reference coords are shifted back to its original
       *  position (the for-cycle above), we need to shift trajectory
       *  coordinates by the same amount (i.e. CM of the reference) 
       *  to get them overlapped with the reference. The actual
       *  translation of trajectory coordinates happens in the calling
       *  routine (actions.c::transformRMS() )
       */
      translation[0] = cofmX1;
      translation[1] = cofmY1;
      translation[2] = cofmZ1;
      
      if (fit == 2) {
      /* First apply the rotation (which was calculated for both 
         trajectory and reference coords shifted to their CMs). The
         order (first rotation, then translation) is important.*/
        for (k=0; k < n; k++) {
          VOP_3x3_TIMES_COORDS(rotation, toFitX[k], toFitY[k], toFitZ[k], xtemp, ytemp, ztemp);
          toFitX[k] += cofmX1;
          toFitY[k] += cofmY1;
          toFitZ[k] += cofmZ1;
        }
      } else if (fit == 1){
        /* Nothing. XYZ moved back. ToFitXYZ moved to (0,0,0) */
      } else if (fit == 0){
      /* Or just move them back to their original position */
        for (k=0; k < n; k++) {
          toFitX[k] += cofmX;
          toFitY[k] += cofmY;
          toFitZ[k] += cofmZ;
        }
      } 
    }
  } else
  ierr = -1;

  if (ierr != 0) {
    switch (ierr) {
      case -1:
      err = "Number of atoms less than 2";
      break;
      case -2:  /* ierr is never set to -2 previously ?? */
      err = "Illegal weights";
      break;
      default:
      err = "Unknown error";
      break;
    }
    error("rms", "KRMS_ reported %s\n", err);
  }
  safe_free(weights);
  return (double) rms_return;
}


   double
rmsf_old(int n, int mode, 
double *mass, int *mask,
float *toFitX, float *toFitY, float *toFitZ,
float *X, float *Y, float *Z, 
float rotation[3][3], float translation[3], int fit)
{

  int ierr=0;
  int i, modifiedCount;
  char *err;
  double rms_return;
  double *weights;
  double rot[9], rtr[9];
  int i3, k3, j, k;
  double mwss;
  double b[9], U[9];
  double *Evector[3], Eigenvalue[3], Emat[9];
  double x, y, z, xx, yy, zz;
  double total_mass;
  double sig3;
  double cp[3];
  double cofmX, cofmY, cofmZ;
  double cofmX1, cofmY1, cofmZ1;

  weights = safe_malloc( sizeof(double) * n );
  total_mass = 0.0;
  for (i=0, modifiedCount=n; i < n; i++) {
    if ( mask != NULL && mask[i] == 0 ) {
      weights[i] = 0.0;
      modifiedCount--;
    } else
    if (mass != NULL)
      weights[i] = mass[i];
    else
      weights[i] = 1.0;
    total_mass += weights[i];
  }

  if ( mode ) 
    if ( rotation == NULL || translation == NULL )
      error("rms", "rotation matrix and translation vector are NULL?");

  if ( modifiedCount > 2 ) {
    memset(rot,  0,   sizeof(double) * 9);
    memset(rtr,  0,   sizeof(double) * 9);
    memset(U,    0,   sizeof(double) * 9);

    cofmX =  0.0;
    cofmY =  0.0;
    cofmZ =  0.0;
    cofmX1 = 0.0;
    cofmY1 = 0.0;
    cofmZ1 = 0.0;

    /*
     *  First calculate CMs but don't really shift any coordinates
     */
    for (k=0; k < n; k++) {
      cofmX += weights[k] * toFitX[k];
      cofmY += weights[k] * toFitY[k];
      cofmZ += weights[k] * toFitZ[k];
      cofmX1 += weights[k] * X[k];
      cofmY1 += weights[k] * Y[k];
      cofmZ1 += weights[k] * Z[k];
    }

    cofmX /= total_mass;
    cofmY /= total_mass;
    cofmZ /= total_mass;
    cofmX1 /= total_mass;
    cofmY1 /= total_mass;
    cofmZ1 /= total_mass;

    mwss = 0.0;
    for (k=0; k < n; k++) {

      x  = toFitX[k]-cofmX;
      y  = toFitY[k]-cofmY;
      z  = toFitZ[k]-cofmZ;
      xx = X[k]-cofmX1;
      yy = Y[k]-cofmY1;
      zz = Z[k]-cofmZ1;

      mwss += weights[k] * ( x*x + y*y + z*z + xx*xx + yy*yy + zz*zz );

      /*
       *  calculate the Kabsch matrix: R = (rij) = Sum(wn*yni*xnj) 
       */
      rot[0] += weights[k] * x * xx;
      rot[1] += weights[k] * x * yy;
      rot[2] += weights[k] * x * zz;

      rot[3] += weights[k] * y * xx;
      rot[4] += weights[k] * y * yy;
      rot[5] += weights[k] * y * zz;

      rot[6] += weights[k] * z * xx;
      rot[7] += weights[k] * z * yy;
      rot[8] += weights[k] * z * zz;
    }

    mwss *= 0.5;   /* E0 = 0.5*Sum(wn*(xn^2+yn^2)) */

    /*
     *  calculate Kabsch multiplied by its transpose: RtR 
     */
    rtr[0] = rot[0]*rot[0] + rot[1]*rot[1] + rot[2]*rot[2];
    rtr[1] = rot[0]*rot[3] + rot[1]*rot[4] + rot[2]*rot[5];
    rtr[2] = rot[0]*rot[6] + rot[1]*rot[7] + rot[2]*rot[8];
    rtr[3] = rot[3]*rot[0] + rot[4]*rot[1] + rot[5]*rot[2];
    rtr[4] = rot[3]*rot[3] + rot[4]*rot[4] + rot[5]*rot[5];
    rtr[5] = rot[3]*rot[6] + rot[4]*rot[7] + rot[5]*rot[8];
    rtr[6] = rot[6]*rot[0] + rot[7]*rot[1] + rot[8]*rot[2];
    rtr[7] = rot[6]*rot[3] + rot[7]*rot[4] + rot[8]*rot[5];
    rtr[8] = rot[6]*rot[6] + rot[7]*rot[7] + rot[8]*rot[8];


    if (!diagEsort(rtr, Emat, Evector, Eigenvalue))
      return(0);


    /*
     *  a3 = a1 x a2 
     */
    VOP_3D_COORDS_CROSS_PRODUCT(Evector[2][0], Evector[2][1], Evector[2][2], 
				Evector[0][0], Evector[0][1], Evector[0][2],
				Evector[1][0], Evector[1][1], Evector[1][2]);


    /*
     *  Evector dot transpose rot:  b = R.ak 
     */
    b[0] = Evector[0][0] * rot[0] + 
           Evector[0][1] * rot[3] + 
           Evector[0][2] * rot[6];
    b[1] = Evector[0][0] * rot[1] + 
           Evector[0][1] * rot[4] + 
           Evector[0][2] * rot[7];
    b[2] = Evector[0][0] * rot[2] + 
           Evector[0][1] * rot[5] + 
           Evector[0][2] * rot[8];
    normalize(&b[0]);
    b[3] = Evector[1][0] * rot[0] + 
           Evector[1][1] * rot[3] + 
           Evector[1][2] * rot[6];
    b[4] = Evector[1][0] * rot[1] + 
           Evector[1][1] * rot[4] + 
           Evector[1][2] * rot[7];
    b[5] = Evector[1][0] * rot[2] + 
           Evector[1][1] * rot[5] + 
           Evector[1][2] * rot[8];
    normalize(&b[3]);
    b[6] = Evector[2][0] * rot[0] + 
           Evector[2][1] * rot[3] + 
           Evector[2][2] * rot[6];
    b[7] = Evector[2][0] * rot[1] + 
           Evector[2][1] * rot[4] + 
           Evector[2][2] * rot[7];
    b[8] = Evector[2][0] * rot[2] + 
           Evector[2][1] * rot[5] + 
           Evector[2][2] * rot[8];
    normalize(&b[6]);

    /*
     *  b3 = b1 x b2 
     */
    VOP_3D_COORDS_CROSS_PRODUCT(cp[0], cp[1], cp[2],
				b[0],   b[1],  b[2],
				b[3],   b[4],  b[5]);

    if ( (cp[0] * b[6] + cp[1] * b[7] + cp[2] * b[8]) < 0.0 )
      sig3 = -1.0;
    else
      sig3 = 1.0;

    b[6] = cp[0]; b[7] = cp[1]; b[8] = cp[2];

    /*
     *  U has the best rotation
     */
    for (k=k3=0; k<3; k++,k3+=3)
      for (i=i3=0;i<3; i++,i3+=3)
    for (j=0; j<3; j++) {
      U[i3 + j] += Evector[k][j] * b[k3 + i];
    }

    /*
     *  E = E0 - sqrt(mu1) - sqrt(mu2) - sig3*sqrt(mu3) 
     */
    rms_return = mwss 
      - sqrt(fabs(Eigenvalue[0]))
      - sqrt(fabs(Eigenvalue[1]))
      - sig3 * sqrt(fabs(Eigenvalue[2]));
    if ( rms_return < 0 ) {
      rms_return = 0.0;
    } else {
      rms_return = sqrt( (2.0 * rms_return) / total_mass);
    }

    if (mode) {
      /*
       *  Save rotation matrix which does the best overlap of trajectory
       *  coordinates to reference coordinates when they are both centered
       *  on their CMs. The actual modification (=rotation) of trajectory
       *  coords happens in the calling routine (actions.c::transformRMS())
       */
      rotation[0][0] = U[0];
      rotation[0][1] = U[1];
      rotation[0][2] = U[2];
      rotation[1][0] = U[3];
      rotation[1][1] = U[4];
      rotation[1][2] = U[5];
      rotation[2][0] = U[6];
      rotation[2][1] = U[7];
      rotation[2][2] = U[8];

      /*
       *  Move the reference back so that it stays unchanged. This is
       *  necessary to preserve the meaning of CM shift on next frame
       * iteration.
       */
      for (k=0; k < n; k++) {
        X[k] += cofmX1;
        Y[k] += cofmY1;
        Z[k] += cofmZ1;
      }
      /*
       *  Once the reference coords are shifted back to its original
       *  position (the for-cycle above), we need to shift trajectory
       *  coordinates by the same amount (i.e. CM of the reference) 
       *  to get them overlapped with the reference. The actual
       *  translation of trajectory coordinates happens in the calling
       *  routine (actions.c::transformRMS() )
       */
      translation[0] = cofmX1;
      translation[1] = cofmY1;
      translation[2] = cofmZ1;
    }
  } else
  ierr = -1;

  if (ierr != 0) {
    switch (ierr) {
      case -1:
      err = "Number of atoms less than 2";
      break;
      case -2:  /* ierr is never set to -2 previously ?? */
      err = "Illegal weights";
      break;
      default:
      err = "Unknown error";
      break;
    }
    error("rms", "KRMS_ reported %s\n", err);
  }
  safe_free(weights);
  return (double) rms_return;
}
