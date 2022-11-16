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
 *  $Header: /home/case/cvsroot/amber11/AmberTools/src/ptraj/correlation.c,v 10.0 2008/04/15 23:24:11 case Exp $
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



/** ACTION ROUTINE *************************************************************
 *
 *  transformCorr()   --- perform correlation analysis (Vickie Tsui, Scripps)
 *
 *  Supplementary routines:
 *    compute_corr()
 *
 ******************************************************************************/


void compute_corr(char *outfile, double *x, double *y, double *z, 
		int tcorr, int totalFrames)

/*   x,y,z arrays contain the coordinates of the (unnormalized) vector
           being considered
     tcorr is the maximum length to compute time correlation functions,
           in units of "Frames"
    totalFrames  is the total number of frames being provided              */
     
{
  typedef struct _complex {
    double real;
    double imaginary;
  } complex;

  int i,j, current_frame;
  double rmag0, rmagt;
  double x0, y0, z0, xt, yt, zt;
  double *p2, *corr, *rcorr;
  double r6ave, r3ave, rave, avecrd[4], rrig;
  double dot, y2asy, y20;
  double th0, phi0;
  double rfac, qfac;
  complex y21, y21c, y22, y22c;
  int atind[3], *cfind, npts, ncorr, ntot, ftyp;
  int doit, jmax;
  FILE *ifp;

  /* allocate space */

  p2     = (double *) safe_malloc( sizeof(double) * tcorr );
  corr   = (double *) safe_malloc( sizeof(double) * tcorr );
  rcorr  = (double *) safe_malloc( sizeof(double) * tcorr );

  /* initialize */

  for (i=1; i<=tcorr; ++i)  { corr[i]=p2[i]=rcorr[i]=0.0; }
  r6ave = r3ave = 0.0;
  avecrd[1] = avecrd[2] = avecrd[3] = 0.0;
  rave = y2asy = y20 = 0.0;
  y21.real = y21.imaginary = 0.0;
  y21c.real = y21c.imaginary = 0.0;
  y22.real = y22.imaginary = 0.0;
  y22c.real = y22c.imaginary = 0.0;

/* main loop for calculating correlation functions */

  for( i=0; i<totalFrames; i++ ){

    /* computations for any straighforward averages go here   */

    rmag0= sqrt( x[i]*x[i] + y[i]*y[i] + z[i]*z[i] );

    x0=x[i]/rmag0;
    y0=y[i]/rmag0;
    z0=z[i]/rmag0;

    r6ave=r6ave+1/pow(rmag0,6);
    r3ave=r3ave+1/pow(rmag0,3);
    rave=rave+rmag0;
    avecrd[1]=avecrd[1]+x[i];
    avecrd[2]=avecrd[2]+y[i];
    avecrd[3]=avecrd[3]+z[i];

    th0=acos(z0);
    phi0=atan2(y0,x0);

    y22.real+=sqrt(3.0/4.0)*pow((sin(th0)),2)*(cos(2*phi0))/pow(rmag0,3);
    y22.imaginary+=sqrt(3.0/4.0)*pow((sin(th0)),2)*(sin(2*phi0))/pow(rmag0,3);
    y22c.real+=sqrt(3.0/4.0)*pow((sin(th0)),2)*(cos(2*phi0))/pow(rmag0,3);
    y22c.imaginary+=sqrt(3.0/4.0)*pow((sin(th0)),2)*(-sin(2*phi0))/pow(rmag0,3);
    y21.real+=sqrt(3.0)*cos(th0)*sin(th0)*cos(phi0)/pow(rmag0,3);
    y21.imaginary+=sqrt(3.0)*cos(th0)*sin(th0)*sin(phi0)/pow(rmag0,3);
    y21c.real+=sqrt(3.0)*cos(th0)*sin(th0)*cos(phi0)/pow(rmag0,3);
    y21c.imaginary+=sqrt(3.0)*cos(th0)*sin(th0)*(-sin(phi0))/pow(rmag0,3);
    y20+=(pow((3*(cos(th0))),2)-1)/(2.0*pow(rmag0,3));

    /*  now, loop over frames that follow frame i in the list:   */

    jmax = i+tcorr < totalFrames ? i+tcorr : totalFrames;
    for (j=i; j<jmax; j++)  {

      rmagt=sqrt( x[j]*x[j] + y[j]*y[j] + z[j]*z[j] );
      xt=x[j]/rmagt;
      yt=y[j]/rmagt;
      zt=z[j]/rmagt;
      dot=(3.0*pow((x0*xt+y0*yt+z0*zt),2)-1.0)/2.0;

      corr[j-i] += dot/pow((rmag0*rmagt),3);
      p2[j-i] += dot;
      rcorr[j-i] += 1.0/pow((rmag0*rmagt),3);

    }
  }

  /* normalize simple averages:  */

  r6ave /= totalFrames;
  r3ave /= totalFrames;
  rave /= totalFrames;
  avecrd[1] /= totalFrames;
  avecrd[2] /= totalFrames;
  avecrd[3] /= totalFrames;
  rrig = sqrt( avecrd[1]*avecrd[1] + avecrd[2]*avecrd[2] + avecrd[3]*avecrd[3] );

  y2asy=(y22.real*y22c.real+y21.real*y21c.real)+pow(y20,2);
  y2asy=y2asy/(totalFrames*totalFrames*r6ave);

  /* normalize the correlation functions:  */

  for (i=0; i<tcorr; i++)  {
    npts = i+tcorr < totalFrames ? tcorr : totalFrames - i;
    corr[i] = corr[i]/(npts*r6ave);
    rcorr[i]=rcorr[i]/(npts*r6ave);
    p2[i] /= npts;
  }

  /* output correlation functions */

  ifp=safe_fopen(outfile, "w");
  if (ifp == NULL) {
    warning("ptraj(), correlation: cannot open output file %s\n",
	    outfile);
  } else {
    fprintf(ifp, "# Rrigid= %lf  Rave= %lf \n", rrig, rave);
    fprintf(ifp, "# <1/r^3>= %lf  <1/r^6>= %lf\n", r3ave, r6ave);
    rfac = r6ave*pow(rave,6);
    qfac = y2asy*rfac;
    fprintf(ifp, "#   time     C(t)      P2(t)      R(t)\n");
    for (j=0; j<tcorr; j++)  {
      fprintf(ifp, "%d   %lf   %lf   %lf\n", j,corr[j],p2[j],rcorr[j]);
    }
    safe_fclose(ifp);
  }

  /* deallocate space */

  safe_free(p2);
  safe_free(corr);
  safe_free(rcorr);

}


int transformCorr(actionInformation *action, 
	      double *x, double *y, double *z, double *box, int mode)
{

  stackType **argumentStackPointer;
  char *buffer;
  transformCorrInfo *corrInfo;
  ptrajState *state;
  int i;
  double *principal;
  double cx, cy, cz, total_mass;
  double vx, vy, vz;


  /*
   *  USAGE:
   *
   *     correlation name mask1 mask2 tmin tcorr tmax [out filename]
   *
   *  action argument usage:
   *
   *  carg1:
   *     a transformCorrInfo structure
   */


  if (mode == PTRAJ_SETUP) {
    /*
     *  ACTION: PTRAJ_SETUP
     */

    argumentStackPointer = (stackType **) action->carg1;
    action->carg1 = NULL;

    /*
     *  set up complex argument
     */
    corrInfo = (transformCorrInfo *) safe_malloc(sizeof(transformCorrInfo));
    INITIALIZE_transformCorrInfo(corrInfo);
    corrInfo->totalFrames = -1;

    corrInfo->name = getArgumentString(argumentStackPointer, NULL);

    buffer = getArgumentString(argumentStackPointer, NULL);
    corrInfo->mask = processAtomMask(buffer, action->state);
    safe_free(buffer);

    buffer = getArgumentString(argumentStackPointer, NULL);
    corrInfo->mask2 = processAtomMask(buffer, action->state);
    safe_free(buffer);
    corrInfo->mode = VECTOR_MASK;

    corrInfo->tmin =  getArgumentInteger(argumentStackPointer, 1.0);
    corrInfo->tcorr = getArgumentInteger(argumentStackPointer, 1.0);
    corrInfo->tmax  = getArgumentInteger(argumentStackPointer, 1.0);

    /*
     *  assume "out" may be missing
     */

    corrInfo->filename = argumentStackKeyToString(argumentStackPointer, "out", NULL);
    if (corrInfo->filename == NULL) {
      corrInfo->filename = getArgumentString(argumentStackPointer, NULL);
      if (corrInfo->filename == NULL) {
	error("ptraj()", "correlation, no out file specified\n");
      }
    }
    if (corrInfo->name == NULL || corrInfo->mask == NULL ||
	corrInfo->mask2 == NULL) {
      error("ptraj()", "correlation arguments\n");
    }

    action->carg1 = (void *) corrInfo;

    return 0;
  }


  corrInfo = (transformCorrInfo *) action->carg1;

  if (mode == PTRAJ_STATUS) {

    /*
     *  ACTION: PTRAJ_STATUS
     */
    fprintf(stdout, "  CORRELATION: storage to array named %s",
            corrInfo->name);
    fprintf(stdout, " -- tmin: %i tcorr: %i tmax: %i\n",
	    corrInfo->tmin, corrInfo->tcorr, corrInfo->tmax);
    fprintf(stdout, "      Atom selection 1 is ");
    printAtomMask(corrInfo->mask, action->state);
    fprintf(stdout, "      Atom selection 2 is ");
    printAtomMask(corrInfo->mask2, action->state);

  } else if (mode == PTRAJ_PRINT) {

    /*
     *  ACTION: PTRAJ_PRINT
     */

    fprintf(stdout, "PTRAJ CORRELATION: calculating correlation %s\n",
	    corrInfo->name);
    if (corrInfo != NULL) {
      compute_corr(corrInfo->filename, corrInfo->vx, corrInfo->vy, corrInfo->vz,
		   corrInfo->tcorr, corrInfo->frame );
    }
    return 0;

  } else if (mode == PTRAJ_CLEANUP) {

    /*
     *  ACTION: PTRAJ_CLEANUP
     */

    safe_free(corrInfo->cx);
    safe_free(corrInfo->cy);
    safe_free(corrInfo->cz);
    safe_free(corrInfo->vx);
    safe_free(corrInfo->vy);
    safe_free(corrInfo->vz);
    safe_free(corrInfo->mask);
    safe_free(corrInfo->mask2);
    safe_free(corrInfo->name);
    INITIALIZE_transformCorrInfo(corrInfo);

    safe_free(corrInfo);
  }


  if (mode != PTRAJ_ACTION) return 0;

  /*
   *  ACTION: PTRAJ_ACTION
   */

  state = (ptrajState *) action->state;

  if (corrInfo->totalFrames < 0) {
    corrInfo->totalFrames = state->maxFrames;
    corrInfo->cx = (double *)
      safe_malloc(sizeof(double) * (corrInfo->tmax - corrInfo->tmin) );
    corrInfo->cy = (double *)
      safe_malloc(sizeof(double) * (corrInfo->tmax - corrInfo->tmin) );
    corrInfo->cz = (double *)
      safe_malloc(sizeof(double) * (corrInfo->tmax - corrInfo->tmin) );
    corrInfo->vx = (double *)
      safe_malloc(sizeof(double) * (corrInfo->tmax - corrInfo->tmin) );
    corrInfo->vy = (double *)
      safe_malloc(sizeof(double) * (corrInfo->tmax - corrInfo->tmin) );
    corrInfo->vz = (double *)
      safe_malloc(sizeof(double) * (corrInfo->tmax - corrInfo->tmin) );
  }


  if (corrInfo->frame > corrInfo->totalFrames) {
    warning("transformCorrr()", "Blowing array; too many frames!!\n");
    return 0;
  }

  if( corrInfo->frame < corrInfo->tmin || corrInfo->frame > corrInfo->tmax ){
    corrInfo->frame++;
    return 1;
  }

  corrInfo->mode = CORR_MASK;
  total_mass = 0.0;
  cx = 0.0;
  cy = 0.0;
  cz = 0.0;
  for (i=0; i < state->atoms; i++) {
    if (corrInfo->mask[i]) {
        cx += state->masses[i] * x[i];
        cy += state->masses[i] * y[i];
        cz += state->masses[i] * z[i];
        total_mass += state->masses[i];
    }
  }
  cx = cx / total_mass;
  cy = cy / total_mass;
  cz = cz / total_mass;

  total_mass = 0.0;
  vx = 0.0;
  vy = 0.0;
  vz = 0.0;
  for (i=0; i < state->atoms; i++) {
    if (corrInfo->mask2[i]) {
        vx += state->masses[i] * x[i];
        vy += state->masses[i] * y[i];
        vz += state->masses[i] * z[i];
        total_mass += state->masses[i];
    }
  }
  vx = vx / total_mass;
  vy = vy / total_mass;
  vz = vz / total_mass;

  corrInfo->vx[corrInfo->frame - corrInfo->tmin] = vx - cx;
  corrInfo->vy[corrInfo->frame - corrInfo->tmin] = vy - cy;
  corrInfo->vz[corrInfo->frame - corrInfo->tmin] = vz - cz;
  corrInfo->cx[corrInfo->frame - corrInfo->tmin] = cx;
  corrInfo->cy[corrInfo->frame - corrInfo->tmin] = cy;
  corrInfo->cz[corrInfo->frame - corrInfo->tmin] = cz;

  corrInfo->frame++;
  return 1;
}
