#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrix.h"
#include "vector.h"
#include "macros.h"

/***=======================================================================***/
/*** TRED2: reduces a symmetric matrix to tridiagonal form.  It is a       ***/
/***        modification more suited for C programs from [REF]:            ***/
/***                                                                       ***/
/***        William H. Press, Saul A. Teukolsky, William T. Vetterling,    ***/
/***        and Brian P. Flannery.  Numerical Recipes in C, Second         ***/
/***        Edition.  Cambridge University Press, 1992.                    ***/
/***=======================================================================***/
void TRED2(double** A, int n, double* d, double* e)
{
  int i, j, k, l;
  double scale, hh, h, g, f;
  double* tmp;
  double* tm2p;

  for (i = n-1; i >= 1; i--) {
    l = i - 1;
    h = 0.0;
    scale = 0.0;
    tmp = A[i];
    if (l > 0) {
      for (k = 0; k <= l; k++) {
	scale += fabs(tmp[k]);
      }
      if (scale == 0.0) {
	e[i] = tmp[l];
      }
      else {
	for (k = 0; k <= l; k++) {
	  tmp[k] /= scale;
	  h += tmp[k]*tmp[k];
	}
	f = tmp[l];
	g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
	e[i] = scale*g;
	h -= f*g;
	tmp[l] = f - g;
	f = 0.0;
	for (j = 0; j <= l; j++) {
	  tm2p = A[j];
	  tm2p[i] = tmp[j]/h;
	  g = 0.0;
	  for (k = 0; k <= j; k++) {
	    g += tm2p[k]*tmp[k];
	  }
	  for (k = j+1; k <=l; k++) {
	    g += A[k][j]*tmp[k];
	  }
	  e[j] = g/h;
	  f += e[j]*tmp[j];
	}
	hh = f/(h + h);
	for (j = 0; j <= l; j++) {
	  f = tmp[j];
	  e[j] = g = e[j] - hh*f;
	  tm2p = A[j];
	  for (k = 0; k <= j; k++) {
	    tm2p[k] -= (f*e[k] + g*tmp[k]);
	  }
	}
      }
    }
    else {
      e[i] = tmp[l];
    }
    d[i] = h;
  }
  d[0] = 0.0;
  e[0] = 0.0;
  for (i = 0; i < n; i++) {
    tmp = A[i];
    l = i - 1;
    if (d[i]) {
      for (j = 0; j <= l; j++) {
        g = 0.0;
        for (k = 0; k <= l; k++) {
          g += tmp[k]*A[k][j];
        }
        for (k = 0; k <= l; k++) {
          A[k][j] -= g*A[k][i];
        }
      }
    }
    d[i] = tmp[i];
    tmp[i] = 1.0;
    for (j = 0; j <= l; j++) {
      A[j][i] = tmp[j] = 0.0;
    }
  }
}

/***=======================================================================***/
/*** TQLI: this function is also modified from [REF]:                      ***/
/***                                                                       ***/
/***        William H. Press, Saul A. Teukolsky, William T. Vetterling,    ***/
/***        and Brian P. Flannery.  Numerical Recipes in C, Second         ***/
/***        Edition.  Cambridge University Press, 1992.                    ***/
/***                                                                       ***/
/***        The input is the super- and sub- diagonal vectors d and e.     ***/
/***        The eigenvalues are returned in d and the eigenvectors are     ***/
/***        returned in Z.                                                 ***/
/***=======================================================================***/
void TQLI(double* d, double* e, int n, double** z)
{
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;
  double* tmp;

  for (i = 1; i < n; i++) {
    e[i-1] = e[i];
  }
  e[n-1] = 0.0;
  for (l = 0; l < n; l++) {
    iter = 0;

    /*** Added to initialize m ***/
    m = l - 1;

    while (m != l) {
      for (m = l; m < n - 1; m++) {
	dd = fabs(d[m]) + fabs(d[m+1]);
	if (fabs(e[m]+dd) == dd) break;
      }
      if (m != l) {
	if (iter++ == 30) {
	  printf("TQLI >> Error: Too many iterations in tqli\n");
	}
	g = (d[l+1]-d[l])/(2.0*e[l]);
	r = pythag(g,1.0);
	g = d[m]-d[l]+e[l]/(g+SIGN2(r,g));
	c = 1.0;
	s = 1.0;
	p = 0.0;
	for (i = m-1; i >= l; i--) {
	  f = s*e[i];
	  b = c*e[i];
	  e[i+1] = (r=pythag(f,g));
	  if (r == 0.0) {
	    d[i+1] -= p;
	    e[m] = 0.0;
	    break;
	  }
	  s = f/r;
	  c = g/r;
	  g = d[i+1]-p;
	  r = (d[i] - g)*s + 2.0*c*b;
	  d[i+1] = g + (p = s*r);
	  g = c*r - b;
	  for (k = 0; k < n; k++) {
	    tmp = z[k];
	    f = tmp[i+1];
	    tmp[i+1] = s*tmp[i] + c*f;
	    tmp[i] = c*tmp[i] - s*f;
	  }
	}
	if (r == 0.0 && i >= l) continue;
	d[l] -= p;
	e[l] = g;
	e[m] = 0.0;
      }
    }
  }
}

/***=======================================================================***/
/*** QuatAlign: function for aligning sets of coordinates using a          ***/
/***            quaternion.  Assumes that the sets have been aligned in    ***/
/***            terms of geometric center.  Returns as U the rotation      ***/
/***            which minimizes the RMSD between the sets of coordinates.  ***/
/***=======================================================================***/
void QuatAlign(double* frameI, double* frameII, int num_atoms, double* mass,
	       int m, dmat *U)
{
  int i, itr, maxeigloc;
  double maxeig, a, x, y, z, totmass, tmass;
  double aa, ab, ac, ba, bb, bc, ca, cb, cc, i1, i2, i3, g1, g2, g3;
  double diag[4], sdiag[4];
  dmat R;

  /*** Compute the quaternion matrix ***/
  R = CreateDmat(4, 4);
  totmass = (m == 1) ? 1.0/DSum(mass, num_atoms) : 1.0/num_atoms;
  for (i = 0; i < num_atoms; i++) {
    itr = 3*i;
    i1 = frameII[itr];
    i2 = frameII[itr+1];
    i3 = frameII[itr+2];
    g1 = frameI[itr];
    g2 = frameI[itr+1];
    g3 = frameI[itr+2];
    aa = i1*g1;
    ab = i1*g2;
    ac = i1*g3;
    ba = i2*g1;
    bb = i2*g2;
    bc = i2*g3;
    ca = i3*g1;
    cb = i3*g2;
    cc = i3*g3;
    if (m == 1) {
      tmass = mass[i];
      R.data[0] += tmass*(aa+bb+cc);
      R.data[1] += tmass*(cb-bc);
      R.data[2] += tmass*(ac-ca);
      R.data[3] += tmass*(ba-ab);
      R.data[5] += tmass*(aa-bb-cc);
      R.data[6] += tmass*(ab+ba);
      R.data[7] += tmass*(ca+ac);
      R.data[10] += tmass*(bb-cc-aa);
      R.data[11] += tmass*(bc+cb);
      R.data[15] += tmass*(cc-aa-bb);
    }
    else {
      R.data[0] += aa+bb+cc;
      R.data[1] += cb-bc;
      R.data[2] += ac-ca;
      R.data[3] += ba-ab;
      R.data[5] += aa-bb-cc;
      R.data[6] += ab+ba;
      R.data[7] += ca+ac;
      R.data[10] += bb-cc-aa;
      R.data[11] += bc+cb;
      R.data[15] += cc-aa-bb;
    }
  }
  R.data[4] = R.data[1];
  R.data[8] = R.data[2];
  R.data[12] = R.data[3];
  R.data[9] = R.data[6];
  R.data[13] = R.data[7];
  R.data[14] = R.data[11];
  for (i = 0; i < 16; i++) {
    R.data[i] *= totmass;
  }
  TRED2(R.map, 4, diag, sdiag);
  TQLI(diag, sdiag, 4, R.map);

  maxeig = diag[0];
  maxeigloc = 0;
  for (i = 1; i < 4; i++) {
    if (diag[i] > maxeig) {
      maxeig = diag[i];
      maxeigloc = i;
    }
  }
  a = R.data[maxeigloc];
  x = R.data[maxeigloc+4];
  y = R.data[maxeigloc+8];
  z = R.data[maxeigloc+12];

  /*** Construct the rotation matrix ***/
  U->data[0] = a*a + x*x -y*y - z*z;
  U->data[1] = 2.0*(x*y + a*z);
  U->data[2] = 2.0*(z*x - a*y);
  U->data[3] = 2.0*(x*y - a*z);
  U->data[4] = a*a - x*x + y*y - z*z;
  U->data[5] = 2.0*(y*z + a*x);
  U->data[6] = 2.0*(z*x + a*y);
  U->data[7] = 2.0*(y*z - a*x);
  U->data[8] = a*a - x*x - y*y + z*z;

  DestroyDmat(&R);
}

/***=======================================================================***/
/*** WeightRMSD: function for finding the weighted RMSD between molecules  ***/
/***             mola and molb with weights mass on each equivalent atom.  ***/
/***=======================================================================***/
double WeightRMSD(double* mola, double* molb, double* mass, int m, int natom)
{
  int i;
  double rms_dev, dx, dy, dz, tmp_mass, tot_mass;
  double* tmp;
  double* tm2p;

  rms_dev = 0.0;
  tot_mass = 0.0;
  for (i = 0; i < natom; i++) {
    tmp = &mola[3*i];
    tm2p = &molb[3*i];
    dx = tmp[0] - tm2p[0];
    dy = tmp[1] - tm2p[1];
    dz = tmp[2] - tm2p[2];
    tmp_mass = (m == 1) ? mass[i] : 1.0;
    tot_mass += tmp_mass;
    rms_dev += tmp_mass*SQ_DIST(dx, dy, dz);
  }
  rms_dev /= tot_mass;

  return sqrt(rms_dev);
}
