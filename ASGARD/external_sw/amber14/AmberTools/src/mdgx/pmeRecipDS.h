#ifndef pmeRecipStructs
#define pmeRecipStructs

#ifndef PREP_API
#include "fftw3.h"

#include "GridDS.h"
#include "mleRecipDS.h"
#endif

struct pmeRecipControlData {

  /*** Smooth particle Mesh Ewald ***/
  int ordr[3];
  int* ng;
  double S;

  /*** Multi-Level (Smooth Particle Mesh) Ewald ***/
  int nlev;
  int nslab;
  int nstrip;
  int ggordr;
  int PadYZ[4];
  double cfac[4];
  g2gmap* SPrv;
  g2gmap* SPcv;
  dbook* Urec;
  dbook* QL;
  fftw_plan* forwplan;
  fftw_plan* backplan;
};
typedef struct pmeRecipControlData reccon;

struct BCMeshKit {
  int plans;             // Flag to indicate that FFT plans have been made
                         //   especially for this struct (the plans may be
                         //   borrowed from elsewhere)
  double SelfEcorr;      // Self energy correction for this mesh setup
  double* Bx;            // B mesh prefactors in X, Y, and Z
  double* By;            //
  double* Bz;            //
  double* mx;            // M value (reciprocal space vector index) in X, Y,
  double* my;            //   and Z
  double* mz;            //
  double* mxs;           // Shifted M values in X, Y, or Z
  double* mys;           //
  double* mzs;           //
  double* Ex;            // Exponential tables for X, Y, and Z (used for both
  double* Ey;            //   orthorhombic and non-orthorhombic unit cells)
  double* Ez;            //
  fftw_plan forwplan;    // Forward FFT plan
  fftw_plan backplan;    // Backward FFT plan
};
typedef struct BCMeshKit bckit;

#endif
