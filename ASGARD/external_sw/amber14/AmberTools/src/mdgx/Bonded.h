#ifndef BondedHeadings
#define BondedHeadings

#include "CrdManipDS.h"
#include "TopologyDS.h"
#include "CompFrcDS.h"
#include "TrajectoryDS.h"
#include "CellManipDS.h"

void CellBondedIntr(cell *C, cellgrid *CG, coord *crd, prmtop *tp,
		    FrcTab *Cfrc, FrcTab *Hfrc, Energy *sysUV);

void AttenuatePairVir(prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc, int atmA,
		      int atmB, double dx, double dy, double dz, double fmag,
		      double* afrc, double* bfrc, Energy *sysUV,
		      double elec14fac, double lj14fac);

void AttenuatePairFrcNrg(prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc, int atmA,
                         int atmB, double dx, double dy, double dz,
			 double fmag, double* afrc, double* bfrc,
			 Energy *sysUV, double elec14fac, double lj14fac);

void AttenuatePairFrc(prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc, int atmA,
		      int atmB, double dx, double dy, double dz, double fmag,
		      double* afrc, double* bfrc, double elec14fac,
		      double lj14fac);

void AttenuatePairNrg(prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc, int atmA,
		      int atmB, double dx, double dy, double dz,
		      Energy *sysUV, double elec14fac, double lj14fac);

void BondVir(double* aptr, double *bptr, double* afrc, double* bfrc,
             bondcomm *bcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
	     Energy *sysUV);

void BondFrcNrg(double* aptr, double *bptr, double* afrc, double* bfrc,
                bondcomm *bcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
		Energy *sysUV);

void BondFrc(double* aptr, double *bptr, double* afrc, double* bfrc,
             bondcomm *bcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc);

void BondNrg(double* aptr, double *bptr, bondcomm *bcom, prmtop *tp,
             FrcTab *Cfrc, FrcTab *Hfrc, Energy *sysUV);

void AnglVir(double *aptr, double *bptr, double *ctpr, double *afrc,
             double *bfrc, double *cfrc, anglcomm *acom, prmtop *tp,
	     FrcTab *Cfrc, FrcTab *Hfrc, Energy *sysUV);

void AnglFrcNrg(double *aptr, double *bptr, double *ctpr, double *afrc,
		double *bfrc, double *cfrc, anglcomm *acom, prmtop *tp,
		FrcTab *Cfrc, FrcTab *Hfrc, Energy *sysUV);

void AnglFrc(double *aptr, double *bptr, double *ctpr, double *afrc,
	     double *bfrc, double *cfrc, anglcomm *acom, prmtop *tp,
	     FrcTab *Cfrc, FrcTab *Hfrc);

void AnglNrg(double *aptr, double *bptr, double *ctpr, anglcomm *acom,
	     prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc, Energy *sysUV);

void DiheVir(double *aptr, double *bptr, double *cptr, double *dptr,
             double *afrc, double *bfrc, double *cfrc, double *dfrc,
             dihecomm *hcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
	     Energy *sysUV);

void DiheFrcNrg(double *aptr, double *bptr, double *cptr, double *dptr,
		double *afrc, double *bfrc, double *cfrc, double *dfrc,
		dihecomm *hcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
		Energy *sysUV);

void DiheFrc(double *aptr, double *bptr, double *cptr, double *dptr,
	     double *afrc, double *bfrc, double *cfrc, double *dfrc,
	     dihecomm *hcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc);

void DiheNrg(double *aptr, double *bptr, double *cptr, double *dptr,
	     dihecomm *hcom, prmtop *tp, FrcTab *Cfrc, FrcTab *Hfrc,
	     Energy *sysUV);

#endif
