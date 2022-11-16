#ifndef VirtualSiteHeadings
#define VirtualSiteHeadings

#include "TopologyDS.h"
#include "CrdManipDS.h"
#include "CellManipDS.h"
#include "TrajectoryDS.h"

void ExpandLJTables(prmtop *tp, double sig, double eps);

void ReadEPRuleFile(prmtop *tp);

void Connect1234(prmtop *tp);

void AddNewElimPair(prmtop *tp, int atmH, int atmX, int atmY, int ncon);

void AllocateElimList(prmtop *tp);

void DetermineEPFrames(prmtop *tp);

void XptLocator(double *aloc, double *agloc, double *bloc, double *cloc,
		double *dloc, double *eploc, double *epgloc, expt *tmr);

#ifdef MPI
void ExtraPointLocations(prmtop *tp, coord *crd, cellgrid *CG, trajcon *tj);
#else
void ExtraPointLocations(prmtop *tp, coord *crd, cellgrid *CG);
#endif

void CellXferEPForces(cell *C, prmtop *tp, coord *crd, cellgrid *CG);

void TransferEPForces(prmtop *tp, coord *crd);

#endif
