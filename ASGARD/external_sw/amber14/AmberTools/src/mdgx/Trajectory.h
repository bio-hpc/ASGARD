#ifndef TrajectoryHeadings
#define TrajectoryHeadings

#include "AmberNetcdf.h"

#include "CrdManipDS.h"
#include "CellManipDS.h"
#include "TopologyDS.h"
#include "TrajectoryDS.h"
#include "pmeDirectDS.h"
#include "pmeRecipDS.h"
#include "CompFrcDS.h"
#include "TimingsDS.h"

void UpdateStepNumber(trajcon *tj, int newstep);

coord ReadRst(prmtop *tp, char* source);

void SpliceFileName(trajcon *tj, char* base, char* suff, char* fname,
                    int dprc);

void InitializeEnergy(Energy *sysUV, trajcon *tj, prmtop *tp,
		      int allocBondUdc);

void DestroyEnergyTracker(Energy *sysUV);

void ExtendCoordinates(coord *tc, prmtop *tp);

coord InitCoords(prmtop *tp, trajcon *tj, int n);

void WriteDiagnostics(trajcon *tj, prmtop *tp, dircon *dcinp, reccon *rcinp,
                      FrcTab *Etab, cellgrid *CG, coord *crd, Energy *sysUV,
		      execon *etimers, int n);

void WriteRst(cellgrid *CG, coord *tc, prmtop *tp, trajcon *tj, int n);

void WriteCrd(cellgrid *CG, coord *crd, int cvf, trajcon *tj, prmtop *tp,
	      int n);

void WriteCDF(cellgrid *CG, coord *crd, int cvf, trajcon *tj, cdftrj *Acdf,
	      prmtop *tp, int n);

void DestroyTrajCon(trajcon *tj);

void SynchronizeCoordinates(cellgrid* CG, trajcon *tj);

#endif
