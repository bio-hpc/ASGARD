#ifndef CommandHeadings
#define CommandHeadings

#include "pmeRecip.h"
#include "TopologyDS.h"
#include "TrajectoryDS.h"
#include "CrdManipDS.h"
#include "ChargeFitDS.h"
#include "ParamFitDS.h"
#include "IPolQDS.h"

void ReadCommFile(dircon *dcinp, reccon *rcinp, prmtop* tp, trajcon *tj,
                  fset *myfit, prmset *myparms, ipqcon *ipqinp, char* source);

void SetGridDims(reccon *rcinp, coord *crd);

#endif
