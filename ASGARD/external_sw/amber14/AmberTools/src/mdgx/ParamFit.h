#ifndef ParamFitHeadings
#define ParamFitHeadings

#include "ParamFitDS.h"
#include "TopologyDS.h"
#include "TrajectoryDS.h"

int str4cmp(char* A, char* B);

int TypeCompare(char* T1a, char* T1b, char* T1c, char* T1d, char* T2a,
                char* T2b, char* T2c, char* T2d, int order, int strict);

void FitParams(prmset *mp, trajcon *tj);

#endif
