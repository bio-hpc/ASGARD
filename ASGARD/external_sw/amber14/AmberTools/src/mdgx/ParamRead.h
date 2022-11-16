#ifndef ParamReadFuncs
#define ParamReadFuncs

#include "ParamFitDS.h"
#include "TrajectoryDS.h"

int CrossRefAtomType(prmset *mp, char* aname);

void ReadParmFile(prmset *mp, trajcon *tj);

void ReadFrcmodFile(prmset *mp, trajcon *tj);

void RecastAtomTypes(prmset *mp);

void CleaveAtomTypes(prmset *mp);

#endif
