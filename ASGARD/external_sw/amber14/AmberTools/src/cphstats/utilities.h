#ifndef UTILITIES_H
#define UTILITIES_H

#include <string>
#include "cpout.h"
#include "types.h"

bool fexists(std::string const&);
int sort_remd_files(CpoutList, std::string const&, const bool);

#endif /* UTILITIES_H */
