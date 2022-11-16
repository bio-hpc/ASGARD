#ifndef STRING_MANIP_H
#define STRING_MANIP_H
// Routines for string manipulation

#include <string>
#include <vector>
#include "exceptions.h"

// Various string splitting routines
std::vector<std::string> split(char*, const char*);
std::vector<std::string> split(std::string const&, const char*);
std::vector<std::string> split(std::string const&, std::string const&);
std::vector<std::string> split(std::string const&);
std::vector<std::string> split(char*);

// White-space stripping routines (only strips from either end)
std::string strip(std::string const&);
std::string strip(char*);

// Casts to upper-case
std::string upper(std::string const&);
std::string upper(const char*);

// Converts a string to an integer, double, or float
int StringToInt(std::string const&);
double StringToDouble(std::string const&);
float StringToFloat(std::string const&);

#endif
