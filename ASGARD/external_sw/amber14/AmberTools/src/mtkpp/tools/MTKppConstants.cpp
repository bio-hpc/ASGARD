/*!
   \file MTKppConstants.cpp

   \brief Program for printing the definitions and constants used in MTK++

   \author Martin B. Peters

   $Date: 2010/04/22 22:19:54 $
   $Revision: 1.7 $

   ----------------------------------------------------------------------------

   MTK++ - C++ package of modeling libraries.

   Copyright (C) 2005-2006  (see AUTHORS file for a list of contributors)

   ----------------------------------------------------------------------------
*/
#include "Utils/constants.h"
#include "Utils/printHeader.h"

// - COMMAND LINE OPTIONS
#include "Parsers/commLineOptions.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <vector>
#include <map>
#include <algorithm>

using namespace MTKpp;
/*!
   \brief Program for printing the definitions and constants used in MTK++
*/
int main (int argc, char **argv)
{
    std::string prog_name = "MTKppConstants";
    std::vector<std::string> authors;
    std::string author = "Martin B. Peters";
    authors.push_back(author);
    printHeader(std::cout, prog_name, authors);

    std::cout << "\n\n Number Types used in MTK++ " << std::endl;
    std::cout << "Size of int types is "
              << sizeof(int) << " bytes"
              << '\n';
    std::cout << "Signed int min: "
              << INT_MIN << " max: "
              << INT_MAX << '\n';
    std::cout << "Unsigned int min: 0 max: "
              << UINT_MAX << "\n\n";

#ifdef C99_OK
    std::cout << "Size of long long int type is "
              << sizeof(long long) << " bytes"
              << '\n';
    std::cout << "Signed long long int min: "
              << LONGLONG_MIN << " max: "
              << LONGLONG_MAX << '\n';
    std::cout << "Size of ulong long int type is "
              << sizeof(unsigned long long) << " bytes"
              << '\n';
    std::cout << "Unsigned long long int min: 0 max: "
              << ULONGLONG_MAX << "\n\n";
#endif

#ifndef C99_OK
    std::cout << "Size of long int type is "
              << sizeof(long) << " bytes"
              << '\n';
    std::cout << "Signed long int min: "
              << LONG_MIN << " max: "
              << LONG_MAX << '\n';
    std::cout << "Size of ulong int type is "
              << sizeof(unsigned long) << " bytes"
              << '\n';
    std::cout << "Unsigned long int min: 0 max: "
              << ULONG_MAX << "\n\n";
#endif

    std::cout << "Size of double types is "
              << sizeof(double) << " bytes \n\n";

    std::cout << " Constants used in MTK++ " << std::endl;
    std::cout << "    MAX ATOMS ALLOWED   = " << MAXATOMS  << std::endl;
    std::cout << "    PI                  = " << PI         << std::endl;
    std::cout << "    (PI)2               = " << PIt2       << std::endl;
    //std::cout << "    (PI)3/2             = " << PIto3over2 << std::endl;
    std::cout << "    1/(PI)              = " << INVPI      << std::endl;
    std::cout << "    RADIAN 2 DEGREE     = " << RAD2DEG    << std::endl;
    std::cout << "    DEGREE 2 RADIAN     = " << DEG2RAD    << std::endl;
    std::cout << "    ANGSTROM 2 BOHR     = " << ANG2BOHR   << std::endl;
    std::cout << "    BOHR 2 ANGSTROM     = " << BOHR2ANG   << std::endl;
    std::cout << "    HARTREE 2 KCAL/MOL  = " << H2KCALMOL      << std::endl;
    std::cout << "    ELECTRON 2 KCAL     = " << E2KCAL         << std::endl;
    std::cout << "    HB^2 TO KCAL/MOLA^2 = " << HB2TOKCALMOLA2 << std::endl;

    return 0;
}

