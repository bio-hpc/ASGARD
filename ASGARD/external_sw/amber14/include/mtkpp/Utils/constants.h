/*!
   \file constants.h
   \brief Container for constants
   \author Martin Peters

   $Date: 2010/04/29 19:38:22 $
   $Revision: 1.15 $

   ----------------------------------------------------------------------------

   MTK++ - C++ package of modeling libraries.

   Copyright (C) 2005-2006  (see AUTHORS file for a list of contributors)

   This file is part of MTK++.

   MTK++ is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   MTK++ is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lessser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   ----------------------------------------------------------------------------
*/

#ifndef CONSTANTS_H
#define CONSTANTS_H

#ifdef __INTEL_COMPILER

// remark #177: variable was declared but never referenced
#pragma warning(disable:177)

// remark #181: argument is incompatible with corresponding format string conversion
#pragma warning(disable:181)

// remark #304: access control not specified ("public" by default)
#pragma warning(disable:304)

// remark #383: value copied to temporary, reference to temporary used
#pragma warning(disable:383)

// remark #424: extra ";" ignored
#pragma warning(disable:424)

// remark #593: variable was set but never used
#pragma warning(disable:593)

// remark #810: conversion from "double" to "int" may lose significant bits
#pragma warning(disable:810)

// remark #869: parameter was never referenced
#pragma warning(disable:869)

// remark #981: operands are evaluated in unspecified order
#pragma warning(disable:981)

// warning #1125: virtual function override intended?
#pragma warning(disable:1125)

// remark #1418: external function definition with no prior declaration
#pragma warning(disable:1418)

// remark #1572: floating-point equality and inequality comparisons are unreliable
// disabled -> everyone knows it, the parser passes this problem
//             deliberately to the user
#pragma warning(disable:1572)

// remark #1599: declaration hides variable "t"
#pragma warning(disable:1599)

// remark #2259: non-pointer conversion from "double" to "int" may lose significant bits
#pragma warning(disable:2259)

#endif

#include <math.h>
namespace MTKpp
{
    /////////////////////////////////////////////
    //            MTK++ defines                //
    /////////////////////////////////////////////

// C99 Compliant compiler
#ifdef C99_OK
#define ULONG_KIND unsigned long long
#define LONG_KIND long long
#define ULONGLONG_MAX 0xffffffffffffffffLLU
#define LONGLONG_MAX 0x7fffffffffffffffLL
#define LONGLONG_MIN (-LONGLONG_MAX - 1LL)
#endif

#ifndef C99_OK
#define ULONG_KIND unsigned long
#define LONG_KIND long
#endif

    /////////////////////////////////////////////
    //           Constant numbers              //
    /////////////////////////////////////////////

    //! Maximum number of atoms
#ifdef C99_OK
    const unsigned int MAXATOMS = 100000;
#endif

#ifndef C99_OK
    const unsigned int MAXATOMS = 1000;
#endif

    //! Very Large Number
    const double BIGNUM = 1.e12;

    //! Very Large Number
    const double APPROXZERO = 1.e-6;

    //! PI
    const double PI = 3.14159265359;

    //! (PI)^2
    const double PIt2 = PI * 2;

    //! (PI)^3/2
    //const double PIto3over2 = pow(PI, 1.5);

    //! 1.0/(PI)
    const double INVPI = 1.0/PI;

    //! Radian to Degree
    const double RAD2DEG = 180.0/PI;

    //! Degree to Radian
    const double DEG2RAD = 1.0/RAD2DEG;

    //! Angstrom to Bohr
    const double ANG2BOHR = 0.529177249;

    //! Bohr to Angstrom
    const double BOHR2ANG = 1.0/ANG2BOHR;

    //! Hartree to kcal/mol
    const double H2KCALMOL = 627.5095;

    //! electron to kcal
    const double E2KCAL = 18.2223;

    /*!
       \brief Hartree/Bohr^2 to kcal/mol A^2

       \verbatim
         H      627.5095 kcal/mol           1 B               1 B
        --- . -------------------- .  -------------- .  ---------------- = 2240.87 kcal/mol A^2
        B^2         1 H                0.5291772 A      0.5291772 A
       \endverbatim
    */
    const double HB2TOKCALMOLA2 = H2KCALMOL*BOHR2ANG*BOHR2ANG;

    /*!
       \brief definition of a bond, covRadius1 + covRadius + tolerance
        Bond Types and Hybridization, J. Comp. Chem. 12, 891-898, 1991
    */
    const double BONDTOLERANCE = 0.4;

    //! sp2 versus sp3 angle cutoff
    const double SP2ANGLE = 115.0;

    //! sp versus sp2 angle cutoff
    const double SPANGLE = 160.0;

    //! Angle below which the type of an atom with 2 heavy bonded atom should be reconsidered
    const double AMBIGUOUSANGLE = 122.0;

    //! Generic message type
    const int MESSAGE = 0;

    //! Error message tyoe
    //const int ERROR = 1;
    const int MTK_ERROR = 1;

    //! Warning message type
    const int WARNING = 2;

    //! Debug message type
    const int mDEBUG = 3;

    //! Info message type
    const int INFO = 4;

} // MTKpp namespace

#endif // CONSTANTS_H

