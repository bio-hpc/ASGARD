/*! 
   \file lbfgs.h
   \brief This class contains code for the limited-memory Broyden-Fletcher-Goldfarb-Shanno (LBFGS) algorithm
    for large-scale multidimensional unconstrained minimization problems. This file is a translation of 
    Fortran code written by Jorge Nocedal. It is distributed as part of the RISO project.

    This code is derived from the Fortran program lbfgs.f.
    The C++ translation was effected mostly mechanically, with some manual clean-up.

    Here's some information on the original LBFGS Fortran source code, 
    available at http://www.netlib.org/opt/lbfgs_um.shar.
    This info is taken verbatim from the Netlib blurb on the Fortran source.

    LBFGS minimizer from the original lbfgs.f, converted to c++ by Kenneth Ayers

   $Date: 2009/04/08 10:58:54 $
   $Revision: 1.4 $

   ----------------------------------------------------------------------------
*/

#ifndef LBFGS_H
#define LBFGS_H

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

#endif

namespace MTKpp
{

void lbfgs_(int*, int*, double*, double*, double*, int*, double*, int*, double*, double*, double*, int*);

void mcsrch_(int*, double*, double*, double*, double*, double*, double*, double*, int*, int*, int*, double*);

void mcstep_(double*, double*, double*, double*, double*, double *, double*, double*, double*, int*, double*, double*, int*);

double ddot_(int*, double*, int*, double*, int*);

void daxpy_(int*, double*, double*, int*, double*, int*);

void lb1_(int*, int*, int*, double*, int*, int*, double*, double*, double*, double*, int*);

}

#endif // LBFGS_H
