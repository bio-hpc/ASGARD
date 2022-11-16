/*!
   \file index.h
   \brief Contains indexing functions for pairs, triplets and quadruples
   \author Martin Peters

   $Date: 2010/08/19 11:24:06 $
   $Revision: 1.8 $

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

#ifndef INDEX_H
#define INDEX_H

#include "Diagnostics/MTKException.h"

#include <iostream>
#include <vector>
#include <algorithm>

namespace MTKpp
{

// ============================================================
// Function: index2dTo1d
// ------------------------------------------------------------
// Returns unique index of i and j
/*

  j -> 
    +-                   -+       +-                   -+
    | 1,1 .,. .,. .,. 1,m |       | 0                   |
 i  | 2,1 2,2 .,. .,. 2,m |       | 1  2                |
 |  | 3,1 .,. 3,3 .,. 3,m |  ==>  | 3  4  5             |
\./ | 4,1 .,. .,. 4,4 4,m |       | 6  7  8  9          |
    | n,. .,. .,. .,. n,m |       |                     |
    +-                   -+       +-                   -+
*/
// ============================================================
inline int index2dTo1d(int i, int j)
{
    if ((i < 1) || (j < 1)) {
      return -1;
    }

    int ki = std::max(i, j);
    int kj = std::min(i, j);
    int index = ( ki*(ki-1) )/2 + kj - 1;
    return index;
};

// ============================================================
// Function: indexAB
// ------------------------------------------------------------
// Returns unique index of a and b ( Bond )
// a -- b
// ============================================================
inline int indexAB(int a, int b, unsigned int Nb)
{
    if ((a < 1) || (b < 1)  || (Nb < 1)) {
      return -1;
    }

    if (b < a) {
      int btemp = a;
      a = b;
      b = btemp;
    }
    return ((a-1)*Nb + (b-1));
};

// ============================================================
// Function: indexABC
// ------------------------------------------------------------
// Returns unique index of a, b and c ( Angle )
/*
       b
      / \
     a   c
*/
// ============================================================
inline int indexABC(int a, int b, int c, unsigned int Nb, unsigned int Nc)
{
    if ((a < 1) || (b < 1) || (Nb < 1) || (Nc < 1)) {
      return -1;
    }

    if (c < a) {
      int ctemp = a;
      a = c;
      c = ctemp;
    }
    return ((a-1)*(Nb*Nc) + (b-1)*Nc + (c-1));
};

// ============================================================
// Function: indexABC_LL
// ------------------------------------------------------------
// Returns unique index of a, b and c ( Angle )
/*
       b
      / \
     a   c
*/
// ============================================================
inline ULONG_KIND indexABC_ULL(int a, int b, int c,
                               unsigned int Nb, unsigned int Nc)
{
    if ((a < 1) || (b < 1) || (Nb < 1) || (Nc < 1)) {
      throw MTKException(" Error in indexABC_ULL ... exiting ");
    }

    if (c < a) {
      int ctemp = a;
      a = c;
      c = ctemp;
    }

    ULONG_KIND i = ((a-1)*(Nb*Nc) + (b-1)*Nc + (c-1));

    return i;
};

// ============================================================
// Function: indexABCD
// ------------------------------------------------------------
// Returns unique index of a, b, c and d ( Torsion )
/*     b -- c
      /      \
     a        d
*/
// ============================================================
inline int indexABCD(int a, int b, int c, int d,
           unsigned int Nb, unsigned int Nc, unsigned int Nd)
{
    if ((a < 1) || (b < 1) || (Nb < 1) || (Nc < 1) || (Nd < 1)) {
      return -1;
    }

    if (d < a) {
      int dtemp = a;
      a = d;
      d = dtemp;
      int ctemp = b;
      b = c;
      c = ctemp;
    }
    int i = ((a-1)*(Nb*Nc*Nd) + (b-1)*(Nc*Nd) + (c-1)*Nd + (d-1));
    return i;
};

// ============================================================
// Function: indexABCD_ULL
// ------------------------------------------------------------
// Returns unique index of a, b, c and d ( Torsion )
/*     b -- c
      /      \
     a        d
*/
// ============================================================
inline ULONG_KIND indexABCD_ULL(int a, int b, int c, int d,
                          unsigned int Nb, unsigned int Nc, unsigned int Nd)
{
    if ((a < 1) || (b < 1) || (Nb < 1) || (Nc < 1) || (Nd < 1)) {
      throw MTKException(" Error in indexABCD_ULL ... exiting ");
    }

    if (d < a) {
      int dtemp = a;
      a = d;
      d = dtemp;
      int ctemp = b;
      b = c;
      c = ctemp;
    }

    ULONG_KIND i = ((a-1)*(Nb*Nc*Nd) + (b-1)*(Nc*Nd) + (c-1)*Nd + (d-1));
    return i;
};

// ============================================================
// Function: indexABCD
// ------------------------------------------------------------
// Returns unique index of a, b, c and d ( Torsion )
/*
             e
           /
     a---b
          \
           d
*/
// ============================================================
inline int indexABCDb(int a, int b, int c, int d,
           unsigned int Nb, unsigned int Nc, unsigned int Nd)
{
    if ((a < 1) || (b < 1) || (Nb < 1) || (Nc < 1) || (Nd < 1)) {
      return -1;
    }
    std::vector<int> tempVec;
    tempVec.push_back(a);
    tempVec.push_back(b);
    tempVec.push_back(c);
    tempVec.push_back(d);
    std::sort( tempVec.begin(), tempVec.end() );

    return ((tempVec[0]-1)*(Nb*Nc*Nd) + (tempVec[1]-1)*(Nc*Nd) + (tempVec[2]-1)*Nd + (tempVec[3]-1));
};

} // MTKpp namespace

#endif // INDEX_H

