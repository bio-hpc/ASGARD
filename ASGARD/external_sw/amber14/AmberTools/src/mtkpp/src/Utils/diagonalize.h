/*! 
   \file diagonalize.h
   \brief Inferface function to boost/lapack syev function
   \author Martin Peters

   $Date: 2010/03/29 20:33:22 $
   $Revision: 1.7 $

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

#ifndef DIAGONALIZE_H
#define DIAGONALIZE_H

#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "boost/numeric/bindings/lapack/syev.hpp"
#include "boost/numeric/bindings/lapack/gesvd.hpp"
#include "boost/numeric/bindings/lapack/gesdd.hpp"
#include "boost/numeric/bindings/traits/ublas_matrix.hpp"
#include "boost/numeric/bindings/traits/ublas_vector.hpp"

namespace ublas  = boost::numeric::ublas;
namespace lapack = boost::numeric::bindings::lapack;

namespace MTKpp
{

/*!
   \brief Diagonalize a matrix

   Diag(S):
    Diagonalizing S involves finding U such that:
    transpose(U) * S * U = D
    where D is the diagonal matrix containig the eigenvalues,
    and U is the eigenvector matrix.

    \param eigenvectors Matrix to be diagonalized, it is destroyed, returned containing the eigenvectors
    \param eigenvalues Vector of resulting eigenvalues
    \return success
*/

inline int diagonalize(ublas::matrix<double, ublas::column_major>& eigenvectors,
                       ublas::vector<double>& eigenvalues) {
    int r = lapack::syev( 'V', 'U', eigenvectors, eigenvalues, lapack::minimal_workspace() );
    return r;
}

/*!
   \brief Singular Value Decomposition (SVD)

    From boost documentation
      gesvd --> simple driver
      gesdd --> divide and conquer driver

     gesvd/gesdd computes the singular value decomposition (SVD) of
     M-by-N matrix A, optionally computing the left and/or right
     singular vectors. The SVD is written

          A = U * S * V^T    or    A = U * S * V^H

     where S is an M-by-N matrix which is zero except for its min(m,n)
     diagonal elements, U is an M-by-M orthogonal/unitary matrix, and V
     is an N-by-N orthogonal/unitary matrix. The diagonal elements of S
     are the singular values of A; they are real and non-negative, and
     are returned in descending  order. The first min(m,n) columns of
     U and V are the left and right singular vectors of A. (Note that
     the routine returns V^T or V^H, not V.

     A_nxp = U_nxn S_nxp V^T_pxp

    Other links
    http://en.wikipedia.org/wiki/Singular_value_decomposition
    http://www.netlib.org/lapack/lug/node32.html

    \return success
*/
inline int svd(ublas::matrix<double, ublas::column_major>& H,
               ublas::vector<double>& s,
               ublas::matrix<double, ublas::column_major>& U,
               ublas::matrix<double, ublas::column_major>& VT) {

    //int r = lapack::gesvd('A', 'A', H, s, U, VT);
    int r = lapack::gesdd(H, s, U, VT);
    return r;
}

/*!
   \brief Sorts eigenvalues and applies this ordering on the eigenvector matrix
   \param eigenvectors Eigenvector matrix
   \param eigenvalues Eigenvalue matrix
   \param order Ascending = 0, Desending = 1
*/
inline void eigenValueSort(ublas::matrix<double, ublas::column_major>& eigenvectors,
                           ublas::vector<double>& eigenvalues, int order) {
    int k;
    double p;
    int size = eigenvectors.size1();

    // Ascending
    if (!order) {
      for (int i = 0; i < size-1 ; ++i) {
        k = i;
        p = eigenvalues(i);

        for (int j = i+1; j < size; ++j) {
          if (eigenvalues(j) < p) {
            k = j;
            p = eigenvalues(j);
          }
        }
        if ( k != i ) {
          eigenvalues(k) = eigenvalues(i);
          eigenvalues(i) = p;
          for (int m = 0; m < size; ++m) {
            p = eigenvectors(m,i);
            eigenvectors(m,i) = eigenvectors(m,k);
            eigenvectors(m,k) = p;
          }
        }
      }
    }
    // Descending
    else {
      for (int i = 0; i < size-1 ; ++i) {
        k = i;
        p = eigenvalues(i);

        for (int j = i+1; j < size; ++j) {
          if (eigenvalues(j) > p) {
            k = j;
            p = eigenvalues(j);
          }
        }
        if ( k != i ) {
          eigenvalues(k) = eigenvalues(i);
          eigenvalues(i) = p;
          for (int m = 0; m < size; ++m) {
            p = eigenvectors(m,i);
            eigenvectors(m,i) = eigenvectors(m,k);
            eigenvectors(m,k) = p;
          }
        }
      }
    }
}

} // MTKpp namespace

#endif // DIAGONALIZE_H
