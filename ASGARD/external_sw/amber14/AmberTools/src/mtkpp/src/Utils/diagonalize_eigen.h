/*! 
   \file diagonalize_eigen.h
   \brief Inferface functions to eigen library
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

#ifndef DIAGONALIZEEIGEN_H
#define DIAGONALIZEEIGEN_H

// Eigen
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
using namespace Eigen;

namespace MTKpp
{

/*!
   \brief Sorts eigenvalues and applies this ordering on the eigenvector matrix
   \param eigenvectors Eigenvector matrix
   \param eigenvalues Eigenvalue matrix
   \param order Ascending = 0, Desending = 1
*/
inline void eigenValueSort(MatrixXd& eigenvectors,
                           VectorXd& eigenvalues, int order) {
    int k;
    double p;
//    int size = eigenvectors.size1();
    int size = eigenvectors.rows();

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

#endif // DIAGONALIZEEIGEN_H
