/*!
   \file BaseStats.cpp
   \brief Base class for statistical routines
   \author Martin Peters

   $Date: 2007/11/28 09:22:17 $
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

#include "BaseStats.h"

namespace MTKpp
{

// ============================================================
// Class : BaseStats()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
BaseStats::BaseStats() {}

// ============================================================
// Function : ~BaseStats()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
//BaseStats::~BaseStats() {}

// ============================================================
// Function : meanColumn(double)
// ------------------------------------------------------------
/* Calculates the Mean of a vector of doubles.

          --
       1  \
mean = -  /  y[i]  == <y>
       N  --

where N is the number of observations or the simply the length of y
(NANs are not allowed)
*/
// ============================================================
double BaseStats::meanColumn(ublas::matrix<double>& m, const int& c)
{
    double mean = 0.0;
    if (m.size1() <= 1) {
      return m(0,c);
    }

    for (unsigned int j = 0; j < m.size2(); j++) {
      if (static_cast<int>(j) == c) {
        for (unsigned int i = 0; i < m.size1(); i++) {
          mean += m(i,j);
        }
      }
    }
    mean /= m.size1();
    return mean;
}

// ============================================================
// Function : meanRow(double)
// ------------------------------------------------------------
// Calculates the Mean of a vector of doubles.
// ============================================================
double BaseStats::meanRow(ublas::matrix<double>& m, const int& r)
{
    std::cout << " TEST THIS CODE --> BaseStats::meanRow " << std::endl;
    double mean = 0.0;
    if (m.size2() <= 1) {
      return m(0,r);
    }

    for (unsigned int i = 0; i < m.size1(); i++) {
      if (static_cast<int>(i) == r) {
        for (unsigned int j = 0; j < m.size2(); j++) {
          mean += m(i,j);
        }
      }
    }
    mean /= m.size2();
    return mean;
}

// ============================================================
// Function : sumColumn(double)
// ------------------------------------------------------------
/* Calculates the Sum of a vector of doubles.

        --
        \
 sum =  /  y[i]
        --

*/
// ============================================================
double BaseStats::sumColumn(ublas::matrix<double>& m, const int& c)
{
    double sum = 0.0;
    if (m.size1() <= 1) {
      return m(0,c);
    }

    for (unsigned int j = 0; j < m.size2(); j++) {
      if (static_cast<int>(j) == c) {
        for (unsigned int i = 0; i < m.size1(); i++) {
          sum += m(i,j);
        }
      }
    }
    return sum;
}

// ============================================================
// Function : sumRow(double)
// ------------------------------------------------------------
// Calculates the Sum of a vector of doubles.
// ============================================================
double BaseStats::sumRow(ublas::matrix<double>& m, const int& r)
{
    double sum = 0.0;
    if (m.size2() <= 1) {
      return m(0,r);
    }

    for (unsigned int i = 0; i < m.size1(); i++) {
      if (static_cast<int>(i) == r) {
        for (unsigned int j = 0; j < m.size2(); j++) {
          sum += m(i,j);
        }
      }
    }
    return sum;
}

// ============================================================
// Function : maxColumn(double)
// ------------------------------------------------------------
// Calculates the max value of a vector of doubles.
// ============================================================
double BaseStats::maxColumn(ublas::matrix<double>& m, const int& c)
{
    if (m.size1() <= 1) {
      return 0.0;
    }

    double maxC = m(0,c);
    for (unsigned int j = 0; j < m.size2(); j++) {
      if (static_cast<int>(j) == c) {
        for (unsigned int i = 0; i < m.size1(); i++) {
          if (m(i,j) > maxC) {
            maxC = m(i,j);
          }
        }
      }
    }
    return maxC;
}

// ============================================================
// Function : maxColumn(double)
// ------------------------------------------------------------
// Calculates the max value of a vector of doubles.
// ============================================================
double BaseStats::maxColumn(ublas::matrix<double>& m, const int& c, int& r)
{
    if (m.size1() <= 1) {
      return 0.0;
    }

    r = 0;
    double maxC = m(0,c);
    for (unsigned int j = 0; j < m.size2(); j++) {
      if (static_cast<int>(j) == c) {
        for (unsigned int i = 0; i < m.size1(); i++) {
          if (m(i,j) > maxC) {
            maxC = m(i,j);
            r = i;
          }
        }
      }
    }
    return maxC;
}

// ============================================================
// Function : maxRow(double)
// ------------------------------------------------------------
// Calculates the max value of a vector of doubles.
// ============================================================
double BaseStats::maxRow(ublas::matrix<double>& m, const int& r)
{
    if (m.size2() <= 1) {
      return 0.0;
    }

    double maxR = m(r,0);

    for (unsigned int i = 0; i < m.size1(); i++) {
      if (static_cast<int>(i) == r) {
        for (unsigned int j = 0; j < m.size2(); j++) {
          if (m(i,j) > maxR) {
            maxR = m(i,j);
          }
        }
      }
    }
    return maxR;
}

// ============================================================
// Function : minColumn(double)
// ------------------------------------------------------------
// Calculates the min value of a vector of doubles.
// ============================================================
double BaseStats::minColumn(ublas::matrix<double>& m, const int& c)
{
    if (m.size1() <= 1) {
      return 0.0;
    }

    double minC = m(0,c);
    for (unsigned int j = 0; j < m.size2(); j++) {
      if (static_cast<int>(j) == c) {
        for (unsigned int i = 0; i < m.size1(); i++) {
          if (m(i,j) < minC) {
            minC = m(i,j);
          }
        }
      }
    }
    return minC;
}

// ============================================================
// Function : minColumn(double)
// ------------------------------------------------------------
// Calculates the min value of a vector of doubles.
// ============================================================
double BaseStats::minColumn(ublas::matrix<double>& m, const int& c, int& r)
{
    if (m.size1() < 1) {
      return 0.0;
    }

    if (m.size1() == 1) {
      return m(0,c);
    }

    r = 0;
    double minC = m(0,c);
    for (unsigned int j = 0; j < m.size2(); j++) {
      if (static_cast<int>(j) == c) {
        for (unsigned int i = 0; i < m.size1(); i++) {
          if (m(i,j) < minC) {
            minC = m(i,j);
            r = i;
          }
        }
      }
    }
    return minC;
}

// ============================================================
// Function : minRow(double)
// ------------------------------------------------------------
// Calculates the max value of a vector of doubles.
// ============================================================
double BaseStats::minRow(ublas::matrix<double>& m, const int& r)
{
    if (m.size2() <= 1) {
      return 0.0;
    }

    double minR = m(r,0);

    for (unsigned int i = 0; i < m.size1(); i++) {
      if (static_cast<int>(i) == r) {
        for (unsigned int j = 0; j < m.size2(); j++) {
          if (m(i,j) > minR) {
            minR = m(i,j);
          }
        }
      }
    }
    return minR;
}

// ============================================================
// Function : getColumnCenters(double)
// ------------------------------------------------------------
// Calculate the mean of each column
// ============================================================
int BaseStats::getColumnCenters(ublas::matrix<double>& mat, ublas::matrix<double>& mat_centers)
{
    if (mat.size2() == mat_centers.size2()) {
      for (unsigned int k = 0; k < mat.size2(); k++) {
        mat_centers(0,k) = meanColumn(mat, k);
      }
    }
    else {
      return 1;
    }
    return 0;
}

// ============================================================
// Function : centerColumns(double)
// ------------------------------------------------------------
// Centers the data: element - mean
// ============================================================
void BaseStats::centerColumns(ublas::matrix<double>& old_mat, ublas::matrix<double>& new_mat)
{
    double mean_array[old_mat.size2()];

    for (unsigned int k = 0; k < old_mat.size2(); k++) {
      mean_array[k] = meanColumn(old_mat,k);
    }

    int c = 0;
    for (unsigned int j = 0; j < old_mat.size2(); j++) {
      if (static_cast<int>(j) == c) {
        for (unsigned int i = 0; i < old_mat.size1(); i++) {
          new_mat(i,j) = old_mat(i,j) - mean_array[j];
        }
      }
      c += 1;
    }
}
/*
// ============================================================
// Function : centerRows(double)
// ------------------------------------------------------------
// Centers the data: element - mean
// ============================================================
void BaseStats::centerRows(Matrix &old_mat, Matrix &new_mat)
{
    double mean_array[old_mat.GetDimensionN()];

    for (int k = 0; k < old_mat.GetDimensionN(); k++) {
      mean_array[k] = meanRow(old_mat,k);
    }

    for (int i = 0; i < old_mat.GetDimensionN(); i++) {
      for (int j = 0; j < old_mat.GetDimensionM(); j++) {
         new_mat[i][j] = old_mat[i][j] - mean_array[i];
      }
    }
}
*/

// ============================================================
// Function : zScoreColumn(double)
// ------------------------------------------------------------
// Calculates the Z-Score of a Column of doubles.
// ============================================================
void BaseStats::zScoreColumns(ublas::matrix<double>& old_mat,
                              ublas::matrix<double>& new_mat)
{
    double stdDevColumn = 0.0;
    this->centerColumns(old_mat, new_mat);

    int c = 0;
    for (unsigned int j = 0; j < old_mat.size2(); j++) {
      stdDevColumn = this->standardDeviationColumn(old_mat, j);
      if (static_cast<int>(j) == c) {
        for (unsigned int i = 0; i < old_mat.size1(); i++) {
          new_mat(i,j) /= stdDevColumn;
        }
      }
      stdDevColumn = 0.0;
      c += 1;
    }
}

// ============================================================
// Function : autoScale(double)
// ------------------------------------------------------------
// Autoscale 
// ============================================================
int BaseStats::autoScale(ublas::matrix<double>& old_mat,
                          ublas::matrix<double>& new_mat,
                          ublas::matrix<double>& centers,
                          ublas::matrix<double>& stdDevs)
{
    if (old_mat.size1() != new_mat.size1() or
        old_mat.size2() != new_mat.size2()) return 1;

    if (centers.size2() != old_mat.size2() or
        stdDevs.size2() != old_mat.size2()) return 1;

    int c = 0;
    for (unsigned int j = 0; j < old_mat.size2(); j++) {
      if (static_cast<int>(j) == c) {
        for (unsigned int i = 0; i < old_mat.size1(); i++) {
          new_mat(i,j) = (old_mat(i,j) - centers(0,j)) / stdDevs(0,j);
        }
      }
      c += 1;
    }
    return 0;
}

/*
// ============================================================
// Function : zScoreRow(double)
// ------------------------------------------------------------
// Calculates the Z-Score of a Row of doubles.
// ============================================================
void BaseStats::zScoreRows(Matrix &old_mat, Matrix &new_mat)
{
    double stdDevRow = 0.0;
    centerRows(old_mat, new_mat);

    for (int i = 0; i < old_mat.GetDimensionN(); i++) {
      stdDevRow = standardDeviationRow(old_mat, i);
      for (int j =0; j < old_mat.GetDimensionM(); j++) {
        new_mat[i][j] /= stdDevRow;
      }
    }
}

// ============================================================
// Function : VarianceColumn(double)
// ------------------------------------------------------------
// Calculates the Sample Variance of a Column of doubles.
// This is a measure of spread about the mean
// ============================================================
double BaseStats::varianceColumn(Matrix &Vdoubles, const int& c)
{
    double variance = 0.0;

    int Vdouble_size = Vdoubles.GetDimensionN();
    if (Vdouble_size <= 1) {
      std::cout << "Vdoubles_size <= 1" << std::endl;
      return 0.0;
    }

    double mean = meanColumn(Vdoubles,c);

    for (int j = 0; j < Vdoubles.GetDimensionM(); j++) {
      for (int i = 0; i < Vdoubles.GetDimensionN(); i++) {
        if (j == c) {
          variance += ((Vdoubles[i][j] - mean) * (Vdoubles[i][j] - mean));
        }
      }
    }
    variance /= (Vdouble_size-1) ;
    return variance;
}
*/
// ============================================================
// Function : VarianceColumn(double)
// ------------------------------------------------------------
// Calculates the Sample Variance of a Column of doubles.
// This is a measure of spread about the mean
// ============================================================
double BaseStats::varianceColumn(ublas::matrix<double>& m, const int& c)
{
    double variance = 0.0;

    int nRows = m.size1();
    if (nRows <= 1) {
      return 0.0;
    }

    double mean = meanColumn(m,c);

    for (unsigned int j = 0; j < m.size2(); j++) {
      if (static_cast<int>(j) == c) {
        for (unsigned int i = 0; i < m.size1(); i++) {
          variance += ((m(i,j) - mean) * (m(i,j) - mean));
        }
      }
    }
    variance /= (nRows-1) ;
    return variance;
}

// ============================================================
// Function : VarianceRow(double)
// ------------------------------------------------------------
// Calculates the Sample Variance of a Row of doubles.
// ============================================================
double BaseStats::varianceRow(ublas::matrix<double>& m, const int& c)
{
    double variance = 0.0;

    unsigned int Vdouble_size = m.size2();
    if (Vdouble_size <= 1) {
      std::cout << "Vdoubles_size <= 1" << std::endl;
      return 0.0;
    }

    double mean = meanRow(m, c);

    for (unsigned int i = 0; i < m.size1(); i++) {
      if (i == static_cast<unsigned int>(c)) {
        for (unsigned int j = 0; j < m.size2(); j++) {
          variance+= ((m(i,j) - mean) * (m(i,j) - mean));
        }
      }
    }
    variance /= (Vdouble_size-1) ;
    return variance;
}
/*
// ============================================================
// Function : StandardDeviationColumn(double)
// ------------------------------------------------------------
// Calculates the Standard Deviation of a Column of doubles.
// ============================================================
double BaseStats::standardDeviationColumn(Matrix &Vdoubles, const int& c)
{
    double stddeviation = 0.0;

    int Vdouble_size = Vdoubles.GetDimensionN();
    if (Vdouble_size <= 1) {
      std::cout << "Vdoubles_size <= 1" << std::endl;
      return 0.0;
    }
    stddeviation = sqrt(varianceColumn(Vdoubles, c));
    return stddeviation;
}
*/
// ============================================================
// Function : StandardDeviationColumn(double)
// ------------------------------------------------------------
// Calculates the Standard Deviation of a Column of doubles.
// ============================================================
int BaseStats::getStdDevColumns(ublas::matrix<double>& mat, ublas::matrix<double>& stdDev_mat)
{
    int nRows = mat.size1();
    if (nRows <= 1) {
      return 1;
    }
    if (stdDev_mat.size1() != 1) {
      return 1;
    }

    if (mat.size2() == stdDev_mat.size2()) {
      for (unsigned int k = 0; k < mat.size2(); k++) {
        stdDev_mat(0,k) = sqrt(varianceColumn(mat, k));
      }
    }
    else {
      return 1;
    }
    return 0;
}

// ============================================================
// Function : StandardDeviationColumn(double)
// ------------------------------------------------------------
// Calculates the Standard Deviation of a Column of doubles.
// ============================================================
double BaseStats::standardDeviationColumn(ublas::matrix<double>& m, const int& c)
{
    double stddeviation = 0.0;

    int nRows = m.size1();
    if (nRows <= 1) {
      return 0.0;
    }
    stddeviation = sqrt(varianceColumn(m, c));
    return stddeviation;
}

// ============================================================
// Function : StandardDeviationRow(double)
// ------------------------------------------------------------
// Calculates the Standard Deviation of a Row of doubles.
// ============================================================
double BaseStats::standardDeviationRow(ublas::matrix<double>& m, const int& r)
{
    double stddeviation = 0.0;

    int nCols = m.size2();
    if (nCols <= 1) {
      return 0.0;
    }
    stddeviation = sqrt(varianceRow(m, r));
    return stddeviation;
}

// ============================================================
// Function : Covariance(double)
// ------------------------------------------------------------
// Calculates the Covariance between two column vectors in two
// different matrices.
// ============================================================
double BaseStats::covarianceColumn(ublas::matrix<double>& m1, const int& c,
                                   ublas::matrix<double>& m2, const int& d)
{
    double covariance = 0.0;

    if (m1.size1() <= 1 || m2.size1() <= 1) {
      std::cout << " # rows <= 1" << std::endl;
      return 0.0;
    }

    double mean1 = meanColumn(m1,c);
    double mean2 = meanColumn(m2,d);

    for (unsigned int i = 0; i < m1.size1(); i++) {
      covariance = covariance + (m1(i,c) - mean1) * (m2(i,d) - mean2);
    }

    covariance = covariance / (m1.size1()-1); // excel divides by N
    return covariance;
}

// ============================================================
// Function : CovarianceColumn(double)
// ------------------------------------------------------------
// Calculates the Covariance between two column vectors in the
// same matrix.
// ============================================================
double BaseStats::covarianceColumn(ublas::matrix<double>& matrix, const int& c, const int& d)
{
    double covariance = 0.0;

    double mean1 = meanColumn(matrix,c);
    double mean2 = meanColumn(matrix,d);

    for (unsigned int i = 0; i < matrix.size1(); i++) {
      covariance = covariance + (matrix(i,c) - mean1) * (matrix(i,d) - mean2);
    }

    covariance = covariance / (matrix.size1()-1); // excel divides by N
    return covariance;
}

// ============================================================
// Function : CovarianceMatrix(double)
// ------------------------------------------------------------
// Calculates the Covariance Matrix
//                (or variance, variance-covariance, dispersion)
// The covariance matrix is the matrix of sample variances [i][i]
// and covariances [i][j] of p variables.
// ============================================================
void BaseStats::covarianceMatrix(ublas::matrix<double>& A, ublas::matrix<double>& CovMat)
{
    unsigned int N = A.size1();
    unsigned int M = A.size2();
    double mean = 0.0;
    double mean_k = 0.0;
    double s_ii = 0.0;
    double s_ij = 0.0;

    if (M == CovMat.size1() && M == CovMat.size2()) {
      // - loop over columns - //
      for (unsigned int j = 0; j < M; j++) {
        mean = meanColumn(A,j);
        for (unsigned int i = 0; i < N; i++) {
          s_ii = s_ii + A(i,j) * A(i,j);
        }
        s_ii = (1.0/(N-1.0)) * (s_ii - (N*mean*mean));
        CovMat(j,j) = s_ii;

        for (unsigned int k = j+1; k < M; k++) {
          mean_k = meanColumn(A,k);
          for (unsigned int l = 0; l < N; l++) {
            s_ij = s_ij + A(l,k) * A(l,j);
          }
          s_ij = (1.0/(N-1.0))*(s_ij - (N*mean*mean_k));
          CovMat(j,k) = s_ij;
          s_ij = 0.0;
        }
        s_ii = 0.0;
      }

      // Fill Lower Triangle
      for (unsigned int i = 0; i < M; i++) {
        for (unsigned int j = i+1; j < M; j++) {
          CovMat(j,i) = CovMat(i,j);
        }
      }
    }
    else {
      std::cout << "BaseStats::covarianceMatrix: Not the same size"
                << std::endl;
    }
}

// ============================================================
// Function : CovarianceMatrix(double)
// ------------------------------------------------------------
// Calculates the Covariance Matrix
//                (or variance, variance-covariance, dispersion)
// The covariance matrix is the matrix of sample variances [i][i]
// and covariances [i][j] of p variables.
// ============================================================
void BaseStats::covarianceMatrix(ublas::matrix<double>& A, ublas::matrix<double, ublas::column_major>& CovMat)
{
    unsigned int N = A.size1();
    unsigned int M = A.size2();
    double mean = 0.0;
    double mean_k = 0.0;
    double s_ii = 0.0;
    double s_ij = 0.0;

    if (M == CovMat.size1() && M == CovMat.size2()) {
      // - loop over columns - //
      for (unsigned int j = 0; j < M; j++) {
        mean = meanColumn(A,j);
        for (unsigned int i = 0; i < N; i++) {
          s_ii = s_ii + A(i,j) * A(i,j);
        }
        s_ii = (1.0/(N-1.0)) * (s_ii - (N*mean*mean));
        CovMat(j,j) = s_ii;

        for (unsigned int k = j+1; k < M; k++) {
          mean_k = meanColumn(A,k);
          for (unsigned int l = 0; l < N; l++) {
            s_ij = s_ij + A(l,k) * A(l,j);
          }
          s_ij = (1.0/(N-1.0))*(s_ij - (N*mean*mean_k));
          CovMat(j,k) = s_ij;
          s_ij = 0.0;
        }
        s_ii = 0.0;
      }

      // Fill Lower Triangle
      for (unsigned int i = 0; i < M; i++) {
        for (unsigned int j = i+1; j < M; j++) {
          CovMat(j,i) = CovMat(i,j);
        }
      }
    }
    else {
      std::cout << "BaseStats::covarianceMatrix: Not the same size"
                << std::endl;
    }
}

// ============================================================
// Function : CorrelationCoefficient(double)
// ------------------------------------------------------------
// Calculates the  Correlation Coefficient between two vectors
// in different matrices.
/*
             [covariance(y,y_pred)]
r(Y,X) = ------------------------------
          [stdDev(y) * stdDev(y_pred)]
*/
// ============================================================
double BaseStats::correlationCoefficientColumn(ublas::matrix<double>& m1,
       const int& c1, ublas::matrix<double>& m2, const int& c2)
{
    double correlation = 0.0;
    double covariance  = 0.0;
    //std::cout << m1.size1() << " " << m1.size2() << std::endl;
    //std::cout << m2.size1() << " " << m2.size2() << std::endl;

    if ((m1.size1() <= 1) || (m2.size1() <= 2)) {
      std::cout << " size is 1 " << std::endl;
      return 0.0;
    }

    double mean1 = meanColumn(m1,c1);
    double mean2 = meanColumn(m2,c2);

    for (unsigned int i = 0; i < m1.size1(); i++) {
      covariance += (m1(i,c1) - mean1) * (m2(i,c2) - mean2);
    }
    double SST1 = SumSquaredDeviationsColumn(m1,c1);
    double SST2 = SumSquaredDeviationsColumn(m2,c2);
    correlation = covariance/(sqrt(SST1*SST2));

    return correlation;
}


// ============================================================
// Function : RSquared(double)
// ------------------------------------------------------------
/* Calculates the  R^2 between two vectors.

       SSR         SSE 
r^2 = ----- = 1 - -----
       TSS         TSS 

Or more simply the square of the correlation coefficient
*/
// ============================================================
double BaseStats::rSquaredColumn(ublas::matrix<double>& m1, const int& c1,
                                 ublas::matrix<double>& m2, const int& c2)
{
    double rsquared = 0.0;
    rsquared = correlationCoefficientColumn(m1, c1, m2, c2);
    rsquared = rsquared * rsquared;
    return rsquared;
}

// ============================================================
// Function : AdjustedRSquared(double)
// ------------------------------------------------------------
/*
 Calculates the Adjusted R^2 between two vectors.

             r^2 - ( 1- r^2) * p
r^2(adj) =  ---------------------
                  n - p - 1
where:
 n is the number of observations
 p is the number of independent variables
*/
// ============================================================
double BaseStats::AdjustedRSquaredColumn(ublas::matrix<double>& Ys,
       const int& c, ublas::matrix<double>& Y_Pred, const int& d, const int& p)
{
    double adjusted_rsquared = 0.0;
    double SSE = SumSquaredResidualErrorsColumn(Ys,c,Y_Pred,d);
    double SST = SumSquaredDeviationsColumn(Ys,c);
    int      n = Ys.size1();

    adjusted_rsquared = (1-((SSE/(n-p-1))/(SST/(n-1))));
    return adjusted_rsquared;
}

// ============================================================
// Function : UnsignedError(double)
// ------------------------------------------------------------
/*
  Calculates the Unsigned or Absolute Error.

         --
         \
        /  abs(y[i] - y_pred[i])
  AE =  --
       ------------------------
                N
*/
// ============================================================
double BaseStats::UnsignedError(ublas::matrix<double>& Ys,
       const int& c, ublas::matrix<double>& Y_Pred, const int& d)
{
    double AE = 0.0;
    for (unsigned int i = 0; i < Ys.size1(); i++) {
      AE += std::abs(Ys(i,c) - Y_Pred(i,d));
    }
    return (AE / double(Ys.size1()));
}

// ============================================================
// Function : SignedError(double)
// ------------------------------------------------------------
/*
  Calculates the Signed Error.

         --
         \
        /  y[i] - y_pred[i]
  AE =  --
       ------------------------
                N
*/
// ============================================================
double BaseStats::SignedError(ublas::matrix<double>& Ys,
       const int& c, ublas::matrix<double>& Y_Pred, const int& d)
{
    double AE = 0.0;
    for (unsigned int i = 0; i < Ys.size1(); i++) {
      AE += (Ys(i,c) - Y_Pred(i,d));
    }
    return (AE / double(Ys.size1()));
}

// ============================================================
// Function : RMSE(double)
// ------------------------------------------------------------
/*
  Calculates the Root Mean Squared Error.
            +                          +
            | --                       |
            | \                        |
            | /  (y[i] - y_pred[i])^2  |
RMSE = sqrt | --                       |
            | ------------------------ |
            |           N              |
            +                          +
*/
// ============================================================
double BaseStats::RMSE(ublas::matrix<double>& Ys,
       const int& c, ublas::matrix<double>& Y_Pred, const int& d)
{
    double RMSE = 0.0;
    for (unsigned int i = 0; i < Ys.size1(); i++) {
      RMSE += ((Ys(i,c) - Y_Pred(i,d)) * (Ys(i,c) - Y_Pred(i,d)));
    }
    return sqrt(RMSE / double(Ys.size1()));
}

// ============================================================
// Function : MSE(double)
// ------------------------------------------------------------
/*
  Calculates the Mean Squared Error.
       +                          +
       | --                       |
       | \                        |
       | /  (y[i] - y_pred[i])^2  |
MSE =  | --                       |
       | ------------------------ |
       |           N              |
       +                          +
*/
// ============================================================
double BaseStats::MSE(ublas::matrix<double>& Ys,
       const int& c, ublas::matrix<double>& Y_Pred, const int& d)
{
    double MSE = 0.0;
    for (unsigned int i = 0; i < Ys.size1(); i++) {
      MSE += ((Ys(i,c) - Y_Pred(i,d)) * (Ys(i,c) - Y_Pred(i,d)));
    }
    return (MSE / double(Ys.size1()));
}

// ============================================================
// Function : SumSquaredDeviationsColumn(double)
// ------------------------------------------------------------
/*
   Calculates the Sum of squared deviations
     total sum of squared deviations in Y from its mean (SST) -->
          --
          \
    SST = /  (y[i] - y_mean)^2
          --
*/
//============================================================
double BaseStats::SumSquaredDeviationsColumn(ublas::matrix<double>& Ys, const int& c)
{
    double SST = 0.0;
    double y_mean = meanColumn(Ys,c);

    for (unsigned int i = 0; i < Ys.size1(); i++) {
      SST += ((Ys(i,c) - y_mean) * (Ys(i,c) - y_mean));
    }
    return SST;
}

// ============================================================
// Function : SumSquaredRegressionColumn(double)
// ------------------------------------------------------------
/*
   Calculates the Sum of squared due to regression
     total sum of squares due to regression (SSR) -->

           _
    SSR = \  (y_pred[i] - y_obs_mean)^2
          /_
*/
//============================================================
double BaseStats::SumSquaredRegressionColumn(ublas::matrix<double>& Ys, ublas::matrix<double>& Y_Pred)
{
    double SSR = 0.0;
    double y_obs_mean = meanColumn(Ys,0);

    for (unsigned int i = 0; i < Y_Pred.size1(); i++) {
      SSR += ((Y_Pred(i,0) - y_obs_mean) * (Y_Pred(i,0) - y_obs_mean));
    }
    return SSR;
}

// ============================================================
// Function : SumSquaredResidualErrors(double)
// ------------------------------------------------------------
/*
   Calculates the Sum of squared residuals (Errors)
    residuals --> e[i] = y_obs[i] - y_pred[i]
    sum of squared residuals (SSE) -->
          --
          \
    SSE = /  e[i]^2
          --
*/
//============================================================
double BaseStats::SumSquaredResidualErrorsColumn(ublas::matrix<double>& Ys,
           ublas::matrix<double>& Y_Pred, ublas::matrix<double>& Residuals)
{
    double SSE = 0.0;
    for (unsigned int i = 0; i < Ys.size1(); i++) {
      Residuals(i,0) = Ys(i,0) - Y_Pred(i,0);
    }

    for (unsigned int i = 0; i < Ys.size1(); i++) {
      SSE += (Residuals(i,0) * Residuals(i,0));
    }
    return SSE;
}

// ============================================================
// Function : SumSquaredResidualErrors(double)
// ------------------------------------------------------------
/*
  Calculates the Sum of squared residuals (Errors)
    residuals --> e[i] = y_obs[i] - y_pred[i]
    sum of squared residuals (SSE) -->
          --
          \
    SSE = /  e[i]^2
          --
*/
//============================================================
double BaseStats::SumSquaredResidualErrorsColumn(ublas::matrix<double>& Ys, const int& c,
           ublas::matrix<double>& Y_Pred, const int& d)
{
    double SSE = 0.0;
    double residual;

    for (unsigned int i = 0; i < Ys.size1(); i++) {
      residual = Ys(i,c) - Y_Pred(i,d);
      SSE += (residual * residual);
      residual = 0.0;
    }
    return SSE;
}

} // MTKpp namespace

