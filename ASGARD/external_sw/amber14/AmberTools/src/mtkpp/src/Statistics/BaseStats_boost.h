/*!
   \file BaseStats.h
   \brief Base class for statistical routines
   \author Martin Peters

   $Date: 2010/03/29 20:35:21 $
   $Revision: 1.10 $

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

#ifndef BASESTATS_h
#define BASESTATS_h

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <math.h>

#include "Utils/constants.h"
#include "Utils/object.h"

// - BOOST - //
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/banded.hpp> // for diagonal matrix
#include <boost/numeric/bindings/blas/blas3.hpp>
#include "boost/numeric/bindings/traits/ublas_matrix.hpp"
#include "boost/numeric/bindings/traits/ublas_vector.hpp"

namespace ublas  = boost::numeric::ublas;
namespace blas   = boost::numeric::bindings::blas;

class object;

namespace MTKpp
{
// ============================================================
// Class : BaseStats()
// ------------------------------------------------------------
/*!
   \class BaseStats
   \brief Base class for statistical routines
   \author Martin Peters
   \date 2005
*/
// ============================================================
class BaseStats : public object
{
public:
    /*!
       \brief BaseStats Constructor
    */
    BaseStats();

    //! BaseStats Destructor
    //virtual ~BaseStats();

    /*!
       \brief Get sample mean of column
       \param m Matrix pointer
       \param i column index
       \return sample average
    */
    double         meanColumn(ublas::matrix<double>& m, const int& i);

    /*!
       \brief Get sample mean of row
       \param m Matrix pointer
       \param i row index
       \return sample average
    */
    double         meanRow(ublas::matrix<double>& m, const int& i);

    /*!
       \brief Get sum of column
       \param m Matrix pointer
       \param i column index
       \return sum of values
    */
    double         sumColumn(ublas::matrix<double>& m, const int& i);

    /*!
       \brief Get sum of row
       \param m Matrix pointer
       \param i row index
       \return sum of values
    */
    double         sumRow(ublas::matrix<double>& m, const int& i);

    /*!
       \brief Get max value of column
       \param m Matrix pointer
       \param i column index
       \return max value
    */
    double         maxColumn(ublas::matrix<double>& m, const int& i);

    /*!
       \brief Get max value of column
       \param m Matrix pointer
       \param i column index
       \param r row index [index of max value]
       \return max value
    */
    double         maxColumn(ublas::matrix<double>& m, const int& i, int& r);

    /*!
       \brief Get max value of row
       \param m Matrix pointer
       \param i row index
       \return max value
    */
    double         maxRow(ublas::matrix<double>& m, const int& i);

    /*!
       \brief Get min value of column
       \param m Matrix pointer
       \param i column index
       \return min value
    */
    double         minColumn(ublas::matrix<double>& m, const int& i);

    /*!
       \brief Get min value of column
       \param m Matrix pointer
       \param i column index
       \param r row index [index of min value]
       \return min value
    */
    double         minColumn(ublas::matrix<double>& m, const int& i, int& r);

    /*!
       \brief Get min value of row
       \param m Matrix pointer
       \param i row index
       \return min value
    */
    double         minRow(ublas::matrix<double>& m, const int& i);

    /*!
       \brief Get Mean values for each column in the matrix
       \param mat Matrix pointer
       \param mat_centers Matrix pointer
       \return sucess
    */
    int getColumnCenters(ublas::matrix<double>& mat, ublas::matrix<double>& mat_centers);

    /*!
       \brief Center matrix by column
       \param m1 Matrix pointer
       \param m2 Matrix pointer
    */
    void           centerColumns(ublas::matrix<double>& m1, ublas::matrix<double>& m2);

    /*!
       \brief Center matrix by row
       \param m1 Matrix pointer
       \param m2 Matrix pointer
    */
    void           centerRows(ublas::matrix<double>& m1, ublas::matrix<double>& m2);

    /*!
       \brief Calculate Z-Scores of matrix by column
       \param m1 Matrix pointer
       \param m2 Matrix pointer
    */
    void           zScoreColumns(ublas::matrix<double>& m1, ublas::matrix<double>& m2);

    /*!
       \brief Calculate Z-Scores of matrix by row
       \param m1 Matrix pointer
       \param m2 Matrix pointer
    */
    void           zScoreRows(ublas::matrix<double>& m1, ublas::matrix<double>& m2);

    /*!
       \brief Autoscale using previously computed means and standard deviations
       \param old_mat matrix pointer
       \param new_mat matrix pointer
       \param centers matrix of means
       \param stdDevs matrix of standard deviations
       \return success
    */
    int            autoScale(ublas::matrix<double>& old_mat,
                             ublas::matrix<double>& new_mat,
                             ublas::matrix<double>& centers,
                             ublas::matrix<double>& stdDevs);

    /*!
       \brief Get variance of column
       \param m Matrix pointer
       \param i column index
       \return variance
    */
    double         varianceColumn(ublas::matrix<double>& m, const int& i);

    /*!
       \brief Get variance of row
       \param m Matrix pointer
       \param i row index
       \return variance
    */
    double         varianceRow(ublas::matrix<double>& m, const int& i);

    /*!
       \brief Get standard deviation for each column
       \param mat matrix pointer
       \param stdDev_mat matrix of standard deviations
       \return success
    */
    int getStdDevColumns(ublas::matrix<double>& mat, ublas::matrix<double>& stdDev_mat);

    /*!
       \brief Get standard deviation of column
       \param m Matrix pointer
       \param i column index
       \return standard deviation
    */
    double         standardDeviationColumn(ublas::matrix<double>& m, const int& i);

    /*!
       \brief Get standard deviation of row
       \param m Matrix pointer
       \param i row index
       \return standard deviation
    */
    double         standardDeviationRow(ublas::matrix<double>& m, const int& i);

    // - Covariances - //
    /*!
       \brief Calculate covariance between two column in different matrices
       \param m1 Matrix pointer
       \param i1 column index
       \param m2 Matrix pointer
       \param i2 column index
       \return covariance
    */
    double         covarianceColumn(ublas::matrix<double>& m1, const int& i1,
                                      ublas::matrix<double>& m2, const int& i2);

    /*!
       \brief Calculate covariance between two columns in the same matrix
       \param m Matrix pointer
       \param i column index
       \param j column index
       \return covariance
    */
    double         covarianceColumn(ublas::matrix<double>& m, const int& i, const int& j);

    /*!
       \brief Calculates the Covariance Matrix (or variance, variance-covariance, dispersion)

       The covariance matrix is the matrix of sample variances [i][i]
       and covariances [i][j] of p variables.

       \param A Matrix object
       \param CovMat covariance matrix
    */
    void           covarianceMatrix(ublas::matrix<double>& A, ublas::matrix<double>& CovMat);
    void           covarianceMatrix(ublas::matrix<double>& A, ublas::matrix<double, ublas::column_major>& CovMat);

    /*!
       \brief Calculate correlation coefficient by column
       \param m1 Matrix object
       \param i1 column index
       \param m2 Matrix object
       \param i2 column index
       \return correlation coefficient
    */
    double         correlationCoefficientColumn(ublas::matrix<double>& m1, const int& i1,
                                                ublas::matrix<double>& m2, const int& i2);

    /*!
       \brief Calculate R^2 by column
       \param m1 Matrix object
       \param i1 column index
       \param m2 Matrix object
       \param i2 column index
       \return R-Squared
    */
    double         rSquaredColumn(ublas::matrix<double>& m1, const int& i1,
                                  ublas::matrix<double>& m2, const int& i2);

    /*!
       \brief Calculate Adjusted R^2 by column
       \param Ys Matrix object
       \param i1 column index
       \param Y_Pred Matrix object
       \param i2 column index
       \param i3 column index
       \return Adjusted R-Squared
    */
    double         AdjustedRSquaredColumn(ublas::matrix<double>& Ys, const int& i1,
                       ublas::matrix<double>& Y_Pred, const int& i2, const int& i3);

    /*!
       \brief Calculate root mean squared error
       \param Ys Matrix object
       \param i column index
       \param Y_Pred Matrix object
       \param j column index
       \return rmse
    */
    double         RMSE(ublas::matrix<double>& Ys, const int& i,
                        ublas::matrix<double>& Y_Pred, const int& j);

    /*!
       \brief Calculate mean squared error
       \param Ys Matrix object
       \param i column index
       \param Y_Pred Matrix object
       \param j column index
       \return mse
    */
    double         MSE(ublas::matrix<double>& Ys, const int& i,
                       ublas::matrix<double>& Y_Pred, const int& j);

    /*!
       \brief Calculate unsigned/absolute error
       \param Ys Matrix object
       \param i column index
       \param Y_Pred Matrix object
       \param j column index
       \return unsigned error
    */
    double         UnsignedError(ublas::matrix<double>& Ys, const int& i,
                                 ublas::matrix<double>& Y_Pred, const int& j);

    /*!
       \brief Calculate signed error
       \param Ys Matrix object
       \param i column index
       \param Y_Pred Matrix object
       \param j column index
       \return signed error
    */
    double         SignedError(ublas::matrix<double>& Ys, const int& i,
                               ublas::matrix<double>& Y_Pred, const int& j);

    /*!
       \brief Calculates the Sum of squared deviations

        total sum of squared deviations in Y from its mean (SST)
              --
              \
        SST = /  (y[i] - y_mean)^2
              --

       \param Ys Matrix object
       \param i column index
       \return Sum Squared Deviation
    */
    double         SumSquaredDeviationsColumn(ublas::matrix<double>& Ys, const int& i);

    /*!
       \brief Calculates the Sum of squared due to regression

       total sum of squares due to regression (SSR)

              __
       SSR = \  (y_pred[i] - y_obs_mean)^2
             /__

       \param Ys Matrix object
       \param Y_Pred Matrix object
       \return Sum Squared Deviation
    */
    double         SumSquaredRegressionColumn(ublas::matrix<double>& Ys,
                                              ublas::matrix<double>& Y_Pred);

    /*!
       \brief Calculates the Sum of squared residuals (Errors)

       residuals --> e[i] = y_obs[i] - y_pred[i]
       sum of squared residuals (SSE)
             --
             \
       SSE = /  e[i]^2
             --
       \param Ys
       \param Y_Pred
       \param Residuals
    */
    double         SumSquaredResidualErrorsColumn(ublas::matrix<double> &Ys,
                     ublas::matrix<double> &Y_Pred,
                     ublas::matrix<double>& Residuals);

    /*!
       \brief Calculates the Sum of squared residuals (Errors)

       residuals --> e[i] = y_obs[i] - y_pred[i]
       sum of squared residuals (SSE)
             --
             \
       SSE = /  e[i]^2
             --
       \param Ys Orginial Y matrix
       \param c column in Ys
       \param Y_Pred Predicted Y matrix 
       \param d column in Y_Pred
    */
    double         SumSquaredResidualErrorsColumn(ublas::matrix<double> &Ys,
                     const int& c, ublas::matrix<double> &Y_Pred, const int& d);
};

} // MTKpp namespace

#endif // BASESTATS_H

