/*!
   \file ols.cpp
   \brief Ordinary Least Squares
   \author Martin Peters

   $Date: 2007/11/28 09:22:17 $
   $Revision: 1.4 $

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

#include "ols.h"

namespace MTKpp
{

// ============================================================
// Class : ols()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
ols::ols():BaseStats() {}

// ============================================================
// Function : ~ols()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
ols::~ols() {}

// ============================================================
// Function : LeastSquaresRegressionLine(double)
// ------------------------------------------------------------
/*
Calculates the  Least Squares Regression Line
 y = mx + c
   y is the dependent variable
   x in the independent variable

     sum[(y_pred[i] - <y_pred>)(y[i]-<y>)]
m = ---------------------------------------
       sum[(y_pred[i] - <y_pred>)^2

c = <y> - m<x>

Returns the list [m,c]
*/
// ============================================================
//std::vector<double> ols::LeastSquaresRegressionLine(ublas::matrix<double>& m1,
//                   const int& i1, ublas::matrix<double>& m2, const int& i2)
std::vector<double> ols::LeastSquaresRegressionLine(Eigen::Matrix<double, Dynamic, Dynamic>& m1,
                   const int& i1, Eigen::Matrix<double, Dynamic, Dynamic>& m2, const int& i2)
{
    double m = 0.0;
    double c = 0.0;
    std::vector<double> line;

    double covariance = 0.0;
    double mean1 = meanColumn(m1,i1);
    double mean2 = meanColumn(m2,i2);

    for (unsigned int i = 0; i < m1.rows(); i++) {
      covariance = covariance + (m1(i,i1) - mean1) * (m2(i,i2) - mean2);
    }
    double SST = SumSquaredDeviationsColumn(m2,i2);

    m = covariance/SST;
    line.push_back(m);

    c = mean1 - m*mean2;
    line.push_back(c);
    return line;
}

// ============================================================
// Function : OLS(double)
// ------------------------------------------------------------
// Calculates the  Least Squares Regression Line
// y = mx + c
// ============================================================
void ols::run(Eigen::Matrix<double, Dynamic, Dynamic>& Ys,
              Eigen::Matrix<double, Dynamic, Dynamic>& Xs,
              Eigen::Matrix<double, Dynamic, Dynamic>& Y_Pred)
{
/*    Y    =              X               B    +   E
    +-  -+   +-                      -+ +-  -+   +-  -+
    | Y1 |   | 1 X11 X12  .   .   X1q | | B1 |   | E1 |
    | Y2 |   | 1 X21 X22  .   .   X2q | | B2 |   | E2 |
    | .  | = | 1  .   .   .   .    .  | | .  | + | .  |
    | .  |   | 1  .   .   .   .    .  | | .  |   | .  |
    | Yn |   | 1 Xn1  .   .   .   Xnq | | Bq |   | En |
    +-  -+   +-                      -+ +-  -+   +-  -+

     B = (X'X)^-1X'Y
*/
/*
    ublas::matrix<double> B(Xs.size1(),Xs.size2()+1);

    Matrix X(Xs.GetDimensionN(),Xs.GetDimensionM()+1);
    for (int i = 0; i < X.GetDimensionN(); i++) {
      for (int j = 0; j < X.GetDimensionM(); j++) {
        if ( j == 0.0) {
          X[i][j] = 1.0;
        }
        else {
          X[i][j] = Xs[i][j-1];
        }
      }
    }

    Matrix X_transpose(Xs.GetDimensionM()+1, Xs.GetDimensionN());

    X.transpose(X_transpose);

    SquareMatrix X_tr_x_X(X_transpose.GetDimensionN()); 

    Matrix temp = X_transpose*X;
    X_tr_x_X = temp;

    std::cout << "X_TR_x_X = \n" << X_tr_x_X << std::endl;

    Matrix X_tr_x_Y(X_transpose.GetDimensionN(), Ys.GetDimensionM());

    X_tr_x_Y = X_transpose*Ys;

    std::cout << "X_TR_x_Y = \n" << X_tr_x_Y << std::endl;

    Matrix inv_X_tr_x_X(X_transpose.GetDimensionN(), Xs.GetDimensionM());

    double det;
    inv_X_tr_x_X = X_tr_x_X.inverse(det);
    std::cout << "DET = " << det << std::endl;

    std::cout << "inv_X_TR_x_X = \n" << inv_X_tr_x_X << std::endl;

    Matrix Bs = inv_X_tr_x_X * X_tr_x_Y;

    std::cout << "Bs = \n" << Bs << std::endl;

    Y_Pred = X*Bs;

    X_transpose.DeleteMatrix();
    X_tr_x_X.DeleteMatrix();
    X_tr_x_Y.DeleteMatrix();
    inv_X_tr_x_X.DeleteMatrix();
    Bs.DeleteMatrix();
*/
}

/*

intercept of 0:
R^2 = sum [y - (a+bx)]^2

dR^2        dR^2
---- =  0,  ---- = 0
 da          db

dR^2
---- = 0 = -2sum[y-bx]x
 db

sum(xy) = bsum(x^2)

    sum(xy)
b = -------
     sum(x^2)
*/

} // MTKpp namespace

