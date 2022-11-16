/*!
   \file ols.h
   \brief Ordinary Least Squares
   \author Martin Peters

   $Date: 2010/03/29 20:35:21 $
   $Revision: 1.6 $

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

#ifndef OLS_H
#define OLS_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <math.h>

#include "Utils/constants.h"
#include "BaseStats.h"

/*
// - BOOST - //
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/banded.hpp> // for diagonal matrix
#include <boost/numeric/bindings/blas/blas3.hpp>
#include "boost/numeric/bindings/traits/ublas_matrix.hpp"
#include "boost/numeric/bindings/traits/ublas_vector.hpp"

namespace ublas  = boost::numeric::ublas;
namespace blas   = boost::numeric::bindings::blas;
*/

#include <Eigen/Dense>
using namespace Eigen;

namespace MTKpp
{

// ============================================================
// Class : ols()
// ------------------------------------------------------------
/*! 
   \class ols
   \brief Ordinary Least Squares
   \author Martin Peters
   \version 0.1
   \date 2006
*/
// ============================================================
class ols : public BaseStats
{
public:
    /*!
       \brief ols Constructor
    */
    ols();

    //! ols Destructor
    virtual ~ols();

    /*!
       \brief Calculate least squares regression line
       \param m1 Matrix object
       \param i1 column index
       \param m2 Matrix object
       \param i2 column index
       \return a and b from ax+b equation
    */
    //std::vector<double> LeastSquaresRegressionLine(ublas::matrix<double>& m1, const int& i1,
    //                    ublas::matrix<double>& m2, const int& i2);

    std::vector<double> LeastSquaresRegressionLine(Eigen::Matrix<double, Dynamic, Dynamic>& m1, const int& i1,
                        Eigen::Matrix<double, Dynamic, Dynamic>& m2, const int& i2);

    /*!
       \brief Ordinary Least Squares
       \param Ys Matrix object
       \param Xs Matrix object
       \param Y_Pred matrix object
    */
    //void run(ublas::matrix<double>& Ys,
    //         ublas::matrix<double>& Xs,
    //         ublas::matrix<double>& Y_Pred);

    void run(Eigen::Matrix<double, Dynamic, Dynamic>& Ys,
             Eigen::Matrix<double, Dynamic, Dynamic>& Xs,
             Eigen::Matrix<double, Dynamic, Dynamic>& Y_Pred);
};

} // MTKpp namespace

#endif // OLS_H

