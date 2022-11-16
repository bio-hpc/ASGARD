/*!
   \file pls.h
   \brief Partial Least Squares
   \author Martin Peters

   $Date: 2010/03/29 20:35:21 $
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

#ifndef PLS_H
#define PLS_H

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
#include "BaseStats.h"

// - BOOST - //
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/banded.hpp> // for diagonal matrix
#include <boost/numeric/bindings/blas/blas3.hpp>
#include "boost/numeric/bindings/traits/ublas_matrix.hpp"
#include "boost/numeric/bindings/traits/ublas_vector.hpp"

#include "table.h"

namespace ublas  = boost::numeric::ublas;
namespace blas   = boost::numeric::bindings::blas;

namespace MTKpp
{

class sheet;
//class table;

// ============================================================
// Class : pls()
// ------------------------------------------------------------
/*!
   \class pls
   \brief Partial Least Squares
   \author Martin Peters
*/
// ============================================================
class pls : public BaseStats
{
public:

    /*!
       \brief pls Constructor
    */
    pls();

    /*!
       \brief pls Constructor
       \param Y Y matrix
       \param X X matrix
       \param method method
       \param nlvs number of latent variables
       \param cv cross validation method
       \param bError error boolean
    */
    pls(table<double>* Y, table<double>* X, std::string method, int nlvs, std::string cv,
        bool& bError);

    /*!
       \brief pls Constructor
       \param Y Y matrix
       \param X X matrix
       \param method method
       \param nlvs number of latent variables
       \param cv cross validation method
       \param output sheet pointer
       \param bError error boolean
    */
    pls(table<double>* Y, table<double>* X, std::string method, int nlvs, std::string cv,
        sheet* output, bool& bError);

    /*!
       \brief pls Constructor
       \param Y Y name
       \param X X name
       \param S sheet pointer
       \param method method
       \param nlvs number of latent variables
       \param cv cross validation method
       \param output sheet pointer
       \param bError error boolean
    */
    pls(std::string Y, std::string X, sheet* S, std::string method, int nlvs,
        std::string cv, sheet* output, bool& bError);

    //! pls Destructor
    //virtual ~pls();

    /*!
       \brief Run PLS
       \param bError error boolean
    */
    void run(bool& bError);

    /*!
       \brief Run CV PLS
       \param bError error boolean
    */
    void runCV(bool& bError);

    /*!
       \brief Set X matrix
       \param x table pointer
    */
    void setX(table<double>* x);

    /*!
       \brief Set Y matrix
       \param y table pointer
    */
    void setY(table<double>* y);

    /*!
       \brief Set PLS Algorithm
       \param g PLS Algorithm
    */
    void setMethod(std::string g);

    /*!
       \brief Set Maximum number of iterations for iterative methods
       \param i max number of iterations
    */
    void setMaxIter(int i);

    /*!
       \brief Set Convergence criteria for iterative methods
       \param e Convergence criteria
    */
    void setEpsilon(double e);

    /*!
       \brief Set Cross validation method
       \param c CV method
    */
    void setCV(std::string c);

    /*!
       \brief Set number of samples to consider in RANDOM CV
       \param i number of samples
    */
    void setNITER(int i);

    /*!
       \brief Set size of test set in RANDOM and LNO CV
       \param i size of test set
    */
    void setNTEST(int i);

    /*!
       \brief Set random number generator seed in RANDOM CV
       \param s random number generator seed
    */
    void setSEED(int s);

    /*!
       \brief Set number of Latent Variables to be considered
       \param l number of LVs
    */
    void setNLV(int l);

    /*!
       \brief The sheet where the model is stored
       \param s sheet pointer
    */
    void setOutModel(sheet* s);

protected: // Functions

    /*!
       \brief Kernel Partial Least Squares

                Y[N][M]                   X[N][R]
       +-                  -+   +-                    -+
       | Y11 Y21  .  .  Y1M |   | X11 X12  .   .   X1R |
       | Y21 Y22  .  .  Y2M |   | X21 X22  .   .   X2R |
       |  .    .  .  .   .  |   |  .   .   .   .    .  |
       |  .    .  .  .   .  |   |  .   .   .   .    .  |
       | YN1 Y2N  .  .  YNM |   | XN1  .   .   .   XNR |
       +-                  -+   +-                    -+

       PLS regression searches for a set of components (latent vectors)
       that performs a simultaneous decomposition of X and Y with the
       constraint that these components explain as much as possible of
       the covariance between them.  Then a regression step, where the
       decomposition of X is used to predict Y is performed.

       T -- Score Matrix for X      T[N][MaxComponents]
       U -- Score Matrix for Y      U[N][MaxComponents]
       P -- Loading Matrix for X    P[R][MaxComponents]
       C -- Weight Matrix for Y     C[M][MaxComponents]
       W -- Weighting Matrix        W[R][MaxComponents]
       E -- Residual Matrix for X   E[N][R]
       F -- Residual Matrix for Y   F[N][M]
       B -- Regression Coefficients B[N][1]

       X = TP' + E
       Y = TC' + E
       
       Y_pred = X*B + E
       where:
        B = W * inv(P' * W) * C' 

       -- Reference: Herve Abdi, Partial Least Squares (PLS) Regression,
                     University of Texas at Dallas
    */
    int kernelPLS();

protected: // Data
    /*!
       \brief X matrix
       \code
              X[N][R]
         +-                    -+
         | X11 X12  .   .   X1R |
         | X21 X22  .   .   X2R |
         |  .   .   .   .    .  |
         |  .   .   .   .    .  |
         | XN1  .   .   .   XNR |
         +-                    -+
       \endcode
    */
    table<double>* itsX;

    /*!
       \brief Y matrix
       \code
              Y[N][M]
        +-                  -+
        | Y11 Y21  .  .  Y1M |
        | Y21 Y22  .  .  Y2M |
        |  .    .  .  .   .  |
        |  .    .  .  .   .  |
        | YN1 Y2N  .  .  YNM |
        +-                  -+
       \endcode
    */
    table<double>* itsY;

    /*!
       \brief Number of rows in Y and X
    */
    unsigned int YRows;

    /*!
       \brief Number of Columns in X
    */
    unsigned int XColumns;

    /*!
       \brief PLS Algorithm
        - KERNELPLS
        - SIMPLS
        - NIPLS
    */
    std::string itsMethod;

    /*!
       \brief Maximum number of iterations for iterative methods
    */
    unsigned int maxIter;

    /*!
       \brief Convergence criteria for iterative methods
    */
    double epsilon;

    /*!
       \brief Cross validation parameters
        - NONE : No CV
        - LOO : Leave One Out cross validation
        - LNO : Leave N Out cross validation, set nTEST parameter
        - RANDOM : Leave nTest out nITER times using SEED
    */
    std::string CV;

    /*!
       \brief Number of prediction set in CV
        nEXT = N/nTest
    */
    unsigned int nEXT;

    /*!
       \brief Number of samples to consider in RANDOM CV
    */
    unsigned int nITER;

    /*!
       \brief Size of test set in RANDOM and LNO CV

       if N = 100
       - LOO: nTest = 1
       - LNO: nTest = 10 ==> nExt = 10
       - RANDOM: nTest = 10 ==> nExt = 10
    */
    unsigned int nTEST;

    /*!
       \brief Random number generator seed in RANDOM CV
    */
    int SEED;

    /*!
       \brief Number of Latent Variables to be considered
    */
    unsigned int nLV;

    /*!
       \brief The sheet where the model is stored
    */
    sheet* outModel;
};

} // MTKpp namespace

#endif // PLS_H
