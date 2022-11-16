/*!
   \file pca.h
   \brief Principal Component Analysis
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

#ifndef PCA_H
#define PCA_H

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
#include "table.h"

namespace MTKpp
{

class sheet;
// ============================================================
// Class : pca()
// ------------------------------------------------------------
/*! 
   \class pca
   \brief Principal Component Analysis
   \author Martin Peters
   \date 2006
*/
// ============================================================
class pca : public BaseStats
{
public:
    /*!
       \brief pca Constructor
    */
    pca();

    /*!
       \brief pca Constructor
       \param X X matrix
       \param output sheet pointer
    */
    pca(table<double>* X, sheet* output);

    //! pca Destructor
    //virtual ~pca();

    /*!
       \brief Run PCA 
       \param nKeep Number of components to keep
    */
    int run(int nKeep);

protected:
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
       \brief Number of rows in X
    */
    unsigned int   nRows;

    /*!
       \brief Number of Columns in X
    */
    unsigned int   nColumns;

    /*!
       \brief The sheet where the model is stored
    */
    sheet*         outModel;
};

} // MTKpp namespace

#endif // PCA_H
