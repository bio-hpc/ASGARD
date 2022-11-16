/*!
   \file linearAlgebraTest.cpp
   \brief Tests the Boost functionality in MTK++
   \author Martin Peters

   $Date: 2010/08/11 21:20:18 $
   $Revision: 1.3 $

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
#include "Utils/printHeader.h"
#include "Log/errorHandler.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include "boost/numeric/bindings/traits/ublas_vector.hpp"

#include "boost/numeric/bindings/lapack/syev.hpp"
//#include <boost/numeric/bindings/lapack/gesvd.hpp>

#include <boost/numeric/bindings/lapack/gesv.hpp>

#include <boost/numeric/bindings/traits/ublas_vector2.hpp>

namespace ublas = boost::numeric::ublas;
namespace lapack= boost::numeric::bindings::lapack;

#include "Utils/diagonalize.h"
#include "Parsers/StringManip.h"

#include "Parsers/commLineOptions.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;
using namespace MTKpp;

/*!
   \brief Tests MTK++'s/Boost functionality
   \param argc
   \param argv
   \return success
*/
int main (int argc, char **argv)
{
    std::string prog_name = "linearAlgebraTest";
    std::vector<std::string> authors;
    std::string author = "Martin B. Peters";
    authors.push_back(author);

    // 1. CREATE AN OBJECT
    commLineOptions *clo = new commLineOptions();

    // 2. SET PREFERENCES
    clo->noUsage();

    // 3. SET THE USAGE/HELP
    clo->addUsage( "  linearAlgebraTest                                  \n" );
    clo->addUsage( "    usage:  linearAlgebraTest [flags] [options]      \n" );
    clo->addUsage( "  options:                                             " );
    clo->addUsage( "          -l log file                                \n" );
    clo->addUsage( "          -o out file                                \n" );
    clo->addUsage( "          -b build directory                         \n" );
    clo->addUsage( "    flags:                                             " );
    clo->addUsage( "          -h help                                    \n" );
    clo->addUsage( "   author:                                             " );
    clo->addUsage( "          Martin B. Peters (c) 2010                    " );

    // 4. SET THE OPTION STRINGS/CHARACTERS
    clo->setOption( "logFile", 'l' );
    clo->setOption( "outFile", 'o' );
    clo->setOption( "buildDir", 'b' );

    // 5. PROVIDE THE COMMANDLINE
    clo->processCommandArgs( argc, argv );

    clo->usageOn();

    // 6. GET THE VALUES
    if ( clo->getFlag( "help" ) || clo->getFlag( 'h' ) ) {
      clo->printUsage();
      return 0;
    }

    std::string logFile = "linearAlgebraTest.log";
    std::string outFile = "linearAlgebraTest.out";
    std::string buildDir = "";

    if ( clo->getValue( "l" ) != 0 ) {
      logFile = std::string(clo->getValue( "l" ));
    }
    else if ( clo->getValue( "logFile" ) != 0 ) {
      logFile =  std::string(clo->getValue( "logFile" ));
    }

    if ( clo->getValue( "o" ) != 0 ) {
      outFile = std::string(clo->getValue( "o" ));
    }
    else if ( clo->getValue( "outFile" ) != 0 ) {
      outFile =  std::string(clo->getValue( "outFile" ));
    }

    if ( clo->getValue( "b" ) != 0 ) {
      buildDir = std::string(clo->getValue( "b" ));
    }
    else if ( clo->getValue( "buildDir" ) != 0 ) {
      buildDir =  std::string(clo->getValue( "buildDir" ));
    }
    else {
      std::cout << "\nNo build directory provided ... exiting " << std::endl;
      exit(1);
    }

    // 6. OPEN LOG AND OUT FILES
    std::ofstream oLog;
    oLog.open(logFile.c_str());

    if (!oLog) {
      std::cout << "\nUNABLE TO OPEN LOG FILE"
                << "\nFILENAME = " << logFile << std::endl;
      exit(1);
    }

    std::ofstream oOut;
    oOut.open(outFile.c_str());

    if (!oOut) {
      oLog.close();
      std::cout << "\nUNABLE TO OPEN OUT FILE"
                << "\nFILENAME = " << outFile << std::endl;
      exit(1);
    }

    // Set errorLog stream to the log file
    MTKpp::errorLogger.setStream(&oLog);

    // Print MTK++ copyright message
    printHeader(oLog, prog_name, authors);

    // 7. DONE
    delete clo;

    // Start
    std::string errorMessage = "";

    // Start & end time
    time_t startTime;
    time_t endTime;

    // Get start time
    time (&startTime);

    std::string errMessage = " Real Symmetric Matrix Diagonalization test ";

/*
    ublas::matrix<double,ublas::column_major> A(3,3);
    ublas::vector<double> b(3);

    A(0, 0) = 0.0; A(0, 1) =  1.0; A(0, 2) = 1.0;
    A(1, 0) = 1.0; A(1, 1) =  0.0; A(1, 2) = 1.0;
    A(2, 0) = 1.0; A(2, 1) =  1.0; A(2, 2) = 0.0;

    b(0) = 0.0; b(1) = 0.0; b(2) = 0.0;
*/
    ublas::matrix<double,ublas::column_major> A(5,5);
    ublas::vector<double> b(5);

/*
   http://software.intel.com/sites/products/documentation/hpc/mkl/lapack/mkl_lapack_examples/dsyev_ex.c.htm

   DSYEV Example.
   ==============

   Program computes all eigenvalues and eigenvectors of a real symmetric
   matrix A:

     1.96  -6.49  -0.47  -7.20  -0.65
    -6.49   3.80  -6.39   1.50  -6.34
    -0.47  -6.39   4.17  -1.51   2.67
    -7.20   1.50  -1.51   5.70   1.80
    -0.65  -6.34   2.67   1.80  -7.10

   Eigenvalues
    -11.07  -6.23   0.86   8.87  16.09

   Eigenvectors (stored columnwise)
    -0.30  -0.61   0.40  -0.37   0.49
    -0.51  -0.29  -0.41  -0.36  -0.61
    -0.08  -0.38  -0.66   0.50   0.40
     0.00  -0.45   0.46   0.62  -0.46
    -0.80   0.45   0.17   0.31   0.16
*/

    A(0, 0) =  1.96; A(0, 1) =  -6.49; A(0, 2) = -0.47; A(0, 3) = -7.20; A(0, 4) = -0.65;
    A(1, 0) = -6.49; A(1, 1) =   3.80; A(1, 2) = -6.39; A(1, 3) =  1.50; A(1, 4) = -6.34;
    A(2, 0) = -0.47; A(2, 1) =  -6.39; A(2, 2) =  4.17; A(2, 3) = -1.51; A(2, 4) =  2.67;
    A(3, 0) = -7.20; A(3, 1) =   1.50; A(3, 2) = -1.51; A(3, 3) =  5.70; A(3, 4) =  1.80;
    A(4, 0) = -0.65; A(4, 1) =  -6.34; A(4, 2) =  2.67; A(4, 3) =  1.80; A(4, 4) = -7.10;

    b(0) = 0.0; b(1) = 0.0; b(2) = 0.0; b(3) = 0.0; b(4) = 0.0;

    std::stringstream A_ss;
    A_ss << A;
    errMessage += "\n ### Start A = " + A_ss.str();

    std::stringstream b_ss;
    b_ss << b;
    errMessage += "\n ### Start b = " + b_ss.str() +"\n";

    int r = diagonalize(A, b);
    if (r == 0) {
      //std::cout << " A \n" << A << std::endl;
      //std::cout << " b \n" << b << std::endl;

      std::stringstream A2_ss;
      A2_ss << A;
      errMessage += " ### A = " + A2_ss.str() +"\n";

      std::stringstream b2_ss;
      b2_ss << b;
      errMessage += " ### b = " + b2_ss.str();
    }
    else {
      std::cout << " Error occured during diagonalization" << std::endl;
    }

    oOut << errMessage << std::endl;
    errorLogger.throwError("linearAlgebraTest", errMessage, INFO);
/*
    // ------- SVD -------
    // general tutorial http://www.puffinwarellc.com/index.php/news-and-articles/articles/30-singular-value-decomposition-tutorial.html
    //A_nxp = U_nxn S_nxp V^T_pxp

    errMessage = " Singular Value Decomposition test 1 ";

    // http://web.mit.edu/be.400/www/SVD/Singular_Value_Decomposition.htm
    ublas::matrix<double, ublas::column_major> A1(4, 2);
    ublas::vector<double> s1(8);
    ublas::matrix<double, ublas::column_major> U1 (4, 4);
    ublas::matrix<double, ublas::column_major> VT1 (2, 2);

    A1(0, 0) = 2.0; A1(0, 1) = 4.0;
    A1(1, 0) = 1.0; A1(1, 1) = 3.0;
    A1(2, 0) = 0.0; A1(2, 1) = 0.0;
    A1(3, 0) = 0.0; A1(3, 1) = 0.0;

    //std::cout << " A \n" << A1 << std::endl;

    std::stringstream A1_ss;
    A1_ss << A1;
    errMessage += "\n ### Start A = " + A1_ss.str();

    r = svd(A1, s1, U1, VT1);
    if (r == 0) {
      //std::cout << "\n A \n" << A1 << std::endl;
      //std::cout << "\n U \n" << U1 << std::endl;
      //std::cout << "\n s \n" << s1 << std::endl;
      //std::cout << "\n VT \n" << VT1 << std::endl;

      std::stringstream A12_ss;
      A12_ss << A1;
      errMessage += "\n ### A = " + A12_ss.str();

      std::stringstream U1_ss;
      U1_ss << U1;
      errMessage += "\n ### U = " + U1_ss.str();

      std::stringstream s1_ss;
      s1_ss << s1;
      errMessage += "\n ### s = " + s1_ss.str();

      std::stringstream VT1_ss;
      VT1_ss << VT1;
      errMessage += "\n ### VT = " + VT1_ss.str();

      // Expected results
      // U = ((0.82, -0.58, 0,0), (0.58, 0.82, 0, 0), (0, 0, 1, 0), (0, 0 ,0 1))
      // S = (5.47, 0.37)
      // V = ((0.4, -0.91), (0.91, 0.40))
    }
    else {
      std::cout << " Error occured during SVD" << std::endl;
    }

    errorLogger.throwError("linearAlgebraTest", errMessage, INFO);
    oOut << errMessage << std::endl;

    errMessage = " Singular Value Decomposition test 2 ";

    //std::cout << "  \n\n SVD2 " << std::endl;
    // http://en.wikipedia.org/wiki/Singular_value_decomposition
    ublas::matrix<double, ublas::column_major> A2(4,5);
    ublas::vector<double> S2(4);
    ublas::matrix<double, ublas::column_major> U2 (4, 4);
    ublas::matrix<double, ublas::column_major> VT2 (5, 5);

    A2(0, 0) = 1.0; A2(0, 1) = 0.0; A2(0, 2) = 0.0; A2(0, 3) = 0.0; A2(0, 4) = 2.0;
    A2(1, 0) = 0.0; A2(1, 1) = 0.0; A2(1, 2) = 3.0; A2(1, 3) = 0.0; A2(1, 4) = 0.0;
    A2(2, 0) = 0.0; A2(2, 1) = 0.0; A2(2, 2) = 0.0; A2(2, 3) = 0.0; A2(2, 4) = 0.0;
    A2(3, 0) = 0.0; A2(3, 1) = 4.0; A2(3, 2) = 0.0; A2(3, 3) = 0.0; A2(3, 4) = 0.0;

    //std::cout << " A \n" << A2 << std::endl;
    std::stringstream A2_ss;
    A2_ss << A2;
    errMessage += "\n ### Start A = " + A2_ss.str();

    r = svd(A2, S2, U2, VT2);
    if (r == 0) {
      //std::cout << "\n A \n" << A2 << std::endl;
      //std::cout << "\n U \n" << U2 << std::endl;
      //std::cout << "\n S \n" << S2 << std::endl;
      //std::cout << "\n VT \n" << VT2 << std::endl;

      std::stringstream A22_ss;
      A22_ss << A2;
      errMessage += "\n ### A = " + A22_ss.str();

      std::stringstream U2_ss;
      U2_ss << U2;
      errMessage += "\n ### U = " + U2_ss.str();

      std::stringstream S2_ss;
      S2_ss << S2;
      errMessage += "\n ### s = " + S2_ss.str();

      std::stringstream VT2_ss;
      VT2_ss << VT2;
      errMessage += "\n ### VT = " + VT2_ss.str();

      // Expected results
      // U = ((0, 0, 1, 0), (0, 1, 0, 0), (0, 0, 0, -1), (1, 0 ,0, 0))
      // S = (4, 3, 2.24, 0)
      // VT = ((0, 1, 0, 0, 0), (0, 0, 1, 0, 0), (0.45, 0,0,0, 0.89), (0, 0, 0, 1, 0), (-0.89, 0, 0, 0, 0.45))
    }
    else {
      std::cout << " Error occured during SVD" << std::endl;
    }

    errorLogger.throwError("linearAlgebraTest", errMessage, INFO);
    oOut << errMessage << std::endl;

    errMessage = " Singular Value Decomposition test 3 ";

    //std::cout << "\n\n SVD 3 " << std::endl;
    ublas::matrix<double, ublas::column_major> A3(3, 3);
    // http://osdir.com/ml/lib.boost.ublas/2006-09/msg00020.html
    A3(0,0) = 1.00000; A3(0,1) = 0.50000; A3(0,2) = 0.33333;
    A3(1,0) = 0.50000; A3(1,1) = 0.33333; A3(1,2) = 0.25000;
    A3(2,0) = 0.33333; A3(2,1) = 0.25000; A3(2,2) = 0.20000;

    ublas::matrix<double, ublas::column_major> U3(3, 3);
    ublas::vector<double> S3(3);
    ublas::matrix<double, ublas::column_major> VT3(3, 3);

    //std::cout << " A \n" << A3 << "\n" << std::endl;

    std::stringstream A3_ss;
    A3_ss << A3;
    errMessage += "\n ### Start A = " + A3_ss.str();

    r = svd(A3, S3, U3, VT3);
    if (r == 0) {
      //std::cout << " A \n" << A3 << "\n" << std::endl;
      //std::cout << " U \n" << U3 << "\n" << std::endl;
      //std::cout << " S \n " << S3 << "\n" << std::endl;
      //std::cout << " VT \n " << VT3 << "\n" << std::endl;

      std::stringstream A32_ss;
      A32_ss << A3;
      errMessage += "\n ### A = " + A32_ss.str();

      std::stringstream U3_ss;
      U3_ss << U3;
      errMessage += "\n ### U = " + U3_ss.str();

      std::stringstream S3_ss;
      S3_ss << S3;
      errMessage += "\n ### s = " + S3_ss.str();

      std::stringstream VT3_ss;
      VT3_ss << VT3;
      errMessage += "\n ### VT = " + VT3_ss.str();

      // Expected results
      // U = ((-0.827046, -0.459864, -0.323297), (0.547445, -0.528278, -0.64902), (0.12767, -0.713756, 0.68866))
      // S = (1.40832, 0.122329, 0.00268506)
      // VT = ((-0.827046 -0.459864 -0.323297),(0.547445, -0.528278, -0.64902), (0.12767, -0.713756, 0.68866))
    }

    errorLogger.throwError("linearAlgebraTest", errMessage, INFO);
    oOut << errMessage << std::endl;
*/
    time (&endTime);
    int diffTime = (int) difftime(endTime, startTime);

    errorMessage = " Exited Normally After " + int2String(diffTime) + " Seconds ";
    errorLogger.throwError("linearAlgebraTest", errorMessage, 4);

    oLog.close();
    oOut.close();

    return 0;
}
