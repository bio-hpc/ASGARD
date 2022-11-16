/*!
   \file linearAlgebraTest.cpp
   \brief Tests the eigen functionality in MTK++
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

#include "Utils/diagonalize_eigen.h"
#include "Parsers/StringManip.h"

#include "Parsers/commLineOptions.h"

// Eigen
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
using namespace Eigen;

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
   \brief Tests MTK++'s/eigen functionality
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
    clo->addUsage( "          -l log file                                  " );
    clo->addUsage( "          -o out file                                  " );
    clo->addUsage( "          -b build directory                         \n" );
    clo->addUsage( "    flags:                                             " );
    clo->addUsage( "          -h help                                    \n" );
    clo->addUsage( "   author:                                             " );
    clo->addUsage( "          Martin B. Peters (c) 2011                    " );

    // 4. SET THE OPTION STRINGS/CHARACTERS
    clo->setOption( "logFile", 'l' );
    clo->setOption( "outFile", 'o' );
    clo->setOption( "buildDir", 'b' );
    clo->setFlag( "help",    'h' );

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

    Eigen::MatrixXd A(5, 5);

    A(0,0) =  1.96; A(0,1) = -6.49; A(0,2) = -0.47; A(0,3) = -7.20; A(0,4) = -0.65;
    A(1,0) = -6.49; A(1,1) =  3.80; A(1,2) = -6.39; A(1,3) =  1.50; A(1,4) = -6.34;
    A(2,0) = -0.47; A(2,1) = -6.39; A(2,2) =  4.17; A(2,3) = -1.51; A(2,4) =  2.67;
    A(3,0) = -7.20; A(3,1) =  1.50; A(3,2) = -1.51; A(3,3) =  5.70; A(3,4) =  1.80;
    A(4,0) = -0.65; A(4,1) = -6.34; A(4,2) =  2.67; A(4,3) =  1.80; A(4,4) = -7.10;

    std::string errMessage = " MTKpp/Eigen3 Test\n";

    std::stringstream A_ss;
    A_ss << A;
    errMessage += " ### A = " + A_ss.str() +"\n";

    oOut << errMessage << std::endl;
    errorLogger.throwError("linearAlgebraTest", errMessage, INFO);

    SelfAdjointEigenSolver<MatrixXd> eigensolver(A);

    MatrixXd evectors = eigensolver.eigenvectors();
    VectorXd evalues = eigensolver.eigenvalues();

    std::stringstream Aevalues_ss;
    Aevalues_ss << evalues;
    errMessage = " ### Eigenvalues = \n" + Aevalues_ss.str() +"\n";

    oOut << errMessage << std::endl;
    errorLogger.throwError("linearAlgebraTest", errMessage, INFO);

    std::stringstream Aevectors_ss;
    Aevectors_ss << evectors;
    errMessage = " ### Eigenvectors = \n" + Aevectors_ss.str() +"\n";

    oOut << errMessage << std::endl;
    errorLogger.throwError("linearAlgebraTest", errMessage, INFO);

    eigenValueSort(evectors, evalues, 1);
    //cout << "V * D * V^(-1) = " << endl 
    //     << evectors * evalues * evectors.inverse() << endl;

    std::stringstream Aevalues2_ss;
    Aevalues2_ss << evalues;
    errMessage = " ### Sorted Eigenvalues = \n" + Aevalues2_ss.str() +"\n";

    oOut << errMessage << std::endl;
    errorLogger.throwError("linearAlgebraTest", errMessage, INFO);

    std::stringstream Aevectors2_ss;
    Aevectors2_ss << evectors;
    errMessage = " ### Eigenvectors = \n" + Aevectors2_ss.str() +"\n";

    oOut << errMessage << std::endl;
    errorLogger.throwError("linearAlgebraTest", errMessage, INFO);

/*
    // 
    // From: http://eigen.tuxfamily.org/dox-devel/TutorialLinearAlgebra.html
    //

    // Basic linear solving
    Matrix3f A2;
    Vector3f b2;
    A2 << 1,2,3,  4,5,6,  7,8,10;
    b2 << 3, 3, 4;
    cout << "Here is the matrix A:\n" << A2 << endl;
    cout << "Here is the vector b:\n" << b2 << endl;
    Vector3f x2 = A2.colPivHouseholderQr().solve(b2);
    cout << "The solution is:\n" << x2 << endl;

    //  LDLT decomposition
    Matrix2f A3, b3;
    A3 << 2, -1, -1, 3;
    b3 << 1, 2, 3, 1;
    cout << "Here is the matrix A:\n" << A3 << endl;
    cout << "Here is the right hand side b:\n" << b3 << endl;
    Matrix2f x3 = A3.ldlt().solve(b3);
    cout << "The solution is:\n" << x3 << endl;

    // Checking if a solution really exists
    MatrixXd A4 = MatrixXd::Random(100,100);
    MatrixXd b4 = MatrixXd::Random(100,50);
    MatrixXd x4 = A4.fullPivLu().solve(b4);
    double relative_error4 = (A4*x4 - b4).norm() / b4.norm(); // norm() is L2 norm
    cout << "The relative error is:\n" << relative_error4 << endl;

    // Computing eigenvalues and eigenvectors
    Matrix2f A5;
    A5 << 1, 2, 2, 3;
    cout << "Here is the matrix A:\n" << A5 << endl;
    SelfAdjointEigenSolver<Matrix2f> eigensolver5(A5);
    cout << "The eigenvalues of A are:\n" << eigensolver5.eigenvalues() << endl;
    cout << "Here's a matrix whose columns are eigenvectors of A "
         << "corresponding to these eigenvalues:\n"
         << eigensolver5.eigenvectors() << endl;

    // Computing inverse and determinant
    Matrix3f A6;
    A6 << 1, 2, 1,
          2, 1, 0,
         -1, 1, 2;
    cout << "Here is the matrix A:\n" << A6 << endl;
    cout << "The determinant of A is " << A6.determinant() << endl;
    cout << "The inverse of A is:\n" << A6.inverse() << endl;

    // Least squares solving 
    // http://www.ces.clemson.edu/~petersj/Agents/MatLabNA/MatLabNA004.html
    Eigen::MatrixXf A7(5, 3);
    A7 << 1 ,1, 1, 2, -1, 2, -1, 4, 3, 4, 2, 1, 3, -3, 4;
    Eigen::VectorXf b7(5);
    b7 << 1, 2, -1, 4, 8;

    cout << "\n Least squares solving" <<endl;
    cout << "Here is the matrix A:\n" << A7 << endl;
    cout << "Here is the right hand side b:\n" << b7 << endl;

    JacobiSVD<MatrixXf> svd(A7, ComputeThinU | ComputeThinV);

    cout << "Its singular values are:" << endl << svd.singularValues() << endl;
    cout << "Its left singular vectors are the columns of the thin U matrix:" << endl << svd.matrixU() << endl;
    cout << "Its right singular vectors are the columns of the thin V matrix:" << endl << svd.matrixV() << endl;

    cout << "A least-squares solution of m*x = rhs is:" << endl << svd.solve(b7) << endl;
    cout << " e " << endl << A7 * svd.solve(b7) - b7;

    // Separating the computation from the construction
    Matrix2f A8, b8;
    LLT<Matrix2f> llt8;
    A8 << 2, -1, -1, 3;
    b8 << 1, 2, 3, 1;
    cout << "Here is the matrix A:\n" << A8 << endl;
    cout << "Here is the right hand side b:\n" << b8 << endl;
    cout << "Computing LLT decomposition..." << endl;
    llt8.compute(A8);
    cout << "The solution is:\n" << llt8.solve(b8) << endl;
    A8(1,1)++;
    cout << "The matrix A is now:\n" << A8 << endl;
    cout << "Computing LLT decomposition..." << endl;
    llt8.compute(A8);
    cout << "The solution is now:\n" << llt8.solve(b8) << endl;

    // Rank-revealing decompositions
    Matrix3f A9;
    A9 << 1, 2, 5,
          2, 1, 4,
          3, 0, 3;
    cout << "Here is the matrix A:\n" << A9 << endl;
    FullPivLU<Matrix3f> lu_decomp9(A9);
    cout << "The rank of A is " << lu_decomp9.rank() << endl;
    cout << "Here is a matrix whose columns form a basis of the null-space of A:\n"
         << lu_decomp9.kernel() << endl;
    cout << "Here is a matrix whose columns form a basis of the column-space of A:\n"
         << lu_decomp9.image(A9) << endl; // yes, have to pass the original A

    Matrix2d A10;
    A10 << 2, 1,
           2, 0.9999999999;
    FullPivLU<Matrix2d> lu10(A10);
    cout << "By default, the rank of A is found to be " << lu10.rank() << endl;
    lu10.setThreshold(1e-5);
    cout << "With threshold 1e-5, the rank of A is found to be " << lu10.rank() << endl;
*/
    time (&endTime);
    int diffTime = (int) difftime(endTime, startTime);

    errorMessage = " Exited Normally After " + int2String(diffTime) + " Seconds ";
    errorLogger.throwError("linearAlgebraTest", errorMessage, 4);

    oLog.close();
    oOut.close();

    return 0;
}
