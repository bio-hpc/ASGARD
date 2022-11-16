/*!
   \file statsTest.cpp
   \brief Tests the stats functionality in MTK++
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

#include "Statistics/BaseStats.h"

#include "Statistics/ols.h"
#include "Statistics/pls.h"
#include "Statistics/pca.h"
#include "Statistics/sheet.h"
#include "Statistics/table.h"

#include "Utils/vector3d.h"
#include "Utils/deleteItem.h"

#include "Parsers/dMParser.h"

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
   \brief Tests MTK++'s stats functionality
   \param argc
   \param argv
   \return success
*/
int main (int argc, char **argv)
{
    std::string prog_name = "statsTest";
    std::vector<std::string> authors;
    std::string author = "Martin B. Peters";
    authors.push_back(author);

    // 1. CREATE AN OBJECT
    commLineOptions *clo = new commLineOptions();

    // 2. SET PREFERENCES
    clo->noUsage();

    // 3. SET THE USAGE/HELP
    clo->addUsage( "  statsTest                                          \n" );
    clo->addUsage( "    usage:  statsTest [flags] [options]              \n" );
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

    std::string logFile = "statsTest.log";
    std::string outFile = "statsTest.out";
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

    // Setup
    if (!getenv("AMBERHOME")) {
      std::cout << " Set the AMBERHOME environment variables " << std::endl;
      oLog.close();
      oOut.close();
      exit(0);
    }

    std::string AMBERHOME = getenv("AMBERHOME");
    if (AMBERHOME == "") {
      std::cout << " Set the AMBERHOME environment variables " << std::endl;
      oLog.close();
      oOut.close();
      exit(0);
    }

/*
    Eigen::MatrixXd A(5, 5);

    A(0, 0) =  1.96; A(0, 1) =  -6.49; A(0, 2) = -0.47; A(0, 3) = -7.20; A(0, 4) = -0.65;
    A(1, 0) = -6.49; A(1, 1) =   3.80; A(1, 2) = -6.39; A(1, 3) =  1.50; A(1, 4) = -6.34;
    A(2, 0) = -0.47; A(2, 1) =  -6.39; A(2, 2) =  4.17; A(2, 3) = -1.51; A(2, 4) =  2.67;
    A(3, 0) = -7.20; A(3, 1) =   1.50; A(3, 2) = -1.51; A(3, 3) =  5.70; A(3, 4) =  1.80;
    A(4, 0) = -0.65; A(4, 1) =  -6.34; A(4, 2) =  2.67; A(4, 3) =  1.80; A(4, 4) = -7.10;

*/

    // create descriptor matrix parser
    dMParser* myParser = new dMParser();

    // sheet storage
    std::vector<sheet*> mySheets;

    //! sheet iterator
    typedef std::vector<sheet*>::iterator  sheetIterator;

    sheet* mySheet = new sheet();
    mySheet->setName("validate");
    mySheets.push_back(mySheet);
    myParser->import(mySheets[0], AMBERHOME + "/AmberTools/examples/mtkpp/stats/pca1/validateX.table");
    myParser->import(mySheets[0], AMBERHOME + "/AmberTools/examples/mtkpp/stats/pca1/validateY.table");

    table<double>* tableX = mySheets[0]->getTable("X");
    table<double>* tableY = mySheets[0]->getTable("Y");

    Eigen::Matrix<double, Dynamic, Dynamic> X = tableX->getMatrix();
    Eigen::Matrix<double, Dynamic, Dynamic> Y = tableY->getMatrix();

    tableX->print();
    tableY->print();

    // Y Predicted
    table<double>* tableYPredicted = mySheets[0]->addTable();
    tableYPredicted->setName("Y-Predicted");
    tableYPredicted->setSizes(Y.rows(), Y.cols());
    Eigen::Matrix<double, Dynamic, Dynamic>& YPredicted = tableYPredicted->getMatrix();
    YPredicted << 4.4,5.05,5.1,6.5,6.01,6.9,
                  8.5,5.4,8.6,8.3,9.6,8.2,9.3,
                  5.1,9.1,7.1,6.8,4.0,5.6,3.5,
                  10.1,7.9,5.0,2.3,4.5,6.2,6.9,
                  8.5,6.8,6.9,6.4,6.5,6.8,6.9,
                  6.3,8.5,6.5,7.1,8.2,5.6,6.1,
                  10.0,15.0,10.1,9.1,6.9,9.15,
                  10.1,11.1,8.2,12.6;
    tableYPredicted->print();

    // Y Residuals
    table<double>* tableYResiduals = mySheets[0]->addTable();
    tableYResiduals->setName("Y-Residuals");
    tableYResiduals->setSizes(Y.rows(), Y.cols());
    Eigen::Matrix<double, Dynamic, Dynamic>& YResiduals = tableYResiduals->getMatrix();

    // X Centers
    table<double>* tableXCenters = mySheets[0]->addTable();
    tableXCenters->setName("X-Centers");
    tableXCenters->setSizes(1, X.cols());
    Eigen::Matrix<double, Dynamic, Dynamic>& XCenters = tableXCenters->getMatrix();

    // X Column Centered
    table<double>* tableXColumnCentered = mySheets[0]->addTable();
    tableXColumnCentered->setName("X-ColumnCentered");
    tableXColumnCentered->setSizes(X.rows(), X.cols());
    Eigen::Matrix<double, Dynamic, Dynamic>& XColumnCentered = tableXColumnCentered->getMatrix();

    // X Row Centered
    //table<double>* tableXRowCentered = mySheets[0]->addTable();
    //tableXRowCentered->setName("X-RowCentered");
    //tableXRowCentered->setSizes(X.rows(), X.cols());
    //Eigen::Matrix<double, Dynamic, Dynamic>& XRowCentered = tableXRowCentered->getMatrix();

    // X Column Standard Deviations
    table<double>* tableXColumnStdDev = mySheets[0]->addTable();
    tableXColumnStdDev->setName("X-ColumnStdDev");
    tableXColumnStdDev->setSizes(1, X.cols());
    Eigen::Matrix<double, Dynamic, Dynamic>& XColumnStdDev = tableXColumnStdDev->getMatrix();

    // X Column ZScores
    table<double>* tableXColumnZscores = mySheets[0]->addTable();
    tableXColumnZscores->setName("X-ColumnZScores");
    tableXColumnZscores->setSizes(X.rows(), X.cols());
    Eigen::Matrix<double, Dynamic, Dynamic>& XColumnZScores = tableXColumnZscores->getMatrix();

    // X Autoscaled
    table<double>* tableXAutoScaled = mySheets[0]->addTable();
    tableXAutoScaled->setName("X-AutoScaled");
    tableXAutoScaled->setSizes(X.rows(), X.cols());
    Eigen::Matrix<double, Dynamic, Dynamic>& XAutoScaled = tableXAutoScaled->getMatrix();

    // X Covariance
    table<double>* tableXCovariance = mySheets[0]->addTable();
    tableXCovariance->setName("X-Covariance");
    tableXCovariance->setSizes(X.cols(), X.cols());
    Eigen::Matrix<double, Dynamic, Dynamic>& XCovariance = tableXCovariance->getMatrix();

    table<std::string>* stringTable = mySheets[0]->addStringTable();
    stringTable->setName("String Table");
    stringTable->setSizes(10, 1);
    stringTable->initialize("----");
    stringTable->setCellValue(1, 0, "hello");
    stringTable->setCellValue(2, 0, "world");
    stringTable->print();

    //
    // Start
    //
    BaseStats* bs = new BaseStats();

    double mC0 = bs->meanColumn(X, 0);
    std::cout << " meanColumn " << mC0 << std::endl; // -90.4582

    double mR0 = bs->meanRow(X, 0);
    std::cout << " meanRow " <<  mR0 << std::endl; // 233.814

    double sC0 = bs->sumColumn(X, 0);
    std::cout << " sumColumn " << sC0 << std::endl; // -4613.37

    double sR0 = bs->sumRow(X, 0);
    std::cout << " sumRow " << sR0 << std::endl; // 3273.40

    double maxC0 = bs->maxColumn(X, 0);
    std::cout << " maxColumn " << maxC0 << std::endl; // 208.0

    int r = 0;
    maxC0 = bs->maxColumn(X, 0, r);
    std::cout << " maxColumn " << maxC0 << " " << r << std::endl; // 208 50

    double maxR0 = bs->maxRow(X, 0);
    std::cout << " maxRow " << maxR0 << std::endl; // 909.780

    double minC0 = bs->minColumn(X, 0);
    std::cout << " minColumn " << minC0 << std::endl; // -515.000

    minC0 = bs->minColumn(X, 0, r);
    std::cout << " minColumn " << minC0 << " " << r << std::endl; // -515.000 47

    double minR0 = bs->minRow(X, 0);
    std::cout << " minRow " << minR0 << std::endl; // -54.4600

    double stdDevC0 = bs->standardDeviationColumn(X, 0);
    std::cout << " standardDeviationColumn " << stdDevC0 << std::endl; // 120.709

    double stdDevR0 = bs->standardDeviationRow(X, 0);
    std::cout << " standardDeviationRow " << stdDevR0 << std::endl; // 291.095

    double varC0 = bs->varianceColumn(X, 0);
    std::cout << " varianceColumn " << varC0 << std::endl; // 14570.7

    double varR0 = bs->varianceRow(X, 0);
    std::cout << " varianceRow " << varR0 << std::endl; // 84736.5

    double corrCoefC12 = bs->correlationCoefficientColumn(Y, 0, YPredicted, 0);
    std::cout << " correlationCoefficient " << corrCoefC12 << std::endl; // 0.907736

    double rSqrC12 = bs->rSquaredColumn(Y, 0, YPredicted, 0);
    std::cout << " rSquared " << rSqrC12 << std::endl; // 0.823984

    double adjRSqrC12 = bs->AdjustedRSquaredColumn(Y, 0, YPredicted, 0, 1);
    std::cout << " adjustedRSquared " << adjRSqrC12 << std::endl; // 0.807511

    double rmse = bs->RMSE(Y, 0, YPredicted, 0);
    std::cout << " rmse " << rmse << std::endl; // 0.980104

    double mse = bs->MSE(Y, 0, YPredicted, 0);
    std::cout << " mse " << mse << std::endl; // 0.960604

    double unSigErr = bs->UnsignedError(Y, 0, YPredicted, 0);
    std::cout << " unSignedError " << unSigErr << std::endl; // 0.622353

    double sigErr = bs->SignedError(Y, 0, YPredicted, 0);
    std::cout << " signedError " << sigErr << std::endl; // 0.0811765

    double sumSqDevColumn = bs->SumSquaredDeviationsColumn(Y, 0);
    std::cout << " SumSquaredDeviationsColumn " << sumSqDevColumn << std::endl; // 

    double sumSqRegColumn = bs->SumSquaredRegressionColumn(Y, YPredicted);
    std::cout << " SumSquaredRegressionColumn " << sumSqRegColumn << std::endl; // 

    double sumSqResErrorColumn = bs->SumSquaredResidualErrorsColumn(Y, 0, YPredicted, 0);
    std::cout << " SumSquaredResidualErrorsColumn " << sumSqResErrorColumn << std::endl; // 

    double sumSqResErrColumn = bs->SumSquaredResidualErrorsColumn(Y, YPredicted, YResiduals);
    std::cout << " SumSquaredResidualErrorsColumn " << sumSqResErrColumn << std::endl; // 
    tableYResiduals->print();

    int rValue = bs->getColumnCenters(X, XCenters);
    if (rValue == 0) {
      tableXCenters->print();
    }

    rValue = bs->getStdDevColumns(X, XColumnStdDev);
    if (rValue == 0) {
      tableXColumnStdDev->print();
    }

    bs->centerColumns(X, XColumnCentered);
    tableXColumnCentered->print();

    // todo
    //bs->centerRows(X, XRowCentered);
    //tableXRowCentered->print();

    bs->zScoreColumns(X, XColumnZScores);
    tableXColumnZscores->print();

    // todo
    //bs->zScoreRows(X, XRowZScores);
    //tableXRowZscores->print();

    rValue = bs->autoScale(X, XAutoScaled, XCenters, XColumnStdDev);
    if (rValue == 0) {
      tableXAutoScaled->print();
    }

    // todo
    //double         covarianceColumn(Eigen::Matrix<double, Dynamic, Dynamic>& m1, const int& i1,
    //                                Eigen::Matrix<double, Dynamic, Dynamic>& m2, const int& i2);
    //double         covarianceColumn(Eigen::Matrix<double, Dynamic, Dynamic>& m, const int& i, const int& j);

    bs->covarianceMatrix(X, XCovariance);
    tableXCovariance->print();

    // PCA
    sheet* pcaModelSheet = new sheet();
    pcaModelSheet->setName("PCA");
    mySheets.push_back(pcaModelSheet);
    pca* myPCA = new pca(tableX, pcaModelSheet);
    int f = myPCA->run(5);
    std::cout << " PCA " << f << std::endl;





    time (&endTime);
    int diffTime = (int) difftime(endTime, startTime);

    errorMessage = " Exited Normally After " + int2String(diffTime) + " Seconds ";
    errorLogger.throwError("statsTest", errorMessage, 4);

    oLog.close();
    oOut.close();

    delete bs;
    return 0;
}
