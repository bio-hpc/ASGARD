/*!
   \file stats.cpp

   \brief Simple stats program

   \author Martin B. Peters

   $Date: 2010/04/22 22:19:54 $
   $Revision: 1.7 $

   ----------------------------------------------------------------------------

   MTK++ - C++ package of modeling libraries.

   Copyright (C) 2005-2006  (see AUTHORS file for a list of contributors)

   ----------------------------------------------------------------------------
*/
#include "Utils/printHeader.h"

#include "Statistics/ols.h"
#include "Statistics/pls.h"
#include "Statistics/pca.h"
#include "Statistics/sheet.h"
#include "Statistics/table.h"

#include "Utils/vector3d.h"
#include "Utils/deleteItem.h"

#include "Parsers/dMParser.h"
#include "Parsers/StringManip.h"

#include "Parsers/commLineOptions.h"
#include "Parsers/inputParser.h"

// time
#include "time.h"

// errorHandler
#include "Log/errorHandler.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <vector>
#include <map>
#include <algorithm>

using namespace MTKpp;

int main (int argc, char *argv[])
{
    std::string prog_name = "stats";
    std::vector<std::string> authors;
    std::string author = "Martin B. Peters";
    authors.push_back(author);

    // 1. CREATE AN OBJECT
    commLineOptions *clo = new commLineOptions();

    // 2. SET PREFERENCES
    clo->noUsage();

    // 3. SET THE USAGE/HELP

    clo->addUsage( "  stats: Simple Stats Program                        \n");
    clo->addUsage( "    usage: stats [flags] [options]                   \n" );
    clo->addUsage( "  options:                                             " );
    clo->addUsage( "          -i input file                                " );
    clo->addUsage( "          -l log file                                \n" );
    clo->addUsage( "    flags:                                             " );
    clo->addUsage( "          -h help                                      " );

    // 4. SET THE OPTION STRINGS/CHARACTERS
    clo->setOption(  "input",    'i' );
    clo->setOption(  "log",      'l' );
    clo->setFlag  (  "help",     'h' );

    // 5. PROVIDE THE COMMANDLINE
    clo->processCommandArgs( argc, argv );

    clo->usageOn();

    // 6. GET THE VALUES
    if ( clo->getFlag( "help" ) || clo->getFlag( 'h' ) ) {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      exit(0);
    }

    std::string inputFile = "";
    std::string logFile = "stats.log";

    bool bInput = false;
    bool bLog = false;

    if ( clo->getValue( "i" ) != 0 ) {
      inputFile = clo->getValue( "i" );
      bInput = true;
    }
    else if ( clo->getValue( "input" ) != 0 ) {
      inputFile = clo->getValue( "input" );
      bInput = true;
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide an input file " << std::endl;
      exit(0);
    }

    if ( clo->getValue( "l" ) != 0 ) {
      logFile = clo->getValue( "l" );
      bLog = true;
    }
    else if ( clo->getValue( "log" ) != 0 ) {
      logFile = clo->getValue( "log" );
      bLog = true;
    }

    delete clo;

    // Open log file
    std::ofstream oLog;
    oLog.open(logFile.c_str());

    if (!oLog) {
      std::cout << "\nUNABLE TO OPEN LOG FILE"
                << "\nFILENAME = " << logFile << std::endl;
      exit(1);
    }

    // Set errorLog stream to the log file
    MTKpp::errorLogger.setStream(&oLog);

    // Print MTK++ copyright message
    printHeader(oLog, prog_name, authors);

    // Start & end time
    time_t startTime;
    time_t endTime;

    // Get start time
    time (&startTime);

    // Read input file
    std::vector<std::vector<std::string> > inputFileContents;
    if (bInput) {
      int failure = readInputFile(inputFile, inputFileContents);
      if (failure) {
        MTKpp::errorLogger.throwError("stats", "Failed to open input file ", 1);
        exit(0);
      }
    }

    // create descriptor matrix parser
    dMParser* myParser = new dMParser();

    // sheet storage
    std::vector<sheet*> mySheets;

    //! sheet iterator
    typedef std::vector<sheet*>::iterator  sheetIterator;

    for (unsigned int i = 0; i < inputFileContents.size(); i++) {
      if (inputFileContents[i][0] == "quit") {
        /*!
           \code
            Function: quit

            Description: Exits program

            syntax: quit

           \endcode
        */
        time (&endTime);
        int diffTime = (int) difftime(endTime, startTime);
        std::string errMessage = " stats Exited Normally After " 
                    + int2String(diffTime) + " Seconds ";
        MTKpp::errorLogger.throwError("stats::quit", errMessage, 0);
        oLog.close();
        exit(0);
      }

      else if (inputFileContents[i][0] == "createSheet") {
        /*!
           \code
           Function: createSheet

           Description: Create sheet

           syntax: createSheet sheetName
           \endcode
        */
        if (inputFileContents[i].size() == 2) {
          // create sheet
          sheet* mySheet = new sheet();
          mySheet->setName(inputFileContents[i][1]);
          mySheets.push_back(mySheet);
          //std::cout << " CREATESHEET " << inputFileContents[i][1] << std::endl;
        }
      }

      else if (inputFileContents[i][0] == "import") {
        /*!
           \code
           Function: import

           Description: Import text file into a particular sheet

           syntax: import t.txt into sheetName
           \endcode
        */
        bool success = false;
        if (inputFileContents[i].size() == 4) {
          for (unsigned int s = 0; s < mySheets.size(); s++) {
            if (mySheets[s]->getName() == inputFileContents[i][3]) {
              myParser->import(mySheets[s], inputFileContents[i][1]);
              success = true;
              break;
            }
          }
        }
        else {
          std::cout << " Incorrect use of import " << std::endl;
          std::cout << " e.g. import t.txt into sheetName " << std::endl;
          exit(0);
        }
        if (!success) {
          std::cout << " Incorrect use of import " << std::endl;
          exit(0);
        }
      }

      else if (inputFileContents[i][0] == "read") {
        /*!
           \code
           Function: read

           Description: read xml file

           syntax: read t.xml
           \endcode
        */
        if (inputFileContents[i].size() == 2) {
          sheet* pSheet = new sheet();
          pSheet->setName(inputFileContents[i][1]);
          myParser->read(pSheet, inputFileContents[i][1]);
          mySheets.push_back(pSheet);
        }
      }

      else if (inputFileContents[i][0] == "write") {
        /*!
           \code
           Function: write

           Description: write xml file

           syntax: read t.xml sheetName
           \endcode
        */
        if (inputFileContents[i].size() == 3) {
          for (unsigned int s = 0; s < mySheets.size(); s++) {
            if (mySheets[s]->getName() == inputFileContents[i][2]) {
              myParser->write(mySheets[s], inputFileContents[i][1]);
              break;
            }
          }
        }
      }

      else if (inputFileContents[i][0] == "print") {
        /*!
           \code
           Function: print

           Description: print

           syntax: print sheetName/TableName
           \endcode
        */
        if (inputFileContents[i].size() == 2) {
          std::vector<std::string> words;
          splitString(inputFileContents[i][1], "/", words, 0);
          if (words.size() != 2) {
            std::cout << " Incorrect use of print " << std::endl;
            exit(1);
          }
          for (unsigned int j = 0; j < mySheets.size(); j++) {
            if (mySheets[j]->getName() == words[0]) {
              table<double>* myTable = mySheets[j]->getTable(words[1]);
              if (myTable) {
                myTable->print();
              }
              table<int>* myIntTable = mySheets[j]->getIntTable(words[1]);
              if (myIntTable) {
                //myIntTable->print();
              }
              break;
            }
          }
        }
      }
/*
/////////////////////////////
      else if (inputFileContents[i][0] == "ols") {
        //
           \code
           Function: ols

           Description: Ordinary Least Squares

           syntax: ols sheetName sheetName/Y sheetName/X
           \endcode
        //
        bool bError = false;
        if (inputFileContents[i].size() == 4) {
          std::vector<std::string> words;
          splitString(inputFileContents[i][2], "/", words, 0);
          if (words.size() != 2) {
            std::cout << " Incorrect use of ols " << std::endl;
            exit(1);
          }

          std::vector<std::string> words2;
          splitString(inputFileContents[i][3], "/", words2, 0);
          if (words2.size() != 2) {
            std::cout << " Incorrect use of ols " << std::endl;
            exit(1);
          }

          sheet* outModelSheet = new sheet();
          outModelSheet->setName(inputFileContents[i][1]);
          mySheets.push_back(outModelSheet);

          sheet* ySheet = 0;
          sheet* xSheet = 0;
          table<double>* yTable = 0;
          table<double>* xTable = 0;

          for (unsigned int j = 0; j < mySheets.size(); j++) {
            if (mySheets[j]->getName() == words[0]) {
              ySheet = mySheets[j];
            }
            if (mySheets[j]->getName() == words2[0]) {
              xSheet = mySheets[j];
            }
          }

          if (ySheet and xSheet and outModelSheet) {
            yTable = ySheet->getTable(words[1]);
            xTable = xSheet->getTable(words2[1]);
          }
          else {
            bError = true;
          }

          if (!bError) {
            // create pls
            bool olsError = false;
            ols* myOLS = new ols(yTable, xTable, outModelSheet, olsError);
            if (olsError) {
              std::cout << " Error running OLS ... exiting " << std::endl;
              exit(0);
            }

            myOLS->run(olsError);
            if (olsError) {
              std::cout << " Error running OLS ... exiting " << std::endl;
              exit(0);
            }

            table<double>* yPred = outModelSheet->getTable("Y Pred");
            if (yPred) {
              yPred->print();
            }

            table<double>* coefficients = outModelSheet->getTable("Coefficients");
            if (coefficients) {
              coefficients->print();
            }

            table<double>* residuals = outModelSheet->getTable("Residuals");
            if (residuals) {
              residuals->print();
            }

            table<double>* lx = outModelSheet->getTable("X");
            if (lx) {
              lx->print();
            }

          }
        }
      }
*/
/////////////////////////

      else if (inputFileContents[i][0] == "pca") {
        /*!
           \code
           Function: pca

           Description: Principal Component Analysis

           syntax: pca sheetname sheetName/tableName
           \endcode
        */
        if (inputFileContents[i].size() == 3) {
          std::vector<std::string> words;
          splitString(inputFileContents[i][2], "/", words, 0);
          if (words.size() != 2) {
            std::cout << " Incorrect use of pca " << std::endl;
            exit(1);
          }

          sheet* outModelSheet = new sheet();
          outModelSheet->setName(inputFileContents[i][1]);
          mySheets.push_back(outModelSheet);

          for (unsigned int j = 0; j < mySheets.size(); j++) {
            if (mySheets[j]->getName() == words[0]) {
              table<double>* myTable = mySheets[j]->getTable(words[1]);

              if (outModelSheet and myTable) {
                // Create a PCA object
                pca* myPCA = new pca(myTable, outModelSheet);

                int f = myPCA->run(5);
                if (f) {
                  MTKpp::errorLogger.throwError("stats::pca", " Failed ", 1);
                  exit(0);
                }
                delete myPCA;
              }
            }
          }
        }
      }

      else if (inputFileContents[i][0] == "pls") {
        /*!
           \code
           Function: pls

           Description: Partial Least Squares

           syntax: pls sheetName sheetName/Y sheetName/X nlvs cv
           \endcode
        */
        bool bError = false;
        if (inputFileContents[i].size() > 3) {

          std::vector<std::string> words;
          splitString(inputFileContents[i][2], "/", words, 0);
          if (words.size() != 2) {
            std::cout << " Incorrect use of pls " << std::endl;
            exit(1);
          }

          std::vector<std::string> words2;
          splitString(inputFileContents[i][3], "/", words2, 0);
          if (words2.size() != 2) {
            std::cout << " Incorrect use of pls " << std::endl;
            exit(1);
          }

          sheet* outModelSheet = new sheet();
          outModelSheet->setName(inputFileContents[i][1]);
          mySheets.push_back(outModelSheet);

          sheet* ySheet = 0;
          sheet* xSheet = 0;
          table<double>* yTable = 0;
          table<double>* xTable = 0;

          for (unsigned int j = 0; j < mySheets.size(); j++) {
            if (mySheets[j]->getName() == words[0]) {
              ySheet = mySheets[j];
            }
            if (mySheets[j]->getName() == words2[0]) {
              xSheet = mySheets[j];
            }
          }

          if (ySheet and xSheet and outModelSheet) {
            yTable = ySheet->getTable(words[1]);
            xTable = xSheet->getTable(words2[1]);
          }
          else {
            bError = true;
          }

          if (!bError) {

            int nlvs = 0;
            if (inputFileContents[i].size() > 4) {
              nlvs = string2Int(inputFileContents[i][4]);
            }
            else {
              nlvs = std::min(xTable->getNumRows(), xTable->getNumColumns());
            }

            std::string cv = "NONE";
            if (inputFileContents[i].size() > 5) {
              cv = inputFileContents[i][5];
            }

            // create pls
            bool plsError = false;
            pls* myPLS = new pls(yTable, xTable, "KERNELPLS", nlvs, cv, outModelSheet, plsError);
            myPLS->setNTEST(3);
            myPLS->setNITER(2); //17
            myPLS->setSEED(5);
            myPLS->run(plsError);

            if (plsError) {
              std::cout << " Error running KERNELPLS ... exiting " << std::endl;
              exit(0);
            }

            table<double>* regCoeff = outModelSheet->getTable("Regression Coefficients");
            if (regCoeff) {
              regCoeff->print();
            }

            table<double>* yPred = outModelSheet->getTable("Y-Pred");
            if (yPred) {
              yPred->print();
            }

            table<double>* r2 = outModelSheet->getTable("R2");
            if (r2) {
              r2->print();
            }

            table<double>* aue = outModelSheet->getTable("Average Unsigned Error");
            if (aue) {
              aue->print();
            }

            table<double>* rmse = outModelSheet->getTable("Root Mean Squared Error");
            if (rmse) {
              rmse->print();
            }

            delete myPLS;
          }
          else {
            std::cout << " Incorrect use of pls " << std::endl;
            exit(1);
          }
        }
      }

      else {
        std::cout << " unknown command: " << inputFileContents[i][0] << std::endl;
      }
    }

    std::for_each(mySheets.begin(), mySheets.end(), delete_element<sheet*>);
    //for (sheetIterator s = mySheets.begin(); s != mySheets.end(); s++) {
    //  sheet* pSheet = *s;
    //  delete pSheet;
    //}
    mySheets.clear();
}
