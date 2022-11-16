/*!
   \file pdbSearcher.cpp

   \brief Search a local copy of the pdb for info, assuming that all file are *.Z

   \author Martin B. Peters

   $Date: 2010/05/04 20:04:07 $
   $Revision: 1.14 $

   ----------------------------------------------------------------------------

   MTK++ - C++ package of modeling libraries.

   Copyright (C) 2005-2006  (see AUTHORS file for a list of contributors)

   ----------------------------------------------------------------------------
*/
#include "Utils/constants.h"
#include "Utils/printHeader.h"

// - MOLECULE
#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/submolecule.h"
#include "Molecule/atom.h"
#include "Molecule/element.h"
#include "Molecule/connections.h"
#include "Molecule/stdLibrary.h"
#include "Molecule/metalCenter.h"

// - MTK++ INCLUDE
#include "Utils/vector3d.h"
#include "Utils/constants.h"
#include "Utils/diagonalize.h"

#include "Statistics/sheet.h"
#include "Statistics/table.h"

// - PARSERS
#include "Parsers/elementParser.h"
#include "Parsers/pdbParser.h"
#include "Parsers/StringManip.h"
#include "Parsers/commLineOptions.h"
#include "Parsers/inputParser.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <map>
#include <algorithm>
#include <math.h>
#include <cmath>

using namespace MTKpp;

//! string - double map for 
typedef std::map<std::string, double>::iterator strDbMapIterator;

//! string - int map for 
typedef std::map<std::string, int>::iterator strIntMapIterator;

//! string - string map for 
typedef std::map<std::string, std::string>::iterator strStrMapIterator;

//! int vector for 
typedef std::vector<unsigned int>::iterator intVectorIterator;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*!
   \brief Uncompress a .Z file using zcat
*/
int uncompressFile(std::string filename, std::string newFileName)
{
    int irunUnCompress = -1;
    std::string runUnCompress =  "zcat " + filename + " > " + newFileName;
    irunUnCompress = system(runUnCompress.c_str());
    return irunUnCompress;
}

/*!
   \brief 
*/
/*!
   \brief Searches a local copy of the Protein Data Bank
   \param argc
   \param argv
   \return success
*/
int main (int argc, char **argv)
{
    std::string prog_name = "pdbSearcher";
    std::vector<std::string> authors;
    std::string author = "Martin B. Peters";
    authors.push_back(author);

    // 1. CREATE AN OBJECT
    commLineOptions *clo = new commLineOptions();

    // 2. SET PREFERENCES
    clo->noUsage();

    // 3. SET THE USAGE/HELP
    clo->addUsage( "  pdbSearcher: Searches a local copy of the            " );
    clo->addUsage( "                   Protein Data Bank                 \n" );
    clo->addUsage( "    usage:  pdbSearcher [flags] [options]            \n" );
    clo->addUsage( "  options:                                             " );
    clo->addUsage( "          -i input file                                " );
    clo->addUsage( "          -l log file                                \n" );
    clo->addUsage( "    flags:                                             " );
    clo->addUsage( "          -h help                                    \n" );

    // 4. SET THE OPTION STRINGS/CHARACTERS
    clo->setOption( "inFile",    'i' );
    clo->setOption( "logFile",   'l' );
    clo->setFlag  ( "help",      'h' );

    // 5. PROVIDE THE COMMANDLINE
    clo->processCommandArgs( argc, argv );

    clo->usageOn();

    // 6. GET THE VALUES
    if ( clo->getFlag( "help" ) || clo->getFlag( 'h' ) ) {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      return 0;
    }

    std::string inputFile  = "";
    std::string logFile = "";
    bool bInput = false;
    bool bLog = false;

    if ( clo->getValue( "i" ) != 0 ) {
      inputFile = clo->getValue( "i" );
      bInput = true;
    }
    else if ( clo->getValue( "inFile" ) != 0 ) {
      inputFile =  clo->getValue( "inFile" );
      bInput = true;
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide an input file " << std::endl;
      return 0;
    }

    if ( clo->getValue( "l" ) != 0 ) {
      logFile = clo->getValue( "l" );
      bLog = true;
    }
    else if ( clo->getValue( "logFile" ) != 0 ) {
      logFile =  clo->getValue( "logFile" );
      bLog = true;
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide an output file " << std::endl;
      return 0;
    }

    // 7. DONE
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
    time (&endTime);

    int failure = 0;
    std::string errorMessage = "";

    // 
    std::map<std::string, unsigned int> numberOfMetals;

    // 
    std::vector<std::vector<std::string> > inputFileContents;
    failure = readInputFile(inputFile, inputFileContents);
    if (failure) {
      MTKpp::errorLogger.throwError("pdbSearcher", "Error reading input file ", 1);
      exit(0);
    }

    std::string AMBERHOME = getenv("AMBERHOME");
    collection* pCollection = new collection();
    pCollection->addParameters();
    pCollection->addStdLibrary();

    // READ ELEMENT DATA
    elementParser* pElementParser = new elementParser(pCollection->pElements);
    std::string elementXmlFile = AMBERHOME + "/dat/mtkpp/elements.xml";
    pElementParser->Read(elementXmlFile);

    // CREATE PARSERS
    pdbParser* pPdbParser = 0;
    pPdbParser = new pdbParser();

    std::map<std::string, std::vector<std::string> > pdbLists;
    std::vector<std::string> expTechs;
    double resolution = 0.0;
    bool bResolution = false;
    bool bWriteMetalEnvironment = false;
    bool bMonoMetal = true;

    for (unsigned int i = 0; i < inputFileContents.size(); i++) {
      if (inputFileContents[i][0] == "quit") {
        /*!
           \code
            Function: quit

            Description: Exits program

            syntax: quit

           \endcode
        */
        std::string removeFile =  "rm -rf pdbSearcherTmpFile ";
        int iR = system(removeFile.c_str());
        if (iR != 0 ) {
          std::cout << " pdbSearcher:Error executing remove " << std::endl;
          exit(1);
        }
        exit(0);
      }

      else if (inputFileContents[i][0] == "expTechniques") {
        /*!
           \code
            Function: expTechniques

            Description: 

            syntax: expTechniques X-RAY NMR

           \endcode
        */
        for (unsigned int j = 1; j < inputFileContents[i].size(); j++) {
          expTechs.push_back(inputFileContents[i][j]);
        }
      }

      else if (inputFileContents[i][0] == "resolution") {
        /*!
           \code
            Function: resolution

            Description:

            syntax: resolution 2.5

           \endcode
        */
        if (inputFileContents[i].size() != 2) {
          std::cout << " Incorrect use of resolution " << std::endl;
          exit(1);
        }
        else {
          resolution = strtod(inputFileContents[i][1].c_str(),0);
          bResolution = true;
        }
      }

      else if (inputFileContents[i][0] == "monoMetal") {
        /*!
           \code
            Function: monoMetal

            Description:

            syntax: monoMetal 1

           \endcode
        */
        if (inputFileContents[i].size() != 2) {
          std::cout << " Incorrect use of monoMetal " << std::endl;
          exit(1);
        }
        else {
          int s = string2Int(inputFileContents[i][1]);
          if (s) {
            bMonoMetal = true;
          }
          else {
            bMonoMetal = false;
          }
        }
      }

      else if (inputFileContents[i][0] == "readPDBList") {
        /*!
           \code
           Function: readPDBList

           Description: Read list of pdb files

           syntax: readPDBList list_name file_name
           \endcode
        */
        if (inputFileContents[i].size() != 3) {
          std::cout << " Incorrect use of readPDBList " << std::endl;
          exit(1);
        }
        else {
          std::vector<std::string> pdbList;
          int f = readListFile(inputFileContents[i][2], pdbList);
          if (f) {
            std::cout << " Error in readPDBList ... exiting " << std::endl;
            exit(1);
          }
          pdbLists[inputFileContents[i][1]] = pdbList;
        }
      }

      else if (inputFileContents[i][0] == "writePDBList") {
        /*!
           \code
           Function: writePDBList

           Description: Write list of pdb files

           syntax: writePDBList list_name file_name
           \endcode
        */
        if (inputFileContents[i].size() != 3) {
          std::cout << " Incorrect use of writePDBList " << std::endl;
          exit(1);
        }
        else {
          std::vector<std::string> cList = pdbLists[inputFileContents[i][1]];
          if (!cList.empty()) {
            std::ofstream opdb;
            opdb.open(inputFileContents[i][2].c_str());

            if (!opdb) {
              std::cout << "\nUNABLE TO OPEN FILE"
                        << "\nFILENAME = " << inputFileContents[i][2]
                        << "\nEXITING...\n" << std::endl;
              exit(0);
            }
            for (unsigned int p = 0; p < cList.size(); p++) {
              opdb << cList[p] << std::endl;
            }
            opdb.close();
          }
        }
      }

      else if (inputFileContents[i][0] == "hasMetal") {
        /*!
           \code
           Function: hasMetal

           Description: Does pdb file contain a metal

           syntax: hasMetal list_name new_list_name metal
           \endcode
        */
        if (inputFileContents[i].size() != 4) {
          std::cout << " Incorrect use of hasMetal " << std::endl;
          exit(1);
        }
        else {
          // START TIME
          //time (&startTime);
          std::string metal = toUpper(inputFileContents[i][3]);
          std::vector<std::string> cList1 = pdbLists[inputFileContents[i][1]];
          if (!cList1.empty()) {
            std::vector<std::string> cList2;
            for (unsigned int p = 0; p < cList1.size(); p++) {
              if (extName(cList1[p]) == "Z") {
                int irunUnCompress = uncompressFile(cList1[p], "pdbSearcherTmpFile");
                if (irunUnCompress != 0 ) {
                  std::cout << " pdbSearcher:Error executing zcat for " << cList1[p] << std::endl;
                  continue;
                }
                else {
                  std::ifstream ipdb;
                  ipdb.open("pdbSearcherTmpFile");

                  if (!ipdb) {
                    std::cout << "\nUNABLE TO OPEN PDB FILE"
                              << "\nFILENAME = " << baseName(cList1[p]) + ".pdb"
                              << "\nEXITING...\n" << std::endl;
                    continue;
                  }
                  std::string fileline;
                  while (ipdb) {
                    getline(ipdb,fileline);
                    if (fileline.substr(0,6) == "FORMUL") {
                      std::string::size_type loc = fileline.find(metal);
                      if (loc != std::string::npos) {
                        cList2.push_back(cList1[p]);
                      }
                    }
                  }
                  ipdb.close();
                }
              }
            }
            pdbLists[inputFileContents[i][2]] = cList2;
          }
        }
        // END TIME
        //time (&endTime);
        //int diffTime = (int) difftime(endTime, startTime);
        //std::cout << " pdbSearcher:hasMetal tool " << diffTime << " seconds." << std::endl;
      }

      else if (inputFileContents[i][0] == "getEnvironment") {
        /*!
           \code
           Function: getEnvironment

           Description: Determine ligating residues

           syntax: getEnvironment list_name output_file metal
           \endcode
        */
        if (inputFileContents[i].size() != 4) {
          std::cout << " Incorrect use of getEnvironment " << std::endl;
          exit(1);
        }
        else {
          // START TIME
          //time (&startTime);
          std::string metal = toUpper(inputFileContents[i][3]);

          std::vector<std::string> cList1 = pdbLists[inputFileContents[i][1]];
          if (!cList1.empty()) {

            std::string metalCenterLogFile = inputFileContents[i][2] + ".log.csv";
            std::string metalCenterInfoFile = inputFileContents[i][2] + ".info.csv";

            // Open log file
            std::ofstream oLog;
            oLog.open(metalCenterLogFile.c_str());
            if (!oLog) {
              std::cout << "\nUNABLE TO OPEN FILE"
                        << "\nFILENAME = " << inputFileContents[i][2] + ".log.csv"
                        << "\nEXITING...\n" << std::endl;
              exit(0);
            }
            oLog << "FILE_NAME,RESOLUTION,EXP_TECH,TOTAL_NUM_ATOMS," <<
                    "NUM_METAL_ATOMS,RES_NAME,RES_ID,AT_NAME,AT_ID," <<
                    "B_FACTOR,PRIM_SHELL,SEC_SHELL,GEOM," <<
                    "GEOM_RMS,ERROR" << std::endl;

            // Open info file
            std::ofstream oInfo;
            oInfo.open(metalCenterInfoFile.c_str());
            if (!oInfo) {
              std::cout << "\nUNABLE TO OPEN FILE"
                        << "\nFILENAME = " << inputFileContents[i][2] + ".info.csv"
                        << "\nEXITING...\n" << std::endl;
              exit(0);
            }
            oInfo << "FILE_NAME,METAL_RES_NAME,METAL_RES_ID,METAL_NAME,"
                  << "METAL_ID,RES_NAME,RES_ID,AT_NAME,AT_ID,DISTANCE,"
                  << "TYPE,B_FACTOR,PRIM_SHELL,SEC_SHELL,GEOM1,GEOM1_RMS,"
                  << "GEOM2,GEOM2_RMS"
                  << std::endl;

            for (unsigned int p = 0; p < cList1.size(); p++) {
              std::cout << " Doing " << cList1[p] << std::endl;
              if (extName(cList1[p]) == "Z") {
                int irunUnCompress = uncompressFile(cList1[p], "pdbSearcherTmpFile");
                if (irunUnCompress != 0 ) {
                  std::cout << " pdbSearcher:Error executing zcat for " << cList1[p] << std::endl;
                  continue;
                }
                pPdbParser->Read("pdbSearcherTmpFile", pCollection);
              }
              else if (extName(cList1[p]) == "pdb") {
                pPdbParser->Read(cList1[p], pCollection);
              }

              bool bOK = false;
              for (unsigned int u = 0; u < expTechs.size(); u++) {
                std::string::size_type loc = pPdbParser->itsPdbInfo->expTechnique.find(expTechs[u], 0 );
                if ( loc != std::string::npos ) {
                  bOK = true;
                }
              }
              if (!bOK) {
                pCollection->delAllMolecules();
                std::cout << " experimental techniques error " << std::endl;
                continue;
              }

              if (bResolution) {
                if (pPdbParser->itsPdbInfo->resolution > resolution) {
                  pCollection->delAllMolecules();
                  std::cout << " resolution error " << std::endl;
                  continue;
                }
              }

              // std::cout << " Check to see if collection has a metal atom " << std::endl;
              std::vector<metalCenter*> metCens;
              if (pCollection->hasMetal()) {
                // std::cout << " hasMetal, now findMetals " << std::endl;
                pCollection->findMetals();

                // std::cout << " determineMetalEnvironments " << std::endl;
                pCollection->determineMetalEnvironments();

                // std::cout << " getMetalCenters " << std::endl;
                metCens = pCollection->getMetalCenters();
                // std::cout << " Number of metal centers = " << metCens.size() << std::endl;

                numberOfMetals[cList1[p]] = metCens.size();
                for (unsigned int mc = 0; mc < metCens.size(); mc++) {
                  oLog << cList1[p] << ","
                       << pPdbParser->itsPdbInfo->resolution << ","
                       << pPdbParser->itsPdbInfo->expTechnique << ","
                       << pCollection->getNumAtoms() << ","
                       << numberOfMetals[cList1[p]] << ","
                       << metCens[mc]->getMetalAtom()->getParent()->getName() << ","
                       << metCens[mc]->getMetalAtom()->getParent()->getSubMolId() << ","
                       << metCens[mc]->getMetalAtom()->getName() << ","
                       << metCens[mc]->getMetalAtom()->getFileID() << ","
                       << metCens[mc]->getMetalAtom()->getTempFactor () << ","
                       << metCens[mc]->getPrimaryShell() << ","
                       << metCens[mc]->getSecondaryShell() << ","
                       << metCens[mc]->getGeometry() << ","
                       << metCens[mc]->getGeometryRMS() << ","
                       << metCens[mc]->getError()
                       << std::endl;
                  metCens[mc]->setPdb(baseName(cList1[p]));
                  metCens[mc]->print(oInfo);
                  if (bWriteMetalEnvironment) {
                    std::string metalEnvironmentFileName = baseName(cList1[p]) + "_MetalCenter.pdb";;
                    pPdbParser->Write(metalEnvironmentFileName, metCens[mc]);
                  }
                }
              }
              else {
                std::cout << " No metal found " << std::endl;
              }

              // Clear collection for next pdb file
              pCollection->clear();
            }
            oLog.close();
            oInfo.close();
          }
        }

        // END TIME
        //time (&endTime);
        //int diffTime = (int) difftime(endTime, startTime);
        //std::cout << " pdbSearcher:getEnvironment tool " << diffTime << " seconds." << std::endl;
      }

      else if (inputFileContents[i][0] == "writeMetalEnvironment") {
        /*!
           \code
            Function: writeMetalEnvironment

            Description: Turn on writing of metal environments

            syntax: writeMetalEnvironment 1

           \endcode
        */
        if (inputFileContents[i].size() != 2) {
          std::cout << " Incorrect use of writeMetalEnvironment " << std::endl;
          exit(1);
        }
        else {
          int s = string2Int(inputFileContents[i][1]);
          if (s) {
            bWriteMetalEnvironment = true;
          }
          else {
            bWriteMetalEnvironment = false;
          }
        }
      }

      else {
        /*!
           Unknown command
        */
        std::cout << " unknown command: " << inputFileContents[i][0] << std::endl;
        exit(0);
      }

    } // for loop

    // - Clean up - //
    delete pPdbParser;
    delete pElementParser;
    delete pCollection;
    return 0;
}

