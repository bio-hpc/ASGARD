/*!
   \file sequenceAligner.cpp

   \brief Aligns two sequences

   \author Martin B. Peters

   $Date: 2010/04/22 22:19:54 $
   $Revision: 1.2 $

   ----------------------------------------------------------------------------

   MTK++ - C++ package of modeling libraries.

   Copyright (C) 2005-2007  (see AUTHORS file for a list of contributors)

   ----------------------------------------------------------------------------
*/
// - Utils
#include "Utils/printHeader.h"
#include "Utils/vector3d.h"
#include "Utils/constants.h"

// - MOLECULE
#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/submolecule.h"
#include "Molecule/atom.h"
#include "Molecule/torsion.h"
#include "Molecule/element.h"
#include "Molecule/superimpose.h"
#include "Molecule/seqAlign.h"
#include "Molecule/parameters.h"
#include "Molecule/stdLibrary.h"
#include "Molecule/atomTyper.h"
#include "Molecule/connections.h"
#include "Molecule/complex.h"

// - STATS
#include "Statistics/sheet.h"
#include "Statistics/table.h"

// - PARSERS
#include "Parsers/elementParser.h"
#include "Parsers/molParser.h"
#include "Parsers/pdbParser.h"
#include "Parsers/pamParser.h"
#include "Parsers/dMParser.h"
#include "Parsers/StringManip.h"
#include "Parsers/paramParser.h"
#include "Parsers/stdLibParser.h"
#include "Parsers/inputParser.h"
#include "Parsers/commLineOptions.h"

// - TIME
#include "time.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>

using namespace MTKpp;

/*!
   \brief Aligns two sequences
   \param argc
   \param argv
   \return success
*/
int main (int argc, char **argv) 
{
    std::string prog_name = "sequenceAligner";
    std::vector<std::string> authors;
    std::string author = "Martin B. Peters";
    authors.push_back(author);

    // 1. CREATE AN OBJECT
    commLineOptions *clo = new commLineOptions();

    // 2. SET PREFERENCES
    clo->noUsage();

    // 3. SET THE USAGE/HELP
    clo->addUsage( "  sequenceAligner: Sequence Alignment and Structural   " );
    clo->addUsage( "                   Superimposition                   \n" );
    clo->addUsage( "    usage:  sequenceAligner [flags] [options]        \n" );
    clo->addUsage( "  options:                                             " );
    clo->addUsage( "          -i input file                                " );
    clo->addUsage( "          -o log file                                \n" );
    clo->addUsage( "    flags:                                             " );
    clo->addUsage( "          -h help                                    \n" );

    // 4. SET THE OPTION STRINGS/CHARACTERS
    clo->setOption( "inFile",    'i' );
    clo->setOption( "logFile",   'o' );
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

    std::string inFile  = "";
    std::string logFile = "";
    bool bInput = false;
    bool bLog = false;

    if ( clo->getValue( "i" ) != 0 ) {
      inFile = clo->getValue( "i" );
      bInput = true;
    }
    else if ( clo->getValue( "inFile" ) != 0 ) {
      inFile =  clo->getValue( "inFile" );
      bInput = true;
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide an input file " << std::endl;
      return 0;
    }

    if ( clo->getValue( "o" ) != 0 ) {
      logFile = clo->getValue( "o" );
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

    int failure = 0;
    std::string errorMessage = "";

    // Read input file
    std::vector<std::vector<std::string> > inputFileContents;
    if (bInput) {
      failure = readInputFile(inFile, inputFileContents);
      if (failure) {
        MTKpp::errorLogger.throwError("sequenceAligner", "Failed to open input file ", 1);
        exit(0);
      }
    }

    // Set up collection
    std::string AMBERHOME = getenv("AMBERHOME");
    collection* pCollection = new collection();
    pCollection->addParameters();
    pCollection->addStdLibrary();
    atomTyper* pAtomTyper = new atomTyper(0);

    std::map<std::string, std::string> variableMap;
    typedef std::vector<std::vector<std::string> >::iterator myIterator;

    // If the keyword 'source' is found in the input file then place the
    // commands found in the sourced file into the input file
    bool sourceFound = true;
    while (sourceFound) {
      sourceFound = false;
      int myIndex = 0;
      int sourceIndex = 0;
      for (myIterator i = inputFileContents.begin(); i != inputFileContents.end(); i++) {
        std::vector<std::string> cur = *i;

        if (cur[0] == "source") {
          /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: source

           Description: Source a global file

           syntax: source file_name
           \endcode
           */
          sourceFound = true;
          if ((cur.size() == 2)) {
            sourceIndex = myIndex;
            if (fileExists(cur[1])) {
              std::vector<std::vector<std::string> > sourceInputFileContents;
              failure = readInputFile(cur[1], sourceInputFileContents);
              if (failure) {
                errorMessage =  " Error Reading File " + cur[1];
                MTKpp::errorLogger.throwError("sequenceAligner::source", errorMessage, 1);
                exit(1);
              }
              inputFileContents.insert(i+1, sourceInputFileContents.begin(),
                                            sourceInputFileContents.end());
            }
            else {
              errorMessage =  " Error Reading File " + cur[1];
              MTKpp::errorLogger.throwError("sequenceAligner::source", errorMessage, 1);
              exit(1);
            }
          }
          else {
            errorMessage =  " Improper use of the source command ";
            MTKpp::errorLogger.throwError("sequenceAligner::source", errorMessage, 1);
            exit(1);
          }
        }
        if (sourceFound) {
          inputFileContents[sourceIndex][0] = "# source_read";
          break;
        }
        myIndex++;
      }
    }

    // Template molecule
    molecule* pMoleculeA = 0;

    // Query molecule
    molecule* pMoleculeB = 0;

    // Ligand molecule
    molecule* pLigand = 0;

    // Active site residue indices
    std::vector<int> fullList;

    // Active site molecule indices
    std::vector<int> fullListMols;

    // Create parsers
    pdbParser* pPdbParser = new pdbParser();

    // Create seqAlign object
    seqAlign* pSeqAlign = new seqAlign();

    // Create superimpose object
    superimpose* pSuperimpose = 0;

    // Create results container
    std::map<std::string, sheet*> resultsMap;
    std::string labels [3] = {"x", "y", "z"};

    // Storage the begining and end of each pdb file read in
    std::vector<std::vector<int> > startEnd;
    std::map<std::string, int> mapStartEnd;
    mapStartEnd["start"] = 0;
    mapStartEnd["end"] = 1;
    typedef std::map<std::string, int>::iterator startEndMapIterator;

    // complex pointer
    complex* pComplex = 0;

    std::vector<std::string> ca;
    ca.push_back(" CA ");

    std::vector<std::string> bb;
    bb.push_back(" CA ");
    bb.push_back(" N  ");
    bb.push_back(" C  ");
    bb.push_back(" O  ");

    std::vector<std::string> bbb;
    bbb.push_back(" CA ");
    bbb.push_back(" N  ");
    bbb.push_back(" C  ");
    bbb.push_back(" O  ");
    bbb.push_back(" CB ");

    std::map<std::string, std::vector<std::string> > atomGroups;
    atomGroups["alphaCarbons"] = ca;
    atomGroups["bb"] = bb;
    atomGroups["bbb"] = bbb;

    // Loop over the contents of the input file contents
    for (unsigned int i = 0; i < inputFileContents.size(); i++) {
      if (inputFileContents[i][0] == "# source_read") continue;

      if (inputFileContents[i][0] == "quit") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
            Function: quit

            Description: Exits program

            syntax: quit

           \endcode
        */
        time (&endTime);
        int diffTime = (int) difftime(endTime, startTime);
        std::string errMessage = " Exited Normally After "
                    + int2String(diffTime) + " Seconds ";
        MTKpp::errorLogger.throwError("sequenceAligner::quit", errMessage, MESSAGE);
        oLog.close();
        exit(0);
      }

      else if (inputFileContents[i][0] == "setLoggingLevel") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: setLoggingLevel

           Description: Set the verbosity of error/warning/info messages

           syntax: setLoggingLevel 1

           Values:
             1 - Error
             2 - Warning
             3 - Debug
             4 - Info
           \endcode
        */
        if ((inputFileContents[i].size() == 2)) {
          int l = atoi(inputFileContents[i][1].c_str());
          MTKpp::errorLogger.setLevel(l);
        }
      }

      else if (inputFileContents[i][0] == "set") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: set

           Description: Set variable

           syntax: set variable_name value
           \endcode
        */
        if ((inputFileContents[i].size() == 3)) {
          variableMap[inputFileContents[i][1]] = inputFileContents[i][2];
          for (unsigned int j = i+1; j < inputFileContents.size(); j++) {
            bool bReplaced = false;
            for (unsigned int k = 0; k < inputFileContents[j].size(); k++) {
              if (containsSubStr(inputFileContents[j][k], inputFileContents[i][1])) {
                inputFileContents[j][k] = replaceSubStr(inputFileContents[j][k],
                         inputFileContents[i][1], inputFileContents[i][2]);
                bReplaced = true;
              }
            }
          }
        }
      }

      else if (inputFileContents[i][0] == "loadElements") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: loadElements

           Description: Loads AMBER Parameters into sequenceAligner

           syntax: loadElements ~/MTKpp/data/elements.xml
           \endcode
        */
        if ((inputFileContents[i].size() != 2) or (!pCollection)) {
          MTKpp::errorLogger.throwError("sequenceAligner::loadElements", " Incorrect use ", MTK_ERROR);
          exit(1);
        }
        else {
          if (!fileExists(inputFileContents[i][1])) {
            errorMessage =  " Can't find " + inputFileContents[i][1] + " ... exiting ";
            MTKpp::errorLogger.throwError("sequenceAligner::loadElements", errorMessage, MTK_ERROR);
            exit(1);
          }

          elementParser* pElementParser = new elementParser(pCollection->pElements);
          if (pElementParser) {
            pElementParser->Read(inputFileContents[i][1]);
            if (pElementParser->getError()) {
              errorMessage =  pElementParser->getErrorMessage() + " ... exiting ";
              MTKpp::errorLogger.throwError("sequenceAligner::loadElements", errorMessage, MTK_ERROR);
              exit(1);
            }
            errorMessage =  " Read elements file: " + inputFileContents[i][1];
            MTKpp::errorLogger.throwError("sequenceAligner::loadElements", errorMessage, INFO);
            delete pElementParser;
          }
          else {
            errorMessage = " Error ... exiting ";
            MTKpp::errorLogger.throwError("sequenceAligner::loadElements", errorMessage, MTK_ERROR);
          }
        }
      }

      else if (inputFileContents[i][0] == "loadParam") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: loadParam

           Description: Loads AMBER Parameters into sequenceAligner

           syntax: loadParam ~/MTKpp/data/parm94.xml
           \endcode
        */
        if ((inputFileContents[i].size() != 2) or (!pCollection)) {
          MTKpp::errorLogger.throwError("sequenceAligner::loadParam", " Incorrect use ", MTK_ERROR);
          exit(1);
        }
        else {
          if (!fileExists(inputFileContents[i][1])) {
            errorMessage =  " Can't find " + inputFileContents[i][1] + " ... exiting ";
            MTKpp::errorLogger.throwError("sequenceAligner::loadParam", errorMessage, MTK_ERROR);
            exit(1);
          }
          paramParser* pParamParser = new paramParser(pCollection->getParameters());

          if (pParamParser) {
            if (!pParamParser->getError()) {
              pParamParser->Read(inputFileContents[i][1]);
              if (pParamParser->getError()) {
                errorMessage =  pParamParser->getErrorMessage() + " ... exiting ";
                MTKpp::errorLogger.throwError("sequenceAligner::loadParam", errorMessage, MTK_ERROR);
                exit(1);
              }
              errorMessage =  " Read parameters file: " + inputFileContents[i][1];
              MTKpp::errorLogger.throwError("sequenceAligner::loadParam", errorMessage, INFO);
              delete pParamParser;
            }
            else {
              errorMessage = " Error ... exiting ";
              MTKpp::errorLogger.throwError("sequenceAligner::loadParam", errorMessage, MTK_ERROR);
              exit(1);
            }
          }
          else {
            errorMessage = " Error ... exiting ";
            MTKpp::errorLogger.throwError("sequenceAligner::loadParam", errorMessage, MTK_ERROR);
            exit(1);
          }
        }
      }

      else if (inputFileContents[i][0] == "loadLib") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: loadLib

           Description: Loads AMBER library files into sequenceAligner

           syntax: loadLib ~/MTKpp/data/amino94.xml
           \endcode
        */
        if ((inputFileContents[i].size() != 2) or (!pCollection)) {
          MTKpp::errorLogger.throwError("sequenceAligner::loadLib", " Incorrect use. ", MTK_ERROR);
          exit(1);
        }
        else {
          if (!fileExists(inputFileContents[i][1])) {
            errorMessage =  " Can't find " + inputFileContents[i][1] + " ... exiting ";
            MTKpp::errorLogger.throwError("sequenceAligner::loadLib", errorMessage, MTK_ERROR);
            exit(1);
          }
          stdLibrary* pStdLibrary = pCollection->getStdLibrary();
          if (pStdLibrary) {
            parameters* pParms = pCollection->getParameters();
            if (!pParms) {
              errorMessage =  " Please read in parameters first ... exiting ";
              MTKpp::errorLogger.throwError("sequenceAligner::loadLib", errorMessage, MTK_ERROR);
              exit(1);
            }
            stdLibParser* pStdLibParser = new stdLibParser(pStdLibrary, pParms);
            if (pStdLibParser) {
              if (!pStdLibParser->getError()) {
                pStdLibParser->Read(inputFileContents[i][1]);
                if (pStdLibParser->getError()) {
                  errorMessage =  pStdLibParser->getErrorMessage() + " ... exiting ";
                  MTKpp::errorLogger.throwError("sequenceAligner::loadLib", errorMessage, MTK_ERROR);
                  exit(1);
                }
                errorMessage =  " Read library file: " + inputFileContents[i][1];
                MTKpp::errorLogger.throwError("sequenceAligner::loadLib", errorMessage, INFO);
                delete pStdLibParser;
              }
              else {
                MTKpp::errorLogger.throwError("sequenceAligner::loadLib", " Incorrect use. ", MTK_ERROR);
                delete pStdLibParser;
                exit(1);
              }
            }
          }
          else {
            MTKpp::errorLogger.throwError("sequenceAligner::loadLib", " Incorrect use. ", MTK_ERROR);
          }
        }
      }

      else if (inputFileContents[i][0] == "loadPam") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: loadPam

           Description: Loads PAM file

           syntax: loadPam ~/MTKpp/data/PAM/PAM250
           \endcode
        */
        if ((inputFileContents[i].size() != 2) or (!pSeqAlign)) {
          MTKpp::errorLogger.throwError("sequenceAligner::loadPam", " Incorrect use ", MTK_ERROR);
          exit(1);
        }
        else {
          if (!fileExists(inputFileContents[i][1])) {
            errorMessage =  " Can't find " + inputFileContents[i][1] + " ... exiting ";
            MTKpp::errorLogger.throwError("sequenceAligner::loadPam", errorMessage, MTK_ERROR);
            exit(1);
          }
          pamParser* pPamParser = new pamParser();
          if (pPamParser) {
            if (!pPamParser->getError()) {
              pPamParser->Read(inputFileContents[i][1], pSeqAlign);
              if (pPamParser->getError()) {
                errorMessage =  pPamParser->getErrorMessage() + " ... exiting ";
                MTKpp::errorLogger.throwError("sequenceAligner::loadPam", errorMessage, MTK_ERROR);
                exit(1);
              }
              errorMessage =  " Read PAM file: " + inputFileContents[i][1];
              MTKpp::errorLogger.throwError("sequenceAligner::loadPam", errorMessage, INFO);
              delete pPamParser;
            }
            else {
              errorMessage = " Error ... exiting ";
              MTKpp::errorLogger.throwError("sequenceAligner::loadParam", errorMessage, MTK_ERROR);
              delete pPamParser;
              exit(1);
            }
          }
          else {
            errorMessage = " Error ... exiting ";
            MTKpp::errorLogger.throwError("sequenceAligner::loadPam", errorMessage, MTK_ERROR);
            exit(1);
          }
        }
      }

      else if (inputFileContents[i][0] == "readPdb") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: readPdb

           Description: Opens a PDB file

           syntax: readPdb 1FEE.pdb
           \endcode
        */
        if (inputFileContents[i].size() != 2) {
          MTKpp::errorLogger.throwError("sequenceAligner::readPdb", " Incorrect use of readPdb ", MTK_ERROR);
          exit(1);
        }
        else {
          std::vector<int> stEn;
          int nStart = pCollection->getNumberMolecules()+1;
          pPdbParser->Read(inputFileContents[i][1], pCollection);
          int nEnd = pCollection->getNumberMolecules();

          stEn.push_back(nStart);
          stEn.push_back(nEnd);
          startEnd.push_back(stEn);
        }
      }

      else if (inputFileContents[i][0] == "removeHs") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: removeHs

           Description: 

           syntax: removeHs
           \endcode
        */
        if (inputFileContents[i].size() != 1) {
          MTKpp::errorLogger.throwError("sequenceAligner::removeHs", " Incorrect use ... exiting ", MTK_ERROR);
          exit(1);
        }
        else {
          std::vector<molecule*> molList = pCollection->getMoleculeList();
          for (unsigned int j = 0; j < molList.size(); j++) {
            molList[j]->removeHydrogens();
          }
        }
      }

      else if (inputFileContents[i][0] == "print") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: print

           Description: Print to screen details of structure

           syntax: print
           syntax: print protein
           syntax: print ligand
           syntax: print solvent
           \endcode
        */
        std::vector<molecule*> molList = pCollection->getMoleculeList();
        if (inputFileContents[i].size() == 1) {
          for (unsigned int j = 0; j < molList.size(); j++) {
            std::cout << " MOL: " << j+1 << " " << molList[j]->getName() << " " << std::endl;
            std::vector<submolecule*> smolList = molList[j]->getSubMoleculeList();
            for (unsigned int k = 0; k < smolList.size(); k++) {
              std::cout << smolList[k]->get1LName();
            }
            std::cout << " " << std::endl;
          }
        }
        else if (inputFileContents[i].size() == 2) {
          for (unsigned int j = 0; j < molList.size(); j++) {
            std::vector<submolecule*> smolList = molList[j]->getSubMoleculeList();
            if (molList[j]->getKind() == 1 and inputFileContents[i][1] == "protein") { // protein
              std::cout << " MOL: " << j+1 << " " << molList[j]->getName() << " " << std::endl;
              for (unsigned int k = 0; k < smolList.size(); k++) {
                std::cout << smolList[k]->get1LName();
              }
              std::cout << " " << std::endl;
            }
            else if (molList[j]->getKind() == 2 and inputFileContents[i][1] == "ligand") { // ligand
              std::cout << " MOL: " << j+1 << " " << molList[j]->getName() << " " << std::endl;
              for (unsigned int k = 0; k < smolList.size(); k++) {
                std::cout << smolList[k]->get1LName();
              }
              std::cout << " " << std::endl;
            }
            else if (molList[j]->getKind() == 3 and inputFileContents[i][1] == "solvent") { // solvent
              std::cout << " MOL: " << j+1 << " " << molList[j]->getName() << " " << std::endl;
              for (unsigned int k = 0; k < smolList.size(); k++) {
                std::cout << smolList[k]->get1LName();
              }
              std::cout << " " << std::endl;
            }
          }
        }
      }

      else if (inputFileContents[i][0] == "assignDisulfideBonds") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function assignDisulfideBonds

           Description: Assigns all disulfide bonds (Needs to be carried out before atomtyping)

           syntax assignDisulfideBonds
           \endcode
        */
        if (pCollection) {
          connections* pConnections = new connections(pCollection);
          pConnections->assignDisulfideBonds();
          errorMessage =  " Assigned ";
          MTKpp::errorLogger.throwError("sequenceAligner::assignDisulfideBonds", errorMessage, INFO);
          delete pConnections;
        }
      }

      else if (inputFileContents[i][0] == "atomType") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: atomType

           Description: Assigns atom types in the collection

           syntax: atomType
           \endcode
        */
        if (inputFileContents[i].size() == 1) {
          pAtomTyper->atomTypeByLib(pCollection);
        }
      }

      else if (inputFileContents[i][0] == "assignConnectivity") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: assignConnectivity

           Description: Assigns all bonds, angles, torsions and impropers

           syntax: assignConnectivity
           \endcode
        */
        if (pCollection) {
          connections* pConnections = new connections(pCollection);
          pConnections->run();
          MTKpp::errorLogger.throwError("sequenceAligner::assignConnectivity",
                 " successful ", INFO);
          delete pConnections;
        }
      }

      else if (inputFileContents[i][0] == "assignParameters") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: assignParameters

           Description: Assigns bond/angle/torsion/improper parameters

           syntax: assignParameters /COL/MOL
           \endcode
        */
        connections* pConnections = new connections(pCollection);

        if (inputFileContents[i].size() == 1) {
          std::vector<molecule*> m = pCollection->getMoleculeList();
          for (unsigned int o = 0; o < m.size(); o++) {
            pConnections->assignStd(m[o]);
          }
        }
        delete pConnections;
      }

      else if (inputFileContents[i][0] == "systemSetup") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: systemSetup

           Description: Select Receptor, Solvent, and Ligand

           The receptor is defined as any molecule with a molecular weight of >1000 A.U.
           The receptor can also contain metal atoms.

           The function relies on the solvent molecule being named either HOH or WAT

           The ligand is what remains.

           syntax: systemSetup
           \endcode
        */
        if (inputFileContents[i].size() == 1) {
          pComplex = new complex(pCollection, 1);

          if (failure) {
            MTKpp::errorLogger.throwError("systemSetup", " failure ", MTK_ERROR);
            exit(1);
          }
          else {
            MTKpp::errorLogger.throwError("systemSetup", " successful ", INFO);
          }
        }
        else {
          MTKpp::errorLogger.throwError("systemSetup", " failure ", MTK_ERROR);
          exit(1);
        }
      }

      else if (inputFileContents[i][0] == "setTemplate") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: setTemplate

           Description: set the template molecule

           syntax: setTemplate 1
           syntax: setTemplate 1@start
           \endcode
        */
        int molIndex = 0;
        if (inputFileContents[i].size() == 2) {
          if (containsSubStr(inputFileContents[i][1], "@")) {
            std::vector<std::string> w;
            splitString(inputFileContents[i][1], "@", w, 0);
            if (w.size() == 2) {
              int m = string2Int(w[0]) - 1;
              int m2 = 0;
              startEndMapIterator p = mapStartEnd.find(w[1]);
              if (p != mapStartEnd.end()) {
                m2 = p->second;
                molIndex = startEnd[m][m2];
              }
              else {
                MTKpp::errorLogger.throwError("sequenceAligner::setTemplate",
                       " Incorrect use ", MTK_ERROR);
                exit(1);
              }
            }
            else {
              MTKpp::errorLogger.throwError("sequenceAligner::setTemplate",
                       " Incorrect use ", MTK_ERROR);
              exit(1);
            }
          }
          else {
            molIndex = string2Int(inputFileContents[i][1]);
          }
          std::vector<molecule*> m = pCollection->getMoleculeList();
          int nMols = m.size();
          for (int o = 0; o < nMols; o++) {
            if (molIndex == o+1) {
              pMoleculeA = m[o];
              failure = pSeqAlign->setTemplate(pMoleculeA);
              if (failure) {
                MTKpp::errorLogger.throwError("sequenceAligner::setTemplate",
                       " Incorrect use ", MTK_ERROR);
                exit(1);
              }
            }
          }
        }
      }

      else if (inputFileContents[i][0] == "setQuery") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: setQuery

           Description: set the query molecule

           syntax: setQuery 1
           syntax: setQuery 2@end
           \endcode
        */
        int molIndex = 0;
        if (inputFileContents[i].size() == 2) {
          if (containsSubStr(inputFileContents[i][1], "@")) {
            std::vector<std::string> w;
            splitString(inputFileContents[i][1], "@", w, 0);
            if (w.size() == 2) {
              int m = string2Int(w[0]) - 1;
              int m2 = 0;
              startEndMapIterator p = mapStartEnd.find(w[1]);
              if (p != mapStartEnd.end()) {
                m2 = p->second;
                molIndex = startEnd[m][m2];
              }
              else {
                MTKpp::errorLogger.throwError("sequenceAligner::setQuery",
                       " Incorrect use ", MTK_ERROR);
                exit(1);
              }
            }
            else {
              MTKpp::errorLogger.throwError("sequenceAligner::setQuery",
                       " Incorrect use ", MTK_ERROR);
              exit(1);
            }
          }
          else {
            molIndex = string2Int(inputFileContents[i][1]);
          }

          std::vector<molecule*> m = pCollection->getMoleculeList();
          int nMols = m.size();
          for (int o = 0; o < nMols; o++) {
            if (molIndex == o+1) {
              pMoleculeB = m[o];
              failure = pSeqAlign->setQuery(pMoleculeB);
              if (failure) {
                MTKpp::errorLogger.throwError("sequenceAligner::setQuery",
                       " Incorrect use ", 1);
                exit(1);
              }
            }
          }
        }
      }

      else if (inputFileContents[i][0] == "setLigand") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: setLigand

           Description: 

           syntax: setLigand 3@start
           syntax: setLigand 3
           \endcode
        */
        if (inputFileContents[i].size() == 2) {
          int ligIndex = 0;
          if (containsSubStr(inputFileContents[i][1], "@")) {
            std::vector<std::string> w;
            splitString(inputFileContents[i][1], "@", w, 0);
            if (w.size() == 2) {
              int m = string2Int(w[0]) - 1;
              int m2 = 0;
              startEndMapIterator p = mapStartEnd.find(w[1]);
              if (p != mapStartEnd.end()) {
                m2 = p->second;
                ligIndex = startEnd[m][m2];
              }
              else {
                MTKpp::errorLogger.throwError("sequenceAligner::setLigand",
                       " Incorrect use ", MTK_ERROR);
                exit(1);
              }
            }
            else {
              MTKpp::errorLogger.throwError("sequenceAligner::setLigand",
                       " Incorrect use ", MTK_ERROR);
              exit(1);
            }
          }
          else {
            ligIndex = string2Int(inputFileContents[i][1]);
          }

          std::vector<molecule*> m = pCollection->getMoleculeList();
          int nMols = m.size();
          for (int o = 0; o < nMols; o++) {
            if (ligIndex == o+1) {
              pLigand = m[o];
            }
          }
        }
        else {
          MTKpp::errorLogger.throwError("sequenceAligner::setLigand", " Incorrect use ", MTK_ERROR);
        }
      }

      else if (inputFileContents[i][0] == "runAlign") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: runAlign

           Description: Run sequence alignment

           syntax: runAlign 1 1 10.0 0.5
           \endcode
        */
        if (inputFileContents[i].size() == 5) {
          int algorithmType = string2Int(inputFileContents[i][1]);
          int gapPenaltyType = string2Int(inputFileContents[i][2]);
          double gapOpen = string2Double(inputFileContents[i][3]);
          double gapExtend = string2Double(inputFileContents[i][4]);

          pSeqAlign->setAlgorithmType(algorithmType);
          pSeqAlign->setGapPenaltyType(gapPenaltyType);
          pSeqAlign->setGapOpen(gapOpen);
          pSeqAlign->setGapExtend(gapExtend);
          failure = pSeqAlign->run();
          if (failure) {
            MTKpp::errorLogger.throwError("sequenceAligner::runAlign", " Failed ", MTK_ERROR);
          }
        }
        else {
          MTKpp::errorLogger.throwError("sequenceAligner::runAlign", " Incorrect use ", MTK_ERROR);
        }
      }

      else if (inputFileContents[i][0] == "superimpose") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: superimpose

           Description: Superimpose

           syntax: superimpose alphaCarbons
           syntax: superimpose bb
           syntax: superimpose bbb
           \endcode
        */
        if (inputFileContents[i].size() == 2) {
          int* corrMap = pSeqAlign->getCorrMap();
          int qNRes = pMoleculeB->getNumSubMolecules();

          double centerA[3];
          double centerB[3];
          double rotMat[3][3];
          int index = 0;

          // Determine the number of matches
          int nMatches = 0;
          for (int t = 0; t < qNRes; t++) {
            if (corrMap[t] > 0) {
              nMatches++;
            }
          }

          if (inputFileContents[i][1] == "alphaCarbons") {
            // Align Alpha Carbons
            double coordsA[nMatches][3];
            double coordsB[nMatches][3];

            for (int t = 0; t < qNRes; t++) {
              if (corrMap[t] > 0) {
                submolecule* pSubMolA = pMoleculeA->getSubMolecule(corrMap[t], 1, 0);
                submolecule* pSubMolB = pMoleculeB->getSubMolecule(t+1, 1, 0);

                atom* pAtomA = pSubMolA->getAtom(" CA ");
                atom* pAtomB = pSubMolB->getAtom(" CA ");

                if (!pAtomA or !pAtomB) {
                  std::cout << " Error in sequenceAligner, residue doesn't have a CA atom " << std::endl;
                  exit(1);
                }

                vector3d* pCoordsA = pAtomA->getCoords();
                vector3d* pCoordsB = pAtomB->getCoords();

                coordsA[index][0] = pCoordsA->getX();
                coordsA[index][1] = pCoordsA->getY();
                coordsA[index][2] = pCoordsA->getZ();

                coordsB[index][0] = pCoordsB->getX();
                coordsB[index][1] = pCoordsB->getY();
                coordsB[index][2] = pCoordsB->getZ();
                index++;
              }
            }

            pSuperimpose = new superimpose(nMatches);
            // Calculate the centers
            pSuperimpose->center(coordsA, centerA);
            pSuperimpose->center(coordsB, centerB);

            failure = pSuperimpose->getRotationMatrix(coordsA, coordsB,
                                    centerA, centerB, nMatches, rotMat);

            // Get RMSD Value
            double rmsd = pSuperimpose->fit(coordsA, coordsB, nMatches);

            errorMessage = " RMSD OF ALPHA-CARBONS = " + double2String(rmsd);
            MTKpp::errorLogger.throwError("sequenceAligner::runAlign", errorMessage, 0);
          }

          // Save transformation
          sheet* pSheet = new sheet();
          resultsMap[inputFileContents[i][1]] = pSheet;
          table<double>* pDoubleTable_c1;
          pDoubleTable_c1 = pSheet->addTable();
          pDoubleTable_c1->setName("Template Center");
          pDoubleTable_c1->setSizes(1, 3);

          table<double>* pDoubleTable_c2;
          pDoubleTable_c2 = pSheet->addTable();
          pDoubleTable_c2->setName("Query Center");
          pDoubleTable_c2->setSizes(1, 3);

          for (unsigned int a = 0; a < 3; a++) {
            pDoubleTable_c1->setColumnLabel(a, labels[a]);
            pDoubleTable_c2->setColumnLabel(a, labels[a]);
            pDoubleTable_c1->setCellValue(0, a, centerA[a]);
            pDoubleTable_c2->setCellValue(0, a, centerB[a]);
          }

          table<double>* pDoubleTable_rotMat;
          pDoubleTable_rotMat = pSheet->addTable();
          pDoubleTable_rotMat->setName("Rotation Matrix");
          pDoubleTable_rotMat->setSizes(3, 3);

          for (int b = 0; b < 3; b++) {
            for (unsigned int a = 0; a < 3; a++) {
              pDoubleTable_rotMat->setCellValue(b, a, rotMat[b][a]);
            }
          }
        }
        else {
          MTKpp::errorLogger.throwError("sequenceAligner::superimpose", " Incorrect use ", MTK_ERROR);
        }
      }

      else if (inputFileContents[i][0] == "transform") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: transform

           Description: transform

           syntax: transform 1 140 alphaCarbons
           syntax: transform 2@start 2@end alphaCarbons
           \endcode
        */
        if (inputFileContents[i].size() == 4) {

          int start = 1;
          int end = 1;
          if (containsSubStr(inputFileContents[i][1], "@")) {
            std::vector<std::string> w;
            splitString(inputFileContents[i][1], "@", w, 0);
            if (w.size() == 2) {
              int m = string2Int(w[0]) - 1;
              int m2 = 0;
              startEndMapIterator p = mapStartEnd.find(w[1]);
              if (p != mapStartEnd.end()) {
                m2 = p->second;
                start = startEnd[m][m2];
              }
              else {
                MTKpp::errorLogger.throwError("sequenceAligner::transform",
                       " Incorrect use ", MTK_ERROR);
                exit(1);
              }
            }
            else {
              MTKpp::errorLogger.throwError("sequenceAligner::transform",
                       " Incorrect use ", MTK_ERROR);
              exit(1);
            }
          }
          else {
            start = string2Int(inputFileContents[i][1]);
          }

          if (containsSubStr(inputFileContents[i][2], "@")) {
            std::vector<std::string> w;
            splitString(inputFileContents[i][2], "@", w, 0);
            if (w.size() == 2) {
              int m = string2Int(w[0]) - 1;
              int m2 = 0;
              startEndMapIterator p = mapStartEnd.find(w[1]);
              if (p != mapStartEnd.end()) {
                m2 = p->second;
                end = startEnd[m][m2];
              }
              else {
                MTKpp::errorLogger.throwError("sequenceAligner::transform",
                       " Incorrect use ", MTK_ERROR);
                exit(1);
              }
            }
            else {
              MTKpp::errorLogger.throwError("sequenceAligner::transform",
                       " Incorrect use ", MTK_ERROR);
              exit(1);
            }
          }
          else {
            end = string2Int(inputFileContents[i][2]);
          }
 
          std::vector<molecule*> m = pCollection->getMoleculeList();
          molecule* pMol = 0;
          int nAtoms = 0;
          int lIndex = 0;
          int nMols = m.size();
          for (int o = 0; o < nMols; o++) {
            if (o+1 >= start and o+1 <= end) {
              pMol = m[o];
              nAtoms += pMol->getNumAtoms();
            }
          }

          // Get atom coordinates
          double coords[nAtoms][3];
          for (int o = 0; o < nMols; o++) {
            if (o+1 >= start and o+1 <= end) {
              pMol = m[o];
              std::vector<atom*> atomList = pMol->getAtomList();
              for (unsigned int a = 0; a < atomList.size(); a++) {
                coords[lIndex][0] = atomList[a]->getX();
                coords[lIndex][1] = atomList[a]->getY();
                coords[lIndex][2] = atomList[a]->getZ();
                lIndex++;
              }
            }
          }

          // Transformation vectors and matrices
          double centerA[3];
          double centerB[3];
          double rotMat[3][3];

          // Get transformation data
          sheet* pSheet = resultsMap[inputFileContents[i][3]];

          table<double>* pDoubleTable_c1;
          pDoubleTable_c1 = pSheet->getTable("Template Center");

          table<double>* pDoubleTable_c2;
          pDoubleTable_c2 = pSheet->getTable("Query Center");

          for (unsigned int y = 0; y < 3; y++) {
            centerA[y] = pDoubleTable_c1->getCellValue(0,y);
            centerB[y] = pDoubleTable_c2->getCellValue(0,y);
          }

          table<double>* pDoubleTable_rotMat;
          pDoubleTable_rotMat = pSheet->getTable("Rotation Matrix");

          for (int b = 0; b < 3; b++) {
            for (unsigned int a = 0; a < 3; a++) {
              rotMat[b][a] = pDoubleTable_rotMat->getCellValue(b, a);
            }
          }

          // Transform the coordinates into the template's frame of reference
          pSuperimpose->updateCoords(coords, nAtoms, centerB, centerA, rotMat);

          // Update coordinates
          lIndex = 0;
          for (int o = 0; o < nMols; o++) {
            if (o+1 >= start and o+1 <= end) {
              pMol = m[o];
              std::vector<atom*> atomList = pMol->getAtomList();
              for (unsigned int a = 0; a < atomList.size(); a++) {
                atom* pAtom = atomList[a];
                vector3d* pCoords = pAtom->getCoords();
                pCoords->set(coords[lIndex][0], coords[lIndex][1], coords[lIndex][2]);
                lIndex++;
              }
            }
          }
        }
        else {
          MTKpp::errorLogger.throwError("sequenceAligner::transform", " Incorrect use ", MTK_ERROR);
        }
      }

      else if (inputFileContents[i][0] == "writePdb") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: writePdb

           Description: writePdb

           syntax: writePdb 140 340 file.pdb
           \endcode
        */
        if (inputFileContents[i].size() == 4) {

          int start = 1;
          int end = 1;
          if (containsSubStr(inputFileContents[i][1], "@")) {
            std::vector<std::string> w;
            splitString(inputFileContents[i][1], "@", w, 0);
            if (w.size() == 2) {
              int m = string2Int(w[0]) - 1;
              int m2 = 0;
              startEndMapIterator p = mapStartEnd.find(w[1]);
              if (p != mapStartEnd.end()) {
                m2 = p->second;
                start = startEnd[m][m2];
              }
              else {
                MTKpp::errorLogger.throwError("sequenceAligner::transform",
                       " Incorrect use ", MTK_ERROR);
                exit(1);
              }
            }
            else {
              MTKpp::errorLogger.throwError("sequenceAligner::transform",
                       " Incorrect use ", MTK_ERROR);
              exit(1);
            }
          }
          else {
            start = string2Int(inputFileContents[i][1]);
          }

          if (containsSubStr(inputFileContents[i][2], "@")) {
            std::vector<std::string> w;
            splitString(inputFileContents[i][2], "@", w, 0);
            if (w.size() == 2) {
              int m = string2Int(w[0]) - 1;
              int m2 = 0;
              startEndMapIterator p = mapStartEnd.find(w[1]);
              if (p != mapStartEnd.end()) {
                m2 = p->second;
                end = startEnd[m][m2];
              }
              else {
                MTKpp::errorLogger.throwError("sequenceAligner::transform",
                       " Incorrect use ", MTK_ERROR);
                exit(1);
              }
            }
            else {
              MTKpp::errorLogger.throwError("sequenceAligner::transform",
                       " Incorrect use ", MTK_ERROR);
              exit(1);
            }
          }
          else {
            end = string2Int(inputFileContents[i][2]);
          }

          // Create a vector of molecules
          std::vector<molecule*> molsToWrite;
          std::vector<molecule*> molList = pCollection->getMoleculeList();
          int nMols = molList.size();
          for (int j = 0; j < nMols; j++) {
            if ((j+1 >= start) and (j+1 <= end)) {
              molsToWrite.push_back(molList[j]);
            }
          }
          pPdbParser->Write(inputFileContents[i][3], molsToWrite);
        }
        else {
          MTKpp::errorLogger.throwError("sequenceAligner::writePdb", " Incorrect use ", MTK_ERROR);
        }
      }

      else if (inputFileContents[i][0] == "writeTransformationFile") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: writeTransformationFile

           Description: writeTransformationFile

           syntax: writeTransformationFile alphaCarbons t.xml
           \endcode
        */
        if (inputFileContents[i].size() == 3) {
          dMParser* pDmParser = new dMParser();
          pDmParser->write(resultsMap[inputFileContents[i][1]], inputFileContents[i][2]);
          delete pDmParser;
        }
        else {
          MTKpp::errorLogger.throwError("sequenceAligner::writeTransformationFile", " Incorrect use ", MTK_ERROR);
        }
      }

      else if (inputFileContents[i][0] == "createActiveSite") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: createActiveSite

           Description: Creates an active site in the template around the provided ligand

           syntax: createActiveSite
           \endcode
        */
        if (inputFileContents[i].size() == 1 and pMoleculeA and pLigand) {
          double cutoff = 10.0;
          int gapFill = 3;
          //MTKpp::errorLogger.throwError("sequenceAligner::createActiveSite", " Get all residues not in the ligand ", INFO);

          std::vector<molecule*> molList = pCollection->getMoleculeList();
          std::vector<submolecule*> recSmolList = pMoleculeA->getSubMoleculeList();

          std::vector<atom*> ligAtoms = pLigand->getAtomList();

          atom* pRecAtom = 0;
          atom* pLigAtom = 0;
          vector3d Coord1;
          vector3d Coord2;
          bool keep = false;

          std::vector<submolecule*> activeSiteResiduesSMol;
          std::vector<int> activeSiteResidues;
          std::vector<int> activeSiteMolecules;

          std::vector<int> aSiteResidues;
          std::vector<int> aSiteMolecules;
 
          bool noH = true;

          for (unsigned int r = 0; r < recSmolList.size(); r++) {
            submolecule* pSubMol = recSmolList[r];
            std::vector<atom*> atomList = pSubMol->getAtomList();
            for (unsigned int j = 0; j < atomList.size(); j++) {
              pRecAtom = atomList[j];
              Coord1 = (*pRecAtom->getCoords());
              for (unsigned int k = 0; k < ligAtoms.size(); k++) {
                pLigAtom = ligAtoms[k];

                if (noH and pRecAtom->getElementSymbol() == "H" and pLigAtom->getElementSymbol() == "H") {
                  continue;
                }

                Coord2 = (*pLigAtom->getCoords());

                if (Coord1.dist(Coord2) < cutoff) {
                  keep = true;
                  break;
                }
              }
              if (keep) {
                break;
              }
            }
            if (keep) {
              activeSiteResiduesSMol.push_back(pSubMol);
              aSiteResidues.push_back(pSubMol->getSubMolId());
              aSiteMolecules.push_back(pSubMol->getParent()->getMolId());
              keep = false;
            }
          }

          // -- Add metal center --

          if (aSiteResidues.size() < 1) {
            MTKpp::errorLogger.throwError("sequenceAligner::createActiveSite", " No active site found", MTK_ERROR);
            return 1;
          }

          std::string resInd = " Residue Indices in the active site:\n";
          std::sort(activeSiteResiduesSMol.begin(), activeSiteResiduesSMol.end(), submolecule::less);
          for (unsigned int a = 0; a < activeSiteResiduesSMol.size(); a++) {
            activeSiteResidues.push_back(activeSiteResiduesSMol[a]->getSubMolId());
            activeSiteMolecules.push_back(activeSiteResiduesSMol[a]->getParent()->getMolId());
            resInd += int2String(activeSiteResiduesSMol[a]->getSubMolId()) + " ";
          }
          MTKpp::errorLogger.throwError("sequenceAligner::createActiveSite", resInd, INFO);

          // Create segments and fill in gaps
          int prevRes = activeSiteResidues[0];
          int prevMol = activeSiteMolecules[0];

          submolecule* pSubMol = 0;
          for (unsigned int a = 0; a < activeSiteResidues.size(); a++) {
            if (activeSiteMolecules[a] == prevMol) {
              molecule* pRecMol = pCollection->getMolecule(activeSiteMolecules[a]);
              int diff = activeSiteResidues[a] - prevRes;
              if ((diff <= gapFill) and (diff > 0)) {
                for (int j = prevRes+1; j < activeSiteResidues[a]; j++) {
                  pSubMol = pRecMol->getSubMolecule(j);
                  if (pSubMol) {
                    if ((pSubMol->getName() != "WAT") and (pSubMol->getName() != "HOH")) {
                      fullList.push_back(j);
                      fullListMols.push_back(activeSiteMolecules[a]);
                    }
                  }
                  else {
                    std::cout << " sequenceAligner::Warning, Can't find residue: " << j << std::endl;
                  }
                }
              }
            }
            fullList.push_back(activeSiteResidues[a]);
            fullListMols.push_back(activeSiteMolecules[a]);
            prevRes = activeSiteResidues[a];
            prevMol = activeSiteMolecules[a];
          }

          resInd = " Residue Indices in the active site after gap filling:\n";
          for (unsigned int f = 0; f < fullList.size(); f++) {
            resInd += int2String(fullList[f]) + " ";
          }
          MTKpp::errorLogger.throwError("sequenceAligner::createActiveSite", resInd, INFO);
        }
        else {
          MTKpp::errorLogger.throwError("sequenceAligner::createActiveSite", " Incorrect use ", MTK_ERROR);
        }
      }

      else if (inputFileContents[i][0] == "compareActiveSites") {
        /*!
           @ingroup sequenceAligner_Commands
           \code
           Function: compareActiveSites

           Description: Compares the active site in the template to the query

           syntax: compareActiveSites alphaCarbons
           syntax: compareActiveSites bb
           syntax: compareActiveSites bbb
           \endcode
        */
/*
file:///Applications/moe/html/apps/progeom.htm#MethodRotamer
            General       Proline        Glycine
      C-N   1.330(0.010)  1.335(0.011)   1.328(0.010)
      N-CA  1.459(0.012)  1.464(0.011)   1.453(0.012)
     CA-C   1.523(0.012)  1.523(0.013)   1.515(0.011)
      C-O   1.233(0.012)  1.234(0.012)   1.233(0.012)
     CB-CA  1.532(0.015)  1.532(0.011)

  pCA-pC-N  116.70(1.355) 117.90(1.299) 116.70(1.326)
   pC-N-CA  121.60(1.554) 120.50(2.203) 121.50(2.637)
   N-CA-C   110.80(2.506) 113.00(2.636) 113.30(1.303)
   CA-C-O   120.50(1.055) 120.40(1.311) 120.50(1.303)
   N-CA-CB  110.60(1.420) 103.40(1.000)

pCA-pC-PO-N 180.00(2.840) 180.70(3.210) [planarity]
pCA-pC-N-CA 179.30(6.600) 179.30(6.600) [omega trans]
pCA-pC-N-CA   0.10(4.700)   1.60(4.700) [omega cis]
pC-N-CA-C   [phi]
N-CA-C-N    [psi]
N-CA-CB-C   [C-beta]
*/
        if (inputFileContents[i].size() == 2 and pSeqAlign and pMoleculeB) {

          std::map<std::string, std::vector<std::string> >::iterator result;
          result = atomGroups.find(inputFileContents[i][1]);
          if (result == atomGroups.end()) {
            MTKpp::errorLogger.throwError("sequenceAligner::compareActiveSites", " Incorrect use ", MTK_ERROR);
            exit(1);
          }

          int nAtomsGroup = atomGroups[inputFileContents[i][1]].size();

          int* corrMap = pSeqAlign->getCorrMap();
          int qNRes = pMoleculeB->getNumSubMolecules();

          double centerA[3];
          double centerB[3];
          double rotMat[3][3];
          int index = 0;

          // Determine the number of matches
          int nMatches = 0;
          for (int t = 0; t < qNRes; t++) {
            if (corrMap[t] > 0) {
              std::vector<int>::iterator findInt;
              findInt = std::find(fullList.begin(), fullList.end(), corrMap[t]);

              if (findInt != fullList.end() ) {
                nMatches++;
              }
            }
          }

          int nIdentical = 0;
          int nSimilar = 0;
          int nDissimilar = 0;
          if (nMatches > 0) {
            double coordsA[nMatches * nAtomsGroup][3];
            double coordsB[nMatches * nAtomsGroup][3];
            for (int t = 0; t < qNRes; t++) {
              if (corrMap[t] > 0) {
                std::vector<int>::iterator findInt;
                findInt = std::find(fullList.begin(), fullList.end(), corrMap[t]);

                if (findInt != fullList.end() ) {
                  submolecule* pSubMolA = pMoleculeA->getSubMolecule(corrMap[t], 1, 0);
                  submolecule* pSubMolB = pMoleculeB->getSubMolecule(t+1, 1, 0);

                  std::string lA = pSubMolA->get1LName();
                  std::string lB = pSubMolA->get1LName();
                  if (lA == lB) {
                    nIdentical++;
                  }
                  else  {
                    double s = pSeqAlign->getSimScore(lA, lB);
                    if (s > -0.00000001) {
                      nSimilar++;
                    }
                    else {
                      nDissimilar++;
                    }
                  }

                  for (unsigned int y = 0; y < atomGroups[inputFileContents[i][1]].size(); y++) {
                    atom* pAtomA = pSubMolA->getAtom(atomGroups[inputFileContents[i][1]][y]);
                    atom* pAtomB = pSubMolB->getAtom(atomGroups[inputFileContents[i][1]][y]);

                    if (!pAtomA or !pAtomB) {
                      MTKpp::errorLogger.throwError("sequenceAligner::compareActiveSites", " Incorrect use ", MTK_ERROR);
                      exit(1);
                    }
                    vector3d* pCoordsA = pAtomA->getCoords();
                    vector3d* pCoordsB = pAtomB->getCoords();

                    coordsA[index][0] = pCoordsA->getX();
                    coordsA[index][1] = pCoordsA->getY();
                    coordsA[index][2] = pCoordsA->getZ();

                    coordsB[index][0] = pCoordsB->getX();
                    coordsB[index][1] = pCoordsB->getY();
                    coordsB[index][2] = pCoordsB->getZ();

                    index++;
                  }
                }
              }
            }

            double precentageIdentical = (double(nIdentical)/double(nMatches)) * 100.0;
            double precentageSimilar = (double(nIdentical + nSimilar)/double(nMatches)) * 100.0;
            double precentageDissimilar = (double(nDissimilar)/double(nMatches)) * 100.0;
            std::string asMessage = " \n Active site stats \n  Identicals = " + int2String(nIdentical) +
                    "/" + int2String(nMatches) +
                    " (" + double2String(precentageIdentical) + ") \n";

            asMessage += "  Similars = " + int2String(nIdentical + nSimilar) + "/" + int2String(nMatches) +
                    " (" + double2String(precentageSimilar) + ") \n";

            asMessage += "  Dissimilars = " + int2String(nDissimilar) + "/" + int2String(nMatches) +
                    " (" + double2String(precentageDissimilar) + ") \n";
            MTKpp::errorLogger.throwError("sequenceAligner::compareActiveSites", asMessage, MESSAGE);

            pSuperimpose = new superimpose(nMatches * nAtomsGroup);
            // Calculate the centers
            pSuperimpose->center(coordsA, centerA);
            pSuperimpose->center(coordsB, centerB);

            failure = pSuperimpose->getRotationMatrix(coordsA, coordsB,
                                    centerA, centerB, nMatches * nAtomsGroup, rotMat);

            // Get RMSD Value
            double rmsd = pSuperimpose->fit(coordsA, coordsB, nMatches * nAtomsGroup);

            errorMessage = " rmsd of active site " + inputFileContents[i][1] + " = " + double2String(rmsd);
            MTKpp::errorLogger.throwError("sequenceAligner::compareActiveSites", errorMessage, MESSAGE);

            // Save transformation
            sheet* pSheet = new sheet();
            std::string sheetName = "ActiveSite_" + inputFileContents[i][1];
            resultsMap[sheetName] = pSheet;
            table<double>* pDoubleTable_c1;
            pDoubleTable_c1 = pSheet->addTable();
            pDoubleTable_c1->setName("Template Center");
            pDoubleTable_c1->setSizes(1, 3);

            table<double>* pDoubleTable_c2;
            pDoubleTable_c2 = pSheet->addTable();
            pDoubleTable_c2->setName("Query Center");
            pDoubleTable_c2->setSizes(1, 3);

            for (unsigned int a = 0; a < 3; a++) {
              pDoubleTable_c1->setColumnLabel(a, labels[a]);
              pDoubleTable_c2->setColumnLabel(a, labels[a]);
              pDoubleTable_c1->setCellValue(0, a, centerA[a]);
              pDoubleTable_c2->setCellValue(0, a, centerB[a]);
            }

            table<double>* pDoubleTable_rotMat;
            pDoubleTable_rotMat = pSheet->addTable();
            pDoubleTable_rotMat->setName("Rotation Matrix");
            pDoubleTable_rotMat->setSizes(3, 3);

            for (int b = 0; b < 3; b++) {
              for (unsigned int a = 0; a < 3; a++) {
                pDoubleTable_rotMat->setCellValue(b, a, rotMat[b][a]);
              }
            }
          }
        }
        else {
          MTKpp::errorLogger.throwError("sequenceAligner::compareActiveSites", " Incorrect use ", MTK_ERROR);
        }
      }

      else {
        std::cout << " unknown command: " << inputFileContents[i][0] << std::endl;
      }

    }

    // Clean up
    delete pAtomTyper;
    delete pSeqAlign;
    delete pSuperimpose;
    delete pPdbParser;
    delete pCollection;

    return 0;
}
