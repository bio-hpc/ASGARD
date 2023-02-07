/*!
   \file mmTest.cpp
   \brief Tests the MM functionality in MTK++
   \author Martin Peters

   $Date: 2010/08/11 21:20:18 $
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
#include "Utils/printHeader.h"
#include "Log/errorHandler.h"

#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/atom.h"
#include "Molecule/bond.h"
#include "Molecule/atomTyper.h"
#include "Molecule/connections.h"

#include "Utils/vector3d.h"

// - PARSERS
#include "Parsers/elementParser.h"
#include "Parsers/pdbParser.h"
#include "Parsers/paramParser.h"
#include "Parsers/stdLibParser.h"
#include "Parsers/commLineOptions.h"
#include "Parsers/inputParser.h"

#include "MM/amber.h"

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
   \brief Tests MTK++'s MM functionality
   \param argc
   \param argv
   \return success
*/
int main (int argc, char **argv)
{
    std::string prog_name = "mmTest";
    std::vector<std::string> authors;
    std::string author = "Martin B. Peters";
    authors.push_back(author);

    // 1. CREATE AN OBJECT
    commLineOptions *clo = new commLineOptions();

    // 2. SET PREFERENCES
    clo->noUsage();

    // 3. SET THE USAGE/HELP
    clo->addUsage( "  mmTest                                             \n" );
    clo->addUsage( "    usage:  mmTest [flags] [options]                 \n" );
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

    std::string logFile = "mmTest.log";
    std::string outFile = "mmTest.out";
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

    int calcBond      = 1;
    int calcAngle     = 1;
    int calcTorsion   = 1;
    int calcImproper  = 1;
    int calcNonBonded = 1;
    int calcHBond     = 0;

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

    collection* pCollection = new collection();
    pCollection->addStdLibrary();
    pCollection->addParameters();

    // Read Element Data
    elementParser* pElementParser = new elementParser(pCollection->pElements);
    std::string elementXmlFile = AMBERHOME+"/dat/mtkpp/elements.xml";
    pElementParser->Read(elementXmlFile);

    // Read parameter Data
    paramParser* pParamParser = new paramParser(pCollection->getParameters());
    if (pParamParser) {
      std::string paramXmlFile = AMBERHOME+"/dat/mtkpp/parm94.xml";
      pParamParser->Read(paramXmlFile);

      std::string paramXmlFile2 = AMBERHOME+"/dat/mtkpp/parm_gaff.xml";
      pParamParser->Read(paramXmlFile2);

      std::string paramXmlFile5 = AMBERHOME+"/dat/mtkpp/metals/metalParm.xml";
      pParamParser->Read(paramXmlFile5);

      delete pParamParser;
    }

    // Read standard library
    stdLibrary* pStdLibrary = pCollection->getStdLibrary();
    if (pStdLibrary) {
      stdLibParser* pStdLibParser = new stdLibParser(pCollection, pStdLibrary, pCollection->getParameters());
      if (pStdLibParser) {
        std::string libXmlFile = AMBERHOME+"/dat/mtkpp/amino94.xml";
        pStdLibParser->Read(libXmlFile);

        std::string libXmlFile2 = AMBERHOME+"/dat/mtkpp/aminoct94.xml";
        pStdLibParser->Read(libXmlFile2);

        std::string libXmlFile3 = AMBERHOME+"/dat/mtkpp/aminont94.xml";
        pStdLibParser->Read(libXmlFile3);

        std::string libXmlFile4 = AMBERHOME+"/dat/mtkpp/metals/metals.xml";
        pStdLibParser->Read(libXmlFile4);

        delete pStdLibParser;
      }
    }

    // Read pdb list
    std::string pdbFilePath = AMBERHOME + "/AmberTools/test/mtkpp/mmData/";
    std::string pdbListFile = pdbFilePath + "list.txt";
    std::vector<std::string> pdbList;
    int f = readListFile(pdbListFile, pdbList);
    if (f) {
      std::cout << " Failed to read PDB file list ... exiting " << std::endl;
      exit(1);
    }

    pdbParser* pPdbParser = new pdbParser();
    atomTyper* pAtomTyper;
    pAtomTyper = new atomTyper(0);
    connections* pConnections = new connections();

    for (unsigned int i = 0; i < pdbList.size(); i++) {
      std::string pdbFile = pdbFilePath + pdbList[i];
      //std::cout<< pdbFile << std::endl;

      pPdbParser->Read(pdbFile, pCollection);
      std::string inputBaseName = baseName(pdbFile);

      pAtomTyper->atomTypeByLib(pCollection);
      pConnections->run(pCollection);

      MTKpp::amber* pAmber = new MTKpp::amber();
      pAmber->setPotential(calcBond, calcAngle, calcTorsion, calcImproper,
                           calcNonBonded, calcHBond);

      // INITIALIZE SOME VARIABLES
      int failure = 0;
      int nAtoms = 0;
      int nUniqueTypes = 0;
      int nExcludedAtoms = 0;
      int nExcluded14Atoms = 0;
      int nResidues = 0;

      int nBonds = 0;
      int nAngles = 0;
      int nTorsions = 0;
      int nImpropers = 0;

      // START ALLOCATION //
      // SET NUMBER OF ATOMS
      nAtoms = pCollection->getNumAtoms();
      pAmber->setNumAtoms(nAtoms);

      // SET NUMBER OF UNIQUE ATOM TYPES
      nUniqueTypes = pCollection->getNumUniqueAtomTypes();
      failure = pAmber->setNumTypes(nUniqueTypes);
      if (failure) {
        std::cout << " Error setting number of unique atom types" << std::endl;
        return 1;
      }

      // SET NUMBER OF EXCLUDED ATOMS
      nExcludedAtoms = pCollection->getNumExcludedAtoms();
      nExcluded14Atoms = pCollection->getNumExcluded14Atoms();
      failure = pAmber->setExcludedSize(nExcludedAtoms, nExcluded14Atoms);
      if (failure) {
        std::cout << " Error setting number of excluded atoms" << std::endl;
        return 1;
      }

      // SET NUMBER OF RESIDUES
      nResidues = pCollection->getNumberSubMolecules();
      failure = pAmber->setNumResidues(nResidues);
      if (failure) {
        std::cout << " Error setting number of residues" << std::endl;
        return 1;
      }

      // BONDS
      if (calcBond) {
        nBonds = pCollection->getNumBonds();
        pAmber->setNumBonds(nBonds);
      }

      // ANGLES
      if (calcAngle) {
        nAngles = pCollection->getNumAngles();
        pAmber->setNumAngles(nAngles);
      }

      // TORSIONS
      if (calcTorsion) {
        nTorsions = pCollection->getNumMMTorsions();
        pAmber->setNumTorsions(nTorsions);
      }

      // IMPROPERS
      if (calcImproper) {
        nImpropers = pCollection->getNumMMImpropers();
        pAmber->setNumImpropers(nImpropers);
      }
      // END ALLOCATION //

      // START LOAD DATA //
      // LOAD COORDINATES FROM MOLECULE TO MM
      failure = pCollection->getCoordinates(pAmber->getCoords());
      if (failure) {
        std::cout << " Error loading coordinates " << std::endl;
        return 0;
      }

      // LOAD ATOM CHARGES FROM MOLECULE TO MM
      failure = pCollection->getMMCharges(pAmber->getCharges());
      if (failure) {
        std::cout << " Error loading charges " << std::endl;
        return 0;
      }

      // LOAD ATOM SYMBOLS FROM MOLECULE TO MM
      failure = pCollection->getAtomSymbols(pAmber->getSymbols());
      if (failure) {
        std::cout << " Error loading symbols " << std::endl;
        return 0;
      }

      // LOAD ATOM INTEGER TYPES FROM MOLECULE TO MM
      failure = pCollection->getAtomTypes(pAmber->getIntTypes());
      if (failure) {
        std::cout << " Error loading atom integer types " << std::endl;
        return 0;
      }

      // LOAD ATOM CHARACTER TYPES FROM MOLECULE TO MM
      failure = pCollection->getAtomTypes(pAmber->getCharTypes());
      if (failure) {
        std::cout << " Error loading atom character types " << std::endl;
        return 0;
      }

      // LOAD BONDS
      failure = pCollection->getBonds(pAmber->getBonds());
      if (failure) {
        std::cout << " Error loading bonds " << std::endl;
        exit(0);
      }

      failure = pCollection->getBondParams(pAmber->getBondParams());
      if (failure) {
        std::cout << " Error loading bonds " << std::endl;
        exit(0);
      }

      failure = pCollection->getAngles(pAmber->getAngles());
      if (failure) {
        std::cout << " Error loading angles " << std::endl;
        exit(0);
      }

      failure = pCollection->getAngleParams(pAmber->getAngleParams());
      if (failure) {
        std::cout << " Error loading angles " << std::endl;
        exit(0);
      }

      failure = pCollection->getMMTorsions(pAmber->getTorsions(),
                                           pAmber->getTorsionParams());
      if (failure) {
        std::cout << " Error loading torsions " << std::endl;
        exit(0);
      }

      failure = pCollection->getMMImpropers(pAmber->getImpropers(),
                                            pAmber->getImproperParams());
      if (failure) {
        std::cout << " Error loading impropers " << std::endl;
        exit(0);
      }

      // LOAD L-J PARAMETERS FROM MOLECULE TO MM
      failure = pCollection->getLJParams(pAmber->getR6Params(), pAmber->getR12Params());
      if (failure) {
        std::cout << " Error loading LJ parameters " << std::endl;
        exit(0);
      }

      // LOAD THE NUMBER OF EXCLUDED ATOMS(1-2, 1-3, 1-4) FOR EACH ATOM
      failure = pCollection->getNumExcludedAtoms(pAmber->getNumExcluded());
      if (failure) {
        std::cout << " Error loading number of excluded atoms " << std::endl;
        exit(0);
      }

      // LOAD THE EXCLUDED ATOMS(1-2, 1-3, 1-4) FOR EACH ATOM
      failure = pCollection->getExcludedAtoms(pAmber->getExcluded());
      if (failure) {
        std::cout << " Error loading excluded atoms  " << std::endl;
        exit(0);
      }

      // LOAD THE NUMBER OF EXCLUDED 1-4 ATOMS FROM EACH ATOM
      failure = pCollection->getNumExcluded14Atoms(pAmber->getNumExcluded14());
      if (failure) {
        std::cout << " Error loading number of 1-4 excluded atoms " << std::endl;
        exit(0);
      }

      // LOAD THE EXCLUDED 1-4 ATOMS FROM EACH ATOM
      failure = pCollection->getExcluded14Atoms(pAmber->getExcluded14());
      if (failure) {
        std::cout << " Error loading 1-4 excluded atoms " << std::endl;
        exit(0);
      }
      // END LOAD DATA //

      // Structures xyz coordinates
      double *xyz;
      xyz = pAmber->getCoords();

      double d = 0.0;
      // Test Structure
      for (int i = 0; i < nAtoms; i++) {
        for (int j = 0; j < nAtoms; j++) {
          if (i == j) continue;
          d = sqrt(pow((xyz[i*3  ] - xyz[j*3  ]),2) +
                   pow((xyz[i*3+1] - xyz[j*3+1]),2) +
                   pow((xyz[i*3+2] - xyz[j*3+2]),2) );
          if (d < 0.1) {
            std::cout << " Two Atoms are less than 0.1 Ang Apart ... exiting" << std::endl;
            exit(0);
          }
        }
      }

      pAmber->calcEnergy();

      molecule* pMolecule = pCollection->getLastAddedMolecule();

      pAmber->printEnergy2(oOut, " ### " + pMolecule->getName() + " Energies ", "");

      pCollection->delMolecule(pMolecule);

      pCollection->setAtomIndex(1);
    }

    // - Clean up - //
    delete pAtomTyper;
    delete pConnections;
    delete pPdbParser;

    delete pCollection;

    time (&endTime);
    int diffTime = (int) difftime(endTime, startTime);

    errorMessage = " Exited Normally After " + int2String(diffTime) + " Seconds ";
    errorLogger.throwError("mmTest", errorMessage, 4);

    oLog.close();
    oOut.close();

    return 0;
}
