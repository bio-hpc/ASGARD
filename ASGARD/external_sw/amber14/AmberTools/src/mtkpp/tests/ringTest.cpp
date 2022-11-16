/*!
   \file ringTest.cpp
   \brief Tests the ring detection functionality in MTK++
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
#include "Utils/vector3d.h"

#include "Log/errorHandler.h"

#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/atom.h"
#include "Molecule/bond.h"
#include "Molecule/torsion.h"
#include "Molecule/element.h"
#include "Molecule/connections.h"

#include "Parsers/elementParser.h"
#include "Parsers/molParser.h"
#include "Parsers/sdfParser.h"
#include "Parsers/pdbParser.h"
#include "Parsers/commLineOptions.h"
#include "Parsers/StringManip.h"
#include "Parsers/inputParser.h"

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
   \brief Tests MTK++'s ring detection functionality
   \param argc
   \param argv
   \return success
*/
int main (int argc, char **argv)
{
    std::string prog_name = "ringTest";
    std::vector<std::string> authors;
    std::string author = "Martin B. Peters";
    authors.push_back(author);

    // 1. CREATE AN OBJECT
    commLineOptions *clo = new commLineOptions();

    // 2. SET PREFERENCES
    clo->noUsage();

    // 3. SET THE USAGE/HELP
    clo->addUsage( "  ringTest                                           \n" );
    clo->addUsage( "    usage:  ringTest [flags] [options]               \n" );
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

    std::string logFile = "ringTest.log";
    std::string outFile = "ringTest.out";
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

    collection* pCollection = new collection();

    // Read Element Data
    elementParser* pElementParser = new elementParser(pCollection->pElements);
    std::string elementXmlFile = AMBERHOME + "/dat/mtkpp/elements.xml";
    pElementParser->Read(elementXmlFile);

    // Labute Parameter file
    std::string parametersFile = AMBERHOME + "/dat/mtkpp/hybridize/labute.txt";

    // Create Parsers
    pdbParser* pPdbParser = 0;
    sdfParser* pSdfParser = 0;

    pPdbParser = new pdbParser();
    pSdfParser = new sdfParser();

    // Read ligand list
    std::string pdbFilePath = AMBERHOME + "/AmberTools/test/mtkpp/hybridizeData/";
    std::string pdbListFile = pdbFilePath + "list.txt";
    std::string curPdbFile = "";
    std::vector<std::string> pdbList;
    int f = readListFile(pdbListFile, pdbList);
    if (f) {
      errorMessage = " Failed to read PDB file list ... exiting ";
      errorLogger.throwError("hybridizeTest", errorMessage, 1);
      exit(1);
    }

    molecule* pMolecule = 0;

    for (unsigned int i = 0; i < pdbList.size(); i++) {
      curPdbFile = pdbFilePath + pdbList[i];

      // Read pdb file
      pPdbParser->Read(curPdbFile, pCollection);
      pMolecule = pCollection->getMolecule(1);

      // Assign bonds using the Meng Algorithm
      connections* pConnections = new connections(pCollection);
      pConnections->assignBonds(pMolecule);

      // Determine atom Hybridizations and formal charges,
      // and bond orders using the Labute Algorithm
      pMolecule->determineHybridizations(2, parametersFile);

      // Assign angles, torsions, and impropers
      pConnections->assignAngles(pMolecule);
      pConnections->assignTorsions(pMolecule);
      pConnections->assignImpropers(pMolecule);
      delete pConnections;

      // Determine ring structure
      pMolecule->determineRings();

      // If rings are present test if there is aromatic rings
      pMolecule->kekulizeRings();
      oOut << pMolecule->getName() << std::endl;
      std::string rI = pMolecule->getRingInfo();
      oOut << rI << std::endl;

      pCollection->delMolecule(pMolecule);
      pMolecule = 0;
    }

    time (&endTime);
    int diffTime = (int) difftime(endTime, startTime);

    errorMessage = " Exited Normally After " + int2String(diffTime) + " Seconds ";
    errorLogger.throwError("ringTest", errorMessage, 4);

    oLog.close();
    oOut.close();

    return 0;
}
