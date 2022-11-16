/*!
   \file fileFormatsTest.cpp
   \brief Tests file reading functionality in MTK++
   \author Martin Peters

   $Date: 2010/08/19 20:03:09 $
   $Revision: 1.6 $

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
#include "Molecule/seqAlign.h"

#include "Statistics/sheet.h"
#include "Statistics/table.h"

// - PARSERS
#include "Parsers/StringManip.h"
#include "Parsers/elementParser.h"
#include "Parsers/paramParser.h"
#include "Parsers/dMParser.h"
#include "Parsers/stdLibParser.h"
#include "Parsers/mtkppParser.h"
#include "Parsers/sdfParser.h"
#include "Parsers/pdbParser.h"
#include "Parsers/commLineOptions.h"
#include "Parsers/inputParser.h"

#include "Parsers/mol2Parser.h"
#include "Parsers/molParser.h"
#include "Parsers/pamParser.h"
#include "Parsers/xyzParser.h"

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
   \brief Tests file reading functionality
   \param argc
   \param argv
   \return success
*/
int main (int argc, char **argv)
{
    std::string prog_name = "fileFormatsTest";
    std::vector<std::string> authors;
    std::string author = "Martin B. Peters";
    authors.push_back(author);

    // 1. CREATE AN OBJECT
    commLineOptions *clo = new commLineOptions();

    // 2. SET PREFERENCES
    clo->noUsage();

    // 3. SET THE USAGE/HELP
    clo->addUsage( "  fileFormatsTest                                    \n" );
    clo->addUsage( "    usage:  fileFormatsTest [flags] [options]        \n" );
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

    std::string logFile = "fileFormatsTest.log";
    std::string outFile = "fileFormatsTest.out";
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

      pParamParser->Write("temp_parm94.xml", "parm94");

      std::string paramXmlFile2 = AMBERHOME+"/dat/mtkpp/parm_gaff.xml";
      pParamParser->Read(paramXmlFile2);

      std::string paramXmlFile3 = AMBERHOME+"/AmberTools/test/mtkpp/fileFormatsData/FM6.frcmod.xml";
      pParamParser->Read(paramXmlFile3);

      std::string paramXmlFile4 = AMBERHOME+"/AmberTools/test/mtkpp/fileFormatsData/1A5T_params.xml";
      pParamParser->Read(paramXmlFile4);

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

        pStdLibParser->Write("temp_amino94.xml");

        std::string libXmlFile2 = AMBERHOME+"/AmberTools/test/mtkpp/fileFormatsData/FM6.xml";
        pStdLibParser->Read(libXmlFile2);

        std::string libXmlFile3 = AMBERHOME+"/dat/mtkpp/metals/metals.xml";
        pStdLibParser->Read(libXmlFile3);

        std::string libXmlFile4 = AMBERHOME+"/AmberTools/test/mtkpp/fileFormatsData/1A5T_stdMol.xml";
        pStdLibParser->Read(libXmlFile4);

        pStdLibParser->Write("temp_1A5T_stdMol.xml", "1A5T");

        delete pStdLibParser;
      }
    }

    // Read table Data
    dMParser* pdMParser = new dMParser();
    if (pdMParser) {
      std::string dmXmlFile = AMBERHOME+"/AmberTools/test/mtkpp/fileFormatsData/tables.xml";
      sheet* s = new sheet();
      pdMParser->read(s, dmXmlFile);

      pdMParser->write(s, "temp_table.xml", false);

      delete pdMParser;
    }

    // Read/Write Protein PDB file

    // Read gaussian output file

    // Write AMBER prmtop file

    // Read/Write AMBER ac file

    // Read/Write sd file
    sdfParser* pSdfParser = new sdfParser();
    std::string sdfFile = AMBERHOME+"/AmberTools/test/mtkpp/fileFormatsData/hcaII.sdf";
    pSdfParser->Read(sdfFile, pCollection);

    pSdfParser->Write("temp_hcaII.sdf", pCollection);

    // Read pdb file
    std::string pdbFile = AMBERHOME + "/AmberTools/test/mtkpp/hybridizeData/1A42_lig.pdb";

    pdbParser* pPdbParser = new pdbParser();
    pPdbParser->Read(pdbFile, pCollection);
    molecule* pMolecule = pCollection->getLastAddedMolecule();

    // Assign bonds, angles, torsions, and impropers
    connections* pConnections = new connections(pCollection);
    pConnections->assignBonds(pMolecule);
    pConnections->assignAngles(pMolecule);
    pConnections->assignTorsions(pMolecule);
    pConnections->assignImpropers(pMolecule);
    delete pConnections;
    pPdbParser->Write("temp_1A42_lig.pdb", pMolecule);

    // Write mtk++ state file
    mtkppParser* pMTKppParser = new mtkppParser();
    pMTKppParser->Write("temp_state.xml", pCollection);

    // Read mtk++ state file
    // todo

    // Read mol2 file
    mol2Parser* pMol2Parser = new mol2Parser();
    std::string mol2File = AMBERHOME+"/AmberTools/test/mtkpp/fileFormatsData/benzene.mol2";
    pMol2Parser->Read(mol2File, pCollection);
    pMolecule = pCollection->getLastAddedMolecule();

    // Write mol2 file
    pMol2Parser->Write("temp_benzene.mol2", pMolecule);

    // Read xyz file
    xyzParser* pXyzParser = new xyzParser();
    std::string xyzFile = AMBERHOME+"/AmberTools/test/mtkpp/fileFormatsData/benzene.xyz";
    pXyzParser->Read(xyzFile, pCollection);
    pMolecule = pCollection->getLastAddedMolecule();

    // Write xyz file
    pXyzParser->Write("temp_benzene.xyz", pMolecule);

    // Read mol file
    molParser* pMolParser = new molParser();
    std::string molFile = AMBERHOME+"/AmberTools/test/mtkpp/fileFormatsData/benzene.mol";
    pMolParser->Read(molFile, pCollection);
    pMolecule = pCollection->getLastAddedMolecule();

    // Write mol file
    pMolParser->Write("temp_benzene.mol", pMolecule);

    // Read/Write pam file
    pamParser* pPamParser = new pamParser();
    if (pPamParser) {
      if (!pPamParser->getError()) {
        std::string pamFile = AMBERHOME+"/dat/mtkpp/PAM/PAM250";
        seqAlign* pSeqAlign = new seqAlign();
        pPamParser->Read(pamFile, pSeqAlign);
      }
    }

    time (&endTime);
    int diffTime = (int) difftime(endTime, startTime);

    errorMessage = " Exited Normally After " + int2String(diffTime) + " Seconds ";
    errorLogger.throwError("fileFormatsTest", errorMessage, 4);

    oLog.close();
    oOut.close();

    return 0;
}
