/*!
   \file hybrid.cpp

   \brief Determines the hybridizations of a molecule in a pdb file

   \author Martin B. Peters

   $Date: 2010/04/22 22:19:54 $
   $Revision: 1.8 $

   ----------------------------------------------------------------------------

   MTK++ - C++ package of modeling libraries.

   Copyright (C) 2005-2006  (see AUTHORS file for a list of contributors)

   ----------------------------------------------------------------------------
*/
#include "Utils/printHeader.h"

// - MOLECULE
#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/element.h"
#include "Molecule/connections.h"

// - MTK++ INCLUDE
#include "Utils/vector3d.h"
#include "Utils/constants.h"

// - PARSERS
#include "Parsers/elementParser.h"
#include "Parsers/pdbParser.h"
#include "Parsers/molParser.h"
#include "Parsers/mol2Parser.h"
#include "Parsers/commLineOptions.h"

// - LOG
#include "Log/errorHandler.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>

using namespace MTKpp;

/*!
   \brief Determines hybridizations, bond orders, and formal charges of a small molecule pdb file
   \param argc
   \param argv
   \return success
*/
int main (int argc, char **argv)
{
    std::string prog_name = "hybrid";
    std::vector<std::string> authors;
    std::string author = "Martin B. Peters";
    authors.push_back(author);

    // 1. CREATE AN OBJECT
    commLineOptions *clo = new commLineOptions();

    // 2. SET PREFERENCES
    clo->noUsage();

    // 3. SET THE USAGE/HELP
    clo->addUsage( "  hybrid: Determines bond orders of ligands          \n" );
    clo->addUsage( "    usage:  hybrid [flags] [options]                 \n" );
    clo->addUsage( "  options:                                             " );
    clo->addUsage( "          -i input pdb file                            " );
    clo->addUsage( "          -p parameters file                           " );
    clo->addUsage( "          -o output mol/pdb file                       " );
    clo->addUsage( "          -l log file                                \n" );
    clo->addUsage( "    flags:                                             " );
    clo->addUsage( "          -h help                                    \n" );

    // 4. SET THE OPTION STRINGS/CHARACTERS
    clo->setOption( "inFile",   'i' );
    clo->setOption( "parmFile", 'p' );
    clo->setOption( "logFile",  'l' );
    clo->setOption( "outFile",  'o' );

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
    std::string parametersFile = "";
    std::string outputFile = "";
    std::string logFile = "hybrid.log";

    if ( clo->getValue( "i" ) != 0 ) {
      inputFile = clo->getValue( "i" );
    }
    else if ( clo->getValue( "inFile" ) != 0 ) {
      inputFile =  clo->getValue( "inFile" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide an input file " << std::endl;
      return 0;
    }

    if ( clo->getValue( "p" ) != 0 ) {
      parametersFile = clo->getValue( "p" );
    }
    else if ( clo->getValue( "parmFile" ) != 0 ) {
      parametersFile =  clo->getValue( "parmFile" );
    }

    if ( clo->getValue( "o" ) != 0 ) {
      outputFile = clo->getValue( "o" );
    }
    else if ( clo->getValue( "outFile" ) != 0 ) {
      outputFile =  clo->getValue( "outFile" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide an output file " << std::endl;
      return 0;
    }

    if ( clo->getValue( "l" ) != 0 ) {
      logFile = std::string(clo->getValue( "l" ));
    }
    else if ( clo->getValue( "logFile" ) != 0 ) {
      logFile =  std::string(clo->getValue( "logFile" ));
    }

    // 6. OPEN LOG FILE
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
    std::string AMBERHOME = getenv("AMBERHOME");
    collection* pCollection = new collection();

    // READ ELEMENT DATA
    elementParser* pElementParser = new elementParser(pCollection->pElements);
    std::string elementXmlFile = AMBERHOME+"/dat/mtkpp/elements.xml";
    pElementParser->Read(elementXmlFile);

    if (parametersFile == "") {
      parametersFile = AMBERHOME+"/dat/mtkpp/hybridize/labute.txt";
    }

    // CREATE PARSERS
    pdbParser* pPdbParser = 0;
    molParser* pMolParser = 0;
    mol2Parser* pMol2Parser = 0;

    pPdbParser = new pdbParser();
    pMolParser = new molParser();
    pMol2Parser = new mol2Parser();

    // READ INPUT FILE
    pPdbParser->Read(inputFile, pCollection);

    molecule* pMolecule = pCollection->getMolecule(1);

    // Assign bonds using the Meng Algorithm
    connections* pConnections = new connections(pCollection);
    pConnections->assignBonds(pMolecule);

    // Determine atom Hybridizations and formal charges, and bond orders
    // using the Labute Algorithm
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

    // Set molecule in solution
    pMolecule->setInSolution();

    // Add hydrogen atoms to the molecule
    pMolecule->addHydrogens();

    // WRITE OUTPUT FILE
    if (extName(outputFile) == "mol") {
      pMolParser->Write(outputFile, pMolecule);
    }
    else if (extName(outputFile) == "pdb") {
      pPdbParser->Write(outputFile, pMolecule);
    }
    else if (extName(outputFile) == "mol2") {
      pMol2Parser->Write(outputFile, pMolecule);
    }
    else {
      std::cout << " Unknown type " << extName(outputFile) << std::endl;
    }

    // - Clean up - //
    delete pCollection;
    delete pPdbParser;
    delete pMolParser;
    delete pMol2Parser;

    time (&endTime);
    int diffTime = (int) difftime(endTime, startTime);
    errorMessage = " hybrid Exited Normally After " 
                    + int2String(diffTime) + " Seconds ";
    MTKpp::errorLogger.throwError("hybrid", errorMessage, 0);
    oLog.close();
    return 0;
}
