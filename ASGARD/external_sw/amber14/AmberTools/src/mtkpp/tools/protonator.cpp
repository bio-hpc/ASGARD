/*!
   \file protonator.cpp

   \brief Protonates mol/sdf/pdb files

   \author Martin B. Peters

   $Date: 2010/04/22 22:19:54 $
   $Revision: 1.11 $

   ----------------------------------------------------------------------------

   MTK++ - C++ package of modeling libraries.

   Copyright (C) 2005-2006  (see AUTHORS file for a list of contributors)

   ----------------------------------------------------------------------------
*/
#include "Utils/printHeader.h"

// - MOLECULE
#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/atom.h"
#include "Molecule/bond.h"
#include "Molecule/torsion.h"
#include "Molecule/element.h"
#include "Molecule/connections.h"
#include "Molecule/parameters.h"
#include "Molecule/atomTyper.h"
#include "Molecule/protonate.h"
#include "Utils/vector3d.h"

// - PARSERS
#include "Parsers/elementParser.h"
#include "Parsers/molParser.h"
#include "Parsers/sdfParser.h"
#include "Parsers/pdbParser.h"

// - PARAMETER, ATOM TYPE, STANDARD LIBRARY
#include "Parsers/paramParser.h"
#include "Parsers/stdLibParser.h"
#include "Parsers/commLineOptions.h"

#include "Log/errorHandler.h"

// - TIME
#include "time.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <vector>
#include <map>
#include <algorithm>

using namespace MTKpp;

/*!
   \brief Protonates mol/sdf/pdb files
   \param argc
   \param argv
   \return success
*/
int main (int argc, char **argv)
{
    std::string prog_name = "protonator";
    std::vector<std::string> authors;
    std::string author = "Martin B. Peters";
    authors.push_back(author);

    // 1. CREATE AN OBJECT
    commLineOptions *clo = new commLineOptions();

    // 2. SET PREFERENCES
    clo->noUsage();

    // 3. SET THE USAGE/HELP
    clo->addUsage( "  protonator: Adds Hs to Pro/Lig/Wat Molecules       \n" );
    clo->addUsage( "    usage:  protonator [flags] [options]             \n" );
    clo->addUsage( "  options:                                             " );
    clo->addUsage( "          -i input file                                " );
    clo->addUsage( "          -o output file                               " );
    clo->addUsage( "          -l log file                                \n" );
    clo->addUsage( "    flags:                                             " );
    clo->addUsage( "          -h help                                    \n" );

    // 4. SET THE OPTION STRINGS/CHARACTERS
    clo->setOption( "inFile",   'i' );
    clo->setOption( "outFile",  'o' );
    clo->setOption( "logFile",  'l' );

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
    std::string outputFile = "";
    std::string logFile = "";

    if ( clo->getValue( "i" ) != NULL ) {
      inputFile = clo->getValue( "i" );
    }
    else if ( clo->getValue( "inFile" ) != NULL ) {
      inputFile =  clo->getValue( "inFile" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide an input file " << std::endl;
      return 0;
    }

    if ( clo->getValue( "o" ) != NULL ) {
      outputFile = clo->getValue( "o" );
    }
    else if ( clo->getValue( "outFile" ) != NULL ) {
      outputFile =  clo->getValue( "outFile" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide an output file " << std::endl;
      return 0;
    }

    if ( clo->getValue( "l" ) != 0 ) {
      logFile = clo->getValue( "l" );
    }
    else if ( clo->getValue( "logFile" ) != 0 ) {
      logFile =  clo->getValue( "logFile" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide a log file " << std::endl;
      return 0;
    }

    // 7. DONE
    delete clo;

    // Open log file
    std::ofstream oLog;
    oLog.open(logFile.c_str());

    if (!oLog) {
      std::cout << "\nUnable to open log file"
                << "\nFilename = " << logFile << std::endl;
      exit(1);
    }

    // Set errorLog stream to the log file
    MTKpp::errorLogger.setStream(&oLog);

    // Print MTK++ copyright message
    printHeader(oLog, prog_name, authors);

    std::string errMessage = "";

    errMessage = " Set up";
    MTKpp::errorLogger.throwError("protonator", errMessage, INFO);

    std::string AMBERHOME = getenv("AMBERHOME");
    collection* pCollection = new collection();
    molecule* pMolecule = 0;

    elementParser* pElementParser = new elementParser(pCollection->pElements);
    std::string elementXmlFile = AMBERHOME + "/dat/mtkpp/elements.xml";
    pElementParser->Read(elementXmlFile);
    delete pElementParser;

    std::string inputType = extName(inputFile);
    std::string outputType = extName(outputFile);

    int end   = outputFile.length();
    int slash = outputFile.find_last_of("/");
    std::string outfile_name = outputFile.substr(slash+1,(end-slash-5));

    end   = outfile_name.length();
    int dot = outputFile.find_last_of(".");
    outfile_name = outputFile.substr(0,dot);

    if (inputType == "pdb") {

      errMessage = " Input type equals pdb, Reading parameters";
      MTKpp::errorLogger.throwError("protonator", errMessage, INFO);

      pCollection->addParameters();
      paramParser* pParamParser = new paramParser(pCollection->getParameters());
      std::string parameterXmlFile = AMBERHOME + "/dat/mtkpp/parm94.xml";
      pParamParser->Read(parameterXmlFile);
      parameterXmlFile = AMBERHOME + "/dat/mtkpp/parm_gaff.xml";
      pParamParser->Read(parameterXmlFile);
      parameterXmlFile = AMBERHOME + "/dat/mtkpp/metals/metalParm.xml";
      pParamParser->Read(parameterXmlFile);
      delete pParamParser;

      //pCollection->getParameters()->printAtomTypes();

      errMessage = " Reading library files";
      MTKpp::errorLogger.throwError("protonator", errMessage, INFO);

      pCollection->addStdLibrary();
      stdLibParser* pStdLibParser = new stdLibParser(
                        pCollection->getStdLibrary(),
                        pCollection->getParameters());

      std::string stdLibXmlFile = "";
      stdLibXmlFile = AMBERHOME + "/dat/mtkpp/aminont94.xml";
      pStdLibParser->Read(stdLibXmlFile);

      stdLibXmlFile = AMBERHOME + "/dat/mtkpp/amino94.xml";
      pStdLibParser->Read(stdLibXmlFile);

      stdLibXmlFile = AMBERHOME + "/dat/mtkpp/aminoct94.xml";
      pStdLibParser->Read(stdLibXmlFile);

      stdLibXmlFile = AMBERHOME + "/dat/mtkpp/metals/metals.xml";
      pStdLibParser->Read(stdLibXmlFile);
      delete pStdLibParser;

      errMessage = " Reading pdb files";
      MTKpp::errorLogger.throwError("protonator", errMessage, INFO);

      pdbParser* pPdbParser = new pdbParser();
      pPdbParser->Read(inputFile, pCollection);
      delete pPdbParser;

      errMessage = " Assign disulfide bonds ";
      MTKpp::errorLogger.throwError("protonator", errMessage, INFO);

      connections* pConnections = new connections(pCollection);
      pConnections->assignDisulfideBonds();

      errMessage = " Atom typing ";
      MTKpp::errorLogger.throwError("protonator", errMessage, INFO);

      atomTyper* pAtomTyper = new atomTyper();
      pAtomTyper->atomTypeByLib(pCollection);
      delete pAtomTyper;

      errMessage = " Defining Connections ";
      MTKpp::errorLogger.throwError("protonator", errMessage, INFO);
      pConnections->run();
      delete pConnections;
    }
    else if (inputType == "mol") {
      errMessage = " Input type equals mol";
      MTKpp::errorLogger.throwError("protonator", errMessage, INFO);

      molParser* pMolParser = new molParser();
      pMolParser->Read(inputFile, pCollection);
      delete pMolParser;
    }
    else if (inputType == "sdf") {
      errMessage = " Input type equals sdf";
      MTKpp::errorLogger.throwError("protonator", errMessage, INFO);
      sdfParser* pSdfParser = new sdfParser();
      pSdfParser->Read(inputFile, pCollection);
      delete pSdfParser;
    }
    else {
      errMessage = " Unknown file type " + inputType;
      MTKpp::errorLogger.throwError("protonator", errMessage, MTK_ERROR);
      return 0;
    }

    // Doing pdb file separately because of water molecules
    if (inputType == "pdb") {
      protonate* pProtonate = new protonate(pCollection);
      pProtonate->run();
      delete pProtonate;
    }
    else {
      std::vector<molecule*> molList = pCollection->getMoleculeList();
      std::string molFileName = outputFile;
      for (unsigned int i = 0; i < molList.size(); i++) {
        pMolecule = molList[i];
        std::string parametersFile = AMBERHOME + "/dat/mtkpp/hybridize/labute.txt";
        pMolecule->determineHybridizations(2, parametersFile);
        pMolecule->determineRings();
        pMolecule->addHydrogens();

        if (outputType == "mol") {
          molParser* pMolParser = new molParser();
          if (molList.size() > 1) {
            std::stringstream ss1;
            ss1 << i+1;
            std::string num = ss1.str().c_str(); 
            molFileName = outfile_name + num + ".mol";
          }
          pMolParser->Write(molFileName, pCollection, i+1);
          delete pMolParser;
        }
      }
    }

    if (outputType == "sdf") {
      sdfParser* pSdfParser = new sdfParser();
      pSdfParser->Write(outputFile, pCollection);
      delete pSdfParser;
    }
    else if (outputType == "pdb") {
      pdbParser* pPdbParser = new pdbParser();
      pPdbParser->Write(outputFile, pCollection);
      delete pPdbParser;
    }
    // - Clean up
    oLog.close();
    delete pCollection;

    return 0;
}

