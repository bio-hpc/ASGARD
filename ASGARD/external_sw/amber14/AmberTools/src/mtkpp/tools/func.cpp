/*!
   \file func.cpp

   \brief Determines the functional groups in ligands

   \author Martin B. Peters

   $Date: 2010/05/04 20:04:07 $
   $Revision: 1.8 $

   ----------------------------------------------------------------------------

   MTK++ - C++ package of modeling libraries.

   Copyright (C) 2005-2006  (see AUTHORS file for a list of contributors)

   ----------------------------------------------------------------------------
*/
#include "Utils/printHeader.h"

#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/torsion.h"
#include "Molecule/element.h"
#include "Molecule/stdLibrary.h"

#include "Parsers/stdLibParser.h"
#include "Parsers/commLineOptions.h"
#include "Parsers/elementParser.h"
#include "Parsers/molParser.h"
#include "Parsers/sdfParser.h"
#include "Parsers/paramParser.h"

#include "Log/errorHandler.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>

using namespace MTKpp;

int main (int argc, char **argv)
{
    std::string prog_name = "func";
    std::vector<std::string> authors;
    std::string author = "Martin B. Peters";
    authors.push_back(author);

    // 1. CREATE AN OBJECT
    commLineOptions *clo = new commLineOptions();

    // 2. SET PREFERENCES
    clo->noUsage();

    // 3. SET THE USAGE/HELP
    clo->addUsage( "  func: Determines the functional groups in ligands  \n" );
    clo->addUsage( "    usage: func [flags] [options]                    \n" );
    clo->addUsage( "  options:                                             " );
    clo->addUsage( "          -i input mol file                            " );
    clo->addUsage( "          -o output file                               " );
    clo->addUsage( "          -a log file                                \n" );
    clo->addUsage( "    flags:                                             " );
    clo->addUsage( "          -h help                                    \n" );

    // 4. SET THE OPTION STRINGS/CHARACTERS
    clo->setOption(  "in",     'i' );
    clo->setOption(  "out",    'o' );
    clo->setOption(  "log",    'a' );
    clo->setFlag  (  "help",   'h' );

    // 5. PROVIDE THE COMMANDLINE
    clo->processCommandArgs( argc, argv );

    clo->usageOn();

    // 6. GET THE VALUES
    if ( clo->getFlag( "help" ) || clo->getFlag( 'h' ) ) {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      return 0;
    }

    std::string molFile = "";
    std::string outFile = "";
    std::string logFile = "";

    if ( clo->getValue( "i" ) != 0 ) {
      molFile = clo->getValue( "l" );
    }
    else if ( clo->getValue( "in" ) != 0 ) {
      molFile =  clo->getValue( "in" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide an input mol file " << std::endl;
      return 0;
    }

    if ( clo->getValue( "o" ) != 0 ) {
      outFile = clo->getValue( "o" );
    }
    else if ( clo->getValue( "out" ) != 0 ) {
      outFile =  clo->getValue( "out" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide an output file " << std::endl;
      return 0;
    }

    if ( clo->getValue( "a" ) != 0 ) {
      logFile = clo->getValue( "a" );
    }
    else if ( clo->getValue( "log" ) != 0 ) {
      logFile =  clo->getValue( "log" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide a log file " << std::endl;
      return 0;
    }

    // 7. DONE
    delete clo;

    std::string errMessage = "";

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

    std::string AMBERHOME = getenv("AMBERHOME");
    if (AMBERHOME == "") {
      errMessage = " Can't find AMBERHOME ";
      MTKpp::errorLogger.throwError("superimposer", errMessage, MTK_ERROR);
      oLog.close();
      exit(0);
    }

    errMessage = " Create a New Collection";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);

    collection* pCollection = new collection();

    errMessage = " Read Element Data";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);
    elementParser* pElementParser = new elementParser(pCollection->pElements);
    //std::string elementXmlFile = MTKppDIR+"data/elements.xml";
    std::string elementXmlFile = AMBERHOME+"/dat/mtkpp/elements.xml";

    pElementParser->Read(elementXmlFile);

    errMessage = " Read MM parameters";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);

    pCollection->addParameters();
    paramParser* pParamParser = new paramParser(pCollection->getParameters());
    std::string parameterXmlFile = AMBERHOME + "/dat/mtkpp/parm94.xml";

    pParamParser->Read(parameterXmlFile);
    parameterXmlFile = AMBERHOME + "/dat/mtkpp/parm_gaff.xml";

    pParamParser->Read(parameterXmlFile);
    parameterXmlFile = AMBERHOME + "/dat/mtkpp/metals/metalParm.xml";

    pParamParser->Read(parameterXmlFile);

    errMessage = " Create a Standard Library";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);

    // STANDARD LIBRARIES
    pCollection->addStdLibrary();
    stdLibrary* pStdLibrary = pCollection->getStdLibrary();
    std::string stdLibXmlFile = "";
    stdLibParser* pStdLibParser = new stdLibParser(pStdLibrary, pCollection->getParameters());

    errMessage = " Reading Terminal Fragment Library";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);

    stdLibXmlFile = AMBERHOME+"/dat/mtkpp/fragLib/terminal.xml";
    pStdLibParser->Read(stdLibXmlFile);

    errMessage = " Reading 2PtLinkers Fragment Library";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);

    stdLibXmlFile = AMBERHOME+"/dat/mtkpp/fragLib/2PtLinkers.xml";
    pStdLibParser->Read(stdLibXmlFile);

    errMessage = " Reading 3PtLinkers Fragment Library";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);

    stdLibXmlFile = AMBERHOME+"/dat/mtkpp/fragLib/3PtLinkers.xml";
    pStdLibParser->Read(stdLibXmlFile);

    errMessage = " Reading 4PtLinkers Fragment Library";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);

    stdLibXmlFile = AMBERHOME+"/dat/mtkpp/fragLib/4PtLinkers.xml";
    pStdLibParser->Read(stdLibXmlFile);

    errMessage = " Reading 5MemRings Fragment Library";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);

    stdLibXmlFile = AMBERHOME+"/dat/mtkpp/fragLib/5MemRings.xml";
    pStdLibParser->Read(stdLibXmlFile);

    errMessage = " Reading 6MemRings Fragment Library";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);

    stdLibXmlFile = AMBERHOME+"/dat/mtkpp/fragLib/6MemRings.xml";
    pStdLibParser->Read(stdLibXmlFile);

    errMessage = " Reading 6MemRings Fragment Library";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);

    errMessage = " Reading fusedRings Fragment Library";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);

    stdLibXmlFile = AMBERHOME+"/dat/mtkpp/fragLib/fusedRings.xml";
    pStdLibParser->Read(stdLibXmlFile);

    errMessage = " Reading gt6MemRings Fragment Library";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);

    stdLibXmlFile = AMBERHOME+"/dat/mtkpp/fragLib/gt6MemRings.xml";
    pStdLibParser->Read(stdLibXmlFile);

    errMessage = " Reading hcaII Fragment Library";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);

    stdLibXmlFile = AMBERHOME+"/dat/mtkpp/fragLib/hcaII.xml";
    pStdLibParser->Read(stdLibXmlFile);

    //std::cout << " Total number of fragments in library = " << pStdLibrary->getNumberStdFrag() << std::endl;

    pStdLibrary->generateSimpleFP();
    pStdLibrary->generateAdjMatrices();
    pStdLibrary->generateAtomKinds();

    // CREATE PARSERS
    molParser* pMolParser = 0;
    sdfParser* pSdfParser = 0;

    // READ INPUT FILE
    pMolParser = new molParser();

    errMessage = " Reading Mol File";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);

    pMolParser->Read(molFile,pCollection);
    molecule* pMolecule = pCollection->getMolecule(1);

    std::string parametersFile = AMBERHOME + "/dat/mtkpp/hybridize/labute.txt";
    pMolecule->determineHybridizations(2, parametersFile);

    errMessage = " Determining rings";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);
    pMolecule->determineRings();
    pMolecule->kekulizeRings();

    std::cout << pMolecule->getName() << std::endl;

    std::string rI = pMolecule->getRingInfo();
    std::cout  << rI << std::endl;

    errMessage = " Protonating";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);
    //pMolecule->addHydrogens();

    errMessage = " Generating Simple Fingerprint";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);
    pMolecule->generateSimpleFP();

    errMessage = " Generating Adjacency Matrix";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);
    pMolecule->generateAdjMatrix();

    errMessage = " Determining Functional Groups";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);

    pMolecule->determineFunctionalGroups(pStdLibrary);

    // WRITE SDF FILE
    pSdfParser = new sdfParser();

    errMessage = " Write SDF file";
    MTKpp::errorLogger.throwError("func", errMessage, INFO);

    pSdfParser->Write(outFile, pCollection);

    errMessage = " Molecule fragment fp = " + pMolecule->getFragmentFP();
    MTKpp::errorLogger.throwError("func", errMessage, INFO);

    //pMolecule->generateFeatureDistMatrix();

    // - Clean up - //
    oLog.close();

    delete pCollection;
    delete pElementParser;
    delete pMolParser;
    delete pSdfParser;
    delete pStdLibParser;

    return 0;
}

