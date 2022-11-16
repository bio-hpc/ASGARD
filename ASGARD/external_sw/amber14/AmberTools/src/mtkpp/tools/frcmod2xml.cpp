/*!
   \file frcmod2xml.cpp

   \brief Convert AMBER frcmod file to MTK++ xml parameter file

   \author Martin B. Peters

   $Date: 2010/04/22 22:19:54 $
   $Revision: 1.5 $

   ----------------------------------------------------------------------------

   MTK++ - C++ package of modeling libraries.

   Copyright (C) 2005-2006  (see AUTHORS file for a list of contributors)

   ----------------------------------------------------------------------------
*/
#include "Utils/printHeader.h"

// - MOLECULE
#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/submolecule.h"
#include "Molecule/atom.h"
#include "Molecule/element.h"
#include "Molecule/bond.h"
#include "Molecule/ring.h"
#include "Molecule/connections.h"
#include "Molecule/stdLibrary.h"
#include "Molecule/stdGroup.h"
#include "Molecule/stdFrag.h"
#include "Molecule/functionalize.h"

// - MTK++ INCLUDE
#include "Utils/vector3d.h"

// - PARSERS
#include "Parsers/elementParser.h"
#include "Parsers/pdbParser.h"
#include "Parsers/molParser.h"
#include "Parsers/stdLibParser.h"
#include "Parsers/frcmodParser.h"
#include "Parsers/paramParser.h"
#include "Parsers/commLineOptions.h"

#include "Log/errorHandler.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <vector>
#include <map>
#include <algorithm>

using namespace MTKpp;

/*!
   \brief Convert AMBER frcmod file to MTK++ xml parameter file
*/
int main (int argc, char **argv)
{
    std::string prog_name = "frcmod2xml";
    std::vector<std::string> authors;
    std::string author = "Martin B. Peters";
    authors.push_back(author);
    //printHeader(std::cout, prog_name, authors);

    // 1. CREATE AN OBJECT
    commLineOptions *clo = new commLineOptions();

    // 2. SET PREFERENCES
    clo->noUsage();

    // 3. SET THE USAGE/HELP
    clo->addUsage( "  frcmod2xml: Converts frcmod file to MTK++ xml file \n" );
    clo->addUsage( "    usage:  frcmod2xml [flags] [options]             \n" );
    clo->addUsage( "  options:                                             " );
    clo->addUsage( "          -i frcmod file                               " );
    clo->addUsage( "          -o parameter xml file                        " );
    clo->addUsage( "          -n parameters name                           " );
    clo->addUsage( "          -a log file                                \n" );
    clo->addUsage( "    flags:                                             " );
    clo->addUsage( "          -h help                                      " );

    // 4. SET THE OPTION STRINGS/CHARACTERS
    clo->setOption( "frcmodFile", 'i' );
    clo->setOption( "paramXml",   'o' );
    clo->setOption( "paramName",  'n' );
    clo->setOption( "logFile",    'a' );
    clo->setFlag  ( "help",       'h' );

    // 5. PROVIDE THE COMMANDLINE 
    clo->processCommandArgs( argc, argv );

    clo->usageOn();

    // 6. GET THE VALUES
    if ( clo->getFlag( "help" ) || clo->getFlag( 'h' ) ) {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      return 0;
    }

    std::string frcmodFile   = "";
    std::string paramXml     = "";
    std::string paramName    = "";
    std::string logFile      = "";

    if ( clo->getValue( "i" ) != 0 ) {
      frcmodFile = clo->getValue( "i" );
    }
    else if ( clo->getValue( "frcmodFile" ) != 0 ) {
      frcmodFile =  clo->getValue( "frcmodFile" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide a frcmod file " << std::endl;
      return 0;
    }

    if ( clo->getValue( "o" ) != 0 ) {
      paramXml = clo->getValue( "o" );
    }
    else if ( clo->getValue( "paramXml" ) != 0 ) {
      paramXml =  clo->getValue( "paramXml" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide a param xml file " << std::endl;
      return 0;
    }

    if ( clo->getValue( "n" ) != 0 ) {
      paramName = clo->getValue( "n" );
    }
    else if ( clo->getValue( "paramName" ) != 0 ) {
      paramName =  clo->getValue( "paramName" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide a parameters name " << std::endl;
      return 0;
    }

    if ( clo->getValue( "a" ) != 0 ) {
      logFile = clo->getValue( "a" );
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

    //std::string MTKppDIR = getenv("MTKppDIR");
    collection* pCollection = new collection();

    errMessage = " Read Parameters ";
    MTKpp::errorLogger.throwError("frcmod2xml", errMessage, INFO);
    pCollection->addParameters();
    parameters* pParameters = pCollection->getParameters();
    frcmodParser* pFrcmodParser = new frcmodParser(pParameters, paramName);
    pFrcmodParser->Read(frcmodFile);

    errMessage = " Write Parameters ";
    MTKpp::errorLogger.throwError("frcmod2xml", errMessage, INFO);
    paramParser* pParamParser = new paramParser(pParameters);
    pParamParser->Write(paramXml, paramName);

    // - Clean up
    oLog.close();
    delete pFrcmodParser;
    delete pParamParser;
    delete pCollection;
    return 0;
}

