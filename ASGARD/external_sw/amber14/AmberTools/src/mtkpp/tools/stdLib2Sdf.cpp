/*!
   \file stdLib2Sdf.cpp

   \brief Converts an xml library file into a sd file

   \author Martin B. Peters

   $Date: 2010/04/22 22:19:54 $
   $Revision: 1.10 $

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
#include "Molecule/bond.h"
#include "Molecule/torsion.h"
#include "Molecule/element.h"
#include "Molecule/connections.h"
#include "Molecule/stdLibrary.h"
#include "Molecule/stdGroup.h"
#include "Molecule/stdFrag.h"

#include "Utils/vector3d.h"

// - PARSERS
#include "Parsers/elementParser.h"
#include "Parsers/sdfParser.h"
#include "Parsers/stdLibParser.h"
#include "Parsers/paramParser.h"

// - COMMAND LINE OPTIONS
#include "Parsers/commLineOptions.h"

// - Log
#include "Log/errorHandler.h"

// temp
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
   \brief Converts an xml library file into a sdf file
   \param argc
   \param argv
   \return success
*/
int main (int argc, char *argv[])
{
    std::string prog_name = "stdLib2Sdf";
    std::vector<std::string> authors;
    std::string author = "Martin B. Peters";
    authors.push_back(author);

    // 1. CREATE AN OBJECT
    commLineOptions *clo = new commLineOptions();

    // 2. SET PREFERENCES
    clo->noUsage();

    // 3. SET THE USAGE/HELP
    clo->addUsage( "  stdLib2Sdf: Converts fragment libraries into sdf   \n" ); 
    clo->addUsage( "  usage: stdLib2Sdf  [flags] [options]               \n" );
    clo->addUsage( "  options:                                             " );
    clo->addUsage( "          -i stdLib xml file                           " );
    clo->addUsage( "          -l param xml file                            " );
    clo->addUsage( "          -o sdf file                                  " );
    clo->addUsage( "          -a log file                                \n" );
    clo->addUsage( "    flags:                                             " );
    clo->addUsage( "          -f reference molecules only                \n" );
    clo->addUsage( "          -h help                                    \n" );

    // 4. SET THE OPTION STRINGS/CHARACTERS
    clo->setFlag  ( "help",   'h' );
    clo->setFlag  ( "ref",    'f' );
    clo->setOption( "stdLib", 'i' );
    clo->setOption( "param",  'l' );
    clo->setOption( "sdf",    'o' );
    clo->setOption( "log",    'a' );

    // 5. PROVIDE THE COMMANDLINE
    clo->processCommandArgs( argc, argv );

    clo->usageOn();

    bool refOnly = false;

    // 6. GET THE VALUES
    if ( clo->getFlag( "help" ) || clo->getFlag( 'h' ) ) {
      clo->printUsage();
      return 0;
    }

    if ( clo->getFlag( "ref" ) || clo->getFlag( 'f' ) ) {
      refOnly = true;
    }

    std::string stdLibFile;
    std::string paramFile;
    std::string sdfFile;
    std::string logFile = "stdLib2Sdf.log";

    if ( clo->getValue( "i" ) != NULL ) {
      stdLibFile = clo->getValue( "i" );
    }
    else if ( clo->getValue( "stdLib" ) != NULL ) {
      stdLibFile =  clo->getValue( "stdLib" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      return 0;
    }

    if ( clo->getValue( "l" ) != NULL ) {
      paramFile = clo->getValue( "l" );
    }
    else if ( clo->getValue( "param" ) != NULL ) {
      paramFile =  clo->getValue( "param" );
    }

    if ( clo->getValue( "o" ) != NULL ) {
      sdfFile = clo->getValue( "o" );
    }
    else if ( clo->getValue( "sdf" ) != NULL ) {
      sdfFile =  clo->getValue( "sdf" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
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

    //std::cout << " stdLib File = " << stdLibFile << std::endl;
    //std::cout << " sdf File    = " << sdfFile << std::endl;

    std::string AMBERHOME = getenv("AMBERHOME");

    collection* pCollection = new collection();

    // READ ELEMENT DATA
    elementParser* pElementParser = new elementParser(pCollection->pElements);
    std::string elementXmlFile = AMBERHOME + "/dat/mtkpp/elements.xml";
    pElementParser->Read(elementXmlFile);

    pCollection->addParameters();
    paramParser* pParamParser = new paramParser(pCollection->getParameters());
    if (pParamParser) {
      std::string paramXmlFile = AMBERHOME + "/dat/mtkpp/parm94.xml";
      pParamParser->Read(paramXmlFile);
      paramXmlFile = AMBERHOME + "/dat/mtkpp/parm_gaff.xml";
      pParamParser->Read(paramXmlFile);

      if (paramFile != "") {
        pParamParser->Read(paramFile);
      }
    }
    delete pParamParser;
    parameters* pParms = pCollection->getParameters();

    // STANDARD LIBRARIES
    pCollection->addStdLibrary();
    stdLibrary* pStdLibrary = pCollection->getStdLibrary();

    stdLibParser* pStdLibParser = new stdLibParser(pCollection, pStdLibrary, pParms);
    pStdLibParser->Read(stdLibFile);

    molecule* pMolecule = 0;
    submolecule* pSubMolecule = 0;
    atom* pAtom = 0;
    atom* pBondAtom1 = 0;
    atom* pBondAtom2 = 0;
    Bond* pBond = 0;
    int success = 0;
    vector3d* com = new vector3d(0.0, 0.0, 0.0);

    if (!refOnly) {
      // stdGroup's
      std::vector<stdGroup*> groupList = pStdLibrary->getStdGroupList();
      for (unsigned int i = 0; i < groupList.size(); i++) {
        //std::cout << groupList[i]->getName() << std::endl;

        // stdFrag's
        std::vector<stdFrag*> fragList = groupList[i]->getStdFragList();
        for (unsigned int j = 0; j < fragList.size(); j++) {
          //std::cout << fragList[j]->getName() << std::endl;
          success = fragList[j]->generateCoordinates();
          if (success == 0) {
            std::vector<vector3d*> fragCoords = fragList[j]->getCoordinates();

            pMolecule  = pCollection->addMolecule();
            pMolecule->setName(fragList[j]->getName());
            pMolecule->setMolId(pCollection->getNumberMolecules());

            pSubMolecule = pMolecule->addSubMolecule();
            pSubMolecule->setName(fragList[j]->getName());
            pSubMolecule->setSubMolId(pMolecule->getNumSubMolecules());

            // stdAtom's
            std::vector<stdAtom*> atomList = fragList[j]->getStdAtomList();
            for (unsigned int k = 0; k < atomList.size(); k++) {
              pAtom = pSubMolecule->addAtom();
              pAtom->setElement(pCollection->pElements->getElement(atomList[k]->atSymbol));
              pAtom->setFileID(k+1);
              pAtom->setCoords(fragCoords[k]->getX(),fragCoords[k]->getY(),fragCoords[k]->getZ());
            }

            // stdBond's
            std::vector<stdBond*> bondList = fragList[j]->getStdBondList();
            for (unsigned int k = 0; k < bondList.size(); k++) {
              if ((bondList[k]->atom1 < 0) or (bondList[k]->atom2 < 0)) continue;
              int at1 = bondList[k]->atom1;
              int at2 = bondList[k]->atom2;
              int bondType = bondList[k]->type;
              int bondStereo = bondList[k]->stereo;
              int bondTopology = bondList[k]->topology;

              pBondAtom1 = pMolecule->getAtom(at1, 1, 0);
              pBondAtom2 = pMolecule->getAtom(at2, 1, 0);
              pBond = pMolecule->addBond(pBondAtom1, pBondAtom2, bondType, bondStereo, bondTopology, 0.0);
              if (!pBond) {
                std::cout << " Error in adding bonds in stdLib2Sdf " << std::endl;
                exit(1);
              }
              pBondAtom1->addBondedAtom(pBondAtom2);
              pBondAtom2->addBondedAtom(pBondAtom1);
            }

            // stdLoop's
            std::vector<stdLoop*> loopList = fragList[j]->getStdLoopList();
            for (unsigned int k = 0; k < loopList.size(); k++) {
              if ((loopList[k]->atom1 < 0) or (loopList[k]->atom2 < 0)) continue;
              int at1 = loopList[k]->atom1;
              int at2 = loopList[k]->atom2;
              int bondType = loopList[k]->type;
              int bondStereo = loopList[k]->stereo;
              int bondTopology = 1;

              pBondAtom1 = pMolecule->getAtom(at1, 1, 0);
              pBondAtom2 = pMolecule->getAtom(at2, 1, 0);
              pBond = pMolecule->addBond(pBondAtom1, pBondAtom2, bondType, bondStereo, bondTopology, 0.0);
              if (!pBond) {
                std::cout << " Error in adding (loop) bonds stdLib2Sdf " << std::endl;
                exit(1);
              }
              pBondAtom1->addBondedAtom(pBondAtom2);
              pBondAtom2->addBondedAtom(pBondAtom1);
            }
            pMolecule->moveCenterOfMass(com);
          }
        }
      }
    }
    else {
      std::vector<molecule*> molList = pCollection->getMoleculeList();
      for (unsigned int x = 0; x < molList.size(); x++) {
        molList[x]->moveCenterOfMass(com);
      }
    }

    // CREATE PARSERS
    sdfParser* pSdfParser = 0;

    // READ INPUT FILE
    pSdfParser = new sdfParser();
    pSdfParser->Write(sdfFile, pCollection);

    oLog.close();

    // Clean up
    delete pElementParser;
    delete pStdLibParser;
    delete pSdfParser;
    delete pCollection;
}
