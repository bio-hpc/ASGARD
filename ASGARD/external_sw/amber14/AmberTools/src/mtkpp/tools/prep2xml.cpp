/*!
   \file prep2xml.cpp

   \brief Convert AMBER prep file to MTK++ xml library

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
#include "Molecule/hydrophobize.h"

// - MTK++ INCLUDE
#include "Utils/vector3d.h"

// - PARSERS
#include "Parsers/elementParser.h"
#include "Parsers/stdLibParser.h"
#include "Parsers/prepParser.h"
#include "Parsers/paramParser.h"
#include "Parsers/sdfParser.h"
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
   \brief Convert AMBER prep file to MTK++ xml library
*/
int main (int argc, char **argv)
{
    std::string prog_name = "prep2xml";
    std::vector<std::string> authors;
    std::string author = "Martin B. Peters";
    authors.push_back(author);

    // 1. CREATE AN OBJECT
    commLineOptions *clo = new commLineOptions();

    // 2. SET PREFERENCES
    clo->noUsage();

    // 3. SET THE USAGE/HELP
    clo->addUsage( "  prep2xml: Convert a prepin file to an MTK++ xml file\n" );
    clo->addUsage( "    usage:  prep2xml -i <file> -o <file> -g <string> -f <string> -n <string> \\" );
    clo->addUsage( "                     -a <file> [-l <string>]         \n" );
    clo->addUsage( "    usage:  prep2xml -h                              \n" );
    clo->addUsage( "  options:                                             " );
    clo->addUsage( "          -i input prep file                           " );
    clo->addUsage( "          -o output XML file                           " );
    clo->addUsage( "          -g group name                                " );
    clo->addUsage( "          -f residue 3-letter code                     " );
    clo->addUsage( "          -n residue description                       " );
    clo->addUsage( "          -a log file                                  " );
    clo->addUsage( "          -l hybridize name                          \n" );
    clo->addUsage( "    flags:                                             " );
    clo->addUsage( "          -c convert only                              " );
    clo->addUsage( "          -h help                                    \n" );

    // 4. SET THE OPTION STRINGS/CHARACTERS
    clo->setOption( "prepFile",  'i' );
    clo->setOption( "libXml",    'o' );
    clo->setOption( "group",     'g' );
    clo->setOption( "frag",      'f' );
    clo->setOption( "mol",       'n' );
    clo->setOption( "hyb",       'l' );
    clo->setOption( "log",       'a' );

    clo->setFlag  ( "convert",   'c' );
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

    std::string prepFile       = "";
    std::string libXmlFile     = "";
    std::string groupName      = "";
    std::string fragName       = "";
    std::string molName        = "";
    std::string logFile        = "";
    bool bConvertOnly          = 0;

    if ( clo->getFlag( "convert" ) || clo->getFlag( 'c' ) ) {
      bConvertOnly = 1;
    }

    std::string AMBERHOME = getenv("AMBERHOME");
    std::string parametersFile = AMBERHOME + "/dat/mtkpp/hybridize/labute.txt";

    if ( clo->getValue( "i" ) != 0 ) {
      prepFile = clo->getValue( "i" );
    }
    else if ( clo->getValue( "prepFile" ) != 0 ) {
      prepFile =  clo->getValue( "prepFile" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide a prep file " << std::endl;
      return 0;
    }

    if ( clo->getValue( "o" ) != 0 ) {
      libXmlFile = clo->getValue( "o" );
    }
    else if ( clo->getValue( "libXml" ) != 0 ) {
      libXmlFile =  clo->getValue( "libXml" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide a lib xml file " << std::endl;
      return 0;
    }

    if ( clo->getValue( "g" ) != 0 ) {
      groupName = clo->getValue( "g" );
    }
    else if ( clo->getValue( "group" ) != 0 ) {
      groupName =  clo->getValue( "group" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide a group name " << std::endl;
      return 0;
    }

    if ( clo->getValue( "f" ) != 0 ) {
      fragName = clo->getValue( "f" );
    }
    else if ( clo->getValue( "frag" ) != 0 ) {
      fragName =  clo->getValue( "frag" );
    }
/*    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide a fragment name " << std::endl;
      return 0;
    }
*/
    if ( clo->getValue( "l" ) != 0 ) {
      parametersFile = clo->getValue( "l" );
    }
    else if ( clo->getValue( "hyb" ) != 0 ) {
      parametersFile =  clo->getValue( "hyb" );
    }

    if ( clo->getValue( "n" ) != 0 ) {
      molName = clo->getValue( "n" );
    }
    else if ( clo->getValue( "mol" ) != 0 ) {
      molName =  clo->getValue( "mol" );
    }
/*
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide a molecule name " << std::endl;
      return 0;
    }
*/
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

    std::string errMessage = "";

    collection* pCollection = new collection();

    errMessage = " Read Element Data";
    MTKpp::errorLogger.throwError("prep2xml", errMessage, INFO);
    elementParser* pElementParser = new elementParser(pCollection->pElements);
    std::string elementXmlFile = AMBERHOME + "/dat/mtkpp/elements.xml";
    pElementParser->Read(elementXmlFile);
    delete pElementParser;

    errMessage = " Read Parameters ";
    MTKpp::errorLogger.throwError("prep2xml", errMessage, INFO);
    pCollection->addParameters();
    parameters* pParameters = pCollection->getParameters();
    paramParser* pParamParser = new paramParser(pCollection->getParameters());
    std::string parameterXmlFile = AMBERHOME + "/dat/mtkpp/parm94.xml";
    pParamParser->Read(parameterXmlFile);
    parameterXmlFile = AMBERHOME + "/dat/mtkpp/parm_gaff.xml";
    pParamParser->Read(parameterXmlFile);
    delete pParamParser;

    errMessage = " Read Libraries ";
    MTKpp::errorLogger.throwError("prep2xml", errMessage, INFO);

    pCollection->addStdLibrary();
    stdLibrary* pStdLib = pCollection->getStdLibrary();
    if (!pStdLib) return 1;

    std::string stdLibXmlFile = "";
    stdLibParser* pStdLibParser = new stdLibParser(pStdLib, pParameters);

    stdLibXmlFile = AMBERHOME + "/dat/mtkpp/fragLib/terminal.xml";
    pStdLibParser->Read(stdLibXmlFile);

    stdLibXmlFile = AMBERHOME + "/dat/mtkpp/fragLib/2PtLinkers.xml";
    pStdLibParser->Read(stdLibXmlFile);

    stdLibXmlFile = AMBERHOME + "/dat/mtkpp/fragLib/3PtLinkers.xml";
    pStdLibParser->Read(stdLibXmlFile);

    stdLibXmlFile = AMBERHOME + "/dat/mtkpp/fragLib/4PtLinkers.xml";
    pStdLibParser->Read(stdLibXmlFile);

    stdLibXmlFile = AMBERHOME + "/dat/mtkpp/fragLib/5MemRings.xml";
    pStdLibParser->Read(stdLibXmlFile);

    stdLibXmlFile = AMBERHOME + "/dat/mtkpp/fragLib/6MemRings.xml";
    pStdLibParser->Read(stdLibXmlFile);

    stdLibXmlFile = AMBERHOME + "/dat/mtkpp/fragLib/fusedRings.xml";
    pStdLibParser->Read(stdLibXmlFile);

    stdLibXmlFile = AMBERHOME + "/dat/mtkpp/fragLib/cores.xml";
    pStdLibParser->Read(stdLibXmlFile);

    stdLibXmlFile = AMBERHOME + "/dat/mtkpp/fragLib/hcaII.xml";
    pStdLibParser->Read(stdLibXmlFile);

    pStdLib->generateSimpleFP();
    pStdLib->generateAdjMatrices();
    pStdLib->generateAtomKinds();

    // Create molecule
    molecule* pMolecule = pCollection->addMolecule();
    pMolecule->setName(fragName);
    pMolecule->setMolId(pCollection->getNumberMolecules());

    submolecule* pSubMolecule = pMolecule->addSubMolecule();
    pSubMolecule->setName(fragName);
    pSubMolecule->setSubMolId(pMolecule->getNumSubMolecules());

    //
    // Create new group and fragment
    //
    stdGroup* pStdGroup = pStdLib->addStdGroup();
    pStdGroup->setName(groupName);

    if (!pStdGroup) {
      std::cout << " Incorrect use of prep2xml " << std::endl;
      exit(1);
    }

    /////////
    if (fragName == "") {
      prepParser* pPrepParser = new prepParser();
      pPrepParser->Read(prepFile, pStdGroup);
      delete pPrepParser;

      pStdLibParser->Write(libXmlFile, groupName);
      delete pStdLibParser;

      // - Clean up - //
      delete pCollection;
      return 0;
    }
    ///////////

    stdFrag* pStdFrag = pStdGroup->addStdFrag();
    if (!pStdFrag) {
      std::cout << " Incorrect use of prep2xml " << std::endl;
      exit(1);
    }

    pStdFrag->setSymbol(fragName);       // 3L code
    pStdFrag->setCode("P2XML"+fragName); // 8L code
    pStdFrag->setName(molName);          // long name
    pStdFrag->setType("m");              // fragment type

    prepParser* pPrepParser = new prepParser();
    pPrepParser->Read(prepFile, pStdFrag);
    delete pPrepParser;

    if (bConvertOnly) {
      pStdLibParser->Write(libXmlFile, groupName);
      delete pStdLibParser;

      // - Clean up - //
      delete pCollection;
      return 0;
    }

    int f = pStdFrag->generateCoordinates();
    if (f) {
      std::cout << " Incorrect use of prep2xml " << std::endl;
      exit(1);
    }
    std::vector<vector3d*> stdAtomCoords = pStdFrag->getCoordinates();
    std::vector<stdAtom*> stdAtoms = pStdFrag->getStdAtomList();
    std::vector<stdBond*> stdBonds = pStdFrag->getStdBondList();
    std::vector<stdLoop*> stdLoops = pStdFrag->getStdLoopList();

    // stdAtoms
    for (unsigned int i = 0; i < stdAtoms.size(); i++) {
      atom* pAtom = pSubMolecule->addAtom();
      pAtom->setElement(pCollection->pElements->getElement(stdAtoms[i]->atSymbol));
      pAtom->setFileID(i+1);
      pAtom->setName(stdAtoms[i]->identity);
      pAtom->setCoords(stdAtomCoords[i]->getX(),stdAtomCoords[i]->getY(),stdAtomCoords[i]->getZ());
      pAtom->setStdAtom(stdAtoms[i]);
    }

    // stdBond's
    for (unsigned int i = 0; i < stdBonds.size(); i++) {
      if ((stdBonds[i]->atom1 < 0) or (stdBonds[i]->atom2 < 0)) continue;
      int at1 = stdBonds[i]->atom1;
      int at2 = stdBonds[i]->atom2;
      int bondType = stdBonds[i]->type;
      int bondStereo = stdBonds[i]->stereo;
      int bondTopology = stdBonds[i]->topology;

      atom* pBondAtom1 = pMolecule->getAtom(at1, 1, 0);
      atom* pBondAtom2 = pMolecule->getAtom(at2, 1, 0);
      Bond* pBond = pMolecule->addBond(pBondAtom1, pBondAtom2, bondType, bondStereo, bondTopology, 0.0);
      if (!pBond) {
        std::cout << " Error in prep2xml ... exiting " << std::endl;
        exit(0);
      }
      pBondAtom1->addBondedAtom(pBondAtom2);
      pBondAtom2->addBondedAtom(pBondAtom1);
    }

    // stdLoop's
    for (unsigned int i = 0; i < stdLoops.size(); i++) {
      if ((stdLoops[i]->atom1 < 0) or (stdLoops[i]->atom2 < 0)) continue;
      int at1 = stdLoops[i]->atom1;
      int at2 = stdLoops[i]->atom2;
      int bondType = stdLoops[i]->type;
      int bondStereo = stdLoops[i]->stereo;
      int bondTopology = 1;

      atom* pBondAtom1 = pMolecule->getAtom(at1, 1, 0);
      atom* pBondAtom2 = pMolecule->getAtom(at2, 1, 0);
      Bond* pBond = pMolecule->addBond(pBondAtom1, pBondAtom2, bondType, bondStereo, bondTopology, 0.0);
      if (!pBond) {
        std::cout << " Error in prep2xml ... exiting " << std::endl;
        exit(0);
      }
      pBondAtom1->addBondedAtom(pBondAtom2);
      pBondAtom2->addBondedAtom(pBondAtom1);
    }

    errMessage = " Determine Hybridization And Bond Orders:Labute Algorithm ";
    MTKpp::errorLogger.throwError("prep2xml", errMessage, INFO);
    pMolecule->determineHybridizations(2, parametersFile);

    errMessage = " Assign angles, torsions, and impropers ";
    MTKpp::errorLogger.throwError("prep2xml", errMessage, INFO);
    connections* pConnections = new connections(pCollection);
    pConnections->assignAngles(pMolecule);
    pConnections->assignTorsions(pMolecule);
    pConnections->assignImpropers(pMolecule);
    delete pConnections;

    errMessage = " Determine ring structure ";
    MTKpp::errorLogger.throwError("prep2xml", errMessage, INFO);
    pMolecule->determineRings();

    //puts("  prep2xml::Kekulize Rings");
    //pMolecule->kekulizeRings();

    errMessage = " Set molecule in solution ";
    MTKpp::errorLogger.throwError("prep2xml", errMessage, INFO);
    pMolecule->setInSolution();

    errMessage = " Generate Simple Fingerprint ";
    MTKpp::errorLogger.throwError("prep2xml", errMessage, INFO);
    pMolecule->generateSimpleFP();

    errMessage = " Generate Adjacency Matrix for functional group recognition ";
    MTKpp::errorLogger.throwError("prep2xml", errMessage, INFO);
    pMolecule->generateAdjMatrix();

    errMessage = " Determine Functional Groups ";
    MTKpp::errorLogger.throwError("prep2xml", errMessage, INFO);
    pMolecule->determineFunctionalGroups(pStdLib);

    //sdfParser* pSdfParser = new sdfParser();
    //pSdfParser->Write("testing.sdf", pCollection);
    //delete pSdfParser;

    errMessage = " Determine Hydrophobic Regions";
    MTKpp::errorLogger.throwError("prep2xml", errMessage, INFO);
    f = pMolecule->determineHydrophobicGroups();
    if (f) {
      std::cout << " Error in determining hydrophobic Groups ... exiting " << std::endl;
      exit(0);
    }

    errMessage = " Update Library file";
    MTKpp::errorLogger.throwError("prep2xml", errMessage, INFO);
    std::vector<atom*> atomList = pMolecule->getAtomList();
    std::map<int, Bond*> bondMap = pMolecule->getBondMap();
    typedef std::map<int, Bond*>::iterator BondMapIterator;

    typedef std::map<stdAtom*, atom*> funcGroupAtomsMap;

    std::vector<ring*> ringList = pMolecule->getRings();
    std::vector<funcGroup*> funcGroupList = pMolecule->getFunctionalGroups();
    std::vector<hydrophobe*> hydrophobeList = pMolecule->getHydrophobicGroups();
    // atoms
    for (unsigned int i = 0; i < atomList.size(); i++) {
      stdAtoms[i]->kind = atomList[i]->getType();
      //std::cout << atomList[i]->getIndex() << " " << atomList[i]->getName() << " " << atomList[i]->getType() << std::endl;
    }

    // bonds
    for (BondMapIterator b = bondMap.begin(); b != bondMap.end(); b++) {
      Bond* pBond = b->second;
      stdBond* pStdBond = pStdFrag->getStdBond(pBond->atom1->getStdAtom(), pBond->atom2->getStdAtom());
      if (pStdBond) {
        pStdBond->type = pBond->type;
        pStdBond->stereo = pBond->stereo;
        pStdBond->topology = pBond->topology;

        atom* pAt1 = pBond->atom1;
        atom* pAt2 = pBond->atom2;
        std::string at1Symbol = pAt1->getElementSymbol();
        std::string at2Symbol = pAt2->getElementSymbol();
        if (at1Symbol == "H" or at2Symbol == "H") {
          if (at1Symbol == "N" or at1Symbol == "O" or at1Symbol == "S" or
              at2Symbol == "N" or at2Symbol == "O" or at2Symbol == "S") {
            pStdBond->kind = 1;
          }
          else {
            pStdBond->kind = pBond->kind;
          }
        }
        else {
          pStdBond->kind = pBond->kind;
        }
      }
      else {
        stdLoop* pStdLoop = pStdFrag->getStdLoop(pBond->atom1->getStdAtom(), pBond->atom2->getStdAtom());
        if (pStdLoop) {
          pStdLoop->type = pBond->type;
          pStdLoop->stereo = pBond->stereo;
        }
        else {
          std::cout << " Error in prep2xml ... exiting " << std::endl;
          exit(1);
        }
      }
    }

    // rings
    for (unsigned int i = 0; i < ringList.size(); i++) {
      stdRing* pStdRing = pStdFrag->addStdRing();
      for (unsigned int j = 0; j < ringList[i]->atoms.size(); j++) {
        pStdRing->atoms.push_back(ringList[i]->atoms[j]->getStdAtom()->index);
      }
      pStdRing->size      = ringList[i]->size;
      pStdRing->planar    = ringList[i]->planar;
      pStdRing->aromatic  = ringList[i]->aromatic;
      pStdRing->hetero    = ringList[i]->hetero;
      pStdRing->nHetero   = ringList[i]->nHetero;
      pStdRing->nNitrogen = ringList[i]->nNitrogen;
      pStdRing->nOxygen   = ringList[i]->nOxygen;
      pStdRing->nSulfur   = ringList[i]->nSulfur;
    }

    // functional groups
    for (unsigned int i = 0; i < funcGroupList.size(); i++) {
      stdFuncGroup* pStdFuncGroup = pStdFrag->addStdFuncGroup();

      pStdFuncGroup->groupName = funcGroupList[i]->pStdFrag->getParent()->getName();
      pStdFuncGroup->fragName = funcGroupList[i]->pStdFrag->getSymbol();

      std::vector<stdAtom*> fragAtoms = funcGroupList[i]->pStdFrag->getStdAtomList();
      for (unsigned int j = 0; j < fragAtoms.size(); j++) {
        pStdFuncGroup->atoms.push_back(funcGroupList[i]->atomMap[fragAtoms[j]]->getIndex());
      }
    }

    // Hydrophobic groups
    for (unsigned int i = 0; i < hydrophobeList.size(); i++) {
      stdFeature* pStdFeature = pStdFrag->addStdFeature();
      pStdFeature->name = "HPB";
      for (unsigned int j = 0; j < hydrophobeList[i]->atoms.size(); j++) {
        pStdFeature->atoms.push_back(hydrophobeList[i]->atoms[j]->getStdAtom()->index);
      }
    }

    pStdLibParser->Write(libXmlFile, groupName);
    delete pStdLibParser;

    // - Clean up - //
    delete pCollection;
    return 0;
}

