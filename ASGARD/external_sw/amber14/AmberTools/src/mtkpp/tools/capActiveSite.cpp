/*!
   \file capActiveSite.cpp

   \brief Caps active site using a cutoff with NME/ACE residues

   \author Martin B. Peters

   $Date: 2010/04/22 22:19:54 $
   $Revision: 1.6 $

   ----------------------------------------------------------------------------

   MTK++ - C++ package of modeling libraries.

   Copyright (C) 2005-2006  (see AUTHORS file for a list of contributors)

   ----------------------------------------------------------------------------
*/
#include "Utils/printHeader.h"
#include "Utils/vector3d.h"
#include "Utils/constants.h"

// - MOLECULE
#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/submolecule.h"
#include "Molecule/atom.h"
#include "Molecule/bond.h"
#include "Molecule/torsion.h"
#include "Molecule/element.h"
#include "Molecule/stdLibrary.h"
#include "Molecule/stdFrag.h"
#include "Molecule/connections.h"
#include "Molecule/atomTyper.h"
#include "Molecule/metalCenter.h"

// - PARSERS
#include "Parsers/elementParser.h"
#include "Parsers/pdbParser.h"
#include "Parsers/molParser.h"
#include "Parsers/paramParser.h"
#include "Parsers/stdLibParser.h"

#include "Log/errorHandler.h"

// - COMMAND LINE OPTIONS
#include "Parsers/commLineOptions.h"

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

int main (int argc, char *argv[])
{
    std::string prog_name = "capActiveSite";
    std::vector<std::string> authors;
    std::string author = "Martin B. Peters";
    authors.push_back(author);

    // 1. CREATE AN OBJECT
    commLineOptions *clo = new commLineOptions();

    // 2. SET PREFERENCES
    clo->noUsage();

    // 3. SET THE USAGE/HELP
    clo->addUsage( "  capActiveSite: Caps active site using a cutoff     \n" );
    clo->addUsage( "    usage: capActiveSite [flags] [options]           \n" );
    clo->addUsage( "  options:                                             " );
    clo->addUsage( "          -r receptor pdb file                         " );
    clo->addUsage( "          -l ligand mol/pdb file                       " );
    clo->addUsage( "          -c distance cutoff [10.0]                    " );
    clo->addUsage( "          -o output pdb file                           " );
    clo->addUsage( "          -a log file                                \n" );
    clo->addUsage( "    flags:                                             " );
    clo->addUsage( "          -h help                                    \n" );

    // 4. SET THE OPTION STRINGS/CHARACTERS
    clo->setOption(  "rec",    'r' );
    clo->setOption(  "lig",    'l' );
    clo->setOption(  "cut",    'c' );
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

    std::string recFile = "";
    std::string ligFile  = "";
    double cutoff = 10.0;
    std::string outFile = "";
    std::string logFile = "capActiveSite.log";

    if ( clo->getValue( "l" ) != 0 ) {
      ligFile = clo->getValue( "l" );
    }
    else if ( clo->getValue( "lig" ) != 0 ) {
      ligFile =  clo->getValue( "lig" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide a ligand name " << std::endl;
      return 0;
    }

    if ( clo->getValue( "r" ) != 0 ) {
      recFile = clo->getValue( "r" );
    }
    else if ( clo->getValue( "rec" ) != 0 ) {
      recFile =  clo->getValue( "rec" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide a receptor name " << std::endl;
      return 0;
    }

    if ( clo->getValue( "c" ) != 0 ) {
      cutoff = strtod(clo->getValue( "cut" ), 0);
    }
    else if ( clo->getValue( "cut" ) != 0 ) {
      cutoff =  strtod(clo->getValue( "cut" ), 0);
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
      MTKpp::errorLogger.throwError("capActiveSite", errMessage, MTK_ERROR);
      oLog.close();
      exit(0);
    }

    errMessage = " Create a New Collection";
    MTKpp::errorLogger.throwError("capActiveSite", errMessage, INFO);

    collection* pCollection = new collection();
    molecule* pLigMol = 0;
    molecule* pRecMol = 0;

    errMessage = " Read Element Data";
    MTKpp::errorLogger.throwError("capActiveSite", errMessage, INFO);

    elementParser* pElementParser = new elementParser(pCollection->pElements);
    std::string elementXmlFile = AMBERHOME+"/dat/mtkpp/elements.xml";
    pElementParser->Read(elementXmlFile);

    errMessage = " Read MM parameters";
    MTKpp::errorLogger.throwError("capActiveSite", errMessage, INFO);

    pCollection->addParameters();
    paramParser* pParamParser = new paramParser(pCollection->getParameters());
    std::string parameterXmlFile = AMBERHOME+"/dat/mtkpp/parm94.xml";
    pParamParser->Read(parameterXmlFile);
    parameterXmlFile = AMBERHOME+"/dat/mtkpp/parm_gaff.xml";
    pParamParser->Read(parameterXmlFile);
    parameterXmlFile = AMBERHOME+"/dat/mtkpp/metals/metalParm.xml";
    pParamParser->Read(parameterXmlFile);

    errMessage = " Create a Standard Library";
    MTKpp::errorLogger.throwError("capActiveSite", errMessage, INFO);

    pCollection->addStdLibrary();
    stdLibrary* pStdLibrary = pCollection->getStdLibrary();
    if (!pStdLibrary) return 1;
    stdLibParser* pStdLibParser = new stdLibParser(pStdLibrary, pCollection->getParameters());
    std::string stdLibXmlFile = "";

    errMessage = " Reading N-Terminal Amino Acids Fragment Library";
    MTKpp::errorLogger.throwError("capActiveSite", errMessage, INFO);

    stdLibXmlFile = AMBERHOME+"/dat/mtkpp/aminont94.xml";
    pStdLibParser->Read(stdLibXmlFile);

    errMessage = " Reading Amino Acids Fragment Library";
    MTKpp::errorLogger.throwError("capActiveSite", errMessage, INFO);

    stdLibXmlFile = AMBERHOME+"/dat/mtkpp/amino94.xml";
    pStdLibParser->Read(stdLibXmlFile);

    errMessage = " Reading C-Terminal Amino Acids Fragment Library";
    MTKpp::errorLogger.throwError("capActiveSite", errMessage, INFO);

    stdLibXmlFile = AMBERHOME+"/dat/mtkpp/aminoct94.xml";
    pStdLibParser->Read(stdLibXmlFile);

    errMessage = " Reading Metals Fragment Library";
    MTKpp::errorLogger.throwError("capActiveSite", errMessage, INFO);

    stdLibXmlFile = AMBERHOME+"/dat/mtkpp/metals/metals.xml";
    pStdLibParser->Read(stdLibXmlFile);

    // CREATE PARSERS
    pdbParser* pPdbParser = 0;
    molParser* pMolParser = 0;
    pPdbParser = new pdbParser();
    pMolParser = new molParser();

    errMessage = " Reading receptor pdb file";
    MTKpp::errorLogger.throwError("capActiveSite", errMessage, INFO);

    pPdbParser->Read(recFile, pCollection);

    std::string ligExt = extName(ligFile);
    if (ligExt == "pdb") {
      errMessage = " Reading ligand pdb file";
      pPdbParser->Read(ligFile, pCollection);
    }
    else if (ligExt == "mol") {
      errMessage = " Reading ligand mol file";
      pMolParser->Read(ligFile, pCollection);
    }
    MTKpp::errorLogger.throwError("capActiveSite", errMessage, INFO);

    errMessage = " Assign Disulfide Bonds";
    MTKpp::errorLogger.throwError("capActiveSite", errMessage, INFO);

    connections* pConnections = new connections(pCollection);
    pConnections->assignDisulfideBonds();

    errMessage = " Atom Type";
    MTKpp::errorLogger.throwError("capActiveSite", errMessage, INFO);

    atomTyper* pAtomTyper = new atomTyper();
    pAtomTyper->atomTypeByLib(pCollection);

    errMessage = " Assign Connections (bonds, angles, torsions & impropers)";
    MTKpp::errorLogger.throwError("capActiveSite", errMessage, INFO);
    pConnections->run();

    errMessage = " Determine if its a metalloprotein ";
    MTKpp::errorLogger.throwError("capActiveSite", errMessage, INFO);
    bool bMetalloprotein = pCollection->hasMetal();
    if (bMetalloprotein) {
      pCollection->findMetals();
      pCollection->determineMetalEnvironments();
    }

    errMessage = " Get Ligand molecule and add Hydrogens ";
    MTKpp::errorLogger.throwError("capActiveSite", errMessage, INFO);

    pLigMol = pCollection->getMolecule();
    if (!pLigMol) {
      errMessage = " Cannot find ligand";
      MTKpp::errorLogger.throwError("capActiveSite", errMessage, MTK_ERROR);
      return 0;
    }

    if (ligExt == "mol") {
      pLigMol->addHydrogens();
    }
    std::vector<atom*> ligAtoms = pLigMol->getAtomList();

    pRecMol = pCollection->getMolecule(1);
    if (!pRecMol) {
      errMessage = " Cannot find receptor";
      MTKpp::errorLogger.throwError("capActiveSite", errMessage, MTK_ERROR);
      return 0;
    }

    errMessage = " Get all residues not in the ligand ";
    MTKpp::errorLogger.throwError("capActiveSite", errMessage, INFO);

    std::vector<molecule*> molList = pCollection->getMoleculeList();
    std::vector<submolecule*> recSmolList;
    for (unsigned int x = 0; x < molList.size(); x++) {
      if (molList[x] == pLigMol) continue;
      std::vector<submolecule*> smolList = molList[x]->getSubMoleculeList();
      for (unsigned int x2 = 0; x2 < smolList.size(); x2++) {
        recSmolList.push_back(smolList[x2]);
      }
    }

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

    for (unsigned int i = 0; i < recSmolList.size(); i++) {
      submolecule* pSubMol = recSmolList[i];
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

    // 
    if (bMetalloprotein) {
      // Determine if metal is included in list,
      // if so then add the ligating residues
      std::vector<metalCenter*> metalCenters = pCollection->getMetalCenters();
      std::vector<atom*> primaryShellAtoms;

      for (unsigned int i = 0; i < aSiteResidues.size(); i++) {
        for (unsigned int j = 0; j < metalCenters.size(); j++) {
          atom* pMetal = metalCenters[j]->getMetalAtom();
          if (aSiteResidues[i] == pMetal->getParent()->getSubMolId()) {
            metalCenters[j]->getPrimaryShellAtoms(primaryShellAtoms);
            break;
          }
        }
      }

      for (unsigned int i = 0; i < primaryShellAtoms.size(); i++) {
        std::vector<int>::iterator result;
        submolecule* pRes = primaryShellAtoms[i]->getParent();

        if (pRes->getParent() == pLigMol) continue;

        result = std::find(aSiteResidues.begin(), aSiteResidues.end(), pRes->getSubMolId());

        if (result == aSiteResidues.end()) {
          activeSiteResiduesSMol.push_back(pRes);
          aSiteResidues.push_back(pRes->getSubMolId());
          aSiteMolecules.push_back(pRes->getParent()->getMolId());
        }
      }
    }

    if (aSiteResidues.size() < 1) {
      errMessage = " No active site found";
      MTKpp::errorLogger.throwError("capActiveSite", errMessage, MTK_ERROR);
      return 1;
    }

/*
    oLog << "\n RESIDUE INDICES:" << std::endl;
    std::sort(activeSiteResidues.begin(), activeSiteResidues.end());
    for (unsigned int i = 0; i < activeSiteResidues.size(); i++) {
      oLog << activeSiteResidues[i] << " ";
    }
*/
    oLog << "\n RESIDUE INDICES IN ACTIVE SITE:" << std::endl;
    std::sort(activeSiteResiduesSMol.begin(), activeSiteResiduesSMol.end(), submolecule::less);
    for (unsigned int i = 0; i < activeSiteResiduesSMol.size(); i++) {
      activeSiteResidues.push_back(activeSiteResiduesSMol[i]->getSubMolId());
      activeSiteMolecules.push_back(activeSiteResiduesSMol[i]->getParent()->getMolId());
      oLog << activeSiteResiduesSMol[i]->getSubMolId() << " ";
    }
    //oLog.close();
    //exit(1);

    // Create segments and fill in gaps
    int prevRes = activeSiteResidues[0];
    int prevMol = activeSiteMolecules[0];
    std::vector<int> fullList;
    std::vector<int> fullListMols;

    int gapFill = 3;
    submolecule* pSubMol = 0;
    for (unsigned int i = 0; i < activeSiteResidues.size(); i++) {
      if (activeSiteMolecules[i] == prevMol) {
        pRecMol = pCollection->getMolecule(activeSiteMolecules[i]);
        int diff = activeSiteResidues[i] - prevRes;
        if ((diff <= gapFill) and (diff > 0)) {
          for (int j = prevRes+1; j < activeSiteResidues[i]; j++) {
            pSubMol = pRecMol->getSubMolecule(j);
            if (pSubMol) {
              if ((pSubMol->getName() != "WAT") and (pSubMol->getName() != "HOH")) {
                fullList.push_back(j);
                fullListMols.push_back(activeSiteMolecules[i]);
              }
            }
            else {
              std::cout << " capActiveSite::Warning, Can't find residue: " << j << std::endl;
            }
          }
        }
      }
      fullList.push_back(activeSiteResidues[i]);
      fullListMols.push_back(activeSiteMolecules[i]);
      prevRes = activeSiteResidues[i];
      prevMol = activeSiteMolecules[i];
    }

    //std::sort(fullList.begin(), fullList.end());

    oLog << "\n RESIDUE INDICES IN ACTIVE SITE AFTER GAP FILLING:" << std::endl;
    for (unsigned int i = 0; i < fullList.size(); i++) {
      oLog << fullList[i] << " ";
    }

    oLog << "\n\n  CORRESPONDING MOLECULE INDICES:" << std::endl;
    for (unsigned int i = 0; i < fullListMols.size(); i++) {
      oLog << fullListMols[i] << " ";
    }
    oLog << " " << std::endl;

    // Get segment to cap
    std::vector<std::vector<int> > segments;
    prevRes = fullList[0];
    std::vector<int> seg;
    seg.push_back(fullList[0]);
    for (unsigned int i = 1; i < fullList.size(); i++) {
      if (fullList[i] - prevRes == 1) {
        seg.push_back(fullList[i]);
      }
      else {
        segments.push_back(seg);
        seg.clear();
        seg.push_back(fullList[i]);
      }
      prevRes = fullList[i];
    }
    segments.push_back(seg);

    oLog << "\n SEGMENTS INDICES:" << std::endl;
    for (unsigned int i = 0; i < segments.size(); i++) {
      for (unsigned int j = 0; j < segments[i].size(); j++) {
        oLog << segments[i][j] << " |";
      }
      oLog << " " << std::endl;
    }
    oLog << " " << std::endl;

    atom* pAtom2 = 0;
    submolecule* pSubMol2 = 0;

    int molIndex = 0;
    std::vector<molecule*> addedMolecules;
    for (unsigned int i = 0; i < segments.size(); i++) {
      molecule* molSeg = pCollection->addMolecule();
      addedMolecules.push_back(molSeg);
      for (unsigned int j = 0; j < segments[i].size(); j++) {
        pRecMol = pCollection->getMolecule(fullListMols[molIndex]);
        molIndex++;
        if (pRecMol->getNumSubMolecules() == 1) {
          pSubMol = pRecMol->getSubMolecule(segments[i][j]);
          if (!pSubMol) continue;
          pSubMol2 = molSeg->addSubMolecule();
          //oLog << " Adding " << pSubMol->getName() << std::endl;
          pSubMol2->copy(pSubMol);
          continue;
        }

        // Cap First amino acid in segment with ACE
        if ((j == 0) and (segments[i][j]-1 != 0)) {
          pSubMol = pRecMol->getSubMolecule(segments[i][j]-1);
          if (!pSubMol) continue;
          if ((pSubMol->getName() != "WAT") and (pSubMol->getName() != "HOH")) {
            std::vector<atom*> atomList = pSubMol->getAtomList();
            pSubMol2 = molSeg->addSubMolecule();
            pSubMol2->setName("ACE");
            pSubMol2->setSubMolId(segments[i][j]-1);
            for (unsigned int k = 0; k < atomList.size(); k++) {
              std::string atomName = atomList[k]->getName();
              if (atomName == " C  ") {
                pAtom2 = pSubMol2->addAtom();
                pAtom2->setName(" C  ");
                pAtom2->setElement(pCollection->pElements->getElement("C"));
                pAtom2->setCoords(atomList[k]->getX(),atomList[k]->getY(),atomList[k]->getZ());
                pAtom2->setFileID(atomList[k]->getFileID());
              }
              else if (atomName == " O  ") {
                pAtom2 = pSubMol2->addAtom();
                pAtom2->setName(" O  ");
                pAtom2->setElement(pCollection->pElements->getElement("O"));
                pAtom2->setCoords(atomList[k]->getX(),atomList[k]->getY(),atomList[k]->getZ());
                pAtom2->setFileID(atomList[k]->getFileID());
              }
              else if (atomName == " CA ") {
                pAtom2 = pSubMol2->addAtom();
                pAtom2->setName(" CH3");
                pAtom2->setElement(pCollection->pElements->getElement("C"));
                pAtom2->setCoords(atomList[k]->getX(),atomList[k]->getY(),atomList[k]->getZ());
                pAtom2->setFileID(atomList[k]->getFileID());
              }
            }
          }
        }

        pSubMol = pRecMol->getSubMolecule(segments[i][j]);
        if (!pSubMol) continue;
        pSubMol2 = molSeg->addSubMolecule();
        pSubMol2->copy(pSubMol);

        // Cap Last amino acid in segment with NME
        if ((j == segments[i].size()-1) and (segments[i][j] != recSmolList[recSmolList.size()-1]->getSubMolId())) {
          pSubMol = pRecMol->getSubMolecule(segments[i][j]+1);
          if (!pSubMol) continue;
          if ((pSubMol->getName() != "WAT") and (pSubMol->getName() != "HOH")) {
            std::vector<atom*> atomList = pSubMol->getAtomList();
            pSubMol2 = molSeg->addSubMolecule();
            pSubMol2->setName("NME");
            pSubMol2->setSubMolId(segments[i][j]+1);
            for (unsigned int k = 0; k < atomList.size(); k++) {
              std::string atomName = atomList[k]->getName();
              if (atomName== " N  ") {
                pAtom2 = pSubMol2->addAtom();
                pAtom2->setName(" N  ");
                pAtom2->setElement(pCollection->pElements->getElement("N"));
                pAtom2->setCoords(atomList[k]->getX(),atomList[k]->getY(),atomList[k]->getZ());
                pAtom2->setFileID(atomList[k]->getFileID());
              }
              else if (atomName== " H  ") {
                pAtom2 = pSubMol2->addAtom();
                pAtom2->setName(" H  ");
                pAtom2->setElement(pCollection->pElements->getElement("H"));
                pAtom2->setCoords(atomList[k]->getX(),atomList[k]->getY(),atomList[k]->getZ());
                pAtom2->setFileID(atomList[k]->getFileID());
              }
              else if (atomName== " CA ") {
                pAtom2 = pSubMol2->addAtom();
                pAtom2->setName(" CH3");
                pAtom2->setElement(pCollection->pElements->getElement("C"));
                pAtom2->setCoords(atomList[k]->getX(),atomList[k]->getY(),atomList[k]->getZ());
                pAtom2->setFileID(atomList[k]->getFileID());
              }
            }
          }
        }
      }
    }

    std::vector<molecule*> myMols = pCollection->getMoleculeList();
    std::vector<molecule*> molsToBeDeleted;
    for (unsigned int k = 0; k < myMols.size(); k++) {
      bool d = true;
      for (unsigned int k2 = 0; k2 < addedMolecules.size(); k2++) {
        if (addedMolecules[k2] == myMols[k]) {
          d = false;
          break;
        }
      }
      if (d) {
        molsToBeDeleted.push_back(myMols[k]);
      }
    }
    for (int k = molsToBeDeleted.size()-1; k > -1; k--) {
      pCollection->delMolecule(molsToBeDeleted[k]);
    }

    pPdbParser->Write(outFile, pCollection);

    // Clean up
    oLog.close();
    delete pPdbParser;
    delete pMolParser;
    return 0;
}

/*
    ADD THIS CODE WHEN I GET SOME TIME
	# Check to see if segments are made up of amino acids
	aminoAcids = 0
	doNotCap = []
	for s in segments:
		for a in s:
			smol = col.getSubMolecule(resNum = a)
			if smol.stdFrag is not None:
				try:
					L = amino_acids.LLL2L[smol.stdFrag.symbol]
					aminoAcids+=1
				except:
					pass
		if len(s) != aminoAcids:
			doNotCap.append(segments.index(s))
		aminoAcids = 0

	# Create a non-library residue list
	doNotCapSeg = []
	for s in doNotCap:
		doNotCapSeg.append(segments[s])

	doNotCapAA = []
	for s in doNotCapSeg:
		for r in s:
			doNotCapAA.append(r)

	# Remove non-library residues
	for s in range(len(doNotCap)-1, -1, -1):
		del segments[doNotCap[s]]
*/
