/*!
   \file superimposer.cpp

   \brief Superimpose functions

   \author Martin B. Peters

   $Date: 2010/04/22 22:19:54 $
   $Revision: 1.8 $

   ----------------------------------------------------------------------------

   MTK++ - C++ package of modeling libraries.

   Copyright (C) 2005-2007  (see AUTHORS file for a list of contributors)

   ----------------------------------------------------------------------------
*/
// - Utils
#include "Utils/printHeader.h"
#include "Utils/vector3d.h"
#include "Utils/constants.h"

// - MOLECULE
#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/torsion.h"
#include "Molecule/element.h"
#include "Molecule/superimpose.h"
#include "Molecule/parameters.h"
#include "Molecule/stdLibrary.h"
#include "Molecule/atomTyper.h"
#include "Molecule/connections.h"

// - PARSERS
#include "Parsers/elementParser.h"
#include "Parsers/molParser.h"
#include "Parsers/sdfParser.h"
#include "Parsers/pdbParser.h"
#include "Parsers/StringManip.h"
#include "Parsers/paramParser.h"
#include "Parsers/stdLibParser.h"
#include "Parsers/inputParser.h"

// - TIME
#include "time.h"

// - COMMAND LINE OPTIONS
#include "Parsers/commLineOptions.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <map>
#include <algorithm>

using namespace MTKpp;

/*!
   \brief Superimpose functions
   \param argc
   \param argv
   \return success
*/
int main (int argc, char **argv)
{

    std::string prog_name = "superimposer";
    std::vector<std::string> authors;
    std::string author = "Martin B. Peters";
    authors.push_back(author);

    // 1. CREATE AN OBJECT
    commLineOptions *clo = new commLineOptions();

    // 2. SET PREFERENCES
    clo->noUsage();

    // 3. SET THE USAGE/HELP
    clo->addUsage( "  superimposer: Structural Superimposition            \n" );
    clo->addUsage( "    usage:  superimposer [flags] [options]            \n" );
    clo->addUsage( "  options:                                             " );
    clo->addUsage( "          -t template file                             " );
    clo->addUsage( "          -l template lib xml file                     " );
    clo->addUsage( "          -f flexible file                             " );
    clo->addUsage( "          -e flexible lib xml file                     " );
    clo->addUsage( "          -p aligned file                              " );
    clo->addUsage( "          -k calculation kind                          " );
    clo->addUsage( "            :- 1 coordinate RMSD [no fitting]          " );
    clo->addUsage( "            :- 2 atom type based RMSD [no fitting] default" );
    clo->addUsage( "            :- 3 RMSD                                  " );
    clo->addUsage( "            :- 4 atom type based RMSD                  " );
    clo->addUsage( "            :- 5 RMSD & Molecules Superimposed         " );
    clo->addUsage( "          -v verbose                                   " );
    clo->addUsage( "            :- 0 rmsd is printed [default]             " );
    clo->addUsage( "            :- 1 all output                            " );
    clo->addUsage( "          -a log file                                  " );
    clo->addUsage( "          -o output file                             \n" );
    clo->addUsage( "    flags:                                             " );
    clo->addUsage( "          -h help                                    \n" );

    // 4. SET THE OPTION STRINGS/CHARACTERS
    clo->setOption( "tempFile",   't' );
    clo->setOption( "libTempXml", 'l' );
    clo->setOption( "flexFile",   'f' );
    clo->setOption( "libFlexXml", 'e' );
    clo->setOption( "alignedPdb", 'p' );
    clo->setOption( "kind",       'k' );
    clo->setOption( "verbose",    'v' );
    clo->setOption( "log",        'a' );
    clo->setOption( "output",     'o' );
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

    std::string templateFile = "";
    std::string tempLibFile  = "";
    std::string flexibleFile = "";
    std::string flexLibFile  = "";
    std::string alignedFile  = "";
    int kind                 = 2;
    int verbose              = 0;
    std::string outputFile   = "superimposer.out";
    std::string logFile      = "superimposer.log";

    if ( clo->getValue( "t" ) != 0 ) {
      templateFile = clo->getValue( "t" );
    }
    else if ( clo->getValue( "tempFile" ) != 0 ) {
      templateFile =  clo->getValue( "tempFile" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide a template pdb file " << std::endl;
      return 0;
    }

    if ( clo->getValue( "l" ) != 0 ) {
      tempLibFile = clo->getValue( "l" );
    }
    else if ( clo->getValue( "libTempXml" ) != 0 ) {
      tempLibFile =  clo->getValue( "libTempXml" );
    }

    if ( clo->getValue( "f" ) != 0 ) {
      flexibleFile = clo->getValue( "f" );
    }
    else if ( clo->getValue( "flexFile" ) != 0 ) {
      flexibleFile =  clo->getValue( "flexFile" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide a flexible pdb file " << std::endl;
      return 0;
    }

    if ( clo->getValue( "e" ) != 0 ) {
      flexLibFile = clo->getValue( "e" );
    }
    else if ( clo->getValue( "libFlexXml" ) != 0 ) {
      flexLibFile =  clo->getValue( "libFlexXml" );
    }

    if ( clo->getValue( "p" ) != 0 ) {
      alignedFile = clo->getValue( "p" );
    }
    else if ( clo->getValue( "alignedPdb" ) != 0 ) {
      alignedFile =  clo->getValue( "alignedPdb" );
    }

    if ( clo->getValue( "k" ) != 0 ) {
      kind = atoi(clo->getValue( "k" ));
    }
    else if ( clo->getValue( "kind" ) != 0 ) {
      kind =  atoi(clo->getValue( "kind" ));
    }

    if ( clo->getValue( "v" ) != 0 ) {
      verbose = atoi(clo->getValue( "v" ));
    }
    else if ( clo->getValue( "verbose" ) != 0 ) {
      verbose =  atoi(clo->getValue( "verbose" ));
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

    if ( clo->getValue( "o" ) != 0 ) {
      outputFile = clo->getValue( "o" );
    }
    else if ( clo->getValue( "output" ) != 0 ) {
      outputFile =  clo->getValue( "output" );
    }
    else {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      std::cout << " Please provide an output file " << std::endl;
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

    // Start & end time
    time_t startTime;
    time_t endTime;

    // Get start time
    time (&startTime);

    std::string AMBERHOME = getenv("AMBERHOME");
    if (AMBERHOME == "") {
      errMessage = " Can't find AMBERHOME ";
      MTKpp::errorLogger.throwError("superimposer", errMessage, MTK_ERROR);
      oLog.close();
      exit(0);
    }

    // Create collection
    collection* pCollection = new collection();

    // Read element data
    elementParser* pElementParser = new elementParser(pCollection->pElements);

    //std::string elementXmlFile = MTKppDIR+"data/elements.xml";
    std::string elementXmlFile = AMBERHOME+"/dat/mtkpp/elements.xml";

    pElementParser->Read(elementXmlFile);
    delete pElementParser;

    std::string tempExt = extName(templateFile);
    std::string flexExt = extName(flexibleFile);

    if (tempExt == "pdb" and (flexExt == "pdb" or flexExt == "log")) {
      if (verbose) {
        errMessage = " Read parameter files ";
        MTKpp::errorLogger.throwError("superimposer", errMessage, 4);
      }
      pCollection->addParameters();
      paramParser* pParamParser = new paramParser(pCollection->getParameters());

      //std::string parameterXmlFile = MTKppDIR+"data/parm94.xml";
      std::string parameterXmlFile = AMBERHOME+"/dat/mtkpp/parm94.xml";

      pParamParser->Read(parameterXmlFile);
      //parameterXmlFile = MTKppDIR+"data/parm_gaff.xml";
      parameterXmlFile = AMBERHOME+"/dat/mtkpp/parm_gaff.xml";

      pParamParser->Read(parameterXmlFile);

      if (verbose) {
        errMessage = " Read library files ";
        MTKpp::errorLogger.throwError("superimposer", errMessage, 4);
      }

      pCollection->addStdLibrary();
      stdLibParser* pStdLibParser = new stdLibParser(
                        pCollection->getStdLibrary(),
                        pCollection->getParameters());

      if (tempLibFile != "") {
        pStdLibParser->Read(tempLibFile);
      }
      if (flexLibFile != "") {
        pStdLibParser->Read(flexLibFile);
      }
    }

    connections* pConnections = new connections(pCollection);
    molecule* pMoleculeA = 0;
    int nAtomsA = 0;

    // - Template molecule - //
    std::vector<std::string> tempNames;
    std::vector<std::string> minConformerNames;
    std::vector<double> minConformerRMSDs;

    if (verbose) {
      errMessage = " Read template file ";
      MTKpp::errorLogger.throwError("superimposer", errMessage, 4);
    }

    if (tempExt == "mol") {
      molParser* pMolParser = 0;
      pMolParser = new molParser();
      pMolParser->Read(templateFile, pCollection);
      delete pMolParser;
    }
    else if (tempExt == "pdb") {
      pdbParser* pPdbParser = 0;
      pPdbParser = new pdbParser();
      pPdbParser->Read(templateFile, pCollection);
      delete pPdbParser;
      atomTyper* pAtomTyper = new atomTyper();
      pAtomTyper->atomTypeByLib(pCollection->getMolecule(1));
      delete pAtomTyper;
      pMoleculeA = pCollection->getMolecule(1);

      pConnections->assignBonds(pMoleculeA);
      pConnections->assignAngles(pMoleculeA);
      pConnections->assignStdBondsAngles(pMoleculeA);
      pConnections->assignTorsions(pMoleculeA);
    }
    else if (tempExt == "sdf") {
      sdfParser* pSdfParser = 0;
      pSdfParser = new sdfParser();

      pSdfParser->Read(templateFile, pCollection);
      delete pSdfParser;

      if (pCollection->getNumberMolecules() < 1) {
        errMessage = " No molecules found in template file ";
        MTKpp::errorLogger.throwError("superimposer", errMessage, MTK_ERROR);
        oLog.close();
        return 0;
      }

      std::vector<molecule*> molList = pCollection->getMoleculeList();
      for (unsigned int m = 0; m < molList.size(); m++) {
        std::string tempName = molList[m]->getName();
        std::string tempNameStrip = stripString(tempName, " ");

        //std::cout << " Template name = " << tempNameStrip << std::endl;
        tempNames.push_back(tempNameStrip);

        minConformerNames.push_back("");
        minConformerRMSDs.push_back(BIGNUM);

        // Only doing heavy atom
        molList[m]->removeHydrogens();
      }
    }

    if (tempExt != "sdf") {
      pMoleculeA = pCollection->getMolecule(1);
      nAtomsA = pMoleculeA->getNumAtoms();
    }

    // Open output file
    if (verbose) {
      MTKpp::errorLogger.throwError("superimposer", " Open output file ", 4);
    }

    std::ofstream oOutput;
    oOutput.open(outputFile.c_str());
    if (!oOutput) {
      errMessage = " Unable to open output file ";
      MTKpp::errorLogger.throwError("superimposer", errMessage, MTK_ERROR);
      oLog.close();
      return 0;
    }

    // - Flexible molecule - //
    molecule* pMoleculeB = 0;

    if (verbose) {
      errMessage = " Read flexible files ";
      MTKpp::errorLogger.throwError("superimposer", errMessage, 4);
    }

    if (flexExt == "mol") {
      molParser* pMolParser = 0;
      pMolParser = new molParser();
      pMolParser->Read(flexibleFile,pCollection);
      delete pMolParser;
    }
    else if (flexExt == "pdb") {
      pdbParser* pPdbParser = 0;
      pPdbParser = new pdbParser();
      pPdbParser->Read(flexibleFile, pCollection);
      delete pPdbParser;
      atomTyper* pAtomTyper = new atomTyper();
      pAtomTyper->atomTypeByLib(pCollection->getMolecule(2));
      delete pAtomTyper;

////////////////

      pMoleculeB = pCollection->getMolecule(2);
      pConnections->assignBonds(pMoleculeB);
      pConnections->assignAngles(pMoleculeB);
      pConnections->assignStdBondsAngles(pMoleculeB);
      pConnections->assignTorsions(pMoleculeB);
/*
      std::vector<std::vector<int> > correspondenceMatrices;
      superimpose* pSuperImpose = new superimpose();

    int nAtomsB = pMoleculeB->getNumAtoms();

          int f = pSuperImpose->initializeCorrespondences(pMoleculeA, pMoleculeB, 2, correspondenceMatrices);
          if (f) {
            errMessage = " Unable to initialize correspondences ";
            MTKpp::errorLogger.throwError("superimposer", errMessage, MTK_ERROR);
            oLog.close();
            oOutput.close();
            delete pSuperImpose;
            exit(0);
          }
          std::cout << " Number of correspondences = " << correspondenceMatrices.size() << std::endl;

exit(1);
*/
//////////////////
    }
    else if (flexExt == "log") {
      std::vector<std::vector<std::string> > fileContents;
      int f = readInputFile(flexibleFile, fileContents);
      if (f) {
        errMessage = " Unable to read .log file ";
        MTKpp::errorLogger.throwError("superimposer", errMessage, MTK_ERROR);
        oLog.close();
        oOutput.close();
        exit(0);
      }
      std::vector<std::vector<int> > correspondenceMatrices;
      superimpose* pSuperImpose = new superimpose();

      double dRMSD = 0.0;

      pdbParser* pPdbParser = 0;
      pPdbParser = new pdbParser();
      int nAtomsB = 0;

      for (unsigned int p = 0; p < fileContents.size(); p++) {
        std::string curFile = fileContents[p][0];
        std::vector< vector3d > CoordsB;

        if (p == 0) {
          pPdbParser->Read(curFile+".pdb", pCollection);
          pMoleculeB = pCollection->getLastAddedMolecule();

          atomTyper* pAtomTyper = new atomTyper();
          pAtomTyper->atomTypeByLib(pMoleculeB);
          delete pAtomTyper;

          pConnections->assignBonds(pMoleculeB);
          pConnections->assignAngles(pMoleculeB);
          pConnections->assignStdBondsAngles(pMoleculeB);
          pConnections->assignTorsions(pMoleculeB);

          nAtomsB = pMoleculeB->getNumAtoms();
          int f = pSuperImpose->initializeCorrespondences(pMoleculeA, pMoleculeB, 2, correspondenceMatrices);
          if (f) {
            errMessage = " Unable to initialize correspondences ";
            MTKpp::errorLogger.throwError("superimposer", errMessage, MTK_ERROR);
            oLog.close();
            oOutput.close();
            delete pSuperImpose;
            exit(0);
          }
          std::cout << " Number of correspondences = " << correspondenceMatrices.size() << std::endl;
        }
        else {
          pPdbParser->updateMolCoords(curFile+".pdb", pMoleculeB);
        }

        // Calculate the RMSD based on atom name
        int nFittedAtoms = 0;
        int type = 2;
        dRMSD = pSuperImpose->minRMSD(pMoleculeA, pMoleculeB, type, nFittedAtoms);
        std::cout << " minRMSD = " << dRMSD << std::endl;

        if (nFittedAtoms != nAtomsA) {
          errMessage = " Unable to determine min rmsd ";
          MTKpp::errorLogger.throwError("superimposer", errMessage, MTK_ERROR);
          oLog.close();
          oOutput.close();
          delete pSuperImpose;
          delete pCollection;
          exit(0);
        }

        pMoleculeB->getCoordinates(CoordsB);
        dRMSD = pSuperImpose->rmsdNoFit(pMoleculeA, CoordsB, correspondenceMatrices);
        std::cout << " rmsdNoFit = " << dRMSD << std::endl;

        for (unsigned int d = 0; d < fileContents[p].size(); d++) {
          oOutput << fileContents[p][d] << " ";
        }
        oOutput << dRMSD << std::endl;
      }
      delete pConnections;
      delete pPdbParser;
      delete pSuperImpose;
      oOutput.close();
      oLog.close();
      exit(0);
    }
    else if (flexExt == "sdf") {
      std::vector<molecule*> tempMolList = pCollection->getMoleculeList();

      // Operate on each molecule in the flexible sd file, one at a time
      std::ifstream isdf;
      isdf.open(flexibleFile.c_str());

      if (!isdf) {
        errMessage = " Unable to open flexible sdf file ";
        MTKpp::errorLogger.throwError("superimposer", errMessage, MTK_ERROR);
        oLog.close();
        oOutput.close();
        exit(1);
      }

      // Open output sd file
      std::string out_sdf = baseName(flexibleFile) + ".out.sdf";
      std::ofstream osdf;
      osdf.open(out_sdf.c_str());

      if (!osdf) {
        std::cout << "\nUNABLE TO OPEN SDF FILe"
                  << "\nFILENAME = " << out_sdf << std::endl;
        oLog.close();
        oOutput.close();
        exit(1);
      }

      int n = 1;
      std::string flexMolName = "";
      int failure = 0;
      sdfParser* pSdf = new sdfParser();
/*
      int nFlexibleMolecules = pSdf->numMolecules(flexibleFile);
      std::cout << " Number of Molecules = " << nFlexibleMolecules << "\n";
*/
      while (!isdf.eof()) {
        pMoleculeB = pCollection->addMolecule();

        failure = pSdf->ReadMolecule(isdf, pMoleculeB);
        //std::cout << " read sdf molecule " << failure << std::endl;

        pMoleculeB->removeHydrogens();

        if (failure == 1) {
          //std::cout << " ERROR " << std::endl;
          continue;
        }
        else if (failure) {
          MTKpp::errorLogger.throwError("superimposer", " Unable to read sd file", MTK_ERROR);
          oLog.close();
          oOutput.close();
          delete pSdf;
          exit(1);
        }
        if (pMoleculeB->getNumAtoms() == 0) {
          std::cout << " ERROR2 " << std::endl;
          break;
        }

        flexMolName = pMoleculeB->getName();
        if (flexMolName.size() == 0) {
          MTKpp::errorLogger.throwError("superimposer", " No name present for compound", MTK_ERROR);
          oLog.close();
          oOutput.close();
          delete pSdf;
          exit(1);
        }

        // Split flexMolName based on '_'
        std::vector<std::string> nameParts;
        splitString(flexMolName, "_", nameParts, 0); // was '_'
        if (nameParts.size() < 2) {
          MTKpp::errorLogger.throwError("nameParts", "No _ present in compound name", MTK_ERROR);
          oLog.close();
          oOutput.close();
          delete pSdf;
          exit(1);
        }

        double dRMSD = 0.0;
        for (unsigned int m = 0; m < tempNames.size(); m++) {

          bool bFoundName = false;
          for (unsigned int p = 0; p < nameParts.size(); p++) {
            if (nameParts[p] == tempNames[m]) {
              bFoundName = true;
            }
/*
            else {
              std::cout << "|" << nameParts[p] <<  "| != |" << tempNames[m] << "|" << std::endl;
            }
*/
          }

          //if (nameParts[0] == tempNames[m] or nameParts[nameParts.size()-1] == tempNames[m]) {
          if (bFoundName) {
            pMoleculeA = tempMolList[m];

            std::vector<std::vector<int> > correspondenceMatrices;
            superimpose* pSuperImpose = new superimpose();

            int f = pSuperImpose->initializeCorrespondences(pMoleculeA, pMoleculeB, 3, correspondenceMatrices);
            if (f) {
              errMessage = " Unable to initialize correspondences ";
              MTKpp::errorLogger.throwError("superimposer", errMessage, 1);
              oLog.close();
              oOutput.close();
              delete pSuperImpose;
              exit(0);
            }

            std::vector< vector3d > CoordsB;
            pMoleculeB->getCoordinates(CoordsB);
            int cor = 0;

            dRMSD = pSuperImpose->rmsd(pMoleculeA, CoordsB, correspondenceMatrices, cor);
            oOutput << tempNames[m] << " " << flexMolName << " " << dRMSD << " " << cor << std::endl;

            if (dRMSD < minConformerRMSDs[m]) {
              minConformerNames[m] = flexMolName;
              minConformerRMSDs[m] = dRMSD;

              superimpose* pSuperImpose2 = new superimpose();
              f = pSuperImpose2->fit(pMoleculeA, pMoleculeB,
                                    correspondenceMatrices[cor], 3);
              delete pSuperImpose2;
              // Write file
              molParser* pMolParser = 0;
              pMolParser = new molParser();
              pMolParser->Write(tempNames[m]+"_conf.mol", pMoleculeB);
              delete pMolParser;
            }

            // Write aligned molecule and rmsd value
            pMoleculeB->addProperty("RMSD", double2String(dRMSD));
            pSdf->Write(osdf, pMoleculeB);

            delete pSuperImpose;
          }
          else {
            //std::cout << tempNames[0] << " " << nameParts[0] << std::endl;
            //std::cout << " no name match template = |" << tempNames[m] << "| flex mol = |" << nameParts[0] << "| " << std::endl;
          }
        }

        pCollection->delMolecule(pMoleculeB);
        pMoleculeB = 0;
        n++;
      }

      errMessage = " Lowest RMSD Conformers: \n ";
      errMessage += " Template | Conformer | RMSD \n";
      for (unsigned int m = 0; m < tempNames.size(); m++) {
        errMessage += "  " + tempNames[m] + " | " + minConformerNames[m] + " | " + double2String(minConformerRMSDs[m]) + " \n";
      }
      MTKpp::errorLogger.throwError("superimposer", errMessage, MESSAGE);

      delete pSdf;

      time (&endTime);
      int diffTime = (int) difftime(endTime, startTime);
      errMessage = " Exiting Normally After "
                    + int2String(diffTime) + " Seconds ";
      MTKpp::errorLogger.throwError("superimposer", errMessage, MESSAGE);

      oOutput.close();
      oLog.close();
      exit(0);
    }

    pMoleculeB = pCollection->getMolecule(2);
    int nAtomsB = pMoleculeB->getNumAtoms();

    if (nAtomsA != nAtomsB) {
      std::cout << "  superimposer::Error - Different number of atoms in molecules ... exiting " << std::endl;
      return 1;
    }
    else {
      double rmsd = 0.0;
      std::string itsType = "";
      int type = 1;
      if (tempExt == "pdb" and flexExt == "pdb") {
        type = 2;
      }

      superimpose* pSuperImpose = new superimpose(nAtomsA);
      double coordsA[nAtomsA][3];
      double coordsB[nAtomsB][3];
      pMoleculeA->getCoordinates(coordsA);
      pMoleculeB->getCoordinates(coordsB);

      if (kind == 1) {
        if (verbose) {
          puts("  superimposer::Coordinate RMSD [No Fitting]");
        }
        rmsd = pSuperImpose->calculateRMSD(coordsA, coordsB);
        itsType = "Coordinate RMSD [No Fitting]";
      }
      else if (kind == 2) {
        if (verbose) {
          puts("  superimposer::Atom Type Based RMSD [no fitting] (default)");
        }
        rmsd = 0.0;
        int nFittedAtoms = 0;
        rmsd = pSuperImpose->minRMSD(pMoleculeA, pMoleculeB, type, nFittedAtoms);
        itsType = "Atom Type Based RMSD [no fitting]";
        if (nFittedAtoms != nAtomsA) {
          puts("  superimposer::Error in pSuperImpose->minRMSD ... exiting ");
          delete pSuperImpose;
          delete pCollection;
          exit(0);
        }
      }
      else if (kind == 3) {
        puts("  superimposer::RMSD");
        rmsd = pSuperImpose->rmsd(pMoleculeA, pMoleculeB);
        itsType = "RMSD";
      }
      else if (kind == 4) {
        puts("  superimposer::Atom Type Based RMSD & Molecules Superimposed");
        std::vector<std::vector<int> > correspondenceMatrices;
        int f = pSuperImpose->initializeCorrespondences(pMoleculeA, pMoleculeB, type, correspondenceMatrices);

        if (f) {
          std::cout << " Error in pSuperImposer::initializeCorrespondences " << std::endl;
          delete pSuperImpose;
          delete pCollection;
          exit(0);
        }

        //std::cout << " Number of correspondences = " << correspondenceMatrices.size() << std::endl;

        std::vector< vector3d > tmp_molBCoords;
        pMoleculeB->getCoordinates(tmp_molBCoords);
        int cor = 0;
        rmsd = pSuperImpose->rmsd(pMoleculeA, tmp_molBCoords, correspondenceMatrices, cor);
        itsType = "Atom Type Based RMSD";

        superimpose* pSuperImpose2 = new superimpose();
        f = pSuperImpose2->fit(pMoleculeA, pMoleculeB,
                               correspondenceMatrices[cor], 3);
        delete pSuperImpose2;

        std::string outExt = extName(alignedFile);
        if (outExt == "mol") {
          molParser* pMolParser = 0;
          pMolParser = new molParser();
          pMolParser->Write(alignedFile, pMoleculeB);
          delete pMolParser;
        }
        else if (outExt == "pdb") {
          pdbParser* pPdbParser = 0;
          pPdbParser = new pdbParser();
          pPdbParser->Write(alignedFile, pMoleculeB);
          delete pPdbParser;
        }
      }
      else if (kind == 5) {
        puts("  superimposer::RMSD & Molecules Superimposed");
        rmsd = pSuperImpose->fit(pMoleculeA, pMoleculeB);
        itsType = "RMSD & Molecules Superimposed";
        std::string outExt = extName(alignedFile);
        if (outExt == "mol") {
          molParser* pMolParser = 0;
          pMolParser = new molParser();
          pMolParser->Write(alignedFile, pMoleculeB);
          delete pMolParser;
        }
        else if (outExt == "pdb") {
          pdbParser* pPdbParser = 0;
          pPdbParser = new pdbParser();
          pPdbParser->Write(alignedFile, pMoleculeB);
          delete pPdbParser;
        }
      }
      else {
        std::cout << "  superimposer::Unknown kind " << kind << " ... exiting " << std::endl;
        return 1;
      }

      oOutput << "RMSD = " << rmsd << std::endl;
      oOutput << "TYPE = |" << itsType << "| " << std::endl;
      delete pSuperImpose;
      oOutput.close();
    }

    // - Clean up - //
    if (verbose) puts("  superimposer::Clean up");
    delete pCollection;
    return 0;
}
