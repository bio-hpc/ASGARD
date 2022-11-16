/*!
   \file mmE.cpp

   \brief Determines the MM energy of a pdb file

   \author Martin B. Peters

   $Date: 2010/04/22 22:19:54 $
   $Revision: 1.9 $

   ----------------------------------------------------------------------------

   MTK++ - C++ package of modeling libraries.

   Copyright (C) 2005-2006  (see AUTHORS file for a list of contributors)

   ----------------------------------------------------------------------------
*/
#include "Utils/printHeader.h"
#include "Utils/vector3d.h"

#include "Log/errorHandler.h"

#include "Minimizers/lbfgs.h"

// - MOLECULE
#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/atom.h"
#include "Molecule/bond.h"
#include "Molecule/atomTyper.h"
#include "Molecule/connections.h"

// - PARSERS
#include "Parsers/elementParser.h"
#include "Parsers/pdbParser.h"
#include "Parsers/paramParser.h"
#include "Parsers/stdLibParser.h"
#include "Parsers/commLineOptions.h"

// - AMBER
#include "MM/amber.h"

// - TIME
#include "time.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <vector>

using namespace MTKpp;

int main (int argc, char *argv[])
{
    std::string prog_name = "mmE";
    std::vector<std::string> authors;
    std::string author = "Martin B. Peters";
    authors.push_back(author);

    // 1. CREATE AN OBJECT
    commLineOptions *clo = new commLineOptions();

    // 2. SET PREFERENCES
    clo->noUsage();

    // 3. SET THE USAGE/HELP
    clo->addUsage( "  mmE: Calculates the AMBER Energy/Gradients          \n" );
    clo->addUsage( "  usage: mmE [options] pdbFile                        \n" );
    clo->addUsage( "  options:                                             " );
    clo->addUsage( "          -s Standard Library XML File                 " );
    clo->addUsage( "          -f Frcmod XML File                           " );
    clo->addUsage( "          -c Non Bonded Cutoff [100.0]                 " );
    clo->addUsage( "          -b Calculate Bond Energy [1]                 " );
    clo->addUsage( "          -a Calculate Angle Energy [1]                " );
    clo->addUsage( "          -t Calculate Torsion Energy [1]              " );
    clo->addUsage( "          -i Calculate Improper Energy  [1]            " );
    clo->addUsage( "          -n Calculate Non Bonded Energy [1]           " );
    clo->addUsage( "          -w Calculate H Bond Energy [0]             \n" );
    clo->addUsage( "          -m Minimize [0]                              " );
    clo->addUsage( "            - 0 None                                   " );
    clo->addUsage( "            - 1 Hydrogens Only                         " );
    clo->addUsage( "            - 2 All Atoms                              " );
    clo->addUsage( "          -k Minimize Method [2]                       " );
    clo->addUsage( "            - 0 Steepest Descents                      " );
    clo->addUsage( "            - 1 Conjugate Gradient (non implemented)   " );
    clo->addUsage( "            - 2 lBFGS                                  " );
    clo->addUsage( "          -q Minimize steps [100]                      " );
    clo->addUsage( "          -p Write Output Every N Steps [1]            " );
    clo->addUsage( "          -o Output file                               " );
    clo->addUsage( "          -z Log file                                \n" );
    clo->addUsage( "    flags:                                             " );
    clo->addUsage( "          -v Verbose [1]                               " );
    clo->addUsage( "          -h Help                                      " );

    // 4. SET THE OPTION STRINGS/CHARACTERS
    clo->setOption(  "stdLib",    's' );
    clo->setOption(  "frcmod",    'f' );
    clo->setOption(  "cutoff",    'c' );
    clo->setOption(  "bond",      'b' );
    clo->setOption(  "angle",     'a' );
    clo->setOption(  "torsion",   't' );
    clo->setOption(  "improper",  'i' );
    clo->setOption(  "nonbonded", 'n' );
    clo->setOption(  "hbond",     'w' );
    clo->setOption(  "minimize",  'm' );
    clo->setOption(  "minMethod", 'k' );
    clo->setOption(  "minSteps",  'q' );
    clo->setOption(  "prtSteps",  'p' );
    clo->setOption(  "output",    'o' );
    clo->setOption(  "log",       'z' );
    clo->setFlag  (  "verbose",   'v' );
    clo->setFlag  (  "help",      'h' );

    // 5. PROVIDE THE COMMANDLINE
    clo->processCommandArgs( argc, argv );

    clo->usageOn();

    // 6. GET THE VALUES
    if ( clo->getFlag( "help" ) || clo->getFlag( 'h' ) ) {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      return 0;
    }

    std::string pdbFile = "";
    std::string libXml = "";
    std::string frcmodXml = "";
    std::string outputFile = "mmE.out";
    std::string logFile = "mmE.log";
    double nonbondedCutOff  = 100.0;
    int calcBond      = 1;
    int calcAngle     = 1;
    int calcTorsion   = 1;
    int calcImproper  = 1;
    int calcNonBonded = 1;
    int calcHBond     = 0;
    int verbose       = 1;
    int minimize      = 0;
    int minMethod     = 2;
    int minSteps      = 100;
    int prtSteps      = 1;

    if (clo->getArgc() < 1) {
      printHeader(std::cout, prog_name, authors);
      clo->printUsage();
      return 0;
    }
    else {
      pdbFile = clo->getArgv(0);
    }

    if ( clo->getValue( "s" ) != 0 ) {
      libXml = clo->getValue( "s" );
    }
    else if ( clo->getValue( "stdLib" ) != 0 ) {
      libXml =  clo->getValue( "stdLib" );
    }

    if ( clo->getValue( "f" ) != 0 ) {
      frcmodXml = clo->getValue( "f" );
    }
    else if ( clo->getValue( "frcmod" ) != 0 ) {
      frcmodXml = clo->getValue( "frcmod" );
    }

    if ( clo->getValue( "c" ) != 0 ) {
      nonbondedCutOff = strtod(clo->getValue( "c" ), 0);
    }
    else if ( clo->getValue( "cutoff" ) != 0 ) {
      nonbondedCutOff = strtod(clo->getValue( "cutoff" ), 0);
    }

    if ( clo->getValue( "b" ) != 0 ) {
      calcBond = atoi(clo->getValue( "b" ));
    }
    else if ( clo->getValue( "bond" ) != 0 ) {
      calcBond = atoi(clo->getValue( "bond" ));
    }

    if ( clo->getValue( "a" ) != 0 ) {
      calcAngle = atoi(clo->getValue( "a" ));
    }
    else if ( clo->getValue( "angle" ) != 0 ) {
      calcAngle = atoi(clo->getValue( "angle" ));
    }

    if ( clo->getValue( "t" ) != 0 ) {
      calcTorsion =  atoi(clo->getValue( "t" ));
    }
    else if ( clo->getValue( "torsion" ) != 0 ) {
      calcTorsion =  atoi(clo->getValue( "torsion" ));
    }

    if ( clo->getValue( "i" ) != 0 ) {
      calcImproper =  atoi(clo->getValue( "i" ));
    }
    else if ( clo->getValue( "improper" ) != 0 ) {
      calcImproper =  atoi(clo->getValue( "improper" ));
    }

    if ( clo->getValue( "n" ) != 0 ) {
      calcNonBonded =  atoi(clo->getValue( "n" ));
    }
    else if ( clo->getValue( "nonbonded" ) != 0 ) {
      calcNonBonded =  atoi(clo->getValue( "nonbonded" ));
    }

    if ( clo->getValue( "w" ) != 0 ) {
      calcHBond =  atoi(clo->getValue( "w" ));
    }
    else if ( clo->getValue( "hbond" ) != 0 ) {
      calcHBond =  atoi(clo->getValue( "hbond" ));
    }

    if ( clo->getValue( "m" ) != 0 ) {
      minimize =  atoi(clo->getValue( "m" ));
    }
    else if ( clo->getValue( "minimize" ) != 0 ) {
      minimize =  atoi(clo->getValue( "minimize" ));
    }

    if ( clo->getValue( "k" ) != 0 ) {
      minMethod =  atoi(clo->getValue( "k" ));
    }
    else if ( clo->getValue( "minMethod" ) != 0 ) {
      minMethod =  atoi(clo->getValue( "minMethod" ));
    }

    if ( clo->getValue( "q" ) != 0 ) {
      minSteps =  atoi(clo->getValue( "q" ));
    }
    else if ( clo->getValue( "minSteps" ) != 0 ) {
      minSteps =  atoi(clo->getValue( "minSteps" ));
    }

    if ( clo->getValue( "p" ) != 0 ) {
      prtSteps =  atoi(clo->getValue( "p" ));
    }
    else if ( clo->getValue( "prtSteps" ) != 0 ) {
      prtSteps =  atoi(clo->getValue( "prtSteps" ));
    }

    if ( clo->getValue( "o" ) != 0 ) {
      outputFile = clo->getValue( "o" );
    }
    else if ( clo->getValue( "output" ) != 0 ) {
      outputFile = clo->getValue( "output" );
    }

    if ( clo->getValue( "z" ) != 0 ) {
      logFile = clo->getValue( "z" );
    }
    else if ( clo->getValue( "log" ) != 0 ) {
      logFile = clo->getValue( "log" );
    }

    if ( clo->getValue( "v" ) != 0 ) {
      verbose =  atoi(clo->getValue( "v" ));
    }
    else if ( clo->getValue( "verbose" ) != 0 ) {
      verbose =  atoi(clo->getValue( "verbose" ));
    }

    // OPEN LOG FILE
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

    std::string errMessage = "";

    //puts("  mmE::Open output file");
    std::ofstream oOutput;
    oOutput.open(outputFile.c_str());
    if (!oOutput) {
      std::cout << "\nUNABLE TO OUTPUT FILE" << std::endl;
      return 0;
    }

    printHeader(oOutput, prog_name, authors);

    if (verbose) {
      oOutput << " PDB File                   : " << pdbFile << std::endl;
      oOutput << " Library File               : " << libXml << std::endl;
      oOutput << " Frcmod File                : " << frcmodXml << std::endl;
      oOutput << " Non Bonded Cut Off         : " << nonbondedCutOff << std::endl;
      oOutput << " Calculate Bond Energy      : " << calcBond << std::endl;
      oOutput << " Calculate Angle Energy     : " << calcAngle << std::endl;
      oOutput << " Calculate Torsion Energy   : " << calcTorsion << std::endl;
      oOutput << " Calculate Improper Energy  : " << calcImproper << std::endl;
      oOutput << " Calculate NonBonded Energy : " << calcNonBonded << std::endl;
      oOutput << " Calculate H-Bond Energy    : " << calcHBond << std::endl;
      oOutput << " Minimize                   : " << minimize << std::endl;
    }

    //! START & END TIME
    time_t                   startTime;
    time_t                   endTime;

    // GET TIME
    time (&startTime);

    std::string AMBERHOME = getenv("AMBERHOME");

    collection* pCollection = new collection();

    // ELEMENTS
    errMessage = " Read Element Data";
    MTKpp::errorLogger.throwError("mmE", errMessage, INFO);

    elementParser* pElementParser = new elementParser(pCollection->pElements);
    std::string elementXmlFile = AMBERHOME + "/dat/mtkpp/elements.xml";
    pElementParser->Read(elementXmlFile);

    // PARAMETERS
    errMessage = " Read MM parameters";
    MTKpp::errorLogger.throwError("mmE", errMessage, INFO);

    pCollection->addParameters();
    paramParser* pParamParser = new paramParser(pCollection->getParameters());
    std::string parameterXmlFile = AMBERHOME + "/dat/mtkpp/parm94.xml";
    pParamParser->Read(parameterXmlFile);
    parameterXmlFile = AMBERHOME + "/dat/mtkpp/parm_gaff.xml";
    pParamParser->Read(parameterXmlFile);
    parameterXmlFile = AMBERHOME + "/dat/mtkpp/metals/metalParm.xml";
    pParamParser->Read(parameterXmlFile);

    // FORCEFIELD MOD FILE
    if (frcmodXml != "") {
      pParamParser->Read(frcmodXml);
    }

     // STANDARD LIBRARIES
    errMessage = " Reading MM Library";
    MTKpp::errorLogger.throwError("mmE", errMessage, INFO);

    pCollection->addStdLibrary();
    parameters* pParms = pCollection->getParameters();
    stdLibParser* pStdLibParser = new stdLibParser(pCollection->getStdLibrary(), pParms);
    std::string stdLibXmlFile = AMBERHOME + "/dat/mtkpp/aminont94.xml";
    pStdLibParser->Read(stdLibXmlFile);
    stdLibXmlFile = AMBERHOME + "/dat/mtkpp/amino94.xml";
    pStdLibParser->Read(stdLibXmlFile);
    stdLibXmlFile = AMBERHOME + "/dat/mtkpp/aminoct94.xml";
    pStdLibParser->Read(stdLibXmlFile);
    stdLibXmlFile = AMBERHOME + "/dat/mtkpp/fragLib/terminal.xml";
    pStdLibParser->Read(stdLibXmlFile);

    if (libXml != "") {
      pStdLibParser->Read(libXml);
    }

    // READ PDB FILE
    errMessage = " Reading PDB File";
    MTKpp::errorLogger.throwError("mmE", errMessage, INFO);

    pdbParser* pPdbParser = new pdbParser();
    pPdbParser->Read(pdbFile, pCollection);
    std::string inputBaseName = baseName(pdbFile);

    // ASSIGN ATOM TYPES
    errMessage = " Assign MM Atom Types";
    MTKpp::errorLogger.throwError("mmE", errMessage, INFO);

    atomTyper* pAtomTyper;
    pAtomTyper = new atomTyper(0);
    pAtomTyper->atomTypeByLib(pCollection);

    // ASSIGN CONNECTIONS
    errMessage = " Assign Connections";
    MTKpp::errorLogger.throwError("mmE", errMessage, INFO);

    connections* pConnections;
    pConnections = new connections(pCollection);
    pConnections->run();

    // AMBER POTENTIAL
    errMessage = " Setting up AMBER Potential";
    MTKpp::errorLogger.throwError("mmE", errMessage, INFO);

    MTKpp::amber* pAmber = new MTKpp::amber();
    // turn on gradient calc
    if (minimize) {
      pAmber->calcForces(1);
    }
    pAmber->setPotential(calcBond, calcAngle, calcTorsion, calcImproper,
                         calcNonBonded, calcHBond);

    // INITIALIZE SOME VARIABLES
    int failure = 0;
    int nAtoms = 0;
    int nUniqueTypes = 0;
    int nExcludedAtoms = 0;
    int nExcluded14Atoms = 0;
    int nResidues = 0;

    int nBonds = 0;
    int nAngles = 0;
    int nTorsions = 0;
    int nImpropers = 0;

    // START ALLOCATION //
    // SET NUMBER OF ATOMS
    nAtoms = pCollection->getNumAtoms();
    pAmber->setNumAtoms(nAtoms);
    if (verbose) {
      oOutput << " Number of Atoms     = " << nAtoms << std::endl;
    }

    // SET NUMBER OF UNIQUE ATOM TYPES
    nUniqueTypes = pCollection->getNumUniqueAtomTypes();
    failure = pAmber->setNumTypes(nUniqueTypes);
    if (failure) {
      std::cout << " Error setting number of unique atom types" << std::endl;
      return 1;
    }

    // SET NUMBER OF EXCLUDED ATOMS
    nExcludedAtoms = pCollection->getNumExcludedAtoms();
    nExcluded14Atoms = pCollection->getNumExcluded14Atoms();
    failure = pAmber->setExcludedSize(nExcludedAtoms, nExcluded14Atoms);
    if (failure) {
      std::cout << " Error setting number of excluded atoms" << std::endl;
      return 1;
    }

    // SET NUMBER OF RESIDUES
    nResidues = pCollection->getNumberSubMolecules();
    failure = pAmber->setNumResidues(nResidues);
    if (failure) {
      std::cout << " Error setting number of residues" << std::endl;
      return 1;
    }

    // BONDS
    if (calcBond) {
      nBonds = pCollection->getNumBonds();
      pAmber->setNumBonds(nBonds);
      if (verbose) {
        oOutput << " Number of Bonds     = " << nBonds << std::endl;
      }
    }

    // ANGLES
    if (calcAngle) {
      nAngles = pCollection->getNumAngles();
      pAmber->setNumAngles(nAngles);
      if (verbose) {
        oOutput << " Number of Angles    = " << nAngles <<  std::endl;
      }
    }

    // TORSIONS
    if (calcTorsion) {
      nTorsions = pCollection->getNumMMTorsions();
      pAmber->setNumTorsions(nTorsions);
      if (verbose) {
        oOutput << " Number of Torsions  = " << nTorsions << std::endl;
      }
    }

    // IMPROPERS
    if (calcImproper) {
      nImpropers = pCollection->getNumMMImpropers();
      pAmber->setNumImpropers(nImpropers);
      if (verbose) {
        oOutput << " Number of Impropers = " << nImpropers << std::endl;
      }
    }
    // END ALLOCATION //

    errMessage = " Memory Allocation complete, starting to load data";
    MTKpp::errorLogger.throwError("mmE", errMessage, INFO);

    // START LOAD DATA //
    // LOAD COORDINATES FROM MOLECULE TO MM
    failure = pCollection->getCoordinates(pAmber->getCoords());
    if (failure) {
      std::cout << " Error loading coordinates " << std::endl;
      return 0;
    }

    // LOAD ATOM CHARGES FROM MOLECULE TO MM
    failure = pCollection->getMMCharges(pAmber->getCharges());
    if (failure) {
      std::cout << " Error loading charges " << std::endl;
      return 0;
    }

    // LOAD ATOM SYMBOLS FROM MOLECULE TO MM
    failure = pCollection->getAtomSymbols(pAmber->getSymbols());
    if (failure) {
      std::cout << " Error loading symbols " << std::endl;
      return 0;
    }

    // LOAD ATOM INTEGER TYPES FROM MOLECULE TO MM
    failure = pCollection->getAtomTypes(pAmber->getIntTypes());
    if (failure) {
      std::cout << " Error loading atom integer types " << std::endl;
      return 0;
    }

    // LOAD ATOM CHARACTER TYPES FROM MOLECULE TO MM
    failure = pCollection->getAtomTypes(pAmber->getCharTypes());
    if (failure) {
      std::cout << " Error loading atom character types " << std::endl;
      return 0;
    }

    // LOAD BONDS
    failure = pCollection->getBonds(pAmber->getBonds());
    if (failure) {
      std::cout << " Error loading bonds " << std::endl;
      exit(0);
    }

    failure = pCollection->getBondParams(pAmber->getBondParams());
    if (failure) {
      std::cout << " Error loading bonds " << std::endl;
      exit(0);
    }

    failure = pCollection->getAngles(pAmber->getAngles());
    if (failure) {
      std::cout << " Error loading angles " << std::endl;
      exit(0);
    }

    failure = pCollection->getAngleParams(pAmber->getAngleParams());
    if (failure) {
      std::cout << " Error loading angles " << std::endl;
      exit(0);
    }

    failure = pCollection->getMMTorsions(pAmber->getTorsions(),
                                         pAmber->getTorsionParams());
    if (failure) {
      std::cout << " Error loading torsions " << std::endl;
      exit(0);
    }

    failure = pCollection->getMMImpropers(pAmber->getImpropers(),
                                          pAmber->getImproperParams());
    if (failure) {
      std::cout << " Error loading impropers " << std::endl;
      exit(0);
    }

    // LOAD L-J PARAMETERS FROM MOLECULE TO MM
    failure = pCollection->getLJParams(pAmber->getR6Params(), pAmber->getR12Params());
    if (failure) {
      std::cout << " Error loading LJ parameters " << std::endl;
      exit(0);
    }

    // LOAD THE NUMBER OF EXCLUDED ATOMS(1-2, 1-3, 1-4) FOR EACH ATOM
    failure = pCollection->getNumExcludedAtoms(pAmber->getNumExcluded());
    if (failure) {
      std::cout << " Error loading number of excluded atoms " << std::endl;
      exit(0);
    }

    // LOAD THE EXCLUDED ATOMS(1-2, 1-3, 1-4) FOR EACH ATOM
    failure = pCollection->getExcludedAtoms(pAmber->getExcluded());
    if (failure) {
      std::cout << " Error loading excluded atoms  " << std::endl;
      exit(0);
    }

    // LOAD THE NUMBER OF EXCLUDED 1-4 ATOMS FROM EACH ATOM
    failure = pCollection->getNumExcluded14Atoms(pAmber->getNumExcluded14());
    if (failure) {
      std::cout << " Error loading number of 1-4 excluded atoms " << std::endl;
      exit(0);
    }

    // LOAD THE EXCLUDED 1-4 ATOMS FROM EACH ATOM
    failure = pCollection->getExcluded14Atoms(pAmber->getExcluded14());
    if (failure) {
      std::cout << " Error loading 1-4 excluded atoms " << std::endl;
      exit(0);
    }
    // END LOAD DATA //

    time (&endTime);
    int diffTime = (int) difftime(endTime, startTime);
    if (verbose) {
      oOutput << " MM Setup took " << diffTime << " seconds." << std::endl;
    }

    errMessage = " Checking Structure is OK";
    MTKpp::errorLogger.throwError("mmE", errMessage, INFO);

    // Structures xyz coordinates
    double *xyz;
    xyz = pAmber->getCoords();

    double d = 0.0;
    // Test Structure
    for (int i = 0; i < nAtoms; i++) {
      for (int j = 0; j < nAtoms; j++) {
        if (i == j) continue;
        d = sqrt(pow((xyz[i*3  ] - xyz[j*3  ]),2) +
                 pow((xyz[i*3+1] - xyz[j*3+1]),2) +
                 pow((xyz[i*3+2] - xyz[j*3+2]),2) );
        if (d < 0.1) {
          errMessage = " Two Atoms are less than 0.1 Ang Apart ... exiting";
          MTKpp::errorLogger.throwError("mmE", errMessage, MTK_ERROR);

          exit(0);
        }
      }
    }

    ///////////////////////////////////////////////////////////////////////////
    /////////////////// STEEPEST DESCENTS OPTIMIZATION ////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    if (minimize and (minMethod == 0)) {
      std::cout << " mmE::Steep Descents optimization " << std::endl;
      double bhigh = 0.5;
      double blow = 0.1;
      double alp = 1.0e-4;
      double tol = 1.0e-6;
      int itc = 1;
      double *grad;

      double *xyz_new;
      try {
        xyz_new = new double [nAtoms*3];
      }
      catch (std::bad_alloc) {
        std::cout << " Memory Allocation Failure " << std::endl;
        return 1;
      }

      for (int a = 0; a < nAtoms*3; a++) {
        xyz_new[a] = 0.0;
      }

      pAmber->resetGMatrix();
      pAmber->calcEnergy();
      double eCurrent = pAmber->getTotalEnergy();
      std::cout << " mmE::Current Energy = " << eCurrent << std::endl;

      xyz = pAmber->getCoords();
      grad = pAmber->getGradients();
      pAmber->calcForces(0);

      double gMax = 0.0;
      double gNorm = 0.0;
      int gMaxAtom = 0;
      pAmber->getGradNorm(gMax, gMaxAtom, gNorm);
      std::cout << " mmE::gMax = " << gMax << std::endl;
      std::cout << " mmE::gMaxAtom = " << gMaxAtom << std::endl;
      std::cout << " mmE::gNorm = " << gNorm << std::endl;

      while (gNorm > tol & itc <= minSteps) {
        std::cout << " mmE::iteration: " << itc << std::endl;
        itc++;
        double lambda = std::min(0.0001, 100.0/(1+gNorm));
        std::cout << " mmE::lambda = " << lambda << std::endl;
        for (int a = 0; a < nAtoms*3; a++) {
          xyz_new[a] = xyz[a] - lambda * grad[a];
          std::cout << " old xyz = " << xyz[a] << " new xyz = " << xyz_new[a] << std::endl;
        }

        double gc = 0.0;
        for (int a = 0; a < nAtoms*3; a++) {
          gc += (grad[a] * grad[a]);
        }
        std::cout << " mmE::gc = " << gc << std::endl;

        for (int a = 0; a < nAtoms*3; a++) {
          xyz[a] = xyz_new[a];
        }

        pAmber->calcEnergy();
        double eNew = pAmber->getTotalEnergy();

        double eGoal = eNew - alp * lambda * gc;
        double q0   = eCurrent;
        double qp0  = -gc;
        double lamc = lambda;
        double qc   = eNew;
        double qm   = 0.0;
        double lamm = 0.0;
        int iarm    = 0;

        std::cout << " mmE::New Energy = " << eNew << std::endl;
        std::cout << " mmE::Energy Goal = " << eGoal << std::endl;

        while (eNew > eGoal) {

          double lleft  = lamc * blow;
          double lright = lamc * bhigh; 

          double lplus = - qp0/(2 * lamc*(qc - q0 - qp0) );
          if (lplus < lleft) {
            lplus = lleft;
          }
          else if (lplus > lright) {
            lplus = lright;
          }
          lambda = lplus;

          qm = qc;
          lamm = lamc;
          lamc = lambda;
          std::cout << " lambda = " << lambda << std::endl;
          for (int a = 0; a < nAtoms*3; a++) {
            xyz_new[a] = xyz[a] - lambda * grad[a];
          }
          for (int a = 0; a < nAtoms*3; a++) {
            xyz[a] = xyz_new[a];
          }

          pAmber->calcEnergy();
          eNew = pAmber->getTotalEnergy();

          if (iarm > 10) {
            std::cout << " Error in SD  ... exiting " << std::endl;
            exit(0);
          }
          eGoal = eNew - alp * lambda * gc;
        }

        pAmber->resetGMatrix();
        pAmber->calcForces(1);
        pAmber->calcEnergy();

        pCollection->setCoordinates(xyz);
        std::string outFile = inputBaseName + "_min_" + ".pdb";
        pPdbParser->Write(outFile, pCollection);
      }

/*
      for (int i = 1; i < minSteps; i++) {
        pAmber->resetGMatrix();
        // Calculate Energy
        pAmber->calcEnergy();
        eCurrent = pAmber->getTotalEnergy();
        if (i == 1) {
          deltaE = eCurrent;
          deltaX = 100.0;
        }

        // Only minimize Hs
        if (minimize == 1) {
          pAmber->setNonHGradsToZero();
        }

        if (std::abs(deltaE) < eTest) {
          oOutput << " Energy Test Passed " << std::endl;
          break;
        }

        xyz = pAmber->getCoords();
        grad = pAmber->getGradients();
      }
*/
      pCollection->setCoordinates(xyz);
      std::string outFile = inputBaseName + "_min_" + ".pdb";
      pPdbParser->Write(outFile, pCollection);
    }

    ///////////////////////////////////////////////////////////////////////////
    //////////////////////////// LBFGS OPTIMIZATION ///////////////////////////
    ///////////////////////////////////////////////////////////////////////////
    if (minimize and (minMethod == 2)) {

      static int iflag = 0;
      static int diagco = 0;
      static int iprint[2] = {-1,0};
      static double eps = 1e-5;
      static double xtol = 1e-16;
      static int n = nAtoms*3;
      double *grad;
      int MLBFGS = 5;

      double *diag;
      try {
        diag = new double [n];
      }
      catch (std::bad_alloc) {
        std::cout << " Memory Allocation Failure: diag array for lbfgs" << std::endl;
        return 1;
      }

      double *w;
      try {
        w = new double [nAtoms*3*(2*MLBFGS+1)+2*MLBFGS];
      }
      catch (std::bad_alloc) {
        std::cout << " Memory Allocation Failure: W array for lbfgs" << std::endl;
        return 1;
      }

      int gMaxAtom = 0;
      double gSum  = 0.0;
      double gMax  = 0.0;
      double gNorm = 0.0;

      oOutput << " step bond angle torsion+improper vdw ele vdW14 ele14 eTotal";
      oOutput << " gradMaxAtom gradMax gradNorm " << std::endl;

      double eLast = 0.0;
      double eCurrent = 0.0;

      double deltaE = 0.0;
      double deltaX = 0.0;

      double eTest = 0.002;
      //double xTest = 0.001;
      //double gTest = 0.5;

      int prtCounter = 0;
      std::string outFile = "";
      for (int i = 1; i < minSteps; i++) {
        prtCounter++;
        pAmber->resetGMatrix();

        // Calculate Energy
        pAmber->calcEnergy();
        eCurrent = pAmber->getTotalEnergy();
        if (i == 1) {
          deltaE = eCurrent;
          deltaX = 100.0;
        }

        // Only minimize Hs
        if (minimize == 1) {
          pAmber->setNonHGradsToZero();
        }

        if (std::abs(deltaE) < eTest) {// && deltaE < 0.0) {
          oOutput << " Energy Test Passed " << std::endl;
          break;
        }

        xyz = pAmber->getCoords();
        grad = pAmber->getGradients();

        gMax = 0.0;
        gNorm = 0.0;
        gSum = 0.0;

        //std::cout << " GRADIENT: " << std::endl;
        //int x = 1;
        for (int j = 0; j < nAtoms*3; j++) {
          //std::cout << " " << grad[j];
          if (std::abs(grad[j]) > gMax) {
            gMax = grad[j];
            gSum += (grad[j]*grad[j]);
            gMaxAtom = j;
          }
/*
          if (x == 3) {
            std::cout << " " << std::endl;
            x = 0;
          }
          x++;
*/
        }

        gMaxAtom = (gMaxAtom-1)/3 + 1;
        gNorm = sqrt(gSum/double(nAtoms*3.0));

        std::string gradString = double2String(gMaxAtom) + " " +
                                 double2String(gMax) + " " +
                                 double2String(gNorm);
        pAmber->printEnergy2(oOutput, int2String(i), gradString);

        deltaE = eCurrent - eLast;
        eLast = eCurrent;
        //deltaX = ;

        // Call LBFGS
        lbfgs_(&n, &MLBFGS, xyz, &eCurrent, grad, &diagco, diag, iprint, &eps, &xtol, w, &iflag);

        if (iflag == 0) {
          break;
        }
        if (iflag < 0) {
          std::cout << " Error in lbfgs ... exiting " << std::endl;
          exit(0);
        }
        //pAmber->printCoords();
        xyz = pAmber->getCoords();
        pCollection->setCoordinates(xyz);

        if (prtCounter == prtSteps) {
          outFile = inputBaseName + "_opt_" + int2String(i) + ".pdb";
          pPdbParser->Write(outFile, pCollection);
          prtCounter = 0;
        }
      }
      outFile = inputBaseName + "_min_" + ".pdb";
      pPdbParser->Write(outFile, pCollection);
    }
    else {

      errMessage = " Calculating Energy ";
      MTKpp::errorLogger.throwError("mmE", errMessage, INFO);

      pAmber->calcEnergy();
    }

    if (verbose) {
      errMessage = " Printing Energy (verbose) ";
      MTKpp::errorLogger.throwError("mmE", errMessage, INFO);

      pAmber->printEnergy(oOutput);
    }
    else {
      errMessage = " Printing Energy ";
      MTKpp::errorLogger.throwError("mmE", errMessage, INFO);

      pAmber->printEnergy2(oOutput);
    }

    // - Clean up - //
    delete pCollection;
    delete pPdbParser;
    delete pParamParser;
    delete pStdLibParser;
    delete pAtomTyper;

    time (&endTime);
    int end_Time = (int) difftime(endTime, startTime);
    errMessage = "\n ###  MM Calculation took "
                    + int2String(end_Time) + " Seconds.\n";

    errMessage += " ###  MM Energy evaluation took "
                    + int2String(end_Time - diffTime) + " Seconds. ";

    MTKpp::errorLogger.throwError("mmE", errMessage, MESSAGE);
    oLog.close();

    return 0;
}
