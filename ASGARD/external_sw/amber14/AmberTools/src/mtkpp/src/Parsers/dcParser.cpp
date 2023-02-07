/*!
   \file dcParser.cpp
   \brief Parses Divcon files
   \author Martin Peters

   Reads divcon output files and write divcon input files

   $Date: 2010/03/29 20:39:34 $
   $Revision: 1.24 $

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

#include "dcParser.h"

#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/submolecule.h"
#include "Molecule/atom.h"
#include "Molecule/bond.h"
#include "Molecule/element.h"

#include "Utils/vector3d.h"
#include "StringManip.h"

#include <cctype>

namespace MTKpp
{

// ============================================================
// Function : dcParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
dcParser::dcParser():baseParser() {
    inputFileName  = "";
    outputFileName = "";
    pCollection = 0;
    pMolecule = 0;
    pSubMolecule = 0;
    pAtom = 0;
    coord = "cartesian";
    bStandard = true;
    bCluster = false;
    dBuff1 = 4.0;
    dBuff2 = 2.0;
    bResidue = false;
    bDirect = true;
    hamiltonian = "AM1";
    bChargeModel = false;
    chargeModel = "cm1";
    bCutBond = false;
    cutbond = 0.0;
    bShift = false;
    shift = 10.0;
    bGuess = false;
    bPDump = false;
    pDump = 0;
    bDump = false;
    dump = 0;
    iDouble = 0;
    bDouble = false;
    bMaxOpt = 0;
    maxopt = 0;
    bOpt = false;
    optimizer = "";
    bScrf = false;
    scrfScale=2.0;
    bWater = false;
    bOctanol = false;
    bNoOverlap = false;
    bScreen = false;
    bVdw = false;
    bDipole = false;
    bFreq = false;
    bThermo = false;
    integrals = "talman";
    bNMR = false;
    calnum = 0;
    bPwd = false;
    bPwdAtom = false;
    bPwdResidue = false;
    bPrtSub = false;
    bPrtVec = false;
    bPrtCoords = false;
    bPrtPar = false;
    bPrtVdw = false;
    bAddMM = true;
    bNoMM = false;
    bFullSCF = false;
    bMaxIt = false;
    maxIt = 100;
    bCutRepul = false;
    cutRepul = 0;
    bIntegrals = false;
    bTempK = false;
    tempK = 0;
    bMaxTime = false;
    maxTime = 0;
    bMinR = false;
    minR = 10.0;
    bQMAlign = false;
    bChkRes = false;
    bDOS = false;
    bIP = false;
    bHomoLumo = false;
    bZmake = false;
    imult = 1;
    bUhf = false;
    bExternal = false;
    setError(0);
}

// ============================================================
// Function : dcParser()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
dcParser::~dcParser() {}

// ============================================================================
//
//                           R E A D  F U N C T I O N S
//
// ============================================================================

// ============================================================
// Function : Read
// ------------------------------------------------------------
// Parsers DivCon output files
// ------------------------------------------------------------
void dcParser::Read(const std::string &dcfile, molecule* pMolecule)
{
    setError(0);
    if ((dcfile == "") or (!pMolecule)){
      std::cout << "\ndcParser read error: molecule/file does not exist" << std::endl;
      setError(1);
      return;
    }

    outputFileName = baseName(dcfile);
    std::ifstream idc;
    idc.open(dcfile.c_str());

    if (!idc) {
      std::cout << "\nUNABLE TO OPEN DivCon Output FILE"
                << "\nFILENAME = " << dcfile << std::endl;
      setError(1);
      return;
    }

    std::string fileline;

    // - Read the file - //
    while (idc) {
      getline(idc,fileline);

      if (fileline.substr(0,33) == " DivCon Calculation performed on:") {
        std::vector<std::string> words;
        splitString(fileline, ":", words, 0);
        pMolecule->setName(stripString(words[1]," "));
      }

      //Check for error in divcon output file
      if (fileline.substr(0,68) ==
        " WARNING: NO CONVERGENCE IN SCF CALCULATION; MAXIMUM ITERATION COUNT")
      {
        pMolecule->bError = true;
        setError(1);
        return;
      }

      std::string::size_type loc = fileline.find( "ERROR", 0 );
      if ( loc != std::string::npos ) {
        pMolecule->bError = true;
        setError(1);
        return;
      }

      if (fileline.substr(0,58) == "    NUMBER     SYMBOL            X           Y           Z") {
        this->readCoordinates(idc, pMolecule);
      }
      if (fileline.substr(0,8) == " CHARGE="){
        int totalCharge = 0;
        std::vector<std::string> wordsCHG;
        splitString(fileline, "=", wordsCHG, 0);
        std::vector<std::string> words2;
        splitString(wordsCHG[1], " ", words2, 0);
        totalCharge= (atoi(words2[0].c_str()));
        pMolecule->setTotalCharge(totalCharge);
      }
      if (fileline.substr(0,38) == " #### Q M - Q S A R   O U T P U T ####") {
        this->readQMQSAR(idc, pMolecule);
      }

      if (fileline.substr(0,52) == "****               NUCLEAR MAGNETIC SHIELDING TENSOR") {
        this->readNMR(idc, pMolecule);
      }

      if (fileline.substr(0,24) == " ELECTRONIC ENERGY     =") {
        const  std::string property = "EEnergy";
        double EEnergy = ReadElecEnergy(dcfile);
        const  std::string energy = double2String(EEnergy);
        pMolecule->addProperty(property, energy);
      }
    }
}

// ============================================================
// Function : ReadHOF()
// ------------------------------------------------------------
// Parsers HOF information from DivCon files
// ------------------------------------------------------------
double dcParser::ReadHOF(const std::string &dcfile)
{
    setError(0);
    std::ifstream idc;
    idc.open(dcfile.c_str());

    if (!idc) {
      std::cout << "\nUNABLE TO OPEN DivCon Output FILE"
                << "\nFILENAME = " << dcfile << std::endl;
      setError(1);
      return 0.0;
    }

    std::string fileline;
    std::vector<std::string> words;
    double hof = 0.0;

    // - Read the file - //
    while (idc) {
      getline(idc,fileline);

      //Check for error in divcon output file
      if (fileline.substr(0,68) ==
        " WARNING: NO CONVERGENCE IN SCF CALCULATION; MAXIMUM ITERATION COUNT")
      {
        setError(1);
        setErrorMessage("No Convergence in SCF");
        return hof;
      }

      std::string::size_type loc = fileline.find( "ERROR", 0 );
      if ( loc != std::string::npos ) {
        setError(1);
        setErrorMessage("Error Found");
        return hof;
      }

      if (fileline.substr(0,24) == " HEAT OF FORMATION     =") {
        splitString(fileline, " ", words, 0);
        hof = strtod(words[4].c_str(), 0);
      }
    }
    return hof;
}

// ============================================================
// Function : ReadElecEnergy()
// ------------------------------------------------------------
// Parsers ELECTRONIC ENERGY information from DivCon files
// ------------------------------------------------------------
double dcParser::ReadElecEnergy(const std::string &dcfile)
{
    setError(0);
    std::ifstream idc;
    idc.open(dcfile.c_str());

    if (!idc) {
      std::cout << "\nUNABLE TO OPEN DivCon OUTPUT FILE"
           << "\nFILENAME = " << dcfile << std::endl;
      setError(1);
      setErrorMessage("Unable to open file");
      return 0.0;
    }

    std::string fileline;
    std::vector<std::string> words;
    double EEnergy = 0.0;

    // - Read the file - //
    while (idc) {
      getline(idc,fileline);

      //Check for error in divcon output file
      if (fileline.substr(0,68) ==
        " WARNING: NO CONVERGENCE IN SCF CALCULATION; MAXIMUM ITERATION COUNT")
      {
        setError(1);
        return EEnergy;
      }

      std::string::size_type loc = fileline.find( "ERROR", 0 );
      if ( loc != std::string::npos ) {
        setError(1);
        return EEnergy;
      }

      if (fileline.substr(0,24) == " ELECTRONIC ENERGY     =") {
        splitString(fileline, " ", words, 0);
        EEnergy = strtod(words[3].c_str(), 0);
      }
    }
    return EEnergy;
}

// ============================================================
// Function : readCoordinates()
// ------------------------------------------------------------
// Parsers coordinate information from DivCon files
// ------------------------------------------------------------
void dcParser::readCoordinates(std::ifstream &idc, molecule* pMolecule)
{
    std::vector<std::string> words;
    int FileID;
    std::string element;
    double x,y,z;
    std::string fileline;

    getline(idc,fileline);
    getline(idc,fileline);
    std::vector<atom*> atomList = pMolecule->getAtomList();

    if (atomList.size() == 0){
      pSubMolecule = pMolecule->addSubMolecule();
      pSubMolecule->setName("NULL");
      pSubMolecule->setSubMolId(pMolecule->getNumSubMolecules());
    }

    while(fileline != "") {
      atom* pAtom;
      splitString(fileline, " ", words, 0);
      FileID    = (atoi(words[0].c_str()));
      element   = words[1];
      x         = strtod(words[2].c_str(), 0);
      y         = strtod(words[3].c_str(), 0);
      z         = strtod(words[4].c_str(), 0);
//      if (element.size() == 2) {
//        int b = tolower(element.substr(1,1).c_str());
//      }

      if (atomList.size() == 0){
        pAtom = pSubMolecule->addAtom();
        pAtom->setFileID(FileID);
        pAtom->setElement(pMolecule->getParent()->pElements->getElement(element));
      }
      else{
        pAtom = pMolecule->getAtom(FileID,1,0);
        if (!pAtom){
          std::cout <<  " dcParser:: Atom in molecule does not exist" << std::endl;
        }
      }
      if (pAtom->getElementSymbol() == "XX"){
        pAtom->setElement(pMolecule->getParent()->pElements->getElement(element));
      }
      pAtom->setCoords(x,y,z);
      words.clear();
      getline(idc,fileline);
    }
}

// ============================================================
// Function : readNMR()
// ------------------------------------------------------------
// Parsers NMR information from DivCon files
// ------------------------------------------------------------
void dcParser::readNMR(std::ifstream &idc, molecule* pMolecule)
{
    std::string fileline;
    std::vector<std::string> words;

    getline(idc,fileline);

    int fileID = 1;

    while(fileline != "") {
      splitString(fileline, " ", words, 0);
      const std::string NMRshield = "DivCon NMR shield";

      if (words[1] == "ATOM"){
        fileID = (atoi(words[2].c_str()));
        pAtom = pMolecule->getAtom(fileID, 0, 1);
      }
      if (words[1] == "AVERAGE") {
        double shieldConst = strtod(words[3].c_str(), 0);
        pAtom->addProperty(NMRshield, shieldConst);
      }
      words.clear();
      getline(idc,fileline);
    }
}

// ============================================================
// Function : readQMQSAR()
// ------------------------------------------------------------
// Parsers QMQSAR information from DivCon files
// ------------------------------------------------------------
void dcParser::readQMQSAR(std::ifstream &idc, molecule* pMolecule)
{
/*
    std::string fileline;
    std::vector<std::string> words;
    int natoms,norbitals; //,nelectrons;
    int atomicNum, atOrbs, atPqn; // FileID
    std::string element;
    double x,y,z,atZchg;

    for (int i = 0; i < 4; i++) {
      getline(idc,fileline);
    }
    splitString(fileline," ",words,0);
    natoms = (atoi(words[0].c_str()));
    if (natoms != pMolecule->getNumAtoms() ) {
      return;
    }
    norbitals  = ( atoi(words[1].c_str()) );
    //nelectrons = ( atoi(words[2].c_str()) );

    //pMolecule->setNumOrbitals(norbitals);
    //pMolecule->setNumElectrons(nelectrons);
    //pEigenVectorMatrix = pMolecule->setEigenVectorMatrix(norbitals);
    //pEigenValueMatrix  = pMolecule->setEigenValueMatrix(norbitals);

    words.clear();

    for (int j=0; j < 3; j++) {
      getline(idc,fileline);
    }

    for (int k=0; k < natoms; k++) {
      getline(idc,fileline);
      splitString(fileline," ",words,0);
      //FileID    = (atoi(words[0].c_str()));
      element   = words[1];
      pAtom = pMolecule->getAtom(k+1);
      if ( pAtom->getElement()->symbol != element ) {
        std::cout << "Atom # " << k << " in the divcon file is different from that in the mol file" << std::endl;
        return;
      }

      x = string2Double(words[2]);
      y = string2Double(words[3]);
      z = string2Double(words[4]);
      atomicNum = (atoi(words[5].c_str()));
      atOrbs    = (atoi(words[6].c_str()));
      atPqn     = (atoi(words[7].c_str()));
      atZchg    = string2Double(words[8]);

      pAtom->setCoords(x,y,z);
      pAtom->setAtomicNum(atomicNum);
      pAtom->setPrcpleQtmNum(atPqn);
      pAtom->setZcharge(atZchg);
      pAtom->setOrbitalNum(atOrbs);
      words.clear();
    }

    for (int m = 0; m < norbitals;) {
      getline(idc,fileline);
      getline(idc,fileline);
      getline(idc,fileline);

      splitString(fileline," ",words,0);

      // - Get EigenValues - //
      //int temp = m;
      for (unsigned int p=1; p < words.size(); p++) {
        //(*pEigenValueMatrix)[m][0] = (atof(words[p].c_str()));
        m = m + 1;
      }
      words.clear();

      // - Get EigenVectors - //
      getline(idc,fileline);

      for (int n = 0; n < norbitals; n++) {
        getline(idc,fileline);
        splitString(fileline," ",words,0);
        for (unsigned int w=1; w < words.size(); w++) {
          // -- EigenVectorMatrix:         -- //
          // -- Each Row is an EigenVector -- //
          //(*pEigenVectorMatrix)[temp+w-1][n] = (atof(words[w].c_str()));
        }
        words.clear();
      }
      words.clear();
    }
    // - Read semi-empirical parameter data - //
    getline(idc,fileline);
    getline(idc,fileline);
    getline(idc,fileline);
    splitString(fileline," ",words,0);
    int unique_atoms = (atoi(words[3].c_str()));

    words.clear();
    double atomData[8];

    typedef std::vector<atom*>::iterator     AtomIterator;
    std::vector<atom*>                       Atoms;
    //vector<orbital*>                    Orbitals;
    Atoms                                = pMolecule->getAtomList();
    getline(idc,fileline);

    for (int i = 0; i < unique_atoms; i++) {
      getline(idc,fileline);
      std::string temp = ( fileline.substr(1,2).c_str() );
      splitString(temp," ",words,0);
      std::string atom_symbol = words[0];
      words.clear();
      for (int j = 0; j < 8; j++) {
        getline(idc,fileline);
        splitString(fileline," ",words,0);
        atomData[j] = string2Double(words[4]);
        words.clear();
      }

      vector3d expnts(atomData[0],atomData[1],atomData[2]);
      vector3d al(atomData[3],atomData[4],atomData[5]);

      for (AtomIterator d=Atoms.begin(); d != Atoms.end(); d++) {
        pAtom = *d;
        if ( pAtom->getElement()->symbol == atom_symbol) {
          pAtom->setSEParameters(atomData);
        }
      }
      getline(idc,fileline);
    }
*/
}

// ============================================================================
//
//                          W R I T E  F U N C T I O N S
//
// ============================================================================

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write a DivCon input file.
// ------------------------------------------------------------
// 
// ============================================================
void dcParser::Write(const std::string &dcfile, collection* pCollection, bool &success)
{
    this->pCollection = pCollection;
    this->pMolecule = 0;
    bool ok = this->check();
    if (!ok) {
      success = ok;
      return;
    }
    inputFileName = baseName(dcfile);
    this->openInputFile(dcfile);
    this->writeHead();
    this->writeCoords();
    this->writeTail();
}

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write a DivCon input file.
// ============================================================
void dcParser::Write(const std::string &dcfile, molecule* pMolecule, bool &success)
{
    std::vector< vector3d > coordinates;
    pMolecule->getCoordinates(coordinates);
    this->Write(dcfile, pMolecule, coordinates, success);
}

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write a DivCon input file.
// ============================================================
void dcParser::Write(const std::string &dcfile, molecule* pMolecule,
                     std::vector< vector3d > &coordinates, bool &success)
{
    this->pMolecule = pMolecule;
    this->pCollection = 0;
    bool ok = this->check();
    if (!ok) {
      success = ok;
      return;
    }
    this->openInputFile(dcfile);
    this->writeHead();
    this->writeCoords(coordinates);
    this->writeTail();
}

// ============================================================
// Function : check()
// ------------------------------------------------------------
// Open input file
// ============================================================
bool dcParser::check()
{
    if (!this->pMolecule and !this->pCollection) {
      std::cout << "dcParser::check detected a problem " << std::endl;
      return false;
    }

    std::vector<atom*> Atoms;
    std::vector<molecule*> molList;

    if (this->pCollection) {
      molList = pCollection->getMoleculeList();
    }
    else {
      molList.push_back(this->pMolecule);
    }

    std::vector<std::string> eles;
    vector3d* coord1;
    vector3d* coord2;

    for (molIterator m = molList.begin(); m != molList.end(); m++) {
      molecule* pMol = *m;
      std::vector<submolecule*> submoles1 = pMol->getSubMoleculeList();
      for (sMolIterator d = submoles1.begin(); d != submoles1.end(); d++) {
        submolecule* pSmol = *d;
        std::vector<atom*> atoms1 = pSmol->getAtomList();
        for (AtomIterator x = atoms1.begin(); x != atoms1.end(); x++) {
          atom* pAt1 = *x;
          coord1 = pAt1->getCoords();

          for (molIterator m2 = molList.begin(); m2 != molList.end(); m2++) {
            molecule* pMol2 = *m2;
            std::vector<submolecule*> submoles2 = pMol2->getSubMoleculeList();
            for (sMolIterator d2 = submoles2.begin(); d2 != submoles2.end(); d2++) {
              submolecule* pSmol2 = *d2;
              std::vector<atom*> atoms2 = pSmol2->getAtomList();
              for (AtomIterator x2 = atoms2.begin(); x2 != atoms2.end(); x2++) {
                atom* pAt2 = *x2;
                if (pAt1 == pAt2) continue;
                coord2 = pAt2->getCoords();
                double at1at2Dist = coord1->dist(*coord2);
                if (at1at2Dist < 0.7) {
                  std::cout << pSmol->getName() << pSmol->getSubMolId() << " " << pAt1->getName()
                            << " is too close to "
                            << pSmol2->getName() << pSmol2->getSubMolId() << " " << pAt2->getName()
                            << std::endl;
                  std::cout << " dcParser::Two atoms are less than 0.7 Ang";
                  std::cout << " apart.\n      This calculation won't";
                  std::cout << " converge ... exiting ";
                  return false;
                }
              }
            }
          }
        }
      }
    }

    elements* pEle = 0;
    for (molIterator m = molList.begin(); m != molList.end(); m++) {
      molecule* pMol = *m;
      pEle = pMol->getParent()->pElements;
      std::vector<submolecule*> submolecules = pMol->getSubMoleculeList();

      for (sMolIterator d = submolecules.begin(); d != submolecules.end(); d++) {
        pSubMolecule = *d;
        Atoms = pSubMolecule->getAtomList();
        for (AtomIterator x = Atoms.begin(); x != Atoms.end(); x++) {
          pAtom = *x;
          strIterator a = std::find(eles.begin(), eles.end(), pAtom->getElementSymbol());
          if (a == eles.end()) {
            eles.push_back(pAtom->getElementSymbol());
          }
        }
      }
    }

    for (strIterator e = eles.begin(); e != eles.end(); e++) {
      std::string sElement = *e;
      bool gotIt = pEle->hasSEHamiltonian(sElement, hamiltonian);
      if (!gotIt) {
        std::cout << " The Hamiltonian: " << hamiltonian
                  << " is currently not available for the element "
                  << sElement << std::endl;
        return gotIt;
      }
    }
    return true;
}

// ============================================================
// Function : openInputFile()
// ------------------------------------------------------------
// Open input file
// ============================================================
void dcParser::openInputFile(const std::string &dcfile)
{
    this->inputFileStream.open(dcfile.c_str());

    if (!inputFileStream) {
      std::cout << "\n UNABLE TO OPEN DIVCON INPUT FILE"
                << "\nFILENAME = " << dcfile<< std::endl;
      return;
    }
}

// ============================================================
// Function : writeHead()
// ------------------------------------------------------------
// Write Head of a DivCon Input File
// ============================================================
void dcParser::writeHead()
{
    if (!this->inputFileStream) return;

    this->inputFileStream << hamiltonian << " &" << std::endl;
    if (bChargeModel) this->inputFileStream << chargeModel << " &" << std::endl;

    if (pMolecule) {
      this->inputFileStream << "charge = " << pMolecule->getFormalCharge() << " &" << std::endl;
    }
    else if (pCollection) {
      this->inputFileStream << "charge = " << pCollection->getFormalCharge() << " &" << std::endl;
    }

    if (bResidue) this->inputFileStream << "residue &" << std::endl;
    if (bDirect) this->inputFileStream << "direct &" << std::endl;

    if (!bStandard) {
      if (bCluster) this->inputFileStream << "cluster & " << std::endl;
      if (bCutBond) this->inputFileStream << "cutbond=" << std::showpoint << cutbond << " &" << std::endl;
    }
    if (bShift > 0) this->inputFileStream << "shift=" << shift << " &" << std::endl;

    if (bGuess) this->inputFileStream << "guess &" << std::endl;
    if (bPDump) this->inputFileStream << "pdump= " << pDump << "&" << std::endl;
    if (bDump) this->inputFileStream << "dump= " << dump << " &" << std::endl;

    if (bOpt) {
      this->inputFileStream << "opt=" << optimizer << " &" << std::endl;
      this->inputFileStream << "maxopt=" << maxopt << " &" << std::endl;
    }

    if (bScrf) {
      this->inputFileStream << "scrf &" << std::endl;
      this->inputFileStream << "scale=" <<  std::showpoint << scrfScale << " &" << std::endl;
      if (bWater) this->inputFileStream << "water & " << std::endl;
      if (bOctanol) this->inputFileStream << "octanol &" << std::endl;
    }

    if (bNoOverlap) this->inputFileStream << "no_overlap &" << std::endl;

    if (bPwd) this->inputFileStream << "pwd &" << std::endl;
    if (bPwdAtom) this->inputFileStream << "pwd_atom &" << std::endl;
    if (bPwdResidue) this->inputFileStream << "pwd_residue &" << std::endl;

    if (bScreen) this->inputFileStream << "screen &" << std::endl;
    //if (vdw.size() > 0) this->inputFileStream  << "vdw &" << std::endl;
    //if (iDouble > -1) this->inputFileStream  << "double = " << iDouble << " &" << std::endl;

    if (bStandard) this->inputFileStream << "standard &" << std::endl;

    if (bDipole) this->inputFileStream << "dipole &" << std::endl;

    if (bThermo) this->inputFileStream << "thermo &" << std::endl;

    if (bFreq) this->inputFileStream << "freq &" << std::endl;

    if (bNMR) this->inputFileStream << "nmr DPSCF=1.0D-9 &" << std::endl;

    if (bQMAlign) this->inputFileStream << "qmalign &" << std::endl;

    if (bIntegrals) this->inputFileStream << "intgls=" << integrals << " &" << std::endl;

    if (calnum) this->inputFileStream << "CALNUC=" << calnum << " & " << std::endl;

    if (bUhf) this->inputFileStream << "UHF IMULT=" << imult << " & " << std::endl;

    if (bExternal) this->inputFileStream << "EXTERNAL" <<" & " << std::endl;

    // Let coord be the last one
    if (!pMolecule) {
      this->inputFileStream << coord << " \n" << std::endl; 
    }

    if (pMolecule) {
      this->inputFileStream << coord << std::endl;
      this->inputFileStream << pMolecule->getName() << std::endl;
    }
}

// ============================================================
// Function : writeCoords()
// ------------------------------------------------------------
// Write Coordinate section of a DivCon Input File
// ============================================================
void dcParser::writeCoords()
{
    if (!this->inputFileStream) return;
    if (!this->pMolecule and !this->pCollection) return;

    std::vector<atom*> Atoms;
    //int counter = 1;
    std::vector<molecule*> molList;

    if (this->pCollection) {
      molList = pCollection->getMoleculeList();
    }
    else {
      molList.push_back(this->pMolecule);
    }

    if (bResidue) {
      int counter = 1;
      for (molIterator m = molList.begin(); m != molList.end(); m++) {
        molecule* pMol = *m;

        std::vector<submolecule*> submolecules = pMol->getSubMoleculeList();

        for (sMolIterator d = submolecules.begin(); d != submolecules.end(); d++) {
          pSubMolecule = *d;
          Atoms = pSubMolecule->getAtomList();
          bool firstAtom = true;
          for (AtomIterator x = Atoms.begin(); x != Atoms.end(); x++) {
            pAtom = *x;
            char temp[100];
            if (firstAtom) {
              sprintf(temp,"%6d %2s %18.10f%18.10f%18.10f %5s%3s%s%4s%s",
                      counter, (pAtom->getElement()->symbol.c_str()),
                      pAtom->getX(), pAtom->getY(), pAtom->getZ(),
                      "RES {", pAtom->getParent()->getName().c_str(), ";", pAtom->getName().c_str(), "}");
              firstAtom = false;
            }
            else {
              sprintf(temp,"%6d %2s %18.10f%18.10f%18.10f %5s%3s%s%4s%s",
                      counter, (pAtom->getElement()->symbol.c_str()),
                      pAtom->getX(), pAtom->getY(), pAtom->getZ(),
                      "    {", pAtom->getParent()->getName().c_str(), ";", pAtom->getName().c_str(), "}");
            }
            this->inputFileStream << temp << std::endl;
            counter++;
          }
        }
      }
    }
    this->inputFileStream << "END_COORD" << std::endl;
}

// ============================================================
// Function : writeCoords()
// ------------------------------------------------------------
// Write Coordinate section of a DivCon Input File
// ============================================================
void dcParser::writeCoords(std::vector< vector3d > &coordinates)
{
    if (!this->inputFileStream) return;
    if (!this->pMolecule) return;

    typedef std::vector<atom*>::iterator AtomIterator;
    std::vector<atom*>                   Atoms;
    int counter = 1;

    if (bResidue) {
      typedef std::vector<submolecule*>::iterator sMolIterator;
      std::vector<submolecule*> submolecules = pMolecule->getSubMoleculeList();

      //int counter = 1;
      for (sMolIterator d = submolecules.begin(); d != submolecules.end(); d++) {
        pSubMolecule = *d;
        Atoms = pSubMolecule->getAtomList();
        int firstAtom = true;
        for (AtomIterator x = Atoms.begin(); x != Atoms.end(); x++) {
          pAtom = *x;
          char temp[100];
          if (firstAtom) {
            sprintf(temp,"%6d %2s %18.10f%18.10f%18.10f %5s%3s%s%4s%s",
                    counter, (pAtom->getElement()->symbol.c_str()),
                    pAtom->getX(),pAtom->getY(),pAtom->getZ(),
                    "RES {", pAtom->getParent()->getName().c_str(), ";", pAtom->getName().c_str(), "}");
            firstAtom = false;
          }
          else {
            sprintf(temp,"%6d %2s %18.10f%18.10f%18.10f %5s%3s%s%4s%s",
                    counter, (pAtom->getElement()->symbol.c_str()),
                    pAtom->getX(),pAtom->getY(),pAtom->getZ(),
                    "    {", pAtom->getParent()->getName().c_str(), ";", pAtom->getName().c_str(), "}");
          }
          this->inputFileStream << temp << std::endl;
          counter++;
        }
      }
    }
    else {
      Atoms = pMolecule->getAtomList();
      int atomIndex = 0;
      double curX, curY, curZ;
      for (AtomIterator d = Atoms.begin(); d != Atoms.end(); d++) {
        pAtom = *d;
        char temp[100];
      //  if (coordinates.size() == Atoms.size()) {
          curX = coordinates[atomIndex].getX();
          curY = coordinates[atomIndex].getY();
          curZ = coordinates[atomIndex].getZ();
      /*  }
        else {
          curX = pAtom->getX();
          curY = pAtom->getY();
          curZ = pAtom->getZ();
        } */
        sprintf(temp,"%6d %2s %18.10f%18.10f%18.10f %s%3s%s%4s%s",
                counter, (pAtom->getElement()->symbol.c_str()),
                curX, curY, curZ,
                "{","XXX" , ";", (pAtom->getElement()->symbol.c_str()), "}");
                // pAtom->getParent()->getName().c_str(), pAtom->getName().c_str(), "}");
        this->inputFileStream << temp << std::endl;
        counter++;
        atomIndex++;
      }
    }
    this->inputFileStream << "END_COORD" << std::endl;
}

// ============================================================
// Function : writeTail()
// ------------------------------------------------------------
// Write Tail of a DivCon Input File
// ============================================================
void dcParser::writeTail()
{
    if (!this->inputFileStream) return;

    //if (bVdw) odc << "vdw \n" << vdw << "\nend_vdw \n" << std::endl;
    if (bCluster) {
      this->inputFileStream << "cluster" << std::endl;
      this->inputFileStream << "  dbuff1=4.0" << std::endl;
      this->inputFileStream << "  dbuff2=2.0" << std::endl;
      this->inputFileStream << "  ncore=1" << std::endl;
      this->inputFileStream << "end_cluster" << std::endl;
    }
    if (bNMR) {
      this->inputFileStream << "nmr" << std::endl;
      this->inputFileStream << " atom 1-" << pMolecule->getNumAtoms() << std::endl;
      this->inputFileStream << "end_nmr" << std::endl;
    }
    if (bGuess) {
      std::string dmxFile = outputFileName + ".dmx";
      this->inputFileStream << "guess" << std::endl;
      this->inputFileStream << dmxFile  << std::endl;
      this->inputFileStream << "end_guess" << std::endl;
    }
    this->inputFileStream.close();
}


/*
    ofstream odc;
    odc.open(dcfile.c_str());

    if (!odc or (pMolecule == 0)) {
      std::cout << "\nUNABLE TO OPEN DIVCON FILE"
           << "\nFILENAME = " << dcfile << std::endl;
      return;
    }


    odc << hamiltonian << " &" << std::endl;
    odc << "charge = " << pMolecule->getTotalCharge() << " &" << std::endl;
    if (bChargeModel) odc << chargeModel << " &" << endl;

    if (bResidue) odc << "residue &" << endl;
    if (bDirect) odc << "direct &" << endl;
    if (!bStandard) {
      if (bCluster) odc << "cluster & " << endl;
      if (bCutBond) odc << "cutbond=" << cutbond << " &" << endl;
    }
    if (bShift > 0) odc << "shift=" << shift << " &" << endl;

    if (bGuess) odc << "guess &" << endl;
    if (bPDump) odc << "pdump= " << pDump << "&" << endl;
    if (bDump) odc << "dump= " << dump << " &" << endl;

    if (bOpt) {
      odc << "opt=" << optimizer << " &" << endl;
      odc << "maxopt=" << maxopt << " &" << endl;
    }

    if (bScrf) {
      odc << "scrf divpb scale=2.0 &" << endl;
      if (bWater) odc << "water & " << endl;
      if (bOctanol) odc << "octanol &" << endl;
    }

    if (bNoOverlap) odc << "no_overlap &" << endl;

    if (bPwd) odc << "pwd &" << endl;
    if (bPwdAtom) odc << "pwd_atom &" << endl;
    if (bPwdResidue) odc << "pwd_residue &" << endl;

    if (bScreen) odc << "screen &" << endl;
    //if (vdw.size() > 0) odc << "vdw &" << endl;

    //if (iDouble > -1) odc << "double = " << iDouble << " &" << endl;

    if (bStandard) odc << "standard &" << endl;

    if (bDipole) odc << "dipole &" << endl;

    if (bThermo) odc << "thermo &" << endl;

    if (bFreq) odc << "freq &" << endl;
    if (bNMR) odc << "nmr DPSCF=1.0D-9 &" << endl;

    if (bQMAlign) odc << "qmalign &" << endl;

    if (bIntegrals) odc << "intgls=" << integrals << " &" << endl;

    if (calnum) odc << "CALNUC=" << calnum << " & " << endl;

    // Let coord be the last one
    odc << coord << " \n" << endl;

    //////////////////////////////////////////////////////////
    // - Writes Coordinate Section of a DivCon Input File - //
    //////////////////////////////////////////////////////////
    typedef std::vector<atom*>::iterator AtomIterator;
    std::vector<atom*>                   Atoms;
    int counter = 1;

    if (bResidue) {
      typedef std::vector<submolecule*>::iterator sMolIterator;
      std::vector<submolecule*>            submolecules = pMolecule->getSubMoleculeList();

      int counter = 1;
      for (sMolIterator d = submolecules.begin(); d != submolecules.end(); d++) {
        pSubMolecule = *d;
        Atoms = pSubMolecule->getAtomList();
        for (AtomIterator d = Atoms.begin(); d != Atoms.end(); d++) {
          pAtom = *d;
          char temp[100];
          sprintf(temp,"%6d %2s %18.10f%18.10f%18.10f %5s%3s%s%4s%s",
                  counter, (pAtom->getElement()->symbol.c_str()),
                  pAtom->getX(),pAtom->getY(),pAtom->getZ(),
                  "RES {", pAtom->getParent()->getName().c_str(), ";", pAtom->getName().c_str(), "}");
          odc << temp << endl;
          counter++;
        }
      }
    }
    else {
      Atoms = pMolecule->getAtomList();
      int atomIndex = 0;
      double curX, curY, curZ;
      for (AtomIterator d = Atoms.begin(); d != Atoms.end(); d++) {
        pAtom = *d;
        char temp[100];
      //  if (coordinates.size() == Atoms.size()) {
          curX = coordinates[atomIndex].getX();
          curY = coordinates[atomIndex].getY();
          curZ = coordinates[atomIndex].getZ();
      //  }
      //  else {
      //    curX = pAtom->getX();
      //    curY = pAtom->getY();
      //    curZ = pAtom->getZ();
      //  }
        sprintf(temp,"%6d %2s %18.10f%18.10f%18.10f %s%3s%s%4s%s",
                counter, (pAtom->getElement()->symbol.c_str()),
                curX, curY, curZ,
                "{","XXX" , ";", (pAtom->getElement()->symbol.c_str()), "}");
                // pAtom->getParent()->getName().c_str(), pAtom->getName().c_str(), "}");
        odc << temp << endl;
        counter++;
        atomIndex++;
      }
    }
    odc << "END_COORD" << endl;

    ////////////////////////////////////////////
    // - Writes Tail of a DivCon Input File - //
    ////////////////////////////////////////////
    //if (bVdw) odc << "vdw \n" << vdw << "\nend_vdw \n" << endl;
    if (bCluster) {
      odc << "cluster \n" << endl;
      odc << "  dbuff1=4.0" << "\n" << endl;
      odc << "  dbuff2=2.0" << "\n" << endl;
      odc << "  ncore=1" << "\n" << endl;
      odc << "\nend_cluster \n" << endl;
    }
    if (bNMR) {
      odc << "nmr" << endl;
      odc << " atom 1-" << pMolecule->getNumAtoms() << endl;
      odc << "end_nmr" << endl;
    }
*/

} // MTKpp namespace

