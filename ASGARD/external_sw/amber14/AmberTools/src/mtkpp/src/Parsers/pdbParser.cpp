/*!
   \file pdbParser.cpp
   \brief Parses pdb files
   \author Martin Peters
   \author Duane Williams

   Reads and writes pdb files

   $Date: 2010/08/19 13:48:48 $
   $Revision: 1.31 $

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

#include "pdbParser.h"
#include "parsingException.h"

#include "Molecule/element.h"
#include <vector>
#include <algorithm>
#include <sstream>

#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/submolecule.h"
#include "Molecule/atom.h"
#include "Molecule/metalCenter.h"
#include "Utils/vector3d.h"

#include "Molecule/bond.h"

#include "Log/errorHandler.h"

namespace MTKpp
{

// ============================================================
// Function : pdbParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
pdbParser::pdbParser()
{
    this->itsPdbInfo = 0;
    this->itsPdbInfo = new pdbInfo();

    if (this->itsPdbInfo) {
      this->setInfo();
    }
    else {
      errorLogger.throwError("pdbParser", " exiting ", 1);
      //exit(0);
      errorLogger.flush();
      throw parsingException("pdbParser exiting ");
    }

    // Residues 3L -> 1L
    this->res21l["ALA"] = "A";
    this->res21l["ARG"] = "R";
    this->res21l["ASN"] = "N";
    this->res21l["ASP"] = "D";
    this->res21l["CYS"] = "C";
    this->res21l["CYM"] = "C";
    this->res21l["GLU"] = "E";
    this->res21l["GLN"] = "Q";
    this->res21l["GLY"] = "G";
    this->res21l["HIS"] = "H";
    this->res21l["HIE"] = "H";
    this->res21l["HID"] = "H";
    this->res21l["HIP"] = "H";
    this->res21l["HIN"] = "H";
    this->res21l["ILE"] = "I";
    this->res21l["LEU"] = "L";
    this->res21l["LYS"] = "K";
    this->res21l["MET"] = "M";
    this->res21l["PHE"] = "F";
    this->res21l["PRO"] = "P";
    this->res21l["SER"] = "S";
    this->res21l["THR"] = "T";
    this->res21l["TRP"] = "W";
    this->res21l["TYR"] = "Y";
    this->res21l["VAL"] = "V";
    this->res21l["ACE"] = "X";
    this->res21l["NME"] = "X";
    this->res21l[" ZN"] = "X";
    this->res21l["HOH"] = "X";
    this->res21l["WAT"] = "X";

    // Modified CYS Residues
    this->res21l["CCS"] = "C";
    this->res21l["CME"] = "C";

    // Modified LYS Residues
    this->res21l["KCX"] = "K";
}

// =========================================================
// Function : setInfo
// ---------------------------------------------------------
// Set info
// =========================================================
void pdbParser::setInfo()
{
    if (this->itsPdbInfo) {
      delete this->itsPdbInfo;
      this->itsPdbInfo = 0;
    }

    this->itsPdbInfo = new pdbInfo();
    this->itsPdbInfo->resolution = 0.0;
    this->itsPdbInfo->expTechnique = "UNKNOWN";
}

// =========================================================
// Function : get1LCode
// ---------------------------------------------------------
// Get 1-Letter code
// =========================================================
std::string pdbParser::get1LCode(std::string s)
{
    std::string oneLetterCode;

    nameMapIterator p = this->res21l.find(s);

    if (p != this->res21l.end()) {
      return p->second;
    }
    std::string eM = "Can't determine 1-letter code for " + s + ". Using 'X'.";
    errorLogger.throwError("pdbParser", eM, WARNING);
    return "X";
}

// =========================================================
// Function : ~pdbParser()
// ---------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// =========================================================
pdbParser::~pdbParser()
{
    if (this->itsPdbInfo) delete this->itsPdbInfo;
}

// =========================================================
// Function : Read
// ---------------------------------------------------------
// parsers a pdb file
// =========================================================
void pdbParser::Read(const std::string &pdbfile, collection* pCollection)
{
    setError(0);
    std::ifstream ipdb;
    ipdb.open(pdbfile.c_str());

    int end   = pdbfile.length();
    int slash = pdbfile.find_last_of("/");
    std::string file_name = pdbfile.substr(slash+1,(end-slash-5));

    if (!ipdb) {
      setError(1);
      std::string errorMessage = "Error in pdbParser::Read -- can't read " + pdbfile;
      setErrorMessage(errorMessage);
      errorLogger.throwError("pdbParser::Read", errorMessage, 1);
    }

    if (!pCollection) {
      errorLogger.throwError("pdbParser::Read", " pCollection is 0 ... exiting ", 1);
      errorLogger.flush();
      throw parsingException("pdbParser::Read pCollection is 0 ... exiting ");
    }

    this->setInfo();
    if (!this->itsPdbInfo) {
      errorLogger.throwError("pdbParser::Read", " exiting ", 1);
      errorLogger.flush();
      throw parsingException("pdbParser::Read exiting ");
    }

    std::string fileline = "";
    int serial = 0;

    std::string atName = "";
    std::string atSymbol = "";
    std::string atCharge = "";
    std::string resName = "";
    std::string chainID = "";
    std::string iCode = "";
    std::string altLoc = "";
    std::string prevResName = "";
    std::string prevChainID = "";
    std::string prevType = "";

    int resSeq = 0;
    int prevResSeq = 0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double occ = 0.0;
    double tempFac = 0.0;

    bool bATOM = false;
    bool bHETATM = false;
    int bPrevATOM = 0;
    bool bTER = false;
    bool bFIRST = true;

    pMolecule = pCollection->addMolecule();
    pMolecule->setName(file_name);
    pMolecule->setMolId(pCollection->getNumberMolecules());
    pMolecule->inFileType = "pdb";

    hasTotalChargeRemark = false;

    // - Read the file - //
    while (ipdb) {
      getline(ipdb,fileline);
      if (fileline.substr(0,6) == "REMARK") {
        this->itsPdbInfo->remarks.push_back(fileline);
      }

      atName = "";
      atSymbol = "";

      // - CHARGE - //
      if (fileline.substr(0,6) == "REMARK") {
        std::vector<std::string> splitstr;
        splitString(fileline, " ", splitstr, 0);
        if (splitstr.size() > 1) {
          if (splitstr[1].compare("CHARGE") == 0) {
            totalCharge = atof((char*)splitstr[2].c_str());
            hasTotalChargeRemark = true;
          }
        }
      }

      if (fileline.substr(0,22) == "REMARK   2 RESOLUTION." and
          fileline.substr(0,38) != "REMARK   2 RESOLUTION. NOT APPLICABLE.") {
        std::vector<std::string> splitLine;
        splitString(fileline, " ", splitLine, 0);
        this->itsPdbInfo->resolution = string2Double(splitLine[3]);
      }
      if (fileline.substr(0,6) == "EXPDTA") {
        this->itsPdbInfo->expTechnique = stripString(fileline.substr(6, std::min(21, int(fileline.size()))), " ");
        std::string::size_type loc = this->itsPdbInfo->expTechnique.find("X-RAY", 0);
        if ( loc != std::string::npos ) {
          this->itsPdbInfo->expTechnique = "X-RAY";
        }
        loc = this->itsPdbInfo->expTechnique.find("NMR", 0);
        if ( loc != std::string::npos ) {
          this->itsPdbInfo->expTechnique = "NMR";
        }
      }

      bHETATM = false;
      bATOM = false;
      int bCurATOM = 1;

      // 1 -  6  Record name
      if (fileline.substr(0,6) == "HETATM") bHETATM = true;
      if (fileline.substr(0,6) == "ATOM  ") bATOM = true;
      if (bATOM) bCurATOM = 1;
      if (bHETATM) bCurATOM = 2;

      if (fileline == "" || std::string(fileline.substr(0,1).c_str()) == "#") continue;
/*
      if (fileline.substr(0,3) == "END" && bTER) {
        pMolecule = pCollection->addMolecule();
        bTER = false;
      }
*/
      if (fileline.substr(0,3) == "TER") {
        bTER = true;
      }

      if (bATOM || bHETATM) {
        //  7 - 11  AtmNum    Integer        serial     Atom serial number.
        serial  =  atoi(fileline.substr(6,11).c_str());

        // 13 - 16  AtmName   Atom           name       Atom name.
        atName  =  fileline.substr(12,4);

        // 17                 Character      altLoc     Alternate location indicator.
        altLoc = fileline.substr(16,1);

        // 18 - 20  ResName   Residue name   resName    Residue name.
        resName =  fileline.substr(17,3);

        // 22       ChainID   Character      chainID    Chain identifier.
        chainID =  fileline.substr(21,1);

        //23 - 26            Integer        resSeq     Residue sequence number.
        resSeq  =  atoi(fileline.substr(22,4).c_str());

        //27                 AChar          iCode      Code for insertion of residues.
        iCode =  fileline.substr(26,1);

        // 31 - 38            Real(8.3)      x          Orthogonal coordinates for X.
        x       =  atof(fileline.substr(30,8).c_str());

        // 39 - 46            Real(8.3)      y          Orthogonal coordinates for Y.
        y       =  atof(fileline.substr(38,8).c_str());

        // 47 - 54            Real(8.3)      z          Orthogonal coordinates for Z.
        z       =  atof(fileline.substr(46,8).c_str());

        try {
          // 55 - 60            Real(6.2)      occupancy  Occupancy.
          occ     =  atof(fileline.substr(54,6).c_str());

          // 61 - 66            Real(6.2)      tempFacto  Temperature factor.
          tempFac =  atof(fileline.substr(60,6).c_str());

          // 77 - 78            LString(2)     element    Element symbol; right-justified.
          atSymbol = stripString(fileline.substr(76,2), " ");
          if (atSymbol.size() == 2) {
            std::string t = atSymbol.substr(1,1);
            atSymbol = atSymbol[0] + toLower(t);
          }
        }
        catch (...) {
          occ = 0.0;
          tempFac = 0.0;
          atSymbol = "";
        }
        try {
          // 79 - 80            LString(2)     charge    Formal Charge; right-justified.
          atCharge = stripString(fileline.substr(78,2), " ");
          reverse(atCharge.begin(), atCharge.end());
        }
        catch (...) {
          atCharge = "";
        }

        atSymbol = this->determineElement(atName);
        std::string curType = this->determineType(atSymbol, resName, bATOM, bHETATM);

        if (resName == "WAT") resName = "HOH";

        if (!bFIRST) {
          if (bTER or (curType == "Met") or ((bPrevATOM == bCurATOM) and (curType != prevType))) {
            pMolecule = pCollection->addMolecule();
            pMolecule->setMolId(pCollection->getNumberMolecules());
            pMolecule->setName(resName);
            bTER = false;
          }
        }

        if (pMolecule->getNumAtoms() == 0) {
          prevResName = resName;
          prevResSeq  = resSeq;
          prevChainID = chainID;
          prevType = curType;
          if (bATOM) bPrevATOM = 1;
          if (bHETATM) bPrevATOM = 2;

          pSubMolecule = pMolecule->addSubMolecule();
          pSubMolecule->setName(resName);
          pSubMolecule->set1LName(this->get1LCode(resName));
          pSubMolecule->setSubMolId(resSeq);
          pSubMolecule->setiCode(iCode);
        }
        else {
          if (resName == prevResName) {
            if (resSeq != prevResSeq) {
              if (bHETATM) {
                if ((prevChainID != chainID) or (resName == "HOH") ) {
                  pMolecule = pCollection->addMolecule();
                  pMolecule->setMolId(pCollection->getNumberMolecules());
                  pMolecule->setName(resName);
                }
              }

              pSubMolecule = pMolecule->addSubMolecule();
              pSubMolecule->setName(resName);
              pSubMolecule->set1LName(this->get1LCode(resName));
              pSubMolecule->setSubMolId(resSeq);
              pSubMolecule->setiCode(iCode);
              prevResName = resName;
              prevResSeq  = resSeq;
              prevChainID = chainID;
              prevType = curType;
              bPrevATOM = bCurATOM;
            }
            else {
              if (bHETATM) {
                if ((prevChainID != chainID)) {
                  pMolecule = pCollection->addMolecule();
                  pMolecule->setMolId(pCollection->getNumberMolecules());
                  pMolecule->setName(resName);
                  pSubMolecule = pMolecule->addSubMolecule();
                  pSubMolecule->setName(resName);
                  pSubMolecule->set1LName(this->get1LCode(resName));
                  pSubMolecule->setSubMolId(resSeq);
                  pSubMolecule->setiCode(iCode);
                  prevResName = resName;
                  prevResSeq  = resSeq;
                  prevChainID = chainID;
                  prevType = curType;
                  bPrevATOM = bCurATOM;
                }
              }
            }
          }
          else {
            if (bHETATM) {
              if ((prevChainID != chainID)) {
                pMolecule = pCollection->addMolecule();
                pMolecule->setMolId(pCollection->getNumberMolecules());
                pMolecule->setName(resName);
              }
            }
            if (bATOM) {
              if (prevChainID != chainID) {
                pMolecule = pCollection->addMolecule();
                pMolecule->setMolId(pCollection->getNumberMolecules());
                pMolecule->setName(resName);
              }
            }
            pSubMolecule = pMolecule->addSubMolecule();
            pSubMolecule->setName(resName);
            pSubMolecule->set1LName(this->get1LCode(resName));
            pSubMolecule->setSubMolId(resSeq);
            pSubMolecule->setiCode(iCode);
            prevResName = resName;
            prevResSeq  = resSeq;
            prevChainID = chainID;
            prevType = curType;
            bPrevATOM = bCurATOM;
          }
        }

        if (!pSubMolecule) {
          errorLogger.throwError("pdbParser::Read", " exiting ", 1);
          //exit(1);
          errorLogger.flush();
          throw parsingException("pdbParser::Read exiting ");
        }

        if (!pSubMolecule->getAtom(atName)) {
          bFIRST = false;
          if (bATOM) pMolecule->setKind(1); // biomolecule
          if (bHETATM) pMolecule->setKind(2); // drug/small molecule
          pAtom = pSubMolecule->addAtom();
          if (atSymbol == "") {
            atSymbol = this->determineElement(atName);
            if (!pCollection->pElements->getElement(atSymbol)) {
              std::string errMes = " Can't assign element for the atomName = " + atName + " ... exiting ";
              errorLogger.throwError("pdbParser::Read", errMes, 1);
              //exit(1);
              errorLogger.flush();
              throw parsingException("pdbParser::Read Can't assign element for the atomName = " + atName + " ... exiting ");
            }

            if (atCharge.size() > 0) {
              pAtom->setFormalCharge(atoi(atCharge.c_str()));
            }
            stringstream ss;
            ss <<" "<<this->determineElement(atName)<<" "<<x<<" "<<y<<" "<<z;
            errorLogger.throwError("pdbParser::Read",ss.str(), 1);

            element* ele = pCollection->pElements->getElement(this->determineElement(atName));
            if (ele) {
              pAtom->setElement(ele);
            }
            else {
              std::string errMes = " Can't find element for " + atName;
              MTKpp::errorLogger.throwError("pdbParser::read", errMes, 1);
              //exit(1);
              errorLogger.flush();
              throw parsingException("pdbParser::read  Can't find element for " + atSymbol);
            }
          }
          else {
            element* ele = pCollection->pElements->getElement(atSymbol);
            if (ele) {
              pAtom->setElement(ele);
            }
            else {
              std::string errMes = " Can't find element for " + atSymbol;
              MTKpp::errorLogger.throwError("pdbParser::read", errMes, 1);
              exit(1);
            }
          }
          pAtom->setName(atName);
          pAtom->setCoords(x,y,z);
          pAtom->setFileID(serial);
          pAtom->setOccupancy(occ);
          pAtom->setTempFactor(tempFac);
          pMolecule->setChain(chainID);
          pMolecule->setMaxFileID(serial);
        }
        else {
          std::string errMes = pSubMolecule->getName() + int2String(pSubMolecule->getSubMolId()) + " already contains " + atName;
          MTKpp::errorLogger.throwError("pdbParser::read", errMes, 4);
        }
      }

/*
GET THIS TO WORK
    int at1_serial = 0;
    int at2_serial = 0;
      if (fileline.substr(0,6) == "CONECT") {
        at1_serial  =  atoi(fileline.substr(6,11).c_str());
        at2_serial  =  atoi(fileline.substr(12,17).c_str());
        pAtom1 = pMolecule->getAtom(at1_serial, false, true);  // WRONG
        pAtom2 = pMolecule->getAtom(at1_serial, false, true);  // WRONG
        if (!pAtom1 or !pAtom2) continue;
        pCoord1 = pAtom1->getCoords();
        pCoord2 = pAtom2->getCoords();
        double distance = pCoord1->dist(*pCoord2);
        // TO-DO: FIGURE OUT CAN I SET BOND TYPES, STEREO, AND TOPOLOGY
        pMolecule->addBond(pAtom1, pAtom2, 0, 0, 0, distance);
        pAtom1->addBondedAtom(pAtom2);
        pAtom2->addBondedAtom(pAtom1);
      }
*/
    } // while (ipdb)
    ipdb.close();
}

// =========================================================
// Function : Read
// ---------------------------------------------------------
// parsers a pdb file
// =========================================================
void pdbParser::Read(const std::string &pdbfile, collection* pCollection, molecule* pMolecule)
{
    std::ifstream ipdb;
    ipdb.open(pdbfile.c_str());

    int end   = pdbfile.length();
    int slash = pdbfile.find_last_of("/");
    std::string file_name = pdbfile.substr(slash+1,(end-slash-5));

    if (!ipdb) {
      std::cout << "\nUNABLE TO OPEN PDB FILE"
           << "\nFILENAME = " << pdbfile
           << "\nEXITING...\n" << std::endl;
      return;
    }

    std::string fileline;
    int serial = 0;

    std::string atName = "";
    std::string resName = "";
    std::string chainID = "";
    std::string iCode = "";
    std::string altLoc = "";
    std::string atCharge = "";
    std::string prevResName;
    int resSeq = 0;
    int prevResSeq = 0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    bool bATOM = false;
    bool bHETATM = false;
    bool bTER = false;

    pMolecule->setName(file_name);
    pMolecule->setMolId(pCollection->getNumberMolecules());
    pMolecule->inFileType = "pdb";

    hasTotalChargeRemark = false;
    // - Read the file - //
    while (ipdb) {
      getline(ipdb,fileline);
      bHETATM = false;
      bATOM = false;

      // - CHARGE - //
      if (fileline.substr(0,6) == "REMARK") {
        std::vector<std::string> splitstr;
        splitString(fileline," ", splitstr, 0);
        if (splitstr.size() > 0) {
          if (splitstr[1].compare("CHARGE") == 0) {
            totalCharge = atof((char*)splitstr[2].c_str());
            hasTotalChargeRemark = true;
          }
        }
      }

      if (fileline.substr(0,6) == "HETATM") bHETATM = true;
      if (fileline.substr(0,6) == "ATOM  ") bATOM = true;

      if (fileline == "" || std::string(fileline.substr(0,1).c_str()) == "#") continue;

      if (fileline.substr(0,3) == "TER") {
        bTER = true;
      }

      if (bATOM || bHETATM) {
        //  7 - 11  AtmNum    Integer        serial     Atom serial number.
        serial  =  atoi(fileline.substr(6,11).c_str());

        // 13 - 16  AtmName   Atom           name       Atom name.
        atName  =  fileline.substr(12,4);

        // 17                 Character      altLoc     Alternate location indicator.
        altLoc = fileline.substr(16,1);

        // 18 - 20  ResName   Residue name   resName    Residue name.
        resName =  fileline.substr(17,3);

        // 22       ChainID   Character      chainID    Chain identifier.
        chainID =  fileline.substr(21,1);

        //23 - 26            Integer        resSeq     Residue sequence number.
        resSeq  =  atoi(fileline.substr(22,4).c_str());

        //27                 AChar          iCode      Code for insertion of residues.
        iCode =  fileline.substr(26,1);

        // 31 - 38            Real(8.3)      x          Orthogonal coordinates for X.
        x       =  atof(fileline.substr(30,8).c_str());

        // 39 - 46            Real(8.3)      y          Orthogonal coordinates for Y.
        y       =  atof(fileline.substr(38,8).c_str());

        // 47 - 54            Real(8.3)      z          Orthogonal coordinates for Z.
        z       =  atof(fileline.substr(46,8).c_str());

        try {
          // 79 - 80            LString(2)     charge    Formal Charge; right-justified.
          atCharge = stripString(fileline.substr(78,2), " ");
          reverse(atCharge.begin(), atCharge.end());
        }
        catch (...) {
          atCharge = "";
        }

        if (bTER) {
          pMolecule = pCollection->addMolecule();
          bTER = false;
        }

        if (pMolecule->getNumAtoms() == 0) {
          prevResName = resName;
          prevResSeq  = resSeq;
          pSubMolecule = pMolecule->addSubMolecule();
          pSubMolecule->setName(resName);
          pSubMolecule->set1LName(this->get1LCode(resName));
          pSubMolecule->setSubMolId(resSeq);
          pSubMolecule->setiCode(iCode);
        }

        if (resName == prevResName) {
          if (resSeq != prevResSeq) {
            pSubMolecule = pMolecule->addSubMolecule();
            pSubMolecule->setName(resName);
            pSubMolecule->set1LName(this->get1LCode(resName));
            pSubMolecule->setSubMolId(resSeq);
            pSubMolecule->setiCode(iCode);
            prevResName = resName;
            prevResSeq  = resSeq;
          }
        }
        else {
          pSubMolecule = pMolecule->addSubMolecule();
          pSubMolecule->setName(resName);
          pSubMolecule->set1LName(this->get1LCode(resName));
          pSubMolecule->setSubMolId(resSeq);
          pSubMolecule->setiCode(iCode);
          prevResName = resName;
          prevResSeq  = resSeq;
        }
        pAtom = pSubMolecule->addAtom();
        pAtom->setElement(pCollection->pElements->getElement(this->determineElement(atName)));
        pAtom->setName(atName);
        pAtom->setCoords(x,y,z);
        if (atCharge.size() > 0) {
          pAtom->setFormalCharge(atoi(atCharge.c_str()));
        }
        pAtom->setFileID(serial);
        pMolecule->setChain(chainID);
        }
/*
    int at1_serial = 0;
    int at2_serial = 0;
      if (fileline.substr(0,6) == "CONECT") {
        at1_serial  =  atoi(fileline.substr(6,11).c_str());
        at2_serial  =  atoi(fileline.substr(12,17).c_str());
        pAtom1 = pMolecule->getAtom(at1_serial, false, true);  // WRONG
        pAtom2 = pMolecule->getAtom(at1_serial, false, true);  // WRONG
        pCoord1 = pAtom1->getCoords();
        pCoord2 = pAtom2->getCoords();
        double distance = pCoord1->dist(*pCoord2);
        // TO-DO: FIGURE OUT CAN I SET BOND TYPES, STEREO, AND TOPOLOGY
        pMolecule->addBond(pAtom1, pAtom2, 0, 0, 0, distance);
        pAtom1->addBondedAtom(pAtom2);
        pAtom2->addBondedAtom(pAtom1);
      }
*/
    } // while (ipdb)
    ipdb.close();
}

// =========================================================
// Function : updateMolCoords
// ---------------------------------------------------------
// parsers a pdb file
// =========================================================
void pdbParser::updateMolCoords(const std::string &pdbfile, molecule* pMolecule)
{
    std::ifstream ipdb;
    ipdb.open(pdbfile.c_str());

    //int end   = pdbfile.length();
    //int slash = pdbfile.find_last_of("/");
    //std::string file_name = pdbfile.substr(slash+1,(end-slash-5));

    if (!ipdb) {
      std::cout << "\nUNABLE TO OPEN PDB FILE"
           << "\nFILENAME = " << pdbfile
           << "\nEXITING...\n" << std::endl;
      return;
    }

    std::string fileline;

    std::string atName = "";
    std::string resName = "";
    std::string prevResName;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    bool bATOM = false;
    bool bHETATM = false;
    //bool bTER = false;

    // - Read the file - //
    while (ipdb) {
      getline(ipdb,fileline);
      bHETATM = false;
      bATOM = false;
      if (fileline.substr(0,6) == "HETATM") bHETATM = true;
      if (fileline.substr(0,6) == "ATOM  ") bATOM = true;

      if (fileline == "" || std::string(fileline.substr(0,1).c_str()) == "#") continue;

      //if (fileline.substr(0,3) == "TER") {
      //  bTER = true;
      //}

      if (bATOM || bHETATM) {
        atName  =  fileline.substr(12,4);
        x       =  atof(fileline.substr(30,8).c_str());
        y       =  atof(fileline.substr(38,8).c_str());
        z       =  atof(fileline.substr(46,8).c_str());

        pAtom = pMolecule->getAtom(atName);

        if (pAtom) {
          pAtom->setCoords(x,y,z);
        }
        else {
          std::cout << " Error in updateMolCoords ... exiting " << std::endl;
          //exit(0);
          errorLogger.flush();
          throw parsingException(" Error in updateMolCoords ... exiting ");
        }
      }
    } // while (ipdb)
    ipdb.close();
}

// ============================================================================
// Function : Write
// ----------------------------------------------------------------------------
// writes a pdb file
// -----------------------------------------------------------------------------
/*
Format:
COLUMNS  Variable  DATA TYPE      FIELD      DEFINITION                    done
--------------------------------------------------------------------------------
 1 -  6  Rerord    Record name    "HETATM"
 7 - 11  AtmNum    Integer        serial     Atom serial number.             *
13 - 16  AtmName   Atom           name       Atom name.                      *
17                 Character      altLoc     Alternate location indicator.
18 - 20  ResName   Residue name   resName    Residue name.                   *
22       ChainID   Character      chainID    Chain identifier.
23 - 26            Integer        resSeq     Residue sequence number.        *
27                 AChar          iCode      Code for insertion of residues.
31 - 38            Real(8.3)      x          Orthogonal coordinates for X.   *
39 - 46            Real(8.3)      y          Orthogonal coordinates for Y.   *
47 - 54            Real(8.3)      z          Orthogonal coordinates for Z.   *
55 - 60            Real(6.2)      occupancy  Occupancy.
61 - 66            Real(6.2)      tempFacto  Temperature factor.
73 - 76            LString(4)     segID      Segment identifier;
                                             left-justified.
77 - 78            LString(2)     element    Element symbol; right-justified.
79 - 80            LString(2)     charge     Charge on the atom.
=================================================================================
*/
void pdbParser::Write(const std::string &pdbfile, molecule* pMolecule)
{
#ifdef DEBUG
    std::cout << " pdbParser::Write " << std::endl;
#endif

    std::ofstream opdb;
    opdb.open(pdbfile.c_str());

    if (!opdb or (pMolecule == 0)) {
      std::cout << "\nUNABLE TO OPEN PDB FILE"
                << "\nFILENAME = " << pdbfile << std::endl;
      return;
    }

    int n = 1;
    opdb << "REMARK" << std::endl;

    std::vector<atom*> atomList = pMolecule->getAtomList();

    // - Loop over the Atoms - //
    for (atomIterator a = atomList.begin(); a != atomList.end(); a++) {
      pAtom = *a;
      pSubMolecule = pAtom->getParent();
      //std::cout << pSubMolecule->getColIndex() << std::endl;
      char temp[100];

      char tmp_element[3];
      if (pAtom->getName() == "") {
        if (pAtom->getElement()->symbol != "") {
          sprintf(tmp_element," %s",pAtom->getElement()->symbol.c_str());
        }
        else {
          if (pAtom->getElement()->symbol == "") {
            //sprintf(tmp_element," %s",pAtom->getElement()->symbol.c_str());
          }
        }
        pAtom->setName(tmp_element);
      }

      sprintf(temp,"ATOM  %5d %-4.4s %-3.3s %1s%4d%1s   %8.3f%8.3f%8.3f",
              pAtom->getFileID(),(pAtom->getName().c_str()),
              (pSubMolecule->getName().c_str()),
              (pMolecule->getChain().c_str()),
              (pSubMolecule->getSubMolId()),
              (pSubMolecule->getiCode().c_str()),
              pAtom->getX(),pAtom->getY(),pAtom->getZ());

      opdb << temp << std::endl;

      n=n+1;
    }
    opdb.close();
}

// ============================================================================
// Function : Write
// ----------------------------------------------------------------------------
// Writes a pdb file
// ============================================================================
void pdbParser::Write(const std::string &pdbfile, std::vector<molecule*> molList)
{
    std::ofstream opdb;
    opdb.open(pdbfile.c_str());

    if (!opdb) {
      std::cout << "\nUNABLE TO OPEN PDB FILE"
                << "\nFILENAME = " << pdbfile << std::endl;
      return;
    }

    int n = 1;
    opdb << "REMARK" << std::endl;

    for (unsigned int m = 0; m < molList.size(); m++) {
      pMolecule = molList[m];
      int molKind = pMolecule->getKind();

      std::vector<atom*> atomList = pMolecule->getAtomList();

      // - Loop over the Atoms - //
      for (atomIterator a = atomList.begin(); a != atomList.end(); a++) {
        pAtom = *a;
        pSubMolecule = pAtom->getParent();
        char temp[100];

        char tmp_element[3];
        if (pAtom->getName() == "") {
          if (pAtom->getElement()->symbol != "") {
            sprintf(tmp_element," %s",pAtom->getElement()->symbol.c_str());
          }
          else {
            if (pAtom->getElement()->symbol == "") {
              //sprintf(tmp_element," %s",pAtom->getElement()->symbol.c_str());
            }
          }
          pAtom->setName(tmp_element);
        }

        std::string atT = "";
        if (molKind == 1) {
          atT = "ATOM  ";
        }
        else {
          atT = "HETATM";
        }

        sprintf(temp,"%-6.6s%5d %-4.4s %-3.3s %1s%4d%1s   %8.3f%8.3f%8.3f",
                (atT.c_str()),
                pAtom->getFileID(),(pAtom->getName().c_str()),
                (pSubMolecule->getName().c_str()),
                (pMolecule->getChain().c_str()),
                (pSubMolecule->getSubMolId()),
                (pSubMolecule->getiCode().c_str()),
                pAtom->getX(),pAtom->getY(),pAtom->getZ());

        opdb << temp << std::endl;

        n=n+1;
      }
      // might need to add ter statements here
    }
    opdb.close();
}

// ============================================================================
// Function : Write
// ----------------------------------------------------------------------------
// Writes a pdb file
// ============================================================================
void pdbParser::Write(const std::string &pdbfile, molecule* pMolecule,
                      std::vector< vector3d > &coordinates)
{
    std::ofstream opdb;
    opdb.open(pdbfile.c_str());

    if (!opdb or (pMolecule == 0)) {
      std::cout << "\nUNABLE TO OPEN PDB FILE"
                << "\nFILENAME = " << pdbfile << std::endl;
      return;
    }

    int n = 1;
    std::vector<atom*> atomList = pMolecule->getAtomList();

    double curX = 0.0;
    double curY = 0.0;
    double curZ = 0.0;

    // - Loop over the Atoms - //
    for (atomIterator a = atomList.begin(); a != atomList.end(); a++) {
      pAtom = *a;
      pSubMolecule = pAtom->getParent();
      char temp[100];

      curX = coordinates[n-1].getX();
      curY = coordinates[n-1].getY();
      curZ = coordinates[n-1].getZ();

      char tmp_element[3];
      if (pAtom->getName() == "") {
        if (pAtom->getElement()->symbol != "") {
          sprintf(tmp_element," %s",pAtom->getElement()->symbol.c_str());
        }
        else {
          if (pAtom->getElement()->symbol == "") {
          }
        }
        pAtom->setName(tmp_element);
      }

      sprintf(temp,"ATOM  %5d %-4.4s %-3.3s %1s%4d%1s   %8.3f%8.3f%8.3f",
              pAtom->getFileID(),(pAtom->getName().c_str()),
              (pSubMolecule->getName().c_str()),
              (pMolecule->getChain().c_str()),
              (pSubMolecule->getSubMolId()),
              (pSubMolecule->getiCode().c_str()),
              curX, curY,curZ);

      opdb << temp << std::endl;

      n++;
    }
    opdb.close();
}

// ============================================================================
// Function : Write
// ----------------------------------------------------------------------------
// Writes a pdb file
// ============================================================================
void pdbParser::Write(const std::string &pdbfile, metalCenter* pMetCen)
{
    std::ofstream opdb;
    opdb.open(pdbfile.c_str());

    if (!opdb or (pMetCen == 0)) {
      std::cout << "\nUNABLE TO OPEN PDB FILE"
                << "\nFILENAME = " << pdbfile << std::endl;
      return;
    }

    int n = 1;

    opdb << "REMARK METAL CENTER" << std::endl;

    std::vector<submolecule*> shellSubMols;
    std::vector<submolecule*>::iterator subMolIterator;

    std::vector<atom*> pAts;
    pMetCen->getPrimaryShellAtoms(pAts);

    std::vector<atom*> sAts;
    pMetCen->getSecondaryShellAtoms(sAts);
/*
    for (unsigned int i = 0; i < this->itsBondedAtoms.size(); i++) {
      if (coordinationType[i] == "p" or coordinationType[i] == "s") {
        submolecule* pSM = this->itsBondedAtoms[i]->getParent();
        subMolIterator = std::find(shellSubMols.begin(), shellSubMols.end(), pSM);
        if (subMolIterator == shellSubMols.end()) {
          shellSubMols.push_back(pSM);
        }
      }
    }
*/
    for (unsigned int i = 0; i < pAts.size(); i++) {
      submolecule* pSM = pAts[i]->getParent();
      subMolIterator = std::find(shellSubMols.begin(), shellSubMols.end(), pSM);
      if (subMolIterator == shellSubMols.end()) {
        shellSubMols.push_back(pSM);
      }
    }

    for (unsigned int i = 0; i < sAts.size(); i++) {
      submolecule* pSM = sAts[i]->getParent();
      subMolIterator = std::find(shellSubMols.begin(), shellSubMols.end(), pSM);
      if (subMolIterator == shellSubMols.end()) {
        shellSubMols.push_back(pSM);
      }
    }

    // Add metal ion
    shellSubMols.push_back(pMetCen->getMetalAtom()->getParent());

    for (unsigned int s = 0; s < shellSubMols.size(); s++) {
      std::vector<atom*> atomList = shellSubMols[s]->getAtomList();

      for (atomIterator a = atomList.begin(); a != atomList.end(); a++) {
        pAtom = *a;
        pSubMolecule = pAtom->getParent();
        char temp[100];

        char tmp_element[3];
        if (pAtom->getName() == "") {
          if (pAtom->getElement()->symbol != "") {
            sprintf(tmp_element," %s",pAtom->getElement()->symbol.c_str());
          }
          else {
            if (pAtom->getElement()->symbol == "") {
              //sprintf(tmp_element," %s",pAtom->getElement()->symbol.c_str());
            }
          }
          pAtom->setName(tmp_element);
        }

        sprintf(temp,"ATOM  %5d %-4.4s %-3.3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%10s%2s",
                pAtom->getFileID(),(pAtom->getName().c_str()),
                (pSubMolecule->getName().c_str()),
                (pMolecule->getChain().c_str()),
                (pSubMolecule->getSubMolId()),
                (pSubMolecule->getiCode().c_str()),
                pAtom->getX(),pAtom->getY(),pAtom->getZ(),
                0.0, // Occupancy
                0.0, // Temperature factor
                " ",
                (pAtom->getElement()->symbol.c_str()));
        opdb << temp << std::endl;
        n=n+1;
      }
    }
/*
    molecule* pM = new molecule(pCollection);
    for (unsigned int j = 0; j < shellSubMols.size(); j++) {
      submolecule* pSM = pM->addSubMolecule();
      pSM->setSubMolId(j+1);
      pSM->setName(shellSubMols[j]->getName());
      std::vector<atom*> smAtoms = shellSubMols[j]->getAtomList();
      for (unsigned int k = 0; k < smAtoms.size(); k++) {
        atom* pA = pSM->addAtom(smAtoms[k]);
        if (!pA) {
          std::cout << " Error in pdbParser::Write " << std::endl;
        }
      }
    }
    this->Write(pdbfile, pM);
*/
    opdb.close();
}

// ============================================================================
// Function : Write
// ----------------------------------------------------------------------------
// Writes a pdb file
// ============================================================================
void pdbParser::Write(const std::string &pdbfile, collection* pCollection)
{
    std::ofstream opdb;
    opdb.open(pdbfile.c_str());

    if (!opdb or (pCollection == 0)) {
      std::cout << "\nUNABLE TO OPEN PDB FILE"
                << "\nFILENAME = " << pdbfile << std::endl;
      return;
    }

    int n = 1;

    opdb << "REMARK" << std::endl;
    for (unsigned int i = 0; i < this->itsPdbInfo->remarks.size(); i++) {
      opdb << this->itsPdbInfo->remarks[i] << std::endl;
    }

    std::vector<molecule*> molList = pCollection->getMoleculeList();
    for (moleculeIterator m = molList.begin(); m != molList.end(); m++) {
      pMolecule = *m;
      if (pMolecule->getName() == "Reference") continue;
      std::vector<atom*> atomList = pMolecule->getAtomList();

      for (atomIterator a = atomList.begin(); a != atomList.end(); a++) {
        pAtom = *a;
        pSubMolecule = pAtom->getParent();
        char temp[100];

        char tmp_element[3];
        if (pAtom->getName() == "") {
          if (pAtom->getElement()->symbol != "") {
            sprintf(tmp_element," %s",pAtom->getElement()->symbol.c_str());
          }
          else {
            if (pAtom->getElement()->symbol == "") {
              //sprintf(tmp_element," %s",pAtom->getElement()->symbol.c_str());
            }
          }
          pAtom->setName(tmp_element);
        }

        //sprintf(temp,"ATOM  %5d %-4.4s %-3.3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%11s%-2s",
        sprintf(temp,"ATOM  %5d %-4.4s %-3.3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f%10s%2s",   //Changed by ZXH.
                pAtom->getFileID(),(pAtom->getName().c_str()),
                (pSubMolecule->getName().c_str()),
                (pMolecule->getChain().c_str()),
                (pSubMolecule->getSubMolId()),
                (pSubMolecule->getiCode().c_str()),
                pAtom->getX(),pAtom->getY(),pAtom->getZ(),
                0.0, // Occupancy
                0.0, // Temperature factor
                " ",
                (pAtom->getElement()->symbol.c_str()));
        opdb << temp << std::endl;
        n=n+1;
      }
      opdb << "TER" << std::endl;
    }
    opdb.close();
}

// ============================================================================
// Function : Write
// ----------------------------------------------------------------------------
// Writes a pdb file
// ============================================================================
void pdbParser::Write(std::ostream& opdb, collection* pCollection)
{
    opdb << "MODEL" << std::endl;
    std::vector<molecule*> molList = pCollection->getMoleculeList();
    for (moleculeIterator m = molList.begin(); m != molList.end(); m++) {
      pMolecule = *m;
      std::vector<atom*> atomList = pMolecule->getAtomList();

      for (atomIterator a = atomList.begin(); a != atomList.end(); a++) {
        pAtom = *a;
        pSubMolecule = pAtom->getParent();
        char temp[100];

        char tmp_element[3];
        if (pAtom->getName() == "") {
          if (pAtom->getElement()->symbol != "") {
            sprintf(tmp_element," %s",pAtom->getElement()->symbol.c_str());
          }
          else {
            if (pAtom->getElement()->symbol == "") {
              //sprintf(tmp_element," %s",pAtom->getElement()->symbol.c_str());
            }
          }
          pAtom->setName(tmp_element);
        }
        sprintf(temp,"ATOM  %5d %-4.4s %-3.3s %1s%4d%1s   %8.3f%8.3f%8.3f",
                pAtom->getFileID(),(pAtom->getName().c_str()),
                (pSubMolecule->getName().c_str()),
                (pMolecule->getChain().c_str()),
                (pSubMolecule->getSubMolId()),
                (pSubMolecule->getiCode().c_str()),
                pAtom->getX(),pAtom->getY(),pAtom->getZ());

        opdb << temp << std::endl;
      }
      opdb << "TER" << std::endl;
    }
}

// ============================================================================
// Function : determineType
// ----------------------------------------------------------------------------
//
// ============================================================================
std::string pdbParser::determineType(std::string a, std::string r, bool at, bool b)
{
    if (at) {
      return "Pro";
    }
    if (b) {
      if ((r == "WAT") or (r == "HOH")) {
        return "Wat";
      }
      else if ((a == "Zn") or (a == "Ca")) {
        return "Met";
      }
      else {
        return "Lig";
      }
    }
    return "Pro";
}

} // MTKpp namespace

