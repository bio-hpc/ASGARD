/*!
   \file mol2Parser.cpp
   \brief Parsers mol2 files
   \author Martin Peters

   $Date: 2010/07/22 10:42:44 $
   $Revision: 1.12 $

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
#include <string.h>
#include <vector>

#include "mol2Parser.h"

#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/submolecule.h"
#include "Molecule/atom.h"
#include "Molecule/element.h"
#include "Molecule/connections.h"

#include "Molecule/bond.h"
#include "Molecule/angle.h"
#include "Molecule/torsion.h"
#include "Molecule/improper.h"
#include "Utils/vector3d.h"
#include "Utils/constants.h"

#include "StringManip.h"

//#include <boost/regex.hpp>

namespace MTKpp
{

// ============================================================
// Function : Mol2Parser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
mol2Parser::mol2Parser():baseParser() {
    // Amino Acids
    this->res21l["ALA"] = "A";
    this->res21l["ARG"] = "R";
    this->res21l["ASN"] = "N";
    this->res21l["ASP"] = "D";
    this->res21l["CYS"] = "C";
    this->res21l["GLU"] = "E";
    this->res21l["GLN"] = "Q";
    this->res21l["GLY"] = "G";
    this->res21l["HIS"] = "H";
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
}

// ============================================================
// Function : ~Mol2Parser()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
mol2Parser::~mol2Parser() {}

// =========================================================
// Function : Read
// ---------------------------------------------------------
// parsers a mol2 file
// ---------------------------------------------------------
// Format:
//
//@<TRIPOS>MOLECULE
// -> mol name
// -> num_atoms [num_bons [num_sust [num_feat [num_sets]]]]
// -> mol_type
// -> charge_type
// -> [status_bits
// -> [mol_comments]]
//
// @<TRIPOS>ATOM
// -> atom_id atom_name x y z atom_type [sust_id [subst_name [charge [ status_bit]]]]
//
//@<TRIPOS>BOND
// -> bond_id origin_atom_id target_atom_id bond_type [status_bits]
//
// =========================================================
void mol2Parser::Read(const std::string &mol2file, molecule* pMol)
{
    std::cout << mol2file << " mol2Parser::Read" << std::endl;
    this->pMol=pMol;
    collection* pCollection = pMol->getParent();
    connections* pConnections = new connections(pCollection);
    int end   = mol2file.length();
    int slash = mol2file.find_last_of("/");
    std::string file_name = mol2file.substr(slash+1,(end-slash-6));

    std::ifstream imol2;
    imol2.open(mol2file.c_str());

    if (!imol2) {
      std::cout << "\nUNABLE TO OPEN MOL2 FILE"
                << "\nFILENAME = " << mol2file << std::endl;
      return;
    }

    if (!pCollection) {
            std::cout << "\nUNABLE TO get pCollection"
                      << "\nFILENAME = " << pCollection << std::endl;
      return;
    }

    std::string fileline;
    int natoms = 0, nbonds = 0, at1 = 0, at2 = 0;
    int bondtype = 0;
    int bondStereo = 0;
    int bondTopology = 0;

    //AtomLine* pthisAtom = new AtomLine;
    chargeType = "";
    hasTotalChargeRemark = false;

    pMol->setName(file_name);
    pMol->setMolId(pCollection->getNumberMolecules());
    pMol->inFileType = "mol2";

    //pSmol = pMol->addSubMolecule();

    // - Read the file - //
    while (imol2) {
      std::string buffer(80,'*');
      getline(imol2,fileline);

      // - CHARGE TYPE - //
      if (fileline.substr(0,11) == "NO_CHARGES") {
            chargeType="NO_CHARGES";
      }
      if (fileline.substr(0,6) == "DEL_RE") {
            chargeType="DEL_RE";
      }
      if (fileline.substr(0,9) == "GASTEIGER") {
            chargeType="GASTEIGER";
      }
      if (fileline.substr(0,9) == "GAST_HUCK") {
            chargeType="GAST_HUCK";
      }
      if (fileline.substr(0,6) == "HUCKEL") {
            chargeType="HUCKEL";
      }
      if (fileline.substr(0,7) == "PULLMAN") {
            chargeType="PULLMAN";
      }
      if (fileline.substr(0,15) == "GAUSS80_CHARGES") {
            chargeType="GAUSS80_CHARGES";
      }
      if (fileline.substr(0,13) == "AMPAC_CHARGES") {
            chargeType="AMPAC_CHARGES";
      }
      if (fileline.substr(0,16) == "MULLIKEN_CHARGES") {
            chargeType="MULLIKEN_CHARGES";
      }
      if (fileline.substr(0,12) == "DICT_CHARGES") {
            chargeType="DICT_CHARGES";
      }
      if (fileline.substr(0,14) == "MMFF94_CHARGES") {
            chargeType="MMFF94_CHARGES";
      }
      if (fileline.substr(0,12) == "USER_CHARGES") {
            chargeType="USER_CHARGES";
      }

      // - CHARGE - //
      if (fileline.substr(0,6) == "REMARK") {
        std::vector<std::string> splitstr;
        splitString(fileline," ", splitstr, 0);
        if(splitstr[1].compare("CHARGE")==0)
        {
            totalCharge = atof((char*)splitstr[2].c_str());
            hasTotalChargeRemark = true;
        }
      }

      // - MOLECULE - //
      if (fileline.substr(0,17) == "@<TRIPOS>MOLECULE") {
        getline(imol2,fileline);
        getline(imol2,fileline);

        std::vector<std::string> splitstr;
        splitString(fileline," ", splitstr, 0);

        natoms = atoi((char*)splitstr[0].c_str());
        nbonds = atoi((char*)splitstr[1].c_str());
      }
      // - ATOM - //
      if (fileline.substr(0,13) == "@<TRIPOS>ATOM") {
        std::string residueID="";
        std::vector<AtomLine*> atomLineVector;
        int n = 1;
        for (n = 1; n <= natoms; n ++) {
          getline(imol2,fileline);
          AtomLine* pthisAtom = new AtomLine;

          if (ReadAtomLine(fileline,pthisAtom)) {
            if(n==1)
            {
                residueID=pthisAtom->segID;
            }
            if(residueID.compare(pthisAtom->segID)!=0)
            {
                buildupSubmolecule(pCollection, residueID, n, atomLineVector);
            }
            atomLineVector.push_back(pthisAtom);
            residueID=pthisAtom->segID;
          }
          else {
                buildupSubmolecule(pCollection, residueID, n, atomLineVector);
            return;
          }
/*
          if (ReadAtomLine(fileline,pthisAtom)) {

            pAtom = pSmol->addAtom();
            pAtom->setElement(pCollection->pElements->getElement(pthisAtom->element));

            std::stringstream ss;
            ss << n;
            std::string at_name_n = ss.str().c_str();
            std::string temp_n_name = (pthisAtom->element)+at_name_n;
            pAtom->setName(temp_n_name);

            pAtom->setCoords(pthisAtom->x,pthisAtom->y,pthisAtom->z);
            // pAtom->setNumBonds(0);
            pAtom->setFileID(n);
            pMol->setMaxFileID(n);
          }
          else {
            return;
          }
*/
        }
        buildupSubmolecule(pCollection, residueID, n-1, atomLineVector);
      }

      // - BOND - //
      //      1      1      2     ar
      if (fileline.substr(0,13) == "@<TRIPOS>BOND") {
        for (int n = 1; n <= nbonds; n ++) {
          getline(imol2,fileline);
          if (fileline.length() < 4) {
            std::cout << "NOT A VALID MOL2 FILE\nThe Bond section"
                      << " should contain:\nbond_id origin_atom_id"
                      << " target_atom_id bond_type [status_bits]"
                      << std::endl;
            return;
          }
          else {
            std::vector<std::string> split_bond_str;
            splitString(fileline," ", split_bond_str, 0);

            at1 = atoi((char*)split_bond_str[1].c_str());
            at2 = atoi((char*)split_bond_str[2].c_str());
            if (split_bond_str[3] == "ar") {
              bondtype = 4;
            }
            else {
              bondtype = atoi((char*)split_bond_str[3].c_str());
            }

            pBondAtom1 = pSmol->getAtom(at1);
            pBondAtom2 = pSmol->getAtom(at2);

            if (pBondAtom1 and pBondAtom2) {
              pBond = 0;
              pBond = pMol->addBond(pBondAtom1,pBondAtom2,bondtype,
                      bondStereo, bondTopology, 0.0);
              if (pBond) {
                pBondAtom1->addBondedAtom(pBondAtom2);
                pBondAtom2->addBondedAtom(pBondAtom1);
              }
            }
          }
        }
      }
    }
    // - Clean Up and update - //
    //delete pthisAtom;
    imol2.close();

    // With a mol2 file, I assume that all bonds and bond types are defined correctly
    pConnections->assignAngles(pMol);
    pConnections->assignTorsions(pMol);
    pConnections->assignImpropers(pMol);

    pMol->setKind(2); // drug/small molecule

    pMol->bBondsAssigned = true;
    pMol->bAnglesAssigned = true;
    pMol->bTorsionsAssigned = true;
    pMol->bImpropersAssigned = true;
}

// =========================================================
// Function : Read
// ---------------------------------------------------------
// parsers a mol2 file
// ---------------------------------------------------------
// Format:
//
//@<TRIPOS>MOLECULE
// -> mol name
// -> num_atoms [num_bons [num_sust [num_feat [num_sets]]]]
// -> mol_type
// -> charge_type
// -> [status_bits
// -> [mol_comments]]
//
// @<TRIPOS>ATOM
// -> atom_id atom_name x y z atom_type [sust_id [subst_name [charge [ status_bit]]]]
//
//@<TRIPOS>BOND
// -> bond_id origin_atom_id target_atom_id bond_type [status_bits]
//
// =========================================================
void mol2Parser::Read(const std::string &mol2file, collection* pCollection)
{
    connections* pConnections = new connections(pCollection);
    int end   = mol2file.length();
    int slash = mol2file.find_last_of("/");
    std::string file_name = mol2file.substr(slash+1,(end-slash-6));

    std::ifstream imol2;
    imol2.open(mol2file.c_str());

    if (!imol2) {
      std::cout << "\nUNABLE TO OPEN MOL2 FILE"
                << "\nFILENAME = " << mol2file << std::endl;
      return;
    }

    std::string fileline;
    int natoms = 0, nbonds = 0, at1 = 0, at2 = 0;
    int bondtype = 0;
    int bondStereo = 0;
    int bondTopology = 0;

    //AtomLine* pthisAtom = new AtomLine;
    chargeType="";
    hasTotalChargeRemark = false;

    // - Read the file - //
    while (imol2) {
      std::string buffer(80,'*');
      getline(imol2,fileline);

      // - CHARGE TYPE - //
      if (fileline.substr(0,11) == "NO_CHARGES") {
            chargeType="NO_CHARGES";
      }
      if (fileline.substr(0,6) == "DEL_RE") {
            chargeType="DEL_RE";
      }
      if (fileline.substr(0,9) == "GASTEIGER") {
            chargeType="GASTEIGER";
      }
      if (fileline.substr(0,9) == "GAST_HUCK") {
            chargeType="GAST_HUCK";
      }
      if (fileline.substr(0,6) == "HUCKEL") {
            chargeType="HUCKEL";
      }
      if (fileline.substr(0,7) == "PULLMAN") {
            chargeType="PULLMAN";
      }
      if (fileline.substr(0,15) == "GAUSS80_CHARGES") {
            chargeType="GAUSS80_CHARGES";
      }
      if (fileline.substr(0,13) == "AMPAC_CHARGES") {
            chargeType="AMPAC_CHARGES";
      }
      if (fileline.substr(0,16) == "MULLIKEN_CHARGES") {
            chargeType="MULLIKEN_CHARGES";
      }
      if (fileline.substr(0,12) == "DICT_CHARGES") {
            chargeType="DICT_CHARGES";
      }
      if (fileline.substr(0,14) == "MMFF94_CHARGES") {
            chargeType="MMFF94_CHARGES";
      }
      if (fileline.substr(0,12) == "USER_CHARGES") {
            chargeType="USER_CHARGES";
      }

      // - CHARGE - //
      if (fileline.substr(0,6) == "REMARK") {
        std::vector<std::string> splitstr;
        splitString(fileline," ", splitstr, 0);
        if(splitstr[1].compare("CHARGE")==0)
        {
            totalCharge = atof((char*)splitstr[2].c_str());
            hasTotalChargeRemark = true;
        }
      }

      // - MOLECULE - //
      if (fileline.substr(0,17) == "@<TRIPOS>MOLECULE") {
        pMol = pCollection->addMolecule();

        pMol->setMolId(pCollection->getNumberMolecules());
        pMol->inFileType = "mol2";

        //pSmol = pMol->addSubMolecule();

        getline(imol2,fileline);
        pMol->setName(fileline);
        getline(imol2,fileline);

        std::vector<std::string> splitstr;
        splitString(fileline," ", splitstr, 0);

        natoms = atoi((char*)splitstr[0].c_str());
        nbonds = atoi((char*)splitstr[1].c_str());
      }
      // - ATOM - //
      if (fileline.substr(0,13) == "@<TRIPOS>ATOM") {

        std::string residueID="";
        std::vector<AtomLine*> atomLineVector;
        int n = 1;
        for (n = 1; n <= natoms; n ++) {
          getline(imol2,fileline);

           AtomLine* pthisAtom = new AtomLine;
          if (ReadAtomLine(fileline,pthisAtom)) {
            if(n==1)
            {
                residueID=pthisAtom->segID;
            }
            if(residueID.compare(pthisAtom->segID)!=0)
            {
                buildupSubmolecule(pCollection, residueID, n, atomLineVector);
            }
            atomLineVector.push_back(pthisAtom);
            residueID=pthisAtom->segID;
          }
          else {
                buildupSubmolecule(pCollection, residueID, n, atomLineVector);
            return;
          }
/*
          if (ReadAtomLine(fileline,pthisAtom)) {
            pAtom = pSmol->addAtom();
            pAtom->setElement(pCollection->pElements->getElement(pthisAtom->element));

            std::stringstream ss;
            ss << n;
            std::string at_name_n = ss.str().c_str();
            std::string temp_n_name = (pthisAtom->element)+at_name_n;
            pAtom->setName(temp_n_name);

            pAtom->setCoords(pthisAtom->x,pthisAtom->y,pthisAtom->z);
            // pAtom->setNumBonds(0);
            pAtom->setFileID(n);
            pMol->setMaxFileID(n);
          }
          else {
            return;
          }
*/
        }
        buildupSubmolecule(pCollection, residueID, n-1, atomLineVector);
      }

      // - BOND - //
      //      1      1      2     ar
      if (fileline.substr(0,13) == "@<TRIPOS>BOND") {
        for (int n = 1; n <= nbonds; n ++) {
          getline(imol2,fileline);
          if (fileline.length() < 4) {
            std::cout << "NOT A VALID MOL2 FILE\nThe Bond section"
                      << " should contain:\nbond_id origin_atom_id"
                      << " target_atom_id bond_type [status_bits]"
                      << std::endl;
            return;
          }
          else {
            std::vector<std::string> split_bond_str;
            splitString(fileline," ", split_bond_str, 0);

            at1 = atoi((char*)split_bond_str[1].c_str());
            at2 = atoi((char*)split_bond_str[2].c_str());
            if (split_bond_str[3] == "ar") {
              bondtype = 4;
            }
            else {
              bondtype = atoi((char*)split_bond_str[3].c_str());
            }

            pBondAtom1 = pSmol->getAtom(at1);
            pBondAtom2 = pSmol->getAtom(at2);

            if (pBondAtom1 and pBondAtom2) {
              pBond = 0;
              pBond = pMol->addBond(pBondAtom1,pBondAtom2,bondtype,
                      bondStereo, bondTopology, 0.0);
              if (pBond) {
                pBondAtom1->addBondedAtom(pBondAtom2);
                pBondAtom2->addBondedAtom(pBondAtom1);
              }
            }
          }
        }

        // With a mol2 file, I assume that all bonds and bond types are defined correctly
        pConnections->assignAngles(pMol);
        pConnections->assignTorsions(pMol);
        pConnections->assignImpropers(pMol);

        pMol->setKind(2); // drug/small molecule

        pMol->bBondsAssigned = true;
        pMol->bAnglesAssigned = true;
        pMol->bTorsionsAssigned = true;
        pMol->bImpropersAssigned = true;
      }
    }

    // - Clean Up and update - //
    //delete pthisAtom;
    imol2.close();
}

void mol2Parser::buildupSubmolecule(collection* pCollection, std::string& residueID, unsigned int n, std::vector<AtomLine*>& atomLineVector)
{
    pSmol = pMol->addSubMolecule();
    for (unsigned int index=0; index < atomLineVector.size(); index++) {
      if (index == 0) {

        std::string segID = atomLineVector[index]->segID;

        std::string resName = "";
        std::string subMolIdStr = "";

        int subMolId = 0;

        for (unsigned int c = 0; c < segID.size(); c++) {
          if (isalpha(segID[c])) {
            resName += segID[c];
          }
          else if (isdigit(segID[c])) {
            subMolIdStr += segID[c];
          }
        }

        if (resName.size() > 0) {
          pSmol->setName(resName);
          pSmol->set1LName(this->get1LCode(resName));
        }

        if (subMolIdStr.size() > 0) {
          subMolId = atoi(subMolIdStr.c_str());
        }
/*
        boost::regex expression("(\\D+)(\\d+)");

        boost::cmatch what;
        if (boost::regex_match(segID.c_str(), what, expression)) {
          std::string resName(what[1]);
          int subMolId = atoi(what[2].str().c_str());
          pSmol->setName(resName);
          pSmol->set1LName(this->get1LCode(resName));
          pSmol->setSubMolId(subMolId);
        }
        else {
           //@todo should this throw an exception?
        }

*/
        pSmol->setIndex(atomLineVector[index]->resSeq);
        //@todo find out if mol2 has insertion code pSmol->setiCode(iCode);
      }
      pAtom = pSmol->addAtom();
      pAtom->setElement(pCollection->pElements->getElement(atomLineVector[index]->element));

      std::stringstream ss;
      ss << n;
      std::string at_name_n = ss.str().c_str();
      std::string temp_n_name = (atomLineVector[index]->element)+at_name_n;
      pAtom->setName(temp_n_name);

      pAtom->setCoords(atomLineVector[index]->x,atomLineVector[index]->y,atomLineVector[index]->z);
      if (atomLineVector[index]->charge.size()>0)pAtom->setFormalCharge(atoi(atomLineVector[index]->charge.c_str()));
      // pAtom->setNumBonds(0);
      pAtom->setFileID(n);
    }
    pMol->setMaxFileID(n);
    atomLineVector.clear();
}

// ==============================================================
// Function : Write
// --------------------------------------------------------------
// Write a XYZ file.
// ==============================================================
void mol2Parser::Write(const std::string &mol2file, molecule* pMolecule)
{
    std::ofstream omol2;
    omol2.open(mol2file.c_str());

    if (!omol2 or (pMolecule == 0)) {
      std::cout << "\nUNABLE TO OPEN MOL2 FILE"
                << "\nFILENAME = " << mol2file << std::endl;
      return;
    }

    // ATOMS
    typedef std::vector<atom*>::iterator AtomIterator;
    std::vector<atom*>::iterator atomIterator1;
    std::vector<atom*>::iterator atomIterator2;
    std::vector<atom*>                   Atoms;
    Atoms  = pMolecule->getAtomList();

    // BONDS
    typedef std::map<int, Bond*>::iterator BondMapIterator;
    std::map<int, Bond*>      moleculeBondMap;
    Bond* pBond;
    moleculeBondMap = pMolecule->getBondMap();
    int nBonds = 0;

    // GET NUMBER OF BONDS -- no intermolecular bonds
    if (!moleculeBondMap.empty()) {
      for (BondMapIterator b = moleculeBondMap.begin(); b!=moleculeBondMap.end();b++) {
        pBond = b->second;
        atomIterator1 = std::find(Atoms.begin(), Atoms.end(), pBond->atom1);
        atomIterator2 = std::find(Atoms.begin(), Atoms.end(), pBond->atom2);
        if ((atomIterator1 != Atoms.end()) and (atomIterator2 != Atoms.end())) {
          nBonds++;
        }
      }
    }

    omol2 << "@<TRIPOS>MOLECULE\n" << std::endl;
    omol2 << Atoms.size() << "  " << nBonds << "\n\n\n" << std::endl;
    omol2 << "@<TRIPOS>ATOM" << std::endl;

    // - Loop over atoms - //
    int atomIndex = 0;
    double curX = 0.0;
    double curY = 0.0;
    double curZ = 0.0;
    int n = 1;

    for (AtomIterator d=Atoms.begin(); d != Atoms.end(); d++) {
      pAtom = *d;
      curX = pAtom->getX();
      curY = pAtom->getY();
      curZ = pAtom->getZ();

      omol2 << std::setw(7) << n << std::setw(2) << pAtom->getElement()->symbol
            << " " << std::setiosflags(std::ios::fixed) << std::setprecision(6);
      omol2 << std::setw(12) << curX<< " ";
      omol2 << std::setw(12) << curY << " ";
      omol2 << std::setw(12) << curZ << " ";
      //omol2 << pAtom->getAtomType()->getName(); /// this needs to be fixed
      omol2 << std::endl;

      // - Create a map to ensure that the numbering of the mol file is correct - //
      itsAtomMap[pAtom->getIndex()] = n;
      n=n+1;
      atomIndex++;
    }
    omol2 << "@<TRIPOS>BOND" << std::endl;

    n = 1;
    // - Loop over bonds - //
    if (!moleculeBondMap.empty()) {
      for (BondMapIterator b = moleculeBondMap.begin(); b!=moleculeBondMap.end();b++) {
        pBond = b->second;
        atomIterator1 = std::find(Atoms.begin(), Atoms.end(), pBond->atom1);
        atomIterator2 = std::find(Atoms.begin(), Atoms.end(), pBond->atom2);

        if ((atomIterator1 != Atoms.end()) and (atomIterator2 != Atoms.end())) {
          int bondType = pBond->type;
          if (bondType == 6) bondType = 1;
          if (bondType == 7) bondType = 2;
          omol2 << std::setw(6) << n << std::setw(5) << pBond->atom1->getIndex() << std::setw(5)
                << pBond->atom2->getIndex() << std::setw(2) << bondType;
          omol2 << std::endl;
          n++;
        }
      }
    }
    omol2.close();
}

// ==========================================
// Function : ReadAtomLine
// ------------------------------------------
// Reads ATOM AND HETATOM line and returns data
// ==========================================
bool mol2Parser::ReadAtomLine(std::string &fileline, AtomLine* pthisAtom)
{
    std::vector<std::string> split_atom_str;
    splitString(fileline," ", split_atom_str, 0);

    if ( split_atom_str.size() < 6) {
      std::cout << "The ATOM line should contain:\n atom_id atom_name x y z"
           << " atom_type [sust_id [subst_name [charge [ status_bit]]]]\n"
           << std::endl;
      std::cout << split_atom_str.size() << std::endl;
      std::cout << fileline << std::endl;
      return false;
    }
    else {
      pthisAtom->x       = atof((char*)split_atom_str[2].c_str());
      pthisAtom->y       = atof((char*)split_atom_str[3].c_str());
      pthisAtom->z       = atof((char*)split_atom_str[4].c_str());
      pthisAtom->typ     = split_atom_str[5].c_str();
      //pthisAtom->element = (split_atom_str[5].substr(0,1).c_str());
      pthisAtom->element = char(toupper(split_atom_str[5].substr(0,1).c_str()[0]));
      std::string name   = split_atom_str[5].c_str();

      pthisAtom->charge="";
      if(split_atom_str.size()>6)
      {
        pthisAtom->resSeq=atoi(split_atom_str[6].c_str());
      }
      if(split_atom_str.size()>7)
      {
        pthisAtom->segID=split_atom_str[7];
      }
      if(split_atom_str.size()>8)
      {
        pthisAtom->charge=split_atom_str[8];
      }

      if (!isalpha(name[1])) {
        split_atom_str.clear();
        // - All is GOOD - //
        return true;
      }
      else {
        if (pthisAtom->element == "C") {
          if (GetAlphaChar(name,2) != "0") {
            std::string second_letter = GetAlphaChar(name,2);

            if (second_letter == "L" || second_letter == "l" ) {
              pthisAtom->element  = "Cl";
            }

            if ((second_letter == "A" || second_letter == "a") &&
                (removeCharacter((pthisAtom->name),' ') == "CA")) {
              pthisAtom->element  = "Ca";
            }
            if (second_letter == "U" || second_letter == "u" ) {
              pthisAtom->element  = "Cu";
            }
            if (second_letter == "R" || second_letter == "r" ) {
              pthisAtom->element  = "Cr";
            }
            if (second_letter == "O" || second_letter == "o" ) {
              pthisAtom->element  = "Co";
            }
          }
          else {
            pthisAtom->element  = "C";
          }
        }

        // - HYDROGEN AND HELIUM - //
        if (pthisAtom->element == "H") {
          if (GetAlphaChar(name,2) != "0") {
            std::string second_letter = GetAlphaChar(name,2);

            if (second_letter == "E" || second_letter == "e" ) {
              pthisAtom->element  = "He";
            }
          }
          else {
            pthisAtom->element  = "H";
          }
        }

        // - LITHIUM - //
        if (pthisAtom->element == "L") {
          if (GetAlphaChar(name,2) != "0") {
            std::string second_letter = GetAlphaChar(name,2);

            if (second_letter == "I" || second_letter == "i" ) {
              pthisAtom->element  = "Li";
            }
          }
          else {
            pthisAtom->element  = "L";
            std::cout << "Unknown element symbol : " << pthisAtom->element
                      <<  std::endl;
          }
        }

        // - NITROGEN, SODIUM and NEON - //
        if (pthisAtom->element == "N") {
          if (GetAlphaChar(name,2) != "0") {
            std::string second_letter = GetAlphaChar(name,2);

            if (second_letter == "A" || second_letter == "a" ) {
              pthisAtom->element  = "Na";
            }
            if (second_letter == "E" || second_letter == "e" ) {
              pthisAtom->element  = "Ne";
            }
          }
        }

        // - MAGNESIUM and MANGANESE - //
        if (pthisAtom->element == "M") {
          if (GetAlphaChar(name,2) != "0") {
            std::string second_letter = GetAlphaChar(name,2);

            if (second_letter == "G" || second_letter == "g" ) {
              pthisAtom->element  = "Mg";
            }
            if (second_letter == "N" || second_letter == "n" ) {
              pthisAtom->element  = "Mn";
            }
          }
          else {
            pthisAtom->element  = "M";
            std::cout << "Unknown element symbol : "
                      << pthisAtom->element <<  std::endl;
          }
        }

        // - BROMINE and BORON - //
        if (pthisAtom->element == "B" ) {
          if (GetAlphaChar(name,2) != "0") {
            std::string second_letter = GetAlphaChar(name,2);

            if (second_letter == "R" || second_letter == "r" ) {
              pthisAtom->element  = "Br";
            }
            if (second_letter == "E" || second_letter == "e" ) {
              pthisAtom->element  = "Be";
            }
          }
          else {
            pthisAtom->element  = "B";
          }
        }

        // - FLOURINE and IRON - //
        if (pthisAtom->element == "F" ) {
          if (GetAlphaChar(name,2) != "0") {
            std::string second_letter = GetAlphaChar(name,2);

            if (second_letter == "E" || second_letter == "e" ) {
              pthisAtom->element  = "Fe";
            }
          }
          else {
            pthisAtom->element  = "F";
          }
        }

        // - ZINC  - //
        if (pthisAtom->element == "Z" ) {
          if (GetAlphaChar(name,2) != "0" ) {
            std::string second_letter = GetAlphaChar(name,2);

            if (second_letter == "N" || second_letter == "n" ) {
              pthisAtom->element  = "Zn";
            }
            if (second_letter == "R" || second_letter == "r" ) {
              pthisAtom->element  = "Zr";
            }
          }
        }
      }
    }

    split_atom_str.clear();

     // - All is GOOD - //
     return true;
}

// =========================================================
// Function : get1LCode
// ---------------------------------------------------------
// Get 1-Letter code
// =========================================================
std::string mol2Parser::get1LCode(std::string s)
{
    std::string oneLetterCode;

    nameMapIterator p = this->res21l.find(s);

    if (p != this->res21l.end()) {
      return p->second;
    }
    return " ";
}

// ==========================================
// Function : getAtom
// ------------------------------------------
// returns an atom of given index
// ==========================================
int mol2Parser::getAtom(int index)
{
   std::map<int, int>::iterator iter;
   iter = itsAtomMap.find(index);

   if (iter != itsAtomMap.end()) {
     return itsAtomMap[index];
   }
   return 0;
}

} // MTKpp namespace


