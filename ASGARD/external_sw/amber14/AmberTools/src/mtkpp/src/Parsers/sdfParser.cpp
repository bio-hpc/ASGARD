/*!
   \file sdfParser.cpp
   \brief Parses sd files
   \author Martin Peters

   Reads and writes sd files

   $Date: 2010/03/29 20:39:35 $
   $Revision: 1.20 $

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

#include "sdfParser.h"

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
#include "Utils/constants.h"

#include "StringManip.h"

namespace MTKpp
{

// ============================================================
// Function : sdfParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
sdfParser::sdfParser():baseParser() {
    pMolecule = 0;
    pSubMolecule  = 0;
    pBond = 0;
    pBondAtom1 = 0;
    pBondAtom2 = 0;
    pConnections = 0;
}

// ============================================================
// Function : sdfParser()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
sdfParser::~sdfParser() {
    delete pConnections;
}

// =========================================================
// Function : Read
// ---------------------------------------------------------
// parsers a sdf file
// ---------------------------------------------------------
// Format: mol block
// 1st line       : Compound Name
// 2nd line       : initials, program, date, dimension, scale factors, energy, reg name
// 3rd line       : Comments
// 4rd line       : Number of atoms(N), Number of bonds(B), atlist, obs, chiral, nentries, obs, obs, obs, obs, nprop
//                  -- N cannot exceed 999 atoms --
// 5th line       : x y z atom  massdiff charge stereo hcount stcare valence 0 0 0 0 0 0
//    .
//    .
// (N+4)th line   : x y z atom charge charge stereo hcount stcare valence 0 0 0 0 0 0
// (N+5)th line   : at1 at2 bondtype bondstereo 0  0  0
//    .
//    .
// (N+4+B)th line : at1 at2 bondtype bondstereo 0  0  0
// M  END
// properties
// $$$$
// =========================================================
void sdfParser::Read(const std::string &sdffile, collection* pCollection, const bool &bohr)
{
    pConnections = new connections(pCollection);
    std::ifstream isdf;
    isdf.open(sdffile.c_str());

    if (!isdf) {
      std::cout << "\nUNABLE TO OPEN SDF FILE"
                << "\nFILENAME = " << sdffile << std::endl;
      return;
    }

    std::string fileline = "";
    std::string name = "";
    std::string molecule_name = "";
    std::string initials = "XX", program = "", date = " ", dim = " ";
    std::string comment = "";
    int natoms = 0, nbonds = 0;
    std::string v2 = "V2000";
    std::string atom = "";

    AtomLine* pthisAtom = new AtomLine;

    while (!isdf.eof()) {
      natoms = 0;
      nbonds = 0;

      for (int m = 0; m < 4; m++) {
        if (!isdf.eof()) {
          getline(isdf,fileline);
          if (m == 0) {
            molecule_name = fileline;
          }
          else if (m == 3) {
            if (! (fileline.length() >= 6)) {
              return;
            }
            else {
              natoms = atoi((fileline.substr(0,3)).c_str());
              nbonds = atoi((fileline.substr(3,3)).c_str());
            }
          }
        }
        else {
          break;
        }
      }

      if (isdf.eof()) {
        break;
      }

      pMolecule  = pCollection->addMolecule();
      pMolecule->setName(molecule_name);
      pMolecule->setMolId(pCollection->getNumberMolecules());
      pMolecule->inFileType = "sdf";

      pSubMolecule = pMolecule->addSubMolecule();
      pSubMolecule->setName(molecule_name);
      pSubMolecule->setSubMolId(pMolecule->getNumSubMolecules());

      // - FIFTH LINE: - //
      // - LOOP OVER NATOMS AND GET ATOM, X, Y, Z. - //
      for (int n = 1; n <= natoms; n ++) {
        getline(isdf,fileline);
        if (isdf.eof()) {
          return;
        }
        if (ReadAtomLine(fileline,pthisAtom)) {
          pAtom = pSubMolecule->addAtom();

          pAtom->setElement(pCollection->pElements->getElement(pthisAtom->element));
          pAtom->setFileID(n);

          if (bohr) {
            pAtom->setCoords(pthisAtom->x/ANG2BOHR, pthisAtom->y/ANG2BOHR,
                             pthisAtom->z/ANG2BOHR);
          }
          else {
            pAtom->setCoords(pthisAtom->x,pthisAtom->y,pthisAtom->z);
          }
        }
        else {
          return;
        }
      }

      // - (N+1)th LINE: - //
      // - LOOP OVER NBONDS AND GET AT1 and AT2. - //
      int at1 = 0, at2 = 0, bondType = 0, bondStereo = 0, bondTopology = 0;
      for (int n = 1; n <= nbonds; n++) {
        if (isdf.eof()) {
          return;
        }
        getline(isdf,fileline);
        if (fileline.length() < 9) {
          std::cout << "NOT A VALID SDF FILE" << std::endl;
          std::cout << "line (NATOMS+5) of a sdf file should have at least 9 characters" << std::endl;
          return;
        }
        else {
          at1 = atoi((fileline.substr(0,3)).c_str());
          at2 = atoi((fileline.substr(3,3)).c_str());

          try {
            bondType = atoi((fileline.substr(6,3)).c_str());
          }
          catch (std::out_of_range) {
            bondType = 0;
          }

          try {
            bondStereo = atoi((fileline.substr(9,3)).c_str());
          }
          catch (std::out_of_range) {
            bondStereo = 0;
          }

          try {
            bondTopology = atoi((fileline.substr(15,3)).c_str());
          }
          catch (std::out_of_range) {
            bondTopology = 0;
          }

          pBondAtom1 = pMolecule->getAtom(at1, 1, 0);

          pBondAtom2 = pMolecule->getAtom(at2, 1, 0);

          if (pBondAtom1 and pBondAtom2) {
            pBond = pMolecule->addBond(pBondAtom1, pBondAtom2, bondType, bondStereo, bondTopology, 0.0);

            pBondAtom1->addBondedAtom(pBondAtom2);
            pBondAtom2->addBondedAtom(pBondAtom1);
          }
        }
      }

      int myMolCharge = 0;

      for (int n = 1; n <= 999; n ++) {
        getline(isdf,fileline);
        if ((fileline.substr(0,6)) == "M  END" ) {
          break;
        }

        if ((fileline.substr(0,6)) == "M  CHG") {
            unsigned int nCHG = atoi((fileline.substr(6,3)).c_str());
            //int st = 10;
            int a,c;

            std::vector<std::string> nWords;
            splitString(fileline, " ", nWords, 0);

            if (nWords.size() == (3 + nCHG*2)) {
              unsigned int x = 3;
              unsigned int y = 4;
              for (unsigned int j = 0; j < nCHG; j++) {
                a = string2UInt(nWords[x]);
                c = string2UInt(nWords[y]);

                pAtom = pMolecule->getAtom(a, 1, 0);

                if (pAtom) {
                  pAtom->setFormalCharge(c);
                }
                else {
                  std::cout << " Error in sdf " << std::endl;
                }
                myMolCharge+=c;

                x += 2;
                y += 2;
              }
            }
/*
            for (int j = 0; j < nCHG; j++) {

              a = atoi((fileline.substr(st,4)).c_str());
              c = atoi((fileline.substr(st+4,4)).c_str());
              pAtom = pMolecule->getAtom(a, 1, 0);

              if (pAtom) {
                pAtom->setFormalCharge(c);
              }
              else {
                std::cout << " Error in sdf " << std::endl;
              }
              st+=4;
              myMolCharge+=c;
            }
*/
            //std::cout << " mol charge = " << myMolCharge << std::endl;
            pMolecule->setTotalCharge(myMolCharge);

        }
      }

      while (isdf) {
        getline(isdf,fileline);
        if ((fileline.substr(0,4)) == "$$$$") {
          break;
        }
        if ((fileline.substr(0,1)) == ">") {
          int st = (fileline.substr(0,fileline.size())).find_first_of("<");
          int ed = (fileline.substr(0,fileline.size())).find_last_of(">");
          int end = ed - (st+1);
          std::string propertyName = fileline.substr(st+1,end);
          getline(isdf,fileline);
          pMolecule->addProperty(propertyName, fileline);
          getline(isdf,fileline);
        }
      }
    }
    // - Clean Up and update - //
    delete pthisAtom;
    pConnections->assignAngles();
    pConnections->assignTorsions();
    pConnections->assignImpropers();
    isdf.close();
}

// =========================================================
// Function : Read
// ---------------------------------------------------------
// parsers a molecule from a sd file
// ---------------------------------------------------------
//
// =========================================================
int sdfParser::ReadMolecule(std::ifstream& isdf, molecule* pMolecule, const bool &bohr)
{
    if (!isdf) {
      return 12;
    }
    if (isdf.eof()) {
      return 1;
    }

    std::string fileline = "";
    std::string name = "";
    std::string molecule_name = "";
    std::string initials = "XX", program = "", date = " ", dim = " ";
    std::string comment = "";
    int natoms = 0, nbonds = 0;
    std::string v2 = "V2000";
    std::string atom = "";

    AtomLine* pthisAtom = new AtomLine;

    pCollection = pMolecule->getParent();
    if (!pConnections) {
      pConnections = new connections(pCollection);
    }

    for (int m = 0; m < 4; m++) {
      if (!isdf.eof()) {
        getline(isdf,fileline);
        if (m == 0) {
          molecule_name = fileline;
        }
        else if (m == 3) {
          if (! (fileline.length() >= 6)) {
            return 12;
          }
          else {
            natoms = atoi((fileline.substr(0,3)).c_str());
            nbonds = atoi((fileline.substr(3,3)).c_str());
          }
        }
      }
      else {
        if (m == 1 and isdf.eof()) {
          return 1;
        }
        return 2;
      }
    }

    if (isdf.eof()) {
      return 3;
    }

    pMolecule->setName(molecule_name);
    pMolecule->setMolId(pCollection->getNumberMolecules());
    pMolecule->inFileType = "sdf";

    pSubMolecule = pMolecule->addSubMolecule();
    pSubMolecule->setName(molecule_name);
    pSubMolecule->setSubMolId(pMolecule->getNumSubMolecules());

    // - FIFTH LINE: - //
    // - LOOP OVER NATOMS AND GET ATOM, X, Y, Z. - //
    for (int n = 1; n <= natoms; n ++) {
      getline(isdf,fileline);
      if (isdf.eof()) {
        return 4;
      }
      if (ReadAtomLine(fileline, pthisAtom)) {
        pAtom = pSubMolecule->addAtom();

        pAtom->setElement(pCollection->pElements->getElement(pthisAtom->element));
        pAtom->setFileID(n);

        if (bohr) {
          pAtom->setCoords(pthisAtom->x/ANG2BOHR, pthisAtom->y/ANG2BOHR,
                           pthisAtom->z/ANG2BOHR);
        }
        else {
          pAtom->setCoords(pthisAtom->x,pthisAtom->y,pthisAtom->z);
        }
      }
      else {
        return 5;
      }
    }

    // - (N+1)th LINE: - //
    // - LOOP OVER NBONDS AND GET AT1 and AT2. - //
    int at1 = 0, at2 = 0, bondType = 0, bondStereo = 0, bondTopology = 0;
    for (int n = 1; n <= nbonds; n++) {
      if (isdf.eof()) {
        return 6;
      }
      getline(isdf,fileline);
      if (fileline.length() < 9) {
        std::cout << "NOT A VALID SDF FILE" << std::endl;
        std::cout << "line (NATOMS+5) of a sdf file should have at least 9 characters" << std::endl;
        return 7;
      }
      else {
        unsigned int lineLength = fileline.size();
        at1 = atoi((fileline.substr(0,3)).c_str());
        at2 = atoi((fileline.substr(3,3)).c_str());
        if (lineLength >= 9) bondType = atoi((fileline.substr(6,3)).c_str());
        if (lineLength >= 12) bondStereo = atoi((fileline.substr(9,3)).c_str());
        if (lineLength >= 15) bondTopology = atoi((fileline.substr(15,3)).c_str());

        pBondAtom1 = pMolecule->getAtom(at1, 1, 0);

        pBondAtom2 = pMolecule->getAtom(at2, 1, 0);

        pBond = pMolecule->addBond(pBondAtom1, pBondAtom2, bondType, bondStereo, bondTopology, 0.0);

        pBondAtom1->addBondedAtom(pBondAtom2);
        pBondAtom2->addBondedAtom(pBondAtom1);
        bondType = 0, bondStereo = 0, bondTopology = 0;
      }
    }

    for (int n = 1; n <= 999; n ++) {
      getline(isdf,fileline);
      if ((fileline.substr(0,6)) == "M  END" ) {
        break;
      }
    }
    while (isdf) {
      getline(isdf,fileline);
      if ((fileline.substr(0,4)) == "$$$$") {
        delete pthisAtom;
        return 0;
      }
      if ((fileline.substr(0,1)) == ">") {
        int st = (fileline.substr(0,fileline.size())).find_first_of("<");
        int ed = (fileline.substr(0,fileline.size())).find_last_of(">");
        int end = ed - (st+1);
        std::string propertyName = fileline.substr(st+1,end);
        getline(isdf,fileline);
        pMolecule->addProperty(propertyName, fileline);
        getline(isdf,fileline);
      }
    }

    delete pthisAtom;
    pConnections->assignAngles(pMolecule);
    pConnections->assignTorsions(pMolecule);
    pConnections->assignImpropers(pMolecule);
    return 8;
}

// ============================================================
// Function : numMolecules
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
//
// ============================================================
int sdfParser::numMolecules(const std::string &sdffile)
{
    std::ifstream isdf;
    isdf.open(sdffile.c_str());

    if (!isdf) {
      std::cout << "\nUNABLE TO OPEN SDF FILE"
                << "\nFILENAME = " << sdffile << std::endl;
      return -1;
    }

    std::string fileline = "";

    int nMolecules = 0;
    while (isdf) {
      getline(isdf,fileline);
      if ((fileline.substr(0,4)) == "$$$$") {
        nMolecules++;
      }
    }

    isdf.close();

    return nMolecules;
}

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write a SDF file.
// ------------------------------------------------------------
// Format: mol block
// 1st line       : Compound Name
// 2nd line       :
// 3rd line       :
// 4rd line       : Number of atoms(N), Number of bonds(B)
// 5th line       : x y z atom charge 0  0  0  0  0  0  0  0  0  0  0
//    .
//    .
// (N+4)th line   : x y z atom charge 0  0  0  0  0  0  0  0  0  0  0
// (N+5)th line   : at1 at2 1  0  0  0  0
//    .
//    .
// (N+4+B)th line : at1 at2 1  0  0  0  0
// M  END
// ============================================================
void sdfParser::Write(const std::string &sdffile, collection* pCollection)
{
    std::ofstream osdf;
    osdf.open(sdffile.c_str());

    if (!osdf) {
      std::cout << "\nUNABLE TO OPEN SDF FILe"
                << "\nFILENAME = " << sdffile << std::endl;
      return;
    }
    unsigned int n_atoms = 0;
    unsigned int n_bonds = 0;

    // ==============================================================
    // Write the file,
    // Loop over molecules
    // ==============================================================

    int n = 1;
    std::vector<molecule*> molList = pCollection->getMoleculeList();
    for (moleculeIterator m = molList.begin(); m != molList.end(); m++) {
      pMolecule = *m;
      std::string mol_name = pMolecule->getName();
      n_atoms  = pMolecule->getNumAtoms();
      n_bonds  = pMolecule->getNumBonds();

      osdf << std::left << std::setw(80) << mol_name << std::endl;
      osdf << " Molecule \n" << std::endl;
      osdf << std::setw(3) << n_atoms << std::setw(3) << n_bonds << "  0  0  1  0  0  0  0  0999 V2000" << std::endl;

      std::vector<atom*> atomList = pMolecule->getAtomList();
      // - Loop over the Atoms - //
      for (atomIterator a = atomList.begin(); a != atomList.end(); a++) {
        pAtom = *a;
        char temp[100];
        sprintf(temp,"%10.4f%10.4f%10.4f %-3s%2d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d",
                pAtom->getX(), pAtom->getY(), pAtom->getZ(),
                (pAtom->getElement()->symbol.c_str()),
                0,0,0,0,0,0,0,0,0,0,0,0);
        osdf << temp << std::endl;

        // - Create a map to ensure that the numbering of the mol file is correct - //
        //itsAtomMap[pAtom->getFileID()] = n;
        itsAtomMap[pAtom->getColIndex()] = n;
        n=n+1;
      }

      // - Loop over the Bonds - //
      moleculeBondMap = pMolecule->getBondMap();
      if (!moleculeBondMap.empty()) {
        for (BondMapIterator b = moleculeBondMap.begin(); b != moleculeBondMap.end(); b++) {
          pBond = b->second;
          char temp[100];
          int bondType = pBond->type;
          if (bondType == 6) bondType = 1;
          if (bondType == 7) bondType = 2;
          sprintf(temp,"%3d%3d%3d%3d%3d%3d%3d",
                  //pBond->atom1->getFileID(),pBond->atom2->getFileID(),
                  itsAtomMap[pBond->atom1->getColIndex()], itsAtomMap[pBond->atom2->getColIndex()],
                  bondType,
                  pBond->stereo,
                  pBond->topology,
                  0,0);
          osdf << temp << std::endl;
        }
      }

      // - Write Formal Charge information - //
      int fc = 0;
      for (atomIterator d=atomList.begin(); d != atomList.end(); d++) {
        pAtom = *d;
        fc = pAtom->getFormalCharge();
        if (fc != 0) {
          char temp[100];
          sprintf(temp,"%-6s%3d%4d%4d",
                  //"M  CHG",1,pAtom->getFileID(),fc);
                  "M  CHG", 1, itsAtomMap[pAtom->getColIndex()], fc);
          osdf << temp << std::endl;
        }
      }

      // - Write Molecule Property Information - //
      osdf << "M  END" << std::endl;
      std::map<std::string, std::string> molProperties = pMolecule->getProperties();
      typedef std::map<std::string, std::string>::iterator   PropertyMapIterator;
      if (!molProperties.empty()) {
        for (PropertyMapIterator p = molProperties.begin(); p != molProperties.end(); p++) {
          std::string propName =  p->first;
          std::string propValue = p->second;
          osdf << "> <";
          osdf << propName;
          osdf << ">" << std::endl;
          osdf << propValue << std::endl;
          osdf << " " << std::endl;
        }
      }
      n=1;
      osdf << "$$$$" << std::endl;
    }
    osdf.close();
}

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write a SDF file.
// ------------------------------------------------------------
// ============================================================
void sdfParser::Write(std::ostream& osdf, molecule* pMolecule)
{

    unsigned int n_atoms = 0;
    unsigned int n_bonds = 0;

    int n = 1;

    std::string mol_name = pMolecule->getName();
    n_atoms  = pMolecule->getNumAtoms();
    n_bonds  = pMolecule->getNumBonds();

    osdf << std::left << std::setw(80) << mol_name << "\n Molecule \n" << std::endl;
    osdf << std::setw(3) << n_atoms << std::setw(3) << n_bonds << "  0  0  1  0  0  0  0  0999 V2000" << std::endl;

    std::vector<atom*> atomList = pMolecule->getAtomList();
    // - Loop over the Atoms - //
    for (atomIterator a = atomList.begin(); a != atomList.end(); a++) {
      pAtom = *a;
      char temp[100];
      sprintf(temp,"%10.4f%10.4f%10.4f %-3s%2d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d",
              pAtom->getX(), pAtom->getY(), pAtom->getZ(),
              (pAtom->getElement()->symbol.c_str()),
              0,0,0,0,0,0,0,0,0,0,0,0);
      osdf << temp << std::endl;

      // - Create a map to ensure that the numbering of the mol file is correct - //
      itsAtomMap[pAtom->getColIndex()] = n;
      n = n + 1;
    }

    // - Loop over the Bonds - //
    moleculeBondMap = pMolecule->getBondMap();
    if (!moleculeBondMap.empty()) {
      for (BondMapIterator b = moleculeBondMap.begin(); b != moleculeBondMap.end(); b++) {
        pBond = b->second;
        char temp[100];
        int bondType = pBond->type;
        if (bondType == 6) bondType = 1;
        if (bondType == 7) bondType = 2;
        sprintf(temp,"%3d%3d%3d%3d%3d%3d%3d",
                itsAtomMap[pBond->atom1->getColIndex()], itsAtomMap[pBond->atom2->getColIndex()],
                bondType,
                pBond->stereo,
                pBond->topology,
                0,0);
        osdf << temp << std::endl;
      }
    }

    // - Write Formal Charge information - //
    int fc = 0;
    for (atomIterator d=atomList.begin(); d != atomList.end(); d++) {
      pAtom = *d;
      fc = pAtom->getFormalCharge();
      if (fc != 0) {
        char temp[100];
        sprintf(temp,"%6s%3d%4d%4d",
                "M  CHG", 1, itsAtomMap[pAtom->getColIndex()], fc);
        osdf << temp << std::endl;
      }
    }

    // - Write Molecule Property Information - //
    osdf << "M  END" << std::endl;
    std::map<std::string, std::string> molProperties = pMolecule->getProperties();
    typedef std::map<std::string, std::string>::iterator PropertyMapIterator;
    if (!molProperties.empty()) {
      for (PropertyMapIterator p = molProperties.begin(); p != molProperties.end(); p++) {
        std::string propName =  p->first;
        std::string propValue = p->second;
        osdf << "> <";
        osdf << propName;
        osdf << ">" << std::endl;
        osdf << propValue << std::endl;
        osdf << " " << std::endl;
      }
    }
    n=1;
    osdf << "$$$$" << std::endl;
}

// ============================================================
// Function : ReadAtomLine
// ------------------------------------------------------------
// Reads atom line and returns data
// ============================================================
bool sdfParser::ReadAtomLine(std::string &fileline, AtomLine* pthisAtom)
{
    std::string buffer(80,'*');
    int length = fileline.length();
    if (length < 32) {
      std::cout << "THERE IS A PROBLEM WITH THE COORDINATES\n"
      << "FORMAT:\n" 
      << "1st line       : Compound Name\n"
      << "2nd line       : initials, program, date, dimension, scale factors, energy, reg name\n"
      << "3rd line       : Comments\n"
      << "4rd line       : Number of atoms(N),Number of bonds(B),atlist,obs,chiral,nentries,obs,obs,obs,obs,nprop\n"
      << "      -- N cannot exceed 999 atoms -- \n"
      << "5th line       : x y z atom  massdiff charge stereo hcount stcare valence 0 0 0 0 0 0\n"
      << "    .\n"
      << "    .\n"
      << " (N+4)th line   : x y z atom charge charge stereo hcount stcare valence 0 0 0 0 0 0\n"
      << " (N+5)th line   : at1 at2 bondtype bondstereo 0  0  0\n"
      << "    .\n"
      << "    .\n"
      << " (N+4+B)th line : at1 at2 bondtype bondstereo 0  0  0\n"
      << " M  END\n" << std::endl;
      return false;
    }

    else if (length < 33) {
      fileline=fileline+buffer;
    }

    if (FieldExists(fileline,0,10)) {
      pthisAtom->x      = strtod(fileline.substr(0,10).c_str(), 0);
    }
    else {
      return false;
    }

    if (FieldExists(fileline,10,10)) {
      pthisAtom->y   = strtod(fileline.substr(10,10).c_str(), 0);
    }
    else {
      return false;
    }

    if (FieldExists(fileline,20,10)) {
      pthisAtom->z   = strtod(fileline.substr(20,10).c_str(), 0);
    }
    else {
      return false;
    }

    if (FieldExists(fileline,30,3)) {
      pthisAtom->element  = removeCharacter((fileline.substr(30,3).c_str()),' ');
    }
    else {
      return false;
    }

    // - Not sure what to do with these yet - //
    // - Might not always be in a pdb file - //
/*    pthisAtom->massdiff      = (atoi(fileline.substr(33,2).c_str()));
    pthisAtom->charge        = (atoi(fileline.substr(35,3).c_str()));
    pthisAtom->stereo        = (atoi(fileline.substr(38,3).c_str()));
    pthisAtom->hcount        = (atoi(fileline.substr(41,3).c_str()));
    pthisAtom->strcare       = (atoi(fileline.substr(44,3).c_str()));
    pthisAtom->valence       = (atoi(fileline.substr(47,3).c_str()));
*/
    // - All is GOOD - //
    return true;
}

// ==========================================
// Function : getAtom
// ------------------------------------------
// returns an atom of given index
// ==========================================
int sdfParser::getAtom(int index)
{
    std::map<int, int>::iterator iter;
    iter = itsAtomMap.find(index);

    if (iter != itsAtomMap.end()) {
      return itsAtomMap[index];
    }
    return 0;
}

} // MTKpp namespace

