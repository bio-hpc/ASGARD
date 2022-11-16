/*!
   \file acParser.cpp
   \brief Parses AMBER antechamber files
   \author Martin Peters

   $Date: 2010/03/29 20:39:34 $
   $Revision: 1.10 $

   ----------------------------------------------------------------------------

   MTK++ - C++ package of modeling libraries.

   Copyright (C) 2005-2007  (see AUTHORS file for a list of contributors)

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

#include "acParser.h"

#include "StringManip.h"

#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/submolecule.h"
#include "Molecule/atom.h"
#include "Molecule/bond.h"
#include "Molecule/element.h"
#include "Molecule/stdFrag.h"

#include "Utils/vector3d.h"

namespace MTKpp
{

// ============================================================
// Function : acParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
acParser::acParser():baseParser() {}

// =========================================================
// Function : acParser()
// ---------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// =========================================================
acParser::~acParser() {}

// =========================================================
// Function : Read
// ---------------------------------------------------------
// parsers an ac file
// ---------------------------------------------------------
// Format :
// =========================================================
void acParser::Read(const std::string &acfile, collection* pCollection)
{
    molecule* pMolecule  = pCollection->addMolecule();
    this->Read(acfile, pMolecule);
}

// =========================================================
// Function : Read
// ---------------------------------------------------------
// parsers an ac file
// ---------------------------------------------------------
// Format :
// =========================================================
void acParser::Read(const std::string &acfile, molecule* pMolecule)
{
    std::ifstream iac;
    iac.open(acfile.c_str());

    if (!iac) {
      std::cout << "\nUNABLE TO OPEN AC FILE"
                << "\nFILENAME = " << acfile
                << "\nEXITING...\n" << std::endl;
      return;
    }

    collection* pCollection = pMolecule->getParent();
    pMolecule->setName("Mol");
    pMolecule->setMolId(pCollection->getNumberMolecules());
    pMolecule->inFileType = "pdb";

    submolecule* pSubMolecule = pMolecule->addSubMolecule();
    pSubMolecule->setName("sMol");
    pSubMolecule->setSubMolId(pMolecule->getNumSubMolecules());

    double x,y,z;
    int n = 0;
    atom* pAtom = 0;
    std::string fileline = "";
    bool bATOM = false;
    bool bBOND = false;

    while (iac) {
      std::string buffer(80,'*');
      getline(iac,fileline);
      bATOM = false;
      bBOND = false;
      if (fileline.substr(0,4) == "ATOM") bATOM = true;
      if (fileline.substr(0,4) == "BOND") bBOND = true;

      if (bATOM) {
        n++;
        std::string atomName = fileline.substr(12,4);
        //std::string resName =  fileline.substr(17,3);
        x = strtod(fileline.substr(30,8).c_str(), 0);
        y = strtod(fileline.substr(38,8).c_str(), 0);
        z = strtod(fileline.substr(46,8).c_str(), 0);

        pAtom = pSubMolecule->addAtom();
        pAtom->setElement(pCollection->pElements->getElement(this->determineElement(atomName)));
        pAtom->setName(atomName);
        pAtom->setCoords(x,y,z);
        pAtom->setFileID(n);
        pMolecule->setMaxFileID(n);
      }

      if (bBOND) {
        std::vector<std::string> splitLine;
        splitString(fileline, " ", splitLine, 4);

        int i2  =  atoi(splitLine[1].c_str());
        int i3  =  atoi(splitLine[2].c_str());

        atom* pBondAtom1 = pMolecule->getAtom(i2, 1, 0);
        atom* pBondAtom2 = pMolecule->getAtom(i3, 1, 0);

        double bondDist = pBondAtom1->getCoords()->dist(*pBondAtom2->getCoords());
        Bond* pBond = pMolecule->addBond(pBondAtom1, pBondAtom2, 0, 0, 0, bondDist);
        if (pBond) {
          pBondAtom1->addBondedAtom(pBondAtom2);
          pBondAtom2->addBondedAtom(pBondAtom1);
        }
      }
    }
    iac.close();
}

// ==============================================================
// Function : Write
// --------------------------------------------------------------
// Write an AMBER ac file.
// ==============================================================
void acParser::Write(const std::string &acfile, molecule* pMolecule)
{
    std::ofstream oac;
    oac.open(acfile.c_str());

    if (!oac or (pMolecule == 0)) {
      std::cout << "\nUNABLE TO OPEN AC FILE"
                << "\nFILENAME = " << acfile << std::endl;
      return;
    }

    std::vector<atom*> atomList = pMolecule->getAtomList();
    char temp[80];
    atom* pAtom = 0;
    submolecule* pSubMol = 0;
    stdAtom* pStdAtom = 0;
    std::string element;
    std::string name;

    for (unsigned int a = 0; a < atomList.size(); a++) {
      pAtom = atomList[a];
      pSubMol = pAtom->getParent();
      pStdAtom = pAtom->getStdAtom();
      element = pAtom->getElementSymbol().c_str();
      name = pAtom->getName().c_str();

      // It can only get the element by matching against the
      // name anyway. So we need to make sure that if the
      // element name is two characters, the second letter
      // of the atom name is lowercase.
      if (element.length() == 2) {
        name[1] = tolower(name[1]);
      }

      sprintf(temp,"ATOM  %5d %-4.4s %-3.3s  %4d    %8.3f%8.3f%8.3f %9.6f        %2s",
              pAtom->getIndex(),
              (name.c_str()),
              (pSubMol->getName().c_str()),(pSubMol->getSubMolId()),
              pAtom->getX(),pAtom->getY(),pAtom->getZ(),
              0.0, pStdAtom->type.c_str());
              //pAtom->getX(),pAtom->getY(),pAtom->getZ(),pAtom->getESP(), pStdAtom->type);
      oac << temp << std::endl;
    }
    oac.close();
}

// ==============================================================
// Function : Write
// --------------------------------------------------------------
// Write a AMBER ac file.
// ==============================================================
void acParser::Write(const std::string &acfile, collection* pCollection)
{
    std::ofstream oac;
    oac.open(acfile.c_str());

    if (!oac or (pCollection == 0)) {
      std::cout << "\nUNABLE TO OPEN AC FILE"
                << "\nFILENAME = " << acfile << std::endl;
      return;
    }
    char temp[80];
    double dTemp = 0.0;
    submolecule* pSubMol = 0;
    atom* pAtom = 0;
    stdAtom* pStdAtom = 0;
    Bond* pBond = 0;
    typedef std::map<int, Bond*>::iterator BondMapIterator;
    int bondIndex = 1;

    std::string element = "";
    std::string name = "";

    double charge = 0.0;
    std::vector<molecule*> molList = pCollection->getMoleculeList();
    for (unsigned int m = 0; m < molList.size(); m++) {
      std::vector<submolecule*> subMolList = molList[m]->getSubMoleculeList();
      for (unsigned int s = 0; s < subMolList.size(); s++) {
        if (subMolList[s]->hasStdFrag()) {
          charge += double(subMolList[s]->getStdFrag()->getCharge());
        }
      }
    }
    sprintf(temp,"CHARGE     %3.2f", charge);
    oac << temp << std::endl;

    for (unsigned int m = 0; m < molList.size(); m++) {
      std::vector<atom*> atomList = molList[m]->getAtomList();
      for (unsigned int a = 0; a < atomList.size(); a++) {
        pAtom = atomList[a];
        pSubMol = pAtom->getParent();

        element = pAtom->getElementSymbol().c_str();
        name = pAtom->getName().c_str();

        // It can only get the element by matching against the
        // name anyway. So we need to make sure that if the
        // element name is two characters, the second letter
        // of the atom name is lowercase.
        if (element.length() == 2) {
          name[1] = tolower(name[1]);
        }

        if (pAtom->getStdAtom()) {
          pStdAtom = pAtom->getStdAtom();
          sprintf(temp,"ATOM  %5d %-4.4s %-3.3s  %4d    %8.3f%8.3f%8.3f %9.6f        %2s",
                  pAtom->getColIndex(),
                  (name.c_str()),
                  (pSubMol->getName().c_str()),(pSubMol->getSubMolId()),
                  pAtom->getX(),pAtom->getY(),pAtom->getZ(),dTemp, pStdAtom->type.c_str());
          oac << temp << std::endl;
        }
      }
      std::map<int, Bond*> molBonds = molList[m]->getBondMap();
      if (!molBonds.empty()) {
        for (BondMapIterator b = molBonds.begin(); b != molBonds.end(); b++) {
          pBond = b->second;
          atom* pAt1 = pBond->atom1;
          atom* pAt2 = pBond->atom2;

          sprintf(temp,"BOND %4d %4d %4d %4d  %4s %4s", bondIndex,
            pAt1->getColIndex(), pAt2->getColIndex(), pBond->type,
            pAt1->getName().c_str(), pAt2->getName().c_str());
          oac << temp << std::endl;
          bondIndex++;
        }
      }
      oac.close();
    }
}

} // MTKpp namespace
