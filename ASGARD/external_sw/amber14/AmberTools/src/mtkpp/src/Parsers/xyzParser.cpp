/*!
   \file xyzParser.cpp
   \brief Parses XYZ files
   \author Martin Peters

   $Date: 2010/03/29 20:39:35 $
   $Revision: 1.11 $

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


#include "xyzParser.h"

#include "StringManip.h"

#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/submolecule.h"
#include "Molecule/atom.h"
#include "Molecule/bond.h"
#include "Molecule/element.h"

namespace MTKpp
{

// ============================================================
// Function : XyzParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
xyzParser::xyzParser():baseParser() {}

// =========================================================
// Function : XyzParser()
// ---------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// =========================================================
xyzParser::~xyzParser() {}

// =========================================================
// Function : Read
// ---------------------------------------------------------
// parsers a xyz file
// ---------------------------------------------------------
// Format :
// 1st line : Number of atoms(N)
// 2nd line : Title
// 3rd line : atom x y z
// 4th line : atom x y z
//   .
//   .
// Nth line : atom x y z
// EOF
// =========================================================

void xyzParser::Read(const std::string &xyzfile, collection* pCollection)
{
    std::ifstream ixyz;
    ixyz.open(xyzfile.c_str());

    if (!ixyz) {
      std::cout << "\nUNABLE TO OPEN XYZ FILE"
                << "\nFILENAME = " << xyzfile
                << "\nEXITING...\n" << std::endl;
      return;
    }

    molecule* pMolecule  = pCollection->addMolecule();
    pMolecule->setName("Mol");
    pMolecule->setMolId(pCollection->getNumberMolecules());

    submolecule* pSubMolecule = pMolecule->addSubMolecule();
    pSubMolecule->setName("sMol");
    pSubMolecule->setSubMolId(pMolecule->getNumSubMolecules());

    int end   = xyzfile.length();
    int slash = xyzfile.find_last_of("/");
    std::string file_name = xyzfile.substr(slash+1,(end-slash-5));
    pMolecule->setName(file_name);

    std::string fileline;
    int natoms = 0;
    std::string title;
    std::vector<std::string> splitstring;
    std::string atom;
    double x,y,z;

    // FIRST LINE:
    getline(ixyz,fileline);
    //if (strlen(fileline.c_str()) == 0) {
    if (fileline.length() == 0) {
      return;
    }
    else {
      natoms = atoi(fileline.c_str() );
    }

    if (natoms == 0) {
      return;
    }

    // SECOND LINE:
    getline(ixyz,fileline);
    //if (strlen(fileline.c_str()) == 0)
    if (fileline.length() == 0)
      title = "";
    else
      title = fileline;

    // THIRD LINE:
    // LOOP OVER NATOMS AND GET ATOM, X, Y, Z.
    unsigned int counter = 0;
    for (int n = 1; n <= natoms; n ++) {
      if (getline(ixyz,fileline)) {
        std::string::size_type lastPos = fileline.find_first_not_of(" ",0);
        std::string::size_type pos = fileline.find_first_of(" ",lastPos);

        while (std::string::npos != pos || std::string::npos != lastPos) {
          splitstring.push_back(fileline.substr(lastPos,pos - lastPos));
          lastPos = fileline.find_first_not_of(" ", pos);
          pos = fileline.find_first_of(" ", lastPos);
        }
        if (splitstring.size() == (counter + 4)) {
          atom = splitstring[counter];
          char *endptr;
          x = strtod((char*)splitstring[counter+1].c_str(),&endptr);
          y = strtod((char*)splitstring[counter+2].c_str(),&endptr);
          z = strtod((char*)splitstring[counter+3].c_str(),&endptr);
          pAtom = pSubMolecule->addAtom();
          pAtom->setElement(pCollection->pElements->getElement(atom));
          pAtom->setCoords(x,y,z);
          pAtom->setFileID(n);
          counter=counter+4;
        }
        else {
          return;
        }
      }
      else {
        return;
      }
    }
    ixyz.close();
}

// ==============================================================
// Function : Write
// --------------------------------------------------------------
// Write a XYZ file.
// ==============================================================
void xyzParser::Write(const std::string &xyzfile, molecule* pMolecule)
{
    std::ofstream oxyz;
    oxyz.open(xyzfile.c_str());

    if (!oxyz or (pMolecule == 0)) {
      std::cout << "\nUNABLE TO OPEN XYZ FILE"
                << "\nFILENAME = " << xyzfile << std::endl;
      return;
    }

    typedef std::vector<atom*>::iterator AtomIterator;
    std::vector<atom*> Atoms;
    Atoms  = pMolecule->getAtomList();

    oxyz << Atoms.size() << "\nCOMMENT" << std::endl;
    double curX = 0.0;
    double curY = 0.0;
    double curZ = 0.0;

    for (AtomIterator d=Atoms.begin(); d != Atoms.end(); d++) {
      pAtom = *d;
      //char temp[100];
      curX = pAtom->getX();
      curY = pAtom->getY();
      curZ = pAtom->getZ();

      oxyz << std::setw(2) << pAtom->getElement()->symbol << " " << std::setiosflags(std::ios::fixed) << std::setprecision(6);
      oxyz << std::setw(12) << curX << " ";
      oxyz << std::setw(12) << curY<< " ";
      oxyz << std::setw(12) << curZ;
      oxyz << std::endl;
    }
    oxyz.close();
}

} // MTKpp namespace

