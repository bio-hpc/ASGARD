/*!
   \file molParser.cpp
   \brief Parses mol files
   \author Martin Peters

   Reads and writes mol files

   $Date: 2010/03/29 20:39:35 $
   $Revision: 1.15 $

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

#include "molParser.h"

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

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

// ============================================================
// Function : molParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
molParser::molParser():baseParser() {
    pMolecule = 0;
    pSubMolecule  = 0;
    pBond = 0;
    pBondAtom1 = 0;
    pBondAtom2 = 0;
    pConnections = 0;
}

// ============================================================
// Function : ~molParser()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
molParser::~molParser() {
    delete pConnections;
}

// ============================================================
// Function : Read
// ------------------------------------------------------------
// parsers a mol file
// ------------------------------------------------------------
// Format:
// 1st line       : Compound Name
// 2nd line       : initials, program, date, dimension, scale factors, energy, reg name
//                  IIPPPPPPPPMMDDYYHHmmddSSssssssssssEEEEEEEEEEEERRRRRR
//                  User's first and last initials (l), program name (P),
//                  date/time (M/D/Y,H:m), dimensional codes (d),
//                  scaling factors (S, s), energy (E) if modeling program input,
//                  internal registry number (R) if input through MDL form.
// 3rd line       : Comments
// 4rd line       : Number of atoms(N default limit = 255),Number of bonds(B),
//                  atlist,obs,chiral, nentries,obs,obs,obs,obs,
//                  nprop (default value = 999)
//                  aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv
// 5th line       : x y z atom  massdiff charge stereo hcount stcare valence 0 0 0 0 0 0
//    .
//    .
// (N+4)th line   : x y z atom charge charge stereo hcount stcare valence 0 0 0 0 0 0
// (N+5)th line   : at1 at2 bondtype bondstereo 0  0  0
//    .
//    .
// (N+4+B)th line : at1 at2 bondtype bondstereo 0  0  0
// M  END
// ============================================================
void molParser::Read(const std::string &molfile, molecule* pMolecule, const bool &bohr)
{
    collection* pCollection = pMolecule->getParent();
    pConnections = new connections(pCollection);
    std::ifstream imol;
    imol.open(molfile.c_str());

    if (!imol) {
      std::cout << "\nUNABLE TO OPEN MOL FILE"
                << "\nFILENAME = " << molfile << std::endl;
      return;
    }

    std::string fileline;
    std::string name;
    std::string initials="XX", program="",date=" ",dim=" ";
    std::string comment;
    int natoms=0, nbonds=0;
    std::string v2="V2000";
    std::string atom;

    AtomLine* pthisAtom = new AtomLine;

    pMolecule->setName("Mol");
    pMolecule->setMolId(pCollection->getNumberMolecules());
    pMolecule->inFileType = "mol";

    pSubMolecule = pMolecule->addSubMolecule();
    pSubMolecule->setName("XXX");
    pSubMolecule->setSubMolId(pMolecule->getNumSubMolecules());

    // - FIRST LINE: - //
    getline(imol,fileline);

    // - SECOND LINE: - //
    getline(imol,fileline);

    // - THIRD LINE: - //
    getline(imol,fileline);

    // - FOURTH LINE: - //
    getline(imol,fileline);
    if (fileline.length() == 0) {
      std::cout << "NOT A VALID MOL FILE" << std::endl;
      std::cout << "line 4 of a mol file should have 38 characters" << std::endl;
      std::cout << "Number of atoms(N), Number of bonds(B), atlist, obs, chiral, nentries, obs, obs, obs, obs, nprop"
                << std::endl;
      std::cout << "x y z atom  massdiff charge stereo hcount stcare valence 0 0 0 0 0 0" << std::endl;
      return;
    }
    else if (! (fileline.length() >= 6)) {
      std::cout << "line 4 of a mol file could have 38 characters" << std::endl;
      std::cout << "But must have at least 6, Natoms(3),Nbonds(3)" << std::endl;
    }
    else {
      natoms = atoi((fileline.substr( 0,3)).c_str());
      pSubMolecule->setNumAtoms(natoms);
      nbonds = atoi((fileline.substr( 3,3)).c_str());
      pSubMolecule->setNumBonds(nbonds);
      //atlist = atoi((fileline.substr( 6,3)).c_str());
      //obs = atoi((fileline.substr( 9,3)).c_str());
      //chiral = atoi((fileline.substr( 12,3)).c_str());
      //nentries = atoi((fileline.substr( 15,3)).c_str());
      //nprop = atoi((fileline.substr( 30,3)).c_str());
      //aaa bbb lll fff ccc sss xxxrrrpppiii mmmvvvvvv
    }

    // - FIFTH LINE:
    // - LOOP OVER NATOMS AND GET ATOM, X, Y, Z.
    for (int n = 1; n <= natoms; n++) {
      getline(imol,fileline);
      if (ReadAtomLine(fileline,pthisAtom)) {
        pAtom = pSubMolecule->addAtom();

        pAtom->setElement(pCollection->pElements->getElement(pthisAtom->element));
        pAtom->setFileID(n);
        pMolecule->setMaxFileID(n);

        if (bohr) {
          pAtom->setCoords(pthisAtom->x/ANG2BOHR, pthisAtom->y/ANG2BOHR,
                           pthisAtom->z/ANG2BOHR);
        }
        else {
          pAtom->setCoords(pthisAtom->x,pthisAtom->y,pthisAtom->z);
        }
      }
      else {
        std::cout << "  ERROR READING ATOM LINE IN MOL2 PARSER " << std::endl;
        return;
      }
    }

    // - (N+1)th LINE:
    // - LOOP OVER NBONDS AND GET AT1, AT2 AND BONDTYPE.
    int at1 = 0, at2 = 0, bondType = 0, bondStereo = 0, bondTopology = 0;
    for (int n = 1; n <= nbonds; n ++) {
      getline(imol,fileline);
      if (fileline.length() < 9) {
        std::cout << " NOT A VALID MOL FILE " << std::endl;
        std::cout << "line (NATOMS+5) of a mol file should have at least 9 characters" << std::endl;
        return;
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

        if (pBondAtom1 and pBondAtom2) {
          pBond = 0;
          pBond = pMolecule->addBond(pBondAtom1, pBondAtom2, bondType, bondStereo, bondTopology, 0.0);
          if (pBond) {
            pBondAtom1->addBondedAtom(pBondAtom2);
            pBondAtom2->addBondedAtom(pBondAtom1);
//std::cout << " Adding bond " << pBondAtom1->getIndex() << "-" << pBondAtom2->getIndex() << std::endl;
          }
        }
        bondType = 0, bondStereo = 0, bondTopology = 0;
      }
    }
    for (int n = 1; n <= 999; n++) {
      getline(imol,fileline);
      if ((fileline.substr(0,6)) == "M  END" ) {
        break;
      }
      if ((fileline.substr(0,6)) == "M  CHG") {
        int nCHG = atoi((fileline.substr(6,3)).c_str());
        int st = 10;
        int a,c;
        for (int j = 0; j < nCHG; j++) {
          a = atoi((fileline.substr(st,4)).c_str());
          c = atoi((fileline.substr(st+4,4)).c_str());
          pAtom = pMolecule->getAtom(a, 1, 0);
          pAtom->setFormalCharge(c);
          st+=8; // why is this 8 and in sdfParser its 4???
        }
      }
    }

    // - Clean Up and update - //
    delete pthisAtom;
    imol.close();

    pMolecule->setKind(2); // drug/small molecule

    // With a mol file, I assume that all bonds and bond types are defined correctly
    pConnections->assignAngles(pMolecule);
    pConnections->assignTorsions(pMolecule);
    pConnections->assignImpropers(pMolecule);

    pMolecule->bBondsAssigned = true;
    pMolecule->bAnglesAssigned = true;
    pMolecule->bTorsionsAssigned = true;
    pMolecule->bImpropersAssigned = true;

    //pMolecule->determineHybridizations(2);
    //pMolecule->determineRings(); // Needs bond types.  Sets atom types, hybridizations for ring systems
    //pMolecule->determineValences(); // needs bond types

    //pMolecule->addHydrogens(); // needs atom types, ring atoms, atom hybridizations.  Updates atom valences

    return;
}

// ============================================================
// Function : Read
// ------------------------------------------------------------
// parsers a mol file
// ============================================================
void molParser::Read(const std::string &molfile, collection* pCollection, const bool &bohr)
{
    pConnections = new connections(pCollection);
    std::ifstream imol;
    imol.open(molfile.c_str());

    if (!imol) {
      std::cout << "\nUNABLE TO OPEN MOL FILE"
                << "\nFILENAME = " << molfile << std::endl;
      return;
    }
    imol.close();

    pMolecule  = pCollection->addMolecule();
    this->Read(molfile, pMolecule, bohr);
}

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write a MOL file.
// ------------------------------------------------------------
// Format:
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
void molParser::Write(const std::string &molfile, molecule* pMolecule,
                      std::vector< vector3d > &coordinates)
{
    std::ofstream omol;
    omol.open(molfile.c_str());

    if (!omol or (pMolecule == 0)) {
      std::cout << "\nUNABLE TO OPEN MOL FILE"
                << "\nFILENAME = " << molfile << std::endl;
      return;
    }

    // ATOMS
    typedef std::vector<atom*>::iterator AtomIterator;
    std::vector<atom*>::iterator atomIterator1;
    std::vector<atom*>::iterator atomIterator2;
    std::vector<atom*>                   Atoms;
    Atoms  = pMolecule->getAtomList();
    int n_atoms = Atoms.size();

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

    std::string mol_name = pMolecule->getName();

    if (n_atoms > 999) {
      std::cout << "ERROR WRITING MOL FILE" << std::endl;
      std::cout << "NUMBER OF ATOMS CANNOT EXCEED 999" << std::endl;
      return;
    }

    omol << mol_name <<  "\n Molecule \n" << std::endl;
    omol << std::setw(3) << n_atoms << std::setw(3) << nBonds << "  0  0  1  0  0  0  0  0999 V2000" << std::endl;

    // ============================================================
    // Write the file,
    // Loop over atoms.
    // ============================================================

    int n = 1;

    // - Loop over atoms - //
    int atomIndex = 0;
    double curX = 0.0;
    double curY = 0.0;
    double curZ = 0.0;

    for (AtomIterator d = Atoms.begin(); d != Atoms.end(); d++) {
      pAtom = *d;
      char temp[100];
      curX = coordinates[atomIndex].getX();
      curY = coordinates[atomIndex].getY();
      curZ = coordinates[atomIndex].getZ();

      if (!pAtom->getElement()) {
        std::cout << " Can't find element symbol ... exiting " << std::endl;
        //exit(0);
        throw MTKException(" Can't find element symbol ... exiting ");
      }
      std::string elSymbol = pAtom->getElement()->symbol;

      sprintf(temp,"%10.4f%10.4f%10.4f %-3s%2d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d",
       curX, curY, curZ,(elSymbol.c_str()),
       0,0,0,0,0,0,0,0,0,0,0,0);
       omol << temp << std::endl;
      // - Create a map to ensure that the numbering of the mol file is correct - //
      //itsAtomMap[pAtom->getFileID()] = n;
      itsAtomMap[pAtom->getIndex()] = n;
      n = n+1;
      atomIndex++;
    }

    // - Loop over bonds - //
    if (!moleculeBondMap.empty()) {
      for (BondMapIterator b = moleculeBondMap.begin(); b != moleculeBondMap.end(); b++) {
        pBond = b->second;

        atomIterator1 = std::find(Atoms.begin(), Atoms.end(), pBond->atom1);
        atomIterator2 = std::find(Atoms.begin(), Atoms.end(), pBond->atom2);

        if ((atomIterator1 != Atoms.end()) and (atomIterator2 != Atoms.end())) {
          char temp[100];
          int bondType = pBond->type;
          if (bondType == 6) bondType = 1;
          if (bondType == 7) bondType = 2;
          sprintf(temp,"%3d%3d%3d%3d%3d%3d%3d",
           itsAtomMap[pBond->atom1->getIndex()], itsAtomMap[pBond->atom2->getIndex()],
           bondType,
           pBond->stereo,
           pBond->topology,
           0,0);
           omol << temp << std::endl;
        }
      }
    }
    int fc = 0;
    for (AtomIterator d=Atoms.begin(); d != Atoms.end(); d++) {
      pAtom = *d;
      fc = pAtom->getFormalCharge();
      if (fc != 0) {
        char temp[100];
        sprintf(temp,"%6s%3d%4d%4d", "M  CHG", 1,
        itsAtomMap[pAtom->getIndex()], fc);
        omol << temp << std::endl;
      }
    }
    omol << "M  END" << std::endl;
    omol.close();
}

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write a MOL file.
// ============================================================
void molParser::Write(const std::string &molfile, molecule* pMolecule)
{
    if (pMolecule != 0) {
      std::vector< vector3d > coordinates;
      pMolecule->getCoordinates(coordinates);
      this->Write(molfile, pMolecule, coordinates);
    }
}

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write a MOL file.
// ============================================================
void molParser::Write(const std::string &molfile,collection* pCollection, const int &molId)
{
    molecule* pMolecule = pCollection->getMolecule(molId);

    if (pMolecule != 0) {
      std::vector< vector3d > coordinates;
      pMolecule->getCoordinates(coordinates);
      this->Write(molfile, pMolecule, coordinates);
    }
}

// ============================================================
// Function : ReadAtomLine
// ------------------------------------------------------------
// Reads ATOM AND HETATOM line and returns data
// ============================================================

bool molParser::ReadAtomLine(std::string &fileline, AtomLine* pthisAtom)
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

// ============================================================
// Function : getAtom
// ------------------------------------------------------------
// returns an atom of given index
// ============================================================

int molParser::getAtom(int index)
{
    std::map<int, int>::iterator iter;
    iter = itsAtomMap.find(index);

    if (iter != itsAtomMap.end()) {
         return itsAtomMap[index];
    }
    return 0;
}

} // MTKpp namespace

