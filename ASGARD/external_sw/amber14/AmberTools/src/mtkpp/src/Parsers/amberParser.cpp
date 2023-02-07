/*!
   \file amberParser.cpp
   \brief Parses amber prmtop/crd files
   \author Martin Peters

   Reads and writes amber prmtop/crd files

   $Date: 2010/04/29 19:06:19 $
   $Revision: 1.7 $

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

#include "amberParser.h"
#include "Molecule/element.h"
#include <vector>
#include <algorithm>

#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/submolecule.h"
#include "Molecule/atom.h"
#include "Utils/vector3d.h"
#include "Molecule/parameters.h"
#include "Molecule/bond.h"
#include "Molecule/angle.h"
#include "Molecule/torsion.h"
#include "Molecule/improper.h"
#include "Molecule/metalCenter.h"
#include "Molecule/stdFrag.h"
#include "time.h"

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

// ============================================================
// Function : amberParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
amberParser::amberParser() {}

// =========================================================
// Function : amberParser()
// ---------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// =========================================================
amberParser::~amberParser() {}

// ============================================================================
// Function : Write
// ----------------------------------------------------------------------------
// Writes an amber prmtop file
// ============================================================================
void amberParser::Write(const std::string &inpcrd, const std::string &prmtop, collection* pCollection)
{
    std::ofstream oinpcrd(inpcrd.c_str());
    std::ofstream oprmtop;
    oprmtop.open(prmtop.c_str());

    if (!oinpcrd or (pCollection == 0) or !oprmtop ) {
      std::cout << "\nUNABLE TO OPEN PRMTOP or CRD FILE"
                << "\nFILENAME = " << inpcrd
                << " OR FILENAME = " << prmtop << std::endl;
      return;
    }

    molecule* pMolecule = 0;
    typedef std::vector<molecule*>::iterator moleculeIterator;
    std::vector<molecule*> molList = pCollection->getMoleculeList();
    std::vector<metalCenter*> metalCenterList = pCollection->getMetalCenters();

    std::vector<submolecule*> smolList;
    std::vector<Bond*> bondList;
    std::vector<Angle*> angleList;
    std::vector<Torsion*> torsionList;
    std::vector<Improper*> improperList;

    typedef std::map<int, Bond*>::iterator BondMapIterator;
    typedef std::map<ULONG_KIND, Angle*>::iterator AngleMapIterator;
    typedef std::map<ULONG_KIND, Torsion*>::iterator TorsionMapIterator;
    typedef std::map<int, Improper*>::iterator ImproperMapIterator;
    typedef std::vector<torsionParam*>::iterator torsionParamIterator;
    typedef std::vector<improperParam*>::iterator improperParamIterator;

    // Get number of unique types, bonds, angles and dihedrals
    unsigned int nUniqueAtomTypes = pCollection->getNumUniqueAtomTypes();
    unsigned int nUniqueBondTypes = pCollection->getNumUniqueBondTypes();
    unsigned int nUniqueAngleTypes = pCollection->getNumUniqueAngleTypes();
    unsigned int nUniqueDihedralTypes = pCollection->getNumUniqueDihedralTypes();

    // Get number of bonds, angles and dihedrals
    unsigned int nBondsWithH = pCollection->getNumBondsWithH();
    unsigned int nBondsWithOutH = pCollection->getNumBondsWithOutH();
    unsigned int nAnglesWithH = pCollection->getNumAnglesWithH();
    unsigned int nAnglesWithOutH = pCollection->getNumAnglesWithOutH();
    unsigned int nDihedralsWithH = pCollection->getNumDihedralsWithH();
    unsigned int nDihedralsWithOutH = pCollection->getNumDihedralsWithOutH();

    for (moleculeIterator c = molList.begin(); c != molList.end(); c++) {
      pMolecule = *c;
      std::vector<submolecule*> lList = pMolecule->getSubMoleculeList();
      for (unsigned int s = 0; s < lList.size(); s++) {
        smolList.push_back(lList[s]);
      }

      std::map<int, Bond*> bMap =  pMolecule->getBondMap();
      for (BondMapIterator b = bMap.begin(); b != bMap.end(); b++) {
        Bond* pBond = b->second;
        bondList.push_back(pBond);
      }

      std::map<ULONG_KIND, Angle*> aMap =  pMolecule->getAngleMap();
      for (AngleMapIterator b = aMap.begin(); b != aMap.end(); b++) {
        Angle* pAngle = b->second;
        angleList.push_back(pAngle);
      }

      std::map<ULONG_KIND, Torsion*> tMap =  pMolecule->getTorsionMap();
      for (TorsionMapIterator b = tMap.begin(); b != tMap.end(); b++) {
        Torsion* pTorsion = b->second;
        torsionList.push_back(pTorsion);
      }

      std::map<int, Improper*> iMap =  pMolecule->getImproperMap();
      for (ImproperMapIterator b = iMap.begin(); b != iMap.end(); b++) {
        Improper* pImproper = b->second;
        improperList.push_back(pImproper);
      }
    }

    for (unsigned int m = 0; m < metalCenterList.size(); m++) {
      metalCenter* metCen = metalCenterList[m];
      std::map<int, Bond*> bMap =  metCen->getBondMap();
      for (BondMapIterator b = bMap.begin(); b != bMap.end(); b++) {
        Bond* pBond = b->second;
        bondList.push_back(pBond);
      }

      std::map<ULONG_KIND, Angle*> aMap =  metCen->getAngleMap();

      for (AngleMapIterator b = aMap.begin(); b != aMap.end(); b++) {
        Angle* pAngle = b->second;
        angleList.push_back(pAngle);
        if (!pAngle->pAngleParam) {
          std::cout << " ANGLE PARAM NOT FOUND " << std::endl;
        }
      }
    }

    unsigned int nResidues = smolList.size();

    std::map<atom*, int> localAtomIndex;
    std::vector<atom*> atomList = pCollection->getAtomList();
    unsigned int nAtoms = atomList.size();

    // Calculate the number of excluded atoms for each atom
    int nExcludedAtoms = 0;
    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAt1 = atomList[i];
      localAtomIndex[pAt1] = i + 1;
      int nE = 0;
      for (unsigned int j = i+1; j < nAtoms; j++) {
        atom* pAt2 = atomList[j];
        if (pAt1->hasBondedAtom(pAt2) or
            pAt1->has13BondedAtom(pAt2) or
            pAt1->has14BondedAtom(pAt2)) {
          nE++;
        }
      }
      if (nE == 0) nE = 1;
      nExcludedAtoms += nE;
    }

    //oprmtop << "%VERSION" << std::endl;
    oprmtop << "%VERSION  VERSION_STAMP = V0001.000  DATE = 01/05/06  15:09:43" << std::endl;
    oprmtop << "%FLAG TITLE" << std::endl;
    oprmtop << "%FORMAT(20a4)" << std::endl;
    oprmtop << std::left << std::setw(80) << pCollection->getName() << std::endl;
    oinpcrd << std::left << std::setw(80) << pCollection->getName() << std::endl;
    oprmtop.unsetf(std::ios::left);
    oinpcrd.unsetf(std::ios::left);

    // Write crd file
    /*
      FORMAT(I5,5E15.7) NATOM,TIME
       NATOM  : total number of atoms in coordinate file
       TIME   : option, current time in the simulation (picoseconds)
    */
    oinpcrd << std::setw(5) << nAtoms << std::endl;

    /*
      FORMAT(6F12.7) (X(i), Y(i), Z(i), i = 1,NATOM)
       X,Y,Z  : coordinates
    */
    int index = 1;
    for (unsigned int i = 0; i < nAtoms; i++) {
      vector3d& coords = (*atomList[i]->getCoords());
      for (unsigned int j = 0; j < 3; j++) {
        if (index < 6) {
          oinpcrd << std::fixed << std::setprecision(7) << std::setw(12) << coords[j];
        }
        else {
          oinpcrd << std::fixed << std::setprecision(7) << std::setw(12)
                  << coords[j] << std::endl;
          index = 0;
        }
        index++;
      }
    }
    oprmtop.unsetf(std::ios::fixed);
    oinpcrd.close();

    /*
      NATOM  : total number of atoms
      NTYPES : total number of distinct atom types
      NBONH  : number of bonds containing hydrogen
      MBONA  : number of bonds not containing hydrogen
      NTHETH : number of angles containing hydrogen
      MTHETA : number of angles not containing hydrogen
      NPHIH  : number of dihedrals containing hydrogen
      MPHIA  : number of dihedrals not containing hydrogen
      NHPARM : currently not used
      NPARM  : currently not used
    */
    oprmtop << "%FLAG POINTERS" << std::endl;
    oprmtop << "%FORMAT(10I8)" << std::endl;
    oprmtop << std::setw(8) << nAtoms
            << std::setw(8) << nUniqueAtomTypes
            << std::setw(8) << nBondsWithH
            << std::setw(8) << nBondsWithOutH
            << std::setw(8) << nAnglesWithH
            << std::setw(8) << nAnglesWithOutH
            << std::setw(8) << nDihedralsWithH
            << std::setw(8) << nDihedralsWithOutH
            << std::setw(8) << 0
            << std::setw(8) << 0
            << std::endl;

    /*
      NEXT   : number of excluded atoms
      NRES   : number of residues
      NBONA  : MBONA + number of constraint bonds
      NTHETA : MTHETA + number of constraint angles
      NPHIA  : MPHIA + number of constraint dihedrals
      NUMBND : number of unique bond types
      NUMANG : number of unique angle types
      NPTRA  : number of unique dihedral types
      NATYP  : number of atom types in parameter file, see SOLTY below
      NPHB   : number of distinct 10-12 hydrogen bond pair types
    */
    oprmtop << std::setw(8) << nExcludedAtoms
            << std::setw(8) << nResidues
            << std::setw(8) << (nBondsWithOutH + 0)
            << std::setw(8) << (nAnglesWithOutH + 0)
            << std::setw(8) << (nDihedralsWithOutH + 0)
            << std::setw(8) << nUniqueBondTypes
            << std::setw(8) << nUniqueAngleTypes
            << std::setw(8) << nUniqueDihedralTypes
            << std::setw(8) << nUniqueAtomTypes
            << std::setw(8) << 0 << std::endl;

    /*
      IFPERT : set to 1 if perturbation info is to be read in
      NBPER  : number of bonds to be perturbed
      NGPER  : number of angles to be perturbed
      NDPER  : number of dihedrals to be perturbed
      MBPER  : number of bonds with atoms completely in perturbed group
      MGPER  : number of angles with atoms completely in perturbed group
      MDPER  : number of dihedrals with atoms completely in perturbed groups
      IFBOX  : set to 1 if standard periodic box, 2 when truncated octahedral
      NMXRS  : number of atoms in the largest residue
      IFCAP  : set to 1 if the CAP option from edit was specified
    */
    oprmtop << std::setw(8) << 0
            << std::setw(8) << 0
            << std::setw(8) << 0
            << std::setw(8) << 0
            << std::setw(8) << 0
            << std::setw(8) << 0
            << std::setw(8) << 0
            << std::setw(8) << 0
            << std::setw(8) << pCollection->largestResidueSize()
            << std::setw(8) << 0 << std::endl;

    oprmtop << std::setw(8) << 0 << std::endl;

    /*
      FORMAT(20a4)  (IGRAPH(i), i=1,NATOM)
      IGRAPH : the user atoms names
    */
    oprmtop << "%FLAG ATOM_NAME" << std::endl;
    oprmtop << "%FORMAT(20a4)" << std::endl;

    index = 1;
    std::string lStr = "";
    for (unsigned int k = 0; k < nAtoms; k++) {
      lStr += atomList[k]->getName();
      if (index < 20) {
        if (k == nAtoms-1) {
          oprmtop << std::left << std::setw(80) << lStr << std::endl;
          oprmtop.unsetf(std::ios::left);
          lStr = "";
        }
      }
      else {
        oprmtop << std::setw(80) << lStr << std::endl;
        lStr = "";
        index = 0;
      }
      index++;
    }
    oprmtop.flush();

    /*
      FORMAT(5E16.8)  (CHRG(i), i=1,NATOM)
      CHRG   : the atom charges.  (Divide by 18.2223 to convert to charge
               in units of the electron charge)
    */
    oprmtop << "%FLAG CHARGE" << std::endl;
    oprmtop << "%FORMAT(5E16.8)" << std::endl;
    index = 1;
    for (unsigned int k = 0; k < nAtoms; k++) {
      stdAtom* pStdAtom = atomList[k]->getStdAtom();
      if (pStdAtom) {
        double chg = pStdAtom->atmCharge;
        if (index < 5) {
          oprmtop << std::uppercase << std::scientific << std::setprecision(8)
                  << std::setw(16) << chg * E2KCAL;
        }
        else {
          oprmtop << std::uppercase << std::scientific << std::setw(16)
                  << std::setprecision(8) << chg * E2KCAL << std::endl;
          index = 0;
        }
        index++;
      }
    }
    if (index != 1) {
      oprmtop << " " << std::endl;
      oprmtop.flush();
    }
    oprmtop.unsetf(std::ios::scientific);

    /*
      FORMAT(5E16.8)  (AMASS(i), i=1,NATOM)
      AMASS  : the atom masses
    */
    oprmtop << "%FLAG MASS" << std::endl;
    oprmtop << "%FORMAT(5E16.8)" << std::endl;
    index = 1;
    for (unsigned int k = 0; k < nAtoms; k++) {
      double mass = atomList[k]->getElement()->mass;
      if (index < 5) {
        oprmtop << std::uppercase << std::scientific << std::setprecision(8)
                << std::setw(16) << mass;
      }
      else {
        oprmtop << std::uppercase << std::scientific << std::setw(16)
                << std::setprecision(8) << mass << std::endl;
        index = 0;
      }
      index++;
    }
    if (index != 1) {
      oprmtop << " " << std::endl;
      oprmtop.flush();
    }
    oprmtop.unsetf(std::ios::scientific);

    /*
      FORMAT(12I6)  (IAC(i), i=1,NATOM)
      IAC    : index for the atom types involved in Lennard Jones (6-12)
               interactions.  See ICO below.
    */
    std::map<std::string, int> atomTypesMap;
    std::vector<std::string> atomTypesUsed = pCollection->getUniqueAtomTypes();
    for (unsigned int k = 0; k < atomTypesUsed.size(); k++) {
      atomTypesMap[atomTypesUsed[k]] = k;
    }

    oprmtop << "%FLAG ATOM_TYPE_INDEX" << std::endl;
    oprmtop << "%FORMAT(10I8)" << std::endl;
    index = 1;
    for (unsigned int k = 0; k < atomList.size(); k++) {
      stdAtom* pStdAtom = atomList[k]->getStdAtom();
      if (pStdAtom) {
        std::string t = pStdAtom->type;
        if (index < 10) {
          if (k != nAtoms-1) {
            oprmtop << std::setw(8) << atomTypesMap[t] + 1;
          }
          else {
            oprmtop << std::setw(8) << atomTypesMap[t] + 1 << std::endl;
          }
        }
        else {
          oprmtop << std::setw(8) << atomTypesMap[t] + 1 << std::endl;
          index = 0;
        }
        index++;
      }
    }

    /*
      FORMAT(12I6)  (NUMEX(i), i=1,NATOM)
      NUMEX  : total number of excluded atoms for atom "i".  See
               NATEX below.
    */
    oprmtop << "%FLAG NUMBER_EXCLUDED_ATOMS" << std::endl;
    oprmtop << "%FORMAT(10I8)" << std::endl;
    index = 1;
    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAt1 = atomList[i];
      int nE = 0;
      for (unsigned int j = i+1; j < nAtoms; j++) {
        atom* pAt2 = atomList[j];
        if (pAt1->hasBondedAtom(pAt2) or
            pAt1->has13BondedAtom(pAt2) or
            pAt1->has14BondedAtom(pAt2)) {
          nE++;
        }
      }
      if (nE == 0) {
        nE = 1;
      }
      if (index < 10) {
        oprmtop << std::setw(8) << nE;
      }
      else {
        oprmtop << std::setw(8) << nE << std::endl;
        index = 0;
      }
      index++;
    }
    if (index != 1) {
      oprmtop << " " << std::endl;
      oprmtop.flush();
    }

    /*
      FORMAT(12I6)  (ICO(i), i=1,NTYPES*NTYPES)
      ICO    : provides the index to the nonbon parameter
               arrays CN1, CN2 and ASOL, BSOL.  All possible 6-12
               or 10-12 atoms type interactions are represented.
               NOTE: A particular atom type can have either a 10-12
               or a 6-12 interaction, but not both.  The index is
               calculated as follows:
                index = ICO(NTYPES*(IAC(i)-1)+IAC(j))
               If index is positive, this is an index into the
               6-12 parameter arrays (CN1 and CN2) otherwise it
               is an index into the 10-12 parameter arrays (ASOL
               and BSOL).
    */

    /*
      Lower triangular matrix:
                              j
          1       2       3       4       5       6       7
      |    --------------------------------------------------
  1   |    1       2       4       7      11      16      22
  2   |    2       3       5       8      12      17      23
  3   |    4       5       6       9      13      18      24
i 4   |    7       8       9      10      14      19      25
  5   |   11      12      13      14      15      20      26
  6   |   16      17      18      19      20      21      27
  7   |   22      23      24      25      26      27      28

          diagonal = [n*(n+1)]/2
    */

    oprmtop << "%FLAG NONBONDED_PARM_INDEX" << std::endl;
    oprmtop << "%FORMAT(10I8)" << std::endl;
    index = 1;

    int* nonBondedIndices;
    try {
      nonBondedIndices         = new int [nUniqueAtomTypes*nUniqueAtomTypes];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      //exit(0);
      throw MTKException(" Memory Allocation Failure ");
    }

    for (unsigned int p = 0; p < nUniqueAtomTypes; p++) {
      for (unsigned int pp = 0; pp < nUniqueAtomTypes; pp++) {
        if (p >= pp) {
          nonBondedIndices[p*nUniqueAtomTypes+pp] = index;
          nonBondedIndices[pp*nUniqueAtomTypes+p] = index;

          //std::cout <<  atomTypesUsed[p] << "-" << atomTypesUsed[pp] << " " << index << std::endl;
          index++;
        }
      }
    }

    index = 1;
    for (unsigned int p = 0; p < nUniqueAtomTypes; p++) {
      for (unsigned int pp = 0; pp < nUniqueAtomTypes; pp++) {
        if (index < 10) {
          if (p*nUniqueAtomTypes+pp != nUniqueAtomTypes*nUniqueAtomTypes-1) {
            oprmtop << std::setw(8) << nonBondedIndices[p*nUniqueAtomTypes+pp];
          }
          else {
            oprmtop << std::setw(8) << nonBondedIndices[p*nUniqueAtomTypes+pp] << std::endl;
          }
        }
        else {
          oprmtop << std::setw(8) << nonBondedIndices[p*nUniqueAtomTypes+pp] << std::endl;
          index = 0;
        }
        index++;
      }
    }
    oprmtop.flush();

    /*
      FORMAT(20A4)  (LABRES(i), i=1,NRES)
      LABRES : the residue labels
    */
    oprmtop << "%FLAG RESIDUE_LABEL" << std::endl;
    oprmtop << "%FORMAT(20a4)" << std::endl;
    index = 1;
    lStr = "";
    for (unsigned int s = 0; s < nResidues; s++) {
      lStr += smolList[s]->getName();
      lStr += " ";
      if (index < 20) {
        if (s == nResidues-1) {
          oprmtop << std::left << std::setw(80) << lStr << std::endl;
          oprmtop.unsetf(std::ios::left);
          lStr = "";
        }
      }
      else {
        oprmtop << std::setw(80) << lStr << std::endl;
        lStr = "";
        index = 0;
      }
      index++;
    }

    /*
      FORMAT(12I6)  (IPRES(i), i=1,NRES)
      IPRES  : atoms in each residue are listed for atom "i" in
               IPRES(i) to IPRES(i+1)-1
    */
    oprmtop << "%FLAG RESIDUE_POINTER" << std::endl;
    oprmtop << "%FORMAT(10I8)" << std::endl;
    index = 1;
    int nAts = 1;
    for (unsigned int s = 0; s < nResidues; s++) {
      if (index < 10) {
        if (s != nResidues-1) {
          oprmtop << std::setw(8) << nAts;
        }
        else {
          oprmtop << std::setw(8) << nAts << std::endl;
        }
      }
      else {
        oprmtop << std::setw(8) << nAts << std::endl;
        index = 0;
      }
      index++;
      nAts += smolList[s]->getNumAtoms();
    }
    oprmtop.flush();

    /*
      FORMAT(5E16.8)  (RK(i), i=1,NUMBND)
      RK     : force constant for the bonds of each type, kcal/mol
    */
    std::map<bondParam*, int> bondTypesMap;
    std::vector<bondParam*> bondTypesUsed = pCollection->getUniqueBondTypes();
    unsigned int nBondTypes = bondTypesUsed.size();
    for (unsigned int k = 0; k < bondTypesUsed.size(); k++) {
      bondTypesMap[bondTypesUsed[k]] = k;
    }

    oprmtop << "%FLAG BOND_FORCE_CONSTANT" << std::endl;
    oprmtop << "%FORMAT(5E16.8)" << std::endl;
    index = 1;
    for (unsigned int k = 0; k < nBondTypes; k++) {
      if (index < 5) {
        if (k != nBondTypes-1) {
          oprmtop << std::uppercase << std::scientific << std::setw(16)
                  << std::setprecision(8) << bondTypesUsed[k]->keq;
        }
        else {
          oprmtop << std::uppercase << std::scientific << std::setw(16)
                  << std::setprecision(8) << bondTypesUsed[k]->keq << std::endl;
        }
      }
      else {
        oprmtop << std::uppercase << std::scientific << std::setw(16)
                << std::setprecision(8) << bondTypesUsed[k]->keq << std::endl;
        index = 0;
      }
      index++;
    }
    oprmtop.unsetf(std::ios::scientific);
    oprmtop.flush();

    /*
      FORMAT(5E16.8)  (REQ(i), i=1,NUMBND)
        REQ    : the equilibrium bond length for the bonds of each type, angstroms
    */
    oprmtop << "%FLAG BOND_EQUIL_VALUE" << std::endl;
    oprmtop << "%FORMAT(5E16.8)" << std::endl;
    index = 1;
    for (unsigned int k = 0; k < nBondTypes; k++) {
      if (index < 5) {
        if (k != nBondTypes-1) {
          oprmtop << std::uppercase << std::scientific << std::setw(16)
                  << std::setprecision(8) << bondTypesUsed[k]->req;
        }
        else {
          oprmtop << std::uppercase << std::scientific << std::setw(16)
                  << std::setprecision(8) << bondTypesUsed[k]->req << std::endl;
        }
      }
      else {
        oprmtop << std::uppercase << std::scientific << std::setw(16)
                << std::setprecision(8) << bondTypesUsed[k]->req << std::endl;
        index = 0;
      }
      index++;
    }
    oprmtop.unsetf(std::ios::scientific);
    oprmtop.flush();

    /*
      FORMAT(5E16.8)  (TK(i), i=1,NUMANG)
        TK     : force constant for the angles of each type, kcal/mol A**2
    */

    std::map<angleParam*, int> angleTypesMap;
    std::vector<angleParam*> angleTypesUsed = pCollection->getUniqueAngleTypes();
    unsigned int nAngleTypes = angleTypesUsed.size();
    for (unsigned int k = 0; k < nAngleTypes; k++) {
      angleTypesMap[angleTypesUsed[k]] = k;
    }

    oprmtop << "%FLAG ANGLE_FORCE_CONSTANT" << std::endl;
    oprmtop << "%FORMAT(5E16.8)" << std::endl;
    index = 1;
    for (unsigned int k = 0; k < nAngleTypes; k++) {
      if (index < 5) {
        if (k != nAngleTypes-1) {
          oprmtop << std::uppercase << std::scientific << std::setw(16)
                  << std::setprecision(8) << angleTypesUsed[k]->keq;
        }
        else {
          oprmtop << std::uppercase << std::scientific << std::setw(16)
                  << std::setprecision(8) << angleTypesUsed[k]->keq << std::endl;
        }
      }
      else {
        oprmtop << std::uppercase << std::scientific << std::setw(16)
                << std::setprecision(8) << angleTypesUsed[k]->keq << std::endl;
        index = 0;
      }
      index++;
    }
    oprmtop.unsetf(std::ios::scientific);
    oprmtop.flush();

    /*
      FORMAT(5E16.8)  (TEQ(i), i=1,NUMANG)
        TEQ    : the equilibrium angle for the angles of each type, radians
    */
    oprmtop << "%FLAG ANGLE_EQUIL_VALUE" << std::endl;
    oprmtop << "%FORMAT(5E16.8)" << std::endl;
    index = 1;
    for (unsigned int k = 0; k < nAngleTypes; k++) {
      if (index < 5) {
        if (k != nAngleTypes-1) {
          oprmtop << std::uppercase << std::scientific << std::setw(16)
                  << std::setprecision(8) << angleTypesUsed[k]->req;
        }
        else {
          oprmtop << std::uppercase << std::scientific << std::setw(16)
                  << std::setprecision(8) << angleTypesUsed[k]->req << std::endl;
        }
      }
      else {
        oprmtop << std::uppercase << std::scientific << std::setw(16)
                << std::setprecision(8) << angleTypesUsed[k]->req << std::endl;
        index = 0;
      }
      index++;
    }
    oprmtop.unsetf(std::ios::scientific);
    oprmtop.flush();

    /*
      FORMAT(5E16.8)  (PK(i), i=1,NPTRA)
        PK     : force constant for the dihedrals of each type, kcal/mol
    */
    int nD = 0;
    std::map<torsionParam*, int> torsionTypesMap;
    std::vector<torsionParam*> torsionTypesUsed = pCollection->getUniqueTorsionTypes();
    unsigned int nTorsionTypes = torsionTypesUsed.size();
    for (unsigned int k = 0; k < nTorsionTypes; k++) {
      torsionTypesMap[torsionTypesUsed[k]] = nD;
      nD++;
    }

    std::map<improperParam*, int> improperTypesMap;
    std::vector<improperParam*> improperTypesUsed = pCollection->getUniqueImproperTypes();
    unsigned int nImproperTypes = improperTypesUsed.size();
    for (unsigned int k = 0; k < nImproperTypes; k++) {
      improperTypesMap[improperTypesUsed[k]] = nD;
      nD++;
    }

    /*
       Set Force Constants to Zero if you want no torsion energy
    */
    oprmtop << "%FLAG DIHEDRAL_FORCE_CONSTANT" << std::endl;
    oprmtop << "%FORMAT(5E16.8)" << std::endl;
    index = 1;
    for (unsigned int k = 0; k < nTorsionTypes; k++) {
      if (index < 5) {
        oprmtop << std::uppercase << std::scientific << std::setw(16)
                << std::setprecision(8) << torsionTypesUsed[k]->Vn / torsionTypesUsed[k]->npth;
      }
      else {
        oprmtop << std::uppercase << std::scientific << std::setw(16)
                << std::setprecision(8) << torsionTypesUsed[k]->Vn  / torsionTypesUsed[k]->npth << std::endl;
        index = 0;
      }
      index++;
    }

    for (unsigned int k = 0; k < nImproperTypes; k++) {
      if (index < 5) {
        if (k != nImproperTypes-1) {
          oprmtop << std::uppercase << std::scientific << std::setw(16)
                  << std::setprecision(8) << improperTypesUsed[k]->Vn;
        }
        else {
          oprmtop << std::uppercase << std::scientific << std::setw(16)
                  << std::setprecision(8) << improperTypesUsed[k]->Vn << std::endl;
        }
      }
      else {
        oprmtop << std::uppercase << std::scientific << std::setw(16)
                << std::setprecision(8) << improperTypesUsed[k]->Vn << std::endl;
        index = 0;
      }
      index++;
    }

    if (index != 1) oprmtop << " " << std::endl;
    oprmtop.unsetf(std::ios::scientific);
    oprmtop.flush();

    /*
      FORMAT(5E16.8)  (PN(i), i=1,NPTRA)
        PN     : periodicity of the dihedral of a given type
    */
    oprmtop << "%FLAG DIHEDRAL_PERIODICITY" << std::endl;
    oprmtop << "%FORMAT(5E16.8)" << std::endl;
    index = 1;
    for (unsigned int k = 0; k < nTorsionTypes; k++) {
      if (index < 5) {
        oprmtop << std::uppercase << std::scientific << std::setw(16) << std::setprecision(8)
                << std::abs(torsionTypesUsed[k]->Nt);
      }
      else {
        oprmtop << std::uppercase << std::scientific << std::setw(16) << std::setprecision(8)
                << std::abs(torsionTypesUsed[k]->Nt) << std::endl;
        index = 0;
      }
      index++;
    }

    for (unsigned int k = 0; k < nImproperTypes; k++) {
      if (index < 5) {
        if (k != nImproperTypes-1) {
          oprmtop << std::uppercase << std::scientific << std::setw(16) << std::setprecision(8)
                  << std::abs(improperTypesUsed[k]->Nt);
        }
        else {
          oprmtop << std::uppercase << std::scientific << std::setw(16) << std::setprecision(8)
                  << std::abs(improperTypesUsed[k]->Nt) << std::endl;
        }
      }
      else {
        oprmtop << std::uppercase << std::scientific << std::setw(16) << std::setprecision(8)
                << std::abs(improperTypesUsed[k]->Nt) << std::endl;
        index = 0;
      }
      index++;
    }
    if (index != 1) oprmtop << " " << std::endl;
    oprmtop.unsetf(std::ios::scientific);
    oprmtop.flush();

    /*
      FORMAT(5E16.8)  (PHASE(i), i=1,NPTRA)
        PHASE  : phase of the dihedral of a given type, radians
    */
    oprmtop << "%FLAG DIHEDRAL_PHASE" << std::endl;
    oprmtop << "%FORMAT(5E16.8)" << std::endl;
    index = 1;
    for (unsigned int k = 0; k < nTorsionTypes; k++) {
      if (index < 5) {
        oprmtop << std::uppercase << std::scientific << std::setw(16)
                << std::setprecision(8) << torsionTypesUsed[k]->gamma;
      }
      else {
        oprmtop << std::uppercase << std::scientific << std::setw(16)
                << std::setprecision(8) << torsionTypesUsed[k]->gamma << std::endl;
        index = 0;
      }
      index++;
    }

    for (unsigned int k = 0; k < nImproperTypes; k++) {
      if (index < 5) {
        if (k != nImproperTypes-1) {
          oprmtop << std::uppercase << std::scientific << std::setw(16)
                  << std::setprecision(8) << improperTypesUsed[k]->gamma;
        }
        else {
          oprmtop << std::uppercase << std::scientific << std::setw(16)
                  << std::setprecision(8) << improperTypesUsed[k]->gamma << std::endl;
        }
      }
      else {
        oprmtop << std::uppercase << std::scientific << std::setw(16) << std::setprecision(8)
                << improperTypesUsed[k]->gamma << std::endl;
        index = 0;
      }
      index++;
    }
    if (index != 1) oprmtop << " " << std::endl;

    oprmtop.unsetf(std::ios::scientific);
    oprmtop.flush();

    /*
      FORMAT(5E16.8)  (SOLTY(i), i=1,NATYP)
        SOLTY  : currently unused (reserved for future use)
    */
    oprmtop << "%FLAG SOLTY" << std::endl;
    oprmtop << "%FORMAT(5E16.8)" << std::endl;
    index = 1;
    for (unsigned int p = 0; p < nUniqueAtomTypes; p++) {
      if (index < 5) {
        if (p != nUniqueAtomTypes-1) {
          oprmtop << std::uppercase << std::scientific << std::setprecision(8)
                  << std::setw(16) << 0.0;
        }
        else {
          oprmtop << std::uppercase << std::scientific << std::setprecision(8)
                  << std::setw(16) << 0.0 << std::endl;
        }
      }
      else {
        oprmtop << std::uppercase << std::scientific << std::setprecision(8)
                << std::setw(16) << 0.0 << std::endl;
        index = 0;
      }
      index++;
    }
    oprmtop.unsetf(std::ios::scientific);
    oprmtop.flush();

    /*
      FORMAT(5E16.8)  (CN1(i), i=1,NTYPES*(NTYPES+1)/2)
        CN1    : Lennard Jones r**12 terms for all possible atom type
                 interactions, indexed by ICO and IAC; for atom i and j
                 where i < j, the index into this array is as follows
                 (assuming the value of ICO(index) is positive):
                 CN1(ICO(NTYPES*(IAC(i)-1)+IAC(j))).
      NOTE: Both CN1 and CN2 are UPPER Triangular matrix!
    */
    oprmtop << "%FLAG LENNARD_JONES_ACOEF" << std::endl;
    oprmtop << "%FORMAT(5E16.8)" << std::endl;
    index = 1;
    parameters* pParameters = pCollection->getParameters();
    pParameters->calculateSigmaEpsilon();

    int valuesSize = (nUniqueAtomTypes *(nUniqueAtomTypes+1))/2;
    double aCoefValues[valuesSize];
    double bCoefValues[valuesSize];

    // Set up indices array
    int indices[nUniqueAtomTypes*nUniqueAtomTypes];
    int l = 0;
    for (unsigned int i = 0; i < nUniqueAtomTypes; i++) {
      for (unsigned int j = 0; j < nUniqueAtomTypes; j++) {
        if (i >= j) {
          indices[j * nUniqueAtomTypes + i] = l;
          indices[i * nUniqueAtomTypes + j] = l;
          l++;
        }
      }
    }

    //std::cout << " i-j epsilon sigma acoeff bcoeff"<<std::endl;
    for (unsigned int p = 0; p < nUniqueAtomTypes; p++) {
      for (unsigned int pp = p; pp < nUniqueAtomTypes; pp++) {
        LJ612SE* lj = pParameters->getLJ612SE(atomTypesUsed[p], atomTypesUsed[pp]);
        if (!lj) {
          std::cout << " Error in amberParser ... exiting " << std::endl;
          //exit(0);
          throw MTKException(" Error in amberParser ");
        }
        double r12 = lj->epsilon * lj->sigma * lj->sigma;

        int paramIndex = indices[p * nUniqueAtomTypes + pp];
        aCoefValues[paramIndex] = r12;
        //std::cout << atomTypesUsed[p] << "-" << atomTypesUsed[pp] << " "
        //    << lj->epsilon << " " << lj->sigma << " " << r12 << " " << 2.0 * lj->epsilon * lj->sigma << std::endl;
      }
    }

    for (int p = 0; p < valuesSize; p++) {
      if (index < 5) {
        oprmtop << std::uppercase << std::scientific << std::setprecision(8)
                << std::setw(16) << aCoefValues[p];
      }
      else {
        oprmtop << std::uppercase << std::scientific << std::setprecision(8)
                << std::setw(16) << aCoefValues[p] << std::endl;
        index = 0;
      }
      index++;
    }
    if (index != 1) {
      oprmtop << " " << std::endl;
    }
    oprmtop.flush();

    /*
      FORMAT(5E16.8)  (CN2(i), i=1,NTYPES*(NTYPES+1)/2)
        CN2    : Lennard Jones r**6 terms for all possible atom type
                 interactions.  Indexed like CN1 above.
      NOTE: the atom numbers in the following arrays that describe bonds,
            angles, and dihedrals are coordinate array indexes for runtime
            speed. The true atom number equals the absolute value of the
            number divided by three, plus one. In the case of the dihedrals,
            if the fourth atom is negative, this implies that the dihedral
            is an improper. If the third atom is negative, this implies that
            the end group interations are to be ignored. End group interactions
            are ignored, for example, in dihedrals of various ring systems
            (to prevent double counting of 1-4 interactions) and in multiterm
            dihedrals.
      FORMAT(12I6)  (IBH(i),JBH(i),ICBH(i), i=1,NBONH)
        IBH    : atom involved in bond "i", bond contains hydrogen
        JBH    : atom involved in bond "i", bond contains hydrogen
        ICBH   : index into parameter arrays RK and REQ
    */
    oprmtop << "%FLAG LENNARD_JONES_BCOEF" << std::endl;
    oprmtop << "%FORMAT(5E16.8)" << std::endl;
    index = 1;
    for (unsigned int p = 0; p < nUniqueAtomTypes; p++) {
      for (unsigned int pp = p; pp < nUniqueAtomTypes; pp++) {
        LJ612SE* lj = pParameters->getLJ612SE(atomTypesUsed[p], atomTypesUsed[pp]);
        if (!lj) {
          std::cout << " Error in amberParser ... exiting " << std::endl;
          //exit(0);
          throw MTKException(" Error in amberParser ");
        }
        double r6 = 2.0 * lj->epsilon * lj->sigma;

        int paramIndex = indices[p * nUniqueAtomTypes + pp];
        bCoefValues[paramIndex] = r6;
      }
    }

    for (int p = 0; p < valuesSize; p++) {
      if (index < 5) {
        oprmtop << std::uppercase << std::scientific << std::setprecision(8)
                << std::setw(16) << bCoefValues[p];
      }
      else {
        oprmtop << std::uppercase << std::scientific << std::setprecision(8)
                << std::setw(16) << bCoefValues[p] << std::endl;
        index = 0;
      }
      index++;
    }

    if (index != 1) {
      oprmtop << " " << std::endl;
    }
    oprmtop.flush();

    /*
      FORMAT(12I6)  (IB(i),JB(i),ICB(i), i=1,NBONA)
        IB     : atom involved in bond "i", bond does contain hydrogen
        JB     : atom involved in bond "i", bond does contain hydrogen
        ICB    : index into parameter arrays RK and REQ
    */
    oprmtop << "%FLAG BONDS_INC_HYDROGEN" << std::endl;
    oprmtop << "%FORMAT(10I8)" << std::endl;
    index = 1;
    int bondInfo[3];
    for (unsigned int i = 0; i < bondList.size(); i++) {
      Bond* pBond = bondList[i];
      if ((pBond->atom1->getElementSymbol() == "H") or
          (pBond->atom2->getElementSymbol() == "H")) {
        bondParam* pBondParam = pBond->pBondParam;
        if (!pBondParam) {
          continue;
        }

        bondInfo[0] = (localAtomIndex[pBond->atom1]-1) * 3;
        bondInfo[1] = (localAtomIndex[pBond->atom2]-1) * 3;
        bondInfo[2] = bondTypesMap[pBondParam] + 1;
        for (unsigned int j = 0; j < 3; j++) {
          if (index < 10) {
            oprmtop << std::setw(8) << bondInfo[j];
          }
          else {
            oprmtop << std::setw(8) << bondInfo[j] << std::endl;
            index = 0;
          }
          index++;
        }
      }
    }
    if (index != 1) {
      oprmtop << " " << std::endl;
    }

    /*
      FORMAT(12I6)  (IB(i),JB(i),ICB(i), i=1,NBONA)
        IB     : atom involved in bond "i", bond does not contain hydrogen
        JB     : atom involved in bond "i", bond does not contain hydrogen
        ICB    : index into parameter arrays RK and REQ
    */
    oprmtop << "%FLAG BONDS_WITHOUT_HYDROGEN" << std::endl;
    oprmtop << "%FORMAT(10I8)" << std::endl;
    index = 1;
    for (unsigned int i = 0; i < bondList.size(); i++) {
      Bond* pBond = bondList[i];
      if ((pBond->atom1->getElementSymbol() != "H") and
          (pBond->atom2->getElementSymbol() != "H")) {
        bondParam* pBondParam = pBond->pBondParam;
        if (!pBondParam) {
          continue;
        }
        bondInfo[0] = (localAtomIndex[pBond->atom1]-1) * 3;
        bondInfo[1] = (localAtomIndex[pBond->atom2]-1) * 3;
        bondInfo[2] = bondTypesMap[pBondParam] + 1;

        for (unsigned int j = 0; j < 3; j++) {
          if (index < 10) {
            oprmtop << std::setw(8) << bondInfo[j];
          }
          else {
            oprmtop << std::setw(8) << bondInfo[j] << std::endl;
            index = 0;
          }
          index++;
        }
      }
    }
    if (index != 1) {
      oprmtop << " " << std::endl;
    }

    /*
      FORMAT(12I6)  (ITH(i),JTH(i),KTH(i),ICTH(i), i=1,NTHETH)
        ITH    : atom involved in angle "i", angle contains hydrogen
        JTH    : atom involved in angle "i", angle contains hydrogen
        KTH    : atom involved in angle "i", angle contains hydrogen
        ICTH   : index into parameter arrays TK and TEQ for angle
                 ITH(i)-JTH(i)-KTH(i)
    */
    oprmtop << "%FLAG ANGLES_INC_HYDROGEN" << std::endl;
    oprmtop << "%FORMAT(10I8)" << std::endl;
    index = 1;
    int angleInfo[4];
    for (unsigned int i = 0; i < angleList.size(); i++) {
      Angle* pAngle = angleList[i];
      if ((pAngle->atom1->getElementSymbol() == "H") or
          (pAngle->atom3->getElementSymbol() == "H")) {
        angleParam* pAngleParam = pAngle->pAngleParam;
        if (!pAngleParam) {
          continue;
        }

        angleInfo[0] = (localAtomIndex[pAngle->atom1]-1) * 3;
        angleInfo[1] = (localAtomIndex[pAngle->atom2]-1) * 3;
        angleInfo[2] = (localAtomIndex[pAngle->atom3]-1) * 3;

        angleInfo[3] = angleTypesMap[pAngleParam] + 1;
        for (unsigned int j = 0; j < 4; j++) {
          if (index < 10) {
            oprmtop << std::setw(8) << angleInfo[j];
          }
          else {
            oprmtop << std::setw(8) << angleInfo[j] << std::endl;
            index = 0;
          }
          index++;
        }
      }
    }
    if (index != 1) {
      oprmtop << " " << std::endl;
    }

    /*
      FORMAT(12I6)  (IT(i),JT(i),KT(i),ICT(i), i=1,NTHETA)
        IT     : atom involved in angle "i", angle does not contain hydrogen
        JT     : atom involved in angle "i", angle does not contain hydrogen
        KT     : atom involved in angle "i", angle does not contain hydrogen
        ICT    : index into parameter arrays TK and TEQ for angle
                 IT(i)-JT(i)-KT(i)
    */
    oprmtop << "%FLAG ANGLES_WITHOUT_HYDROGEN" << std::endl;
    oprmtop << "%FORMAT(10I8)" << std::endl;
    index = 1;
    for (unsigned int i = 0; i < angleList.size(); i++) {
      Angle* pAngle = angleList[i];
      if ((pAngle->atom1->getElementSymbol() != "H") and
          (pAngle->atom3->getElementSymbol() != "H")) {
        angleParam* pAngleParam = pAngle->pAngleParam;
        if (!pAngleParam) {
          continue;
        }

        angleInfo[0] = (localAtomIndex[pAngle->atom1]-1) * 3;
        angleInfo[1] = (localAtomIndex[pAngle->atom2]-1) * 3;
        angleInfo[2] = (localAtomIndex[pAngle->atom3]-1) * 3;

        angleInfo[3] = angleTypesMap[pAngleParam] + 1;
        for (unsigned int j = 0; j < 4; j++) {
          if (index < 10) {
            oprmtop << std::setw(8) << angleInfo[j];
          }
          else {
            oprmtop << std::setw(8) << angleInfo[j] << std::endl;
            index = 0;
          }
          index++;
        }
      }
    }
    if (index != 1) {
      oprmtop << " " << std::endl;
    }

    /*
      FORMAT(12I6)  (IPH(i),JPH(i),KPH(i),LPH(i),ICPH(i), i=1,NPHIH)
        IPH    : atom involved in dihedral "i", dihedral contains hydrogen
        JPH    : atom involved in dihedral "i", dihedral contains hydrogen
        KPH    : atom involved in dihedral "i", dihedral contains hydrogen
        LPH    : atom involved in dihedral "i", dihedral contains hydrogen
        ICPH   : index into parameter arrays PK, PN, and PHASE for
                 dihedral IPH(i)-JPH(i)-KPH(i)-LPH(i)
    */
    oprmtop << "%FLAG DIHEDRALS_INC_HYDROGEN" << std::endl;
    oprmtop << "%FORMAT(10I8)" << std::endl;
    index = 1;
    int dihedralInfo[5];
    torsionParam* pTorsionParam = 0;
    for (unsigned int i = 0; i < torsionList.size(); i++) {
      Torsion* pTorsion = torsionList[i];
      if ((pTorsion->atom1->getElementSymbol() == "H") or
          (pTorsion->atom4->getElementSymbol() == "H")) {

        std::vector<torsionParam*> torsionParamList =
        pParameters->getTorsionParamList(
        pTorsion->atom1->getStdAtom()->type, pTorsion->atom2->getStdAtom()->type,
        pTorsion->atom3->getStdAtom()->type, pTorsion->atom4->getStdAtom()->type);

        dihedralInfo[0] = (localAtomIndex[pTorsion->atom1]-1) * 3;
        dihedralInfo[1] = (localAtomIndex[pTorsion->atom2]-1) * 3;
        dihedralInfo[2] = (localAtomIndex[pTorsion->atom3]-1) * 3;
        dihedralInfo[3] = (localAtomIndex[pTorsion->atom4]-1) * 3;

        if (!torsionParamList.empty()) {
          int torCounter = 0;
          for (torsionParamIterator c = torsionParamList.begin();
               c != torsionParamList.end(); c++) {
            pTorsionParam = *c;

            if (torCounter > 0) dihedralInfo[2] = -dihedralInfo[2];
            dihedralInfo[4] = torsionTypesMap[pTorsionParam] + 1;
            for (unsigned int j = 0; j < 5; j++) {
              if (index < 10) {
                oprmtop << std::setw(8) << dihedralInfo[j];
              }
              else {
                oprmtop << std::setw(8) << dihedralInfo[j] << std::endl;
                index = 0;
              }
              index++;
            }
            torCounter++;
          }
        }
      }
    }

    improperParam* pImproperParam = 0;
    for (unsigned int i = 0; i < improperList.size(); i++) {
      Improper* pImproper = improperList[i];
      if ((pImproper->atom1->getElementSymbol() == "H") or
          (pImproper->atom2->getElementSymbol() == "H") or
          (pImproper->atom3->getElementSymbol() == "H") or
          (pImproper->atom4->getElementSymbol() == "H")) {

        std::vector<std::vector<int> > order;
        std::vector<improperParam*> improperParamList =
        pParameters->getImproperParamList(
        pImproper->atom1->getStdAtom()->type, pImproper->atom2->getStdAtom()->type,
        pImproper->atom3->getStdAtom()->type, pImproper->atom4->getStdAtom()->type, order);

        if (!improperParamList.empty()) {
          atom* tempImproper[4] = {pImproper->atom1, pImproper->atom2, pImproper->atom3, pImproper->atom4};
          int ii = 0;
          for (improperParamIterator c = improperParamList.begin();
               c != improperParamList.end(); c++) {
            pImproperParam = *c;

            dihedralInfo[0] =  (localAtomIndex[tempImproper[order[ii][0]]]-1) * 3;
            dihedralInfo[1] =  (localAtomIndex[tempImproper[order[ii][1]]]-1) * 3;
            dihedralInfo[2] = -(localAtomIndex[tempImproper[order[ii][2]]]-1) * 3;
            dihedralInfo[3] = -(localAtomIndex[tempImproper[order[ii][3]]]-1) * 3;

            ii++;
            dihedralInfo[4] = improperTypesMap[pImproperParam] + 1;
            for (unsigned int j = 0; j < 5; j++) {
              if (index < 10) {
                oprmtop << std::setw(8) << dihedralInfo[j];
              }
              else {
                oprmtop << std::setw(8) << dihedralInfo[j] << std::endl;
                index = 0;
              }
              index++;
            }
          }
        }
      }
    }
    if (index != 1) {
      oprmtop << " " << std::endl;
    }

    /*
      FORMAT(12I6)  (IP(i),JP(i),KP(i),LP(i),ICP(i), i=1,NPHIA)
        IP     : atom involved in dihedral "i", dihedral does not contain hydrogen
        JP     : atom involved in dihedral "i", dihedral does not contain hydrogen
        KP     : atom involved in dihedral "i", dihedral does not contain hydrogen
        LP     : atom involved in dihedral "i", dihedral does not contain hydrogen
        ICP    : index into parameter arrays PK, PN, and PHASE for
                 dihedral IPH(i)-JPH(i)-KPH(i)-LPH(i).  Note, if the
                 periodicity is negative, this implies the following entry
                 in the PK, PN, and PHASE arrays is another term in a
                 multitermed dihedral.
    */
    oprmtop << "%FLAG DIHEDRALS_WITHOUT_HYDROGEN" << std::endl;
    oprmtop << "%FORMAT(10I8)" << std::endl;
    index = 1;
    int numPrinted = 0;
    for (unsigned int i = 0; i < torsionList.size(); i++) {
      Torsion* pTorsion = torsionList[i];
      if ((pTorsion->atom1->getElementSymbol() != "H") and
          (pTorsion->atom4->getElementSymbol() != "H")) {

        std::vector<torsionParam*> torsionParamList =
        pParameters->getTorsionParamList(
        pTorsion->atom1->getStdAtom()->type, pTorsion->atom2->getStdAtom()->type,
        pTorsion->atom3->getStdAtom()->type, pTorsion->atom4->getStdAtom()->type);

        dihedralInfo[0] = (localAtomIndex[pTorsion->atom1]-1) * 3;
        dihedralInfo[1] = (localAtomIndex[pTorsion->atom2]-1) * 3;
        dihedralInfo[2] = (localAtomIndex[pTorsion->atom3]-1) * 3;
        dihedralInfo[3] = (localAtomIndex[pTorsion->atom4]-1) * 3;

        if (!torsionParamList.empty()) {
          if (pTorsion->atom1->getParent()->getName() == "PRO") {
            pParameters->removeProlineTorsion(torsionParamList);
          }
          int torCounter = 0;
          for (torsionParamIterator c = torsionParamList.begin();
               c != torsionParamList.end(); c++) {
            pTorsionParam = *c;
            if (torCounter > 0) dihedralInfo[2] = -(abs(dihedralInfo[2]));
            dihedralInfo[4] = torsionTypesMap[pTorsionParam] + 1;
            for (unsigned int j = 0; j < 5; j++) {
              if (index < 10) {
                oprmtop << std::setw(8) << dihedralInfo[j];
              }
              else {
                oprmtop << std::setw(8) << dihedralInfo[j] << std::endl;
                index = 0;
              }
              index++;
              numPrinted++;
            }
            torCounter++;
          }
        }
      }
    }

    for (unsigned int i = 0; i < improperList.size(); i++) {
      Improper* pImproper = improperList[i];
      if ((pImproper->atom1->getElementSymbol() != "H") and
          (pImproper->atom2->getElementSymbol() != "H") and
          (pImproper->atom3->getElementSymbol() != "H") and
          (pImproper->atom4->getElementSymbol() != "H")) {

        std::vector<std::vector<int> > order;
        std::vector<improperParam*> improperParamList =
        pParameters->getImproperParamList(
        pImproper->atom1->getStdAtom()->type, pImproper->atom2->getStdAtom()->type,
        pImproper->atom3->getStdAtom()->type, pImproper->atom4->getStdAtom()->type, order);

        if (!improperParamList.empty()) {
          atom* tempImproper[4] = {pImproper->atom1, pImproper->atom2, pImproper->atom3, pImproper->atom4};
          int ii = 0;
          for (improperParamIterator c = improperParamList.begin();
               c != improperParamList.end(); c++) {
            pImproperParam = *c;

            dihedralInfo[0] =  (localAtomIndex[tempImproper[order[ii][0]]]-1) * 3;
            dihedralInfo[1] =  (localAtomIndex[tempImproper[order[ii][1]]]-1) * 3;
            dihedralInfo[2] = -(localAtomIndex[tempImproper[order[ii][2]]]-1) * 3;
            dihedralInfo[3] = -(localAtomIndex[tempImproper[order[ii][3]]]-1) * 3;

            ii++;
            dihedralInfo[4] = improperTypesMap[pImproperParam] + 1;
            for (unsigned int j = 0; j < 5; j++) {
              if (index < 10) {
                oprmtop << std::setw(8) << dihedralInfo[j];
              }
              else {
                oprmtop << std::setw(8) << dihedralInfo[j] << std::endl;
                index = 0;
              }
              index++;
              numPrinted++;
            }
          }
        }
      }
    }
    if (index != 1) {
      oprmtop << " " << std::endl;
    }

    if (numPrinted == 0) {
      oprmtop << " " << std::endl;
    }

    /*
      FORMAT(12I6)  (NATEX(i), i=1,NEXT)
      NATEX  : the excluded atom list.  To get the excluded list for atom
               "i" you need to traverse the NUMEX list, adding up all
               the previous NUMEX values, since NUMEX(i) holds the number
               of excluded atoms for atom "i", not the index into the
               NATEX list.  Let IEXCL = SUM(NUMEX(j), j=1,i-1), then
               excluded atoms are NATEX(IEXCL) to NATEX(IEXCL+NUMEX(i)).
    */
    oprmtop << "%FLAG EXCLUDED_ATOMS_LIST" << std::endl;
    oprmtop << "%FORMAT(10I8)" << std::endl;
    index = 1;
    for (unsigned int i = 0; i < nAtoms; i++) {
      atom* pAt1 = atomList[i];
      bool g = 0;
      for (unsigned int j = i+1; j < nAtoms; j++) {
        atom* pAt2 = atomList[j];
        if (pAt1->hasBondedAtom(pAt2) or
            pAt1->has13BondedAtom(pAt2) or
            pAt1->has14BondedAtom(pAt2)) {
          g = 1;
          if (index < 10) {
            oprmtop << std::setw(8) << j + 1;
          }
          else {
            oprmtop << std::setw(8) << j + 1 << std::endl;
            index = 0;
          }
          index++;
        }
      }
      if (g == 0) {
        if (index < 10) {
          oprmtop << std::setw(8) << 0;
        }
        else {
          oprmtop << std::setw(8) << 0 << std::endl;
          index = 0;
        }
        index++;
      }
    }
    if (index != 1) {
      oprmtop << " " << std::endl;
    }

    /*
      FORMAT(5E16.8)  (ASOL(i), i=1,NPHB)
        ASOL   : the value for the r**12 term for hydrogen bonds of all
                 possible types.  Index into these arrays is equivalent
                 to the CN1 and CN2 arrays, however the index is negative.
                 For example, for atoms i and j, with i < j, the index is
                 -ICO(NTYPES*(IAC(i)-1+IAC(j)).
    */
    oprmtop << "%FLAG HBOND_ACOEF" << std::endl;
    oprmtop << "%FORMAT(5E16.8)" << std::endl;
    oprmtop << " " << std::endl;

    /*
      FORMAT(5E16.8)  (BSOL(i), i=1,NPHB)
        BSOL   : the value for the r**10 term for hydrogen bonds of all
                 possible types.  Indexed like ASOL.
    */
    oprmtop << "%FLAG HBOND_BCOEF" << std::endl;
    oprmtop << "%FORMAT(5E16.8)" << std::endl;
    oprmtop << " " << std::endl;

    /*
      FORMAT(5E16.8)  (HBCUT(i), i=1,NPHB)
        HBCUT  : no longer in use
    */
    oprmtop << "%FLAG HBCUT" << std::endl;
    oprmtop << "%FORMAT(5E16.8)" << std::endl;
    oprmtop << " " << std::endl;

    /*
      FORMAT(20A4)  (ISYMBL(i), i=1,NATOM)
        ISYMBL : the AMBER atom types for each atom
    */
    oprmtop << "%FLAG AMBER_ATOM_TYPE" << std::endl;
    oprmtop << "%FORMAT(20a4)" << std::endl;

    lStr = "";
    index = 1;
    for (unsigned int k = 0; k < nAtoms; k++) {
      stdAtom* pStdAtom = atomList[k]->getStdAtom();
      if (!pStdAtom) {
        //exit(0);
        throw MTKException(" amberParser::null pointer ");
      }

      lStr += pStdAtom->type;
      if (pStdAtom->type.size() == 1) {
        lStr += "   ";
      }
      else if (pStdAtom->type.size() == 2) {
        lStr += "  ";
      }
      if (index < 20) {
        if (k == nAtoms-1) {
          oprmtop << std::left << std::setw(80) << lStr << std::endl;
          oprmtop.unsetf(std::ios::left);
          lStr = "";
        }
      }
      else {
        oprmtop << std::setw(80) << lStr << std::endl;
        lStr = "";
        index = 0;
      }
      index++;
    }
    oprmtop.flush();

    /*
      FORMAT(20A4)  (ITREE(i), i=1,NATOM)
        ITREE  : the list of tree joining information, classified into five
                 types.  M -- main chain, S -- side chain, B -- branch point,
                 3 -- branch into three chains, E -- end of the chain
    */
    oprmtop << "%FLAG TREE_CHAIN_CLASSIFICATION" << std::endl;
    oprmtop << "%FORMAT(20a4)" << std::endl;
    lStr = "";
    index = 1;
    for (unsigned int k = 0; k < nAtoms; k++) {
      stdAtom* pStdAtom = atomList[k]->getStdAtom();
      if (!pStdAtom) {
        //exit(0);
        throw MTKException(" amberParser:: null pointer ");
      }

      lStr += pStdAtom->chain;
      lStr += "   ";
      if (index < 20) {
        if (k == nAtoms-1) {
          oprmtop << std::left << std::setw(80) << lStr << std::endl;
          oprmtop.unsetf(std::ios::left);
          lStr = "";
        }
      }
      else {
        oprmtop << std::setw(80) << lStr << std::endl;
        lStr = "";
        index = 0;
      }
      index++;
    }
    oprmtop.flush();

    /*
      FORMAT(12I6)  (JOIN(i), i=1,NATOM)
        JOIN   : tree joining information, potentially used in ancient
                 analysis programs.  Currently unused in sander or gibbs.
    */
    oprmtop << "%FLAG JOIN_ARRAY" << std::endl;
    oprmtop << "%FORMAT(10i8)" << std::endl;
    index = 1;
    for (unsigned int k = 0; k < nAtoms; k++) {
      if (index < 10) {
        oprmtop << std::setw(8) << 0;
      }
      else {
        oprmtop << std::setw(8) << 0 << std::endl;
        index = 0;
      }
      index++;
    }
    if (index != 1) {
      oprmtop << " " << std::endl;
      oprmtop.flush();
    }

    oprmtop << "%FLAG IROTAT" << std::endl;
    oprmtop << "%FORMAT(10i8)" << std::endl;
    index = 1;
    for (unsigned int k = 0; k < nAtoms; k++) {
      if (index < 10) {
        oprmtop << std::setw(8) << 0;
      }
      else {
        oprmtop << std::setw(8) << 0 << std::endl;
        index = 0;
      }
      index++;
    }
    if (index != 1) {
      oprmtop << " " << std::endl;
      oprmtop.flush();
    }

    oprmtop << "%FLAG RADII" << std::endl;
    oprmtop << "%FORMAT(5E16.8)" << std::endl;
    index = 1;
    for (unsigned int k = 0; k < nAtoms; k++) {
      if (index < 5) {
        oprmtop << std::uppercase << std::scientific << std::setw(16)
                << std::setprecision(8) << 0.0;
      }
      else {
        oprmtop << std::uppercase << std::scientific << std::setw(16)
                << std::setprecision(8) << 0.0 << std::endl;
        index = 0;
      }
      index++;
    }
    if (index != 1) {
      oprmtop << " " << std::endl;
      oprmtop.flush();
    }

    oprmtop << "%FLAG SCREEN" << std::endl;
    oprmtop << "%FORMAT(5E16.8)" << std::endl;
    index = 1;
    for (unsigned int k = 0; k < nAtoms; k++) {
      if (index < 5) {
        oprmtop << std::uppercase << std::scientific << std::setw(16)
                << std::setprecision(8) << 0.0;
      }
      else {
        oprmtop << std::uppercase << std::scientific << std::setw(16)
                << std::setprecision(8) << 0.0 << std::endl;
        index = 0;
      }
      index++;
    }
    if (index != 1) {
      oprmtop << " " << std::endl;
      oprmtop.flush();
    }

    oprmtop.close();
/*
    int seconds = 1000;
    clock_t endwait;
    endwait = clock () + seconds * CLOCKS_PER_SEC ;
    while (clock() < endwait) {}
*/
}

} // MTK++ namespace


