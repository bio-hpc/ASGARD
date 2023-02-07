/*!
   \file mmPotential.cpp
   \brief Base class for MM potentials
   \author Martin Peters

   $Date: 2010/03/29 20:28:34 $
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

#include "mmPotential.h"
#include <stdlib.h>
#include <math.h>
#include <cmath>

#include "Diagnostics/MTKException.h"
#include "Log/errorHandler.h"

namespace MTKpp
{

// ============================================================
// Function : mmPotential()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
mmPotential::mmPotential()
{
    this->nAtoms             = 0;
    this->nDistinctTypes     = 0;
    this->nBondsWithH        = 0;
    this->nBondsWithOutH     = 0;
    this->nAnglesWithH       = 0;
    this->nAnglesWithOutH    = 0;
    this->nDihedralsWithH    = 0;
    this->nDihedralsWithOutH = 0;
    this->nResidues          = 0;
    this->nUniqueBonds       = 0;
    this->nUniqueAngles      = 0;
    this->nUniqueDihedrals   = 0;
    this->nAtomsBigRes       = 0;

    this->xyz                = 0;
    this->gradients          = 0;

    this->nBonds             = 0;
    this->bonds              = 0;
    this->bondParams         = 0;

    this->nAngles            = 0;
    this->angles             = 0;
    this->angleParams        = 0;

    this->nTorsions          = 0;
    this->torsions           = 0;
    this->torsionParams      = 0;

    this->nImpropers         = 0;
    this->impropers          = 0;
    this->improperParams     = 0;

    this->nNonBonded         = 0;
    this->nonBonded          = 0;
    this->nonBondedParams    = 0;

    this->nonBondedParameterIndex = 0;
    this->r6Params           = 0;
    this->r12Params          = 0;
    this->nExcluded          = 0;
    this->numExcluded        = 0;
    this->excluded           = 0;
    this->nExcluded14        = 0;
    this->numExcluded14      = 0;
    this->excluded14         = 0;
    this->atomFlags          = 0;
    this->outFlags           = 0;
    this->iTypes             = 0;
    this->charges            = 0;
    this->masses             = 0;
    this->resPointers        = 0;
    this->initialize();
}

// ============================================================
// Function : ~mmPotential()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
mmPotential::~mmPotential()
{
    delete [] xyz;
    delete [] bonds;
    delete [] bondParams;
    delete [] angles;
    delete [] angleParams;
    delete [] torsions;
    delete [] torsionParams;
    delete [] impropers;
    delete [] improperParams;
    delete [] nonBonded;
    delete [] nonBondedPtrs;
    delete [] nonBonded14Ptrs;
    delete [] nonBondedParameterIndex;
    delete [] charges;
    delete [] symbols;
    delete [] names;
    delete [] masses;
    delete [] iTypes;
    delete [] cTypes;
    delete [] r6Params;
    delete [] r12Params;
    delete [] numExcluded;
    delete [] numExcluded14;
    delete [] atomFlags;
    delete [] outFlags;
    delete [] excluded;
    delete [] excluded14;
    delete [] resNames;
    delete [] resPointers;
    if (this->bGradient) {
      delete [] gradients;
    }
}

// ============================================================
// Function : initialize()
// ------------------------------------------------------------
//
// ============================================================
void mmPotential::initialize()
{
    this->bEnergy    = 1;
    this->bGradient  = 0;

    this->bBond      = 0;
    this->bAngle     = 0;
    this->bTorsion   = 0;
    this->bImproper  = 0;
    this->bNonBonded = 0;
    this->bHBond     = 0;

    this->bBondDecompose = 0;
    this->bAngleDecompose = 0;
    this->bTorsionDecompose = 0;
    this->bImproperDecompose = 0;
    this->bNonBondedDecompose = 0;
    this->bHBondDecompose = 0;
    this->pwdFile   =  "";
    this->pwdGZFile = "";

    bondEnergy      = 0.0;
    angleEnergy     = 0.0;
    torsionEnergy   = 0.0;
    improperEnergy  = 0.0;

    vdWEnergy       = 0.0;
    vdW14Energy     = 0.0;

    R6              = 0.0;
    R12             = 0.0;

    eleEnergy       = 0.0;
    ele14Energy     = 0.0;

    nonBondedEnergy = 0.0;

    hBondEnergy     = 0.0;
    totalEnergy     = 0.0;
}

// ============================================================
// Function : setNumAtoms()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::setNumAtoms(int n)
{
    this->nAtoms = n;
    try {
      this->xyz         = new double [this->nAtoms*3]; // ATOM_COORD
      this->charges     = new double [this->nAtoms  ]; // CHARGE
      this->symbols     = new char   [this->nAtoms*2]; // ATOM_ELEMENT
      this->names       = new char   [this->nAtoms*4]; // ATOM_NAME
      this->masses      = new double [this->nAtoms  ]; // MASS
      this->iTypes      = new int    [this->nAtoms  ]; // ATOM_TYPE_INDEX
      this->cTypes      = new char   [this->nAtoms*2]; // ATOM_TYPE_CHAR
      this->numExcluded = new int    [this->nAtoms  ]; // NUMBER_EXCLUDED_ATOMS
      this->numExcluded14 = new int  [this->nAtoms  ]; // NUMBER_14EXCLUDED_ATOMS
      this->atomFlags   = new int    [this->nAtoms  ]; //
      this->outFlags    = new int    [this->nAtoms  ]; // OUTPUT_FLAGS_DECOMPOSE
      if (this->bGradient) {
        this->gradients = new double [this->nAtoms*3]; //
      }
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function : setNumTypes()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::setNumTypes(int n)
{
    this->nDistinctTypes = n;

    try {
      this->nonBondedParameterIndex = new int [this->nDistinctTypes * this->nDistinctTypes]; // NONBONDED_PARM_INDEX
      this->r12Params = new double [(this->nDistinctTypes * (this->nDistinctTypes+1))/2];    // LENNARD_JONES_ACOEF
      this->r6Params  = new double [(this->nDistinctTypes * (this->nDistinctTypes+1))/2];    // LENNARD_JONES_BCOEF
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }

    int l = 0;
    for (int i = 0; i < this->nDistinctTypes; i++) {
      for (int j = 0; j < this->nDistinctTypes; j++) {
        if (i >= j) {
          this->nonBondedParameterIndex[j * this->nDistinctTypes + i] = l;
          this->nonBondedParameterIndex[i * this->nDistinctTypes + j] = l;
          l++;
        }
      }
    }
    return 0;
}

// ============================================================
// Function : setNumResidues()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::setNumResidues(int n)
{
    this->nResidues = n; // NRES
    try {
      this->resNames    = new char   [this->nResidues*3]; // RESIDUE_LABEL
      this->resPointers = new int    [this->nResidues];   // RESIDUE_POINTER
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function : setExcludedSize()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::setExcludedSize(int n, int n2)
{
    this->nExcluded = n; // NEXT
    this->nExcluded14 = n2;
    try {
      this->excluded = new int    [this->nExcluded];   // EXCLUDED_ATOMS_LIST
      this->excluded14 = new int    [this->nExcluded14];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getNumAtoms()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::getNumAtoms()
{
    return this->nAtoms;
}

// ============================================================
// Function : getResidueNames()
// ------------------------------------------------------------
//
// ============================================================
char* mmPotential::getResidueNames()
{
    return this->resNames;
}

// ============================================================
// Function : getResiduePointers()
// ------------------------------------------------------------
//
// ============================================================
int* mmPotential::getResiduePointers()
{
    return this->resPointers;
}

// ============================================================
// Function : setNumAtomsBigRes()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::setNumAtomsBigRes(int n)
{
    this->nAtomsBigRes = n;
    return 0;
}

// ============================================================
// Function : getCoords()
// ------------------------------------------------------------
//
// ============================================================
double* mmPotential::getCoords()
{
    return this->xyz;
}

// ============================================================
// Function : getCharges()
// ------------------------------------------------------------
//
// ============================================================
double* mmPotential::getCharges()
{
    return this->charges;
}

// ============================================================
// Function : getSymbols()
// ------------------------------------------------------------
//
// ============================================================
char* mmPotential::getSymbols()
{
    return this->symbols;
}

// ============================================================
// Function : getAtomNames()
// ------------------------------------------------------------
//
// ============================================================
char* mmPotential::getAtomNames()
{
    return this->names;
}

// ============================================================
// Function : getNumTypes()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::getNumTypes()
{
    return this->nDistinctTypes;
}

// ============================================================
// Function : getIntTypes()
// ------------------------------------------------------------
//
// ============================================================
int* mmPotential::getIntTypes()
{
    return this->iTypes;
}

// ============================================================
// Function : getCharTypes()
// ------------------------------------------------------------
//
// ============================================================
char* mmPotential::getCharTypes()
{
    return this->cTypes;
}

// ============================================================
// Function : getNumExcluded()
// ------------------------------------------------------------
//
// ============================================================
int* mmPotential::getNumExcluded()
{
    return this->numExcluded;
}

// ============================================================
// Function : getExcluded()
// ------------------------------------------------------------
//
// ============================================================
int* mmPotential::getExcluded()
{
    return this->excluded;
}

// ============================================================
// Function : getNumExcluded14()
// ------------------------------------------------------------
//
// ============================================================
int* mmPotential::getNumExcluded14()
{
    return this->numExcluded14;
}

// ============================================================
// Function : getExcluded14()
// ------------------------------------------------------------
//
// ============================================================
int* mmPotential::getExcluded14()
{
    return this->excluded14;
}

// ============================================================
// Function : getAtomFlags()
// ------------------------------------------------------------
//
// ============================================================
int* mmPotential::getAtomFlags()
{
    return this->atomFlags;
}

// ============================================================
// Function : getGradients()
// ------------------------------------------------------------
//
// ============================================================
double* mmPotential::getGradients()
{
    return this->gradients;
}

// ============================================================
// Function : setNumBonds()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::setNumBonds(int n)
{
    this->nBonds = n;
    try {
      this->bonds       = new int [this->nBonds*2];
      this->bondParams  = new double [this->nBonds*2];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }
    for (int i = 0; i > this->nBonds*2; i++) {
      this->bonds[i] = -1;
      this->bondParams[i] = 0.0;
    }
    return 0;
}

// ============================================================
// Function : getNumBonds()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::getNumBonds()
{
    return this->nBonds;
}

// ============================================================
// Function : setNumBondsWithH()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::setNumBondsWithH(int n)
{
    this->nBondsWithH = n;
    return 0;
}

// ============================================================
// Function : setNumBondsWithOutH()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::setNumBondsWithOutH(int n)
{
    this->nBondsWithOutH = n;
    return 0;
}

// ============================================================
// Function : setNumUniqueBonds()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::setNumUniqueBonds(int n)
{
    this->nUniqueBonds = n;
    return 0;
}

// ============================================================
// Function : getBonds()
// ------------------------------------------------------------
//
// ============================================================
int* mmPotential::getBonds()
{
    return this->bonds;
}

// ============================================================
// Function : getBondParams()
// ------------------------------------------------------------
//
// ============================================================
double* mmPotential::getBondParams()
{
    return this->bondParams;
}

// ============================================================
// Function : setNumAngles()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::setNumAngles(int n)
{
    this->nAngles = n;
    try {
      this->angles      = new int [nAngles*3];
      this->angleParams = new double [nAngles*2];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }
    for (int i = 0; i > this->nAngles*3; i++) {
      this->angles[i] = -1;
    }
    for (int i = 0; i > this->nAngles*2; i++) {
      this->angleParams[i] = 0.0;
    }
    return 0;
}

// ============================================================
// Function : getNumAngles()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::getNumAngles()
{
    return this->nAngles;
}

// ============================================================
// Function : setNumAnglesWithH()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::setNumAnglesWithH(int n)
{
    this->nAnglesWithH = n;
    return 0;
}

// ============================================================
// Function : setNumAnglesWithOutH()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::setNumAnglesWithOutH(int n)
{
    this->nAnglesWithOutH = n;
    return 0;
}

// ============================================================
// Function : setNumUniqueAngles()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::setNumUniqueAngles(int n)
{
    this->nUniqueAngles = n;
    return 0;
}

// ============================================================
// Function : getAngles()
// ------------------------------------------------------------
//
// ============================================================
int* mmPotential::getAngles()
{
    return this->angles;
}

// ============================================================
// Function : getAngleParams()
// ------------------------------------------------------------
//
// ============================================================
double* mmPotential::getAngleParams()
{
    return this->angleParams;
}

// ============================================================
// Function : setNumTorsions()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::setNumTorsions(int n)
{
    this->nTorsions = n;
    try {
      this->torsions    = new int [this->nTorsions*4];
      this->torsionParams = new double [this->nTorsions*3];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }
    for (int i = 0; i > this->nTorsions*4; i++) {
      this->torsions[i] = -1;
    }
    for (int i = 0; i > this->nTorsions*3; i++) {
      this->torsionParams[i] = 0.0;
    }
    return 0;
}

// ============================================================
// Function : getNumTorsions()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::getNumTorsions()
{
    return this->nTorsions;
}

// ============================================================
// Function : getTorsions()
// ------------------------------------------------------------
//
// ============================================================
int* mmPotential::getTorsions()
{
    return this->torsions;
}

// ============================================================
// Function : getTorsionParams()
// ------------------------------------------------------------
//
// ============================================================
double* mmPotential::getTorsionParams()
{
    return this->torsionParams;
}

/*
// ============================================================
// Function : getNTorsionParams()
// ------------------------------------------------------------
//
// ============================================================
int* mmPotential::getNTorsionParams()
{
    return this->nTorsionParams;
}

// ============================================================
// Function : setNumTorsionParams()
// ------------------------------------------------------------
//
// ============================================================

void mmPotential::setNumTorsionParams(int n)
{
    this->totalNumTorsionParams = n;
    try {
      this->torsionParams = new double [n*4];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      //return 1; // make it return an int
    }
}
*/

// ============================================================
// Function : setNumImpropers()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::setNumImpropers(int n)
{
    this->nImpropers = n;
    try {
      this->impropers   = new int [nImpropers*4];
      this->improperParams = new double [nImpropers*3];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getNumImpropers()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::getNumImpropers()
{
    return this->nImpropers;
}

// ============================================================
// Function : getImpropers()
// ------------------------------------------------------------
//
// ============================================================
int* mmPotential::getImpropers()
{
    return this->impropers;
}

// ============================================================
// Function : getImproperParams()
// ------------------------------------------------------------
//
// ============================================================
double* mmPotential::getImproperParams()
{
    return this->improperParams;
}

// ============================================================
// Function : setNumDihedralsWithH()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::setNumDihedralsWithH(int n)
{
    this->nDihedralsWithH = n;
    return 0;
}

// ============================================================
// Function : setNumDihedralsWithOutH()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::setNumDihedralsWithOutH(int n)
{
    this->nDihedralsWithOutH = n;
    return 0;
}

// ============================================================
// Function : setNumUniqueDihedrals()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::setNumUniqueDihedrals(int n)
{
    this->nUniqueDihedrals = n;
    return 0;
}

// ============================================================
// Function : setNumNonBonded()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::setNumNonBonded(int n)
{
    this->nNonBonded = n;
    try {
      this->nonBonded   = new int [this->nNonBonded];
      this->nonBondedPtrs = new int [this->nAtoms];
      this->nonBondedParams = new double [this->nNonBonded*2];
      this->nonBonded14Ptrs = new int [this->nNonBonded];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }
    return 0;
}

// ============================================================
// Function : getNumNonBonded()
// ------------------------------------------------------------
//
// ============================================================
int mmPotential::getNumNonBonded()
{
    return this->nNonBonded;
}

// ============================================================
// Function : getNonBonded()
// ------------------------------------------------------------
//
// ============================================================
int* mmPotential::getNonBonded()
{
    return this->nonBonded;
}

// ============================================================
// Function : getNonBondedParams()
// ------------------------------------------------------------
//
// ============================================================
double* mmPotential::getNonBondedParams()
{
    return this->nonBondedParams;
}

// ============================================================
// Function : getNonBondedPtrs()
// ------------------------------------------------------------
//
// ============================================================
int* mmPotential::getNonBondedPtrs()
{
    return this->nonBondedPtrs;
}

// ============================================================
// Function : getNonBonded14Ptrs()
// ------------------------------------------------------------
//
// ============================================================
int* mmPotential::getNonBonded14Ptrs()
{
    return this->nonBonded14Ptrs;
}

// ============================================================
// Function : getNonBondedParameterIndex()
// ------------------------------------------------------------
//
// ============================================================
int* mmPotential::getNonBondedParameterIndex()
{
    return this->nonBondedParameterIndex;
}

// ============================================================
// Function : getR6Params()
// ------------------------------------------------------------
//
// ============================================================
double* mmPotential::getR6Params()
{
    return this->r6Params;
}

// ============================================================
// Function : getR12Params()
// ------------------------------------------------------------
//
// ============================================================
double* mmPotential::getR12Params()
{
    return this->r12Params;
}

// ============================================================
// Function : calcForces()
// ------------------------------------------------------------
//
// ============================================================
void mmPotential::calcForces(int i)
{
    if (i) {
      this->bGradient = 1;
    }
    else {
      this->bGradient = 0;
    }
}

// ============================================================
// Function : setPotential()
// ------------------------------------------------------------
//
// ============================================================
void mmPotential::setPotential(bool bBond, bool bAngle, bool bTorsion,
                               bool bImproper, bool bNonBonded, bool bHBond)
{
    this->bBond      = bBond;
    this->bAngle     = bAngle;
    this->bTorsion   = bTorsion;
    this->bImproper  = bImproper;
    this->bNonBonded = bNonBonded;
    this->bHBond     = bHBond;
}

// ============================================================
// Function : getPotential()
// ------------------------------------------------------------
//
// ============================================================
void mmPotential::getPotential(bool& bBond, bool& bAngle, bool& bTorsion,
                               bool& bImproper, bool& bNonBonded, bool& bHBond)
{
    bBond = this->bBond;
    bAngle = this->bAngle;
    bTorsion = this->bTorsion;
    bImproper = this->bImproper;
    bNonBonded = this->bNonBonded;
    bHBond = this->bHBond;
}

// ============================================================
// Function : calcEnergy()
// ------------------------------------------------------------
//
// ============================================================
void mmPotential::calcEnergy()
{
#ifdef HAVE_ZLIB
    if (this->bBondDecompose or this->bAngleDecompose or
        this->bTorsionDecompose or this->bImproperDecompose or
        this->bNonBondedDecompose or this->bHBondDecompose) {
      this->pwdGZFileStream = gzopen(pwdGZFile.c_str(), "wb9");
      if (this->pwdGZFileStream == 0) {
        std::string errMessage = "Unable to open MM pwd.gz file: " + pwdGZFile;
        errorLogger.throwError(__FUNCTION__, errMessage, 1);
        throw MTKException(errMessage);
      }
    }

    this->totalEnergy = 0;

    if (this->bBond) {
      try {
        this->totalEnergy += calcBondEnergy();
      }
      catch (std::exception& e) {
        std::cout << "Standard exception: " << e.what() << std::endl;
      }
    }

    if (this->bAngle) {
      this->totalEnergy += calcAngleEnergy();
    }

    if (this->bTorsion) {
      this->totalEnergy += calcTorsionEnergy();
    }

    if (this->bImproper) {
      this->totalEnergy += calcImproperEnergy();
    }

    if (this->bNonBonded) {
      this->totalEnergy += calcNonBondedEnergy();
    }

    if (this->bHBond) {
      this->totalEnergy += calcHBondEnergy();
    }

    if (this->pwdFileStream.is_open()) {
      this->pwdFileStream.close();
    }

    if (this->bBondDecompose or this->bAngleDecompose or
        this->bTorsionDecompose or this->bImproperDecompose or
        this->bNonBondedDecompose or this->bHBondDecompose) {
      gzclose(this->pwdGZFileStream);
    }
#endif
}

// ============================================================
// Function : decompose()
// ------------------------------------------------------------
//
// ============================================================
void mmPotential::decomposeEnergy(bool bBond, bool bAngle, bool bTorsion,
                  bool bImproper, bool bNonBonded, bool bHBond,
                  std::string f)
{
    this->bBondDecompose = bBond;
    this->bAngleDecompose = bAngle;
    this->bTorsionDecompose = bTorsion;
    this->bImproperDecompose = bImproper;
    this->bNonBondedDecompose = bNonBonded;
    this->bHBondDecompose = bHBond;
    this->pwdGZFile = f;
}

// ============================================================
// Function : getOutputFlags()
// ------------------------------------------------------------
//
// ============================================================
int* mmPotential::getOutputFlags()
{
    return this->outFlags;
}

// ============================================================
// Function : printEnergy()
// ------------------------------------------------------------
// Print Individual terms after a potential energy calculation
// ============================================================
void mmPotential::printEnergy()
{
#ifdef DEBUG
        std::cout << "  mmPotential::printEnergy " << std::endl;
#endif

    std::cout << "  (1)  BOND ENERGY               = " << this->bondEnergy << std::endl;
    std::cout << "  (2)  ANGLE ENERGY              = " << this->angleEnergy << std::endl;
    std::cout << "  (3)  TORSION ENERGY            = " << this->torsionEnergy << std::endl;
    std::cout << "  (4)  IMPROPER ENERGY           = " << this->improperEnergy << std::endl;
    std::cout << "  (5)  DIHEDRAL ENERGY (3+4)     = " << this->torsionEnergy+this->improperEnergy << std::endl;
    std::cout << "  (6)  1-4 VDW ENERGY            = " << this->vdW14Energy << std::endl;
    std::cout << "  (7)  VDW ENERGY                = " << this->vdWEnergy << std::endl;
    std::cout << "  (8)  TOTAL VDW ENERGY (6+7)    = " << this->vdW14Energy + this->vdWEnergy << std::endl;
    std::cout << "  (9)  ELE ENERGY                = " << this->eleEnergy << std::endl;
    std::cout << " (10)  1-4 ELE ENERGY            = " << this->ele14Energy << std::endl;
    std::cout << " (11)  TOTAL ELE ENERGY (9+10)   = " << this->ele14Energy + this->eleEnergy << std::endl;
    std::cout << " (12)  NON-BONDED ENERGY (8+11)  = " << this->nonBondedEnergy << std::endl;
    std::cout << " (13)  H-BOND ENERGY             = " << this->hBondEnergy << std::endl;
    std::cout << " --------------------------------------------------" << std::endl;
    std::cout << "       TOTAL ENERGY (1+2+5+12)   = " << this->totalEnergy << std::endl;
    std::cout << "       R^6 LENNARD JONES ENERGY  = " << this->R6 << std::endl;
    std::cout << "       R^12 LENNARD JONES ENERGY = " << this->R12 << std::endl;
}

// ============================================================
// Function : printEnergy()
// ------------------------------------------------------------
// Print Individual terms after a potential energy calculation
// ============================================================
void mmPotential::printEnergy(std::ostream& os)
{
    os << "  (1)  BOND ENERGY               = " << this->bondEnergy << std::endl;
    os << "  (2)  ANGLE ENERGY              = " << this->angleEnergy << std::endl;
    os << "  (3)  TORSION ENERGY            = " << this->torsionEnergy << std::endl;
    os << "  (4)  IMPROPER ENERGY           = " << this->improperEnergy << std::endl;
    os << "  (5)  DIHEDRAL ENERGY (3+4)     = " << this->torsionEnergy+this->improperEnergy << std::endl;
    os << "  (6)  1-4 VDW ENERGY            = " << this->vdW14Energy << std::endl;
    os << "  (7)  VDW ENERGY                = " << this->vdWEnergy << std::endl;
    os << "  (8)  TOTAL VDW ENERGY (6+7)    = " << this->vdW14Energy + this->vdWEnergy << std::endl;
    os << "  (9)  ELE ENERGY                = " << this->eleEnergy << std::endl;
    os << " (10)  1-4 ELE ENERGY            = " << this->ele14Energy << std::endl;
    os << " (11)  TOTAL ELE ENERGY (9+10)   = " << this->ele14Energy + this->eleEnergy << std::endl;
    os << " (12)  NON-BONDED ENERGY (8+11)  = " << this->nonBondedEnergy << std::endl;
    os << " (13)  H-BOND ENERGY             = " << this->hBondEnergy << std::endl;
    os << " --------------------------------------------------" << std::endl;
    os << "       TOTAL ENERGY (1+2+5+12)   = " << this->totalEnergy << std::endl;
    os << "       R^6 LENNARD JONES ENERGY  = " << this->R6 << std::endl;
    os << "       R^12 LENNARD JONES ENERGY = " << this->R12 << std::endl;
}

// ============================================================
// Function : printEnergy2()
// ------------------------------------------------------------
// Print Individual terms after a potential energy calculation
// ============================================================
void mmPotential::printEnergy2()
{
    std::cout << this->bondEnergy << " " << this->angleEnergy << " "
              << this->torsionEnergy + this->improperEnergy << " "
              << this->vdWEnergy << " " << this->eleEnergy << " "
              << this->vdW14Energy << " " << this->ele14Energy << " "
              << this->totalEnergy << std::endl;
}

// ============================================================
// Function : printEnergy2()
// ------------------------------------------------------------
// Print Individual terms after a potential energy calculation
// ============================================================
void mmPotential::printEnergy2(std::ostream& os)
{
    os << this->bondEnergy << " " << this->angleEnergy << " "
       << this->torsionEnergy + this->improperEnergy << " "
       << this->vdWEnergy << " " << this->eleEnergy << " "
       << this->vdW14Energy << " " << this->ele14Energy << " "
       << this->totalEnergy << std::endl;
}

// ============================================================
// Function : printEnergy2()
// ------------------------------------------------------------
// Print Individual terms after a potential energy calculation
// ============================================================
void mmPotential::printEnergy2(std::string t)
{
    std::cout << t << " " << this->bondEnergy << " " << this->angleEnergy << " "
              << this->torsionEnergy + this->improperEnergy << " "
              << this->vdWEnergy << " " << this->eleEnergy << " "
              << this->vdW14Energy << " " << this->ele14Energy << " "
              << this->totalEnergy << std::endl;
}

// ============================================================
// Function : printEnergy2()
// ------------------------------------------------------------
// Print Individual terms after a potential energy calculation
// ============================================================
void mmPotential::printEnergy2(std::string t, std::string t2)
{
    std::cout << t << " " << this->bondEnergy << " " << this->angleEnergy << " "
              << this->torsionEnergy + this->improperEnergy << " "
              << this->vdWEnergy << " " << this->eleEnergy << " "
              << this->vdW14Energy << " " << this->ele14Energy << " "
              << this->totalEnergy << " " << t2 << std::endl;
}

// ============================================================
// Function : printEnergy2()
// ------------------------------------------------------------
// Print Individual terms after a potential energy calculation
// ============================================================
void mmPotential::printEnergy2(std::ostream& os, std::string t, std::string t2)
{
    os << t << " " << this->bondEnergy << " " << this->angleEnergy << " "
              << this->torsionEnergy + this->improperEnergy << " "
              << this->vdWEnergy << " " << this->eleEnergy << " "
              << this->vdW14Energy << " " << this->ele14Energy << " "
              << this->totalEnergy << " " << t2 << std::endl;
}

// ============================================================
// Function : printCoords()
// ------------------------------------------------------------
// Print Individual terms after a potential energy calculation
// ============================================================
void mmPotential::printCoords()
{
    for (int i = 0; i < nAtoms; i++) {
      for (int j = 0; j < 3; j++) {
        std::cout << " " << this->xyz[i*nAtoms+j];
      }
      std::cout << " " << std::endl;
    }
}

// ============================================================
// Function : calcBondEnergy()
// ------------------------------------------------------------
//
// ============================================================
double mmPotential::calcBondEnergy()
{
#ifdef DEBUG
        std::cout << "  mmPotential::calcBondEnergy " << std::endl;
#endif
    this->bondEnergy = 0.0;
    return 0.0;
}

// ============================================================
// Function : calcAngleEnergy()
// ------------------------------------------------------------
//
// ============================================================
double mmPotential::calcAngleEnergy()
{
#ifdef DEBUG
        std::cout << "  mmPotential::calcAngleEnergy " << std::endl;
#endif
    this->angleEnergy = 0.0;
    return 0.0;
}

// ============================================================
// Function : calcTorsionEnergy()
// ------------------------------------------------------------
//
// ============================================================
double mmPotential::calcTorsionEnergy()
{
#ifdef DEBUG
        std::cout << "  mmPotential::calcTorsionEnergy " << std::endl;
#endif
    this->torsionEnergy = 0.0;
    return 0.0;
}

// ============================================================
// Function : calcImproperEnergy()
// ------------------------------------------------------------
//
// ============================================================
double mmPotential::calcImproperEnergy()
{
#ifdef DEBUG
        std::cout << "  mmPotential::calcImproperEnergy " << std::endl;
#endif
    this->improperEnergy = 0.0;
    return 0.0;
}

// ============================================================
// Function : calcNonBondedEnergy()
// ------------------------------------------------------------
//
// ============================================================
double mmPotential::calcNonBondedEnergy()
{
#ifdef DEBUG
        std::cout << "  mmPotential::calcNonBondedEnergy " << std::endl;
#endif
    this->nonBondedEnergy = 0.0;
    return 0.0;
}

// ============================================================
// Function : calcHBondEnergy()
// ------------------------------------------------------------
//
// ============================================================
double mmPotential::calcHBondEnergy()
{
#ifdef DEBUG
        std::cout << "  mmPotential::calcHBondEnergy " << std::endl;
#endif
    this->hBondEnergy = 0.0;
    return 0.0;
}

// ============================================================
// Function : getTotalEnergy()
// ------------------------------------------------------------
//
// ============================================================
double mmPotential::getTotalEnergy()
{
    return this->totalEnergy;
}

// ============================================================
// Function : getTotalVDWEnergy()
// ------------------------------------------------------------
//
// ============================================================
double mmPotential::getTotalVDWEnergy()
{
    return this->vdWEnergy + this->vdW14Energy;
}

// ============================================================
// Function : getTotalEleEnergy()
// ------------------------------------------------------------
//
// ============================================================
double mmPotential::getTotalEleEnergy()
{
    return this->eleEnergy + this->ele14Energy;
}

// ============================================================
// Function : getTotalEleEnergy()
// ------------------------------------------------------------
//
// ============================================================
void mmPotential::resetGMatrix()
{
    for (int i = 0; i < nAtoms*3-3; i++) {
      this->gradients[i] = 0;
    }
}

// ============================================================
// Function : setNonHGradsToZero()
// ------------------------------------------------------------
//
// ============================================================
void mmPotential::setNonHGradsToZero()
{
    for (int i = 0; i < nAtoms; i++) {
      if (this->symbols[i*2] != 'H') {
        for (int j = 0; j < 3; j++) {
          this->gradients[i*3+j] = 0;
        }
      }
    }
}

// ============================================================
// Function : getGradNorm()
// ------------------------------------------------------------
//
// ============================================================
void mmPotential::getGradNorm(double& gMax, int& gMaxAtom, double& gNorm)
{
    double gSum = 0;
    for (int j = 0; j < nAtoms*3; j++) {
      if (std::abs(gradients[j]) > gMax) {
        gMax = gradients[j];
        gSum += (gradients[j]*gradients[j]);
        gMaxAtom = j;
      }
    }
    gMaxAtom = (gMaxAtom-1)/3 + 1;
    gNorm = sqrt(gSum/double(nAtoms*3.0));
}

} // MTKpp namespace
