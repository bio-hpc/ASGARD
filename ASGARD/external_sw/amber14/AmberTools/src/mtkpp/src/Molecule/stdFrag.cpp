/*!
   \file stdFrag.cpp
   \brief Container for standard atoms, bonds, angles, torsions and impropers
   \author Martin Peters

   $Date: 2010/05/03 18:21:38 $
   $Revision: 1.21 $

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
#include <sstream>
#include "stdio.h"
#include "stdFrag.h"
#include "stdGroup.h"
#include "Log/errorHandler.h"
#include "Diagnostics/MTKException.h"

#include "Utils/vector3d.h"

namespace MTKpp
{

// ============================================================
// Function : stdFrag()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
stdFrag::stdFrag(stdGroup *parent):pParent(parent)
{
    this->itsName         = "";
    this->itsSymbol       = "";
    this->itsCharacter    = "";
    this->itsCode         = "";
    this->itsType         = "l";
    this->itsSymmetry     = "";
    this->pStdAtom        = 0;
    this->pStdBond        = 0;
    this->pStdLoop        = 0;
    this->pStdAlias       = 0;
    this->pStdImproper    = 0;
    this->pStdRing        = 0;
    this->pStdFeature     = 0;
    this->pStdFuncGroup   = 0;
    this->pStdConnTorsion = 0;
    this->pStdRotBond     = 0;
    this->adjMatrix       = 0;
    this->adjMatrixSize   = 0;
    this->atomSymbols     = 0;
    this->atomKinds       = 0;
    this->heavyAtomAdjMatrix = 0;
    this->heavyAtomAdjMatrixSize = 0;
    this->heavyAtomSymbols = 0;
    this->heavyAtomKinds   = 0;
}

// ============================================================
// Function : stdFrag()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
stdFrag::stdFrag(stdFrag *sf, stdGroup *parent)
{
    this->pParent         = parent;
    this->pStdAtom        = 0;
    this->pStdBond        = 0;
    this->pStdLoop        = 0;
    this->pStdAlias       = 0;
    this->pStdImproper    = 0;
    this->pStdRing        = 0;
    this->pStdFeature     = 0;
    this->pStdFuncGroup   = 0;
    this->pStdConnTorsion = 0;
    this->pStdRotBond     = 0;
    this->adjMatrix       = 0;
    this->adjMatrixSize   = 0;
    this->atomSymbols     = 0;
    this->atomKinds       = 0;
    this->heavyAtomAdjMatrix = 0;
    this->heavyAtomAdjMatrixSize = 0;
    this->heavyAtomSymbols = 0;
    this->heavyAtomKinds   = 0;

    if (sf) {
      this->itsName       = sf->getName();
      this->itsSymbol     = sf->getSymbol();
      this->itsCharacter  = sf->getCharacter();
      this->itsCode       = sf->getCode();
      this->itsType       = sf->getType();
      this->itsSymmetry   = sf->getSymmetry();

      this->setSubGraphs(sf->getSubGraphs());

      std::vector<stdAtom*> vStdAtoms = sf->getStdAtomList();
      for (stdAtomIterator c = vStdAtoms.begin(); c != vStdAtoms.end(); c++) {
        pStdAtom = *c;
        stdAtom* newStdAtom = this->addStdAtom(pStdAtom);
        if (!newStdAtom) {
          std::cout << " Error creating stdAtom in stdFrag " << std::endl;
        }
      }

      std::vector<stdBond*> vStdBonds = sf->getStdBondList();
      for (stdBondIterator c = vStdBonds.begin(); c != vStdBonds.end(); c++) {
        pStdBond = *c;
        stdBond* newStdBond = this->addStdBond(pStdBond);
        if (!newStdBond) {
          std::cout << " Error creating stdBond in stdFrag " << std::endl;
        }
      }

      std::vector<stdLoop*> vStdLoops = sf->getStdLoopList();
      for (stdLoopIterator c = vStdLoops.begin(); c != vStdLoops.end(); c++) {
        pStdLoop = *c;
        stdLoop* newStdLoop = this->addStdLoop(pStdLoop);
        if (!newStdLoop) {
          std::cout << " Error creating stdLoop in stdFrag " << std::endl;
        }
      }

      std::vector<stdAlias*> vStdAliass = sf->getStdAliasList();
      for (stdAliasIterator c = vStdAliass.begin(); c != vStdAliass.end(); c++) {
        pStdAlias = *c;
        stdAlias* newStdAlias = this->addStdAlias(pStdAlias);
        if (!newStdAlias) {
          std::cout << " Error creating stdAlias in stdFrag " << std::endl;
        }
      }

      std::vector<stdImproper*> vStdImpropers = sf->getStdImproperList();
      for (stdImproperIterator c = vStdImpropers.begin(); c != vStdImpropers.end(); c++) {
        pStdImproper = *c;
        stdImproper* newStdImproper = this->addStdImproper(pStdImproper);
        if (!newStdImproper) {
          std::cout << " Error creating stdImproper in stdFrag " << std::endl;
        }
      }

      std::vector<stdRing*> vStdRings = sf->getStdRingList();
      for (stdRingIterator c = vStdRings.begin(); c != vStdRings.end(); c++) {
        pStdRing = *c;
        stdRing* newStdRing = this->addStdRing(pStdRing);
        if (!newStdRing) {
          std::cout << " Error creating stdRing in stdFrag " << std::endl;
        }
      }

      std::vector<stdFeature*> vStdFeatures = sf->getStdFeatureList();
      for (stdFeatureIterator c = vStdFeatures.begin(); c != vStdFeatures.end(); c++) {
        pStdFeature = *c;
        stdFeature* newStdFeature = this->addStdFeature(pStdFeature);
        if (!newStdFeature) {
          std::cout << " Error creating stdFeature in stdFrag " << std::endl;
        }
      }

      std::vector<stdFuncGroup*> vStdFuncGroups = sf->getStdFuncGroupList();
      for (stdFuncGroupIterator c = vStdFuncGroups.begin(); c != vStdFuncGroups.end(); c++) {
        pStdFuncGroup = *c;
        stdFuncGroup* newStdFuncGroup = this->addStdFuncGroup(pStdFuncGroup);
        if (!newStdFuncGroup) {
          std::cout << " Error creating stdFuncGroup in stdFrag " << std::endl;
        }
      }

      std::vector<int> cPL = sf->getStdConnPtsList();
      for (unsigned int e = 0; e < cPL.size(); e++) {
        this->itsStdConnPtsList.push_back(cPL[e]);
      }
      //this->generateSimpleFP(); // fix getting atom number
    }
}

// ============================================================
// Function : ~stdFrag()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
stdFrag::~stdFrag()
{
    // stdAtoms
    for (stdAtomIterator c = this->itsStdAtomList.begin(); c != this->itsStdAtomList.end(); c++) {
      pStdAtom = *c;
      delete pStdAtom;
    }
    this->itsStdAtomList.clear();

    // stdBonds
    for (stdBondIterator c = this->itsStdBondList.begin(); c != this->itsStdBondList.end(); c++) {
      pStdBond = *c;
      delete pStdBond;
    }

    // stdLoops
    for (stdLoopIterator c = this->itsStdLoopList.begin(); c != this->itsStdLoopList.end(); c++) {
      pStdLoop = *c;
      delete pStdLoop;
    }

    // stdAlias
    for (stdAliasIterator c = this->itsStdAliasList.begin(); c != this->itsStdAliasList.end(); c++) {
      pStdAlias = *c;
      delete pStdAlias;
    }

    // stdImproper
    for (stdImproperIterator c = this->itsStdImproperList.begin(); c != this->itsStdImproperList.end(); c++) {
      pStdImproper = *c;
      delete pStdImproper;
    }

    // stdImproper
    for (stdRingIterator c = this->itsStdRingList.begin(); c != this->itsStdRingList.end(); c++) {
      pStdRing = *c;
      delete pStdRing;
    }

    // stdFeature
    for (stdFeatureIterator c = this->itsStdFeatureList.begin(); c != this->itsStdFeatureList.end(); c++) {
      pStdFeature = *c;
      delete pStdFeature;
    }

    // stdFuncGroup
    for (stdFuncGroupIterator c = this->itsStdFuncGroupList.begin(); c != this->itsStdFuncGroupList.end(); c++) {
      pStdFuncGroup = *c;
      delete pStdFuncGroup;
    }

    // stdConnTorsion
    for (stdConnTorsionIterator c = this->itsStdConnTorsionList.begin(); c != this->itsStdConnTorsionList.end(); c++) {
      pStdConnTorsion = *c;
      delete pStdConnTorsion;
    }

    // stdConnTorsion
    for (stdRotBondIterator c = this->itsStdRotBondList.begin(); c != this->itsStdRotBondList.end(); c++) {
      pStdRotBond = *c;
      delete pStdRotBond;
    }

    delete [] adjMatrix;
    delete [] atomSymbols;
}

// ============================================================
// Function : getParent()
// ------------------------------------------------------------
//
// ============================================================
stdGroup* stdFrag::getParent()
{
    return pParent;
}

// ============================================================
// Function : setName()
// ------------------------------------------------------------
//
// ============================================================
void stdFrag::setName(const std::string &name)
{
    this->itsName = name;
}

// ============================================================
// Function : getName()
// ------------------------------------------------------------
//
// ============================================================
std::string stdFrag::getName()
{
    return this->itsName;
}

// ============================================================
// Function : setSymbol()
// ------------------------------------------------------------
//
// ============================================================
void stdFrag::setSymbol(const std::string &symbol)
{
    if (symbol.size() == 3) {
      this->itsSymbol = symbol;
    }
    else {
      std::cout << " Length of Fragment Symbol (" << symbol 
                << ") = " << symbol.size()
                << " --> Must be 3." << std::endl;
    }
}

// ============================================================
// Function : getSymbol()
// ------------------------------------------------------------
//
// ============================================================
std::string stdFrag::getSymbol()
{
    return this->itsSymbol;
}

// ============================================================
// Function : setCode()
// ------------------------------------------------------------
//
// ============================================================
void stdFrag::setCode(const std::string &code)
{
    if (code.size() == 8) {
      this->itsCode = code;
    }
    else {
      std::cout << " stdFrag::setCode for " << code
                << ". Length of Fragment Code = " << code.size()
                << " --> Must be 8." << std::endl;
    }
}

// ============================================================
// Function : getCode()
// ------------------------------------------------------------
//
// ============================================================
std::string stdFrag::getCode()
{
    return this->itsCode;
}

// ============================================================
// Function : setCharacter()
// ------------------------------------------------------------
//
// ============================================================
void stdFrag::setCharacter(const std::string &c)
{
    if (c.size() == 1) {
      this->itsCharacter = c;
    }
    else {
      std::cout << " stdFrag::setCharacter for |" << this->itsSymbol
                << "| Length of Fragment Character = " << c.size()
                << " --> Must be 1." << std::endl;
      exit(0);
    }
}

// ============================================================
// Function : getCharacter()
// ------------------------------------------------------------
//
// ============================================================
std::string stdFrag::getCharacter()
{
    return this->itsCharacter;
}

// ============================================================
// Function : setType()
// ------------------------------------------------------------
//
// ============================================================
void stdFrag::setType(const std::string &t)
{
    if (t.size() == 1) {
      this->itsType = t;
    }
    else {
      std::cout << " Length of fragment type = " << t.size()
                << " --> Must be 1." << std::endl;
    }
}

// ============================================================
// Function : getType()
// ------------------------------------------------------------
//
// ============================================================
std::string stdFrag::getType()
{
    return this->itsType;
}

// ============================================================
// Function : setSymmetry()
// ------------------------------------------------------------
//
// ============================================================
void stdFrag::setSymmetry(const std::string &t)
{
    this->itsSymmetry = t;
}

// ============================================================
// Function : getSymmetry()
// ------------------------------------------------------------
//
// ============================================================
std::string stdFrag::getSymmetry()
{
    return this->itsSymmetry;
}

// ============================================================
// Function : setSubGraphs()
// ------------------------------------------------------------
//
// ============================================================
void stdFrag::setSubGraphs(const std::vector<std::string> &g)
{
    for (unsigned int i = 0; i < g.size(); i++) {
      if (g[i].size() == 3) {
        this->itsSubGraphs.push_back(g[i]);
      }
      else {
        std::cout << " User must use 3L code in subGraphs tag " << std::endl;
      }
    }
}

// ============================================================
// Function : getSubGraphs()
// ------------------------------------------------------------
//
// ============================================================
std::vector<std::string> stdFrag::getSubGraphs()
{
    return this->itsSubGraphs;
}

// ============================================================
// Function : getSubGraphStr()
// ------------------------------------------------------------
//
// ============================================================
std::string stdFrag::getSubGraphStr()
{
    std::string sGstr = "";
    for (unsigned int i = 0; i < this->itsSubGraphs.size(); i++) {
      sGstr+=itsSubGraphs[i];
    }
    return sGstr;
}

// ============================================================
// Function : addStdAtom()
// ------------------------------------------------------------
//
// ============================================================
stdAtom* stdFrag::addStdAtom()
{
    pStdAtom = new stdAtom();

    pStdAtom->kind = 0;

    this->itsStdAtomList.push_back(pStdAtom);
    return pStdAtom;
}

// ============================================================
// Function : addStdAtom(stdAtom*)
// ------------------------------------------------------------
//
// ============================================================
stdAtom* stdFrag::addStdAtom(stdAtom* s)
{
    pStdAtom = new stdAtom();
    this->itsStdAtomList.push_back(pStdAtom);

    pStdAtom->identity = s->identity;
    pStdAtom->index = s->index;
    pStdAtom->type = s->type;
    pStdAtom->chain = s->chain;
    pStdAtom->atmCharge = s->atmCharge;
    pStdAtom->bond12 = s->bond12;
    pStdAtom->bondLength = s->bondLength;
    pStdAtom->bond13 = s->bond13;
    pStdAtom->bondAngle = s->bondAngle;
    pStdAtom->bond14 = s->bond14;
    pStdAtom->bondTorsion = s->bondTorsion;
    //pStdAtom->atNum = s->atNum;
    //pStdAtom->atSymbol = s->atSymbol;
    //pStdAtom->hybridization = s->hybridization;
    pStdAtom->kind = s->kind;

    return pStdAtom;
}

// ============================================================
// Function : addStdBond()
// ------------------------------------------------------------
//
// ============================================================
stdBond* stdFrag::addStdBond()
{
    pStdBond = new stdBond();
    this->itsStdBondList.push_back(pStdBond);
    return pStdBond;
}

// ============================================================
// Function : addStdBond(stdBond*)
// ------------------------------------------------------------
//
// ============================================================
stdBond* stdFrag::addStdBond(stdBond* b)
{
    pStdBond = new stdBond();
    this->itsStdBondList.push_back(pStdBond);

    pStdBond->atom1 = b->atom1;
    pStdBond->atom2 = b->atom2;
    pStdBond->type = b->type;
    pStdBond->stereo = b->stereo;
    pStdBond->topology = b->topology;
    pStdBond->length = b->length;
    return pStdBond;
}

// ============================================================
// Function : addStdImproper()
// ------------------------------------------------------------
//
// ============================================================
stdImproper* stdFrag::addStdImproper()
{
    pStdImproper = new stdImproper();
    this->itsStdImproperList.push_back(pStdImproper);
    return pStdImproper;
}

// ============================================================
// Function : addStdImproper(stdImproper*)
// ------------------------------------------------------------
//
// ============================================================
stdImproper* stdFrag::addStdImproper(stdImproper* i)
{
    pStdImproper = new stdImproper();
    this->itsStdImproperList.push_back(pStdImproper);

    pStdImproper->atom1 = i->atom1;
    pStdImproper->atom2 = i->atom2;
    pStdImproper->atom3 = i->atom3;
    pStdImproper->atom4 = i->atom4;
    return pStdImproper;
}

// ============================================================
// Function : addStdLoop()
// ------------------------------------------------------------
//
// ============================================================
stdLoop* stdFrag::addStdLoop()
{
    pStdLoop = new stdLoop();
    this->itsStdLoopList.push_back(pStdLoop);
    return pStdLoop;
}

// ============================================================
// Function : addStdLoop()
// ------------------------------------------------------------
//
// ============================================================
stdLoop* stdFrag::addStdLoop(stdLoop* l)
{
    pStdLoop = new stdLoop();
    this->itsStdLoopList.push_back(pStdLoop);
    pStdLoop->atom1 = l->atom1;
    pStdLoop->atom2 = l->atom2;
    pStdLoop->type = l->type;
    pStdLoop->stereo = l->stereo;
    return pStdLoop;
}

// ============================================================
// Function : getStdLoop()
// ------------------------------------------------------------
//
// ============================================================
stdLoop* stdFrag::getStdLoop(const int &index1, const int &index2)
{
    for (stdLoopIterator c = this->itsStdLoopList.begin(); c != this->itsStdLoopList.end(); c++) {
      pStdLoop = *c;
      if ((pStdLoop->atom1 == index1 && pStdLoop->atom2 == index2) ||
         (pStdLoop->atom1 == index2 && pStdLoop->atom2 == index1) ) {
        return pStdLoop;
      }
    }
    return 0;
}

// ============================================================
// Function : getStdLoop()
// ------------------------------------------------------------
//
// ============================================================
stdLoop* stdFrag::getStdLoop(stdAtom* pAt1, stdAtom* pAt2)
{
    for (stdLoopIterator c = this->itsStdLoopList.begin(); c != this->itsStdLoopList.end(); c++) {
      pStdLoop = *c;
      if ((pStdLoop->atom1 == pAt1->index && pStdLoop->atom2 == pAt2->index) ||
          (pStdLoop->atom1 == pAt2->index && pStdLoop->atom2 == pAt1->index) ) {
        return pStdLoop;
      }
    }
    return 0;
}

// ============================================================
// Function : addStdAlias()
// ------------------------------------------------------------
//
// ============================================================
stdAlias* stdFrag::addStdAlias()
{
    pStdAlias = new stdAlias();
    this->itsStdAliasList.push_back(pStdAlias);
    return pStdAlias;
}

// ============================================================
// Function : addStdAlias(stdAlias*)
// ------------------------------------------------------------
//
// ============================================================
stdAlias* stdFrag::addStdAlias(stdAlias* a)
{
    pStdAlias = new stdAlias();
    this->itsStdAliasList.push_back(pStdAlias);
    pStdAlias->atom1 = a->atom1;
    pStdAlias->atom2 = a->atom2;
    return pStdAlias;
}

// ============================================================
// Function : addStdRing()
// ------------------------------------------------------------
//
// ============================================================
stdRing* stdFrag::addStdRing()
{
    pStdRing = new stdRing();
    this->itsStdRingList.push_back(pStdRing);
    return pStdRing;
}

// ============================================================
// Function : addStdRing(stdRing*)
// ------------------------------------------------------------
//
// ============================================================
stdRing* stdFrag::addStdRing(stdRing* r)
{
    pStdRing = new stdRing();
    this->itsStdRingList.push_back(pStdRing);

    for (unsigned int i = 0; i < r->atoms.size(); i++) {
      pStdRing->atoms.push_back(r->atoms[i]);
    }
    pStdRing->size = r->size;
    pStdRing->planar = r->planar;
    pStdRing->aromatic = r->aromatic;
    pStdRing->hetero = r->hetero;
    pStdRing->nHetero = r->nHetero;
    pStdRing->nNitrogen = r->nNitrogen;
    pStdRing->nOxygen = r->nOxygen;
    pStdRing->nSulfur = r->nSulfur;
    return pStdRing;
}

// ============================================================
// Function : addStdFeature()
// ------------------------------------------------------------
//
// ============================================================
stdFeature* stdFrag::addStdFeature()
{
    pStdFeature = new stdFeature();
    this->itsStdFeatureList.push_back(pStdFeature);
    return pStdFeature;
}

// ============================================================
// Function : addStdFeature(stdFeature*)
// ------------------------------------------------------------
//
// ============================================================
stdFeature* stdFrag::addStdFeature(stdFeature* f)
{
    pStdFeature = new stdFeature();
    this->itsStdFeatureList.push_back(pStdFeature);
    pStdFeature->name = f->name;
    for (unsigned int i = 0; i < f->atoms.size(); i++) {
      pStdFeature->atoms.push_back(f->atoms[i]);
    }
    return pStdFeature;
}

// ============================================================
// Function : hasStdFeature()
// ------------------------------------------------------------
//
// ============================================================
bool stdFrag::hasStdFeature(stdAtom* a, std::string f)
{
    for (unsigned int i = 0; i < this->itsStdFeatureList.size(); i++) {
      if (this->itsStdFeatureList[i]->name == f) {
        for (unsigned int j = 0; j < this->itsStdFeatureList[i]->atoms.size();
             j++) {
          if (a->index == this->itsStdFeatureList[i]->atoms[j]) {
            return 1;
          }
        }
      }
    }
    return 0;
}

// ============================================================
// Function : addStdFuncGroup()
// ------------------------------------------------------------
//
// ============================================================
stdFuncGroup* stdFrag::addStdFuncGroup()
{
    pStdFuncGroup = new stdFuncGroup();
    this->itsStdFuncGroupList.push_back(pStdFuncGroup);
    return pStdFuncGroup;
}

// ============================================================
// Function : addStdFuncGroup(stdFuncGroup*)
// ------------------------------------------------------------
//
// ============================================================
stdFuncGroup* stdFrag::addStdFuncGroup(stdFuncGroup* f)
{
    pStdFuncGroup = new stdFuncGroup();
    this->itsStdFuncGroupList.push_back(pStdFuncGroup);
    pStdFuncGroup->groupName = f->groupName;
    pStdFuncGroup->fragName = f->fragName;
    for (unsigned int i = 0; i < f->atoms.size(); i++) {
      pStdFuncGroup->atoms.push_back(f->atoms[i]);
    }
    return pStdFuncGroup;
}

// ============================================================
// Function : addStdConnPts()
// ------------------------------------------------------------
//
// ============================================================
void stdFrag::addStdConnPts(std::vector<int> &v)
{
    for (unsigned int i = 0; i < v.size(); i++) {
      this->itsStdConnPtsList.push_back(v[i]);
    }
}

// ============================================================
// Function : addStdConnTorsion()
// ------------------------------------------------------------
//
// ============================================================
stdConnTorsion* stdFrag::addStdConnTorsion()
{
    pStdConnTorsion = new stdConnTorsion();
    this->itsStdConnTorsionList.push_back(pStdConnTorsion);
    return pStdConnTorsion;
}

// ============================================================
// Function : addStdRotBond()
// ------------------------------------------------------------
//
// ============================================================
stdRotBond* stdFrag::addStdRotBond()
{
    pStdRotBond = new stdRotBond();
    this->itsStdRotBondList.push_back(pStdRotBond);
    return pStdRotBond;
}

// ============================================================
// Function : getStdAtomList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<stdAtom*> stdFrag::getStdAtomList()
{
    return this->itsStdAtomList;
}

// ============================================================
// Function : numStdAtoms()
// ------------------------------------------------------------
//
// ============================================================
int stdFrag::numStdAtoms()
{
    return this->itsStdAtomList.size();
}

// ============================================================
// Function : numStdHeavyAtoms()
// ------------------------------------------------------------
//
// ============================================================
int stdFrag::numStdHeavyAtoms()
{
    int nHeavies = 0;
    unsigned int nAts = this->itsStdAtomList.size();
    for (unsigned int i = 0; i < nAts; i++) {
      if (this->itsStdAtomList[i]->atSymbol != "H") {
        nHeavies++;
      }
    }
    return nHeavies;
}

// ============================================================
// Function : getStdBondList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<stdBond*> stdFrag::getStdBondList()
{
    return this->itsStdBondList;
}

// ============================================================
// Function : getStdLoopList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<stdLoop*> stdFrag::getStdLoopList()
{
    return this->itsStdLoopList;
}

// ============================================================
// Function : getStdAliasList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<stdAlias*> stdFrag::getStdAliasList()
{
    return this->itsStdAliasList;
}

// ============================================================
// Function : getStdImproperList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<stdImproper*> stdFrag::getStdImproperList()
{
    return this->itsStdImproperList;
}

// ============================================================
// Function : getStdRingList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<stdRing*> stdFrag::getStdRingList()
{
    return this->itsStdRingList;
}

// ============================================================
// Function : getStdFeatureList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<stdFeature*> stdFrag::getStdFeatureList()
{
    return this->itsStdFeatureList;
}

// ============================================================
// Function : getStdFuncGroupList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<stdFuncGroup*> stdFrag::getStdFuncGroupList()
{
    return this->itsStdFuncGroupList;
}

// ============================================================
// Function : getStdConnPtsList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<int> stdFrag::getStdConnPtsList()
{
    return this->itsStdConnPtsList;
}

// ============================================================
// Function : getNumStdConnPtsList()
// ------------------------------------------------------------
//
// ============================================================
int stdFrag::getNumStdConnPtsList()
{
    return this->itsStdConnPtsList.size();
}

// ============================================================
// Function : getStdConnTorList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<stdConnTorsion*> stdFrag::getStdConnTorList()
{
    return this->itsStdConnTorsionList;
}

// ============================================================
// Function : getStdRotBondList()
// ------------------------------------------------------------
//
// ============================================================
std::vector<stdRotBond*> stdFrag::getStdRotBondList()
{
    return this->itsStdRotBondList;
}

// ============================================================
// Function : getStdAtom()
// ------------------------------------------------------------
//
// ============================================================
stdAtom* stdFrag::getStdAtom(const int &index)
{
    if (index > 0) {
      for (stdAtomIterator c = this->itsStdAtomList.begin();
           c != this->itsStdAtomList.end(); c++) {
        pStdAtom = *c;
        if (pStdAtom->index == index) {
          return pStdAtom;
        }
      }
    }
    else if (index == -1) {
      std::cout << " stdFrag::getStdAtom -1 " << std::endl;
      return 0;
    }
    else if (index == -4) {
      std::cout << " stdFrag::getStdAtom -4 " << std::endl;
      return 0;
    }
    else {
      std::vector<stdAtom*> v2(itsStdAtomList.size());
      std::reverse_copy(itsStdAtomList.begin(),
                        this->itsStdAtomList.end(), v2.begin());
      int counter = 0;
      for (stdAtomIterator c = v2.begin(); c != v2.end(); c++) {
        pStdAtom = *c;
        if (pStdAtom->chain == "M") {
          counter++;
        }
        if (counter == abs(index)) {
          return pStdAtom;
        }
      }
    }
    return 0;
}

// ============================================================
// Function : getStdAtom()
// ------------------------------------------------------------
//
// ============================================================
stdAtom* stdFrag::getStdAtom(const std::string &name)
{
    for (stdAtomIterator c = this->itsStdAtomList.begin();
         c != this->itsStdAtomList.end(); c++) {
      pStdAtom = *c;
      if (pStdAtom->identity == name) {
        return pStdAtom;
      }
    }

    std::string al = "";
    al = this->getAlias(name);
    if (al != "") {
      for (stdAtomIterator c = this->itsStdAtomList.begin();
           c != this->itsStdAtomList.end(); c++) {
        pStdAtom = *c;
        if (pStdAtom->identity == al) {
          return pStdAtom;
        }
      }
    }
    return 0;
}

// ============================================================
// Function : hasStdAtom()
// ------------------------------------------------------------
//
// ============================================================
bool stdFrag::hasStdAtom(const std::string &name)
{
    if (name.size() != 4) {
/*
      std::cout << " Error in stdFrag::hasStdAtom " << name
                << " stdAtom name must have 4 characters ... exiting "
                << std::endl;
      exit(0);
*/
      std::stringstream ss;
      ss << " Error in stdFrag::hasStdAtom " << name
                << " stdAtom name must have 4 characters ... exiting "
                << std::endl;
      std::cout << ss.str();
      throw MTKException(ss.str());
    }
    for (stdAtomIterator c = this->itsStdAtomList.begin();
         c != this->itsStdAtomList.end(); c++) {
      pStdAtom = *c;
      if (pStdAtom->identity == name) {
        return true;
      }
    }
    std::string al = "";
    al = this->getAlias(name);
    if (al != "") {
      for (stdAtomIterator c = this->itsStdAtomList.begin();
           c != this->itsStdAtomList.end(); c++) {
        pStdAtom = *c;
        if (pStdAtom->identity == al) {
          return true;
        }
      }
    }
    return false;
}

// ============================================================
// Function : getBondedStdAtoms()
// ------------------------------------------------------------
//
// ============================================================
std::vector<stdAtom*> stdFrag::getBondedStdAtoms(stdAtom* sAt)
{
    std::vector<stdAtom*> bondedAtoms;
    for (stdAtomIterator c = this->itsStdAtomList.begin();
         c != this->itsStdAtomList.end(); c++) {
      pStdAtom = *c;
      if (this->hasStdBond(sAt->index, pStdAtom->index)) {
        bondedAtoms.push_back(pStdAtom);
      }
    }
    return bondedAtoms;
}

// ============================================================
// Function : getStdBond()
// ------------------------------------------------------------
//
// ============================================================
stdBond* stdFrag::getStdBond(const int &index1, const int &index2)
{
    for (stdBondIterator c = this->itsStdBondList.begin();
         c != this->itsStdBondList.end(); c++) {
      pStdBond = *c;
      if ((pStdBond->atom1 == index1 &&
           pStdBond->atom2 == index2) ||
          (pStdBond->atom1 == index2 &&
           pStdBond->atom2 == index1) ) {
        return pStdBond;
      }
    }
    return 0;
}

// ============================================================
// Function : getStdBond()
// ------------------------------------------------------------
//
// ============================================================
stdBond* stdFrag::getStdBond(const std::string &at1, const std::string &at2)
{
    stdAtom* pAt1 = getStdAtom(at1);
    stdAtom* pAt2 = getStdAtom(at2);

    for (stdBondIterator c = this->itsStdBondList.begin();
         c != this->itsStdBondList.end(); c++) {
      pStdBond = *c;

      if ((pStdBond->atom1 == pAt1->index &&
           pStdBond->atom2 == pAt2->index) ||
          (pStdBond->atom1 == pAt2->index &&
           pStdBond->atom2 == pAt1->index) ) {
        return pStdBond;
      }
    }
    return 0;
}

// ============================================================
// Function : getStdBond()
// ------------------------------------------------------------
//
// ============================================================
stdBond* stdFrag::getStdBond(stdAtom* pAt1, stdAtom* pAt2)
{
    for (stdBondIterator c = this->itsStdBondList.begin();
         c != this->itsStdBondList.end(); c++) {
      pStdBond = *c;
      if ((pStdBond->atom1 == pAt1->index && pStdBond->atom2 == pAt2->index) ||
          (pStdBond->atom1 == pAt2->index && pStdBond->atom2 == pAt1->index) ) {
        return pStdBond;
      }
    }
    return 0;
}

// ============================================================
// Function : hasStdBond()
// ------------------------------------------------------------
//
// ============================================================
bool stdFrag::hasStdBond(const int &index1, const int &index2)
{
    for (stdBondIterator c = this->itsStdBondList.begin();
         c != this->itsStdBondList.end(); c++) {
      pStdBond = *c;
      if ((pStdBond->atom1 == index1 && pStdBond->atom2 == index2) ||
          (pStdBond->atom1 == index2 && pStdBond->atom2 == index1) ) {
        return true;
      }
    }
    for (stdLoopIterator c = this->itsStdLoopList.begin();
         c != this->itsStdLoopList.end(); c++) {
      pStdLoop = *c;
      if ((pStdLoop->atom1 == index1 && pStdLoop->atom2 == index2) ||
          (pStdLoop->atom1 == index2 && pStdLoop->atom2 == index1) ) {
        return true;
      }
    }
    return false;
}

// ============================================================
// Function : numStdBonds()
// ------------------------------------------------------------
//
// ============================================================
int stdFrag::numStdBonds()
{
    return this->itsStdBondList.size();
}

// ============================================================
// Function : numStdBonds()
// ------------------------------------------------------------
//
// ============================================================
int stdFrag::numStdBonds(stdAtom* pStdAt)
{
    int i = pStdAt->index;
    int counter = 0;
    for (stdBondIterator c = this->itsStdBondList.begin();
         c != this->itsStdBondList.end(); c++) {
      pStdBond = *c;
      if ((pStdBond->atom1 == i) or (pStdBond->atom2 == i)) {
        counter++;
      }
    }
    for (stdLoopIterator c = this->itsStdLoopList.begin();
         c != this->itsStdLoopList.end(); c++) {
      pStdLoop = *c;
      if ((pStdLoop->atom1 == i) or (pStdLoop->atom2 == i)) {
        counter++;
      }
    }
    return counter;
}

// ============================================================
// Function : getAlias()
// ------------------------------------------------------------
//
// ============================================================
std::string stdFrag::getAlias(const std::string &name)
{
    for (stdAliasIterator c = this->itsStdAliasList.begin();
         c != this->itsStdAliasList.end(); c++) {
      pStdAlias = *c;
      if (pStdAlias->atom2 == name) {
        return pStdAlias->atom1;
      }
    }
    return "";
}

// ============================================================
// Function : generateCoordinates()
// ------------------------------------------------------------
// Computes cartesian coordinates from internal coordinates.
// ============================================================
int stdFrag::generateCoordinates()
{
    std::string errorMessage = "\n";

    int nAtoms = this->itsStdAtomList.size();

    vector3d* dum1 = new vector3d(0.0);
    vector3d* dum2 = new vector3d(0.0);
    dum2->setX(1.449);
    vector3d* dum3 = new vector3d(0.0);
    dum3->setX(1.449);
    double a = 111.1*DEG2RAD - 90.0*DEG2RAD;
    double xref = 0.0;
    dum3->setX(xref + 1.522*sin(a));
    dum3->setY(1.522*cos(a));

/*
    // Coordinates for atoms 1, 2, and 3 are generated separately.
    if (nAtoms > 0) {
      // atom 1
      pCoords = new vector3d(0.0);
      this->itsCoords.push_back(pCoords);
      if (nAtoms == 1) return 0;

      // atom 2
      pStdAtom = this->getStdAtom(2);
      pCoords = new vector3d(0.0);
      pCoords->setX(pStdAtom->bondLength);
      this->itsCoords.push_back(pCoords);
      if (nAtoms == 2) return 0;

      // atom 3
      pCoords = new vector3d(0.0);
      pStdAtom = this->getStdAtom(3);
      int bondPartner = pStdAtom->bond12;
      int anglePartner = pStdAtom->bond13;
      double bondLength = pStdAtom->bondLength;
      double ang = pStdAtom->bondAngle;
      double a = ang*DEG2RAD - 90.0*DEG2RAD;
      double xref = this->itsCoords[1]->getX();

      // The 3-1-2 numbering scheme has a different convention:
      if ((bondPartner == 1) and (anglePartner == 2)) {
        a = -a;
        xref = 0.0;
      }

      pCoords->setX(xref + bondLength*sin(a));
      pCoords->setY(bondLength*cos(a));
      this->itsCoords.push_back(pCoords);
      if (nAtoms == 3) return 0;
    }
*/

    // Now begin general procedure for determining cartesian coordinates
    // from internal coordinates.
    for (int i = 0; i < nAtoms; i++) {
      pCoords = new vector3d(0.0);
      pStdAtom = this->getStdAtom(i+1);

      int bondPartner    = pStdAtom->bond12;
      vector3d* bondPartnerCoords = 0;
      if (bondPartner < 0) {
        if (bondPartner == -1) bondPartnerCoords = dum3;
        if (bondPartner == -2) bondPartnerCoords = dum2;
        if (bondPartner == -3) bondPartnerCoords = dum1;
      }
      else {
        bondPartnerCoords = this->itsCoords[bondPartner-1];
      }

      int anglePartner   = pStdAtom->bond13;
      vector3d* anglePartnerCoords = 0;
      if (anglePartner < 0) {
        if (anglePartner == -1) anglePartnerCoords = dum3;
        if (anglePartner == -2) anglePartnerCoords = dum2;
        if (anglePartner == -3) anglePartnerCoords = dum1;
      }
      else {
        anglePartnerCoords = this->itsCoords[anglePartner-1];
      }

      int torsionPartner = pStdAtom->bond14;
      vector3d* torsionPartnerCoords = 0;
      if (torsionPartner < 0) {
        if (torsionPartner == -1) torsionPartnerCoords = dum3;
        if (torsionPartner == -2) torsionPartnerCoords = dum2;
        if (torsionPartner == -3) torsionPartnerCoords = dum1;
      }
      else {
        torsionPartnerCoords = this->itsCoords[torsionPartner-1];
      }

      double bondLength  = pStdAtom->bondLength;
      double angleSize   = pStdAtom->bondAngle;
      double torsionSize = pStdAtom->bondTorsion;

      buildCoord((*pCoords), (*bondPartnerCoords),
                 (*anglePartnerCoords), (*torsionPartnerCoords),
                 bondLength, angleSize*DEG2RAD, torsionSize*DEG2RAD);

      this->itsCoords.push_back(pCoords);

      char temp[100];
      sprintf(temp,"ATOM  %5d %-4.4s %-3.3s %1s%4d%1s   %8.3f%8.3f%8.3f\n",
              i+1,(pStdAtom->identity.c_str()),
              "XXX",
              "A",
              (1),
              "A",
              this->itsCoords[i]->getX(),this->itsCoords[i]->getY(),this->itsCoords[i]->getZ());

      errorMessage += temp;
    }
    errorLogger.throwError("stdFrag::generateCoordinates", errorMessage, INFO);

    return 0;
}

/////////
// ============================================================
// Function : generateCoordinates()
// ------------------------------------------------------------
// Computes cartesian coordinates from internal coordinates.
// ============================================================
int stdFrag::generateCoordinates(vector3d* bd, vector3d* ag,
                                 vector3d* tr, const int &forward)
{
    int nAtoms = this->itsStdAtomList.size();
    std::vector<vector3d*> coords(nAtoms);
    int connTable[nAtoms][3];
    double bdAgTors[nAtoms][3];

    int negTors = 0;
    for (unsigned int i = 0; i < this->itsStdAtomList.size(); i++) {
      if (itsStdAtomList[i]->bond14 < 0) negTors++;
    }

    double negativeTorsions[negTors-3];
    if (negTors == 4) {
      negativeTorsions[0] = 0.0;
    }
    else if (negTors == 5) {
      negativeTorsions[0] =  60.0;
      negativeTorsions[1] = 300.0;
    }

    // If we want to start from the end of the standard residue, we need
    // to make a new connection table
    std::vector<stdAtom*> newOrdering;
    //std::vector<stdAtom*>::iterator stdAtomIter;
    int b = -1;
    int a = -2;
    int t = -3;
    int ang1 = 1;
    int angM1 = 1;
    int tor1 = 0;
    int torM1 = 1;
    int torM2 = 1;
    int cT = 0;

    if (!forward) {
      // setup bonds, angles, and torsions to previous 3 atoms
      bdAgTors[0][0] = 1.335; bdAgTors[0][1] = 116.6; bdAgTors[0][2] = 180.0;
                              bdAgTors[1][1] = 120.0; bdAgTors[1][2] = 180.0;
                                                      bdAgTors[2][2] = 180.0;

      int atm = 0;
      for (int i = (itsStdAtomList.size()-1); i > -1; --i) {
        if (itsStdAtomList[i]->chain == "M") {
          if (b == 0) b++;
          if (a == 0) a++;
          if (t == 0) t++;
          newOrdering.push_back(itsStdAtomList[i]);
          connTable[atm][0] = b;
          connTable[atm][1] = a;
          connTable[atm][2] = t;
          atm++;
          b++;
          a++;
          t++;
          if (itsStdAtomList[i]->bond12 > 0) {
            bdAgTors[atm  ][0] = this->itsStdAtomList[i]->bondLength;
          }
          if (itsStdAtomList[i]->bond13 > 0) {
            bdAgTors[atm+1][1] = this->itsStdAtomList[i]->bondAngle;
          }
          if (itsStdAtomList[i]->bond14 > 0) {
            bdAgTors[atm+2][2] = this->itsStdAtomList[i]->bondTorsion;
          }
        }
      }

      for (int i = (itsStdAtomList.size()-1); i > -1; --i) {
        if (itsStdAtomList[i]->chain != "M") {
          if (b == 0) b++;
          if (a == 0) a++;
          if (t == 0) t++;
          newOrdering.push_back(itsStdAtomList[i]);

          stdAtom* bondAtom = 0;
          if (itsStdAtomList[i]->bond12 > 0) {
            bondAtom = this->getStdAtom(itsStdAtomList[i]->bond12);
            bdAgTors[atm][0] = this->itsStdAtomList[i]->bondLength;
          }
          else {
            b = this->itsStdAtomList[i]->bond12;
          }
          for (unsigned int x = 0; x < newOrdering.size(); x++) {
            if (newOrdering[x] == bondAtom) {
              b = x+1;
            }
          }

          stdAtom* angleAtom = 0;
          if (b > 2) {
            if (itsStdAtomList[i]->bond13 > 0) {
              angleAtom = this->getStdAtom(itsStdAtomList[i]->bond13);
              bdAgTors[atm][1] = this->itsStdAtomList[i]->bondAngle;
            }
            else {
              a = this->itsStdAtomList[i]->bond13;
            }
            for (unsigned int x = 0; x < newOrdering.size(); x++) {
              if (newOrdering[x] == angleAtom) {
                a = x+1;
              }
            }
          }
          else {
            if (b == 1) {
              a = -1;
              angM1++;
            }
            if (b == 2) {
              a =  1;
              ang1++;
            }
          }

          stdAtom* torsionAtom = 0;
          if (a > 2) {
            if (itsStdAtomList[i]->bond14 > 0) {
              torsionAtom = this->getStdAtom(itsStdAtomList[i]->bond14);
              bdAgTors[atm][2] = this->itsStdAtomList[i]->bondTorsion;
            }
            else {
              t = this->itsStdAtomList[i]->bond14;
              bdAgTors[atm][2] = negativeTorsions[cT];
              cT++;
            }
            for (unsigned int x = 0; x < newOrdering.size(); x++) {
              if (newOrdering[x] == torsionAtom) {
                t = x+1;
              }
            }
          }
          else {
            if (a == -1) {
              t = -2;
              torM2++;
            }
            if (a ==  1) {
              t = -1;
              torM1++;
            }
            if (a ==  2) {
              t =  1;
              tor1++;
            }
          }

          connTable[atm][0] = b;
          connTable[atm][1] = a;
          connTable[atm][2] = t;
          atm++;
        }
      }

      double sp2Angles[2] = {0.0, 180.0};
      double sp3Angles[3] = {180.0, 60.0, 300.0};
      double sp2Angle = 120.0;

      int cSp3Angles1 = 0;
      int cSp3Angles2 = 0;
      int cSp2Angles = 0;
      for (int i = 2; i < nAtoms; i++) {
        if (connTable[i][0] == 2) {
          if (connTable[i][1] == 1) {
            if (ang1 == 3) { // means atom 2 is sp3
              bdAgTors[i][1] = 109.5;
              bdAgTors[i][2] = sp3Angles[cSp3Angles2];
              cSp3Angles2++;
            }
            else if (ang1 == 2) { // means atom 2 is sp2
              bdAgTors[i][1] = sp2Angle;
              bdAgTors[i][2] = sp2Angles[cSp2Angles];
              cSp2Angles++;
            }
          }
        }
        else if (connTable[i][0] == 1) {
          if (connTable[i][1] == -1) {
            if (angM1 == 3) { // means atom1 is sp3
              bdAgTors[i][1] = 109.5;
              bdAgTors[i][2] = sp3Angles[cSp3Angles1];
              cSp3Angles1++;
            }
            else { // means atom 1 is sp2
              bdAgTors[i][1] = sp2Angle;
              if (torM2 == 2) {
                bdAgTors[i][2] = 0.0;
              }
            }
          }
        }
      }

/*
      for (int i = 0; i < nAtoms; i++) {
        std::cout << this->itsStdAtomList[i]->identity << " "
                  << this->itsStdAtomList[i]->bond12 << " " 
                  << this->itsStdAtomList[i]->bond13 << " "
                  << this->itsStdAtomList[i]->bond14 << " "
                  << this->itsStdAtomList[i]->bondLength << " "
                  << this->itsStdAtomList[i]->bondAngle << " "
                  << this->itsStdAtomList[i]->bondTorsion << " \n" << std::endl;
      }

      for (int i = 0; i < nAtoms; i++) {
        std::cout << newOrdering[i]->identity << " ";
        for (int j = 0; j < 3; j++) {
          std::cout << connTable[i][j] << " ";
        }
        for (int j = 0; j < 3; j++) {
          std::cout << bdAgTors[i][j] << " ";
        }
        std::cout << " " << std::endl;
      }
*/
      std::cout << " exit in stdFrag ... not sure why? " << std::endl;
      throw MTKException(" exit in stdFrag ... not sure why? ");

      //exit(0);
///////////////////////////////
    }
    // forward
    else {
      for (unsigned int i = 0; i < this->itsStdAtomList.size(); i++) {
        newOrdering.push_back(itsStdAtomList[i]);
        connTable[i][0] = this->itsStdAtomList[i]->bond12;
        connTable[i][1] = this->itsStdAtomList[i]->bond13;
        connTable[i][2] = this->itsStdAtomList[i]->bond14;

        bdAgTors[i][0] = this->itsStdAtomList[i]->bondLength;
        bdAgTors[i][1] = this->itsStdAtomList[i]->bondAngle;
        bdAgTors[i][2] = this->itsStdAtomList[i]->bondTorsion;
      }
    }

    // Now begin general procedure for determining cartesian coordinates
    // from internal coordinates.

    for (int i = 0; i < nAtoms; i++) {

      pCoords = new vector3d(0.0);
      pStdAtom = newOrdering[i];

      int bondPartner = connTable[i][0];
      vector3d* bondPartnerCoords = 0;

      if (bondPartner < 0) {
        if (bondPartner == -1) bondPartnerCoords = bd;
        if (bondPartner == -2) bondPartnerCoords = ag;
        if (bondPartner == -3) bondPartnerCoords = tr;
      }
      else {
        bondPartnerCoords = coords[bondPartner-1];
      }

      int anglePartner   = connTable[i][1];
      vector3d* anglePartnerCoords = 0;

      if (anglePartner < 0) {
        if (anglePartner == -1) anglePartnerCoords = bd;
        if (anglePartner == -2) anglePartnerCoords = ag;
        if (anglePartner == -3) anglePartnerCoords = tr;
      }
      else {
        anglePartnerCoords = coords[anglePartner-1];
      }

      int torsionPartner = connTable[i][2];
      vector3d* torsionPartnerCoords = 0;

      if (torsionPartner < 0) {
        if (torsionPartner == -1) torsionPartnerCoords = bd;
        if (torsionPartner == -2) torsionPartnerCoords = ag;
        if (torsionPartner == -3) torsionPartnerCoords = tr;
      }
      else {
        torsionPartnerCoords = coords[torsionPartner-1];
      }

      double bondLength  = bdAgTors[i][0];
      double angleSize   = bdAgTors[i][1];
      double torsionSize = bdAgTors[i][2];

      //std::cout << bondLength << " " << angleSize << " " << torsionSize << std::endl;

      buildCoord((*pCoords), (*bondPartnerCoords),
                 (*anglePartnerCoords), (*torsionPartnerCoords),
                 bondLength, angleSize*DEG2RAD, torsionSize*DEG2RAD);

      coords[i] = pCoords;
    }

    this->itsCoords.resize(nAtoms);
    for (int i = 0; i < nAtoms; i++) {
      int index = this->getStdAtomIndex(newOrdering[i]);
      this->itsCoords[index] = coords[i];
    }
    return 0;
}

// ============================================================
// Function : getCoordinates()
// ------------------------------------------------------------
//
// ============================================================
std::vector<vector3d*> stdFrag::getCoordinates()
{
    return this->itsCoords;
}

// ============================================================
// Function : print()
// ------------------------------------------------------------
//
// ============================================================
void stdFrag::print()
{
    unsigned int n = this->itsStdAtomList.size();
    for (unsigned int i = 0; i < n; i++) {
      std::cout << i+1 << " " << this->itsStdAtomList[i]->identity << " "
                << this->itsStdAtomList[i]->bond12 << " "
                << this->itsStdAtomList[i]->bond13 << " "
                << this->itsStdAtomList[i]->bond14 << " "
                << this->itsStdAtomList[i]->bondLength << " "
                << this->itsStdAtomList[i]->bondAngle << " "
                << this->itsStdAtomList[i]->bondTorsion << std::endl;
    }

}

//----------------//
// FINGERPRINTS  -//
//----------------//

// ============================================================
// Function : generateSimpleFP()
// ------------------------------------------------------------
// Generate Simple Fingerprint
// ============================================================
void stdFrag::generateSimpleFP()
{
/*
#ifdef DEBUG
      std::cout << " stdFrag::generateSimpleFP " << this->itsName <<  std::endl;
#endif
*/
    unsigned int nRings = this->itsStdRingList.size();

    // Initialize vector
    for (unsigned int i = 0; i < 106; i++) {
      this->itsSimpleFP.push_back(0);
    }

    // Atoms
    int el = 0;
    for (stdAtomIterator c = this->itsStdAtomList.begin();
         c != this->itsStdAtomList.end(); c++) {
      pStdAtom = *c;
      el = pStdAtom->atNum;
      if (el < 54) this->itsSimpleFP[el - 1]++;
    }

    // Bonds
    std::string at1 = "";
    std::string at2 = "";
    std::string at3 = "";
    for (stdBondIterator c = this->itsStdBondList.begin();
         c != this->itsStdBondList.end(); c++) {
      pStdBond = *c;
      if (pStdBond->atom1 < 0) continue;
      if (pStdBond->atom2 < 0) continue;

      pStdAtom1 = this->getStdAtom(pStdBond->atom1);
      pStdAtom2 = this->getStdAtom(pStdBond->atom2);

      at1 = pStdAtom1->atSymbol;
      at2 = pStdAtom2->atSymbol;

      if (at1 == "H" or at2 == "H") {
        if (at1 == "B" or at2 == "B" ) this->itsSimpleFP[53]++;
        if (at1 == "C" or at2 == "C" ) this->itsSimpleFP[54]++;
        if (at1 == "N" or at2 == "N" ) this->itsSimpleFP[55]++;
        if (at1 == "O" or at2 == "O" ) this->itsSimpleFP[56]++;
        if (at1 == "S" or at2 == "S" ) this->itsSimpleFP[57]++;
        continue;
      }

      if (at1 == "B" or at2 == "B") {
        if (at1 == "B") at3 = at2;
        if (at2 == "B") at3 = at1;
        if (at3 == "C") {
          if (pStdBond->type == 1) this->itsSimpleFP[58]++;
          if (pStdBond->type == 2) this->itsSimpleFP[59]++;
        }
        if (at3 == "O")  this->itsSimpleFP[60]++;
        if (at3 == "N")  this->itsSimpleFP[61]++;
        if (at3 == "F")  this->itsSimpleFP[62]++;
        if (at3 == "S")  this->itsSimpleFP[63]++;
        if (at3 == "Cl") this->itsSimpleFP[64]++;
        if (at3 == "Br") this->itsSimpleFP[65]++;
        if (at3 == "I")  this->itsSimpleFP[66]++;
      }

      if (at1 == "C" and at2 == "C") {
        if (pStdBond->type == 1) this->itsSimpleFP[67]++;
        if (pStdBond->type == 2) this->itsSimpleFP[68]++;
        if (pStdBond->type == 3) this->itsSimpleFP[69]++;
      }

      if (at1 == "N" and at2 == "N") {
        if (pStdBond->type == 1) this->itsSimpleFP[70]++;
        if (pStdBond->type == 2) this->itsSimpleFP[71]++;
      }

      if (at1 == "C" or at2 == "C") {
        if (at1 == "C") at3 = at2;
        if (at2 == "C") at3 = at1;

        if (at3 == "N") {
          if (pStdBond->type == 1) this->itsSimpleFP[72]++;
          if (pStdBond->type == 2) this->itsSimpleFP[73]++;
          if (pStdBond->type == 3) this->itsSimpleFP[74]++;
        }
        if (at3 == "F")  this->itsSimpleFP[88]++;
        if (at3 == "Cl") this->itsSimpleFP[89]++;
        if (at3 == "Br") this->itsSimpleFP[90]++;
        if (at3 == "I")  this->itsSimpleFP[91]++;
        if (at3 == "S") {
          if (pStdBond->type == 1) this->itsSimpleFP[92]++;
          if (pStdBond->type == 2) this->itsSimpleFP[93]++;
        }
        if (at3 == "P") this->itsSimpleFP[94]++;
        if (at3 == "Se") {
          if (pStdBond->type == 1) this->itsSimpleFP[95]++;
          if (pStdBond->type == 2) this->itsSimpleFP[96]++;
        }
      }

      if (at1 == "N" or at2 == "N") {
        if (at1 == "N") at3 = at2;
        if (at2 == "N") at3 = at1;
        if (at3 == "O") {
          if (pStdBond->type == 1) this->itsSimpleFP[75]++;
          if (pStdBond->type == 2) this->itsSimpleFP[76]++;
        }
        if (at3 == "P") this->itsSimpleFP[77]++;

        if (at3 == "Se") {
          if (pStdBond->type == 1) this->itsSimpleFP[78]++;
          if (pStdBond->type == 2) this->itsSimpleFP[79]++;
        }
      }

      if (at1 == "O" and at2 == "O") this->itsSimpleFP[80]++;

      if (at1 == "O" or at2 == "O") {
        if (at1 == "O") at3 = at2;
        if (at2 == "O") at3 = at1;
        if (at3 == "C") {
          if (pStdBond->type == 1) this->itsSimpleFP[81]++;
          if (pStdBond->type == 2) this->itsSimpleFP[82]++;
        }

        if (at3 == "Si") this->itsSimpleFP[83]++;

        if (at3 == "S") {
          if (pStdBond->type == 1) this->itsSimpleFP[84]++;
          if (pStdBond->type == 2) this->itsSimpleFP[85]++;
        }

        if (at3 == "Se") {
          if (pStdBond->type == 1) this->itsSimpleFP[86]++;
          if (pStdBond->type == 2) this->itsSimpleFP[87]++;
        }
      }

      if (at1 == "S" and at2 == "S") this->itsSimpleFP[97]++;

      if (at1 == "S" or at2 == "S") {
        if (at1 == "S") at3 = at2;
        if (at2 == "S") at3 = at1;
        if (at3 == "N") this->itsSimpleFP[98]++;
        if (at3 == "P") this->itsSimpleFP[99]++;
      }

      if (at1 == "P" and at2 == "P") this->itsSimpleFP[100]++;

      if (at1 == "P" or at2 == "P") {
        if (at1 == "P") at3 = at2;
        if (at2 == "P") at3 = at1;

        if (at3 == "O") {
          if (pStdBond->type == 1) this->itsSimpleFP[101]++;
          if (pStdBond->type == 2) this->itsSimpleFP[102]++;
        }
        if (at3 == "Se") this->itsSimpleFP[103]++;
      }

      if (at1 == "Se" and at2 == "Se") this->itsSimpleFP[104]++;
    } // bonds

    // Loops
    at1 = "";
    at2 = "";
    at3 = "";
    for (stdLoopIterator c = this->itsStdLoopList.begin();
         c != this->itsStdLoopList.end(); c++) {
      pStdLoop = *c;
      pStdAtom1 = this->getStdAtom(pStdLoop->atom1);
      pStdAtom2 = this->getStdAtom(pStdLoop->atom2);

      if (!(pStdAtom1 and pStdAtom2)) return;

      at1 = pStdAtom1->atSymbol;
      at2 = pStdAtom2->atSymbol;

      if (at1 == "H" or at2 == "H") {
        if (at1 == "B" or at2 == "B" ) this->itsSimpleFP[53]++;
        if (at1 == "C" or at2 == "C" ) this->itsSimpleFP[54]++;
        if (at1 == "N" or at2 == "N" ) this->itsSimpleFP[55]++;
        if (at1 == "O" or at2 == "O" ) this->itsSimpleFP[56]++;
        if (at1 == "S" or at2 == "S" ) this->itsSimpleFP[57]++;
        continue;
      }

      if (at1 == "B" or at2 == "B") {
        if (at1 == "B") at3 = at2;
        if (at2 == "B") at3 = at1;
        if (at3 == "C") {
          if (pStdLoop->type == 1) this->itsSimpleFP[58]++;
          if (pStdLoop->type == 2) this->itsSimpleFP[59]++;
        }
        if (at3 == "O")  this->itsSimpleFP[60]++;
        if (at3 == "N")  this->itsSimpleFP[61]++;
        if (at3 == "F")  this->itsSimpleFP[62]++;
        if (at3 == "S")  this->itsSimpleFP[63]++;
        if (at3 == "Cl") this->itsSimpleFP[64]++;
        if (at3 == "Br") this->itsSimpleFP[65]++;
        if (at3 == "I")  this->itsSimpleFP[66]++;
      }

      if (at1 == "C" and at2 == "C") {
        if (pStdLoop->type == 1) this->itsSimpleFP[67]++;
        if (pStdLoop->type == 2) this->itsSimpleFP[68]++;
        if (pStdLoop->type == 3) this->itsSimpleFP[69]++;
      }

      if (at1 == "N" and at2 == "N") {
        if (pStdLoop->type == 1) this->itsSimpleFP[70]++;
        if (pStdLoop->type == 2) this->itsSimpleFP[71]++;
      }

      if (at1 == "C" or at2 == "C") {
        if (at1 == "C") at3 = at2;
        if (at2 == "C") at3 = at1;

        if (at3 == "N") {
          if (pStdLoop->type == 1) this->itsSimpleFP[72]++;
          if (pStdLoop->type == 2) this->itsSimpleFP[73]++;
          if (pStdLoop->type == 3) this->itsSimpleFP[74]++;
        }
        if (at3 == "F")  this->itsSimpleFP[88]++;
        if (at3 == "Cl") this->itsSimpleFP[89]++;
        if (at3 == "Br") this->itsSimpleFP[90]++;
        if (at3 == "I")  this->itsSimpleFP[91]++;
        if (at3 == "S") {
          if (pStdLoop->type == 1) this->itsSimpleFP[92]++;
          if (pStdLoop->type == 2) this->itsSimpleFP[93]++;
        }
        if (at3 == "P") this->itsSimpleFP[94]++;
        if (at3 == "Se") {
          if (pStdLoop->type == 1) this->itsSimpleFP[95]++;
          if (pStdLoop->type == 2) this->itsSimpleFP[96]++;
        }
      }

      if (at1 == "N" or at2 == "N") {
        if (at1 == "O") at3 = at2;
        if (at2 == "O") at3 = at1;
        if (at3 == "O") {
          if (pStdLoop->type == 1) this->itsSimpleFP[75]++;
          if (pStdLoop->type == 2) this->itsSimpleFP[76]++;
        }
        if (at1 == "P" or at2 == "P") this->itsSimpleFP[77]++;

        if (at1 == "Se") at3 = at2;
        if (at2 == "Se") at3 = at1;
        if (at3 == "Se") {
          if (pStdLoop->type == 1) this->itsSimpleFP[78]++;
          if (pStdLoop->type == 2) this->itsSimpleFP[79]++;
        }
      }

      if (at1 == "O" and at2 == "O") this->itsSimpleFP[80]++;

      if (at1 == "O" or at2 == "O") {
        if (at1 == "O") at3 = at2;
        if (at2 == "O") at3 = at1;
        if (at3 == "C") {
          if (pStdLoop->type == 1) this->itsSimpleFP[81]++;
          if (pStdLoop->type == 2) this->itsSimpleFP[82]++;
        }

        if (at1 == "Si" or at2 == "Si") this->itsSimpleFP[83]++;

        if (at3 == "S") {
          if (pStdLoop->type == 1) this->itsSimpleFP[84]++;
          if (pStdLoop->type == 2) this->itsSimpleFP[85]++;
        }

        if (at3 == "Se") {
          if (pStdLoop->type == 1) this->itsSimpleFP[86]++;
          if (pStdLoop->type == 2) this->itsSimpleFP[87]++;
        }
      }

      if (at1 == "S" and at2 == "S") this->itsSimpleFP[97]++;

      if (at1 == "S" or at2 == "S") {
        if (at1 == "S") at3 = at2;
        if (at2 == "S") at3 = at1;
        if (at3 == "N") this->itsSimpleFP[98]++;
        if (at3 == "P") this->itsSimpleFP[99]++;
      }

      if (at1 == "P" and at2 == "P") this->itsSimpleFP[100]++;

      if (at1 == "P" or at2 == "P") {
        if (at1 == "P") at3 = at2;
        if (at2 == "P") at3 = at1;

        if (at3 == "O") {
          if (pStdLoop->type == 1) this->itsSimpleFP[101]++;
          if (pStdLoop->type == 2) this->itsSimpleFP[102]++;
        }
        if (at3 == "Se") this->itsSimpleFP[103]++;
      }

      if (at1 == "Se" and at2 == "Se") this->itsSimpleFP[104]++;
    }

    // Rings
    if (nRings > 0) {
      this->itsSimpleFP[105] = nRings;

      for (stdRingIterator c = this->itsStdRingList.begin();
           c != this->itsStdRingList.end(); c++) {
        pStdRing = *c;
        this->itsSimpleFP.push_back(pStdRing->size);
        this->itsSimpleFP.push_back(pStdRing->planar);
        this->itsSimpleFP.push_back(pStdRing->aromatic);
        this->itsSimpleFP.push_back(pStdRing->hetero);
        this->itsSimpleFP.push_back(pStdRing->nNitrogen);
        this->itsSimpleFP.push_back(pStdRing->nOxygen);
        this->itsSimpleFP.push_back(pStdRing->nSulfur);
      }
    }
/*
#ifdef DEBUG
    for (unsigned int i = 0; i < this->itsSimpleFP.size(); i++) {
      std::cout << this->itsSimpleFP[i] << " ";
    }
    std::cout << " " << std::endl;
#endif
*/
}

// ============================================================
// Function : getSimpleFP()
// ------------------------------------------------------------
// Get Simple Fingerprint
// ============================================================
std::vector<unsigned int> stdFrag::getSimpleFP()
{
    return this->itsSimpleFP;
}

// ============================================================
// Function : generateAdjMatrix()
// ------------------------------------------------------------
// Generate Adjacency Matrix
// ============================================================
int stdFrag::generateAdjMatrix()
{
    int nAtoms = this->numStdAtoms();
    this->adjMatrixSize = nAtoms*nAtoms;

    try {
      this->adjMatrix   = new int [this->adjMatrixSize];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }

    for (int i = 0; i < nAtoms; i++) {
      for (int j = 0; j < nAtoms; j++) {
        pStdBond = this->getStdBond(i+1, j+1);
        pStdLoop = this->getStdLoop(i+1, j+1);
        if (pStdBond) {
          this->adjMatrix[i*nAtoms+j] = pStdBond->type;
        }
        else if (pStdLoop) {
          this->adjMatrix[i*nAtoms+j] = pStdLoop->type;
        }
        else {
          this->adjMatrix[i*nAtoms+j] = 0;
        }
      }
    }
/*
#ifdef DEBUG
    std::cout << " stdFrag::generateAdjMatrix : "
              << " fragAdjMatrix = " << std::endl;
    for (int i = 0; i < nAtoms; i++) {
      for (int j = 0; j < nAtoms; j++) {
        std::cout << this->adjMatrix[i*nAtoms+j] << " ";
      }
      std::cout << " " << std::endl;
    }
    std::cout << " " << std::endl;
#endif
*/
    int r = this->generateAtomSymbols();
    return r;
}

// ============================================================
// Function : generateHeavyAdjMatrix()
// ------------------------------------------------------------
// Generate Heavy Atom Adjacency Matrix
// ============================================================
int stdFrag::generateHeavyAdjMatrix()
{
    int nAtoms = this->numStdAtoms();

    int nHeavyAtoms = this->numStdHeavyAtoms();
    this->heavyAtomAdjMatrixSize = nHeavyAtoms*nHeavyAtoms;

    try {
      this->heavyAtomAdjMatrix   = new int [this->heavyAtomAdjMatrixSize];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }

    int iIndex = -1;
    int jIndex = -1;
    for (int i = 0; i < nAtoms; i++) {
      if (this->itsStdAtomList[i]->atSymbol == "H") {
        continue;
      }
      else {
        iIndex++;
      }
      for (int j = 0; j < nAtoms; j++) {
        if (this->itsStdAtomList[j]->atSymbol == "H") {
          continue;
        }
        else {
          jIndex++;
        }
        pStdBond = this->getStdBond(i+1, j+1);
        pStdLoop = this->getStdLoop(i+1, j+1);
        if (pStdBond) {
          this->heavyAtomAdjMatrix[iIndex*nHeavyAtoms+jIndex] = pStdBond->type;
        }
        else if (pStdLoop) {
          this->heavyAtomAdjMatrix[iIndex*nHeavyAtoms+jIndex] = pStdLoop->type;
        }
        else {
          this->heavyAtomAdjMatrix[iIndex*nHeavyAtoms+jIndex] = 0;
        }
      }
      jIndex = -1;
    }
/*
#ifdef DEBUG
    std::cout << " stdFrag::generateAdjMatrix : "
              << " fragAdjMatrix = " << std::endl;
    for (int i = 0; i < nHeavyAtoms; i++) {
      for (int j = 0; j < nHeavyAtoms; j++) {
        std::cout << this->adjMatrix[i*nHeavyAtoms+j] << " ";
      }
      std::cout << " " << std::endl;
    }
    std::cout << " " << std::endl;
#endif
*/
    return 0;
}

// ============================================================
// Function : getAdjMatrix()
// ------------------------------------------------------------
// Get Adjacency Matrix Size
// ============================================================
int* stdFrag::getAdjMatrix()
{
    return this->adjMatrix;
}

// ============================================================
// Function : getHeavyAdjMatrix()
// ------------------------------------------------------------
// Get Heavy Atom Adjacency Matrix Size
// ============================================================
int* stdFrag::getHeavyAdjMatrix()
{
    return this->heavyAtomAdjMatrix;
}

// ============================================================
// Function : getAdjMatrixSize()
// ------------------------------------------------------------
// Get Adjacency Matrix Size
// ============================================================
int stdFrag::getAdjMatrixSize()
{
    return this->adjMatrixSize;
}

// ============================================================
// Function : getHeavyAdjMatrixSize()
// ------------------------------------------------------------
// Get Heavy Atom Adjacency Matrix Size
// ============================================================
int stdFrag::getHeavyAdjMatrixSize()
{
    return this->heavyAtomAdjMatrixSize;
}

// ============================================================
// Function : generateAtomSymbols()
// ------------------------------------------------------------
// Generate atom symbols array
// ============================================================
int stdFrag::generateAtomSymbols()
{
    int nAtoms = this->numStdAtoms();

    try {
      this->atomSymbols   = new char [nAtoms*2];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }

    int index = 0;
    for (stdAtomIterator c = this->itsStdAtomList.begin();
                         c != this->itsStdAtomList.end(); c++) {
      pStdAtom = *c;
      this->atomSymbols[index  ] = pStdAtom->atSymbol[0];
      this->atomSymbols[index+1] = pStdAtom->atSymbol[1];
      index+=2;
    }
    return 0;
}

// ============================================================
// Function : generateHeavyAtomSymbols()
// ------------------------------------------------------------
// Generate heavy atom symbols array
// ============================================================
int stdFrag::generateHeavyAtomSymbols()
{
    int nAtoms = this->numStdHeavyAtoms();

    try {
      this->heavyAtomSymbols   = new char [nAtoms*2];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }

    int index = 0;
    for (stdAtomIterator c = this->itsStdAtomList.begin();
                         c != this->itsStdAtomList.end(); c++) {
      pStdAtom = *c;
      if (pStdAtom->atSymbol != "H") {
        this->heavyAtomSymbols[index  ] = pStdAtom->atSymbol[0];
        this->heavyAtomSymbols[index+1] = pStdAtom->atSymbol[1];
        index+=2;
      }
    }
    return 0;
}

// ============================================================
// Function : getAtomSymbols()
// ------------------------------------------------------------
// Get atom symbols array
// ============================================================
char* stdFrag::getAtomSymbols()
{
    return this->atomSymbols;
}

// ============================================================
// Function : getHeavyAtomSymbols()
// ------------------------------------------------------------
// Get heavy atom symbols array
// ============================================================
char* stdFrag::getHeavyAtomSymbols()
{
    int f = this->generateHeavyAtomSymbols();
    if (f) {
      std::cout << "   Error in stdFrag::getHeavyAtomSymbols ... " << std::endl;
      //exit(0);
      throw MTKException( "   Error in stdFrag::getHeavyAtomSymbols ... ");
    }
    return this->heavyAtomSymbols;
}

// ============================================================
// Function : generateAtomKinds()
// ------------------------------------------------------------
// Generate atom kinds array
// ============================================================
int stdFrag::generateAtomKinds()
{
    int nAtoms = this->numStdAtoms();

    try {
      this->atomKinds   = new int [nAtoms];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }
    int index = 0;
    for (stdAtomIterator c = this->itsStdAtomList.begin();
         c != this->itsStdAtomList.end(); c++) {
      pStdAtom = *c;
      this->atomKinds[index] = pStdAtom->kind;
      index++;
    }
    return 0;
}

// ============================================================
// Function : generateHeavyAtomKinds()
// ------------------------------------------------------------
// Generate heavy atom kinds array
// ============================================================
int stdFrag::generateHeavyAtomKinds()
{
    int nAtoms = this->numStdHeavyAtoms();

    try {
      this->heavyAtomKinds   = new int [nAtoms];
    }
    catch (std::bad_alloc) {
      std::cout << " Memory Allocation Failure " << std::endl;
      return 1;
    }
    int index = 0;
    for (stdAtomIterator c = this->itsStdAtomList.begin();
         c != this->itsStdAtomList.end(); c++) {
      pStdAtom = *c;
      if (pStdAtom->atSymbol != "H") {
        this->heavyAtomKinds[index] = pStdAtom->kind;
        index++;
      }
    }
    return 0;
}

// ============================================================
// Function : getAtomKinds()
// ------------------------------------------------------------
// Get atom kinds array
// ============================================================
int* stdFrag::getAtomKinds()
{
    return this->atomKinds;
}

// ============================================================
// Function : getHeavyAtomKinds()
// ------------------------------------------------------------
// Get heavy atom kinds array
// ============================================================
int* stdFrag::getHeavyAtomKinds()
{
    int f = generateHeavyAtomKinds();
    if (f) {
      std::cout << "   Error in stdFrag::generateHeavyAtomKinds ... exiting \n";
      throw MTKException("   Error in stdFrag::generateHeavyAtomKinds ... exiting \n");
      //exit(0);
    }
    return this->heavyAtomKinds;
}

// ============================================================
// Function : getCharge()
// ------------------------------------------------------------
// Get fragment charge
// ============================================================
double stdFrag::getCharge()
{
    double charge = 0.0;
    for (stdAtomIterator c = this->itsStdAtomList.begin();
         c != this->itsStdAtomList.end(); c++) {
      pStdAtom = *c;
      charge += pStdAtom->atmCharge;
    }
    return charge;
}

// ============================================================
// Function : getStdAtomIndex()
// ------------------------------------------------------------
//
// ============================================================
int stdFrag::getStdAtomIndex(stdAtom* a)
{
    int i = 0;
    for (stdAtomIterator c = this->itsStdAtomList.begin();
         c != this->itsStdAtomList.end(); c++) {
      if (a == *c) return i;
      i++;
    }
    return -1;
}

} // MTKpp namespace


