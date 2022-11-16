/*!
   \file selection.cpp
   \brief Allows for the selection of certain parts of the collection
   \author Martin Peters

   $Date: 2010/08/19 11:33:30 $
   $Revision: 1.13 $

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

#include "collection.h"
#include "selection.h"
#include "molecule.h"
#include "submolecule.h"
#include "atom.h"
#include "element.h"
#include "utility.h"
#include "Log/errorHandler.h"

namespace MTKpp
{

// ============================================================
// Function : selection()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
selection::selection(collection *parent):pParent(parent)
{
    selectionType = 0;
    itsName  = "";
    selnMol  = 0;
    selnSMol = 0;
    selnAtom = 0;
}

// ============================================================
// Function : ~selection()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
selection::~selection() {}

// ============================================================
// Function : getParent()
// ------------------------------------------------------------
//
// ============================================================
collection* selection::getParent()
{
    return pParent;
}

// ============================================================
// Function : parse()
// ------------------------------------------------------------
//
// ============================================================
int selection::parse(std::string seln)
{
    /*
      Hierarchy
      /col/mol/submol/atom
       ^    ^    ^     ^
       |    |    |     |++ atom name, number, or name-number
       |    |    |++ submolecule name, number, or name-number
       |    |++ molecule name, number, or name-number
       |++ collection name

       Selection possibilities:
       1) /col/mol

       2) /col/mol/smol
       3) /col/   /smol

       4) /col/mol/smol/atom
       5) /col/mol/    /atom
       6) /col/   /smol/atom
       7) /col/   /    /atom
    */

    std::string errMessage = " selection string = " + seln;
    errorLogger.throwError("selection::parse", errMessage, 4);

////
/// remove trailing /
////

    while (seln.find("//", 0)) {
      if (seln.find("//", 0) > seln.size()) break;
      std::string::size_type pos = seln.find("//", 0);
      seln.insert(pos+1, " ");
    }

    std::vector<std::string> selnSplit;
    this->splitString(seln, "/", selnSplit, 0);

    std::vector<int> types;
    int nAlpha = 0;
    int nNum = 0;
    int nPunct = 0;

    for (unsigned int i = 0; i < selnSplit.size(); i++) {
      nAlpha = 0;
      nNum = 0;
      nPunct = 0;
      for (unsigned int j = 0; j < selnSplit[i].size(); j++) {
        if (selnSplit[i][j] == '.') selnSplit[i][j] = ' ';
        if (isalpha(selnSplit[i][j])) {
          nAlpha++;
        }
        if (isdigit(selnSplit[i][j])) {
          nNum++;
        }
        if (ispunct(selnSplit[i][j])) {
          nPunct++;
        }
      }
      if (nNum == int(selnSplit[i].size())) {
        types.push_back(2);
      }
      else if (nPunct > 0) {
        types.push_back(3);
      }
      else {
        types.push_back(1);
      }
    }

    //std::vector<atom*> selnAtomList;

    const char sep = '/';
    // Check to see if first character is a '/'
    if (seln[0] == sep) {
      int selnSplitSize = selnSplit.size();
      if (selnSplitSize > 4) {
        MTKpp::errorLogger.throwError("selection::parse", " selection has more than 5 terms ... exiting ", 1);
        return 1;
      }
      // molecule
      else if (selnSplitSize == 2) { // 1. /col/mol
        if (selnSplit[1] != " ") {
          if (types[1] == 1) {
            this->selnMol = pParent->getMolecule(selnSplit[1]);
            if (this->selnMol) {
              this->itsAtoms = this->selnMol->getAtomList();
              selectionType = 1;
            }
            else {
              MTKpp::errorLogger.throwError("selection::parse", " /col/mol selection error (1) ... exiting ", 1);
              return 1;
            }
          }
          else if (types[1] > 1) {
            int molnum = 0;
            std::string molname = "";
            if (types[1] == 3) {
              std::vector<std::string> nameNum;
              this->splitString(selnSplit[1], "-", nameNum, 0);
              if (nameNum.size() == 1) {
                this->selnMol = pParent->getMolecule(selnSplit[1]);
                if (this->selnMol) {
                  this->itsAtoms = this->selnMol->getAtomList();
                  selectionType = 1;
                }
                else {
                  MTKpp::errorLogger.throwError("selection::parse", " /col/mol selection error (1) ... exiting ", 1);
                  return 1;
                }
                return 0;
              }
              molname = nameNum[0];
              molnum = atoi(nameNum[1].c_str());
            }
            else {
              molnum = atoi(selnSplit[1].c_str());
            }
            this->selnMol = pParent->getMolecule(molnum);

            if (this->selnMol) {
              if ((types[1] == 3) and (this->selnMol->getName() != molname)) {
                std::string strTemp = " User supplied " + selnSplit[1] + ", however molecule " + i2s(molnum)
                          + " has a name of " + this->selnMol->getName();
                MTKpp::errorLogger.throwError("selection::parse", strTemp, 1);

                return 1;
              }
              this->itsAtoms = this->selnMol->getAtomList();
              selectionType = 1;
            }
            else {
              MTKpp::errorLogger.throwError("selection::parse", " /col/mol selection error (2) ... exiting ", 1);
              return 1;
            }
          }
        }
      }
      // submolecule
      else if (selnSplitSize == 3) { // /col/mol/smol
        selectionType = 2;
        if (selnSplit[1] != " ") { //  /col/mol
          int molnum = 0;
          std::string molname = "";
          if (types[1] == 1) {
            this->selnMol = pParent->getMolecule(selnSplit[1]);
          }
          else if (types[1] > 1) {
            if (types[1] == 3) {
              std::vector<std::string> nameNum;
              this->splitString(selnSplit[1], "-", nameNum, 0);
              molname = nameNum[0];
              molnum = atoi(nameNum[1].c_str());
            }
            else {
              molnum = atoi(selnSplit[1].c_str());
            }
            this->selnMol = pParent->getMolecule(molnum);
          }

          if (this->selnMol) {
            if ((types[1] == 3) and (this->selnMol->getName() != molname)) {
              std::string strTemp = " User supplied " + selnSplit[1] + ", however molecule " + i2s(molnum)
                          + " has a name of " + this->selnMol->getName();
              MTKpp::errorLogger.throwError("selection::parse", strTemp, 1);
              return 1;
            }
          }
          else {
            // add error message
            MTKpp::errorLogger.throwError("selection::parse", " /col/mol/smol selection error (3) ... exiting ", 1);
            return 1;
          }

          if (types[2] == 1) {
            std::vector<submolecule*> sMolList = this->selnMol->getSubMoleculeList(selnSplit[2]);
            for (unsigned int s = 0; s < sMolList.size(); s++) { // 2. /col/mol/smol
              this->selnSMol = sMolList[s];
              std::vector<atom*> selnAtomList = sMolList[s]->getAtomList();
              for (unsigned int a = 0; a < selnAtomList.size(); a++) {
                this->itsAtoms.push_back(selnAtomList[a]);
                this->selnAtom = selnAtomList[a];
              }
            }
          }
          else if (types[2] > 1) {
            int resnum = 0;
            std::string resname = "";
            if (types[2] == 3) {
              std::vector<std::string> nameNum;
              this->splitString(selnSplit[2], "-", nameNum, 0);
              resname = nameNum[0];
              resnum = atoi(nameNum[1].c_str());
            }
            else {
              resnum = atoi(selnSplit[2].c_str());
            }
            this->selnSMol = this->selnMol->getSubMolecule(resnum);

            if (this->selnSMol) {
              if ((types[2] == 3) and (this->selnSMol->getName() != resname)) {
                std::string strTemp = " User supplied " + selnSplit[2] + ", however residue " + i2s(resnum)
                          + " has a name of " + this->selnSMol->getName();
                MTKpp::errorLogger.throwError("selection::parse", strTemp, 1);
                return 1;
              }
            }
            else {
              std::string strTemp = " Can't find: " + resname + "-" + i2s(resnum) + " in molecule: "
                                  + this->selnMol->getName() + " /col/mol/smol selection error (4)";

              MTKpp::errorLogger.throwError("selection::parse", strTemp, 1);
              return 1;
            }

            std::vector<atom*> selnAtomList = this->selnSMol->getAtomList();
            for (unsigned int a = 0; a < selnAtomList.size(); a++) {
              this->itsAtoms.push_back(selnAtomList[a]);
              this->selnAtom = selnAtomList[a];
            }
          }
        }
        else {

          //std::cout << "\n\n Selection: " << seln << std::endl;
          //std::cout << " type = " << types[2] << std::endl;

          std::vector<molecule*> molList = pParent->getMoleculeList();
          for (unsigned int m = 0; m < molList.size(); m++) { // 3. /col/ /smol
            //this->selnMol = molList[m];
            std::string smolname;
            int smolnum;
            std::vector<submolecule*> sMolList;

            if (types[2] == 1) {
              sMolList = molList[m]->getSubMoleculeList(selnSplit[2]);
            }
            else {
              if (types[2] == 3) {
                std::vector<std::string> smolNameNum;
                this->splitString(selnSplit[2], "-", smolNameNum, 0);
                smolname = smolNameNum[0];
                smolnum = atoi(smolNameNum[1].c_str());
                sMolList = molList[m]->getSubMoleculeList(smolname, smolnum);
                //std::cout << smolname << " " << smolnum << std::endl;
                //std::cout << molList[m]->getName() << std::endl;
              }
              else {
                smolnum = atoi(selnSplit[2].c_str());
                sMolList = molList[m]->getSubMoleculeList(smolnum);
              }
            }
            for (unsigned int s = 0; s < sMolList.size(); s++) {
              this->selnSMol = sMolList[s];
              this->selnMol = this->selnSMol->getParent();
              std::vector<atom*> selnAtomList = sMolList[s]->getAtomList();
              for (unsigned int a = 0; a < selnAtomList.size(); a++) {
                this->itsAtoms.push_back(selnAtomList[a]);
              }
            }
          }

          //std::cout << "result:" << std::endl;
          //std::cout << this->selnMol->getName() << std::endl;
          //std::cout << " " << std::endl;
        }
      }
      // atom
      else if (selnSplitSize == 4) {
        int atomnum = 0;
        std::string atomname = "";
        if (types[3] == 1) { // atom name
            atomname = selnSplit[3];
        }
        else if (types[3] > 1) {
          if (types[3] == 3) {
            std::vector<std::string> nameNum;
            this->splitString(selnSplit[3], "-", nameNum, 0);
            atomname = nameNum[0];
            atomnum = atoi(nameNum[1].c_str());
          }
          else {
            atomnum = atoi(selnSplit[3].c_str());
          }
        }
        if (selnSplit[1] != " ") { //  /col/mol
          int molnum = 0;
          std::string molname = "";
          if (types[1] == 1) { // mol name
            this->selnMol = pParent->getMolecule(selnSplit[1]);
          }
          else if (types[1] > 1) {
            if (types[1] == 3) {
              std::vector<std::string> nameNum;
              this->splitString(selnSplit[1], "-", nameNum, 0);
              molname = nameNum[0];
              molnum = atoi(nameNum[1].c_str());
            }
            else {
              molnum = atoi(selnSplit[1].c_str());
            }
            this->selnMol = pParent->getMolecule(molnum);
          }
          if (this->selnMol) {
            if ((types[1] == 3) and (this->selnMol->getName() != molname)) {
              std::string strTemp = " User supplied " + selnSplit[1] + ", however molecule " + i2s(molnum)
                          + " has a name of " + this->selnMol->getName();
              MTKpp::errorLogger.throwError("selection::parse", strTemp, 1);
              return 1;
            }
          }
          else {
            MTKpp::errorLogger.throwError("selection::parse", " /col/mol selection error (5) ... exiting ", 1);
            return 1;
          }

          if (selnSplit[2] != " ") { //  /col/mol/smol/atom
            int resnum = 0;
            std::string resname = "";
            if (types[2] == 1) {
              this->selnSMol = this->selnMol->getSubMolecule(selnSplit[2]);
            }
            else if (types[2] > 1) {
              if (types[2] == 3) {
                std::vector<std::string> nameNum;
                this->splitString(selnSplit[2], "-", nameNum, 0);
                resname = nameNum[0];
                resnum = atoi(nameNum[1].c_str());
              }
              else {
                resnum = atoi(selnSplit[2].c_str());
              }
              this->selnSMol = this->selnMol->getSubMolecule(resnum);
            }
            if (this->selnSMol) {
              if ((types[2] == 3) and (this->selnSMol->getName() != resname)) {
                std::string strTemp = " User supplied " + selnSplit[1] + ", however residue " + i2s(resnum)
                            + " has a name of " + this->selnSMol->getName();
                MTKpp::errorLogger.throwError("selection::parse", strTemp, 1);
                return 1;
              }
            }
            else {
              MTKpp::errorLogger.throwError("selection::parse", " /col/mol/smol/atom selection error (6) ... exiting ", 1);
              return 1;
            }
            if (types[3] == 1) {
              this->selnAtom = this->selnSMol->getAtom(atomname); // 4. /col/mol/smol/atom
            }
            else if (types[3] == 2) {
              this->selnAtom = this->selnSMol->getAtom(atomnum); // 4. /col/mol/smol/atom
            }
            else {
              this->selnAtom = this->selnSMol->getAtom(atomnum); // 4. /col/mol/smol/atom
              if (!this->selnAtom) return 1;
              if (this->selnAtom->getName() != atomname) {
                std::string strTemp = " User supplied " + selnSplit[3] + ", however atom " + i2s(atomnum)
                            + " has a name of " + this->selnAtom->getName();
                MTKpp::errorLogger.throwError("selection::parse", strTemp, 1);
                return 1;
              }
            }
            if (!this->selnAtom) {
              MTKpp::errorLogger.throwError("selection::parse", " /col/mol/smol/atom selection error (7) ... exiting ", 1);
              return 1;
            }
            itsAtoms.push_back(this->selnAtom);
            selectionType = 3;
          }
          else { // 5. /col/mol/ /atom
            std::vector<submolecule*> smolList = this->selnMol->getSubMoleculeList();
            for (unsigned int s = 0; s < smolList.size(); s++) {
              this->selnSMol = smolList[s];
              if (types[3] == 1) {
                this->selnAtom = smolList[s]->getAtom(atomname);
              }
              else if (types[3] == 2) {
                this->selnAtom = smolList[s]->getAtom(atomnum);
              }
              else {
                this->selnAtom = smolList[s]->getAtom(atomnum);
                if (!this->selnAtom) return 1;
                if (this->selnAtom->getName() != atomname) {
                  std::string strTemp = " User supplied " + selnSplit[3] + ", however atom " + i2s(atomnum)
                              + " has a name of " + this->selnAtom->getName();
                  MTKpp::errorLogger.throwError("selection::parse", strTemp, 1);
                  return 1;
                }
              }
              if (!this->selnAtom) continue;
              itsAtoms.push_back(this->selnAtom);
              selectionType = 3;
            }
          }
        }
        else { //  /col/ /
          std::vector<molecule*> molList = pParent->getMoleculeList();
          for (unsigned int m = 0; m < molList.size(); m++) {
            this->selnMol = molList[m];
            if (selnSplit[2] != " ") { // 6. /col/ /smol/atom
              int resnum = 0;
              std::string resname = "";
              if (types[2] == 1) {
                this->selnSMol = molList[m]->getSubMolecule(selnSplit[2]);
              }
              else if (types[2] > 1) {
                if (types[2] == 3) {
                  std::vector<std::string> nameNum;
                  this->splitString(selnSplit[2], "-", nameNum, 0);
                  resname = nameNum[0];
                  resnum = atoi(nameNum[1].c_str());
                }
                else {
                  resnum = atoi(selnSplit[2].c_str());
                }
                this->selnSMol = molList[m]->getSubMolecule(resnum);
              }
              if (this->selnSMol) {
                if ((types[2] == 3) and (this->selnSMol->getName() != resname)) {
                  std::string strTemp = " User supplied " + selnSplit[2] + ", however residue " + i2s(resnum)
                              + " has a name of " + this->selnSMol->getName();
                  MTKpp::errorLogger.throwError("selection::parse", strTemp, 1);
                  return 1;
                }
                if (types[3] == 1) {
                  this->selnAtom = this->selnSMol->getAtom(atomname);
                }
                else if (types[3] > 1) {
                  this->selnAtom  = this->selnSMol->getAtom(atomnum);
                }
                if (this->selnAtom) {
                  if ((types[3] == 3) and (this->selnAtom->getName() != atomname)) {
                    std::string strTemp = " User supplied " + selnSplit[3] + ", however atom " + i2s(atomnum)
                                + " has a name of " + this->selnAtom->getName();
                    MTKpp::errorLogger.throwError("selection::parse", strTemp, 1);
                    return 1;
                  }
                  itsAtoms.push_back(this->selnAtom);
                  selectionType = 3;
                }
              }
            }
            else { // 7. /col/ / /atom
              std::vector<submolecule*> smolList = molList[m]->getSubMoleculeList();
              for (unsigned int s = 0; s < smolList.size(); s++) {
                if (types[3] == 1) {
                  this->selnAtom = smolList[s]->getAtom(atomname);
                }
                else if (types[3] > 1) {
                  this->selnAtom = smolList[s]->getAtom(atomnum);
                }
                if (this->selnAtom) {
                  if ((types[3] == 3) and (this->selnAtom->getName() != atomname)) {
                    std::string strTemp = " User supplied " + selnSplit[3] + ", however atom " + i2s(atomnum)
                                + " has a name of " + this->selnAtom->getName();
                    MTKpp::errorLogger.throwError("selection::parse", strTemp, 1);
                    return 1;
                  }
                }
                else {
                  MTKpp::errorLogger.throwError("selection::parse", " /col/ / /atom selection error (10) ... exiting ", 1);
                  return 1;
                }
                itsAtoms.push_back(this->selnAtom);
                selectionType = 3;
              }
            }
          }
        }
      }
    }
    else {
      MTKpp::errorLogger.throwError("selection::parse", " Need to implement the filling of the hierarchy from the bottom ... exiting ", 1);
      return 1;
    }
/*
#ifdef DEBUG
    std::string errorMessage = " Selection Size: " + itsAtoms.size();
    MTKpp::errorLogger.throwError("selection::parse", errorMessage, 4);

    std::stringstream number;
    for (unsigned int x = 0; x < itsAtoms.size(); x++) {
      number << x+1;
      std::string errorMessage = "  seln " + number.str() + ": " +
           itsAtoms[x]->getFileID() + " " + itsAtoms[x]->getName();
      MTKpp::errorLogger.throwError("selection::parse", errorMessage, 4);

//      MTKpp::errorLogger.getStream() << "  seln " << x+1 << ": "
//              << itsAtoms[x]->getFileID()
//              << " " << itsAtoms[x]->getName() << std::endl;
    }
#endif
*/
    return 0;
}

//=============================================================
// Function : getSelectionType
// ------------------------------------------------------------
//
// ============================================================
int selection::getSelectionType()
{
    return this->selectionType;
}

//=============================================================
// Function : getAtoms
// ------------------------------------------------------------
//
// ============================================================
std::vector<atom*> selection::getAtoms()
{
    return this->itsAtoms;
}

//=============================================================
// Function : setMol
// ------------------------------------------------------------
//
// ============================================================
void selection::setMol(molecule* m)
{
    this->selnMol = m;
}

//=============================================================
// Function : getMol
// ------------------------------------------------------------
//
// ============================================================
molecule* selection::getMol()
{
    return this->selnMol;
}

//=============================================================
// Function : setSMol
// ------------------------------------------------------------
//
// ============================================================
void selection::setSMol(submolecule* s)
{
    this->selnSMol = s;
}

//=============================================================
// Function : getSMol
// ------------------------------------------------------------
//
// ============================================================
submolecule* selection::getSMol()
{
    return this->selnSMol;
}

//=============================================================
// Function : setAtom
// ------------------------------------------------------------
//
// ============================================================
void selection::setAtom(atom* a)
{
    this->selnAtom = a;
}

//=============================================================
// Function : getAtom
// ------------------------------------------------------------
//
// ============================================================
atom* selection::getAtom()
{
    return this->selnAtom;
}

//=============================================================
// Function : splitString
// ------------------------------------------------------------
// splits up a string based on a separator and returns a vector
// ============================================================
void selection::splitString(std::string &text, const std::string separator,
                 std::vector<std::string> &words, int mystart)
{
    int n = text.length();
    int start, stop;

    start = text.find_first_not_of(separator, mystart);

    while ( (start >=0) && (start < n)) {
      stop = text.find_first_of(separator, start);
      if ((stop < 0) || (stop > n)) {
        stop = n;
      }

      words.push_back(text.substr(start, stop-start));
      start = text.find_first_not_of(separator, stop+1);
    }
}

} // MTKpp namespace
