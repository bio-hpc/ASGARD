/*!
   \file prepParser.cpp
   \brief Parses AMBER prep files
   \author Martin Peters

   $Date: 2010/03/29 20:39:35 $
   $Revision: 1.13 $

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
#include "prepParser.h"

#include "StringManip.h"

#include "Molecule/collection.h"
#include "Molecule/molecule.h"
#include "Molecule/submolecule.h"
#include "Molecule/atom.h"
#include "Molecule/bond.h"
#include "Molecule/element.h"
#include "Molecule/stdLibrary.h"
#include "Molecule/stdGroup.h"
#include "Molecule/stdFrag.h"

#include "Log/errorHandler.h"

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

// ============================================================
// Function : prepParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
prepParser::prepParser():baseParser() {}

// =========================================================
// Function : prepParser()
// ---------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// =========================================================
prepParser::~prepParser() {}

// =========================================================
// Function : openFile
// ---------------------------------------------------------
// open a prep file
// ---------------------------------------------------------
// Format :
// =========================================================
int prepParser::openFile(const std::string &prepfile)
{
    iprep.open(prepfile.c_str());

    //std::cout << "\n\n  prepParser::openFile " << prepfile << std::endl;
    if (!iprep) {
      setError(1);
      std::string errorMessage = "  Error, Can't Find " + prepfile;
      setErrorMessage(errorMessage);
      errorLogger.throwError("prepParser::Read", errorMessage, 1);
      return 1;
    }

    std::string fileline;

    // SKIP FIRST TWO LINES
    getline(iprep,fileline);
    getline(iprep,fileline);

    return 0;
}

// =========================================================
// Function : read fragment header section
// ---------------------------------------------------------
// 
// ---------------------------------------------------------
// Format :
// =========================================================
int prepParser::readHeader(std::string &name, std::string &symbol, double &charge)
{
    std::string fileline;

    getline(iprep,fileline);
    if (stripString(fileline, " ") == "STOP") return 1;

    name = stripString(fileline," ");

    // blank line
    getline(iprep,fileline);

    // symbol
    getline(iprep,fileline);
    std::vector<std::string> splitLine;
    splitString(fileline, " ", splitLine, 0);
    symbol = splitLine[0];

    // 
    getline(iprep,fileline);

    // charge
    getline(iprep,fileline);
    charge = string2Double(stripString(fileline, " "));

    return 0;
}

// =========================================================
// Function : read fragment main section
// ---------------------------------------------------------
// 
// ---------------------------------------------------------
// Format :
// =========================================================
void prepParser::readFragment(stdFrag* pStdFrag)
{
    std::string fileline;

    std::vector<stdBond*> todos;

    while (iprep) {
      getline(iprep,fileline);

      std::string filelineTest = stripString(fileline, " ");
      if (filelineTest == "") break;

//std::cout << "|" << fileline << "|" << std::endl;

      std::vector<std::string> splitstring;
      splitString(fileline, " ", splitstring, 0);
      if (atoi(splitstring[0].c_str()) > 3) {
        stdAtom* pStdAtom = pStdFrag->addStdAtom();
        int atIndex = atoi(splitstring[0].c_str()) - 3;

////////
        //std::string atomName = fileline.substr(6,4);
        std::string atomName = splitstring[1];
        if (atomName.size() == 1) atomName = " " + atomName + "  ";
        else if (atomName.size() == 2) atomName = " " + atomName + " ";
        else if (atomName.size() == 3) atomName = " " + atomName;

//std::cout << "|" << atomName << "|" << std::endl;
/////////

        std::string atomType = splitstring[2];
        std::string atomChain = splitstring[3];

        int b12 = atoi(splitstring[4].c_str()) - 3;
        int b13 = atoi(splitstring[5].c_str()) - 3;
        int b14 = atoi(splitstring[6].c_str()) - 3;

        if (atIndex <= 0) atIndex++;
        if (b12 <= 0) b12--;
        if (b13 <= 0) b13--;
        if (b14 <= 0) b14--;

        double bondLength = strtod(splitstring[7].c_str(), 0);
        double bondAngle = strtod(splitstring[8].c_str(), 0);
        double bondTorsion = strtod(splitstring[9].c_str(), 0);

        double atomCharge = 0.0;
        if (splitstring.size() > 10) {
          atomCharge = strtod(splitstring[10].c_str(), 0);
        }

        pStdAtom->index = atIndex;
        pStdAtom->identity = atomName;
        pStdAtom->atSymbol = determineElement(atomName);
        pStdAtom->type = atomType;
        pStdAtom->chain = atomChain;
        pStdAtom->kind = 0;

        pStdAtom->bond12 = b12;
        pStdAtom->bond13 = b13;
        pStdAtom->bond14 = b14;
        pStdAtom->bondLength = bondLength;
        pStdAtom->bondAngle = bondAngle;
        pStdAtom->bondTorsion = bondTorsion;
        pStdAtom->atmCharge = atomCharge;
        stdBond* pStdBond = pStdFrag->addStdBond();
        pStdBond->atom1 = pStdAtom->index;
        pStdBond->atom2 = pStdAtom->bond12;
        pStdBond->type = 1;

        if (pStdAtom->atSymbol == "H") {
          if (pStdAtom->bond12 > 0) {
            stdAtom* pStdBondedAtom = pStdFrag->getStdAtom(pStdAtom->bond12);

            if (pStdBondedAtom) {
              std::string bondedAtomSymbol = pStdBondedAtom->atSymbol;
              if (bondedAtomSymbol == "N" or bondedAtomSymbol == "O" or bondedAtomSymbol == "S") {
                pStdBond->kind = 1;
              }
            }
            else {
              todos.push_back(pStdBond);
            }
          }
        }
        else {
          pStdBond->kind = 0;
        }

        pStdBond->topology = 2;
        pStdBond->stereo = 0;
        pStdBond->length = pStdAtom->bondLength;
      }
    }

    for (unsigned int x = 0; x < todos.size(); x++) {
      stdAtom* pStdBondedAtom1 = pStdFrag->getStdAtom(todos[x]->atom1);
      stdAtom* pStdBondedAtom2 = pStdFrag->getStdAtom(todos[x]->atom2);
      if (pStdBondedAtom1 and pStdBondedAtom2) {
        std::string bondedAtomSymbol = pStdBondedAtom2->atSymbol;
        if (bondedAtomSymbol == "N" or bondedAtomSymbol == "O" or bondedAtomSymbol == "S") {
          todos[x]->kind = 1;
        }
      }
    }

    // Read charge, alias, loop, and improper data
    while (iprep) {
      getline(iprep,fileline);
//std::cout << "!" << fileline << "!" << std::endl;
      if (fileline != "") {
        if (fileline.substr(0,6) == "CHARGE") {
          while (iprep) {
            getline(iprep,fileline);

            std::string filelineTest = stripString(fileline, " ");
            if (filelineTest == "") break;
          }
        }
        else if (fileline.substr(0,5) == "ALIAS")  {
          while (iprep) {
            getline(iprep,fileline);
            if (fileline == "") break;

            std::vector<std::string> atomNames;
            splitString(fileline, " ", atomNames, 0);

            std::vector< std::vector<std::string> > possibleNames;
            for (int d = 0; d < 2; d++) {
              int atNameSize = atomNames[d].size();
              std::string name = "";
              std::vector<std::string> atns;
              if (atNameSize == 1) {
                name = " " + atomNames[d] + "  ";
                atns.push_back(name);
                name = atomNames[d] + "   ";
                atns.push_back(name);
              }
              else if (atNameSize == 2) {
                name = " " + atomNames[d] + " ";
                atns.push_back(name);
                name = atomNames[d] + "  ";
                atns.push_back(name);
              }
              else if (atNameSize == 3) {
                name = " " + atomNames[d];
                atns.push_back(name);
                name = atomNames[d] + " ";
                atns.push_back(name);
              }
              possibleNames.push_back(atns);
            }

            std::vector<stdAtom*> stdAtoms;
            for (int d = 0; d < 2; d++) {
              for (unsigned int d2 = 0; d2 < possibleNames[d].size(); d2++) {
                if (pStdFrag->getStdAtom(possibleNames[d][d2])) {
                  stdAtoms.push_back(pStdFrag->getStdAtom(possibleNames[d][d2]));
                  break;
                }
              }
            }

            if (stdAtoms.size() == 1) {
              for (unsigned int d = 0; d < possibleNames[1].size(); d++) {
                stdAlias* pStdAlias = pStdFrag->addStdAlias();
                pStdAlias->atom1 = stdAtoms[0]->identity;
                pStdAlias->atom2 = possibleNames[1][d];
              }
            }
            else {
              std::string errorMessage = "  ALIAS Tag Error ... exiting ";
              errorLogger.throwError("prepParser::Read", errorMessage, 1);

              std::stringstream ss;
              ss << "prepParser::Read"<< errorMessage;
              throw MTKException(ss.str());
            }

          }
        }
        else if (fileline.substr(0,4) == "LOOP")  {
          while (iprep) {
            getline(iprep,fileline);

            std::string filelineTest = stripString(fileline, " ");
            if (filelineTest == "") break;

            //if (fileline == "") break;

            std::vector<std::string> atomNames;
            splitString(fileline, " ", atomNames, 0);
/*
for (unsigned int x = 0; x <atomNames.size();x++){
  std::cout << atomNames[x] <<std::endl;
}
pStdFrag->print();
*/
            std::vector< std::vector<std::string> > possibleNames;
            for (int d = 0; d < 2; d++) {
              int atNameSize = atomNames[d].size();
              std::string name = "";
              std::vector<std::string> atns;
              if (atNameSize == 1) {
                name = " " + atomNames[d] + "  ";
                atns.push_back(name);
                name = atomNames[d] + "   ";
                atns.push_back(name);
              }
              else if (atNameSize == 2) {
                name = " " + atomNames[d] + " ";
                atns.push_back(name);
                name = atomNames[d] + "  ";
                atns.push_back(name);
              }
              else if (atNameSize == 3) {
                name = " " + atomNames[d];
                atns.push_back(name);
                name = atomNames[d] + " ";
                atns.push_back(name);
              }
              possibleNames.push_back(atns);
            }

            std::vector<stdAtom*> stdAtoms;
            for (int d = 0; d < 2; d++) {
              for (unsigned int d2 = 0; d2 < possibleNames[d].size(); d2++) {
                if (pStdFrag->getStdAtom(possibleNames[d][d2])) {
                  stdAtoms.push_back(pStdFrag->getStdAtom(possibleNames[d][d2]));
                  break;
                }
              }
            }

            if (stdAtoms.size() == 2) {
              stdLoop* pStdLoop = pStdFrag->addStdLoop();
              pStdLoop->atom1 = stdAtoms[0]->index;
              pStdLoop->atom2 = stdAtoms[1]->index;

              pStdLoop->type = 0;
              pStdLoop->stereo = 0;
            }
            else {
              std::string errorMessage = "  LOOP Tag Error ... exiting ";
              errorLogger.throwError("prepParser::Read", errorMessage, 1);
              std::stringstream ss;
              ss << "prepParser::Read"<< errorMessage;
              throw MTKException(ss.str());
            }
          }
        }
        else if (fileline.substr(0,4) == "DONE")  {
          break;
        }
        else if (fileline.substr(0,8) == "IMPROPER")  {
          while (iprep) {
            getline(iprep,fileline);

            std::string filelineTest = stripString(fileline, " ");
            if (filelineTest == "") break;

            if (fileline == "") break;

            std::vector<std::string> atomNames;
            splitString(fileline, " ", atomNames, 0);

            std::vector< std::vector<std::string> > possibleNames;
            for (int d = 0; d < 4; d++) {
              int atNameSize = atomNames[d].size();
              std::string name = "";
              std::vector<std::string> atns;
              if (atNameSize == 1) {
                name = " " + atomNames[d] + "  ";
                atns.push_back(name);
                name = atomNames[d] + "   ";
                atns.push_back(name);
              }
              else if (atNameSize == 2) {
                name = " " + atomNames[d] + " ";
                atns.push_back(name);
                name = atomNames[d] + "  ";
                atns.push_back(name);
              }
              else if (atNameSize == 3) {
                name = " " + atomNames[d];
                atns.push_back(name);
                name = atomNames[d] + " ";
                atns.push_back(name);
              }
              possibleNames.push_back(atns);
            }

            std::vector<stdAtom*> stdAtoms;
            std::vector<int> atIndex;
            for (int d = 0; d < 4; d++) {

              bool bGotIt = false;
              for (unsigned int d2 = 0; d2 < possibleNames[d].size(); d2++) {
                if (d == 0) {
                  if (containsSubStr(possibleNames[d][d2], "-M")) {
                    atIndex.push_back(-1);
                    stdAtoms.push_back(0);
                    bGotIt = true;
                    break;
                  }
                }

                if (d == 0) {
                  if (containsSubStr(possibleNames[d][d2], "+M")) {
                    atIndex.push_back(-4);
                    stdAtoms.push_back(0);
                    bGotIt = true;
                    break;
                  }
                }

                if (d == 1) {
                  if (containsSubStr(possibleNames[d][d2], "+M") or
                      containsSubStr(possibleNames[d][d2], "-M")
                     ) {
                    atIndex.push_back(-5);
                    stdAtoms.push_back(0);
                    bGotIt = true;
                    break;
                  }
                }

                if (pStdFrag->getStdAtom(possibleNames[d][d2])) {
                  stdAtoms.push_back(pStdFrag->getStdAtom(possibleNames[d][d2]));
                  atIndex.push_back(pStdFrag->getStdAtom(possibleNames[d][d2])->index);
                  bGotIt = true;
                  break;
                }
              }
              if (!bGotIt) {
                stdAtoms.push_back(0);
                atIndex.push_back(-5);
              }
            }

            bool bOK = true;
            for (int d = 0; d < 4; d++) {
              if (atIndex[d] == -5) {
                bOK = false;
              }
            }

            if (bOK) {
              stdImproper* pStdImproper = pStdFrag->addStdImproper();

              if (pStdImproper) {
                pStdImproper->atom1 = atIndex[0];
                pStdImproper->atom2 = atIndex[1];
                pStdImproper->atom3 = atIndex[2];
                pStdImproper->atom4 = atIndex[3];
              }
              else {
                std::string errorMessage = "  IMPROPER Tag Error ... exiting ";
                errorLogger.throwError("prepParser::Read", errorMessage, 1);

                std::stringstream ss;
                ss << "prepParser::Read" << errorMessage;
                throw MTKException(ss.str());
              }
            }
          }
        }
      }
    }

    std::string eMessage = " Atom Name Aliases \n";

    // Add possible atom name aliases
    std::vector<stdAtom*> stdAtoms = pStdFrag->getStdAtomList();
    for (unsigned int i = 0; i < stdAtoms.size(); i++) {
      std::string itsType = "";
      std::string itsName = stdAtoms[i]->identity;
      for (unsigned int j = 0; j < itsName.size(); j++) {
        if (isalpha(itsName[j])) {
          itsType += "C"; // Character
        }
        else if (itsName[j] == ' ') {
          itsType += "B"; // Blank
        }
        else {
          itsType += "N"; // Number
        }
      }
      eMessage += " |" + itsName + "| " + itsType + " \n";

      std::string blankString = " ";
      std::string aliasName = "";

      if (itsType == "BCNB") {
        continue;
      }
      else if (itsType == "CNBB") {
        aliasName = "    ";
        aliasName[1] = itsName[0];
        aliasName[2] = itsName[1];
        if (!pStdFrag->hasStdAtom(aliasName)) {
          stdAlias* pStdAlias = pStdFrag->addStdAlias();
          pStdAlias->atom1 = std::string(itsName);
          pStdAlias->atom2 = std::string(aliasName);
          eMessage += "   |" + aliasName + "|\n";
        }
      }
      else if (itsType == "CNNB") {
        aliasName = "    ";
        aliasName[1] = itsName[0];
        aliasName[2] = itsName[1];
        aliasName[3] = itsName[2];
        if (!pStdFrag->hasStdAtom(aliasName)) {
          stdAlias* pStdAlias = pStdFrag->addStdAlias();
          pStdAlias->atom1 = std::string(itsName);
          pStdAlias->atom2 = std::string(aliasName);
          eMessage += "   |" + aliasName + "|\n";
        }
      }
      else if (itsType == "BCNN") {
        aliasName = "    ";
        aliasName[0] = itsName[1];
        aliasName[1] = itsName[2];
        aliasName[2] = itsName[3];
        if (!pStdFrag->hasStdAtom(aliasName)) {
          stdAlias* pStdAlias = pStdFrag->addStdAlias();
          pStdAlias->atom1 = std::string(itsName);
          pStdAlias->atom2 = std::string(aliasName);
          eMessage += "   |" + aliasName + "|\n";
        }

        aliasName = "    ";
        aliasName[0] = itsName[3];
        aliasName[1] = itsName[1];
        aliasName[2] = itsName[2];
        if (!pStdFrag->hasStdAtom(aliasName)) {
          stdAlias* pStdAlias = pStdFrag->addStdAlias();
          pStdAlias->atom1 = std::string(itsName);
          pStdAlias->atom2 = std::string(aliasName);
          eMessage += "   |" + aliasName + "|\n";
        }

        aliasName = "    ";
        aliasName[1] = itsName[3];
        aliasName[2] = itsName[1];
        aliasName[3] = itsName[2];
        if (!pStdFrag->hasStdAtom(aliasName)) {
          stdAlias* pStdAlias = pStdFrag->addStdAlias();
          pStdAlias->atom1 = std::string(itsName);
          pStdAlias->atom2 = std::string(aliasName);
          eMessage += "   |" + aliasName + "|\n";
        }
      }
      else if (itsType == "BCCN") {
        aliasName = "    ";
        aliasName[0] = itsName[1];
        aliasName[1] = itsName[2];
        aliasName[2] = itsName[3];
        if (!pStdFrag->hasStdAtom(aliasName)) {
          stdAlias* pStdAlias = pStdFrag->addStdAlias();
          pStdAlias->atom1 = std::string(itsName);
          pStdAlias->atom2 = std::string(aliasName);
          eMessage += "   |" + aliasName + "|\n";
        }

        aliasName = "    ";
        aliasName[0] = itsName[3];
        aliasName[1] = itsName[1];
        aliasName[2] = itsName[2];
        if (!pStdFrag->hasStdAtom(aliasName)) {
          stdAlias* pStdAlias = pStdFrag->addStdAlias();
          pStdAlias->atom1 = std::string(itsName);
          pStdAlias->atom2 = std::string(aliasName);
          eMessage += "   |" + aliasName + "|\n";
        }

        aliasName = "    ";
        aliasName[1] = itsName[3];
        aliasName[2] = itsName[1];
        aliasName[3] = itsName[2];
        if (!pStdFrag->hasStdAtom(aliasName)) {
          stdAlias* pStdAlias = pStdFrag->addStdAlias();
          pStdAlias->atom1 = std::string(itsName);
          pStdAlias->atom2 = std::string(aliasName);
          eMessage += "   |" + aliasName + "|\n";
        }
      }
      else if (itsType == "CCNN") {
        aliasName = "    ";
        aliasName[0] = itsName[2];
        aliasName[1] = itsName[0];
        aliasName[2] = itsName[1];
        aliasName[3] = itsName[3];
        if (!pStdFrag->hasStdAtom(aliasName)) {
          stdAlias* pStdAlias = pStdFrag->addStdAlias();
          pStdAlias->atom1 = std::string(itsName);
          pStdAlias->atom2 = std::string(aliasName);
          eMessage += "   |" + aliasName + "|\n";
        }

        aliasName = "    ";
        aliasName[0] = itsName[3];
        aliasName[1] = itsName[0];
        aliasName[2] = itsName[1];
        aliasName[3] = itsName[2];
        if (!pStdFrag->hasStdAtom(aliasName)) {
          stdAlias* pStdAlias = pStdFrag->addStdAlias();
          pStdAlias->atom1 = std::string(itsName);
          pStdAlias->atom2 = std::string(aliasName);
          eMessage += "   |" + aliasName + "|\n";
        }
      }
      else if (itsType == "CNNN") {
        aliasName = "    ";
        aliasName[0] = itsName[3];
        aliasName[1] = itsName[0];
        aliasName[2] = itsName[1];
        aliasName[3] = itsName[2];
        if (!pStdFrag->hasStdAtom(aliasName)) {
          stdAlias* pStdAlias = pStdFrag->addStdAlias();
          pStdAlias->atom1 = std::string(itsName);
          pStdAlias->atom2 = std::string(aliasName);
          eMessage += "   |" + aliasName + "|\n";
        }
      }
// QBio start
      else if (itsType == "CNCN") {
        aliasName = "    ";
        aliasName[0] = itsName[3];
        aliasName[1] = itsName[0];
        aliasName[2] = itsName[1];
        aliasName[3] = itsName[2];
        if (!pStdFrag->hasStdAtom(aliasName)) {
          stdAlias* pStdAlias = pStdFrag->addStdAlias();
          pStdAlias->atom1 = std::string(itsName);
          pStdAlias->atom2 = std::string(aliasName);
        }
      }
      else if (itsType == "CCCN") {
        aliasName = "    ";
        aliasName[0] = itsName[3];
        aliasName[1] = itsName[0];
        aliasName[2] = itsName[1];
        aliasName[3] = itsName[2];
        if (!pStdFrag->hasStdAtom(aliasName)) {
          stdAlias* pStdAlias = pStdFrag->addStdAlias();
          pStdAlias->atom1 = std::string(itsName);
          pStdAlias->atom2 = std::string(aliasName);
        }
      }
      else if (itsType == "CCNB") {
        aliasName = "    ";
        aliasName[0] = itsName[2];
        aliasName[1] = itsName[0];
        aliasName[2] = itsName[1];
        aliasName[3] = itsName[3];
        if (!pStdFrag->hasStdAtom(aliasName)) {
          stdAlias* pStdAlias = pStdFrag->addStdAlias();
          pStdAlias->atom1 = std::string(itsName);
          pStdAlias->atom2 = std::string(aliasName);
        }
      }
      else if (itsType == "CNNB") {
        aliasName = "    ";
        aliasName[0] = itsName[2];
        aliasName[1] = itsName[0];
        aliasName[2] = itsName[1];
        aliasName[3] = itsName[3];
        if (!pStdFrag->hasStdAtom(aliasName)) {
          stdAlias* pStdAlias = pStdFrag->addStdAlias();
          pStdAlias->atom1 = std::string(itsName);
          pStdAlias->atom2 = std::string(aliasName);
        }
      }
// QBio end
    }
    errorLogger.throwError("prepParser::Read", eMessage, INFO);
}

// =========================================================
// Function : Read
// ---------------------------------------------------------
// parsers a prep file
// ---------------------------------------------------------
// Format :
// =========================================================
void prepParser::Read(const std::string &prepfile, stdGroup* pStdGroup)
{
    int r = openFile(prepfile);
    if (r == 1) return;

    while(iprep) {
      std::string name = "";
      std::string symbol = "";
      double charge = 0.0;
      int h = readHeader(name, symbol, charge);
      if (h == 1) break;

      stdFrag* pStdFrag = pStdGroup->addStdFrag();
      if (!pStdFrag) {
        setError(1);
        std::string errorMessage = "  Error creating fragment ";
        setErrorMessage(errorMessage);
        errorLogger.throwError("prepParser::Read", errorMessage, 1);
        return;
      }

      if (symbol.size() == 1) symbol = symbol + "  ";
      if (symbol.size() == 2) symbol = symbol + " ";
      if (symbol.size() > 3) symbol = symbol.substr(1,3);

      pStdFrag->setSymbol(symbol);       // 3L code
      //pStdFrag->setCode("P2XML"+name); // 8L code
      pStdFrag->setName(name);           // long name
      pStdFrag->setType("m");            // fragment type

      //std::cout << " prepParser::Read " << name << "  charge " << charge << " symbol " << symbol << std::endl;
      readFragment(pStdFrag);
    }

    iprep.close();

    return;
}

// =========================================================
// Function : Read
// ---------------------------------------------------------
// parsers a prep file
// ---------------------------------------------------------
// Format :
// =========================================================
void prepParser::Read(const std::string &prepfile, stdFrag* pStdFrag)
{
    errorLogger.throwError("prepParser", "Read", INFO);
/*
    std::ifstream iprep;
    iprep.open(prepfile.c_str());

    if (!iprep) {
      setError(1);
      std::string errorMessage = "  Error, Can't Find " + prepfile;
      setErrorMessage(errorMessage);
      errorLogger.throwError("prepParser::Read", errorMessage, 1);
      return;
    }
*/

    int r = openFile(prepfile);
    if (r == 1) return;

    std::string fileline;
    std::string title;
    std::string atom;
    //std::string buffer(80,'*');

    // SKIP FIRST TWO LINES
    //getline(iprep,fileline);
    //getline(iprep,fileline);

/*
    getline(iprep,fileline);
    //std::string name = stripString(fileline," ");
    //pStdFrag->setName(fileline);
    getline(iprep,fileline);
    getline(iprep,fileline);
    std::vector<std::string> splitLine;
    splitString(fileline, " ", splitLine, 0);
    pStdFrag->setSymbol(splitLine[0]);
    getline(iprep,fileline);
    getline(iprep,fileline);
    //pStdFrag->setCharge(splitLine[0]);
*/

    std::string name = "";
    std::string symbol = "";
    double charge = 0.0;
    readHeader(name, symbol, charge);

    pStdFrag->setSymbol(symbol);

    readFragment(pStdFrag);

    iprep.close();
}

// ==============================================================
// Function : Write
// --------------------------------------------------------------
// Write a AMBER prep file.
// ==============================================================
void prepParser::Write(const std::string &prepfile, stdFrag* stdF)
{
    std::ofstream oprep;
    oprep.open(prepfile.c_str());

    if (!oprep or (stdF == 0)) {
      setError(1);
      std::string errorMessage = "  Error, Can't Open " + prepfile;
      setErrorMessage(errorMessage);
      errorLogger.throwError("prepParser::Write", errorMessage, 1);
      return;
    }

    oprep.close();
}

// ==============================================================
// Function : Write
// --------------------------------------------------------------
// Write a AMBER prep file.
// ==============================================================
void prepParser::Write(const std::string &prepfile, stdGroup* stdG)
{
    std::ofstream oprep;
    oprep.open(prepfile.c_str());

    if (!oprep or (stdG == 0)) {
      setError(1);
      std::string errorMessage = "  Error, Can't Open " + prepfile;
      setErrorMessage(errorMessage);
      errorLogger.throwError("prepParser::Write", errorMessage, 1);
      return;
    }

    std::vector<stdFrag*> frags = stdG->getStdFragList();

    oprep << "    1    1    2" << std::endl;

    // Title
    oprep << stdG->getName() << " fragment set created by MTK++/MCPB " << std::endl;

    for (unsigned int i = 0; i < frags.size(); i++) {
      stdFrag* pStdFrag = frags[i];

      // Descriptive header for the residue
      oprep << pStdFrag->getName() << "\n" << std::endl;
      oprep << std::right << std::setw(4) << pStdFrag->getSymbol() << " INT 1" << std::endl;

      oprep << " CORR OMIT DU   BEG " << std::endl;
      oprep << "0.00000" << std::endl;

      std::vector<stdAtom*> atoms = pStdFrag->getStdAtomList();
      oprep << "   1  DUMM  DU    M    0  -1  -2     0.000      .0        .0      .00000" << std::endl;
      oprep << "   2  DUMM  DU    M    1   0  -1     1.449      .0        .0      .00000" << std::endl;
      oprep << "   3  DUMM  DU    M    2   1   0     1.522   111.1        .0      .00000" << std::endl;

      for (unsigned int j = 0; j < atoms.size(); j++) {
        stdAtom* pStdAtom = atoms[j];
        int b12 = 0;
        int b13 = 0;
        int b14 = 0;

        if (pStdAtom->bond12 < 0) {
          b12 = pStdAtom->bond12 + 4;
        }
        else {
          b12 = pStdAtom->bond12 + 3;
        }

        if (pStdAtom->bond13 < 0) {
          b13 = pStdAtom->bond13 + 4;
        }
        else {
          b13 = pStdAtom->bond13 + 3;
        }

        if (pStdAtom->bond14 < 0) {
          b14 = pStdAtom->bond14 + 4;
        }
        else {
          b14 = pStdAtom->bond14 + 3;
        }

        oprep << std::right << std::setw(4) << j+4 << " " << pStdAtom->identity << "   "
              << std::left << std::setw(2) << pStdAtom->type << "    "
              << pStdAtom->chain << " "
              << std::right
              << std::setw(4) << b12
              << std::setw(4) << b13
              << std::setw(4) << b14
              << std::setiosflags(std::ios::fixed) << std::setprecision(3) << std::setw(10) << pStdAtom->bondLength
              << std::setiosflags(std::ios::fixed) << std::setprecision(3) << std::setw(10) << pStdAtom->bondAngle
              << std::setiosflags(std::ios::fixed) << std::setprecision(3) << std::setw(10) << pStdAtom->bondTorsion
             << std::setiosflags(std::ios::fixed) << std::setprecision(5) << std::setw(10) << pStdAtom->atmCharge
             << std::endl;
      }

      oprep << " " << std::endl;

      std::vector<stdLoop*> bonds = pStdFrag->getStdLoopList();
      if (bonds.size() > 0) {
        oprep << "LOOP" << std::endl;
        for (unsigned int j = 0; j < bonds.size(); j++) {
          int at1 = bonds[j]->atom1;
          int at2 = bonds[j]->atom2;

          stdAtom* pStdAtom1 = pStdFrag->getStdAtom(at1);
          stdAtom* pStdAtom2 = pStdFrag->getStdAtom(at2);

          oprep << std::setw(4) << pStdAtom1->identity 
          << std::setw(5) << pStdAtom2->identity
          << std::endl;
        }
        oprep << " " << std::endl;
      }

      std::vector<stdImproper*> imps = pStdFrag->getStdImproperList();
      if (imps.size() > 0) {
        oprep << "IMPROPER" << std::endl;
        for (unsigned int j = 0; j < imps.size(); j++) {
          int at1 = imps[j]->atom1;
          int at2 = imps[j]->atom2;
          int at3 = imps[j]->atom3;
          int at4 = imps[j]->atom4;
/*
std::cout << " prepParser: Improper " << at1 << " " << at2 << " "
          << at3 << " " << at4 << std::endl;
*/
          stdAtom* pStdAtom1 = 0;
          stdAtom* pStdAtom2 = 0;
          stdAtom* pStdAtom3 = 0;
          stdAtom* pStdAtom4 = 0;

          std::string at1Name = "";
          std::string at2Name = "";
          std::string at3Name = "";
          std::string at4Name = "";

          if (at1 == -1) {
            at1Name = " -M ";
          }
          else {
            pStdAtom1 = pStdFrag->getStdAtom(at1);
            if (pStdAtom1) {
              at1Name = pStdAtom1->identity;
            }
          }

          if (at2 == -4) {
            at2Name = " +M ";
          }
          else {
            pStdAtom2 = pStdFrag->getStdAtom(at2);
            if (pStdAtom2) {
              at2Name = pStdAtom2->identity;
            }
          }

          pStdAtom3 = pStdFrag->getStdAtom(at3);
          if (pStdAtom3) {
            at3Name = pStdAtom3->identity;
          }

          pStdAtom4 = pStdFrag->getStdAtom(at4);
          if (pStdAtom4) {
            at4Name = pStdAtom4->identity;
          }
/*
std::cout << "    " 
          << at1Name << " "
          << at2Name << " "
          << at3Name << " "
          << at4Name
          << std::endl;
*/
          oprep << std::setw(4) << at1Name
                << std::setw(5) << at2Name
                << std::setw(5) << at3Name
                << std::setw(5) << at4Name << std::endl;
        }
        oprep << " " << std::endl;
      }

      oprep << "DONE" << std::endl;
    }
    oprep << "STOP" << std::endl;

    oprep.close();
}

} // MTKpp namespace

