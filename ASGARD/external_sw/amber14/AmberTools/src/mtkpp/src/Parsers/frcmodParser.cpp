/*!
   \file frcmodParser.cpp
   \brief Parses AMBER frcmod files
   \author Martin Peters

   $Date: 2010/03/29 20:39:34 $
   $Revision: 1.7 $

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

#include "frcmodParser.h"
#include "Molecule/parameters.h"
#include "StringManip.h"
#include "Utils/constants.h"

#ifdef USE_XERCES
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#endif

#include "Diagnostics/MTKException.h"

namespace MTKpp
{

// ============================================================
// Function : frcmodParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
frcmodParser::frcmodParser(parameters *c, std::string grpName)
{
    this->pParameters = c;
    this->groupName = grpName;
    if (!pParameters) {
      std::cout << " Error in frcmodParser ... exiting " << std::endl;
      //exit(0);
      throw MTKException(" Error in frcmodParser ... exiting ");
    }
}

// =========================================================
// Function : frcmodParser()
// ---------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// =========================================================
frcmodParser::~frcmodParser() {}

// =========================================================
// Function : Read
// ---------------------------------------------------------
// parsers an AMBER frcmod file
// ---------------------------------------------------------
// Format :
// =========================================================
void frcmodParser::Read(const std::string &frcmodfile)
{
    std::ifstream ifrcmod;
    ifrcmod.open(frcmodfile.c_str());

    if (!ifrcmod) {
      std::cout << "\nUnable to open frcmod file"
                << "\nFilename = " << frcmodfile
                << "\nExiting ...\n" << std::endl;
      return;
    }

    bool bMASS = false;
    bool bBOND = false;
    bool bANGLE = false;
    bool bDIHE = false;
    bool bIMPROPER = false;
    bool bNONBON = false;
    std::string fileline = "";
    while (ifrcmod) {
      std::string buffer(80,'*');
      getline(ifrcmod,fileline);
//std::cout << "|" << fileline << "|" << std::endl;
      bMASS = false;
      bBOND = false;
      bANGLE = false;
      bDIHE = false;
      bIMPROPER = false;
      bNONBON = false;

      if ( (fileline.size() > 0) and (stripString(fileline, " ") != "") ) {
        std::vector<std::string> splitstring;
        splitString(fileline, " ", splitstring, 0);

        if (splitstring[0] == "MASS") bMASS = true;
        if (splitstring[0] == "BOND") bBOND = true;
        if (splitstring[0] == "ANGL" or splitstring[0] == "ANGLE") bANGLE = true;
        if (splitstring[0] == "DIHE" or splitstring[0] == "DIHEDRAL") bDIHE = true;
        if (splitstring[0] == "IMPROP" or  splitstring[0] == "IMPROPER") bIMPROPER = true;
        if (splitstring[0] == "NONB" or splitstring[0] == "NONBON") bNONBON = true;
      }
/*
      try {
        if (fileline.substr(0,4) == "MASS") bMASS = true;
      }
      catch (std::out_of_range) {
         std::cout << " out_of_range MAS " << std::endl;
      }

      try {
        if (fileline.substr(0,4) == "BOND") bBOND = true;
      }
      catch (std::out_of_range) {
         std::cout << " out_of_range BOND " << std::endl;
      }

      try {
        if (fileline.substr(0,5) == "ANGLE") bANGLE = true;
      }
      catch (std::out_of_range) {
         std::cout << " out_of_range ANG" << std::endl;
      }

      try {
        if (fileline.substr(0,4) == "DIHE") bDIHE = true;
      }
      catch (std::out_of_range) {
         std::cout << " out_of_range DIHE " << std::endl;
      }

      try {
        if (fileline.substr(0,8) == "IMPROPER") bIMPROPER = true;
      }
      catch (std::out_of_range) {
         std::cout << " out_of_range IMP " << std::endl;
      }

      try {
        if (fileline.substr(0,6) == "NONBON") bNONBON = true;
      }
      catch (std::out_of_range) {
         std::cout << " out_of_range NONB" << std::endl;
      }
*/

      if (bMASS) {
        getline(ifrcmod, fileline);
        while (stripString(fileline, " ") != "") {
        //while (GetAlphaChar(fileline, 0) != "0") {
          std::vector<std::string> splitLine;
          splitString(fileline, " ", splitLine, 0);
          if (splitLine.size() >= 2) {
            for (unsigned int i = 0; i < 2; i++) {
              splitLine[i] = stripString(splitLine[i], " ");
            }
            pAtomType = pParameters->addAtomType();
            pAtomType->name = splitLine[0];
            pAtomType->groupName = this->groupName;
            pAtomType->mass = string2Double(splitLine[1]); //strtod(fileline.substr(3,5).c_str(), 0);

            if (splitLine.size() > 2) {
              char *end;
              double ret_double = strtod(splitLine[2].c_str(), &end);
              if (!*end) {
                pAtomType->atomPolarizability = ret_double;
              }
//std::cout << splitLine[2] << " " << ret_double << std::endl;
            }

            getline(ifrcmod,fileline);
          }
        }
      }

      if (bBOND) {
        getline(ifrcmod,fileline);
        //while (fileline != "") {
        while (stripString(fileline, " ") != "") {
          //std::string atomNames = fileline.substr(0,5);
          std::vector<std::string> splitLine;
          splitString(fileline, "-", splitLine, 0);

          std::vector<std::string> splitLine2;
          splitString(splitLine[1], " ", splitLine2, 0);

          //if (splitLine.size() == 2) {
            //for (unsigned int i = 0; i < 2; i++) {
            //  splitLine[i] = stripString(splitLine[i]," ");
            //}
            pBondParam = pParameters->addBondParam();
            pBondParam->groupName = this->groupName;
            pBondParam->atomType1 = stripString(splitLine[0], " ");
            pBondParam->atomType2 = stripString(splitLine2[0], " ");
//std::cout << "|" << pBondParam->atomType1 << "|" << pBondParam->atomType2 << "|" << std::endl;
            pBondParam->keq = string2Double(splitLine2[1]);
            pBondParam->req = string2Double(splitLine2[2]);

//            pBondParam->keq = strtod(fileline.substr(7,5).c_str(), 0);
//            pBondParam->req = strtod(fileline.substr(16,5).c_str(), 0);

//std::cout << "\n " << fileline << std::endl;
//std::cout << " BOND " << pBondParam->atomType1 << " " << pBondParam->atomType2 << " " << pBondParam->keq << " " << pBondParam->req << std::endl;
//std::cout << " BOND " << fileline.substr(7,5) << " "  << fileline.substr(16,5) <<std::endl;

            getline(ifrcmod,fileline);
          //}
        }
      }

      if (bANGLE) {
        getline(ifrcmod,fileline);
        //while (fileline != "") {
        while (stripString(fileline, " ") != "") {
//          std::string atomNames = fileline.substr(0,8);

          std::vector<std::string> splitLine;
          splitString(fileline, "-", splitLine, 0);

          std::vector<std::string> splitLine2;
          splitString(splitLine[2], " ", splitLine2, 0);

          //std::vector<std::string> splitLine;
          //splitString(atomNames, "-", splitLine, 0);

          //std::vector<std::string> splitLine2;
          //splitString(fileline, " ", splitLine2, 0);

          //if (splitLine.size() == 3) {
            //for (unsigned int i = 0; i < 3; i++) {
            //  splitLine[i] = stripString(splitLine[i]," ");
            //}
            pAngleParam = pParameters->addAngleParam();
            pAngleParam->groupName = this->groupName;
            pAngleParam->atomType1 = stripString(splitLine[0], " ");
            pAngleParam->atomType2 = stripString(splitLine[1], " ");
            pAngleParam->atomType3 = stripString(splitLine2[0], " ");
//std::cout << "|" << pAngleParam->atomType1 << "|" << pAngleParam->atomType2 << "|" << pAngleParam->atomType3 << "|" << std::endl;

            pAngleParam->keq = string2Double(splitLine2[1]);
            pAngleParam->req = string2Double(splitLine2[2]) * DEG2RAD;

            //pAngleParam->keq = strtod(fileline.substr(11,7).c_str(), 0);
            //pAngleParam->req = strtod(fileline.substr(20,8).c_str(), 0) * DEG2RAD;

//std::cout << "\n " << fileline << std::endl;
//std::cout << " ANGLE " << pAngleParam->atomType1 << " " << pAngleParam->atomType2 << " " << pAngleParam->atomType3 << " " << pAngleParam->keq << " " << string2Double(splitLine2[2]) << std::endl;

            getline(ifrcmod,fileline);
          //}
        }
      }

      if (bDIHE) {
        getline(ifrcmod,fileline);
        //while (fileline != "") {
        while (stripString(fileline, " ") != "") {

//std::cout << "\n " << fileline << std::endl;

//std::cout << fileline.substr(0,12) << std::endl;

// snip = str.substr(0,12); 
//(3) erase the substring of 0..where str.erase(0,where); 

          std::string atTypes = fileline.substr(0,12);
          std::vector<std::string> splitLine;
          splitString(atTypes, "-", splitLine, 0);

/*for (unsigned int i = 0; i < splitLine.size(); i++) {
  std::cout << splitLine[i] << std::endl;
}*/

          fileline.erase(0,12);

          std::vector<std::string> splitLine2;
          splitString(fileline, " ", splitLine2, 0);
/*
for (unsigned int i = 0; i < splitLine2.size(); i++) {
  std::cout << splitLine2[i] << std::endl;
}*/

          //std::string atomNames = fileline.substr(0,11);
          //std::vector<std::string> splitLine;
          //splitString(atomNames, "-", splitLine, 0);

          //if (splitLine.size() == 4) {
            //for (unsigned int i = 0; i < 4; i++) {
              //splitLine[i] = stripString(splitLine[i], " ");
            //}
            pTorsionParam = pParameters->addTorsionParam();
            pTorsionParam->groupName = this->groupName;
            pTorsionParam->atomType1 = stripString(splitLine[0], " ");
            pTorsionParam->atomType2 = stripString(splitLine[1], " ");
            pTorsionParam->atomType3 = stripString(splitLine[2], " ");
            pTorsionParam->atomType4 = stripString(splitLine[3], " ");

            pTorsionParam->npth = string2Int(splitLine2[0]);
            pTorsionParam->Vn = string2Double(splitLine2[1]);
            pTorsionParam->gamma = string2Double(splitLine2[2]) * DEG2RAD;
            pTorsionParam->Nt = string2Double(splitLine2[3]);

            //pTorsionParam->npth = atoi(fileline.substr(11,4).c_str());
            //pTorsionParam->Vn = strtod(fileline.substr(18,8).c_str(), 0);
            //pTorsionParam->gamma = strtod(fileline.substr(31,8).c_str(), 0) * DEG2RAD;
            //pTorsionParam->Nt = strtod(fileline.substr(42,12).c_str(), 0);

//std::cout << " DIHE " << pTorsionParam->atomType1 << " " << pTorsionParam->atomType2 << " " << pTorsionParam->atomType3 << " " << pTorsionParam->atomType4 << " " << pTorsionParam->npth << " " <<  pTorsionParam->Vn << " " << string2Double(splitLine2[2]) << " " << pTorsionParam->Nt << std::endl;
 
            getline(ifrcmod,fileline);
//std::cout << " DIHE LINE: " << fileline << std::endl;
          //}
        }
      }

      if (bIMPROPER) {
        getline(ifrcmod,fileline);
        //while (fileline != "") {
        while (stripString(fileline, " ") != "") {

          std::string atTypes = fileline.substr(0,12);
          std::vector<std::string> splitLine;
          splitString(atTypes, "-", splitLine, 0);

          fileline.erase(0,12);

          std::vector<std::string> splitLine2;
          splitString(fileline, " ", splitLine2, 0);

//          std::vector<std::string> splitLine;
//          splitString(fileline, "-", splitLine, 0);

//          std::vector<std::string> splitLine2;
//          splitString(splitLine[3], " ", splitLine2, 0);

          //std::string atomNames = fileline.substr(0,11);
          //std::vector<std::string> splitLine;
          //splitString(atomNames, "-", splitLine, 0);
          //if (splitLine.size() == 4) {
            //for (unsigned int i = 0; i < 4; i++) {
              //splitLine[i] = stripString(splitLine[i]," ");
            //}
            pImproperParam = pParameters->addImproperParam();
            pImproperParam->groupName = this->groupName;
            pImproperParam->atomType1 = stripString(splitLine[0], " ");
            pImproperParam->atomType2 = stripString(splitLine[1], " ");
            pImproperParam->atomType3 = stripString(splitLine[2], " ");
            pImproperParam->atomType4 = stripString(splitLine[3], " ");

            pImproperParam->Vn = string2Double(splitLine2[0]);
            pImproperParam->gamma = string2Double(splitLine2[1]) * DEG2RAD;
            pImproperParam->Nt = string2Double(splitLine2[2]);

            //pImproperParam->Vn = strtod(fileline.substr(18,8).c_str(), 0);
            //pImproperParam->gamma = strtod(fileline.substr(31,8).c_str(), 0) * DEG2RAD;
            //pImproperParam->Nt = strtod(fileline.substr(42,12).c_str(), 0);

//std::cout << "\n " << fileline << std::endl;
//std::cout << " IMPR " << pImproperParam->atomType1 << " " << pImproperParam->atomType2 << " " << pImproperParam->atomType3 << " " << pImproperParam->atomType4 << " " << pImproperParam->Vn << " " << string2Double(splitLine2[1]) << " " << pImproperParam->Nt << std::endl;

            getline(ifrcmod,fileline);
          //}
        }
      }
      if (bNONBON) {
        getline(ifrcmod,fileline);
//std::cout << fileline << std::endl;
        //while (fileline != "") {
        while (stripString(fileline, " ") != "") {
          std::vector<std::string> splitLine;
          splitString(fileline, " ", splitLine, 0);
          if (splitLine.size() >= 3) {
            for (unsigned int i = 0; i < 3; i++) {
              splitLine[i] = stripString(splitLine[i], " ");
            }
//std::cout << " |" << splitLine[0] << "|" << std::endl;

//std::cout << "\n " << fileline << std::endl;

            pAtomType = pParameters->getAtomType(stripString(splitLine[0], " "));
            if (pAtomType) {
              pAtomType->rvalue = string2Double(splitLine[1]);
              pAtomType->evalue = string2Double(splitLine[2]);

//std::cout << pAtomType->name << " " << pAtomType->rvalue << " " << pAtomType->evalue << std::endl;

            }
            else {
              std::cout << " Can't find atom type " << std::endl;
            }
            getline(ifrcmod,fileline);
          }
        }
      }
    }

    ifrcmod.close();
}

// =========================================================
// Function : Write
// ---------------------------------------------------------
// Write an AMBER frcmod file
// ---------------------------------------------------------
// Format :
// =========================================================
void frcmodParser::Write(const std::string &frcmodfile, const std::string &ps)
{
    std::ofstream ofrcmod;
    ofrcmod.open(frcmodfile.c_str());

    if (!ofrcmod or (!this->pParameters)) {
      std::cout << "\nUNABLE TO OPEN FRCMOD FILE"
                << "\nFILENAME = " << frcmodfile << std::endl;
      return;
    }

    // Title
    ofrcmod << ps << " parameter set created by MTK++/MCPB " << std::endl;

    ofrcmod << "\nMASS" << std::endl;

    std::vector<atomType*> atomTypes = this->pParameters->getAtomTypes();
    for (unsigned int a = 0; a < atomTypes.size(); a++) {
      if (atomTypes[a]->groupName == ps) {
        ofrcmod << std::left << std::setw(2) << atomTypes[a]->name 
                << std::right << std::setiosflags(std::ios::fixed) << std::setprecision(3)
                << std::setw(8) << atomTypes[a]->mass << "\n";
      }
    }

    std::vector<bondParam*> bondTypes = this->pParameters->getBondParams();
    if (bondTypes.size() > 0) {
      ofrcmod << "\nBOND" << std::endl;
    }

    for (unsigned int b = 0; b < bondTypes.size(); b++) {
      if (bondTypes[b]->groupName == ps) {
        ofrcmod << std::left << std::setw(2) << bondTypes[b]->atomType1 << "-" << std::setw(2) << bondTypes[b]->atomType2
                << std::right
                << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(8) << bondTypes[b]->keq
                << std::setiosflags(std::ios::fixed) << std::setprecision(3) << std::setw(8) << bondTypes[b]->req << "\n";
      }
    }

    std::vector<angleParam*> angleTypes = this->pParameters->getAngleParams();
    if (angleTypes.size() > 0) {
      ofrcmod << "\nANGL" << std::endl;
    }

    for (unsigned int a = 0; a < angleTypes.size(); a++) {
      if (angleTypes[a]->groupName == ps) {
        ofrcmod << std::left << std::setw(2) << angleTypes[a]->atomType1 << "-" << std::setw(2) << angleTypes[a]->atomType2
                << "-" << std::setw(2) << angleTypes[a]->atomType3
                << std::right
                << std::setiosflags(std::ios::fixed) << std::setprecision(3) << std::setw(9) << angleTypes[a]->keq
                << std::setiosflags(std::ios::fixed) << std::setprecision(3) << std::setw(12) << angleTypes[a]->req / DEG2RAD << "\n";
      }
    }

    std::vector<torsionParam*> diheTypes = this->pParameters->getTorsionParams();
    if (diheTypes.size() > 0) {
      ofrcmod << "\nDIHE" << std::endl;
    }

    for (unsigned int a = 0; a < diheTypes.size(); a++) {
      if (diheTypes[a]->groupName == ps) {
        ofrcmod << std::left << std::setw(2) << diheTypes[a]->atomType1 << "-" << std::setw(2) << diheTypes[a]->atomType2
                << "-" << std::setw(2) << diheTypes[a]->atomType3 << "-" << std::setw(2) << diheTypes[a]->atomType4
                << std::right
                << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(4) << diheTypes[a]->npth
                << std::setiosflags(std::ios::fixed) << std::setprecision(3) << std::setw(9) << diheTypes[a]->Vn 
                << std::setiosflags(std::ios::fixed) << std::setprecision(3) << std::setw(14) << diheTypes[a]->gamma / DEG2RAD
                << std::setiosflags(std::ios::fixed) << std::setprecision(3) << std::setw(16) << diheTypes[a]->Nt
                << "\n";
      }
    }

    std::vector<improperParam*> imprTypes = this->pParameters->getImproperParams();
    if (imprTypes.size() > 0) {
      ofrcmod << "\nIMPR" << std::endl;
    }

    for (unsigned int a = 0; a < imprTypes.size(); a++) {
      if (imprTypes[a]->groupName == ps) {
        ofrcmod << std::left << std::setw(2) << imprTypes[a]->atomType1 << "-" << std::setw(2) << imprTypes[a]->atomType2
                << "-" << std::setw(2) << imprTypes[a]->atomType3 << "-" << std::setw(2) << imprTypes[a]->atomType4
                << std::right
                << std::setiosflags(std::ios::fixed) << std::setprecision(3) << std::setw(12) << imprTypes[a]->Vn 
                << std::setiosflags(std::ios::fixed) << std::setprecision(3) << std::setw(15) << imprTypes[a]->gamma / DEG2RAD
                << std::setiosflags(std::ios::fixed) << std::setprecision(3) << std::setw(12) << imprTypes[a]->Nt
                << "\n";
      }
    }

    ofrcmod << "\nNONB" << std::endl;
    for (unsigned int a = 0; a < atomTypes.size(); a++) {
      if (atomTypes[a]->groupName == ps) {
        ofrcmod << std::left << std::setw(2) << atomTypes[a]->name 
                << std::right 
                << std::setiosflags(std::ios::fixed) << std::setprecision(2) << std::setw(9) << atomTypes[a]->rvalue
                << std::setiosflags(std::ios::fixed) << std::setprecision(3) << std::setw(9) << atomTypes[a]->evalue
                << "\n";
      }
    }

    ofrcmod.close();
}

} // MTKpp namespace

