/*!
   \file stdLibParser.cpp
   \brief Parses standard library xml files using XERCES-C and Trolltech's Qt
   \author Martin Peters

   $Date: 2010/08/19 13:48:48 $
   $Revision: 1.26 $

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
#include "stdLibParser.h"
#include "mtkppParser.h"

#include "StringManip.h"

#include "Molecule/collection.h"
#include "Molecule/stdLibrary.h"
#include "Molecule/stdGroup.h"
#include "Molecule/stdFrag.h"
#include "Molecule/parameters.h"

#include "Diagnostics/MTKException.h"
#include "Log/errorHandler.h"

#include "config.h"

#ifdef USE_XERCES
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#endif // USE_XERCES

namespace MTKpp
{

// ============================================================
// Function : stdLibParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
stdLibParser::stdLibParser(collection* c, stdLibrary* s, parameters* p)
{
    setError(0);

    if (c and s and p) {
      pCollection = c;
      pStdLibrary = s;
      pParameters = p;
    }
    else {
      setError(1);
      std::string errorMessage = "stdLibParser::Error, can't find collection or parameters or stdLibrary objects";
      setErrorMessage(errorMessage);
    }
}

// ============================================================
// Function : stdLibParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
stdLibParser::stdLibParser(stdLibrary* s, parameters* p)
{
    setError(0);

    if (s and p) {
      pStdLibrary = s;
      pParameters = p;
      pCollection = 0;
    }
    else {
      setError(1);
      std::string errorMessage = "stdLibParser::Error, can't find parameters or stdLibrary objects";
      setErrorMessage(errorMessage);
    }
}

// ============================================================
// Function : stdLibParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
stdLibParser::stdLibParser(stdLibrary* s)
{
    if (s) {
      pStdLibrary = s;
      pParameters = 0;
      pCollection = 0;
    }
    else {
      //errorLogger.throwError("stdLibParser", " Can't find pointer ", ERROR);
      errorLogger.throwError("stdLibParser", " Can't find pointer ", MTK_ERROR);
      //exit(1);
      std::stringstream ss;
      ss <<"stdLibParser"<< " Can't find pointer ";
      throw MTKException(ss.str());
    }
    setError(0);
}

// ============================================================
// Function : stdLibParser()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
stdLibParser::~stdLibParser() {}

    // ---------------------------- //
    // -    QT READ FUNCTIONS     - //
    // ---------------------------- //

#ifdef USE_TINYXML
// ============================================================
// Function : Read
// ------------------------------------------------------------
// Parses standard library xml files using tinyxml
// ------------------------------------------------------------
void stdLibParser::Read(std::string fileName)
{
    std::string errMessage = " Reading " + fileName;
    errorLogger.throwError("stdLibParser", errMessage, INFO);

    TiXmlDocument doc(fileName);
    bool loadOkay = doc.LoadFile();

    if (!loadOkay) {
      //std::cout << doc.ErrorDesc();
      errMessage = " Error loading " + fileName;
      errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);

      std::stringstream ss;
      ss <<"stdLibParser"<< " Can't read " << fileName;
      throw MTKException(ss.str());
    }
    else {
      TiXmlNode* node = 0;
      TiXmlElement* element;

      node = doc.RootElement();
      element = node->ToElement();

      //TiXmlNode* groupsNode = node->FirstChild("group");

      for (TiXmlNode* groupNode = node->FirstChild("group"); groupNode; groupNode = groupNode->NextSibling("group")) {
        TiXmlElement* groupElement = groupNode->ToElement();
 
        pStdGroup = pStdLibrary->addStdGroup();
        if (groupElement->Attribute("identity")) {
          //std::cout << "    group " << groupElement->Attribute("identity") << std::endl;
          pStdGroup->setName(groupElement->Attribute("identity"));
        }

        if (groupElement->Attribute("info")) {
          //std::cout << "    group info " << groupElement->Attribute("info") << std::endl;
          pStdGroup->setInfo(groupElement->Attribute("info"));
        }

        for (TiXmlNode* molNode = groupNode->FirstChild("molecule"); molNode; molNode = molNode->NextSibling("molecule")) {
          if (pCollection) {
            mtkppParser* pMTKppParser = new mtkppParser;
            molecule* pStdMolecule = pCollection->addMolecule();
            pMTKppParser->readMolecule(molNode, pStdMolecule);
            pStdGroup->setStdMolecule(pStdMolecule);
          }
          else {
            std::stringstream ss;
            ss << "stdLibParser" << " Can't find pointer ";
            throw MTKException(ss.str());
          }
        }

        for (TiXmlNode* fragmentNode = groupNode->FirstChild("fragment"); fragmentNode; fragmentNode = fragmentNode->NextSibling("fragment")) {
          TiXmlElement* fragmentElement = fragmentNode->ToElement();
          pStdFrag = pStdGroup->addStdFrag();
          std::string att;
          std::string lll = "";
          std::string  l = "";

          if (fragmentElement->Attribute("identity")) {
            //std::cout << "    fragment " << fragmentElement->Attribute("identity") << std::endl;
            pStdFrag->setName(fragmentElement->Attribute("identity") );
          }

          if (fragmentElement->Attribute("symbol")) {
            lll = fragmentElement->Attribute("symbol");
            pStdFrag->setSymbol(lll);
          }

          if (fragmentElement->Attribute("character")) {
            l = fragmentElement->Attribute("character");
            pStdFrag->setCharacter(l);
          }

          if (fragmentElement->Attribute("code")) {
            pStdFrag->setCode(fragmentElement->Attribute("code"));
          }

          if (fragmentElement->Attribute("type")) {
            pStdFrag->setType(fragmentElement->Attribute("type"));
          }

          if (fragmentElement->Attribute("symmetry")) {
            pStdFrag->setSymmetry(fragmentElement->Attribute("symmetry"));
          }

          if (fragmentElement->Attribute("subGraphs")) {
            att = fragmentElement->Attribute("subGraphs");
            std::string sGs = att;
            std::vector<std::string> vsCodes;
            splitString(sGs, " ", vsCodes, 0);
            pStdFrag->setSubGraphs(vsCodes);
          }

          if (lll != "" and l != "") {
            pStdLibrary->setL(lll, l);
          }

          for (TiXmlNode* atomNode = fragmentNode->FirstChild("atom"); atomNode; atomNode = atomNode->NextSibling("atom")) {
            TiXmlElement* atomElement = atomNode->ToElement();

            pStdAtom = pStdFrag->addStdAtom();
            int bType = 1;
            int bKind = 0;
            int bTop = 0;
            std::string errMessage;
            std::string att;

            if (atomElement->Attribute("identity")) {
              //std::cout << "      atom |" << atomElement->Attribute("identity") << "|" << std::endl;
              att = atomElement->Attribute("identity");
              if (att.size() != 4) {
                errMessage = "|" + att + "|\n";
                errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
                //exit(1);
                std::stringstream ss;
                ss <<"stdLibParser"<< " Can't find pointer ";
                throw MTKException(ss.str());
              }
              pStdAtom->identity = att;
            }

            if (atomElement->Attribute("index")) {
              att = atomElement->Attribute("index");
              pStdAtom->index = string2Int(att);
            }

            if (atomElement->Attribute("type")) {
              pStdAtom->type = atomElement->Attribute("type");
            }

            if (atomElement->Attribute("chain")) {
              pStdAtom->chain = atomElement->Attribute("chain");
            }

            if (atomElement->Attribute("atmCharge")) {
              att = atomElement->Attribute("atmCharge");
              pStdAtom->atmCharge = string2Double(att);
            }

            if (atomElement->Attribute("bond12")) {
              att = atomElement->Attribute("bond12");
              pStdAtom->bond12 = string2Int(att);
            }

            if (atomElement->Attribute("bType")) {
              att = atomElement->Attribute("bType");
              bType = string2Int(att);
            }

            if (atomElement->Attribute("bKind")) {
              att = atomElement->Attribute("bKind");
              bKind = string2Int(att);
            }

            if (atomElement->Attribute("bTop")) {
              att = atomElement->Attribute("bTop");
              bTop = string2Int(att);
            }

            if (atomElement->Attribute("bond13")) {
              att = atomElement->Attribute("bond13");
              pStdAtom->bond13 = string2Int(att);
            }

            if (atomElement->Attribute("bond14")) {
              att = atomElement->Attribute("bond14");
              pStdAtom->bond14 = string2Int(att);
            }

            if (atomElement->Attribute("bondLength")) {
              att = atomElement->Attribute("bondLength");
              pStdAtom->bondLength = string2Double(att);
            }

            if (atomElement->Attribute("bondAngle")) {
              att = atomElement->Attribute("bondAngle");
              pStdAtom->bondAngle = string2Double(att);
            }

            if (atomElement->Attribute("bondTorsion")) {
              att = atomElement->Attribute("bondTorsion");
              pStdAtom->bondTorsion = string2Double(att);
            }

            if (atomElement->Attribute("kind")) {
              att = atomElement->Attribute("kind");
              pStdAtom->kind = string2Int(att);
            }

            if (pStdAtom->type != "" and pParameters) {
              std::string atSymbol = pParameters->getAtomTypeSymbol(pStdAtom->type);
              if (atSymbol == "") {
                errMessage = " Cannot find element symbol for " + pStdAtom->type;
                errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);

                std::stringstream ss;
                ss << errMessage;
                throw MTKException(ss.str());
              }
              pStdAtom->atSymbol = atSymbol;
              pStdAtom->atNum = pParameters->getAtomicNum(pStdAtom->type);
            }
            else {
              errMessage = " Cannot find element symbol because the standard atom does not have a type or the parameters were not loaded ";
              errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);

              std::stringstream ss;
              ss << errMessage;
              throw MTKException(ss.str());
            }

            if (pStdAtom->index && pStdAtom->bondLength) {
              pStdBond = pStdFrag->addStdBond();
              pStdBond->atom1 = pStdAtom->index;
              pStdBond->atom2 = pStdAtom->bond12;
              pStdBond->type = bType;
              pStdBond->kind = bKind;
              pStdBond->topology = bTop;
              pStdBond->stereo = 0;
              bType = 1;
              bKind = 0;
              bTop = 0;
              pStdBond->length = pStdAtom->bondLength;
            }
          }

          for (TiXmlNode* aliasNode = fragmentNode->FirstChild("alias"); aliasNode; aliasNode = aliasNode->NextSibling("alias")) {
            TiXmlElement* aliasElement = aliasNode->ToElement();
            std::string errMessage;
            std::string strAlias;
            std::string strOriginal;
            std::string att;

            if (aliasElement->Attribute("alias") and aliasElement->Attribute("original")) {
              strAlias = aliasElement->Attribute("alias");
              strOriginal = aliasElement->Attribute("original");
              if ((strAlias != "") and (strOriginal != "")) {
                pStdAlias = pStdFrag->addStdAlias();
                pStdAlias->atom1 = strOriginal;
                pStdAlias->atom2 = strAlias;
              }
            }
          }

          for (TiXmlNode* loopNode = fragmentNode->FirstChild("loop"); loopNode; loopNode = loopNode->NextSibling("loop")) {
            TiXmlElement* loopElement = loopNode->ToElement();

            std::string errMessage;
            std::string att;

            int iAtom1 = 0;
            int iAtom2 = 0;
            int bType = 1;

            if (loopElement->Attribute("bType")) {
              bType = string2Int(loopElement->Attribute("bType"));
            }

            if (loopElement->Attribute("Atom1") and loopElement->Attribute("Atom2")) {
              iAtom1 = string2Int(loopElement->Attribute("Atom1"));
              iAtom2 = string2Int(loopElement->Attribute("Atom2"));

              pStdLoop = pStdFrag->addStdLoop();
              pStdLoop->atom1 = iAtom1;
              pStdLoop->atom2 = iAtom2;
              pStdLoop->type = bType;
              bType = 1;
            }
          }

          for (TiXmlNode* impNode = fragmentNode->FirstChild("improper"); impNode; impNode = impNode->NextSibling("improper")) {
            TiXmlElement* impElement = impNode->ToElement();

            std::string errMessage;

            std::string strAtom1 = "";
            std::string strAtom2 = "";
            std::string strAtom3 = "";
            std::string strAtom4 = "";
            std::string att;

            if (impElement->Attribute("Atom1") and impElement->Attribute("Atom2") and
                impElement->Attribute("Atom3") and impElement->Attribute("Atom4")) {
              strAtom1 = impElement->Attribute("Atom1");
              strAtom2 = impElement->Attribute("Atom2");
              strAtom3 = impElement->Attribute("Atom3");
              strAtom4 = impElement->Attribute("Atom4");

              if ((strAtom1 != "") and (strAtom2 != "") and
                  (strAtom3 != "") and (strAtom4 != "")) {
                pStdImproper = pStdFrag->addStdImproper();
                pStdImproper->atom1 = string2Int(strAtom1);
                pStdImproper->atom2 = string2Int(strAtom2);
                pStdImproper->atom3 = string2Int(strAtom3);
                pStdImproper->atom4 = string2Int(strAtom4);
              }
            }
          }

          for (TiXmlNode* ringNode = fragmentNode->FirstChild("ring"); ringNode; ringNode = ringNode->NextSibling("ring")) {
            TiXmlElement* ringElement = ringNode->ToElement();

            std::string errMessage;
            std::string att;

            int rSize = 0;
            std::string rAtoms = "";
            int rPlanar = 0;
            int rAromatic = 0;
            int rHetero = 0;
            int rNHetero = 0;
            int rNNitrogen = 0;
            int rNOxygen = 0;
            int rNSulfur = 0;

            if (ringElement->Attribute("size")) {
              rSize = string2Int(ringElement->Attribute("size"));
            }

            if (ringElement->Attribute("atoms")) {
              rAtoms = ringElement->Attribute("atoms");
            }

            if (ringElement->Attribute("planar")) {
              rPlanar = string2Int(ringElement->Attribute("planar"));
            }

            if (ringElement->Attribute("aromatic")) {
              rAromatic = string2Int(ringElement->Attribute("aromatic"));
            }

            if (ringElement->Attribute("hetero")) {
              rHetero = string2Int(ringElement->Attribute("hetero"));
            }

            if (ringElement->Attribute("nHetero")) {
              rNHetero = string2Int(ringElement->Attribute("nHetero"));
            }

            if (ringElement->Attribute("nHetero")) {
              rNHetero = string2Int(ringElement->Attribute("nHetero"));
            }

            if (ringElement->Attribute("nNitrogen")) {
              rNNitrogen = string2Int(ringElement->Attribute("nNitrogen"));
            }

            if (ringElement->Attribute("nOxygen")) {
              rNOxygen = string2Int(ringElement->Attribute("nOxygen"));
            }

            if (ringElement->Attribute("nSulfur")) {
              rNSulfur = string2Int(ringElement->Attribute("nSulfur"));
            }

            if (rAtoms != "") {
              pStdRing = pStdFrag->addStdRing();
              std::vector<std::string> vsAtoms;
              splitString(rAtoms, " ", vsAtoms, 0);
              for (unsigned int x = 0; x < vsAtoms.size(); x++) {
                pStdRing->atoms.push_back(atoi(vsAtoms[x].c_str()));
              }

              if (rSize) {
                if (rSize == static_cast<int>(vsAtoms.size())) {
                  pStdRing->size = rSize;
                }
                else {
                  // put error here
                }
              }

              pStdRing->planar = rPlanar;
              pStdRing->aromatic = rAromatic;
              pStdRing->hetero = rHetero;
              pStdRing->nHetero = rNHetero;
              pStdRing->nNitrogen = rNNitrogen;
              pStdRing->nOxygen = rNOxygen;
              pStdRing->nSulfur = rNSulfur;
            }
          }

          for (TiXmlNode* ftNode = fragmentNode->FirstChild("feature"); ftNode; ftNode = ftNode->NextSibling("feature")) {
            TiXmlElement* ftElement = ftNode->ToElement();

            std::string strAtI = "";
            std::string strAtoms = "";

            if (ftElement->Attribute("identity")) {
              strAtI = ftElement->Attribute("identity");
            }

            if (ftElement->Attribute("atoms")) {
              strAtoms = ftElement->Attribute("atoms");
            }

            if (strAtI != "" and strAtoms != "") {
              pStdFeature = pStdFrag->addStdFeature();
              pStdFeature->name = strAtI;
              std::vector<std::string> vsAtoms;
              splitString(strAtoms, " ", vsAtoms, 0);
              for (unsigned int x = 0; x < vsAtoms.size(); x++) {
                pStdFeature->atoms.push_back(atoi(vsAtoms[x].c_str()));
              }
            }
          }

          for (TiXmlNode* fgNode = fragmentNode->FirstChild("funcGroup"); fgNode; fgNode = fgNode->NextSibling("funcGroup")) {
            TiXmlElement* fgElement = fgNode->ToElement();

            std::string strGroup = "";
            std::string strFrag = "";
            std::string strAtoms = "";
            std::string att = "";

            if (fgElement->Attribute("group")) {
              strGroup = fgElement->Attribute("group");
            }

            if (fgElement->Attribute("fragment")) {
              strFrag = fgElement->Attribute("fragment");
            }

            if (fgElement->Attribute("atoms")) {
              strAtoms = fgElement->Attribute("atoms");
            }

            if (strGroup != "" and strFrag != "" and strAtoms != "") {
              pStdFuncGroup = pStdFrag->addStdFuncGroup();
              pStdFuncGroup->groupName = strGroup;
              pStdFuncGroup->fragName = strFrag;
              std::vector<std::string> vsAtoms;
              splitString(strAtoms, " ", vsAtoms, 0);
              for (unsigned int x = 0; x < vsAtoms.size(); x++) {
                pStdFuncGroup->atoms.push_back(atoi(vsAtoms[x].c_str()));
              }
            }
          }

          for (TiXmlNode* conPtsNode = fragmentNode->FirstChild("connectionPoints"); conPtsNode; conPtsNode = conPtsNode->NextSibling("connectionPoints")) {
            TiXmlElement* conPtsElement = conPtsNode->ToElement();

            std::string strAtoms = "";

            if (conPtsElement->Attribute("atoms")) {
              strAtoms = conPtsElement->Attribute("atoms");
            }

            if (strAtoms != "") {
              std::vector<std::string> vsAtoms;
              std::vector<int> ivsAtoms;
              splitString(strAtoms, " ", vsAtoms, 0);
              for (unsigned int x = 0; x < vsAtoms.size(); x++) {
                ivsAtoms.push_back(atoi(vsAtoms[x].c_str()));
              }
              pStdFrag->addStdConnPts(ivsAtoms);
            }
          }

          for (TiXmlNode* conTorNode = fragmentNode->FirstChild("connectionTorsion"); conTorNode; conTorNode = conTorNode->NextSibling("connectionTorsion")) {
            TiXmlElement* conTorElement = conTorNode->ToElement();
            std::string errMessage;
            int bd = -1;
            int ag = -1;
            int tr = -1;
            double tor = 0.0;

            if (conTorElement->Attribute("bd")) {
              bd = string2Int(conTorElement->Attribute("bd"));
            }

            if (conTorElement->Attribute("ag")) {
              ag = string2Int(conTorElement->Attribute("ag"));
            }

            if (conTorElement->Attribute("tr")) {
              tr = string2Int(conTorElement->Attribute("tr"));
            }

            if (conTorElement->Attribute("value")) {
              tor = string2Double(conTorElement->Attribute("value"));
            }

            if (bd > -1 and ag > -1 and tr > -1) {
              stdConnTorsion* pStdConnTorsion = pStdFrag->addStdConnTorsion();
              pStdConnTorsion->bondAtom = bd;
              pStdConnTorsion->angleAtom = ag;
              pStdConnTorsion->torsionAtom = tr;
              pStdConnTorsion->torsion = tor;
            }
          }

          for (TiXmlNode* rotBondNode = fragmentNode->FirstChild("rotatableBond"); rotBondNode; rotBondNode = rotBondNode->NextSibling("rotatableBond")) {
            TiXmlElement* rotBondElement = rotBondNode->ToElement();

            std::string errMessage;
            int atom1 = -1;
            int atom2 = -1;
            int atom3 = -1;
            int atom4 = -1;
            std::string values = "";

            if (rotBondElement->Attribute("atom1")) {
              atom1 = string2Int(rotBondElement->Attribute("atom1"));
            }

            if (rotBondElement->Attribute("atom2")) {
              atom2 = string2Int(rotBondElement->Attribute("atom2"));
            }

            if (rotBondElement->Attribute("atom3")) {
              atom3 = string2Int(rotBondElement->Attribute("atom3"));
            }

            if (rotBondElement->Attribute("atom4")) {
              atom4 = string2Int(rotBondElement->Attribute("atom4"));
            }

            if (rotBondElement->Attribute("values")) {
              values = string2Int(rotBondElement->Attribute("values"));
            }

            if (atom1 > -1 and atom2 > -1 and atom3 > -1 and atom4 > -1) {
              stdRotBond* pStdRotBond = pStdFrag->addStdRotBond();
              pStdRotBond->atom1 = atom1;
              pStdRotBond->atom2 = atom2;
              pStdRotBond->atom3 = atom3;
              pStdRotBond->atom4 = atom4;
              if (values != "") {
                std::vector<std::string> vsValues;
                splitString(values, " ", vsValues, 0);
                for (unsigned int x = 0; x < vsValues.size(); x++) {
                  pStdRotBond->values.push_back(strtod(vsValues[x].c_str(),0));
                }
              }
              else {
                // error
              }
            }
          }
        }
      }
    }
}
#endif

#ifdef USE_QT
// ============================================================
// Function : Read
// ------------------------------------------------------------
// Parses standard library xml files using Qt
// ------------------------------------------------------------
void stdLibParser::Read(std::string fileName)
{
    QString qFileName = QString::fromStdString(fileName);
    std::string errMessage = " Reading " + fileName;
    errorLogger.throwError("stdLibParser", errMessage, INFO);

    if (!pParameters or !pStdLibrary) {
      errMessage = " Can't read lib file:  " + fileName;
      errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
      return;
    }

    // Read elements file using DOM into memory
    QString mtkppInfo = PACKAGE_TARNAME;
    QDomDocument doc(mtkppInfo);

    QFile file(qFileName);
    if (!file.open(QIODevice::ReadOnly)) {
      errMessage = " Reading " + fileName;
      errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
    }

    QString domErrorMessage = "";
    int domErrorLine = 0;
    int domErrorColumn = 0;

    if (!doc.setContent(&file, &domErrorMessage, &domErrorLine, &domErrorColumn)) {
      errMessage = " Reading " + fileName + " " + domErrorMessage.toStdString() + " " + int2String(domErrorLine) + " " + int2String(domErrorColumn);
      errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
      file.close();
     // return 1;
    }
    file.close();

    // Parse the contents of the elements file
    QString att;
    QDomElement docElem = doc.documentElement();
    if (docElem.tagName() == "stdLib") {

      QDomNode groupNode = docElem.firstChild();
      while (!groupNode.isNull()) {
        if (groupNode.nodeName() == "group" and groupNode.hasAttributes()) {
          groupFiller(groupNode);
        }
        groupNode = groupNode.nextSibling();
      }
    }
}

// ============================================================
// Function : groupfiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void stdLibParser::groupFiller(QDomNode g)
{
    QDomElement groupElement = g.toElement();

    if (!groupElement.isNull()) {
      QString att;
      pStdGroup = pStdLibrary->addStdGroup();
      if (groupElement.hasAttribute("identity")) {
        att = groupElement.attribute("identity");
        pStdGroup->setName(att.toStdString());
      }

      if (groupElement.hasAttribute("info")) {
        att = groupElement.attribute("info");
        pStdGroup->setInfo(att.toStdString());
      }

      QDomNode curNode = g.firstChild();

      while (!curNode.isNull()) {
        if (curNode.nodeName() == "fragment" and curNode.hasAttributes()) {
          fragmentFiller(curNode);
        }
        if (curNode.nodeName() == "molecule" and curNode.hasAttributes()) {
          if (pCollection) {
            mtkppParser* pMTKppParser = new mtkppParser;
            molecule* pStdMolecule = pCollection->addMolecule();
            pMTKppParser->readMolecule(curNode, pStdMolecule);
            pStdGroup->setStdMolecule(pStdMolecule);
          }
        }
        curNode = curNode.nextSibling();
      }
    }
}

// ============================================================
// Function : fragmentfiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void stdLibParser::fragmentFiller(QDomNode f)
{
    QDomElement fragmentElement = f.toElement();

    if (!fragmentElement.isNull()) {
      QString att;
      QString lll = "";
      QString l = "";
      pStdFrag = pStdGroup->addStdFrag();

      if (fragmentElement.hasAttribute("identity")) {
        att = fragmentElement.attribute("identity");
        pStdFrag->setName(att.toStdString());
      }

      if (fragmentElement.hasAttribute("symbol")) {
        lll = fragmentElement.attribute("symbol");
        pStdFrag->setSymbol(lll.toStdString());
      }

      if (fragmentElement.hasAttribute("character")) {
        l = fragmentElement.attribute("character");
        pStdFrag->setCharacter(l.toStdString());
      }

      if (fragmentElement.hasAttribute("code")) {
        att = fragmentElement.attribute("code");
        pStdFrag->setCode(att.toStdString());
      }

      if (fragmentElement.hasAttribute("type")) {
        att = fragmentElement.attribute("type");
        pStdFrag->setType(att.toStdString());
      }

      if (fragmentElement.hasAttribute("symmetry")) {
        att = fragmentElement.attribute("symmetry");
        pStdFrag->setSymmetry(att.toStdString());
      }

      if (fragmentElement.hasAttribute("subGraphs")) {
        att = fragmentElement.attribute("subGraphs");
        std::string sGs = att.toStdString();
        std::vector<std::string> vsCodes;
        splitString(sGs, " ", vsCodes, 0);
        pStdFrag->setSubGraphs(vsCodes);
      }

      QDomNode node = f.firstChild();
      while (!node.isNull()) {
        if (node.nodeName() == "atom" and node.hasAttributes()) {
          atomFiller(node);
        }
        if (node.nodeName() == "alias" and node.hasAttributes()) {
          aliasFiller(node);
        }
        if (node.nodeName() == "improper" and node.hasAttributes()) {
          improperFiller(node);
        }
        if (node.nodeName() == "loop" and node.hasAttributes()) {
          loopFiller(node);
        }
        if (node.nodeName() == "ring" and node.hasAttributes()) {
          ringFiller(node);
        }
        if (node.nodeName() == "feature" and node.hasAttributes()) {
          featureFiller(node);
        }
        if (node.nodeName() == "funcGroup" and node.hasAttributes()) {
          funcGroupFiller(node);
        }
        if (node.nodeName() == "connectionPoints" and node.hasAttributes()) {
          connPtsFiller(node);
        }
        if (node.nodeName() == "connectionTorsion" and node.hasAttributes()) {
          connTorFiller(node);
        }
        if (node.nodeName() == "rotatableBond" and node.hasAttributes()) {
          rotBondFiller(node);
        }
        node = node.nextSibling();
      }

      if (lll != "" and l != "") {
        pStdLibrary->setL(lll.toStdString(), l.toStdString());
      }
    }
}

// ============================================================
// Function : atomFiller
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::atomFiller(QDomNode a)
{
    QDomElement atomElement = a.toElement();
    std::string errMessage;

    int bType = 1;
    int bKind = 0;
    int bTop = 0;

    if (!atomElement.isNull()) {
      QString att;
      bool bConvertOk = true;
      pStdAtom = pStdFrag->addStdAtom();

      if (atomElement.hasAttribute("identity")) {
        att = atomElement.attribute("identity");

        if (att.toStdString().size() != 4) {
          errMessage = "|" + att.toStdString() + "|\n";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss <<"stdLibParser"<<errMessage;
          throw MTKException(ss.str());
        }
        pStdAtom->identity = att.toStdString();
        //std::cout << pStdFrag->getSymbol() << "@|" << pStdAtom->identity  << "|" << std::endl;
      }

      if (atomElement.hasAttribute("index")) {
        att = atomElement.attribute("index");
        int index = att.toInt(&bConvertOk, 10);
        if (bConvertOk) {
          pStdAtom->index = index;
        }
        else {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (atomElement.hasAttribute("type")) {
        att = atomElement.attribute("type");
        pStdAtom->type = att.toStdString();
      }

      if (atomElement.hasAttribute("chain")) {
        att = atomElement.attribute("chain");
        pStdAtom->chain = att.toStdString();
      }

      if (atomElement.hasAttribute("atmCharge")) {
        att = atomElement.attribute("atmCharge");
        double atmCharge = att.toDouble(&bConvertOk);
        if (bConvertOk) {
          pStdAtom->atmCharge = atmCharge; // * 18.2223 ; // electron to kcal
        }
        else {
          errMessage = " Error converting string to double ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (atomElement.hasAttribute("bond12")) {
        att = atomElement.attribute("bond12");
        int bond12 = att.toInt(&bConvertOk, 10);
        if (bConvertOk) {
          pStdAtom->bond12 = bond12;
        }
        else {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (atomElement.hasAttribute("bType")) {
        att = atomElement.attribute("bType");
        bType = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (atomElement.hasAttribute("bKind")) {
        att = atomElement.attribute("bKind");
        bKind = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (atomElement.hasAttribute("bTop")) {
        att = atomElement.attribute("bTop");
        bTop = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (atomElement.hasAttribute("bond13")) {
        att = atomElement.attribute("bond13");
        int bond13 = att.toInt(&bConvertOk, 10);
        if (bConvertOk) {
          pStdAtom->bond13 = bond13;
        }
        else {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (atomElement.hasAttribute("bond14")) {
        att = atomElement.attribute("bond14");
        int bond14 = att.toInt(&bConvertOk, 10);
        if (bConvertOk) {
          pStdAtom->bond14 = bond14;
        }
        else {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (atomElement.hasAttribute("bondLength")) {
        att = atomElement.attribute("bondLength");
        double bondLength = att.toDouble(&bConvertOk);
        if (bConvertOk) {
          pStdAtom->bondLength = bondLength;
        }
        else {
          errMessage = " Error converting string to double ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (atomElement.hasAttribute("bondAngle")) {
        att = atomElement.attribute("bondAngle");
        double bondAngle = att.toDouble(&bConvertOk);
        if (bConvertOk) {
          pStdAtom->bondAngle = bondAngle;
        }
        else {
          errMessage = " Error converting string to double ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (atomElement.hasAttribute("bondTorsion")) {
        att = atomElement.attribute("bondTorsion");
        double bondTorsion = att.toDouble(&bConvertOk);
        if (bConvertOk) {
          pStdAtom->bondTorsion = bondTorsion;
        }
        else {
          errMessage = " Error converting string to double ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (atomElement.hasAttribute("kind")) {
        att = atomElement.attribute("kind");
        int kind = att.toInt(&bConvertOk, 10);
        if (bConvertOk) {
          pStdAtom->kind = kind;
        }
        else {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (pStdAtom->type != "" and pParameters) {
        std::string atSymbol = pParameters->getAtomTypeSymbol(pStdAtom->type);
        if (atSymbol == "") {
          errMessage = " Cannot find element symbol for " + pStdAtom->type;
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
        pStdAtom->atSymbol = atSymbol;
        pStdAtom->atNum = pParameters->getAtomicNum(pStdAtom->type);
      }
      else {
        errMessage = " Cannot find element symbol ";
        errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
        //exit(1);

        std::stringstream ss;
        ss << "stdLibParser" << errMessage;
        throw MTKException(ss.str());
      }

      if (pStdAtom->index && pStdAtom->bondLength) {
        pStdBond = pStdFrag->addStdBond();
        pStdBond->atom1 = pStdAtom->index;
        pStdBond->atom2 = pStdAtom->bond12;
        pStdBond->type = bType;
        pStdBond->kind = bKind;
        pStdBond->topology = bTop;
        pStdBond->stereo = 0;
        bType = 1;
        bKind = 0;
        bTop = 0;
        pStdBond->length = pStdAtom->bondLength;
      }
    }
}

// ============================================================
// Function : aliasFiller
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::aliasFiller(QDomNode a)
{
    QDomElement aliasElement = a.toElement();
    std::string errMessage;

    std::string strAlias;
    std::string strOriginal;

    if (!aliasElement.isNull()) {
      QString att;

      if (aliasElement.hasAttribute("alias")) {
        att = aliasElement.attribute("alias");
        strAlias = att.toStdString();
      }

      if (aliasElement.hasAttribute("original")) {
        att = aliasElement.attribute("original");
        strOriginal = att.toStdString();
      }
    }

    if ((strAlias != "") and (strOriginal != "")) {
      pStdAlias = pStdFrag->addStdAlias();
      pStdAlias->atom1 = strOriginal;
      pStdAlias->atom2 = strAlias;
    }
}

// ============================================================
// Function : improperFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void stdLibParser::improperFiller(QDomNode a)
{
    QDomElement impElement = a.toElement();
    std::string errMessage;

    std::string strAtom1 = "";
    std::string strAtom2 = "";
    std::string strAtom3 = "";
    std::string strAtom4 = "";

    if (!impElement.isNull()) {
      QString att;

      if (impElement.hasAttribute("Atom1")) {
        att = impElement.attribute("Atom1");
        strAtom1 = att.toStdString();
      }

      if (impElement.hasAttribute("Atom2")) {
        att = impElement.attribute("Atom2");
        strAtom2 = att.toStdString();
      }

      if (impElement.hasAttribute("Atom3")) {
        att = impElement.attribute("Atom3");
        strAtom3 = att.toStdString();
      }

      if (impElement.hasAttribute("Atom4")) {
        att = impElement.attribute("Atom4");
        strAtom4 = att.toStdString();
      }
    }

    if ((strAtom1 != "") and (strAtom2 != "") and
        (strAtom3 != "") and (strAtom4 != "")) {
      pStdImproper = pStdFrag->addStdImproper();

/*
      if (std::string(strAtom1) == "p1") {
        pStdImproper->atom1 = -1;
      }
      else if (std::string(strAtom1) == "n1") {
        pStdImproper->atom1 = -4;
      }
      else {
        pStdImproper->atom1 = atoi(strAtom1.c_str());
      }
*/

      pStdImproper->atom1 = atoi(strAtom1.c_str());
      pStdImproper->atom2 = atoi(strAtom2.c_str());
      pStdImproper->atom3 = atoi(strAtom3.c_str());
      pStdImproper->atom4 = atoi(strAtom4.c_str());
    }
}

// ============================================================
// Function : loopFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void stdLibParser::loopFiller(QDomNode l)
{
    QDomElement loopElement = l.toElement();
    std::string errMessage;

    int iAtom1 = 0;
    int iAtom2 = 0;
    int bType = 1;
    bool bConvertOk = true;

    if (!loopElement.isNull()) {
      QString att;

      if (loopElement.hasAttribute("Atom1")) {
        att = loopElement.attribute("Atom1");
        iAtom1 = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (loopElement.hasAttribute("Atom2")) {
        att = loopElement.attribute("Atom2");
        iAtom2 = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (loopElement.hasAttribute("bType")) {
        att = loopElement.attribute("bType");
        bType = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }
    }
    else {
      errMessage = " Error reading loop in standard library file ";
      errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
      //exit(1);

      std::stringstream ss;
      ss << "stdLibParser" << errMessage;
      throw MTKException(ss.str());
    }
    pStdLoop = pStdFrag->addStdLoop();
    pStdLoop->atom1 = iAtom1;
    pStdLoop->atom2 = iAtom2;
    pStdLoop->type = bType;
    bType = 1;
}

// ============================================================
// Function : ringFiller
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void stdLibParser::ringFiller(QDomNode r)
{
    QDomElement ringElement = r.toElement();
    std::string errMessage;

    int rSize = 0;
    std::string rAtoms = "";
    int rPlanar = 0;
    int rAromatic = 0;
    int rHetero = 0;
    int rNHetero = 0;
    int rNNitrogen = 0;
    int rNOxygen = 0;
    int rNSulfur = 0;

    bool bConvertOk = true;

    if (!ringElement.isNull()) {
      QString att;

      if (ringElement.hasAttribute("size")) {
        att = ringElement.attribute("size");
        rSize = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (ringElement.hasAttribute("atoms")) {
        att = ringElement.attribute("atoms");
        rAtoms = att.toStdString();
      }

      if (ringElement.hasAttribute("planar")) {
        att = ringElement.attribute("planar");
        rPlanar = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (ringElement.hasAttribute("aromatic")) {
        att = ringElement.attribute("aromatic");
        rAromatic = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (ringElement.hasAttribute("hetero")) {
        att = ringElement.attribute("hetero");
        rHetero = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (ringElement.hasAttribute("nHetero")) {
        att = ringElement.attribute("nHetero");
        rNHetero = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (ringElement.hasAttribute("nNitrogen")) {
        att = ringElement.attribute("nNitrogen");
        rNNitrogen = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (ringElement.hasAttribute("nOxygen")) {
        att = ringElement.attribute("nOxygen");
        rNOxygen = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (ringElement.hasAttribute("nSulfur")) {
        att = ringElement.attribute("nSulfur");
        rNSulfur = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (rAtoms != "") {
        pStdRing = pStdFrag->addStdRing();
        std::vector<std::string> vsAtoms;
        splitString(rAtoms, " ", vsAtoms, 0);
        for (unsigned int x = 0; x < vsAtoms.size(); x++) {
          pStdRing->atoms.push_back(atoi(vsAtoms[x].c_str()));
        }

        if (rSize) {
          if (rSize == static_cast<int>(vsAtoms.size())) {
            pStdRing->size = rSize;
          }
          else {
            // put error here
          }
        }

        pStdRing->planar = rPlanar;
        pStdRing->aromatic = rAromatic;
        pStdRing->hetero = rHetero;
        pStdRing->nHetero = rNHetero;
        pStdRing->nNitrogen = rNNitrogen;
        pStdRing->nOxygen = rNOxygen;
        pStdRing->nSulfur = rNSulfur;
      }
    }
}

// ============================================================
// Function : featureFiller
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::featureFiller(QDomNode f)
{
    QDomElement featureElement = f.toElement();
    std::string errMessage;

    std::string strAtI = "";
    std::string strAtoms = "";

    if (!featureElement.isNull()) {
      QString att;

      if (featureElement.hasAttribute("identity")) {
        att = featureElement.attribute("identity");
        strAtI = att.toStdString();
      }

      if (featureElement.hasAttribute("atoms")) {
        att = featureElement.attribute("atoms");
        strAtoms = att.toStdString();
      }
    }

    if (strAtI != "" and strAtoms != "") {
      pStdFeature = pStdFrag->addStdFeature();
      pStdFeature->name = strAtI;
      std::vector<std::string> vsAtoms;
      splitString(strAtoms, " ", vsAtoms, 0);
      for (unsigned int x = 0; x < vsAtoms.size(); x++) {
        pStdFeature->atoms.push_back(atoi(vsAtoms[x].c_str()));
      }
    }
}

// ============================================================
// Function : funcGroupFiller
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::funcGroupFiller(QDomNode f)
{
    QDomElement funcGroupElement = f.toElement();
    std::string errMessage;

    std::string strGroup = "";
    std::string strFrag = "";
    std::string strAtoms = "";

    if (!funcGroupElement.isNull()) {
      QString att;

      if (funcGroupElement.hasAttribute("group")) {
        att = funcGroupElement.attribute("group");
        strGroup = att.toStdString();
      }

      if (funcGroupElement.hasAttribute("fragment")) {
        att = funcGroupElement.attribute("fragment");
        strFrag = att.toStdString();
      }

      if (funcGroupElement.hasAttribute("atoms")) {
        att = funcGroupElement.attribute("atoms");
        strAtoms = att.toStdString();
      }
    }

    if (strGroup != "" and strFrag != "" and strAtoms != "") {
      pStdFuncGroup = pStdFrag->addStdFuncGroup();
      pStdFuncGroup->groupName = strGroup;
      pStdFuncGroup->fragName = strFrag;
      std::vector<std::string> vsAtoms;
      splitString(strAtoms, " ", vsAtoms, 0);
      for (unsigned int x = 0; x < vsAtoms.size(); x++) {
        pStdFuncGroup->atoms.push_back(atoi(vsAtoms[x].c_str()));
      }
    }
}

// ============================================================
// Function : connPtsFiller
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::connPtsFiller(QDomNode f)
{
    QDomElement connPtsElement = f.toElement();
    std::string errMessage;
    std::string strAtoms = "";

    if (!connPtsElement.isNull()) {
      QString att;

      if (connPtsElement.hasAttribute("atoms")) {
        att = connPtsElement.attribute("atoms");
        strAtoms = att.toStdString();
      }
    }

    if (strAtoms != "") {
      std::vector<std::string> vsAtoms;
      std::vector<int> ivsAtoms;
      splitString(strAtoms, " ", vsAtoms, 0);
      for (unsigned int x = 0; x < vsAtoms.size(); x++) {
        ivsAtoms.push_back(atoi(vsAtoms[x].c_str()));
      }
      pStdFrag->addStdConnPts(ivsAtoms);
    }
}

// ============================================================
// Function : connTorFiller
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::connTorFiller(QDomNode f)
{
    QDomElement connTorElement = f.toElement();
    std::string errMessage;
    int bd = -1;
    int ag = -1;
    int tr = -1;
    double tor = 0.0;
    bool bConvertOk = true;

    if (!connTorElement.isNull()) {
      QString att;

      if (connTorElement.hasAttribute("bd")) {
        att = connTorElement.attribute("bd");
        bd = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (connTorElement.hasAttribute("ag")) {
        att = connTorElement.attribute("ag");
        ag = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (connTorElement.hasAttribute("tr")) {
        att = connTorElement.attribute("tr");
        tr = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (connTorElement.hasAttribute("value")) {
        att = connTorElement.attribute("value");
        tor = att.toDouble(&bConvertOk);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }
    }

    if (bd > -1 and ag > -1 and tr > -1) {
      stdConnTorsion* pStdConnTorsion = pStdFrag->addStdConnTorsion();
      pStdConnTorsion->bondAtom = bd;
      pStdConnTorsion->angleAtom = ag;
      pStdConnTorsion->torsionAtom = tr;
      pStdConnTorsion->torsion = tor;
    }
}

// ============================================================
// Function : rotBondFiller
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::rotBondFiller(QDomNode a)
{
    QDomElement rotBondElement = a.toElement();
    std::string errMessage;
    int atom1 = -1;
    int atom2 = -1;
    int atom3 = -1;
    int atom4 = -1;
    std::string values = "";
    bool bConvertOk = true;

    if (!rotBondElement.isNull()) {
      QString att;

      if (rotBondElement.hasAttribute("atom1")) {
        att = rotBondElement.attribute("atom1");
        atom1 = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (rotBondElement.hasAttribute("atom2")) {
        att = rotBondElement.attribute("atom2");
        atom2 = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (rotBondElement.hasAttribute("atom3")) {
        att = rotBondElement.attribute("atom3");
        atom3 = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (rotBondElement.hasAttribute("atom4")) {
        att = rotBondElement.attribute("atom4");
        atom4 = att.toInt(&bConvertOk, 10);
        if (!bConvertOk) {
          errMessage = " Error converting string to integer ";
          errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
          //exit(1);

          std::stringstream ss;
          ss << "stdLibParser" << errMessage;
          throw MTKException(ss.str());
        }
      }

      if (rotBondElement.hasAttribute("values")) {
        att = rotBondElement.attribute("values");
        values = att.toStdString();
      }
    }

    if (atom1 > -1 and atom2 > -1 and atom3 > -1 and atom4 > -1) {
      stdRotBond* pStdRotBond = pStdFrag->addStdRotBond();
      pStdRotBond->atom1 = atom1;
      pStdRotBond->atom2 = atom2;
      pStdRotBond->atom3 = atom3;
      pStdRotBond->atom4 = atom4;
      if (values != "") {
        std::vector<std::string> vsValues;
        splitString(values, " ", vsValues, 0);
        for (unsigned int x = 0; x < vsValues.size(); x++) {
          pStdRotBond->values.push_back(strtod(vsValues[x].c_str(),0));
        }
      }
      else {
        // error
      }
    }
}

#endif // USE_QT

    // ---------------------------------- //
    // -    XERCES-C READ FUNCTIONS     - //
    // ---------------------------------- //

#ifdef USE_XERCES
// ============================================================
// Function : Read
// ------------------------------------------------------------
// Parses standard library xml files using xercesc
// ------------------------------------------------------------
void stdLibParser::Read(std::string fileName)
{
    std::string errMessage = " Reading " + fileName;
    errorLogger.throwError("stdLibParser", errMessage, INFO);

    errMessage = "";
    setError(0);
    if (!fileExists(fileName)) {
      setError(1);
      std::string errorMessage = "  stdLibParser::Read Error, Can't Find " + fileName;
      setErrorMessage(errorMessage);
      return;
    }

    if (!pParameters or !pStdLibrary) {
      errMessage = " Can't read lib file:  " + fileName;
      errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
      return;
    }

    try {
       XMLPlatformUtils::Initialize();
    }
    catch (const XMLException& toCatch) {
       char* message = XMLString::transcode(toCatch.getMessage());
       XMLString::release(&message);
       setError(1);
       std::string errorMessage = "  stdLibParser::Error during initialization " + std::string(message);
       setErrorMessage(errorMessage);
       return;
    }

    XMLCh tempStr[100];
    XMLString::transcode("LS", tempStr, 99);
    DOMImplementation *impl = DOMImplementationRegistry::getDOMImplementation(tempStr);

#if (XERCES_VERSION_MAJOR == 3)
    DOMLSParser* parser=((DOMImplementationLS*)impl)->createLSParser(DOMImplementationLS::MODE_SYNCHRONOUS,0, XMLPlatformUtils::fgMemoryManager, NULL);
#else
    DOMBuilder* parser=((DOMImplementationLS*)impl)->createDOMBuilder(DOMImplementationLS::MODE_SYNCHRONOUS,0);

    // optionally you can set some features on this builder
    if (parser->canSetFeature(XMLUni::fgDOMValidation, true))
       parser->setFeature(XMLUni::fgDOMValidation, true);
    if (parser->canSetFeature(XMLUni::fgDOMNamespaces, true))
       parser->setFeature(XMLUni::fgDOMNamespaces, true);
    if (parser->canSetFeature(XMLUni::fgDOMDatatypeNormalization, true))
       parser->setFeature(XMLUni::fgDOMDatatypeNormalization, true);
#endif

    char* xmlFile = (char*)(fileName.c_str());
    XERCES_CPP_NAMESPACE::DOMDocument *doc = 0;

    //ErrorHandler* errHandler = (ErrorHandler*) new HandlerBase();
    //parser->setErrorHandler(errHandler);
    //MyDOMErrorHandler* errHandler = new MyDOMErrorHandler();
    //parser->setErrorHandler(errHandler);
/*
    try {
       doc = parser->parseURI(xmlFile);
       stdLibFiller((DOMNode*)doc->getDocumentElement());
    }
    catch (const XMLException& toCatch) {
       char* message = XMLString::transcode(toCatch.getMessage());
       std::cout << "Exception message is: \n"
                 << message << "\n";
       XMLString::release(&message);
       setError(1);
       return;
    }
    catch (const DOMException& toCatch) {
       char* message = XMLString::transcode(toCatch.msg);
       std::cout << "Exception message is: \n"
                 << message << "\n";
       XMLString::release(&message);
       setError(1);
       return;
    }
    catch (...) {
       std::cout << "Unexpected Exception \n" ;
       setError(1);
       return;
    }
*/
    doc = parser->parseURI(xmlFile);
    if (doc == NULL) {
      throw MTKException("stdLibParser: parse effort turned up nothing for "+fileName+".");
    }

    stdLibFiller((DOMNode*)doc->getDocumentElement());

    parser->release();
    XMLPlatformUtils::Terminate();
    setError(0);
    return;
}

// ============================================================
// Function : stdLibFiller
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::stdLibFiller(DOMNode *rootnode)
{
    char *name = XMLString::transcode(rootnode->getNodeName());

    if (std::string(name) == std::string("stdLib")) {
       groupFiller(rootnode);
    }
    delete name;
}

// ============================================================
// Function : groupfiller
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::groupFiller(DOMNode *rootnode)
{
    DOMNode *c;
    for (c = rootnode->getFirstChild(); c != 0; c = c->getNextSibling()) {
      if (c->getNodeType() == 1) {
        if ((std::string)(XC(c->getNodeName())) == "group") {
          pStdGroup = pStdLibrary->addStdGroup();
          if (c->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = c->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("identity")) {
                char *cGroupName = XMLString::transcode(pAttributeNode->getNodeValue());
                pStdGroup->setName(std::string(cGroupName));
                delete cGroupName;
              }
              if (std::string(name) == std::string("info")) {
                char *cGroupInfo = XMLString::transcode(pAttributeNode->getNodeValue());
                pStdGroup->setInfo(std::string(cGroupInfo));
                delete cGroupInfo;
              }
              delete name;
            }
          }
          fragmentFiller(c);
        }
        else {
        }
      }
      else {
      }
    }
}

// ============================================================
// Function : fragmentFiller
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::fragmentFiller(DOMNode *rootnode)
{
    DOMNode *child;
    DOMNode *a;
    for (child = rootnode->getFirstChild(); child != 0; child = child->getNextSibling()) {
      if (child->getNodeType() == 1) {
        if ((std::string)(XC(child->getNodeName())) == "fragment") {
          pStdFrag = pStdGroup->addStdFrag();
          if (child->hasAttributes()) {
            // get all the attributes of the node
            DOMNamedNodeMap *pAttributes = child->getAttributes();
            int nSize = pAttributes->getLength();
            for (int i = 0; i < nSize; ++i) {
              DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
              // get attribute name
              char *name = XMLString::transcode(pAttributeNode->getName());
              if (std::string(name) == std::string("identity")) {
                char *cFragName = XMLString::transcode(pAttributeNode->getNodeValue());
                pStdFrag->setName(std::string(cFragName)); // long name
                delete cFragName;
              }
              if (std::string(name) == std::string("symbol")) {
                char *cSymbol = XMLString::transcode(pAttributeNode->getNodeValue());
                pStdFrag->setSymbol(std::string(cSymbol)); // 3L code
                delete cSymbol;
              }
              if (std::string(name) == std::string("character")) {
                char *cChar = XMLString::transcode(pAttributeNode->getNodeValue());
                pStdFrag->setCharacter(std::string(cChar)); // 1L code
                delete cChar;
              }
              if (std::string(name) == std::string("code")) {
                char *cCode = XMLString::transcode(pAttributeNode->getNodeValue());
                pStdFrag->setCode(std::string(cCode)); // 8L code
                delete cCode;
              }
              if (std::string(name) == std::string("type")) {
                char *cType = XMLString::transcode(pAttributeNode->getNodeValue());
                pStdFrag->setType(std::string(cType));
                delete cType;
              }
              if (std::string(name) == std::string("symmetry")) {
                char *cSymmetry = XMLString::transcode(pAttributeNode->getNodeValue());
                pStdFrag->setSymmetry(std::string(cSymmetry));
                delete cSymmetry;
              }
              if (std::string(name) == std::string("subGraphs")) {
                char *cGs = XMLString::transcode(pAttributeNode->getNodeValue());
                std::string sGs = std::string(cGs);
                std::vector<std::string> vsCodes;
                splitString(sGs, " ", vsCodes, 0);
                pStdFrag->setSubGraphs(vsCodes);
                delete cGs;
              }
              delete name;
            }
          }
          for (a = child->getFirstChild(); a != 0; a = a->getNextSibling()) {
            if (a->getNodeType() == 1) {
              if ((std::string)(XC(a->getNodeName())) == "atom") {
                pStdAtom = pStdFrag->addStdAtom();
                atomFiller(a);
              }
              else if ((std::string)(XC(a->getNodeName())) == "alias") {
                aliasFiller(a);
              }
              else if ((std::string)(XC(a->getNodeName())) == "improper") {
                improperFiller(a);
              }
              else if ((std::string)(XC(a->getNodeName())) == "loop") {
                loopFiller(a);
              }
              else if ((std::string)(XC(a->getNodeName())) == "ring") {
                ringFiller(a);
              }
              else if ((std::string)(XC(a->getNodeName())) == "feature") {
                featureFiller(a);
              }
              else if ((std::string)(XC(a->getNodeName())) == "funcGroup") {
                funcGroupFiller(a);
              }
              else if ((std::string)(XC(a->getNodeName())) == "connectionPoints") {
                connPtsFiller(a);
              }
              else if ((std::string)(XC(a->getNodeName())) == "connectionTorsion") {
                connTorFiller(a);
              }
              else if ((std::string)(XC(a->getNodeName())) == "rotatableBond") {
                rotBondFiller(a);
              }
            }
          }
        }
        else if ((std::string)(XC(child->getNodeName())) == "molecule") {

          if (pCollection) {
            mtkppParser* pMTKppParser = new mtkppParser;
            molecule* pStdMolecule = pCollection->addMolecule();
            pMTKppParser->readMolecule(child, pStdMolecule);
            pStdGroup->setStdMolecule(pStdMolecule);
          }

        }
        else {

        }
      }
      else {
      }
    }
}

// ============================================================
// Function : atomFiller
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::atomFiller(DOMNode *rootnode)
{
    if (rootnode->hasAttributes()) {
      // get all the attributes of the node
      DOMNamedNodeMap *pAttributes = rootnode->getAttributes();
      int nSize = pAttributes->getLength();
      int bType = 1;
      int bKind = 0;
      int bTop = 0;

      for (int i = 0; i < nSize; ++i) {
        DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
        char *name = XMLString::transcode(pAttributeNode->getName());
        if (std::string(name) == std::string("identity")) {
          char *cAtomName = XMLString::transcode(pAttributeNode->getNodeValue());
          if (std::string(cAtomName).size() != 4) {
            std::cout << "|" << std::string(cAtomName) << "|" << std::endl;
            std::cout << " Error in stdLib xml file ... exiting " << std::endl;
            //exit(0);
            throw MTKException(" Error in stdLib xml file ... exiting ");
          }
          pStdAtom->identity = std::string(cAtomName);
          delete cAtomName;
        }
        if (std::string(name) == std::string("index")) {
          char *cIndex = XMLString::transcode(pAttributeNode->getNodeValue());
          pStdAtom->index = atoi(cIndex);
          delete cIndex;
        }
        if (std::string(name) == std::string("type")) {
          char *cType = XMLString::transcode(pAttributeNode->getNodeValue());
          pStdAtom->type = std::string(cType);
          delete cType;
        }
        if (std::string(name) == std::string("chain")) {
          char *cChain = XMLString::transcode(pAttributeNode->getNodeValue());
          pStdAtom->chain = std::string(cChain);
          delete cChain;
        }
        if (std::string(name) == std::string("atmCharge")) {
          double atmCharge = strtod(XC(pAttributeNode->getNodeValue()), 0);
          pStdAtom->atmCharge = atmCharge; // * 18.2223 ; // electron to kcal
        }
        if (std::string(name) == std::string("bond12")) {
          char *cBond12 = XMLString::transcode(pAttributeNode->getNodeValue());
          pStdAtom->bond12 = atoi(cBond12);
          delete cBond12;
        }
        if (std::string(name) == std::string("bType")) {
          char *cBType = XMLString::transcode(pAttributeNode->getNodeValue());
          bType = atoi(cBType);
          delete cBType;
        }
        if (std::string(name) == std::string("bKind")) {
          char *cBKind = XMLString::transcode(pAttributeNode->getNodeValue());
          bKind = atoi(cBKind);
          delete cBKind;
        }
        if (std::string(name) == std::string("bTop")) {
          char *cBTop = XMLString::transcode(pAttributeNode->getNodeValue());
          bTop = atoi(cBTop);
          delete cBTop;
        }
        if (std::string(name) == std::string("bond13")) {
          char *cBond13 = XMLString::transcode(pAttributeNode->getNodeValue());
          pStdAtom->bond13 = atoi(cBond13);
          delete cBond13;
        }
        if (std::string(name) == std::string("bond14")) {
          char *cBond14 = XMLString::transcode(pAttributeNode->getNodeValue());
          pStdAtom->bond14 = atoi(cBond14);
          delete cBond14;
        }
        if (std::string(name) == std::string("bondLength")) {
          double bondLength = strtod(XC(pAttributeNode->getNodeValue()), 0);
          pStdAtom->bondLength = bondLength;
        }
        if (std::string(name) == std::string("bondAngle")) {
          double bondAngle = strtod(XC(pAttributeNode->getNodeValue()), 0);
          pStdAtom->bondAngle = bondAngle;
        }
        if (std::string(name) == std::string("bondTorsion")) {
          double bondTorsion = strtod(XC(pAttributeNode->getNodeValue()), 0);
          pStdAtom->bondTorsion = bondTorsion;
        }
        //if (std::string(name) == std::string("atNum")) {
        //  char *cAtNum = XMLString::transcode(pAttributeNode->getNodeValue());
        //  pStdAtom->atNum = atoi(cAtNum);
        //}
        //if (std::string(name) == std::string("symbol")) {
        //  char *cSymbol = XMLString::transcode(pAttributeNode->getNodeValue());
        //  pStdAtom->atSymbol = std::string(cSymbol);
        //}
        //if (std::string(name) == std::string("hybridization")) {
        //  char *cHybrid = XMLString::transcode(pAttributeNode->getNodeValue());
        //  pStdAtom->hybridization = std::string(cHybrid);
        //}
        if (std::string(name) == std::string("kind")) {
          char *cKind = XMLString::transcode(pAttributeNode->getNodeValue());
          pStdAtom->kind = atoi(cKind);
          delete cKind;
        }
        delete name;
      }
      if (pStdAtom->type != "") {
        std::string atSymbol = pParameters->getAtomTypeSymbol(pStdAtom->type);
        if (atSymbol == "") {
          std::cout << " Cannot find element symbol for " << pStdAtom->type << std::endl;
          std::cout << " Error in Library file ... exiting " << std::endl;
          //exit(0);
          throw MTKException(" Error in stdLib xml file ... exiting ");
        }
        pStdAtom->atSymbol = atSymbol;
        pStdAtom->atNum = pParameters->getAtomicNum(pStdAtom->type);
      }
      if (pStdAtom->index && pStdAtom->bondLength) {
        pStdBond = pStdFrag->addStdBond();
        pStdBond->atom1 = pStdAtom->index;
        pStdBond->atom2 = pStdAtom->bond12;
        pStdBond->type = bType;
        pStdBond->kind = bKind;
        pStdBond->topology = bTop;
        pStdBond->stereo = 0;
        bType = 1;
        bKind = 0;
        bTop = 0;
        pStdBond->length = pStdAtom->bondLength;
      }
    }
}

// ============================================================
// Function : aliasFiller
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::aliasFiller(DOMNode *rootnode)
{
    if (rootnode->hasAttributes()) {
      char *cAlias = 0;
      char *cOriginal = 0;
      // get all the attributes of the node
      DOMNamedNodeMap *pAttributes = rootnode->getAttributes();
      int nSize = pAttributes->getLength();

      for (int i = 0; i < nSize; ++i) {
        DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
        char *name = XMLString::transcode(pAttributeNode->getName());
        if (std::string(name) == std::string("alias")) {
          cAlias = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("original")) {
          cOriginal = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        delete name;
      }
      if (cAlias && cOriginal) {
        pStdAlias = pStdFrag->addStdAlias();
        pStdAlias->atom1 = std::string(cOriginal);
        pStdAlias->atom2 = std::string(cAlias);
      }
      else {
      }
      delete cAlias;
      delete cOriginal;
    }
}

// ============================================================
// Function : improperFiller
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::improperFiller(DOMNode *rootnode)
{
    if (rootnode->hasAttributes()) {
      char *cAtom1 = 0;
      char *cAtom2 = 0;
      char *cAtom3 = 0;
      char *cAtom4 = 0;
      // get all the attributes of the node
      DOMNamedNodeMap *pAttributes = rootnode->getAttributes();
      int nSize = pAttributes->getLength();

      for (int i = 0; i < nSize; ++i) {
        DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
        char *name = XMLString::transcode(pAttributeNode->getName());
        if (std::string(name) == std::string("Atom1")) {
          cAtom1 = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("Atom2")) {
          cAtom2 = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("Atom3")) {
          cAtom3 = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("Atom4")) {
          cAtom4 = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        delete name;
      }

      if (cAtom1 && cAtom2 && cAtom3 && cAtom4) {
        pStdImproper = pStdFrag->addStdImproper();

/*
        if (std::string(cAtom1) == "p1") {
          pStdImproper->atom1 = -1;
        }
        else if (std::string(cAtom1) == "n1") {
          pStdImproper->atom1 = -4;
        }
        else {
          pStdImproper->atom1 = atoi(cAtom1);
        }
*/
        pStdImproper->atom1 = atoi(cAtom1);
        pStdImproper->atom2 = atoi(cAtom2);
        pStdImproper->atom3 = atoi(cAtom3);
        pStdImproper->atom4 = atoi(cAtom4);
      }
      else {
      }

      delete cAtom1;
      delete cAtom2;
      delete cAtom3;
      delete cAtom4;
    }
}

// ============================================================
// Function : loopFiller
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::loopFiller(DOMNode *rootnode)
{
    if (rootnode->hasAttributes()) {
      char *cAtom1 = 0;
      char *cAtom2 = 0;
      int bType = 1;
      // get all the attributes of the node
      DOMNamedNodeMap *pAttributes = rootnode->getAttributes();
      int nSize = pAttributes->getLength();

      for (int i = 0; i < nSize; ++i) {
        DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
        char *name = XMLString::transcode(pAttributeNode->getName());
        if (std::string(name) == std::string("Atom1")) {
          cAtom1 = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("Atom2")) {
          cAtom2 = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("bType")) {
          char *cBType = XMLString::transcode(pAttributeNode->getNodeValue());
          bType = atoi(cBType);
          delete cBType;
        }
        delete name;
      }
      if (cAtom1 && cAtom2) {
        pStdLoop = pStdFrag->addStdLoop();
        pStdLoop->atom1 = atoi(cAtom1);
        pStdLoop->atom2 = atoi(cAtom2);
        pStdLoop->type = bType;
        bType = 1;
      }
      else {
      }
      delete cAtom1;
      delete cAtom2;
    }
}

// ============================================================
// Function : ringFiller
// ------------------------------------------------------------
// <ring size="6" atoms="1 2 3 4 5 14" planar="0" aromatic="0"
//  hetero="0" nHetero="0" nNitrogen="0" nOxygen="0" nSulfur="0"/>
// ------------------------------------------------------------
void stdLibParser::ringFiller(DOMNode *rootnode)
{
    if (rootnode->hasAttributes()) {
      char *cSize = 0;
      char *cAtoms = 0;
      char *cPlanar = 0;
      char *cAromatic = 0;
      char *cHetero = 0;
      char *cNHetero = 0;
      char *cNNitrogen = 0;
      char *cNOxygen = 0;
      char *cNSulfur = 0;

      // get all the attributes of the node
      DOMNamedNodeMap *pAttributes = rootnode->getAttributes();
      int nSize = pAttributes->getLength();
      for (int i = 0; i < nSize; ++i) {
        DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
        char *name = XMLString::transcode(pAttributeNode->getName());
        if (std::string(name) == std::string("size")) {
          cSize = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("atoms")) {
          cAtoms = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("planar")) {
          cPlanar = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("aromatic")) {
          cAromatic = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("hetero")) {
          cHetero = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("nHetero")) {
          cNHetero = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("nNitrogen")) {
          cNNitrogen = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("nOxygen")) {
          cNOxygen = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("nSulfur")) {
          cNSulfur = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        delete name;
      }
      if (cAtoms) {
        pStdRing = pStdFrag->addStdRing();
        std::string sAtoms = std::string(cAtoms);
        std::vector<std::string> vsAtoms;
        splitString(sAtoms, " ", vsAtoms, 0);
        for (unsigned int x = 0; x < vsAtoms.size(); x++) {
          pStdRing->atoms.push_back(atoi(vsAtoms[x].c_str()));
        }
        if (cSize) {
          if (atoi(cSize) == static_cast<int>(vsAtoms.size())) {
            pStdRing->size = atoi(cSize);
          }
        }
        if (cPlanar) {
          pStdRing->planar = atoi(cPlanar);
        }
        if (cAromatic) {
          pStdRing->aromatic = atoi(cAromatic);
        }
        if (cHetero) {
          pStdRing->hetero = atoi(cHetero);
        }
        if (cNHetero) {
          pStdRing->nHetero = atoi(cNHetero);
        }
        if (cNNitrogen) {
          pStdRing->nNitrogen = atoi(cNNitrogen);
        }
        if (cNOxygen) {
          pStdRing->nOxygen = atoi(cNOxygen);
        }
        if (cNSulfur) {
          pStdRing->nSulfur = atoi(cNSulfur);
        }
      }
      delete cSize;
      delete cAtoms;
      delete cPlanar;
      delete cAromatic;
      delete cHetero;
      delete cNHetero;
      delete cNNitrogen;
      delete cNOxygen;
      delete cNSulfur;
    }
}

// ============================================================
// Function : featureFiller
// ------------------------------------------------------------
// <feature identity="PIC" atoms="1 2 3 4 5 14"/>
// ------------------------------------------------------------
void stdLibParser::featureFiller(DOMNode *rootnode)
{
    if (rootnode->hasAttributes()) {
      char *cIdentity = 0;
      char *cAtoms = 0;

      // get all the attributes of the node
      DOMNamedNodeMap *pAttributes = rootnode->getAttributes();
      int nSize = pAttributes->getLength();
      for (int i = 0; i < nSize; ++i) {
        DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
        char *name = XMLString::transcode(pAttributeNode->getName());
        if (std::string(name) == std::string("identity")) {
          cIdentity = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("atoms")) {
          cAtoms = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        delete name;
      }
      if (cAtoms and cIdentity) {
        pStdFeature = pStdFrag->addStdFeature();
        pStdFeature->name = std::string(cIdentity);
        std::string sAtoms = std::string(cAtoms);
        std::vector<std::string> vsAtoms;
        splitString(sAtoms, " ", vsAtoms, 0);
        for (unsigned int x = 0; x < vsAtoms.size(); x++) {
          pStdFeature->atoms.push_back(atoi(vsAtoms[x].c_str()));
        }
      }
      delete cIdentity;
      delete cAtoms;
    }
}

// ============================================================
// Function : funcGroupFiller
// ------------------------------------------------------------
// <funcGroup group="terminal" fragment="CH3" atoms="1 2 3 4"/>
// ------------------------------------------------------------
void stdLibParser::funcGroupFiller(DOMNode *rootnode)
{
    if (rootnode->hasAttributes()) {
      char *cGroup = 0;
      char *cFrag = 0;
      char *cAtoms = 0;

      // get all the attributes of the node
      DOMNamedNodeMap *pAttributes = rootnode->getAttributes();
      int nSize = pAttributes->getLength();
      for (int i = 0; i < nSize; ++i) {
        DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
        char *name = XMLString::transcode(pAttributeNode->getName());
        if (std::string(name) == std::string("group")) {
          cGroup = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("fragment")) {
          cFrag = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("atoms")) {
          cAtoms = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        delete name;
      }
      if (cAtoms and cGroup and cFrag) {
        pStdFuncGroup = pStdFrag->addStdFuncGroup();
        pStdFuncGroup->groupName = std::string(cGroup);
        pStdFuncGroup->fragName = std::string(cFrag);
        std::string sAtoms = std::string(cAtoms);
        std::vector<std::string> vsAtoms;
        splitString(sAtoms, " ", vsAtoms, 0);
        for (unsigned int x = 0; x < vsAtoms.size(); x++) {
          pStdFuncGroup->atoms.push_back(atoi(vsAtoms[x].c_str()));
        }
      }
      delete cGroup;
      delete cFrag;
      delete cAtoms;
    }
}

// ============================================================
// Function : connPtsFiller
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::connPtsFiller(DOMNode *rootnode)
{
    if (rootnode->hasAttributes()) {
      char *cAtoms = 0;
      // get all the attributes of the node
      DOMNamedNodeMap *pAttributes = rootnode->getAttributes();
      int nSize = pAttributes->getLength();
      for (int i = 0; i < nSize; ++i) {
        DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
        char *name = XMLString::transcode(pAttributeNode->getName());
        if (std::string(name) == std::string("atoms")) {
          cAtoms = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        delete name;
      }
      if (cAtoms) {
        std::string sAtoms = std::string(cAtoms);
        std::vector<std::string> vsAtoms;
        std::vector<int> ivsAtoms;
        splitString(sAtoms, " ", vsAtoms, 0);
        for (unsigned int x = 0; x < vsAtoms.size(); x++) {
          ivsAtoms.push_back(atoi(vsAtoms[x].c_str()));
        }
        pStdFrag->addStdConnPts(ivsAtoms);
      }
      else {
      }
      delete cAtoms;
    }
}

// ============================================================
// Function : connTorFiller
// ------------------------------------------------------------
// <connectionTorsion bd="1" ag="4" tr="5" value="94.54"/>
// ------------------------------------------------------------
void stdLibParser::connTorFiller(DOMNode *rootnode)
{
    if (rootnode->hasAttributes()) {
      char *cBd = 0;
      char *cAg = 0;
      char *cTr = 0;
      char *cValue = 0;
      // get all the attributes of the node
      DOMNamedNodeMap *pAttributes = rootnode->getAttributes();
      int nSize = pAttributes->getLength();
      for (int i = 0; i < nSize; ++i) {
        DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
        char *name = XMLString::transcode(pAttributeNode->getName());
        if (std::string(name) == std::string("bd")) {
          cBd = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("ag")) {
          cAg = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("tr")) {
          cTr = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("value")) {
          cValue = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        delete name;
      }
      // Need to do some error checking
      if (cBd and cAg and cTr and cValue) {
        stdConnTorsion* pStdConnTorsion = pStdFrag->addStdConnTorsion();
        pStdConnTorsion->bondAtom = atoi(cBd);
        pStdConnTorsion->angleAtom = atoi(cAg);
        pStdConnTorsion->torsionAtom = atoi(cTr);
        pStdConnTorsion->torsion = strtod(cValue,0);
      }
      else {
        std::cout << " Error in stdLibParser::connTorFiller ... exiting " << std::endl;
        exit(0);
      }
      delete cBd;
      delete cAg;
      delete cTr;
      delete cValue;
    }
}

// ============================================================
// Function : rotBondFiller
// ------------------------------------------------------------
// <rotatableBond atom1="14" atom2="13" atom3="8" atom4="9" values="15.00"/>
// ------------------------------------------------------------
void stdLibParser::rotBondFiller(DOMNode *rootnode)
{
    if (rootnode->hasAttributes()) {
      char *cAt1 = 0;
      char *cAt2 = 0;
      char *cAt3 = 0;
      char *cAt4 = 0;
      char *cValues = 0;
      // get all the attributes of the node
      DOMNamedNodeMap *pAttributes = rootnode->getAttributes();
      int nSize = pAttributes->getLength();
      for (int i = 0; i < nSize; ++i) {
        DOMAttr *pAttributeNode = (DOMAttr*) pAttributes->item(i);
        char *name = XMLString::transcode(pAttributeNode->getName());
        if (std::string(name) == std::string("atom1")) {
          cAt1 = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("atom2")) {
          cAt2 = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("atom3")) {
          cAt3 = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("atom4")) {
          cAt4 = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        if (std::string(name) == std::string("values")) {
          cValues = XMLString::transcode(pAttributeNode->getNodeValue());
        }
        delete name;
      }
      // Need to do some error checking
      if (cAt1 and cAt2 and cAt3 and cAt4) {
        stdRotBond* pStdRotBond = pStdFrag->addStdRotBond();
        pStdRotBond->atom1 = atoi(cAt1);
        pStdRotBond->atom2 = atoi(cAt2);
        pStdRotBond->atom3 = atoi(cAt3);
        pStdRotBond->atom4 = atoi(cAt4);
        if (cValues) {
          std::string sValues = std::string(cValues);
          std::vector<std::string> vsValues;
          splitString(sValues, " ", vsValues, 0);
          for (unsigned int x = 0; x < vsValues.size(); x++) {
            pStdRotBond->values.push_back(strtod(vsValues[x].c_str(),0));
          }
        }
      }
      else {
        std::cout << " Error in stdLibParser::connTorFiller ... exiting " << std::endl;
        //exit(0);
        throw MTKException(" Error in stdLibParser::connTorFiller ... exiting ");
      }
      delete cAt1;
      delete cAt2;
      delete cAt3;
      delete cAt4;
      delete cValues;
    }
}
#endif // USE_XERCES

    // ----------------------------- //
    // -    Qt WRITE FUNCTIONS     - //
    // ----------------------------- //
#ifdef USE_QT
// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write library xml files using Qt
// ------------------------------------------------------------
void stdLibParser::Write(std::string fileName)
{
    QString qFileName = QString::fromStdString(fileName);
    std::string errMessage = " Writing " + fileName;
    errorLogger.throwError("stdLibParser", errMessage, INFO);

    // Write library file using DOM into memory
    QString mtkppInfo = PACKAGE_TARNAME;
    QDomDocument doc(mtkppInfo);

    QDomElement rootElem = doc.createElement("stdLib");
    doc.appendChild(rootElem);

    // stdGroup's
    std::vector<stdGroup*> groupList = pStdLibrary->getStdGroupList();
    for (unsigned int i = 0; i < groupList.size(); i++) {
      QDomElement  groupElem = doc.createElement("group");
      rootElem.appendChild(groupElem);

      groupElem.setAttribute("identity", string2QString(groupList[i]->getName()));

      // stdFrag's
      std::vector<stdFrag*> fragList = groupList[i]->getStdFragList();
      for (unsigned int j = 0; j < fragList.size(); j++) {
        QDomElement  fragElem = doc.createElement("fragment");
        groupElem.appendChild(fragElem);

        fragElem.setAttribute("identity", string2QString(fragList[j]->getName()));
        fragElem.setAttribute("symbol" , string2QString(fragList[j]->getSymbol()));
        fragElem.setAttribute("code"   , string2QString(fragList[j]->getCode()));
        fragElem.setAttribute("type"   , string2QString(fragList[j]->getType()));

        std::string sGstr = fragList[j]->getSubGraphStr();
        if (sGstr != "") fragElem.setAttribute("subGraphs", string2QString(sGstr));

        // stdAtom's
        std::vector<stdAtom*> atomList = fragList[j]->getStdAtomList();
        for (unsigned int k = 0; k < atomList.size(); k++) {
          QDomElement  atomElem = doc.createElement("atom");
          fragElem.appendChild(atomElem);

          atomElem.setAttribute("identity"   , string2QString(atomList[k]->identity));
          atomElem.setAttribute("index"      , int2QString(atomList[k]->index));
          atomElem.setAttribute("type"       , string2QString(atomList[k]->type));
          atomElem.setAttribute("chain"      , string2QString(atomList[k]->chain));
          atomElem.setAttribute("atmCharge"  , double2QString(atomList[k]->atmCharge));
          atomElem.setAttribute("bond12"     , int2QString(atomList[k]->bond12));
          atomElem.setAttribute("bondLength" , double2QString(atomList[k]->bondLength));
          atomElem.setAttribute("bond13"     , int2QString(atomList[k]->bond13));
          atomElem.setAttribute("bondAngle"  , double2QString(atomList[k]->bondAngle));
          atomElem.setAttribute("bond14"     , int2QString(atomList[k]->bond14));
          atomElem.setAttribute("bondTorsion", double2QString(atomList[k]->bondTorsion));

          atomElem.setAttribute("kind", int2QString(atomList[k]->kind));
          pStdBond = fragList[j]->getStdBond(atomList[k]->index,atomList[k]->bond12);
          if (pStdBond) {
            atomElem.setAttribute("bType", int2QString(pStdBond->type));
            atomElem.setAttribute("bTop", int2QString(pStdBond->topology));
            atomElem.setAttribute("bKind", int2QString(pStdBond->kind));
            atomElem.setAttribute("bondLength", double2QString(pStdBond->length));
          }
          else {
            atomElem.setAttribute("bType", int2QString(1));
            atomElem.setAttribute("bTop", int2QString(0));
            atomElem.setAttribute("bKind", int2QString(0));
          }
        }

        // stdLoop's
        std::vector<stdLoop*> loopList = fragList[j]->getStdLoopList();
        for (unsigned int k = 0; k < loopList.size(); k++) {
          QDomElement loopElem = doc.createElement("loop");
          fragElem.appendChild(loopElem);

          loopElem.setAttribute("Atom1"  , int2QString(loopList[k]->atom1));
          loopElem.setAttribute("Atom2"  , int2QString(loopList[k]->atom2));
          loopElem.setAttribute("bType"  , int2QString(loopList[k]->type));
          loopElem.setAttribute("bStereo", int2QString(loopList[k]->stereo));
        }

        // stdAlias's
        std::vector<stdAlias*> aliasList = fragList[j]->getStdAliasList();
        for (unsigned int k = 0; k < aliasList.size(); k++) {
          QDomElement aliasElem = doc.createElement("alias");
          fragElem.appendChild(aliasElem);

          aliasElem.setAttribute("original", string2QString(aliasList[k]->atom1));
          aliasElem.setAttribute("alias"   , string2QString(aliasList[k]->atom2));
        }

        // stdImproper's
        std::vector<stdImproper*> improperList = fragList[j]->getStdImproperList();
        for (unsigned int k = 0; k < improperList.size(); k++) {
          QDomElement impElem = doc.createElement("improper");
          fragElem.appendChild(impElem);

          if (improperList[k]->atom1 == -1) {
            impElem.setAttribute("Atom1", "p1");
          }
          else if (improperList[k]->atom1 == -4) {
            impElem.setAttribute("Atom1", "n1");
          }
          else {
            impElem.setAttribute("Atom1", int2QString(improperList[k]->atom1));
          }
          impElem.setAttribute("Atom2", int2QString(improperList[k]->atom2));
          impElem.setAttribute("Atom3", int2QString(improperList[k]->atom3));
          impElem.setAttribute("Atom4", int2QString(improperList[k]->atom4));
        }

        // stdRing's
        std::vector<stdRing*> ringList = fragList[j]->getStdRingList();
        for (unsigned int k = 0; k < ringList.size(); k++) {
          QDomElement ringElem = doc.createElement("ring");
          fragElem.appendChild(ringElem);

          ringElem.setAttribute("size"      , int2QString(ringList[k]->size));
          ringElem.setAttribute("planar"    , int2QString(ringList[k]->planar));
          ringElem.setAttribute("aromatic"  , int2QString(ringList[k]->aromatic));
          ringElem.setAttribute("hetero"    , int2QString(ringList[k]->hetero));
          ringElem.setAttribute("nHetero"   , int2QString(ringList[k]->nHetero));
          ringElem.setAttribute("nNitrogen" , int2QString(ringList[k]->nNitrogen));
          ringElem.setAttribute("nOxygen"   , int2QString(ringList[k]->nOxygen));
          ringElem.setAttribute("nSulfur"   , int2QString(ringList[k]->nSulfur));

          std::string ringAtoms = "";
          for (unsigned int x = 0; x < ringList[k]->atoms.size(); x++) {
            ringAtoms+= (int2String(ringList[k]->atoms[x]) + " ");
          }
          ringElem.setAttribute("atoms", string2QString(ringAtoms));
        }

        // stdFeature's
        std::vector<stdFeature*> featureList = fragList[j]->getStdFeatureList();
        for (unsigned int k = 0; k < featureList.size(); k++) {
          QDomElement featureElem = doc.createElement("feature");
          fragElem.appendChild(featureElem);

          featureElem.setAttribute("identity", string2QString(featureList[k]->name));

          std::string featureAtoms = "";
          for (unsigned int x = 0; x < featureList[k]->atoms.size(); x++) {
            featureAtoms+= (int2String(featureList[k]->atoms[x]) + " ");
          }
          featureElem.setAttribute("atoms", string2QString(featureAtoms));
        }

        // stdFuncGroup's
        std::vector<stdFuncGroup*> funcGroupList = fragList[j]->getStdFuncGroupList();
        for (unsigned int k = 0; k < funcGroupList.size(); k++) {
          QDomElement funcGroupElem = doc.createElement("funcGroup");
          fragElem.appendChild(funcGroupElem);

          funcGroupElem.setAttribute("group"    , string2QString(funcGroupList[k]->groupName));
          funcGroupElem.setAttribute("fragment" , string2QString(funcGroupList[k]->fragName));

          std::string fgroupAtoms = "";
          for (unsigned int x = 0; x < funcGroupList[k]->atoms.size(); x++) {
            fgroupAtoms += (int2String(funcGroupList[k]->atoms[x]) + " ");
          }
          funcGroupElem.setAttribute("atoms", string2QString(fgroupAtoms));
        }

        // stdConnPts
        std::vector<int> connPtsList = fragList[j]->getStdConnPtsList();
        if (connPtsList.size() > 0) {
          QDomElement connPtsElem = doc.createElement("connectionPoints");
          fragElem.appendChild(connPtsElem);

          std::string connPtsAtoms = "";
          for (unsigned int x = 0; x < connPtsList.size(); x++) {
            connPtsAtoms+= (int2String(connPtsList[x]) + " ");
          }
          connPtsElem.setAttribute("atoms", string2QString(connPtsAtoms));
        }
      } // fragment
    } // group

    // Write dom xml to file
    QString domXml = doc.toString();

    QFile file(qFileName);
    if (!file.open(QIODevice::WriteOnly)) {
      errMessage = " Writing " + fileName;
      errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
      return;
    }

    QTextStream out(&file);
    out << domXml;
    file.close();
}

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write library xml files using Qt
// ------------------------------------------------------------
void stdLibParser::Write(std::string fileName, std::string groupName)
{
    QString qFileName = QString::fromStdString(fileName);
    std::string errMessage = " Writing " + fileName;
    errorLogger.throwError("stdLibParser", errMessage, INFO);

    QString mtkppInfo = PACKAGE_TARNAME;
    QDomDocument doc(mtkppInfo);

    QDomElement rootElem = doc.createElement("stdLib");
    doc.appendChild(rootElem);

    std::vector<stdGroup*> groupList = pStdLibrary->getStdGroupList();
    for (unsigned int i = 0; i < groupList.size(); i++) {
      if (groupName == groupList[i]->getName()) {
        this->writeGroup(doc, groupList[i]);
      }
    }

    // Write dom xml to file
    QString domXml = doc.toString();

    QFile file(qFileName);
    if (!file.open(QIODevice::WriteOnly)) {
      errMessage = " Writing " + fileName;
      errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
      return;
    }

    QTextStream out(&file);
    out << domXml;
    file.close();
}

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write library xml files using Qt
// ------------------------------------------------------------
void stdLibParser::Write(std::string fileName, std::string groupName, std::string fragName)
{
    QString qFileName = QString::fromStdString(fileName);
    std::string errMessage = " Writing " + fileName;
    errorLogger.throwError("stdLibParser", errMessage, INFO);

    // Read elements file using DOM into memory
    QString mtkppInfo = PACKAGE_TARNAME;
    QDomDocument doc(mtkppInfo);

    QDomElement rootElem = doc.createElement("stdLib");
    doc.appendChild(rootElem);

    std::vector<stdGroup*> groupList = pStdLibrary->getStdGroupList();
    for (unsigned int i = 0; i < groupList.size(); i++) {
      if (groupName == groupList[i]->getName()) {
        this->writeGroup(doc, groupList[i], fragName);
      }
    }

    // Write dom xml to file
    QString domXml = doc.toString();

    QFile file(qFileName);
    if (!file.open(QIODevice::WriteOnly)) {
      errMessage = " Writing " + fileName;
      errorLogger.throwError("stdLibParser", errMessage, MTK_ERROR);
      return;
    }

    QTextStream out(&file);
    out << domXml;
    file.close();
}

// ============================================================
// Function : writeGroup
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeGroup(QDomDocument doc, stdGroup* stdGp)
{
    QDomElement rootElem = doc.documentElement();
    QDomElement groupElem = doc.createElement("group");
    rootElem.appendChild(groupElem);
    groupElem.setAttribute("identity", string2QString(stdGp->getName()));

    // stdFrag's
    std::vector<stdFrag*> fragList = stdGp->getStdFragList();
    for (unsigned int j = 0; j < fragList.size(); j++) {
      this->writeFrag(doc, groupElem, fragList[j]);
    }

    // standard structure
    if (stdGp->hasStdMolecule()) {
      mtkppParser* pMTKppParser = new mtkppParser;
      pMTKppParser->writeMolecule(doc, groupElem, stdGp->getStdMolecule());
    }
}

// ============================================================
// Function : writeGroup
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeGroup(QDomDocument doc, stdGroup* stdGp, std::string frag3L)
{
    QDomElement rootElem = doc.documentElement();
    QDomElement groupElem = doc.createElement("group");
    rootElem.appendChild(groupElem);
    groupElem.setAttribute("identity", string2QString(stdGp->getName()));

    // stdFrag's
    std::vector<stdFrag*> fragList = stdGp->getStdFragList();
    for (unsigned int j = 0; j < fragList.size(); j++) {
      if (fragList[j]->getSymbol() == frag3L) {
        this->writeFrag(doc, groupElem, fragList[j]);
      }
    }
}

// ============================================================
// Function : writeFrag
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeFrag(QDomDocument doc, QDomElement groupElem, stdFrag* stdFg)
{
    QDomElement fragElem = doc.createElement("fragment");
    groupElem.appendChild(fragElem);

    fragElem.setAttribute("identity", string2QString(stdFg->getName()));
    fragElem.setAttribute("symbol", string2QString(stdFg->getSymbol()));
    fragElem.setAttribute("character", string2QString(stdFg->getCharacter()));
    fragElem.setAttribute("code", string2QString(stdFg->getCode()));
    fragElem.setAttribute("type", string2QString(stdFg->getType()));

    std::string sGstr = stdFg->getSubGraphStr();
    if (sGstr != "") fragElem.setAttribute("subGraphs", string2QString(sGstr));

    // stdAtom's
    std::vector<stdAtom*> atomList = stdFg->getStdAtomList();
    for (unsigned int k = 0; k < atomList.size(); k++) {
      this->writeAtom(doc, fragElem, stdFg, atomList[k]);
    }

    // stdLoop's
    std::vector<stdLoop*> loopList = stdFg->getStdLoopList();
    for (unsigned int k = 0; k < loopList.size(); k++) {
      this->writeLoop(doc, fragElem, loopList[k]);
    }

    // stdAlias's
    std::vector<stdAlias*> aliasList = stdFg->getStdAliasList();
    for (unsigned int k = 0; k < aliasList.size(); k++) {
      this->writeAlias(doc, fragElem, aliasList[k]);
    }

    // stdImproper's
    std::vector<stdImproper*> improperList = stdFg->getStdImproperList();
    for (unsigned int k = 0; k < improperList.size(); k++) {
      this->writeImproper(doc, fragElem, improperList[k]);
    }

    // stdRing's
    std::vector<stdRing*> ringList = stdFg->getStdRingList();
    for (unsigned int k = 0; k < ringList.size(); k++) {
      this->writeRing(doc, fragElem, ringList[k]);
    }

    // stdFeature's
    std::vector<stdFeature*> featureList = stdFg->getStdFeatureList();
    for (unsigned int k = 0; k < featureList.size(); k++) {
      this->writeFeature(doc, fragElem, featureList[k]);
    }

    // stdFuncGroup's
    std::vector<stdFuncGroup*> funcGroupList = stdFg->getStdFuncGroupList();
    for (unsigned int k = 0; k < funcGroupList.size(); k++) {
      this->writeFuncGroup(doc, fragElem, funcGroupList[k]);
    }

    // stdConnPts
    this->writeConnPts(doc, fragElem, stdFg);
}

// ============================================================
// Function : writeAtom
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeAtom(QDomDocument doc, QDomElement  fragElem, stdFrag* stdFg, stdAtom* stdAt)
{
    // stdAtom's
    QDomElement atomElem = doc.createElement("atom");
    fragElem.appendChild(atomElem);

    atomElem.setAttribute("identity"   , string2QString(stdAt->identity));
    atomElem.setAttribute("index"      , int2QString(stdAt->index));
    atomElem.setAttribute("type"       , string2QString(stdAt->type));
    atomElem.setAttribute("chain"      , string2QString(stdAt->chain));
    atomElem.setAttribute("atmCharge"  , double2QString(stdAt->atmCharge));
    atomElem.setAttribute("bond12"     , int2QString(stdAt->bond12));
    atomElem.setAttribute("bondLength" , double2QString(stdAt->bondLength));
    atomElem.setAttribute("bond13"     , int2QString(stdAt->bond13));
    atomElem.setAttribute("bondAngle"  , double2QString(stdAt->bondAngle));
    atomElem.setAttribute("bond14"     , int2QString(stdAt->bond14));
    atomElem.setAttribute("bondTorsion", double2QString(stdAt->bondTorsion));

    atomElem.setAttribute("kind", int2QString(stdAt->kind));
    pStdBond = stdFg->getStdBond(stdAt->index, stdAt->bond12);
    if (pStdBond) {
      atomElem.setAttribute("bType", int2QString(pStdBond->type));
      atomElem.setAttribute("bTop", int2QString(pStdBond->topology));
      atomElem.setAttribute("bKind", int2QString(pStdBond->kind));
      atomElem.setAttribute("bondLength", double2QString(pStdBond->length));
    }
    else {
      atomElem.setAttribute("bType", int2QString(1));
      atomElem.setAttribute("bTop", int2QString(0));
      atomElem.setAttribute("bKind", int2QString(0));
    }
}

// ============================================================
// Function : writeLoop
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeLoop(QDomDocument doc, QDomElement fragElem, stdLoop* stdLp)
{
    QDomElement loopElem = doc.createElement("loop");
    fragElem.appendChild(loopElem);

    loopElem.setAttribute("Atom1"  , int2QString(stdLp->atom1));
    loopElem.setAttribute("Atom2"  , int2QString(stdLp->atom2));
    loopElem.setAttribute("bType"  , int2QString(stdLp->type));
    loopElem.setAttribute("bStereo", int2QString(stdLp->stereo));
}

// ============================================================
// Function : writeAlias
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeAlias(QDomDocument doc, QDomElement fragElem, stdAlias* stdAl)
{
    QDomElement aliasElem = doc.createElement("alias");
    fragElem.appendChild(aliasElem);

    aliasElem.setAttribute("original", string2QString(stdAl->atom1));
    aliasElem.setAttribute("alias"   , string2QString(stdAl->atom2));
}

// ============================================================
// Function : writeImproper
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeImproper(QDomDocument doc, QDomElement fragElem, stdImproper* stdIm)
{
    QDomElement impElem = doc.createElement("improper");
    fragElem.appendChild(impElem);
/*
    if (stdIm->atom1 == -1) {
      impElem.setAttribute("Atom1", "p1");
    }
    else if (stdIm->atom1 == -4) {
      impElem.setAttribute("Atom1", "n1");
    }
    else {
      impElem.setAttribute("Atom1", int2QString(stdIm->atom1));
    }
*/

    impElem.setAttribute("Atom1", int2QString(stdIm->atom1));
    impElem.setAttribute("Atom2", int2QString(stdIm->atom2));
    impElem.setAttribute("Atom3", int2QString(stdIm->atom3));
    impElem.setAttribute("Atom4", int2QString(stdIm->atom4));
}

// ============================================================
// Function : writeRing
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeRing(QDomDocument doc, QDomElement fragElem, stdRing* stdRg)
{
    QDomElement ringElem = doc.createElement("ring");
    fragElem.appendChild(ringElem);

    ringElem.setAttribute("size"      , int2QString(stdRg->size));
    ringElem.setAttribute("planar"    , int2QString(stdRg->planar));
    ringElem.setAttribute("aromatic"  , int2QString(stdRg->aromatic));
    ringElem.setAttribute("hetero"    , int2QString(stdRg->hetero));
    ringElem.setAttribute("nHetero"   , int2QString(stdRg->nHetero));
    ringElem.setAttribute("nNitrogen" , int2QString(stdRg->nNitrogen));
    ringElem.setAttribute("nOxygen"   , int2QString(stdRg->nOxygen));
    ringElem.setAttribute("nSulfur"   , int2QString(stdRg->nSulfur));

    std::string ringAtoms = "";
    for (unsigned int x = 0; x < stdRg->atoms.size(); x++) {
      ringAtoms+= (int2String(stdRg->atoms[x]) + " ");
    }
    ringElem.setAttribute("atoms", string2QString(ringAtoms));
}

// ============================================================
// Function : writeFeature
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeFeature(QDomDocument doc, QDomElement fragElem, stdFeature* stdFt)
{
    QDomElement featureElem = doc.createElement("feature");
    fragElem.appendChild(featureElem);

    featureElem.setAttribute("identity", string2QString(stdFt->name));

    std::string featureAtoms = "";
    for (unsigned int x = 0; x < stdFt->atoms.size(); x++) {
      featureAtoms+= (int2String(stdFt->atoms[x]) + " ");
    }
    featureElem.setAttribute("atoms", string2QString(featureAtoms));
}

// ============================================================
// Function : writeFuncGroup
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeFuncGroup(QDomDocument doc, QDomElement fragElem, stdFuncGroup* stdFg)
{
    QDomElement funcGroupElem = doc.createElement("funcGroup");
    fragElem.appendChild(funcGroupElem);

    funcGroupElem.setAttribute("group"    , string2QString(stdFg->groupName));
    funcGroupElem.setAttribute("fragment" , string2QString(stdFg->fragName));

    std::string fgroupAtoms = "";
    for (unsigned int x = 0; x < stdFg->atoms.size(); x++) {
      fgroupAtoms += (int2String(stdFg->atoms[x]) + " ");
    }
    funcGroupElem.setAttribute("atoms", string2QString(fgroupAtoms));
}

// ============================================================
// Function : writeConnPts
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeConnPts(QDomDocument doc, QDomElement fragElem, stdFrag* stdFg)
{
    QDomElement connPtsElem = doc.createElement("connectionPoints");
    fragElem.appendChild(connPtsElem);

    std::vector<int> connPtsList = stdFg->getStdConnPtsList();
    std::string connPtsAtoms = "";
    for (unsigned int x = 0; x < connPtsList.size(); x++) {
      connPtsAtoms+= (int2String(connPtsList[x]) + " ");
    }
    connPtsElem.setAttribute("atoms", string2QString(connPtsAtoms));
}
#endif // USE_QT

    // ----------------------------------- //
    // -     TINYXML WRITE FUNCTIONS     - //
    // ----------------------------------- //
#ifdef USE_TINYXML
// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write xml files using tinyxml
// ------------------------------------------------------------
void stdLibParser::Write(std::string fileName) // TinyXML
{
    std::string errMessage = " Writing " + fileName;
    errorLogger.throwError("stdLibParser", errMessage, INFO);

    // Write library file using DOM into memory
    //QString mtkppInfo = PACKAGE_TARNAME;
    //QDomDocument doc(mtkppInfo);

    TiXmlDocument doc;
    TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "", "" );
    doc.LinkEndChild( decl );

    TiXmlElement* root = new TiXmlElement( "stdLib" );
    doc.LinkEndChild( root );

    // stdGroup's
    std::vector<stdGroup*> groupList = pStdLibrary->getStdGroupList();
    for (unsigned int i = 0; i < groupList.size(); i++) {

      TiXmlElement* group = new TiXmlElement("group");
      root->LinkEndChild(group);

      group->SetAttribute("identity", groupList[i]->getName());

      // stdFrag's
      std::vector<stdFrag*> fragList = groupList[i]->getStdFragList();
      for (unsigned int j = 0; j < fragList.size(); j++) {
        this->writeFrag(group, fragList[j]);
      } // fragment

    } // group

    doc.SaveFile(fileName);

    //dump_to_stdout(&doc);
}

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write xml files using tinyxml
// ------------------------------------------------------------
void stdLibParser::Write(std::string fileName, std::string groupName)
{
    std::string errMessage = " Writing " + fileName;
    errorLogger.throwError("stdLibParser", errMessage, INFO);

    // Write library file using DOM into memory
    //QString mtkppInfo = PACKAGE_TARNAME;
    //QDomDocument doc(mtkppInfo);

    TiXmlDocument doc;
    TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "", "" );
    doc.LinkEndChild( decl );

    TiXmlElement* root = new TiXmlElement( "stdLib" );
    doc.LinkEndChild(root);

    // stdGroup's
    std::vector<stdGroup*> groupList = pStdLibrary->getStdGroupList();
    for (unsigned int i = 0; i < groupList.size(); i++) {
      if (groupName == groupList[i]->getName()) {
        this->writeGroup(root, groupList[i]);
      }
    }

    doc.SaveFile(fileName);

    //dump_to_stdout(&doc);
}

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write xml files using tinyxml
// ------------------------------------------------------------
void stdLibParser::Write(std::string fileName, std::string groupName, std::string fragName)
{
    std::string errMessage = " Writing " + fileName;
    errorLogger.throwError("stdLibParser", errMessage, INFO);

    // Write library file using DOM into memory
    //QString mtkppInfo = PACKAGE_TARNAME;
    //QDomDocument doc(mtkppInfo);

    TiXmlDocument doc;
    TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "", "" );
    doc.LinkEndChild( decl );

    TiXmlElement* root = new TiXmlElement( "stdLib" );
    doc.LinkEndChild(root);

    // stdGroup's
    std::vector<stdGroup*> groupList = pStdLibrary->getStdGroupList();
    for (unsigned int i = 0; i < groupList.size(); i++) {
      if (groupName == groupList[i]->getName()) {
        this->writeGroup(root, groupList[i], fragName);
      }
    }

    doc.SaveFile(fileName);

    //dump_to_stdout(&doc);
}

// ============================================================
// Function : writeGroup
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeGroup(TiXmlElement* root, stdGroup* stdGp)
{
    TiXmlElement* group = new TiXmlElement("group");
    root->LinkEndChild(group);
    group->SetAttribute("identity", stdGp->getName());

    // stdFrag's
    std::vector<stdFrag*> fragList = stdGp->getStdFragList();
    for (unsigned int j = 0; j < fragList.size(); j++) {
      this->writeFrag(group, fragList[j]);
    }

    // standard structure
    if (stdGp->hasStdMolecule()) {
      mtkppParser* pMTKppParser = new mtkppParser;
      pMTKppParser->writeMolecule(group, stdGp->getStdMolecule());
    }
}

// ============================================================
// Function : writeGroup
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeGroup(TiXmlElement* root, stdGroup* stdGp, std::string frag3L)
{
    TiXmlElement* group = new TiXmlElement("group");
    root->LinkEndChild(group);
    group->SetAttribute("identity", stdGp->getName());

    // stdFrag's
    std::vector<stdFrag*> fragList = stdGp->getStdFragList();
    for (unsigned int j = 0; j < fragList.size(); j++) {
      if (fragList[j]->getSymbol() == frag3L) {
        this->writeFrag(group, fragList[j]);
      }
    }
}

// ============================================================
// Function : writeFrag
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeFrag(TiXmlElement* root, stdFrag* stdFg)
{
    TiXmlElement* frag = new TiXmlElement("fragment");
    root->LinkEndChild(frag);

    frag->SetAttribute("identity", stdFg->getName());
    frag->SetAttribute("symbol" , stdFg->getSymbol());
    frag->SetAttribute("code"   , stdFg->getCode());
    frag->SetAttribute("type"   , stdFg->getType());

    std::string sGstr = stdFg->getSubGraphStr();
    if (sGstr != "") {
      frag->SetAttribute("subGraphs", sGstr);
    }

    // stdAtom's
    std::vector<stdAtom*> atomList = stdFg->getStdAtomList();
    for (unsigned int k = 0; k < atomList.size(); k++) {
      TiXmlElement* atomElem = new TiXmlElement("atom");
      frag->LinkEndChild(atomElem);

      atomElem->SetAttribute("identity"   , atomList[k]->identity);
      atomElem->SetAttribute("index"      , int2String(atomList[k]->index));
      atomElem->SetAttribute("type"       , atomList[k]->type);
      atomElem->SetAttribute("chain"      , atomList[k]->chain);
      atomElem->SetAttribute("atmCharge"  , double2String(atomList[k]->atmCharge));
      atomElem->SetAttribute("bond12"     , int2String(atomList[k]->bond12));
      atomElem->SetAttribute("bondLength" , double2String(atomList[k]->bondLength));
      atomElem->SetAttribute("bond13"     , int2String(atomList[k]->bond13));
      atomElem->SetAttribute("bondAngle"  , double2String(atomList[k]->bondAngle));
      atomElem->SetAttribute("bond14"     , int2String(atomList[k]->bond14));
      atomElem->SetAttribute("bondTorsion", double2String(atomList[k]->bondTorsion));

      atomElem->SetAttribute("kind", int2String(atomList[k]->kind));

      pStdBond = stdFg->getStdBond(atomList[k]->index,atomList[k]->bond12);
      if (pStdBond) {
        atomElem->SetAttribute("bType", int2String(pStdBond->type));
        atomElem->SetAttribute("bTop", int2String(pStdBond->topology));
        atomElem->SetAttribute("bKind", int2String(pStdBond->kind));
        atomElem->SetAttribute("bondLength", double2String(pStdBond->length));
      }
      else {
        atomElem->SetAttribute("bType", int2String(1));
        atomElem->SetAttribute("bTop", int2String(0));
        atomElem->SetAttribute("bKind", int2String(0));
      }
    }

    // stdLoop's
    std::vector<stdLoop*> loopList = stdFg->getStdLoopList();
    for (unsigned int k = 0; k < loopList.size(); k++) {
      TiXmlElement* loopElem = new TiXmlElement("loop");
      frag->LinkEndChild(loopElem);

      loopElem->SetAttribute("Atom1"  , int2String(loopList[k]->atom1));
      loopElem->SetAttribute("Atom2"  , int2String(loopList[k]->atom2));
      loopElem->SetAttribute("bType"  , int2String(loopList[k]->type));
      loopElem->SetAttribute("bStereo", int2String(loopList[k]->stereo));
    }

    // stdAlias's
    std::vector<stdAlias*> aliasList = stdFg->getStdAliasList();
    for (unsigned int k = 0; k < aliasList.size(); k++) {
      TiXmlElement* aliasElem = new TiXmlElement("alias");
      frag->LinkEndChild(aliasElem);

      aliasElem->SetAttribute("original", aliasList[k]->atom1);
      aliasElem->SetAttribute("alias"   , aliasList[k]->atom2);
    }

    // stdImproper's
    std::vector<stdImproper*> improperList = stdFg->getStdImproperList();
    for (unsigned int k = 0; k < improperList.size(); k++) {
      TiXmlElement* impElem = new TiXmlElement("improper");
      frag->LinkEndChild(impElem);

      /*
      if (improperList[k]->atom1 == -1) {
        impElem->SetAttribute("Atom1", "p1");
      }
      else if (improperList[k]->atom1 == -4) {
        impElem->SetAttribute("Atom1", "n1");
      }
      else {
        impElem->SetAttribute("Atom1", int2String(improperList[k]->atom1));
      }*/
      impElem->SetAttribute("Atom1", int2String(improperList[k]->atom1));
      impElem->SetAttribute("Atom2", int2String(improperList[k]->atom2));
      impElem->SetAttribute("Atom3", int2String(improperList[k]->atom3));
      impElem->SetAttribute("Atom4", int2String(improperList[k]->atom4));
    }

    // stdRing's
    std::vector<stdRing*> ringList = stdFg->getStdRingList();
    for (unsigned int k = 0; k < ringList.size(); k++) {
      TiXmlElement* ringElem = new TiXmlElement("ring");
      frag->LinkEndChild(ringElem);

      ringElem->SetAttribute("size", int2String(ringList[k]->size));
      ringElem->SetAttribute("planar", int2String(ringList[k]->planar));
      ringElem->SetAttribute("aromatic", int2String(ringList[k]->aromatic));
      ringElem->SetAttribute("hetero", int2String(ringList[k]->hetero));
      ringElem->SetAttribute("nHetero", int2String(ringList[k]->nHetero));
      ringElem->SetAttribute("nNitrogen", int2String(ringList[k]->nNitrogen));
      ringElem->SetAttribute("nOxygen", int2String(ringList[k]->nOxygen));
      ringElem->SetAttribute("nSulfur", int2String(ringList[k]->nSulfur));

      std::string ringAtoms = "";
      for (unsigned int x = 0; x < ringList[k]->atoms.size(); x++) {
        ringAtoms+= (int2String(ringList[k]->atoms[x]) + " ");
      }
      ringElem->SetAttribute("atoms", ringAtoms);
    }

    // stdFeature's
    std::vector<stdFeature*> featureList = stdFg->getStdFeatureList();
    for (unsigned int k = 0; k < featureList.size(); k++) {
      TiXmlElement* featureElem = new TiXmlElement("feature");
      frag->LinkEndChild(featureElem);

      featureElem->SetAttribute("identity", featureList[k]->name);

      std::string featureAtoms = "";
      for (unsigned int x = 0; x < featureList[k]->atoms.size(); x++) {
        featureAtoms+= (int2String(featureList[k]->atoms[x]) + " ");
      }
      featureElem->SetAttribute("atoms", featureAtoms);
    }

    // stdFuncGroup's
    std::vector<stdFuncGroup*> funcGroupList = stdFg->getStdFuncGroupList();
    for (unsigned int k = 0; k < funcGroupList.size(); k++) {
      TiXmlElement* funcGroupElem = new TiXmlElement("funcGroup");
      frag->LinkEndChild(funcGroupElem);

      funcGroupElem->SetAttribute("group"    , funcGroupList[k]->groupName);
      funcGroupElem->SetAttribute("fragment" , funcGroupList[k]->fragName);

      std::string fgroupAtoms = "";
      for (unsigned int x = 0; x < funcGroupList[k]->atoms.size(); x++) {
        fgroupAtoms += (int2String(funcGroupList[k]->atoms[x]) + " ");
      }
      funcGroupElem->SetAttribute("atoms", fgroupAtoms);
    }

    // stdConnPts
    std::vector<int> connPtsList = stdFg->getStdConnPtsList();
    if (connPtsList.size() > 0) {
      TiXmlElement* connPtsElem = new TiXmlElement("connectionPoints");
      frag->LinkEndChild(connPtsElem);

      std::string connPtsAtoms = "";
      for (unsigned int x = 0; x < connPtsList.size(); x++) {
        connPtsAtoms+= (int2String(connPtsList[x]) + " ");
      }
      connPtsElem->SetAttribute("atoms", connPtsAtoms);
    }
}

//////////////////////////////////////////////////////////////////////////////

#endif // USE_TINYXML

    // ----------------------------------- //
    // -    XERCES-C WRITE FUNCTIONS     - //
    // ----------------------------------- //
#ifdef USE_XERCES
// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write xml files using xercesc
// ------------------------------------------------------------
void stdLibParser::Write(std::string fileName, std::string groupName)
{
    char* xmlFile = (char*)(fileName.c_str());

    try {
      XMLPlatformUtils::Initialize();
    }
    catch (const XMLException& toCatch) {
      char* message = XMLString::transcode(toCatch.getMessage());
      std::cout << "Error during initialization! :\n"
                << message << "\n";
      XMLString::release(&message);
      setError(1);
      return;
    }

    int errorCode = 0;

    XMLCh tempStr[100];
    XMLString::transcode("LS", tempStr, 99);
    DOMImplementation *impl = DOMImplementationRegistry::getDOMImplementation(tempStr);

    // xml settings
    bool gDoNamespaces          = false;
    bool gDoSchema              = false;
    bool gSchemaFullChecking    = false;
    bool gDoCreate              = false;

#if (XERCES_VERSION_MAJOR == 3)
#else
    bool gSplitCdataSections    = true;
    bool gDiscardDefaultContent = true;
    bool gFormatPrettyPrint     = true;
    bool gWriteBOM              = false;
#endif

    XercesDOMParser::ValSchemes gValScheme  = XercesDOMParser::Val_Auto;

    //  Create our parser, then attach an error handler to the parser.
    //  The parser will call back to methods of the ErrorHandler if it
    //  discovers errors during the course of parsing the XML document.
    XercesDOMParser *parser = new XercesDOMParser;
    parser->setValidationScheme(gValScheme);
    parser->setDoNamespaces(gDoNamespaces);
    parser->setDoSchema(gDoSchema);
    parser->setValidationSchemaFullChecking(gSchemaFullChecking);
    parser->setCreateEntityReferenceNodes(gDoCreate);

    //DOMTreeErrorReporter *errReporter = new DOMTreeErrorReporter();
    //parser->setErrorHandler(errReporter);

    {
      if (impl != 0) {
        try {
#if (XERCES_VERSION_MAJOR == 3)
          DOMLSSerializer         *theSerializer = ((DOMImplementationLS*)impl)->createLSSerializer(XMLPlatformUtils::fgMemoryManager);
#else
          DOMWriter         *theSerializer = ((DOMImplementationLS*)impl)->createDOMWriter();
#endif

          //DOMErrorHandler *myErrorHandler = new DOMPrintErrorHandler();
          //theSerializer->setErrorHandler(myErrorHandler);

          XERCES_CPP_NAMESPACE::DOMDocument* doc = impl->createDocument(
            0,           // root element namespace URI.
            X("stdLib"), // root element name
            0);          // document type object (DTD).

          std::vector<stdGroup*> groupList = pStdLibrary->getStdGroupList();
          for (unsigned int i = 0; i < groupList.size(); i++) {
            if (groupName == groupList[i]->getName()) {
              this->writeGroup(doc, groupList[i]);
            }
          }

#if (XERCES_VERSION_MAJOR == 3)
#else
          // set feature if the serializer supports the feature/mode
          if (theSerializer->canSetFeature(XMLUni::fgDOMWRTSplitCdataSections, gSplitCdataSections))
              theSerializer->setFeature(XMLUni::fgDOMWRTSplitCdataSections, gSplitCdataSections);

          if (theSerializer->canSetFeature(XMLUni::fgDOMWRTDiscardDefaultContent, gDiscardDefaultContent))
              theSerializer->setFeature(XMLUni::fgDOMWRTDiscardDefaultContent, gDiscardDefaultContent);

          if (theSerializer->canSetFeature(XMLUni::fgDOMWRTFormatPrettyPrint, gFormatPrettyPrint))
              theSerializer->setFeature(XMLUni::fgDOMWRTFormatPrettyPrint, gFormatPrettyPrint);

          if (theSerializer->canSetFeature(XMLUni::fgDOMWRTBOM, gWriteBOM))
              theSerializer->setFeature(XMLUni::fgDOMWRTBOM, gWriteBOM);
#endif

          XMLFormatTarget *myFormTarget;
          myFormTarget = new LocalFileFormatTarget(xmlFile);

#if (XERCES_VERSION_MAJOR == 3)
          DOMLSOutput         *lsOutput = ((DOMImplementationLS*)impl)->createLSOutput(XMLPlatformUtils::fgMemoryManager);
          lsOutput->setByteStream(myFormTarget);

          // do the serialization to xml file.
          theSerializer->write(doc, lsOutput);

#else
          // do the serialization through DOMWriter::writeNode();
          theSerializer->writeNode(myFormTarget, *doc);
#endif

          delete theSerializer;
          delete myFormTarget;
          //delete myErrorHandler;

          doc->release();
        }
        catch (const OutOfMemoryException&) {
          XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
          errorCode = 5;
        }
        catch (const DOMException& e) {
          XERCES_STD_QUALIFIER cerr << "DOMException code is:  " << e.code << XERCES_STD_QUALIFIER endl;
          errorCode = 2;
        }
        catch (...) {
          XERCES_STD_QUALIFIER cerr << "An error occurred creating the document" << XERCES_STD_QUALIFIER endl;
          std::cout << " stdLibParser::Write(" << fileName << ", " << groupName << std::endl;
          errorCode = 3;
        }
      }
      else {
        XERCES_STD_QUALIFIER cerr << "Requested implementation is not supported" << XERCES_STD_QUALIFIER endl;
        errorCode = 4;
      }
    }

    XMLPlatformUtils::Terminate();
    setError(errorCode);
    return;
}



// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write xml files using xercesc
// ------------------------------------------------------------
void stdLibParser::Write(std::string fileName, std::string groupName, std::string fragName)
{
    char* xmlFile = (char*)(fileName.c_str());

    try {
      XMLPlatformUtils::Initialize();
    }
    catch (const XMLException& toCatch) {
      char* message = XMLString::transcode(toCatch.getMessage());
      std::cout << "Error during initialization! :\n"
                << message << "\n";
      XMLString::release(&message);
      setError(1);
      return;
    }

    int errorCode = 0;

    XMLCh tempStr[100];
    XMLString::transcode("LS", tempStr, 99);
    DOMImplementation *impl = DOMImplementationRegistry::getDOMImplementation(tempStr);

    // xml settings
    bool gDoNamespaces          = false;
    bool gDoSchema              = false;
    bool gSchemaFullChecking    = false;
    bool gDoCreate              = false;

#if (XERCES_VERSION_MAJOR == 3)
#else
    bool gSplitCdataSections    = true;
    bool gDiscardDefaultContent = true;
    bool gFormatPrettyPrint     = true;
    bool gWriteBOM              = false;
#endif

    XercesDOMParser::ValSchemes gValScheme  = XercesDOMParser::Val_Auto;

    //  Create our parser, then attach an error handler to the parser.
    //  The parser will call back to methods of the ErrorHandler if it
    //  discovers errors during the course of parsing the XML document.
    XercesDOMParser *parser = new XercesDOMParser;
    parser->setValidationScheme(gValScheme);
    parser->setDoNamespaces(gDoNamespaces);
    parser->setDoSchema(gDoSchema);
    parser->setValidationSchemaFullChecking(gSchemaFullChecking);
    parser->setCreateEntityReferenceNodes(gDoCreate);

    //DOMTreeErrorReporter *errReporter = new DOMTreeErrorReporter();
    //parser->setErrorHandler(errReporter);

    {
      if (impl != 0) {
        try {
#if (XERCES_VERSION_MAJOR == 3)
          DOMLSSerializer         *theSerializer = ((DOMImplementationLS*)impl)->createLSSerializer(XMLPlatformUtils::fgMemoryManager);
#else
          DOMWriter         *theSerializer = ((DOMImplementationLS*)impl)->createDOMWriter();
#endif

          //DOMErrorHandler *myErrorHandler = new DOMPrintErrorHandler();
          //theSerializer->setErrorHandler(myErrorHandler);

          XERCES_CPP_NAMESPACE::DOMDocument* doc = impl->createDocument(
            0,           // root element namespace URI.
            X("stdLib"), // root element name
            0);          // document type object (DTD).

          std::vector<stdGroup*> groupList = pStdLibrary->getStdGroupList();
          for (unsigned int i = 0; i < groupList.size(); i++) {
            if (groupName == groupList[i]->getName()) {
              //std::cout << "    stdLibParser::Write Group " << groupList[i] << std::endl;
              this->writeGroup(doc, groupList[i], fragName);
            }
          }

#if (XERCES_VERSION_MAJOR == 3)
#else
          // set feature if the serializer supports the feature/mode
          if (theSerializer->canSetFeature(XMLUni::fgDOMWRTSplitCdataSections, gSplitCdataSections))
              theSerializer->setFeature(XMLUni::fgDOMWRTSplitCdataSections, gSplitCdataSections);

          if (theSerializer->canSetFeature(XMLUni::fgDOMWRTDiscardDefaultContent, gDiscardDefaultContent))
              theSerializer->setFeature(XMLUni::fgDOMWRTDiscardDefaultContent, gDiscardDefaultContent);

          if (theSerializer->canSetFeature(XMLUni::fgDOMWRTFormatPrettyPrint, gFormatPrettyPrint))
              theSerializer->setFeature(XMLUni::fgDOMWRTFormatPrettyPrint, gFormatPrettyPrint);

          if (theSerializer->canSetFeature(XMLUni::fgDOMWRTBOM, gWriteBOM))
              theSerializer->setFeature(XMLUni::fgDOMWRTBOM, gWriteBOM);
#endif
          XMLFormatTarget *myFormTarget;
          myFormTarget = new LocalFileFormatTarget(xmlFile);

#if (XERCES_VERSION_MAJOR == 3)
          DOMLSOutput         *lsOutput = ((DOMImplementationLS*)impl)->createLSOutput(XMLPlatformUtils::fgMemoryManager);
          lsOutput->setByteStream(myFormTarget);

          // do the serialization to xml file.
          theSerializer->write(doc, lsOutput);

#else
          // do the serialization through DOMWriter::writeNode();
          theSerializer->writeNode(myFormTarget, *doc);
#endif

          delete theSerializer;
          delete myFormTarget;
          //delete myErrorHandler;

          doc->release();
        }
        catch (const OutOfMemoryException&) {
          XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
          errorCode = 5;
        }
        catch (const DOMException& e) {
          XERCES_STD_QUALIFIER cerr << "DOMException code is:  " << e.code << XERCES_STD_QUALIFIER endl;
          errorCode = 2;
        }
        catch (...) {
          XERCES_STD_QUALIFIER cerr << "An error occurred creating the document" << XERCES_STD_QUALIFIER endl;
          std::cout << " stdLibParser::Write(" << fileName << ", " << groupName << std::endl;
          errorCode = 3;
        }
      }
      else {
        XERCES_STD_QUALIFIER cerr << "Requested implementation is not supported" << XERCES_STD_QUALIFIER endl;
        errorCode = 4;
      }
    }

    XMLPlatformUtils::Terminate();
    setError(errorCode);
    return;
}

// ============================================================
// Function : Write
// ------------------------------------------------------------
// Write xml files using xercesc
// ------------------------------------------------------------
void stdLibParser::Write(std::string fileName)
{
    char* xmlFile = (char*)(fileName.c_str());

    try {
      XMLPlatformUtils::Initialize();
    }
    catch (const XMLException& toCatch) {
      char* message = XMLString::transcode(toCatch.getMessage());
      std::cout << "Error during initialization! :\n"
                << message << "\n";
      XMLString::release(&message);
      setError(1);
      return;
    }

    int errorCode = 0;

    XMLCh tempStr[100];
    XMLString::transcode("LS", tempStr, 99);
    DOMImplementation *impl = DOMImplementationRegistry::getDOMImplementation(tempStr);

    // xml settings
    bool gDoNamespaces          = false;
    bool gDoSchema              = false;
    bool gSchemaFullChecking    = false;
    bool gDoCreate              = false;

#if (XERCES_VERSION_MAJOR == 3)
#else
    bool gSplitCdataSections    = true;
    bool gDiscardDefaultContent = true;
    bool gFormatPrettyPrint     = true;
    bool gWriteBOM              = false;
#endif

    XercesDOMParser::ValSchemes gValScheme  = XercesDOMParser::Val_Auto;

    //  Create our parser, then attach an error handler to the parser.
    //  The parser will call back to methods of the ErrorHandler if it
    //  discovers errors during the course of parsing the XML document.
    XercesDOMParser *parser = new XercesDOMParser;
    parser->setValidationScheme(gValScheme);
    parser->setDoNamespaces(gDoNamespaces);
    parser->setDoSchema(gDoSchema);
    parser->setValidationSchemaFullChecking(gSchemaFullChecking);
    parser->setCreateEntityReferenceNodes(gDoCreate);

    //DOMTreeErrorReporter *errReporter = new DOMTreeErrorReporter();
    //parser->setErrorHandler(errReporter);

    {
      if (impl != 0) {
        try {
#if (XERCES_VERSION_MAJOR == 3)
          DOMLSSerializer *theSerializer = ((DOMImplementationLS*)impl)->createLSSerializer(XMLPlatformUtils::fgMemoryManager);
#else
          DOMWriter *theSerializer = ((DOMImplementationLS*)impl)->createDOMWriter();
#endif

          //DOMErrorHandler *myErrorHandler = new DOMPrintErrorHandler();
          //theSerializer->setErrorHandler(myErrorHandler);

          XERCES_CPP_NAMESPACE::DOMDocument* doc = impl->createDocument(
            0,           // root element namespace URI.
            X("stdLib"), // root element name
            0);          // document type object (DTD).

          XERCES_CPP_NAMESPACE::DOMElement* rootElem = doc->getDocumentElement();

          // stdGroup's
          std::vector<stdGroup*> groupList = pStdLibrary->getStdGroupList();
          for (unsigned int i = 0; i < groupList.size(); i++) {
            //std::cout << groupList[i]->getName() << endl;

            DOMElement*  groupElem = doc->createElement(X("group"));
            rootElem->appendChild(groupElem);

            groupElem->setAttribute(X("identity"), X(groupList[i]->getName().c_str()));

            // stdFrag's
            std::vector<stdFrag*> fragList = groupList[i]->getStdFragList();
            for (unsigned int j = 0; j < fragList.size(); j++) {
              //std::cout << fragList[j]->getName() << endl;
              DOMElement*  fragElem = doc->createElement(X("fragment"));
              groupElem->appendChild(fragElem);

              fragElem->setAttribute(X("identity"), X(fragList[j]->getName().c_str()));
              fragElem->setAttribute(X("symbol"), X(fragList[j]->getSymbol().c_str()));
              fragElem->setAttribute(X("character"), X(fragList[j]->getCharacter().c_str()));
              fragElem->setAttribute(X("code"), X(fragList[j]->getCode().c_str()));
              fragElem->setAttribute(X("type"), X(fragList[j]->getType().c_str()));

              std::string sGstr = fragList[j]->getSubGraphStr();
              if (sGstr != "") fragElem->setAttribute(X("subGraphs"), X(sGstr.c_str()));

              // stdAtom's
              std::vector<stdAtom*> atomList = fragList[j]->getStdAtomList();
              for (unsigned int k = 0; k < atomList.size(); k++) {
                DOMElement*  atomElem = doc->createElement(X("atom"));
                fragElem->appendChild(atomElem);

                atomElem->setAttribute(X("identity"), X(atomList[k]->identity.c_str()));
                atomElem->setAttribute(X("index"), X(int2String(atomList[k]->index).c_str()));
                atomElem->setAttribute(X("type"), X(atomList[k]->type.c_str()));
                atomElem->setAttribute(X("chain"), X(atomList[k]->chain.c_str()));
                atomElem->setAttribute(X("atmCharge"), X(double2String(atomList[k]->atmCharge).c_str())); //18.2223
                atomElem->setAttribute(X("bond12"), X(int2String(atomList[k]->bond12).c_str()));
                atomElem->setAttribute(X("bondLength"), X(double2String(atomList[k]->bondLength).c_str()));
                atomElem->setAttribute(X("bond13"), X(int2String(atomList[k]->bond13).c_str()));
                atomElem->setAttribute(X("bondAngle"), X(double2String(atomList[k]->bondAngle).c_str()));
                atomElem->setAttribute(X("bond14"), X(int2String(atomList[k]->bond14).c_str()));
                atomElem->setAttribute(X("bondTorsion"), X(double2String(atomList[k]->bondTorsion).c_str()));
                //atomElem->setAttribute(X("atNum"), X(int2String(atomList[k]->atNum).c_str()));
                //atomElem->setAttribute(X("symbol"), X(atomList[k]->atSymbol.c_str()));
                //atomElem->setAttribute(X("hybridization"), X(atomList[k]->hybridization.c_str()));
                atomElem->setAttribute(X("kind"), X(int2String(atomList[k]->kind).c_str()));
                pStdBond = fragList[j]->getStdBond(atomList[k]->index,atomList[k]->bond12);
                if (pStdBond) {
                  atomElem->setAttribute(X("bType"), X(int2String(pStdBond->type).c_str()));
                  atomElem->setAttribute(X("bTop"), X(int2String(pStdBond->topology).c_str()));
                  atomElem->setAttribute(X("bKind"), X(int2String(pStdBond->kind).c_str()));
                  atomElem->setAttribute(X("bondLength"), X(double2String(pStdBond->length).c_str()));
                }
                else {
                  atomElem->setAttribute(X("bType"), X(int2String(1).c_str()));
                  atomElem->setAttribute(X("bTop"), X(int2String(0).c_str()));
                  atomElem->setAttribute(X("bKind"), X(int2String(0).c_str()));
                }
              }

              // stdLoop's
              std::vector<stdLoop*> loopList = fragList[j]->getStdLoopList();
              for (unsigned int k = 0; k < loopList.size(); k++) {
                DOMElement*  loopElem = doc->createElement(X("loop"));
                fragElem->appendChild(loopElem);

                loopElem->setAttribute(X("Atom1"), X(int2String(loopList[k]->atom1).c_str()));
                loopElem->setAttribute(X("Atom2"), X(int2String(loopList[k]->atom2).c_str()));
                loopElem->setAttribute(X("bType"), X(int2String(loopList[k]->type).c_str()));
                loopElem->setAttribute(X("bStereo"), X(int2String(loopList[k]->stereo).c_str()));
              }

              // stdAlias's
              std::vector<stdAlias*> aliasList = fragList[j]->getStdAliasList();
              for (unsigned int k = 0; k < aliasList.size(); k++) {
                DOMElement*  aliasElem = doc->createElement(X("alias"));
                fragElem->appendChild(aliasElem);

                aliasElem->setAttribute(X("original"), X( aliasList[k]->atom1.c_str()));
                aliasElem->setAttribute(X("alias"),    X( aliasList[k]->atom2.c_str()));
              }

              // stdImproper's
              std::vector<stdImproper*> improperList = fragList[j]->getStdImproperList();
              for (unsigned int k = 0; k < improperList.size(); k++) {
                DOMElement*  impElem = doc->createElement(X("improper"));
                fragElem->appendChild(impElem);

                if (improperList[k]->atom1 == -1) {
                  impElem->setAttribute(X("Atom1"), X("p1"));
                }
                else if (improperList[k]->atom1 == -4) {
                  impElem->setAttribute(X("Atom1"), X("n1"));
                }
                else {
                  impElem->setAttribute(X("Atom1"), X(int2String(improperList[k]->atom1).c_str()));
                }
                impElem->setAttribute(X("Atom2"), X(int2String(improperList[k]->atom2).c_str()));
                impElem->setAttribute(X("Atom3"), X(int2String(improperList[k]->atom3).c_str()));
                impElem->setAttribute(X("Atom4"), X(int2String(improperList[k]->atom4).c_str()));
              }

              // stdRing's
              std::vector<stdRing*> ringList = fragList[j]->getStdRingList();
              for (unsigned int k = 0; k < ringList.size(); k++) {
                DOMElement*  ringElem = doc->createElement(X("ring"));
                fragElem->appendChild(ringElem);

                ringElem->setAttribute(X("size"),      X(int2String(ringList[k]->size).c_str()));
                ringElem->setAttribute(X("planar"),    X(int2String(ringList[k]->planar).c_str()));
                ringElem->setAttribute(X("aromatic"),  X(int2String(ringList[k]->aromatic).c_str()));
                ringElem->setAttribute(X("hetero"),    X(int2String(ringList[k]->hetero).c_str()));
                ringElem->setAttribute(X("nHetero"),   X(int2String(ringList[k]->nHetero).c_str()));
                ringElem->setAttribute(X("nNitrogen"), X(int2String(ringList[k]->nNitrogen).c_str()));
                ringElem->setAttribute(X("nOxygen"),   X(int2String(ringList[k]->nOxygen).c_str()));
                ringElem->setAttribute(X("nSulfur"),   X(int2String(ringList[k]->nSulfur).c_str()));

                std::string ringAtoms = "";
                for (unsigned int x = 0; x < ringList[k]->atoms.size(); x++) {
                  ringAtoms+= (int2String(ringList[k]->atoms[x]) + " ");
                }
                ringElem->setAttribute(X("atoms"), X(ringAtoms.c_str()));
              }

              // stdFeature's
              std::vector<stdFeature*> featureList = fragList[j]->getStdFeatureList();
              for (unsigned int k = 0; k < featureList.size(); k++) {
                DOMElement*  featureElem = doc->createElement(X("feature"));
                fragElem->appendChild(featureElem);

                featureElem->setAttribute(X("identity"), X(featureList[k]->name.c_str()));

                std::string featureAtoms = "";
                for (unsigned int x = 0; x < featureList[k]->atoms.size(); x++) {
                  featureAtoms+= (int2String(featureList[k]->atoms[x]) + " ");
                }
                featureElem->setAttribute(X("atoms"), X(featureAtoms.c_str()));
              }

              // stdFuncGroup's
              std::vector<stdFuncGroup*> funcGroupList = fragList[j]->getStdFuncGroupList();
              for (unsigned int k = 0; k < funcGroupList.size(); k++) {
                DOMElement* funcGroupElem = doc->createElement(X("funcGroup"));
                fragElem->appendChild(funcGroupElem);

                funcGroupElem->setAttribute(X("group"), X(funcGroupList[k]->groupName.c_str()));
                funcGroupElem->setAttribute(X("fragment"), X(funcGroupList[k]->fragName.c_str()));

                std::string fgroupAtoms = "";
                for (unsigned int x = 0; x < funcGroupList[k]->atoms.size(); x++) {
                  fgroupAtoms += (int2String(funcGroupList[k]->atoms[x]) + " ");
                }
                funcGroupElem->setAttribute(X("atoms"), X(fgroupAtoms.c_str()));
              }

              // stdConnPts
              std::vector<int> connPtsList = fragList[j]->getStdConnPtsList();
              if (connPtsList.size() > 0) {
                DOMElement*  connPtsElem = doc->createElement(X("connectionPoints"));
                fragElem->appendChild(connPtsElem);

                std::string connPtsAtoms = "";
                for (unsigned int x = 0; x < connPtsList.size(); x++) {
                  connPtsAtoms+= (int2String(connPtsList[x]) + " ");
                }
                connPtsElem->setAttribute(X("atoms"), X(connPtsAtoms.c_str()));
              }
            }
          }

          // Now count the number of elements in the above DOM tree.
          //unsigned int elementCount = doc->getElementsByTagName(X("*"))->getLength();
          //XERCES_STD_QUALIFIER cout << "The tree just created contains: " << elementCount
          //                          << " elements." << XERCES_STD_QUALIFIER endl;

#if (XERCES_VERSION_MAJOR == 3)
#else
          // set feature if the serializer supports the feature/mode
          if (theSerializer->canSetFeature(XMLUni::fgDOMWRTSplitCdataSections, gSplitCdataSections))
              theSerializer->setFeature(XMLUni::fgDOMWRTSplitCdataSections, gSplitCdataSections);

          if (theSerializer->canSetFeature(XMLUni::fgDOMWRTDiscardDefaultContent, gDiscardDefaultContent))
              theSerializer->setFeature(XMLUni::fgDOMWRTDiscardDefaultContent, gDiscardDefaultContent);

          if (theSerializer->canSetFeature(XMLUni::fgDOMWRTFormatPrettyPrint, gFormatPrettyPrint))
              theSerializer->setFeature(XMLUni::fgDOMWRTFormatPrettyPrint, gFormatPrettyPrint);

          if (theSerializer->canSetFeature(XMLUni::fgDOMWRTBOM, gWriteBOM))
              theSerializer->setFeature(XMLUni::fgDOMWRTBOM, gWriteBOM);
#endif

          XMLFormatTarget *myFormTarget;
          myFormTarget = new LocalFileFormatTarget(xmlFile);

#if (XERCES_VERSION_MAJOR == 3)
          DOMLSOutput         *lsOutput = ((DOMImplementationLS*)impl)->createLSOutput(XMLPlatformUtils::fgMemoryManager);
          lsOutput->setByteStream(myFormTarget);

          // do the serialization to xml file.
          theSerializer->write(doc, lsOutput);

#else
          // do the serialization through DOMWriter::writeNode();
          theSerializer->writeNode(myFormTarget, *doc);
#endif

          delete theSerializer;
          delete myFormTarget;
          //delete myErrorHandler;

          doc->release();
        }
        catch (const OutOfMemoryException&) {
          XERCES_STD_QUALIFIER cerr << "OutOfMemoryException" << XERCES_STD_QUALIFIER endl;
          errorCode = 5;
        }
        catch (const DOMException& e) {
          XERCES_STD_QUALIFIER cerr << "DOMException code is:  " << e.code << XERCES_STD_QUALIFIER endl;
          errorCode = 2;
        }
        catch (...) {
          XERCES_STD_QUALIFIER cerr << "An error occurred creating the document" << XERCES_STD_QUALIFIER endl;
          std::cout << " stdLibParser::Write(" << fileName << ", " << std::endl;

          errorCode = 3;
        }
      }
      else {
        XERCES_STD_QUALIFIER cerr << "Requested implementation is not supported" << XERCES_STD_QUALIFIER endl;
        errorCode = 4;
      }
    }

    XMLPlatformUtils::Terminate();
    setError(errorCode);
    return;
}

// ============================================================
// Function : writeGroup
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeGroup(XERCES_CPP_NAMESPACE::DOMDocument* doc, stdGroup* stdGp)
{
    DOMElement* rootElem = doc->getDocumentElement();
    DOMElement* groupElem = doc->createElement(X("group"));
    rootElem->appendChild(groupElem);
    groupElem->setAttribute(X("identity"), X(stdGp->getName().c_str()));

    // stdFrag's
    std::vector<stdFrag*> fragList = stdGp->getStdFragList();
    for (unsigned int j = 0; j < fragList.size(); j++) {
      this->writeFrag(doc, groupElem, fragList[j]);
    }

    // standard structure
    if (stdGp->hasStdMolecule()) {
      mtkppParser* pMTKppParser = new mtkppParser;
      pMTKppParser->writeMolecule(doc, stdGp->getStdMolecule());
    }
}

// ============================================================
// Function : writeGroup
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeGroup(XERCES_CPP_NAMESPACE::DOMDocument* doc, stdGroup* stdGp, std::string frag3L)
{
    DOMElement* rootElem = doc->getDocumentElement();
    DOMElement* groupElem = doc->createElement(X("group"));
    rootElem->appendChild(groupElem);
    groupElem->setAttribute(X("identity"), X(stdGp->getName().c_str()));

    // stdFrag's
    std::vector<stdFrag*> fragList = stdGp->getStdFragList();
    for (unsigned int j = 0; j < fragList.size(); j++) {
      if (fragList[j]->getSymbol() == frag3L) {
        this->writeFrag(doc, groupElem, fragList[j]);
      }
    }
}

// ============================================================
// Function : writeFrag
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeFrag(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* groupElem, stdFrag* stdFg)
{
    DOMElement* fragElem = doc->createElement(X("fragment"));

    groupElem->appendChild(fragElem);

    fragElem->setAttribute(X("identity"), X(stdFg->getName().c_str()));
    fragElem->setAttribute(X("symbol"), X(stdFg->getSymbol().c_str()));
    fragElem->setAttribute(X("character"), X(stdFg->getCharacter().c_str()));
    fragElem->setAttribute(X("code"), X(stdFg->getCode().c_str()));
    fragElem->setAttribute(X("type"), X(stdFg->getType().c_str()));

    std::string sGstr = stdFg->getSubGraphStr();
    if (sGstr != "") fragElem->setAttribute(X("subGraphs"), X(sGstr.c_str()));

    // stdAtom's
    std::vector<stdAtom*> atomList = stdFg->getStdAtomList();
    for (unsigned int k = 0; k < atomList.size(); k++) {
      this->writeAtom(doc, fragElem, stdFg, atomList[k]);
    }

    // stdLoop's
    std::vector<stdLoop*> loopList = stdFg->getStdLoopList();
    for (unsigned int k = 0; k < loopList.size(); k++) {
      this->writeLoop(doc, fragElem, loopList[k]);
    }

    // stdAlias's
    std::vector<stdAlias*> aliasList = stdFg->getStdAliasList();
    for (unsigned int k = 0; k < aliasList.size(); k++) {
      this->writeAlias(doc, fragElem, aliasList[k]);
    }

    // stdImproper's
    std::vector<stdImproper*> improperList = stdFg->getStdImproperList();
    for (unsigned int k = 0; k < improperList.size(); k++) {
      this->writeImproper(doc, fragElem, improperList[k]);
    }

    // stdRing's
    std::vector<stdRing*> ringList = stdFg->getStdRingList();
    for (unsigned int k = 0; k < ringList.size(); k++) {
      this->writeRing(doc, fragElem, ringList[k]);
    }

    // stdFeature's
    std::vector<stdFeature*> featureList = stdFg->getStdFeatureList();
    for (unsigned int k = 0; k < featureList.size(); k++) {
      this->writeFeature(doc, fragElem, featureList[k]);
    }

    // stdFuncGroup's
    std::vector<stdFuncGroup*> funcGroupList = stdFg->getStdFuncGroupList();
    for (unsigned int k = 0; k < funcGroupList.size(); k++) {
      this->writeFuncGroup(doc, fragElem, funcGroupList[k]);
    }

    // stdConnPts
    this->writeConnPts(doc, fragElem, stdFg);
}

// ============================================================
// Function : writeAtom
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeAtom(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement*  fragElem, stdFrag* stdFg, stdAtom* stdAt)
{
    DOMElement* atomElem = doc->createElement(X("atom"));
    fragElem->appendChild(atomElem);

    atomElem->setAttribute(X("identity"), X(stdAt->identity.c_str()));
    atomElem->setAttribute(X("index"), X(int2String(stdAt->index).c_str()));
    atomElem->setAttribute(X("type"), X(stdAt->type.c_str()));
    atomElem->setAttribute(X("chain"), X(stdAt->chain.c_str()));
    atomElem->setAttribute(X("atmCharge"), X(double2String(stdAt->atmCharge).c_str())); // 18.2223
    atomElem->setAttribute(X("bond12"), X(int2String(stdAt->bond12).c_str()));
    atomElem->setAttribute(X("bondLength"), X(double2String(stdAt->bondLength).c_str()));
    atomElem->setAttribute(X("bond13"), X(int2String(stdAt->bond13).c_str()));
    atomElem->setAttribute(X("bondAngle"), X(double2String(stdAt->bondAngle).c_str()));
    atomElem->setAttribute(X("bond14"), X(int2String(stdAt->bond14).c_str()));
    atomElem->setAttribute(X("bondTorsion"), X(double2String(stdAt->bondTorsion).c_str()));
    //atomElem->setAttribute(X("atNum"), X(int2String(stdAt->atNum).c_str()));
    //atomElem->setAttribute(X("symbol"), X(stdAt->atSymbol.c_str()));
    //atomElem->setAttribute(X("hybridization"), X(stdAt->hybridization.c_str()));
    atomElem->setAttribute(X("kind"), X(int2String(stdAt->kind).c_str()));
    pStdBond = stdFg->getStdBond(stdAt->index,stdAt->bond12);
    if (pStdBond) {
      atomElem->setAttribute(X("bType"), X(int2String(pStdBond->type).c_str()));
      atomElem->setAttribute(X("bTop"), X(int2String(pStdBond->topology).c_str()));
      atomElem->setAttribute(X("bKind"), X(int2String(pStdBond->kind).c_str()));
      atomElem->setAttribute(X("bondLength"), X(double2String(pStdBond->length).c_str()));
    }
    else {
      atomElem->setAttribute(X("bType"), X(int2String(1).c_str()));
      atomElem->setAttribute(X("bTop"), X(int2String(0).c_str()));
      atomElem->setAttribute(X("bKind"), X(int2String(0).c_str()));
    }
}

// ============================================================
// Function : writeLoop
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeLoop(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* fragElem, stdLoop* stdLp)
{
    DOMElement*  loopElem = doc->createElement(X("loop"));
    fragElem->appendChild(loopElem);

    loopElem->setAttribute(X("Atom1"), X(int2String(stdLp->atom1).c_str()));
    loopElem->setAttribute(X("Atom2"), X(int2String(stdLp->atom2).c_str()));
    loopElem->setAttribute(X("bType"), X(int2String(stdLp->type).c_str()));
    loopElem->setAttribute(X("bStereo"), X(int2String(stdLp->stereo).c_str()));
}

// ============================================================
// Function : writeAlias
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeAlias(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* fragElem, stdAlias* stdAl)
{

    DOMElement* aliasElem = doc->createElement(X("alias"));
    fragElem->appendChild(aliasElem);

    aliasElem->setAttribute(X("original"), X( stdAl->atom1.c_str()));
    aliasElem->setAttribute(X("alias"),    X( stdAl->atom2.c_str()));
}

// ============================================================
// Function : writeImproper
// ------------------------------------------------------------
//
// ------------------------------------------------------------
void stdLibParser::writeImproper(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* fragElem, stdImproper* stdIm)
{
    DOMElement* impElem = doc->createElement(X("improper"));
    fragElem->appendChild(impElem);

/*
    if (stdIm->atom1 == -1) {
      impElem->setAttribute(X("Atom1"), X("p1"));
    }
    else if (stdIm->atom1 == -4) {
      impElem->setAttribute(X("Atom1"), X("n1"));
    }
    else {
      impElem->setAttribute(X("Atom1"), X(int2String(stdIm->atom1).c_str()));
    }
*/

    impElem->setAttribute(X("Atom1"), X(int2String(stdIm->atom1).c_str()));
    impElem->setAttribute(X("Atom2"), X(int2String(stdIm->atom2).c_str()));
    impElem->setAttribute(X("Atom3"), X(int2String(stdIm->atom3).c_str()));
    impElem->setAttribute(X("Atom4"), X(int2String(stdIm->atom4).c_str()));
}

// ============================================================
// Function : writeRing
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void stdLibParser::writeRing(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* fragElem, stdRing* stdRg)
{
    DOMElement*  ringElem = doc->createElement(X("ring"));
    fragElem->appendChild(ringElem);

    ringElem->setAttribute(X("size"),      X(int2String(stdRg->size).c_str()));
    ringElem->setAttribute(X("planar"),    X(int2String(stdRg->planar).c_str()));
    ringElem->setAttribute(X("aromatic"),  X(int2String(stdRg->aromatic).c_str()));
    ringElem->setAttribute(X("hetero"),    X(int2String(stdRg->hetero).c_str()));
    ringElem->setAttribute(X("nHetero"),   X(int2String(stdRg->nHetero).c_str()));
    ringElem->setAttribute(X("nNitrogen"), X(int2String(stdRg->nNitrogen).c_str()));
    ringElem->setAttribute(X("nOxygen"),   X(int2String(stdRg->nOxygen).c_str()));
    ringElem->setAttribute(X("nSulfur"),   X(int2String(stdRg->nSulfur).c_str()));

    std::string ringAtoms = "";
    for (unsigned int x = 0; x < stdRg->atoms.size(); x++) {
      ringAtoms += (int2String(stdRg->atoms[x]) + " ");
    }
    ringElem->setAttribute(X("atoms"), X(ringAtoms.c_str()));
}

// ============================================================
// Function : writeFeature
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void stdLibParser::writeFeature(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* fragElem, stdFeature* stdFt)
{
    DOMElement*  featureElem = doc->createElement(X("feature"));
    fragElem->appendChild(featureElem);

    featureElem->setAttribute(X("identity"), X(stdFt->name.c_str()));

    std::string featureAtoms = "";
    for (unsigned int x = 0; x < stdFt->atoms.size(); x++) {
      featureAtoms+= (int2String(stdFt->atoms[x]) + " ");
    }
    featureElem->setAttribute(X("atoms"), X(featureAtoms.c_str()));
}

// ============================================================
// Function : writeFuncGroup
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void stdLibParser::writeFuncGroup(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* fragElem, stdFuncGroup* stdFg)
{
    DOMElement* funcGroupElem = doc->createElement(X("funcGroup"));
    fragElem->appendChild(funcGroupElem);

    funcGroupElem->setAttribute(X("group"), X(stdFg->groupName.c_str()));
    funcGroupElem->setAttribute(X("fragment"), X(stdFg->fragName.c_str()));

    std::string fgroupAtoms = "";
    for (unsigned int x = 0; x < stdFg->atoms.size(); x++) {
      fgroupAtoms+= (int2String(stdFg->atoms[x]) + " ");
    }
    funcGroupElem->setAttribute(X("atoms"), X(fgroupAtoms.c_str()));
}

// ============================================================
// Function : writeConnPts
// ------------------------------------------------------------
// 
// ------------------------------------------------------------
void stdLibParser::writeConnPts(XERCES_CPP_NAMESPACE::DOMDocument* doc, DOMElement* fragElem, stdFrag* stdFg)
{
    std::vector<int> connPtsList = stdFg->getStdConnPtsList();
    if (connPtsList.size() > 0) {
      DOMElement*  connPtsElem = doc->createElement(X("connectionPoints"));
      fragElem->appendChild(connPtsElem);

      std::string connPtsAtoms = "";
      for (unsigned int x = 0; x < connPtsList.size(); x++) {
        connPtsAtoms+= (int2String(connPtsList[x]) + " ");
      }
      connPtsElem->setAttribute(X("atoms"), X(connPtsAtoms.c_str()));
    }
}
#endif // USE_XERCES

} // MTKpp namespace

