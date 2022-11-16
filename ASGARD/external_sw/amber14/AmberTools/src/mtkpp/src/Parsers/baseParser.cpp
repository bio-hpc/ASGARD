/*!
   \file baseParser.cpp
   \brief Base parser class
   \author Martin Peters

   Base Class for all parsers

   $Date: 2010/04/29 19:06:19 $
   $Revision: 1.10 $

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

#include "baseParser.h"

#include <string>
#include <sstream>

#include "StringManip.h"

#include "parsingException.h"

namespace MTKpp
{

// ============================================================
// Function : baseParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
baseParser::baseParser() {
    this->bError = false;
    this->errorMessage = "";
    this->NUM_INDENTS_PER_SPACE=2;
}

// ============================================================
// Function : ~baseParser()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
baseParser::~baseParser() {
}

// ============================================================
// Function : Read()
// ------------------------------------------------------------
//
// ============================================================
void baseParser::Read() {
}

// ============================================================
// Function : Write()
// ------------------------------------------------------------
//
// ============================================================
void baseParser::Write() {
}

// ============================================================
// Function : OpenFile()
// ------------------------------------------------------------
//
// ============================================================
std::ofstream& baseParser::OpenFile(std::string fileName)
{
    this->outputFileStream.open(fileName.c_str());

    if (!this->outputFileStream) {
      std::cout << "\n UNABLE TO OPEN FILE"
                << "\nFILENAME = " << fileName << std::endl;
    }
    return this->outputFileStream;
}

// =========================================================
// Function : determineElement
// ---------------------------------------------------------
// Preceive element symbol from the atom name 
// =========================================================
std::string baseParser::determineElement(std::string &name) {
    // HARD CODING -- TO GET ELEMENT SYMBOL FROM ATOM NAME
    //   H                                                  He
    //   Li Be                               B  S  N  O  F  Ne
    //   Na Mg                               Al Si P  S  Cl Ar
    //   K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
    //   Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
    //   Cs Ba Ln Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
    //   Fr Ra An Rf Db Sg Bh Hs Mt Ds Rg Cn

    //   Ln: La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu
    //   An: Ac Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No Lr

    // - Check for the first alphabetic character and then the second - //
    std::string element  = GetAlphaChar(name,1);

    // AC, AG, AL, AM, AR, AS, AT AND AU //
    if (element == "A" ) {
      if (GetAlphaChar(name,2) != "0") {
        std::string second_letter = GetAlphaChar(name,2);
        if (second_letter == "C" || second_letter == "c" ) {
          element  = "Ac";
        }
        if (second_letter == "G" || second_letter == "g" ) {
          element  = "Ag";
        }
        if (second_letter == "L" || second_letter == "l" ) {
          element  = "Al";
        }
        if (second_letter == "M" || second_letter == "m" ) {
          element  = "Am";
        }
        if (second_letter == "R" || second_letter == "r" ) {
          element  = "Ar";
        }
        if (second_letter == "S" || second_letter == "s" ) {
          element  = "As";
        }
        if (second_letter == "T" || second_letter == "t" ) {
          element  = "At";
        }
        if (second_letter == "U" || second_letter == "u" ) {
          element  = "Au";
        }
      }
    }

    // B, BA, BE, BH, BI, BK AND BR //
    if (element == "B" ) {
      if (GetAlphaChar(name,2) != "0") {
        std::string second_letter = GetAlphaChar(name,2);
        if (second_letter == "A" || second_letter == "a" ) {
          element  = "Ba";
        }
        if (second_letter == "E" || second_letter == "e" ) {
          element  = "Be";
        }
        if (second_letter == "H" || second_letter == "h" ) {
          element  = "Bh";
        }
        if (second_letter == "I" || second_letter == "i" ) {
          element  = "Bi";
        }
        if (second_letter == "K" || second_letter == "k" ) {
          element  = "Bk";
        }
        if (second_letter == "R" || second_letter == "r" ) {
          element  = "Br";
        }
      }
      else {
        element  = "B";
      }
    }

    // C, CA, CD, CE, CF, CL, CM, CN, CO, CR, CS AND CU //
    if (element == "C") {
      if (isalnum(name[0])) {
        if (GetAlphaChar(name,2) != "0") {
          std::string second_letter = GetAlphaChar(name,2);

          if ((second_letter == "A" || second_letter == "a") &&
              (removeCharacter((name),' ') == "CA")) { 
            element  = "Ca";
          }
          if (second_letter == "D" || second_letter == "d" ) {
            element  = "Cd";
          }
          if (second_letter == "E" || second_letter == "e" ) {
            element  = "Ce";
          }
          if (second_letter == "F" || second_letter == "f" ) {
            element  = "Cf";
          }
          if (second_letter == "L" || second_letter == "l" ) {
            element  = "Cl";
          }
          if (second_letter == "M" || second_letter == "m" ) {
            element  = "Cm";
          }
          if (second_letter == "N" || second_letter == "n" ) {
            element  = "Cn";
          }
          if (second_letter == "O" || second_letter == "o" ) {
            element  = "Co";
          }
          if (second_letter == "R" || second_letter == "r" ) {
            element  = "Cr";
          }
          if (second_letter == "S" || second_letter == "s" ) {
            element  = "Cs";
          }
          if (second_letter == "U" || second_letter == "u" ) {
            element  = "Cu";
          }
        }
      }
      else {
        std::string second_letter = GetAlphaChar(name,2);
        if (second_letter == "L" || second_letter == "l" ) {
          element  = "Cl";
        }
        else {
          element  = "C";
        }
      }
    }

    // DB, DS AND DY //
    if (element == "D" ) {
      if (GetAlphaChar(name,2) != "0") {
        std::string second_letter = GetAlphaChar(name,2);
        if (second_letter == "B" || second_letter == "b" ) {
          element  = "Db";
        }
        if (second_letter == "S" || second_letter == "s" ) {
          element  = "Ds";
        }
        if (second_letter == "Y" || second_letter == "y" ) {
          element  = "Dy";
        }
      }
    }

    // ER, ES AND EU //
    if (element == "E" ) {
      if (GetAlphaChar(name,2) != "0") {
        std::string second_letter = GetAlphaChar(name,2);
        if (second_letter == "R" || second_letter == "r" ) {
          element  = "Er";
        }
        if (second_letter == "S" || second_letter == "s" ) {
          element  = "Es";
        }
        if (second_letter == "U" || second_letter == "u" ) {
          element  = "Eu";
        }
      }
    }

    // F, FE, FM AND FR //
    if (element == "F" ) {
      if (GetAlphaChar(name,2) != "0") {
        std::string second_letter = GetAlphaChar(name,2);
        if (second_letter == "E" || second_letter == "e" ) {
          element  = "Fe";
        }
        if (second_letter == "M" || second_letter == "m" ) {
          element  = "Fm";
        }
        if (second_letter == "R" || second_letter == "r" ) {
          element  = "Fr";
        }
      }
      else {
        element  = "F";
      }
    }

    // GA, GD AND GE //
    if (element == "G" ) {
      if (GetAlphaChar(name,2) != "0") {
        std::string second_letter = GetAlphaChar(name,2);
        if (second_letter == "A" || second_letter == "a" ) {
          element  = "Ga";
        }
        if (second_letter == "D" || second_letter == "d" ) {
          element  = "Gd";
        }
        if (second_letter == "E" || second_letter == "e" ) {
          element  = "Ge";
        }
      }
    }

    // H, HE, HF, HG, HO AND HS //
    if (element == "H") {
      if (isalnum(name[0])) {
        if (GetAlphaChar(name,2) != "0") {
          std::string second_letter = GetAlphaChar(name,2);
          if (second_letter == "E" || second_letter == "e" ) {
            element  = "He";
          }
          if (second_letter == "F" || second_letter == "f" ) {
            element  = "Hf";
          }
          if (second_letter == "G" || second_letter == "g" ) {
            element  = "Hg";
          }
          if (second_letter == "O" || second_letter == "o" ) {
            element  = "Ho";
          }
          if (second_letter == "S" || second_letter == "s" ) {
            element  = "Hs";
          }
        }
      }
      else {
        element  = "H"; 
      }
    }

    // I, IN AND IR //
    if (element == "I" ) {
      if (GetAlphaChar(name,2) != "0") {
        std::string second_letter = GetAlphaChar(name,2);
        if (second_letter == "N" || second_letter == "n" ) {
          element  = "In";
        }
        if (second_letter == "R" || second_letter == "r" ) {
          element  = "Ir";
        }
      }
      else {
        element  = "I";
      }
    }

    // K AND KR //
    if (element == "K" ) {
      if (GetAlphaChar(name,2) != "0") {
        std::string second_letter = GetAlphaChar(name,2);
        if (second_letter == "R" || second_letter == "r" ) {
          element  = "Kr";
        }
      }
      else {
        element  = "K";
      }
    }

    // LA, LI, LR AND LU //
    if (element == "L") {
      if (GetAlphaChar(name,2) != "0") {
        std::string second_letter = GetAlphaChar(name,2);

        if (second_letter == "A" || second_letter == "a" ) {
          element  = "La";
        }
        if (second_letter == "I" || second_letter == "i" ) {
          element  = "Li";
        }
        if (second_letter == "R" || second_letter == "r" ) {
          element  = "Lr";
        }
        if (second_letter == "U" || second_letter == "u" ) {
          element  = "Lu";
        }
      }
      else {
        element  = "L";
        std::stringstream ss;
        ss << "Unknown element symbol : " << element
                  << " found in baseParser ... " << std::endl;
         std::cout << ss.str();
        throw parsingException(ss.str());
      }
    }

    // MD, MG, MN, MO AND MT //
    if (element == "M") {
      if (GetAlphaChar(name,2) != "0") {
        std::string second_letter = GetAlphaChar(name,2);
        if (second_letter == "D" || second_letter == "d" ) {
          element  = "Md";
        }
        if (second_letter == "G" || second_letter == "g" ) {
          element  = "Mg";
        }
        if (second_letter == "N" || second_letter == "n" ) {
          element  = "Mn";
        }
        if (second_letter == "O" || second_letter == "o" ) {
          element  = "Mo";
        }
        if (second_letter == "T" || second_letter == "t" ) {
          element  = "Mt";
        }
      }
      else {
        element  = "M";
        std::stringstream ss;
        ss << "Unknown element symbol : " << element
                  << " found in baseParser ... " << std::endl;
         std::cout << ss.str();
        throw parsingException(ss.str());
      }
    }

    // N, NA, NB, ND, NE, NI, NO AND NP //
    if (element == "N") {
      if (isalnum(name[0])) {
        if (GetAlphaChar(name,2) != "0") {
          std::string second_letter = GetAlphaChar(name,2);

          if (second_letter == "A" || second_letter == "a" ) {
            element  = "Na";
          }
          if (second_letter == "B" || second_letter == "b" ) {
            element  = "Nb";
          }
          if (second_letter == "D" || second_letter == "d" ) {
            element  = "Nd";
          }
          if (second_letter == "E" || second_letter == "e" ) {
            element  = "Ne";
          }
          if (second_letter == "I" || second_letter == "i" ) {
            element  = "Ni";
          }
          if (second_letter == "O" || second_letter == "o" ) {
            element  = "No";
          }
          if (second_letter == "P" || second_letter == "p" ) {
            element  = "Np";
          }
        }
      }
    }

    // O AND OS //
    if (element == "O") {
      if (isalnum(name[0])) {
        if (GetAlphaChar(name,2) != "0") {
          std::string second_letter = GetAlphaChar(name,2);
          if (second_letter == "S" || second_letter == "s" ) {
            element  = "Os";
          }
        }
      }
      else {
        element  = "O";
      }
    }

    // P, PA, PB, PD, PM, PO, PR, PT AND PU //
    if (element == "P") {
      if (isalnum(name[0])) {
        if (GetAlphaChar(name,2) != "0") {
          std::string second_letter = GetAlphaChar(name,2);
          if (second_letter == "A" || second_letter == "a" ) {
            element  = "Pa";
          }
          if (second_letter == "B" || second_letter == "b" ) {
            element  = "Pb";
          }
          if (second_letter == "D" || second_letter == "d" ) {
            element  = "Pd";
          }
          if (second_letter == "M" || second_letter == "m" ) {
            element  = "Pm";
          }
          if (second_letter == "O" || second_letter == "o" ) {
            element  = "Po";
          }
          if (second_letter == "R" || second_letter == "r" ) {
            element  = "Pr";
          }
          if (second_letter == "T" || second_letter == "t" ) {
            element  = "Pt";
          }
          if (second_letter == "U" || second_letter == "u" ) {
            element  = "Pu";
          }
        }
      }
      else {
        element  = "P";
      }
    }

    // RA, RB, RE, RF, RG, RH, RN AND RU //
    if (element == "R" ) {
      if (GetAlphaChar(name,2) != "0" ) {
        std::string second_letter = GetAlphaChar(name,2);
        if (second_letter == "A" || second_letter == "a" ) {
          element  = "Ra";
        }
        if (second_letter == "B" || second_letter == "b" ) {
          element  = "Rb";
        }
        if (second_letter == "E" || second_letter == "e" ) {
          element  = "Re";
        }
        if (second_letter == "F" || second_letter == "f" ) {
          element  = "Rf";
        }
        if (second_letter == "G" || second_letter == "g" ) {
          element  = "Rg";
        }
        if (second_letter == "H" || second_letter == "h" ) {
          element  = "Rh";
        }
        if (second_letter == "N" || second_letter == "n" ) {
          element  = "Rn";
        }
        if (second_letter == "U" || second_letter == "u" ) {
          element  = "Ru";
        }
      }
    }

    // S, SB, SC, SE, SG, SI, SM, SN AND SR //
    if (element == "S") {
      if (isalnum(name[0])) {
        if (GetAlphaChar(name,2) != "0") {
          std::string second_letter = GetAlphaChar(name,2);
          if (second_letter == "B" || second_letter == "b" ) {
            element  = "Sb";
          }
          if (second_letter == "C" || second_letter == "c" ) {
            element  = "Sc";
          }
          if (second_letter == "E" || second_letter == "e" ) {
            element  = "Se";
          }
          if (second_letter == "G" || second_letter == "g" ) {
            element  = "Sg";
          }
          if (second_letter == "I" || second_letter == "i" ) {
            element  = "Si";
          }
          if (second_letter == "M" || second_letter == "m" ) {
            element  = "Sm";
          }
          if (second_letter == "N" || second_letter == "n" ) {
            element  = "Sn";
          }
          if (second_letter == "R" || second_letter == "r" ) {
            element  = "Sr";
          }
        }
      }
      else {
        element  = "S";
      }
    }

    // TA, TB, TC, TE, TH, TI, TL AND TM //
    if (element == "T" ) {
      if (GetAlphaChar(name,2) != "0" ) {
        std::string second_letter = GetAlphaChar(name,2);
        if (second_letter == "A" || second_letter == "a" ) {
          element  = "Ta";
        }
        if (second_letter == "B" || second_letter == "b" ) {
          element  = "Tb";
        }
        if (second_letter == "C" || second_letter == "c" ) {
          element  = "Tc";
        }
        if (second_letter == "E" || second_letter == "e" ) {
          element  = "Te";
        }
        if (second_letter == "H" || second_letter == "h" ) {
          element  = "Th";
        }
        if (second_letter == "I" || second_letter == "i" ) {
          element  = "Ti";
        }
        if (second_letter == "L" || second_letter == "l" ) {
          element  = "Tl";
        }
        if (second_letter == "M" || second_letter == "m" ) {
          element  = "Tm";
        }
      }
    }

    // Y AND YB //
    if (element == "Y") {
      if (isalnum(name[0])) {
        if (GetAlphaChar(name,2) != "0") {
          std::string second_letter = GetAlphaChar(name,2);
          if (second_letter == "B" || second_letter == "b" ) {
            element  = "Yb";
          }
        }
      }
      else {
        element  = "Y";
      }
    }

    // ZN AND ZR //
    if (element == "Z" ) {
      if (GetAlphaChar(name,2) != "0" ) {
        std::string second_letter = GetAlphaChar(name,2);
        if (second_letter == "N" || second_letter == "n" ) {
          element  = "Zn";
        }
        if (second_letter == "R" || second_letter == "r" ) {
          element  = "Zr";
        }
      }
    }
    //std::cout << "   baseParser::determineElement " << name << " --> " << element << std::endl;

    return element;
}

#ifdef USE_QT

// ==========================================
// Function : string2QString
// ------------------------------------------
// Returns a QString of the std::string
// ==========================================
QString baseParser::string2QString(std::string s)
{
    return QString::fromStdString(s);
}

// ==========================================
// Function : int2QString
// ------------------------------------------
// Returns a QString of the int
// ==========================================
QString baseParser::int2QString(int i)
{
   return QString::number(i, 10);
}

// ==========================================
// Function : double2QString
// ------------------------------------------
// Returns a QString of the double
// ==========================================
QString baseParser::double2QString(double d)
{
   return QString::number(d);
}
#endif // USE_QT

#ifdef USE_TINYXML

// =========================================================
// Function : getIndent
// ---------------------------------------------------------
//
// =========================================================
const char* baseParser::getIndent( unsigned int numIndents )
{
     static const char * pINDENT="                                      + ";
     static const unsigned int LENGTH=strlen( pINDENT );
     unsigned int n=numIndents*NUM_INDENTS_PER_SPACE;
     if ( n > LENGTH ) n = LENGTH;

     return &pINDENT[ LENGTH-n ];
}

// =========================================================
// Function : getIndentAlt
// ---------------------------------------------------------
// same as getIndent but no "+" at the end
// =========================================================
const char* baseParser::getIndentAlt( unsigned int numIndents )
{
     static const char * pINDENT="                                        ";
     static const unsigned int LENGTH=strlen( pINDENT );
     unsigned int n=numIndents*NUM_INDENTS_PER_SPACE;
     if ( n > LENGTH ) n = LENGTH;

     return &pINDENT[ LENGTH-n ];
}

// =========================================================
// Function : dump_attribs_to_stdout
// ---------------------------------------------------------
//
// =========================================================
int baseParser::dump_attribs_to_stdout(TiXmlElement* pElement, unsigned int indent)
{
     if ( !pElement ) return 0;

     TiXmlAttribute* pAttrib=pElement->FirstAttribute();
     int i=0;
     int ival;
     double dval;
     const char* pIndent=getIndent(indent);
     printf("\n");
     while (pAttrib)
     {
          printf( "%s%s: value=[%s]", pIndent, pAttrib->Name(), pAttrib->Value());

          if (pAttrib->QueryIntValue(&ival)==TIXML_SUCCESS)    printf( " int=%d", ival);
          if (pAttrib->QueryDoubleValue(&dval)==TIXML_SUCCESS) printf( " d=%1.1f", dval);
          printf( "\n" );
          i++;
          pAttrib=pAttrib->Next();
     }
     return i;
}

// =========================================================
// Function : dump_to_stdout
// ---------------------------------------------------------
//
// =========================================================
void baseParser::dump_to_stdout( TiXmlNode* pParent, unsigned int indent)
{
     if ( !pParent ) return;

     TiXmlNode* pChild;
     TiXmlText* pText;
     int t = pParent->Type();
     printf( "%s", getIndent(indent));
     int num;

     switch ( t )
     {
     case TiXmlNode::DOCUMENT:
          printf( "Document" );
          break;

     case TiXmlNode::ELEMENT:
          printf( "Element [%s]", pParent->Value() );
          num=dump_attribs_to_stdout(pParent->ToElement(), indent+1);
          switch(num)
          {
               case 0:  printf( " (No attributes)"); break;
               case 1:  printf( "%s1 attribute", getIndentAlt(indent)); break;
               default: printf( "%s%d attributes", getIndentAlt(indent), num); break;
          }
          break;

     case TiXmlNode::COMMENT:
          printf( "Comment: [%s]", pParent->Value());
          break;

     case TiXmlNode::UNKNOWN:
          printf( "Unknown" );
          break;

     case TiXmlNode::TEXT:
          pText = pParent->ToText();
          printf( "Text: [%s]", pText->Value() );
          break;

     case TiXmlNode::DECLARATION:
          printf( "Declaration" );
          break;
     default:
          break;
     }
     printf( "\n" );
     for ( pChild = pParent->FirstChild(); pChild != 0; pChild = pChild->NextSibling()) 
     {
          dump_to_stdout( pChild, indent+1 );
     }
}

// =========================================================
// Function : dump_to_stdout
// ---------------------------------------------------------
// load the named file and dump its structure to STDOUT
// =========================================================
void baseParser::dump_to_stdout(const char* pFilename)
{
     TiXmlDocument doc(pFilename);
     bool loadOkay = doc.LoadFile();
     if (loadOkay)
     {
          printf("\n%s:\n", pFilename);
          dump_to_stdout( &doc ); // defined later in the tutorial
     }
     else
     {
          printf("Failed to load file \"%s\"\n", pFilename);
     }
}
#endif // USE_TINYXML

} // MTKpp namespace

