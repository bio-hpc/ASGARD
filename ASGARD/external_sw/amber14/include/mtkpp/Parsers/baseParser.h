/*!
   \file baseParser.h
   \brief Base parser class
   \author Martin Peters

   Base Class for all parsers

   $Date: 2010/08/11 21:11:00 $
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

#ifndef BASEPARSER_H
#define BASEPARSER_H

#ifdef __INTEL_COMPILER

// remark #177: variable was declared but never referenced
#pragma warning(disable:177)

// remark #181: argument is incompatible with corresponding format string conversion
#pragma warning(disable:181)

// remark #304: access control not specified ("public" by default)
#pragma warning(disable:304)

// remark #383: value copied to temporary, reference to temporary used
#pragma warning(disable:383)

// remark #424: extra ";" ignored
#pragma warning(disable:424)

// remark #593: variable was set but never used
#pragma warning(disable:593)

// remark #810: conversion from "double" to "int" may lose significant bits
#pragma warning(disable:810)

// remark #869: parameter was never referenced
#pragma warning(disable:869)

// remark #981: operands are evaluated in unspecified order
#pragma warning(disable:981)

// warning #1125: virtual function override intended?
#pragma warning(disable:1125)

// remark #1418: external function definition with no prior declaration
#pragma warning(disable:1418)

// remark #1572: floating-point equality and inequality comparisons are unreliable
// disabled -> everyone knows it, the parser passes this problem
//             deliberately to the user
#pragma warning(disable:1572)

// remark #1599: declaration hides variable "t"
#pragma warning(disable:1599)

// remark #2259: non-pointer conversion from "double" to "int" may lose significant bits
#pragma warning(disable:2259)

#endif

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <stdlib.h>
#include <stdio.h>

#ifdef USE_QT
#include <QtCore>
#include <QtXml>
#elif defined(USE_TINYXML)
#include <sstream>
#include "tinyxml/tinyxml.h"
#endif

namespace MTKpp
{

// ============================================================
// Struct : AtomLine
// ------------------------------------------------------------
/*! 
   \struct AtomLine
   \brief Temporary container for atom info.
   \author Martin Peters
   \date 2005
*/
// ============================================================
struct AtomLine
{
     std::string typ;
     int serial;
     std::string name;
     std::string altLoc;
     std::string resName;
     int resSeq, massdiff,stereo, hcount,strcare,valence;
     std::string iCode;
     double x;
     double y;
     double z;
     std::string segID; 
     std::string element;
     std::string charge;
};

// ============================================================
// Class : baseParser()
// ------------------------------------------------------------
/*! 
   \class baseParser
   \brief base class to all parsers
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================
class baseParser
{
public:

    /*!
       \brief baseParser Constructor
    */
    baseParser();

    //! baseParser Destructor
    virtual ~baseParser();

    /*!
       \brief Read function
    */
    virtual void             Read();

    /*!
       \brief Write function
    */
    virtual void             Write();

    /*!
       \brief Write function
    */
    std::ofstream&            OpenFile(std::string fileName);

    /*!
       \brief Preceive element symbol from the atom name
    */
    std::string               determineElement(std::string &name);

    /*!
       \brief Error handling within Parsers
       \param i error
    */
    void setError(int i) {
      if (i) {
        this->bError = true;
      }
      else {
        this->bError = false;
      }
    }

    /*!
       \brief Error handling within dcParser
       \return error
    */
    bool getError() {
      return this->bError;
    }

    /*!
       \brief Error handling within Parsers
       \param s error message
    */
    void setErrorMessage(std::string s) {
      this->errorMessage = s;
    }

    /*!
       \brief Error handling within Parsers
       \return error message
    */
    std::string getErrorMessage() {
      return this->errorMessage;
    }

#ifdef USE_QT
protected:
    /*!
       \brief Convert an int to a QString
       \param s std::string to be converted
       \return string
    */
    QString string2QString(std::string s);

    /*!
       \brief Convert an int to a QString
       \param i int to be converted
       \return string
    */
    QString int2QString(int i);

    /*!
       \brief Convert a double to a QString
       \param d double to be converted
       \return string
    */
    QString double2QString(double d);

#endif // USE_QT

#ifdef USE_TINYXML
protected:

    /*!
       \brief 
       \return 
    */
    const char* getIndent(unsigned int numIndents);

    /*!
       \brief 
       \return 
    */
    const char* getIndentAlt(unsigned int numIndents);

    /*!
       \brief 
       \return 
    */
    int dump_attribs_to_stdout(TiXmlElement* pElement, unsigned int indent);

    /*!
       \brief 
       \return 
    */
    void dump_to_stdout(TiXmlNode* pParent, unsigned int indent = 0);

    /*!
       \brief 
       \return 
    */
    void dump_to_stdout(const char* pFilename);

#endif // USE_TINYXML

protected:

    //! Output File Stream
    std::ofstream outputFileStream;

    //! Error occured
    bool bError;

    //! Error message
    std::string errorMessage;

    //!
    unsigned int NUM_INDENTS_PER_SPACE;
};

} // MTKpp namespace

#endif // BASEPARSER_H
