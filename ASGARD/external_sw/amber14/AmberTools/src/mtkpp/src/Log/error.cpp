/*! 
   \file error.cpp
   \brief Error object used in errorHandler
   \author Martin Peters

   $Date: 2007/11/19 11:54:57 $
   $Revision: 1.2 $

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

#include "error.h"

namespace MTKpp
{

// ============================================================
// Function : error()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
error::error(const std::string &function,
             const std::string &message,
             const int type)
{
    this->function = function;
    this->message = message;
    this->type = type;
}

// ============================================================
// Function : ~error()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
error::~error() {}

// ============================================================
// Function : getMessage()
// ------------------------------------------------------------
// 
// ============================================================
std::string error::getFormattedMessage() const
{
    std::string lMessage = "\n ### ### ### ### \n";
    if (this->type == 0) {
      lMessage += " ### MTK++ ### \n";
    }
    if (this->type == 1) {
      lMessage += " ### MTK++ Error ### \n";
    }
    else if (this->type == 2) {
      lMessage += " ### MTK++ Warning ### \n";
    }
    else if (this->type == 3) {
      lMessage += " ### MTK++ Debug ### \n";
    }
    else if (this->type == 4) {
      lMessage += " ### MTK++ Info ### \n";
    }
    lMessage += " ### Function: " + this->function + " ### \n";
    lMessage += " ### Message: ";
    lMessage += this->message;
    lMessage += "\n ### ### ### ### \n";
    return lMessage;
}

} // MTKpp namespace
