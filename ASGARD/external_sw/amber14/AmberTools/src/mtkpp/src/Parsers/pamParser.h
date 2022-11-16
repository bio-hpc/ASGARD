/*!
   \file pamParser.h
   \brief Parses pam files
   \author Martin Peters

   Reads pam files

   $Date: 2010/02/20 01:33:16 $
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

#ifndef PAMPARSER_H
#define PAMPARSER_H

#include "StringManip.h"
#include "baseParser.h"

namespace MTKpp
{

class seqAlign;

// ============================================================
// Class : pamParser()
// ------------------------------------------------------------
/*!
   \class pamParser
   \brief Reads pam format files
   \author Martin Peters
   \date 2007
*/
// ============================================================
class pamParser : public baseParser
{

public:
    /*!
       \brief pamParser Constructor
    */
    pamParser();

    //! pamParser Destructor
    ~pamParser();

    /*!
       \brief Read PAM file
       \param i pam file
       \param pSeqAlign seqAlign pointer
    */
    void           Read(const std::string &i, seqAlign* pSeqAlign);

};

} // MTKpp namespace

#endif // PAMPARSER_H
