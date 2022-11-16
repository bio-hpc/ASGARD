/*!
   \file pamParser.cpp
   \brief Parses pam files
   \author Martin Peters

   Reads pam files

   $Date: 2010/03/29 20:39:35 $
   $Revision: 1.3 $

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

#include "pamParser.h"
#include <vector>
#include <algorithm>

#include "Molecule/seqAlign.h"
#include "Log/errorHandler.h"

namespace MTKpp
{

// ============================================================
// Function : pamParser()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
pamParser::pamParser()
{

}

// =========================================================
// Function : ~pamParser()
// ---------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// =========================================================
pamParser::~pamParser()
{

}

// =========================================================
// Function : Read
// ---------------------------------------------------------
// parsers a pam file
/*

#
# This matrix was produced by "pam" Version 1.0.8 [12-Oct-03]
# using Dayhoff et al. (1978) mutability data.
#
# PAM 250 substitution matrix, scale = ln(2)/3 = 0.231049
#
# Expected score = -0.844, Entropy = 0.354 bits
#
# Lowest score = -8, Highest score = 17
#
  A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  2 -2  0  0 -2  0  0  1 -1 -1 -2 -1 -1 -3  1  1  1 -6 -3  0  0  0  0 -8
R -2  6  0 -1 -4  1 -1 -3  2 -2 -3  3  0 -4  0  0 -1  2 -4 -2 -1  0 -1 -8
N  0  0  2  2 -4  1  1  0  2 -2 -3  1 -2 -3  0  1  0 -4 -2 -2  2  1  0 -8
D  0 -1  2  4 -5  2  3  1  1 -2 -4  0 -3 -6 -1  0  0 -7 -4 -2  3  3 -1 -8
C -2 -4 -4 -5 12 -5 -5 -3 -3 -2 -6 -5 -5 -4 -3  0 -2 -8  0 -2 -4 -5 -3 -8
Q  0  1  1  2 -5  4  2 -1  3 -2 -2  1 -1 -5  0 -1 -1 -5 -4 -2  1  3 -1 -8
E  0 -1  1  3 -5  2  4  0  1 -2 -3  0 -2 -5 -1  0  0 -7 -4 -2  3  3 -1 -8
G  1 -3  0  1 -3 -1  0  5 -2 -3 -4 -2 -3 -5  0  1  0 -7 -5 -1  0  0 -1 -8
H -1  2  2  1 -3  3  1 -2  6 -2 -2  0 -2 -2  0 -1 -1 -3  0 -2  1  2 -1 -8
I -1 -2 -2 -2 -2 -2 -2 -3 -2  5  2 -2  2  1 -2 -1  0 -5 -1  4 -2 -2 -1 -8
L -2 -3 -3 -4 -6 -2 -3 -4 -2  2  6 -3  4  2 -3 -3 -2 -2 -1  2 -3 -3 -1 -8
K -1  3  1  0 -5  1  0 -2  0 -2 -3  5  0 -5 -1  0  0 -3 -4 -2  1  0 -1 -8
M -1  0 -2 -3 -5 -1 -2 -3 -2  2  4  0  6  0 -2 -2 -1 -4 -2  2 -2 -2 -1 -8
F -3 -4 -3 -6 -4 -5 -5 -5 -2  1  2 -5  0  9 -5 -3 -3  0  7 -1 -4 -5 -2 -8
P  1  0  0 -1 -3  0 -1  0  0 -2 -3 -1 -2 -5  6  1  0 -6 -5 -1 -1  0 -1 -8
S  1  0  1  0  0 -1  0  1 -1 -1 -3  0 -2 -3  1  2  1 -2 -3 -1  0  0  0 -8
T  1 -1  0  0 -2 -1  0  0 -1  0 -2  0 -1 -3  0  1  3 -5 -3  0  0 -1  0 -8
W -6  2 -4 -7 -8 -5 -7 -7 -3 -5 -2 -3 -4  0 -6 -2 -5 17  0 -6 -5 -6 -4 -8
Y -3 -4 -2 -4  0 -4 -4 -5  0 -1 -1 -4 -2  7 -5 -3 -3  0 10 -2 -3 -4 -2 -8
V  0 -2 -2 -2 -2 -2 -2 -1 -2  4  2 -2  2 -1 -1 -1  0 -6 -2  4 -2 -2 -1 -8
B  0 -1  2  3 -4  1  3  0  1 -2 -3  1 -2 -4 -1  0  0 -5 -3 -2  3  2 -1 -8
Z  0  0  1  3 -5  3  3  0  2 -2 -3  0 -2 -5  0  0 -1 -6 -4 -2  2  3 -1 -8
X  0 -1  0 -1 -3 -1 -1 -1 -1 -1 -1 -1 -1 -2 -1  0  0 -4 -2 -1 -1 -1 -1 -8
* -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8 -8  1
*/
// =========================================================
void pamParser::Read(const std::string &pamfile, seqAlign* pSeqAlign)
{
    setError(0);
    std::ifstream ipam;
    ipam.open(pamfile.c_str());

    int end   = pamfile.length();
    int slash = pamfile.find_last_of("/");
    std::string file_name = pamfile.substr(slash+1,(end-slash-5));

    if (!ipam) {
      std::string errorMessage = " Error, Can't Find " + pamfile;
      setErrorMessage(errorMessage);
      errorLogger.throwError( "seqAlign::Read", errorMessage, 1);
    }

    if (pSeqAlign) {
      pSeqAlign->setPAMFile(pamfile);
    }
    else {
      std::string errorMessage = " Error ";
      setErrorMessage(errorMessage);
      errorLogger.throwError( "seqAlign::Read", errorMessage, 1);
    }

    std::vector<std::string> characters;
    int i = -1;
    std::string fileline = "";
    while (ipam) {
      getline(ipam,fileline);
      if (fileline.substr(0,1) == "#") {
        continue;
      }
      if (i == -1) { // first line
        splitString(fileline, " ", characters, 0);
        pSeqAlign->setPAMSize(static_cast<int>(characters.size()));
        for (unsigned int i = 0; i < characters.size(); i++) {
          pSeqAlign->setPAMType(i, characters[i]);
        }
      }
      else {
        std::vector<std::string> words;
        splitString(fileline, " ", words, 0);
        for (unsigned int j = 1; j < words.size(); j++) {
          double v = string2Double(words[j]);
          pSeqAlign->setPAMValue(i, j-1, v);
        }
      }
      i++;
    }
    ipam.close();
}

} // MTKpp namespace
