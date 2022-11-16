/*!
   \file inputParser.h
   \brief Parser an input file
   \author Martin Peters

   $Date: 2010/03/29 20:39:35 $
   $Revision: 1.5 $

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

#ifndef INPUTPARSER_H
#define INPUTPARSER_H

#include "baseParser.h"
#include "Log/errorHandler.h"
#include "Utils/constants.h"
namespace MTKpp
{

/*!
   \brief Read generic input file

    \param inputFile
    \param inputFileContents list of options from file

    example:
    \code
    # Read a PDB file
    readPdb file.pdb
    addHs /col/mol # Add Hydrogen atoms
    # Write out MOL file
    writeMol file.mol
    \endcode
    \return success
*/
inline int readInputFile(std::string inputFile, std::vector<std::vector<std::string> >& inputFileContents)
{
    std::ifstream iFile;
    iFile.open(inputFile.c_str());

    if (!iFile) {
      std::string errMess = "UNABLE TO OPEN INPUT FILE FILENAME = " + inputFile + " ... exiting ";
      errorLogger.throwError("readInputFile", errMess, MTK_ERROR);
      return 1;
    }

    std::string errorMessage = "Commands read: \n";

    std::string fileline = "";
    const char comment = '#';

    while (iFile) {
      std::string buffer(80,'*');
      std::string preFileLine;
      getline(iFile, preFileLine);

      // Strip any spaces at the begin and end of the line
      std::string fileline = stripString(preFileLine, " ");

      // If line is not blank or starts with "#"
      if ((fileline != "") and (fileline[0] != comment)) {
        errorMessage += fileline + "\n";
        // User may have a comment at end of line
        std::vector<std::string> splitLine1;
        splitString(fileline, "#", splitLine1, 0);
        fileline = splitLine1[0];

        std::vector<std::string> splitLine2;
        splitString(fileline, " ", splitLine2, 0);
        inputFileContents.push_back(splitLine2);
      }
    } // while (iFile)
    iFile.close();

    errorLogger.throwError("readInputFile", errorMessage, INFO);

    return 0;
};

/*!
   \brief Read generic list file

    \param listFile
    \param listFileContents list file
    \return success
*/
inline int readListFile(std::string listFile, std::vector<std::string>& listFileContents)
{
    std::ifstream iFile;
    iFile.open(listFile.c_str());

    if (!iFile) {
      std::cout << "\nUNABLE TO OPEN LIST FILE" << "\nFILENAME = "
                << listFile << "\nEXITING...\n" << std::endl;
      return 1;
    }
    std::string fileline = "";

    while (iFile) {
      std::string buffer(80,'*');
      getline(iFile, fileline);
      if (fileline != "") {
        listFileContents.push_back(fileline);
      }
    } // while (iFile)
    iFile.close();
    return 0;
};

} // MTKpp namespace

#endif // INPUTPARSER_H

