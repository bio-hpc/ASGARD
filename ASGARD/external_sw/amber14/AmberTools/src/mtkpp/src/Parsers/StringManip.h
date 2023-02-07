/*!
   \file StringManip.h
   \brief Miscellaneous string manipulation functions
   \author Andrew Wollacott
   \author Martin Peters

   $Date: 2010/03/29 20:39:34 $
   $Revision: 1.17 $

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

#ifndef STRINGMANIP_H
#define STRINGMANIP_H

#include <stdlib.h>
#include <string.h>
#include <vector>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include "Utils/constants.h"

namespace MTKpp
{

/*!
   \brief Splits string
   \param s string to work on
   \param s2 separator
   \param v vector of string that gets returned
   \param i starting point
*/
void splitString(std::string& s, const std::string s2, std::vector<std::string>& v, int i);

/*!
   \brief Returns a string with leading/trailing characters of a set stripped
   \param str string to work on
   \param sepSet separator
   \return string
*/
std::string stripString(std::string const& str, char const* sepSet);

/*!
   \brief Remove character from string
   \param s string to work on
   \param c character to remove
   \return new string
*/
std::string removeCharacter(const std::string& s, const char& c);

/*!
   \brief Replace character in string
   \param s string to work on
   \param c character to be replaced
   \param n new character
   \return new string
*/
std::string replaceCharacter(const std::string& s, const char& c, const char& n);

/*!
   \brief Add character to the string
   \param s string to work on
   \param c character to add
   \param i number of times to add
   \return new string
*/
std::string addCharacter(const std::string& s, char c,int i);

/*!
   \brief A check to see if a file field exist or not
   \param s string to check
   \param i start point
   \param j end point
   \return boolean yes/no
*/
bool FieldExists(std::string& s, int i, int j);

/*!
   \brief Get the first alpha character in the string
   \param s string to check
   \param i start point
   \return an alphabetic letter in a string at position number.
*/
std::string GetAlphaChar(std::string& s, int i);

/*!
   \brief Convert an int to a string
   \param i int to be converted
   \return string
*/
std::string int2String(int i);

/*!
   \brief Convert an int to a string
   \param i int to be converted
   \return string
*/
std::string uLongKind2String(ULONG_KIND i);

/*!
   \brief Convert a string to an unsigned integer
   \param s string to be converted
   \return unsigned int
*/
unsigned int string2UInt(std::string s);

/*!
   \brief Convert a string to an integer
   \param s string to be converted
   \return int
*/
int string2Int(std::string s);

/*!
   \brief Convert a double to a string
   \param d Double to be converted
   \return string
*/
std::string double2String(double d);

/*!
   \brief Convert a double to a string
   \param d Double to be converted
   \param precision Number precision
   \return string
*/
std::string double2String(double d, int precision);

/*!
   \brief Convert a string to a double
   \param s string to be converted
   \return double
*/
double string2Double(std::string s);

/*!
   \brief Converts to upper case
   \param s string to convert
   \return string
*/
std::string toUpper(std::string s);

/*!
   \brief Converts to lower case
   \param s string to convert
   \return string
*/
std::string toLower(std::string s);

/*!
   \brief  Returns a base of the string, e.g "asdf" gets returned from "asdf.mol"
   \param name string to check
   \return string
*/
std::string baseName(std::string name);

/*!
   \brief  Returns a file extension, e.g "mol" gets returned from "asdf.mol"
   \param e string to check
   \return string
*/
std::string extName(std::string e);

/*!
   \brief Returns a boolean if a file exists or not
   \param fileName file name
   \return boolean
*/
bool fileExists(const std::string& fileName);

/*!
   \brief Returns a boolean if a string contains a certain sub-string
   \param s main string
   \param ss sub-string
   \return boolean
*/
bool containsSubStr(const std::string& s, const std::string& ss);

/*!
   \brief Replace part of string with another string
   \param s1 main string
   \param s2 sub-string to be replaced
   \param s3 string that will be entered into the new string
   \return new string
*/
std::string replaceSubStr(const std::string& s1, const std::string& s2, const std::string& s3);

} // MTKpp namespace

#endif // STRINGMANIP_H

