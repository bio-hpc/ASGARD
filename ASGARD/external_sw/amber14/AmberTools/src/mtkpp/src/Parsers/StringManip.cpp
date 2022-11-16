/*!
   \file StringManip.cpp
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

#include "StringManip.h"
#include <sstream>
#include <cmath>

namespace MTKpp
{

//=============================================================
// Function : splitString
// ------------------------------------------------------------
// splits up a string based on a separator and returns a vector
// ============================================================
void splitString(std::string &text, const std::string separator,
                 std::vector<std::string> &words, int mystart)
{
    int n = text.length();
    int start, stop;

    start = text.find_first_not_of(separator, mystart);

    while ((start >= 0) && (start < n)) {
      stop = text.find_first_of(separator, start);
      if ((stop < 0) || (stop > n)) {
        stop = n;
      }

      words.push_back(text.substr(start, stop-start));
      start = text.find_first_not_of(separator, stop+1);
    }
}

//=============================================================
// Function : stripString
// ------------------------------------------------------------
// Returns a string with leading/trailing characters of a set stripped
// ============================================================
std::string stripString(std::string const& str, char const* sepSet)
{
    std::string::size_type const first = str.find_first_not_of(sepSet);
    return ( first == std::string::npos )
    ? std::string()
    : str.substr(first, str.find_last_not_of(sepSet)-first+1);
}

//=============================================================
// Function : removeCharacter
// ------------------------------------------------------------
// removes a character from a string
// ============================================================
std::string removeCharacter(const std::string &original, const char &remove)
{
   std::string output = "";
   int n = original.length();
   for (int i = 0; i < n; i++) {
     if (original[i] != remove) {
       output += original[i];
     }
   }
   return output;
}

//=============================================================
// Function : replaceCharacter
// ------------------------------------------------------------
// removes a character from a string
// ============================================================
std::string replaceCharacter(const std::string &original, const char &remove, const char& replace)
{
   std::string output = "";
   int n = original.length();
   for (int i = 0; i < n; i++) {
     if (original[i] == remove) {
       output += replace;
     }
     else {
       output += original[i];
     }
   }
   return output;
}

//=============================================================
// Function : addCharacter
// ------------------------------------------------------------
// adds a character to a string
// ============================================================
std::string addCharacter(const std::string &original, char c_add, int numchar)
{
    std::string output = original;
    for (int i = 0; i < numchar; i++) {
      output = output + c_add;
    }
    return output;
}

// ==========================================
// Function : FieldExists
// ------------------------------------------
// Check to see if field in file exists or not
// ==========================================
bool FieldExists(std::string &fileline, int start, int size)
{
    int numBlanks = (fileline.substr(start,size)).find_first_not_of(" ");
    if (numBlanks < 0) {
      return false;
    }
    return true;
}

// ==========================================
// Function : GetAlphaChar
// ------------------------------------------
// Returns the an alphabetic letter in a string at position number.
// ==========================================
std::string GetAlphaChar(std::string &original, int number)
{
    std::string st = original;
    int n=1;
    for (unsigned int i = 0; i < st.length(); i++) {
      if ( isalpha(st[i]) ) {
        if (n == number) {
          return st.substr(i,1);
        }
        else {
          n=n+1;
        }
      }
    }
    return "0";
}

// ==========================================
// Function : int2String
// ------------------------------------------
// Returns a string of the int
// ==========================================
std::string int2String(int i)
{
    std::stringstream number;
    number << i;
    return number.str();
}

// ==========================================
// Function : uLongKind2String
// ------------------------------------------
// Returns a string of the ulong_kind
// ==========================================
std::string uLongKind2String(ULONG_KIND i)
{
    std::stringstream number;
    number << i;
    return number.str();
}

// ==========================================
// Function : string2UInt
// ------------------------------------------
// Returns a string of the int
// ==========================================
unsigned int string2UInt(std::string s)
{
    //return static_cast<unsigned int>(atoi(s.c_str()));
    unsigned int temp;
    std::stringstream ssout(s);
    ssout>>temp;
    return static_cast<unsigned int>(temp);
}

// ==========================================
// Function : string2Int
// ------------------------------------------
// Returns a string of the int
// ==========================================
int string2Int(std::string s)
{
    //return atoi(s.c_str());
    int temp;
    std::stringstream ssout(s);
    ssout>>temp;
    return temp;
}

// ==========================================
// Function : double2String
// ------------------------------------------
// Returns a string of the double
// ==========================================
std::string double2String(double d)
{
    std::stringstream number;
    number << d;
    return number.str();
}

// ==========================================
// Function : double2String
// ------------------------------------------
// Returns a string of the double
// ==========================================
std::string double2String(double d, int p)
{
    std::stringstream number;
    number.setf(std::ios_base::fixed);
    number.setf(std::ios_base::showpoint);
    number.precision(p);
    number << d;
    return number.str();
}

// ==========================================
// Function : string2Double
// ------------------------------------------
// Returns a string of the double
// ==========================================
double string2Double(std::string s)
{
    //return strtod(s.c_str(), 0);
    double temp;
    std::stringstream ssout(s);
    ssout>>temp;
    return temp;
}

// ==========================================
// Function : toUpper
// ------------------------------------------
// Converts string to uppercase
// ==========================================
std::string toUpper(std::string s)
{
    for (unsigned int j=0; j < s.size(); ++j) {
      s[j]=toupper(s[j]);
    }
    return s;
}

// ==========================================
// Function : toLower
// ------------------------------------------
// Converts string to lowercase
// ==========================================
std::string toLower(std::string s)
{
    for (unsigned int j=0; j < s.size(); ++j) {
      s[j]=tolower(s[j]);
    }
    return s;
}

// ==========================================
// Function : baseName
// ------------------------------------------
// Returns a base of the string
// e.g "asdf" gets returned from "asdf.mol"
// ==========================================
std::string baseName(std::string name)
{
    int end   = name.length();
    int slash = name.find_last_of("/");
    std::string outfile_name = name.substr(slash+1, (end-slash-5));
    end   = outfile_name.length();
    int dot = outfile_name.find_last_of(".");
    outfile_name = outfile_name.substr(0,dot);
    return outfile_name;
}

// ==========================================
// Function : extName
// ------------------------------------------
// Returns file extension
// e.g "mol" gets returned from "asdf.mol"
// ==========================================
std::string extName(std::string fileName)
{
    int end   = fileName.length();
    int dot = fileName.find_last_of(".");
    fileName = fileName.substr(dot+1,end);
    return fileName;
}

// ==========================================
// Function : fileExists
// ------------------------------------------
// Check to see if a file exists or not
// ==========================================
bool fileExists(const std::string& fileName)
{
    std::fstream fin;
    fin.open(fileName.c_str(), std::ios::in);

    if ( fin.is_open() ) {
      fin.close();
      return true;
    }
    fin.close();
    return false;
}

// ==========================================
// Function : containsSubStr
// ------------------------------------------
// Check to see if a file exists or not
// ==========================================
bool containsSubStr(const std::string& s, const std::string& ss)
{
    std::string::size_type loc = s.find(ss, 0);
    if (loc != std::string::npos) {
      return true;
    } 
    return false;
}


// ==========================================
// Function : replaceSubStr
// ------------------------------------------
//
// ==========================================
std::string replaceSubStr(const std::string& s1, const std::string& s2, const std::string& s3)
{
    std::string nString = s1;
    bool done = false;

    while (!done) {
      std::string::size_type loc = nString.find(s2, 0);
      if (loc == std::string::npos) {
        done = true;
      }
      else {
        nString.erase( loc, s2.size() );
        nString.insert( loc, s3 );
      }
    }
    return nString;
}

} // MTKpp namespace

