/*!
   \file printHeader.h

   \brief Print header

   \author Martin B. Peters

   $Date: 2010/03/29 20:33:22 $
   $Revision: 1.4 $

   ----------------------------------------------------------------------------

   MTK++ - C++ package of modeling libraries.

   Copyright (C) 2005-2006  (see AUTHORS file for a list of contributors)

   ----------------------------------------------------------------------------
*/

#ifndef PRINT_HEADER_H
#define PRINT_HEADER_H

#include "copyright.h"
#include <vector>
namespace MTKpp
{

/*!
   \brief Print MTK++ and Program Headers
*/
void printHeader(std::ostream& os, std::string prog_name, std::vector<std::string> authors) {
    copyright(os);
    os << "|                                                                              |" << std::endl;
    std::string outline = "                                                                              ";

    if (prog_name.size() > outline.size()+2) {
      os << " Error in printHeader: Program name is too long ... exiting " << std::endl;
      exit(0);
    }

    for (unsigned int i = 0; i < authors.size(); i++) {
      if (authors[i].size() > outline.size()+2) {
        os << " Error in printHeader: Author name is too long ... exiting " << std::endl;
        exit(0);
      }
    }

    unsigned int s = outline.size() - prog_name.size();
    int odd = 0;
    if (s % 2) {
      odd = 1;
      s--;
    }
    unsigned int w = s/2;
    std::string sub1 = outline.substr(0,w);
    std::string sub2 = outline.substr(0,w+odd);
    std::string sub = "|" + sub1 + prog_name + sub2 + "|";
    if (sub.size() > 80) {
      os << " Error in printHeader ... exiting " << std::endl;
      exit(0);
    }
    os << sub << std::endl;

    for (unsigned int i = 0; i < authors.size(); i++) {
      s = 0;
      w = 0;
      odd = 0;
      sub1 = "";
      sub2 = "";
      sub = "";
      s = outline.size() - authors[i].size();
      if (s % 2) {
        odd = 1;
        s--;
      }
      w = s/2;
      sub1 = outline.substr(0,w);
      sub2 = outline.substr(0,w+odd);
      sub = "|" + sub1 + authors[i] + sub2 + "|";
      if (sub.size() > 80) {
        os << " Error in printHeader ... exiting " << std::endl;
        exit(0);
      }
      os << sub << std::endl;
    }
    os << "|                            Copyright (c) 2005-2011                           |" << std::endl;
    os << "|                              All Rights Reserved.                            |" << std::endl;
    os << "|                                                                              |" << std::endl;
    os << "+------------------------------------------------------------------------------+" << std::endl;
};

} // MTKpp namespace

#endif // PRINT_HEADER_H
