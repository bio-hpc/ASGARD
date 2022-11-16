/*!
   \file copyright.h

   \brief MTK++ copyright

   \author Martin B. Peters

   $Date: 2010/03/29 20:33:22 $
   $Revision: 1.7 $

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

#ifndef COPYRIGHT_H
#define COPYRIGHT_H
#include "Utils/constants.h"

#include <stdlib.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include "config.h"

namespace MTKpp
{
void copyright(std::ostream& os)
{
    std::string outline = "                                                                              ";
    std::string MTK_NAME = PACKAGE_STRING;

    if (MTK_NAME.size() > outline.size()+2) {
      os << " Error in copyright: Program name is too long ... exiting " << std::endl;
      exit(0);
    }

    unsigned int s = outline.size() - MTK_NAME.size();
    int odd = 0;
    if (s % 2) {
      odd = 1;
      s--;
    }
    unsigned int w = s/2;
    std::string sub1 = outline.substr(0,w);
    std::string sub2 = outline.substr(0,w+odd);
    std::string sub = "|" + sub1 + MTK_NAME + sub2 + "|";
    if (sub.size() > 80) {
      os << " Error in copyright ... exiting " << std::endl;
      exit(0);
    }

    os << "+------------------------------------------------------------------------------+"
       << std::endl;
    os << "|____/\\/\\______/\\/\\__/\\/\\/\\/\\/\\/\\__/\\/\\____/\\/\\________________________________|"
       << std::endl;
    os << "|____/\\/\\/\\__/\\/\\/\\______/\\/\\______/\\/\\__/\\/\\________/\\/\\__________/\\/\\________|"
       << std::endl;
    os << "|____/\\/\\/\\/\\/\\/\\/\\______/\\/\\______/\\/\\/\\/\\______/\\/\\/\\/\\/\\/\\__/\\/\\/\\/\\/\\/\\____|"
       << std::endl;
    os << "|____/\\/\\__/\\__/\\/\\______/\\/\\______/\\/\\__/\\/\\________/\\/\\__________/\\/\\________|"
       << std::endl;
    os << "|____/\\/\\______/\\/\\______/\\/\\______/\\/\\____/\\/\\________________________________|"
       << std::endl;
    os << "+------------------------------------------------------------------------------+" << std::endl;
    os << "|                                                                              |" << std::endl;
    os << "|                                 ____    ____                                 |" << std::endl;
    os << "|                                /    \\  /    \\                                |" << std::endl;
    os << "|                               /      \\/      \\                               |" << std::endl;
    os << "|                               \\              /                               |" << std::endl;
    os << "|                                \\            /                                |" << std::endl;
    os << "|                                /            \\                                |" << std::endl;
    os << "|                               /              \\                               |" << std::endl;
    os << "|                               \\      /\\      /                               |" << std::endl;
    os << "|                                \\____/  \\____/                                |" << std::endl;
    os << "|                                                                              |" << std::endl;
    os << sub << std::endl;
    //os << "| MTK++ - C++ package of modeling libraries.                                   |" << std::endl;
    os << "|                                                                              |" << std::endl;
    os << "| Copyright (C) 2005-2010  (see AUTHORS file for a list of contributors)       |" << std::endl;
    os << "|                                                                              |" << std::endl;
    os << "| MTK++ is free software; you can redistribute it and/or modify                |" << std::endl;
    os << "| it under the terms of the GNU Lesser General Public License as published by  |" << std::endl;
    os << "| the Free Software Foundation; either version 3 of the License, or            |" << std::endl;
    os << "| (at your option) any later version.                                          |" << std::endl;
    os << "|                                                                              |" << std::endl;
    os << "| MTK++ is distributed in the hope that it will be useful,                     |" << std::endl;
    os << "| but WITHOUT ANY WARRANTY; without even the implied warranty of               |" << std::endl;
    os << "| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                |" << std::endl;
    os << "| GNU Lessser General Public License for more details.                         |" << std::endl;
    os << "|                                                                              |" << std::endl;
    os << "| You should have received a copy of the GNU Lesser General Public License     |" << std::endl;
    os << "| along with this program.  If not, see <http://www.gnu.org/licenses/>.        |" << std::endl;
    os << "|                                                                              |" << std::endl;
    os << "+------------------------------------------------------------------------------+" << std::endl;
    os << "|                                                                              |" << std::endl;
    os << "| This program uses the MTK++ package.                                         |" << std::endl;
    os << "|                                                                              |" << std::endl;
    os << "+------------------------------------------------------------------------------+" << std::endl;
};

} // namespace MTKpp

#endif // COPYRIGHT_H
