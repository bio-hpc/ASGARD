/*!
   \file sheet.cpp
   \brief 
   \author Martin Peters

   $Date: 2008/07/22 08:08:29 $
   $Revision: 1.4 $

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

#include "sheet.h"

namespace MTKpp
{
// ============================================================
// Class : sheet()
// ------------------------------------------------------------
// Constructor for the class.
// ============================================================
sheet::sheet()
{
    this->itsName = "";
    this->pTable = 0;
}

// ============================================================
// Function : ~sheet()
// ------------------------------------------------------------
// Destructor for the class
// All data is destroyed.
// ============================================================
sheet::~sheet()
{
    this->clear();
}

// ============================================================
// Function : setName
// ------------------------------------------------------------
// 
// ============================================================
void sheet::setName(std::string n)
{
    itsName = n;
}

// ============================================================
// Function : getName
// ------------------------------------------------------------
// 
// ============================================================
std::string sheet::getName()
{
    return this->itsName;
}

// ============================================================
// Function : addTable
// ------------------------------------------------------------
// 
// ============================================================
table<double>* sheet::addTable()
{
    table<double>* myTable = new table<double>();
    this->tableList.push_back(myTable);
    return myTable;
}

// ============================================================
// Function : addIntTable
// ------------------------------------------------------------
// 
// ============================================================
table<int>* sheet::addIntTable()
{
    table<int>* myTable = new table<int>();
    this->intTableList.push_back(myTable);
    return myTable;
}

// ============================================================
// Function : addStringTable
// ------------------------------------------------------------
// 
// ============================================================
table<std::string>* sheet::addStringTable()
{
    table<std::string>* myTable = new table<std::string>();
    this->stringTableList.push_back(myTable);
    return myTable;
}

// ============================================================
// Function : getTable
// ------------------------------------------------------------
// 
// ============================================================
table<double>* sheet::getTable(std::string name)
{
    for (unsigned int i = 0; i < this->tableList.size(); i++) {
      if (this->tableList[i]->getName() == name) {
        return this->tableList[i];
      }
    }
    return NULL;
}

// ============================================================
// Function : getIntTable
// ------------------------------------------------------------
// 
// ============================================================
table<int>* sheet::getIntTable(std::string name)
{
    for (unsigned int i = 0; i < this->intTableList.size(); i++) {
      if (this->intTableList[i]->getName() == name) {
        return this->intTableList[i];
      }
    }
    return NULL;
}

// ============================================================
// Function : getTables
// ------------------------------------------------------------
// 
// ============================================================
std::vector<table<double>*> sheet::getTables()
{
    return this->tableList;
}

// ============================================================
// Function : getIntTables
// ------------------------------------------------------------
// 
// ============================================================
std::vector<table<int>*> sheet::getIntTables()
{
    return this->intTableList;
}

// ============================================================
// Function : clear
// ------------------------------------------------------------
// 
// ============================================================
void sheet::clear()
{
    for (tableIterator t = this->tableList.begin(); t != this->tableList.end(); t++) {
      pTable = *t;
      delete pTable;
    }
    this->tableList.clear();

    for (intTableIterator t = this->intTableList.begin(); t != this->intTableList.end(); t++) {
      pIntTable = *t;
      delete pIntTable;
    }
    this->intTableList.clear();
}

}
