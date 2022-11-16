/*!
   \file sheet.h
   \brief 
   \author Martin Peters

   $Date: 2010/03/29 20:35:21 $
   $Revision: 1.6 $

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

#ifndef SHEET_H
#define SHEET_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <math.h>

#include "Utils/constants.h"
#include "table.h"

namespace MTKpp
{

//template <class T> class table;

// ============================================================
// Class : sheet()
// ------------------------------------------------------------
/*!
   \class sheet
   \brief 
   \author Martin Peters
   \version 0.1
   \date 2007
*/
// ============================================================
class sheet
{
public:
    /*!
       \brief sheet Constructor
    */
    sheet();

    //! sheet Destructor
    virtual ~sheet();

    /*!
       \brief Set the name of the sheet
       \param n sheet name
    */
    void setName(std::string n);

    /*!
       \brief Get the name of the sheet
       \return sheet name
    */
    std::string getName();

    /*!
       \brief Add a table to the sheet
       \return table pointer
    */
    table<double>* addTable();

    /*!
       \brief Add an integer table to the sheet
       \return table pointer
    */
    table<int>* addIntTable();

    /*!
       \brief Add a string table to the sheet
       \return table pointer
    */
    table<std::string>* addStringTable();

    /*!
       \brief Get a table by name
       \param n table name
       \return table name
    */
    table<double>* getTable(std::string n);

    /*!
       \brief Get an integer table by name
       \param n table name
       \return table name
    */
    table<int>* getIntTable(std::string n);

    /*!
       \brief Get tables
    */
    std::vector<table<double>*> getTables();

    /*!
       \brief Get the integer tables
    */
    std::vector<table<int>*> getIntTables();

    /*!
       \brief Delete all tables
    */
    void clear();

protected:
    //! sheet name
    std::string itsName;

    //! table pointer
    table<double>* pTable;

    //! table pointer
    table<int>* pIntTable;

    //! list of table pointers
    std::vector<table<double>*> tableList;

    //! table iterator
    typedef std::vector<table<double>*>::iterator tableIterator;

    //! list of int table pointers
    std::vector<table<int>*> intTableList;

    //! int table iterator
    typedef std::vector<table<int>*>::iterator intTableIterator;

    //! list of string table pointers
    std::vector<table<std::string>*> stringTableList;

    //! string table iterator
    typedef std::vector<table<std::string>*>::iterator stringTableIterator;
};

}

#endif // SHEET_H
