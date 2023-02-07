/*!
   \file table.h
   \brief Extension of ublas::matrix to store labels
   \author Martin Peters

   $Date: 2010/03/29 20:35:21 $
   $Revision: 1.9 $

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

#ifndef TABLE_H
#define TABLE_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>
#include <math.h>

#include "Diagnostics/MTKException.h"

#include <Eigen/Dense>
using namespace Eigen;

namespace MTKpp
{
// ============================================================
// Class : table()
// ------------------------------------------------------------
/*! 
   \class table
   \brief Extension of eigen::matrix to store labels
   \author Martin Peters
   \date 2007
*/
// ============================================================
template <class T> class table
{
public:
    /*!
       \brief table Constructor
    */
    table() {
      this->itsName = "";
      this->nColumns = 0;
      this->nRows = 0;
    }

    //! table Destructor
    virtual ~table() {}

    /*!
       \brief Set the name of the table
       \param n table name
    */
    void setName(std::string n) {
      this->itsName = n;
    }

    /*!
       \brief Get the name of the table
       \return table name
    */
    std::string getName() {
      return this->itsName;
    }

    /*!
       \brief Set the type of the table
       \param n table name
    */
    void setType(std::string n) {
      this->itsType = n;
    }

    /*!
       \brief Get the type of the table
       \return table name
    */
    std::string getType() {
      return this->itsType;
    }

    /*!
       \brief Set the number of columns in the table
       \param i number of columns
    */
    void setNumColumns(const int& i) {
      this->nColumns = i;
    }

    /*!
       \brief Get the number of columns in the table
       \return number of columns
    */
    int getNumColumns() {
      return this->nColumns;
    }

    /*!
       \brief Set a column label in the table
       \param i column index
       \param s label
    */
    void setColumnLabel(const int& i, std::string s) {
      if (i > this->nColumns-1) {
        std::stringstream ss;
        ss << " Out of bounds error: i = " << i+1
                  << " > " << this->nColumns << std::endl;
        std::cout << ss.str();
        throw MTKException(ss.str());
      }
      this->columnLabels[i] = s;
    }

    /*!
       \brief Set column labels
       \param v column labels
    */
    void setColumnLabels(std::vector<std::string> v) {
      if (v.size() >= this->columnLabels.size()) {
        for (unsigned int i = 0; i < columnLabels.size(); i++) {
          this->columnLabels[i] = v[i];
        }
      }
    }

    /*!
       \brief Get column label
       \param i column index
       \return column label
    */
    std::string getColumnLabel(const int& i) {
      if (i < static_cast<int>(this->columnLabels.size())) {
        return this->columnLabels[i];
      }
      return "";
    }

    /*!
       \brief Get column labels
       \return column labels
    */
    std::vector<std::string> getColumnLabels() {
      return this->columnLabels;
    }

    /*!
       \brief Set row labels
       \param v row labels
    */
    void setRowLabels(std::vector<std::string> v) {
      if (v.size() == this->rowLabels.size()) {
        for (unsigned int i = 0; i < v.size(); i++) {
          this->rowLabels[i] = v[i];
        }
      }
    }

    /*!
       \brief Get row label
       \param i row index
       \return row labels
    */
    std::string getRowLabel(const int& i) {
      if (i < static_cast<int>(this->rowLabels.size())) {
        return this->rowLabels[i];
      }
      return "";
    }

    /*!
       \brief Get row labels
       \return row labels
    */
    std::vector<std::string> getRowLabels() {
      return this->rowLabels;
    }

    /*!
       \brief Set the number of rows in the table
       \param i number of rows
    */
    void setNumRows(const int& i) {
      this->nRows = i;
    }

    /*!
       \brief Get the number of rows in the table
       \return number of rows
    */
    int getNumRows() {
      return this->nRows;
    }

    /*!
       \brief Set a row label in the table
       \param i row index
       \param s label
    */
    void setRowLabel(const int& i, std::string s) {
      if (i > this->nRows-1) {
        std::stringstream ss;
        ss << " Out of bounds error: i = " << i+1
                  << " > " << this->nRows << std::endl;
        std::cout << ss.str();
        throw MTKException(ss.str());
      }
      this->rowLabels[i] = s;
    }

    /*!
       \brief Set the number of columns and rows in the table
       \param i number of rows
       \param j number of columns
    */
    void setSizes(const int& i, const int& j) {
      this->nRows = i;
      this->nColumns = j;
      this->setSize(i, j);
    }

    /*!
       \brief table setup
    */
    void setup() {
      if (this->nColumns > 0 and this->nRows > 0) {
        this->setSize(this->nRows, this->nColumns);
      }
      else {
        std::cout << " ERROR " << std::endl;
        throw MTKException(" table::setup ERROR ");
      }
    }

    /*!
       \brief Set the value of a cell in the table
       \param i row index
       \param j column index
       \param v cell value
    */
    void setCellValue(const int& i, const int& j, T v) {
      if ((i > this->nRows-1) or (j > this->nColumns-1)) {
        std::stringstream ss;
        ss << " Out of bounds error: i = " << i << " j = " << j << " value = " << v
                  << " " << this->nRows << " " << this->nColumns << std::endl;
        std::cout << ss.str();
        throw MTKException(ss.str());
      }
      else {
        this->itsMatrix(i,j) = v;
      }
    }

    /*!
       \brief Get the value of a cell in the table
       \param i row index
       \param j column index
       \return cell value
    */
    T getCellValue(const int& i, const int& j) {
      if ((i > this->nRows-1) or (j > this->nColumns-1)) {
        std::stringstream ss;
        ss << " Out of bounds error: i = " << i << " j = " << j
                  << " " << this->nRows << " " << this->nColumns << std::endl;
        std::cout << ss.str();
        throw MTKException(ss.str());
      }
      else {
        return this->itsMatrix(i,j);
      }
      return 0;
    }

    /*!
       \brief Print table to the screen
    */
    void print() {
      std::cout << " \n " << this->getName() << std::endl;
      std::cout << "- ";
      for (unsigned j = 0; j < this->itsMatrix.cols(); ++ j) {
        std::cout << this->columnLabels[j] << " ";
      }
      std::cout << " " << std::endl;

      for (unsigned i = 0; i < this->itsMatrix.rows(); ++ i) {
        std::cout << this->rowLabels[i] << " ";
        for (unsigned j = 0; j < this->itsMatrix.cols(); ++ j) {
          std::cout << std::showpoint << this->itsMatrix(i,j) << " ";
        }
        std::cout << " " << std::endl;
      }
    }

    /*!
       \brief Print row of table to the screen
    */
    void printRow(const int& r) {
      for (unsigned j = 0; j < this->itsMatrix.cols(); ++ j) {
        std::cout << this->columnLabels[j] << " ";
      }
      std::cout << " " << std::endl;
      if (r > static_cast<int>(this->itsMatrix.rows())) {
        std::cout << " Error in table " << std::endl;
      }
      std::cout << this->rowLabels[r] << " ";
      for (unsigned j = 0; j < this->itsMatrix.cols(); ++ j) {
        std::cout << std::showpoint << this->itsMatrix(r,j) << " ";
      }
      std::cout << " " << std::endl;
    }

    /*!
       \brief Print table matrix to the screen
    */
    void printMatrix() {
      for (unsigned i = 0; i < this->itsMatrix.rows(); ++ i) {
        for (unsigned j = 0; j < this->itsMatrix.cols(); ++ j) {
          std::cout << std::showpoint << this->itsMatrix(i,j) << " ";
        }
        std::cout << " " << std::endl;
      }
    }

    /*!
       \brief Get the matrix of type M which the table stores
    */
    Eigen::Matrix<T, Dynamic, Dynamic>& getMatrix() {
      return this->itsMatrix;
    }

    /*!
       \brief Print table matrix to the screen
       \param j column to be sorted
       \param order Ascending = 0, Desending = 1
    */
    void sortByColumn(int j, int order) {

      Eigen::Matrix<T, Dynamic, Dynamic> bkupMatrix;
      bkupMatrix.resize(this->nRows, this->nColumns);

      for (unsigned i = 0; i < this->itsMatrix.rows(); ++i) {
        for (unsigned k = 0; k < this->itsMatrix.cols(); ++k) {
          bkupMatrix(i,k) = this->itsMatrix(i,k);
        }
      }

      int size = nRows;
      //int size2 = nColumns;

//      Eigen::Vector<T> column;
//      Eigen::Vector<int> columnIndices;

      Eigen::Matrix<T, Dynamic, 1> column;
      Eigen::Matrix<int, Dynamic, 1> columnIndices;

      column.resize(size,1);
      columnIndices.resize(size, 1);

      int k = 0;
      double p = 0.0;
      std::string label = "";

      for (int i = 0; i < size; i++) {
        column(i, 0) = itsMatrix(i,j);
        columnIndices(i, 0) = i;
      }

      for (int i = 0; i < size; ++i) {
        k = i;
        p = column(i, 0);
        label = rowLabels[i];

        for (int j = i+1; j < size; ++j) {
          if (!order) { // Ascending
            if (column(j, 0) < p) {
              k = j;
              p = column(j, 0);
              label = rowLabels[j];
            }
          }
          else { // Descending
            if (column(j,0) > p) {
              k = j;
              p = column(j,0);
              label = rowLabels[j];
            }
          }
        }
        if ( k != i ) {
          column(k,0) = column(i,0);
          column(i,0) = p;
          int temK = columnIndices(k,0);
          columnIndices(k,0) = columnIndices(i,0);
          columnIndices(i,0) = temK;
          rowLabels[k] = rowLabels[i];
          rowLabels[i] = label;
/*
          for (int m = 0; m < size2; ++m) {
            p = itsMatrix(i,m);
            itsMatrix(i,m) = itsMatrix(k,m);
            itsMatrix(k,m) = p;
          }
*/
        }
      }
      for (unsigned i = 0; i < this->itsMatrix.rows(); ++i) {
        for (unsigned k = 0; k < this->itsMatrix.cols(); ++k) {
          this->itsMatrix(i,k) = bkupMatrix(columnIndices(i,0), k);
        }
      }
    }

    /*!
       \brief Initialize
    */
    void initialize(double t) {
      for (unsigned i = 0; i < this->itsMatrix.rows(); ++ i) {
        for (unsigned j = 0; j < this->itsMatrix.cols(); ++ j) {
          this->itsMatrix(i,j) = t;
        }
      }
    }

    /*!
       \brief Initialize
    */
    void initialize(int t) {
      for (unsigned i = 0; i < this->itsMatrix.rows(); ++ i) {
        for (unsigned j = 0; j < this->itsMatrix.cols(); ++ j) {
          this->itsMatrix(i,j) = t;
        }
      }
    }

    /*!
       \brief Initialize
    */
    void initialize(std::string t) {
      for (unsigned i = 0; i < this->itsMatrix.rows(); ++ i) {
        for (unsigned j = 0; j < this->itsMatrix.cols(); ++ j) {
          this->itsMatrix(i,j) = t;
        }
      }
    };

protected: // Functions
    /*!
       \brief Set the dimensions of the table
    */
    void setSize(const int& i, const int& j) {
//      try {
        this->itsMatrix.resize(i,j);
        this->rowLabels.resize(i);
        this->columnLabels.resize(j);
        this->initializeAll();
/*
      }
      catch (boost::numeric::ublas::bad_size) {
        std::cout << " ERROR IN TABLE ... EXITING " << std::endl;
        throw MTKException(" ERROR IN TABLE ... EXITING ");
      }
*/
    }

    void initializeAll() {
      for (unsigned i = 0; i < this->itsMatrix.rows(); ++ i) {
        this->rowLabels[i] = "";
      }
      for (unsigned j = 0; j < this->itsMatrix.cols(); ++ j) {
        this->columnLabels[j] = "";
      }
    }

protected: // Data
    //! Name
    std::string itsName;

    //! Type
    std::string itsType;

    //! Number of columns
    int nColumns;

    //! Columns labels
    std::vector<std::string> columnLabels;

    //! Number of rows
    int nRows;

    //! Row labels
    std::vector<std::string> rowLabels;

    //! matrix of type T
    Eigen::Matrix<T, Dynamic, Dynamic> itsMatrix;
};

}

#endif // TABLE_H
