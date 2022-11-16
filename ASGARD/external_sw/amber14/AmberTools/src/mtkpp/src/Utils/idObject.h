/*!
   \file idObject.h
   \brief 
   \author Martin B. Peters

   $Date: 2010/03/29 20:33:22 $
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

#ifndef IDOBJECT_H
#define IDOBJECT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iomanip>

namespace MTKpp
{

// ============================================================
// Class : idObject()
// ------------------------------------------------------------
/*! 
   \class idObject
   \brief
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================
class idObject
{
public:

    /*!
       \brief idObject Constructor
    */
    idObject() {
      this->i = 0;
      this->d = 0.0;
    }

    /*!
       \brief idObject Constructor
    */
    idObject(int i, double d) {
      this->i = i;
      this->d = d;
    }

    /*!
       \brief idObject Constructor
    */
    idObject(int i, double d, std::string n) {
      this->i = i;
      this->d = d;
      this->n = n;
    }

    /*!
       \brief idObject Destructor
    */
    virtual ~idObject() {};

    /*!
       \brief Compares two idObject based on d
       \param lhs first idObject
       \param rhs Second idObject
       \return boolean

       After this function is defined the STL sort() function can be used
    */
    static bool less(const idObject *lhs, const idObject *rhs) {
        return lhs->d < rhs->d;
    }

    /*!
       \brief Compares two idObject based on d
       \param lhs first idObject
       \param rhs Second idObject
       \return boolean

       After this function is defined the STL sort() function can be used
    */
    static bool greater(const idObject *lhs, const idObject *rhs) {
        return lhs->d > rhs->d;
    }

    void setI(int i) {
      this->i = i;
    }

    void setD(double d) {
      this->d = d;
    }

    void setN(std::string n) {
      this->n = n;
    }

    int getI() {
      return this->i;
    }

    double getD() {
      return this->d;
    }

    std::string getN() {
      return this->n;
    }

protected:

   //! index
   int i;

   //! numerical value
   double d;

   //! object name
   std::string n;
};

// ============================================================
// Class : idObjectList()
// ------------------------------------------------------------
/*! 
   \class idObjectList
   \brief
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================
class idObjectList
{
public:

    //! index
    int index;

    //! list of idObject
    std::vector<idObject*> objList;

public:

    /*!
       \brief idObjectList Constructor
    */
    idObjectList() {
      this->index = 0;
    }

    /*!
       \brief idObject Destructor
    */
    virtual ~idObjectList() {};

    void addObject(std::string name, double value) {
      idObject* idO = new idObject(index, value, name);
      objList.push_back(idO);
      index++;
    }

    void sortList(int d) {
      if (d > 0) { // ascending
        std::sort(objList.begin(), objList.end(), idObject::less);
      }
      else { // descending
        std::sort(objList.begin(), objList.end(), idObject::greater);
      }
    }

    int getRank(std::string n) {
      for (unsigned int i = 0; i < objList.size(); i++) {
        if (objList[i]->getN() == n) {
          return i;
        }
      }
      return -1;
    }

    idObject* getObject(std::string n) {
      if (index == 0) return 0;
      for (unsigned int i = 0; i < objList.size(); i++) {
        if (objList[i]->getN() == n) {
          return objList[i];
        }
      }
      return 0;
    }

    idObject* getFirst() {
      return objList[0];
    }

    void print() {
      std::cout << " Number of objects = " << objList.size() << std::endl;
      for (unsigned int i = 0; i < objList.size(); i++) {
        std::cout << " Name = " << objList[i]->getN() << " Value = " << objList[i]->getD() << std::endl;
      }
    }

    double lMin() {
      double m = objList[0]->getD();
      for (unsigned int i = 1; i < objList.size(); i++) {
        if (objList[i]->getD() < m) {
          m = objList[i]->getD();
        }
      }
      return m;
    }

    double lMax() {
      double m = objList[0]->getD();
      for (unsigned int i = 1; i < objList.size(); i++) {
        if (objList[i]->getD() > m) {
          m = objList[i]->getD();
        }
      }
      return m;
    }

    double mean() {
      double m = 0.0;
      for (unsigned int i = 0; i < objList.size(); i++) {
        m += objList[i]->getD();
      }
      m /= double(objList.size());
      return m;
    }

    double variance() {
      if (index == 0) return 0.0;
      double variance = 0.0;
      double m = this->mean();

      for (unsigned int i = 0; i < objList.size(); i++) {
        double x = objList[i]->getD() - m;
        variance += (x * x);
      }
      variance /= double(objList.size() - 1);
      return variance;
    }

    double standardDeviation() {
      if (index == 0) return 0.0;
      return sqrt(this->variance());
    }

    void normalize() {
      double mn = lMin();
      double mx = lMax();
      double mxmn = mx - mn;

      for (unsigned int i = 0; i < objList.size(); i++) {
        objList[i]->setD(  (objList[i]->getD() - mn) / mxmn  );
      }
    }

    void autoScale() {
      double m = this->mean();
      double s = this->standardDeviation();

      for (unsigned int i = 0; i < objList.size(); i++) {
        objList[i]->setD(  (objList[i]->getD() - m) / s  );
      }
    }

};

} // MTKpp namespace

#endif // IDOBJECT_H

