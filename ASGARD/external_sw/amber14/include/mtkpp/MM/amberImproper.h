/*! 
   \file amberImproper.h
   \brief AMBER improper energy and gradient
   \author Martin Peters

   $Date: 2010/03/29 20:28:34 $
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

#ifndef AMBERIMPROPER_H
#define AMBERIMPROPER_H

namespace MTKpp
{

class amber;

// ============================================================
// Class : amberImproper()
// ------------------------------------------------------------
/*! 
   \class amberImproper

   \brief AMBER improper energy and gradient

   \author Martin Peters

   \version 0.1

   \date 2006

   \section AMBERImproperEnG AMBER Torsion Energy And Gradient

*/
// ============================================================
class amberImproper
{
public:
    /*!
       \brief amberImproper Constructor
    */
    amberImproper();

    /*!
       \brief amberImproper Constructor
       \param pAmber parent potential
    */
    amberImproper(amber* pAmber);

    /*!
       \brief amberImproper Destructor
    */
    virtual ~amberImproper();

    /*!
       \brief Calculate angle Energy
       \return angle energy
    */
    double calculateE();

    /*!
       \brief Calculate angle Energy and Gradient
       \return angle energy
    */
    double calculateG();

protected: // FUNCTIONS
    /*!
       \brief Determine the improper between atoms 1, 2, 3, and 4
       \param x1 atom 1 x coordinate
       \param y1 atom 1 y coordinate
       \param z1 atom 1 z coordinate
       \param x2 atom 2 x coordinate
       \param y2 atom 2 y coordinate
       \param z2 atom 2 z coordinate
       \param x3 atom 3 x coordinate
       \param y3 atom 3 y coordinate
       \param z3 atom 3 z coordinate
       \param x4 atom 4 x coordinate
       \param y4 atom 4 y coordinate
       \param z4 atom 4 z coordinate
       \return improper between atoms 1, 2, 3, and 4
    */
    double improper(double x1, double y1, double z1, double x2, double y2, double z2,
                    double x3, double y3, double z3, double x4, double y4, double z4);

protected: // DATA
    //! improper energy
    double         energy;

    //! parent potential
    amber*         pAmber;
};

} // MTKpp namespace

    //double calculateE(collection*);
    //double calculateE(molecule*);
    //double calculateG(collection*);
    //double calculateG(molecule*);
    //double calculateG(submolecule*);
    //typedef std::map<int, Improper*>::iterator ImproperMapIterator;
    //std::map<int, Improper*>      moleculeImproperMap;
    //std::vector<molecule*>    moleculeList;

    //collection*    pCollection;
    //molecule*      pMolecule;
    //submolecule*   pSubMolecule;
    //Improper*      pImproper;
    //parameters*    pParameters;
    //improperParam* pImproperParam;
    //std::vector<improperParam*> improperParamList;
    //typedef std::vector<improperParam*>::iterator improperParamIterator;

#endif // amberImproper_H

