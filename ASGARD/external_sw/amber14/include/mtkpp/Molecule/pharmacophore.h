/*!
   \file pharmacophore.h
   \brief Determines the Maximum Common Pharmacophore (MCP) between two molecules
   \author Martin Peters

   $Date: 2010/03/29 20:44:27 $
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

#ifndef PHARMACOPHORE_H
#define PHARMACOPHORE_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include "Utils/vector3d.h"

namespace MTKpp
{

class molecule;
class atom;
//class vector3d;

// ============================================================
// Class : clique
// ------------------------------------------------------------
/*!
   \struct clique
   \brief Container for clique info
*/
// ============================================================
class clique
{
public:
    //! pharmacophore Constructor
    clique();

    //! pharmacophore Destructor
    virtual ~clique();

    /*!
       \brief Compares two cliques based on score
       \param lhs first clique
       \param rhs Second clique
       \return boolean

       After this function is defined the STL sort() function can be used
    */
    static bool less(const clique *lhs, const clique *rhs) {
      return lhs->score < rhs->score;
    }

    /*!
       \brief Compares two cliques based on score
       \param lhs first clique
       \param rhs Second clique
       \return boolean

       After this function is defined the STL sort() function can be used
    */
    static bool greater(const clique *lhs, const clique *rhs) {
      return lhs->score > rhs->score;
    }

    /*!
       \brief Checks equality of two vectors
       \param lhs Left hand side vector3d object
       \param rhs Right hand side vector3d object
    */
    inline friend bool operator==(clique &lhs, clique &rhs) {
      unsigned int nHits = 0;
      unsigned int counter = 0;
      unsigned int lhsSize = lhs.indicesA.size();
      unsigned int rhsSize = rhs.indicesA.size();
      if (lhsSize != rhsSize) return false;

      // Not sure about this test, may remove some useful cliques
      if (std::abs(lhs.score - rhs.score) < 0.00001) {
        return true;
      }

      for (unsigned int i = 0; i < lhsSize; i++) {
        for (unsigned int j = 0; j < rhsSize; j++) {
          if (lhs.indicesA[i] == rhs.indicesA[j]) {
            if (lhs.indicesB[i] == rhs.indicesB[j]) {
              counter++;
            }
          }
        }
      }

      for (unsigned int i = 0; i < lhsSize; i++) {
        for (unsigned int j = 0; j < rhsSize; j++) {
          double coordDiff = lhs.coords[i].dist(rhs.coords[j]);
          if (std::abs(coordDiff) < 0.0001) {
            nHits++;
            break;
          }
        }
      }

      if ((counter == rhsSize) or (nHits == rhsSize)) {
        return true;
      }
      return false;
    }

    /*!
       \brief Compares two cliques
       \return bool
    */
    bool compare(clique* rhs);

    //! score
    double score;

    //! score
    double score2;

    //! dMax
    double dMax;

    //! feature coordinates
    std::vector<vector3d> coords;

    //! feature indices A
    std::vector<int> indicesA;

    //! feature indices B
    std::vector<int> indicesB;
};

// ============================================================
// Class : pharmacophore()
// ------------------------------------------------------------
/*!
   \class pharmacophore

   \brief Determines the Maximum Common Pharmacophore (MCP) between two molecules
*/
// ============================================================
class pharmacophore
{
public:
    /*!
       \brief pharmacophore Constructor
       \param parent molecule pointer
       \param distMax Feature distance parameter
    */
    pharmacophore(molecule *parent, double distMax = 0.5);

    //! pharmacophore Destructor
    virtual ~pharmacophore();

    /*!
       \brief Determines the Maximum Common Pharmacophore (MCP) between two molecules
       \param mol molecule pointer
       \param mcp Maximum Common Pharmacophore
       \param mcpCoords mcp coordinates
       \return success
    */
    int run(molecule* mol,
            std::vector<std::vector<unsigned int> > &mcp,
            std::vector<vector3d> &mcpCoords);

    /*!
       \brief Determines all possible Pharmacophores between two molecules
       \param mol molecule pointer
       \param cliqueList all possible pharmacophores
       \return success
    */
    int getCliques(molecule* mol, std::vector<clique*> &cliqueList);

protected: // DATA

    /*!
       \brief Generate the correspondence matrix between molecule A and B
       \param nFeaturesA Number of feature in mol A
       \param nFeaturesB Number of feature in mol B
       \param featureLabelsA Feature labels in mol A
       \param featureLabelsB Feature labels in mol B
       \param crpdeMatrix Correspondence matrix
       \todo determine the doxygen error here
       \return success
    */
    int generateCorrMatrix(unsigned int nFeaturesA, unsigned int nFeaturesB,
                           std::vector<std::string> featureLabelsA,
                           std::vector<std::string> featureLabelsB,
                           int crpdeMatrix[]);
//                           int* crpdeMatrix);

    /*!
       \brief Get the c good value
       \param cGood c good value
       \param nFeaturesA Number of feature in mol A
       \param nFeaturesB Number of feature in mol B
       \param featureLabelsA Feature labels in mol A
       \param featureDistMatrixA
       \param featureDistMatrixB
       \param crpdeMatrix Correspondence matrix
       \return success
    */
    int getCGood(double &cGood, unsigned int nFeaturesA, unsigned int nFeaturesB,
                 std::vector<std::string> featureLabelsA,
                 double* featureDistMatrixA,
                 double* featureDistMatrixB,
                 int crpdeMatrix[]);

    /*!
       \brief Get clique
       \param iA ?
       \param jB ?
       \param nFeaturesA Number of feature in mol A
       \param nFeaturesB Number of feature in mol B
       \param featureLabelsA Feature labels in mol A
       \param featureDistMatrixA
       \param featureDistMatrixB
       \param crpdeMatrix Correspondence matrix
       \param cGood ?
       \param clique ?
       \param score ?
       \param score2 ?
       \param dMax ?
       \return success
    */
    int getClique(unsigned int iA, unsigned int jB,
                   unsigned int nFeaturesA, unsigned int nFeaturesB,
                   std::vector<std::string> featureLabelsA,
                   double* featureDistMatrixA, double* featureDistMatrixB,
                   int crpdeMatrix[], double cGood,
                   std::vector<std::vector<unsigned int> > &clique,
                   double &score, double &score2, double &dMax);

    //---------------//
    // -  POINTERS - //
    //---------------//

    //! molecule pointer
    molecule*                pParent;

    //! atom pointer
    atom*                    pAtom;

    //! ddMax
    double                   deltaDistMax;
};

} // MTKpp namespace

#endif // PHARMACOPHORE_H
