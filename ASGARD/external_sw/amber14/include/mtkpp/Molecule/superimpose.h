/*!
   \file superimpose.h
   \brief Superimposes molecules
   \author Martin Peters

   $Date: 2010/03/29 20:45:26 $
   $Revision: 1.14 $

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

#ifndef SUPERIMPOSE_H
#define SUPERIMPOSE_H

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include "Utils/constants.h"

// - BOOST - //
/*
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/banded.hpp> // for diagonal matrix
#include <boost/numeric/bindings/blas/blas3.hpp>
#include "boost/numeric/bindings/traits/ublas_matrix.hpp"
#include "boost/numeric/bindings/traits/ublas_vector.hpp"

namespace ublas  = boost::numeric::ublas;
namespace blas   = boost::numeric::bindings::blas;
*/
#include <Eigen/Dense>
using namespace Eigen;

namespace MTKpp
{

class molecule;
class vector3d;

// ============================================================
// Class : superimpose()
// ------------------------------------------------------------
/*!
   \class superimpose
   \brief Functions to align molecules
   \author Martin Peters
   \date 2006
*/
// ============================================================
class superimpose
{
public:

    //! superimpose Constructor
    superimpose();

    /*!
       \brief superimpose Constructor
       \param i number of atoms
    */
    superimpose(int i);

    //! superimpose Destructor.
    virtual ~superimpose();

    /*!
      \brief Superimposes molecule B onto molecule A and reports an rmsd
      \param pMoleculeA Fixed molecule
      \param pMoleculeB Molecule to be moved
      \return rmsd value
    */
    double    fit(molecule* pMoleculeA, molecule* pMoleculeB);

    /*!
      \brief Superimposes coordinates B onto coordinates A and reports an rmsd
      \param coordsA Fixed coordinates
      \param coordsB Coordinates to be moved
      \param n Number of coordinates
      \return rmsd value
    */
    double    fit(double coordsA[][3], double coordsB[][3], int n);

    /*!
      \brief Superimposes molecule B onto molecule A using a correspondence vector
      \param pMoleculeA Molecule A
      \param pMoleculeB Molecule B
      \param cor correspondence
      \param type of matching
          - 0 match based on atom type
          - 1 match based on atom symbol
          - 2 match based on heavy atom type
          - 3 match based on heavy atom symbol
      \return success
    */
    int       fit(molecule* pMoleculeA, molecule* pMoleculeB,
                  std::vector<int> cor, int type);

    /*!
      \brief Computes the rmsd between molecule A and molecule B however
             does not transform the coordinates after the fit
      \param pMoleculeA Fixed molecule
      \param pMoleculeB Molecule to be moved
      \return rmsd value
    */
    double    rmsd(molecule* pMoleculeA, molecule* pMoleculeB);

    /*!
      \brief Computes the lowest rmsd between molecule A and molecule B
             using the correspondence matrices and does not transform the
             coordinates after the fit
      \param pMoleculeA Fixed molecule
      \param pMoleculeB Molecule to be moved
      \param correspondenceMatrices Correspondence Matrices
      \return rmsd value
    */
    double    rmsd(molecule* pMoleculeA, molecule* pMoleculeB,
              std::vector<std::vector<int> > &correspondenceMatrices);

    /*!
      \brief Computes the lowest rmsd between molecule A and molecule B
             using the correspondence matrices and does not transform the
             coordinates after the fit
      \param pMoleculeA Fixed molecule
      \param molBCoords Molecule to be moved coordinates
      \param correspondenceMatrices Correspondence Matrices
      \param cor final correspondence
      \return rmsd value
    */
    double    rmsd(molecule* pMoleculeA, std::vector< vector3d > &molBCoords,
              std::vector<std::vector<int> > &correspondenceMatrices, int &cor);

    /*!
      \brief Computes the lowest rmsd between molecule A and molecule B
             using the correspondence matrices
      \param pMoleculeA Fixed molecule
      \param molBCoords Molecule to be moved coordinates
      \param correspondenceMatrices Correspondence Matrices
      \return rmsd value
    */
    double    rmsdNoFit(molecule* pMoleculeA,
              std::vector< vector3d > &molBCoords,
              std::vector<std::vector<int> > &correspondenceMatrices);

    /*!
      \brief Get Rotation matrix and centers of A and B
      \param centerA Center of A
      \param centerB Center of B
      \param rotMat Rotation Matrix
    */
    void       getRotationMatrix(double centerA[3], double centerB[3],
                                 double rotMat[][3]);

    /*!
      \brief Get Rotation matrix which superimpose B onto A
      \param coordsA Fixed coordinates
      \param coordsB Coordinates to be moved
      \param n Number of coordinates
      \param rotMat Rotation Matrix
      \return success
    */
    int       getRotationMatrix(double coordsA[][3], double coordsB[][3],
                                     int n, double rotMat[][3]);

    /*!
      \brief Get Rotation matrix which superimpose B onto A
      \param coordsA Fixed coordinates
      \param coordsB Coordinates to be moved
      \param centerA center of B
      \param centerB center of A
      \param n Number of coordinates
      \param rotMat Rotation Matrix
      \return success
    */
    int       getRotationMatrix(double coordsA[][3], double coordsB[][3],
                                double centerA[3], double centerB[3],
                                int n, double rotMat[][3]);

    /*!
      \brief Calculate rmsd between two sets of coordinates
      \param coordsA Coordinate set A
      \param coordsB Coordinate set B
      \return rmsd between coordsA and coordsB
    */
    double    calculateRMSD(double coordsA[][3], double coordsB[][3]);

    /*!
      \brief Calculate rmsd between two sets of coordinates
      \param pMoleculeA molecule A
      \param pMoleculeB molecule B
      \param type of matching
          -0 match based on atom type
          -1 match based on atom symbol
      \param nFittedAtoms number of atoms used to determine the rmsd
      \return minimum rmsd between two different molecule
    */
    double    minRMSD(molecule* pMoleculeA, molecule* pMoleculeB, int type,
                      int &nFittedAtoms);

    /*!
       \brief Initialize the correspondences between the atoms in a molecule
       \param pMol molecule pointer
       \param type of matching
          - 0 match based on atom type
          - 1 match based on atom symbol
          - 2 match based on heavy atom type
          - 3 match based on heavy atom symbol
       \param correspondenceMatrices Correspondence Matrices
       \return success
    */
    int initializeCorrespondences(molecule* pMol,
                 int type,
                 std::vector<std::vector<int> > &correspondenceMatrices);

    int initializeCorrespondences(molecule* pMolB,
                 molecule* pMolA,
                 int type,
                 std::vector<std::vector<int> > &correspondenceMatrices);

    /*!
       \brief Calculate the center of mass
       \param coords Coordinate set
       \param center Center of mass
    */
    void      center(double coords[][3], double center[3]);

    /*!
       \brief Update coordinates
       \param coords Coordinate set
       \param n number of atoms
       \param c1 first center
       \param c2 second center
       \param rotMat rotation matrix
    */
    void      updateCoords(double coords[][3], int n, double c1[3],
                           double c2[3], double rotMat[][3]);

protected: // functions

    /*!
      \brief
      \param coordsA Coordinate set A
      \param coordsB Coordinate set B
      \param itsDXM Work array
      \param itsDXP Work array
    */
    void      coordinateDiff(double coordsA[][3], double coordsB[][3],
                             double itsDXM[][3], double itsDXP[][3]);

    /*!
      \brief Constructs the quaternion matrix
      \param quaternion Quaternion matrix
      \param itsDXM Work array
      \param itsDXP Work array
    */
/*
    void      buildQuaternion(ublas::matrix<double, ublas::column_major>& quaternion,
                              double itsDXM[][3], double itsDXP[][3]);
*/
    void      buildQuaternion(Eigen::MatrixXd& quaternion,
                              double itsDXM[][3], double itsDXP[][3]);

    /*!
      \brief Builds rotation matrix, t, from the decomposed Quaternion matrix
      \param t Rotation matrix
      \param eigenvectors Decomposed Quaternion matrix
    */
/*
    void      buildRotation(double t[][3],
              ublas::matrix<double, ublas::column_major> eigenvectors);
*/
    void      buildRotation(double t[][3],
              Eigen::MatrixXd eigenvectors);

protected: // data

    //! molecule pointer
    molecule* pMoleculeA;

    //! molecule pointer
    molecule* pMoleculeB;

    //! Number of atoms
    int       nAtoms;

    //! Center of mass, MoleculeA
    double    centerA[3];

    //! Center of mass, MoleculeB
    double    centerB[3];

    //! Rotation matrix
    double    rotMat[3][3];

    //! Root mean squared deviation
    double    dRMSD;

    //! match matrix
    int* genMatchMatrix;

    /*!
       Correspondence Type
        - 0
        - 1
        - 2
        - 3
    */
    int correspondenceType;

    //! Number of Heavy Atoms
    int                           nHeavyAtomsA;

    //! Heavy Atom indices
    int*                          heavyAtomIndicesA;

    //! Number of Heavy Atoms
    int                           nHeavyAtomsB;

    //! Heavy Atom indices
    int*                          heavyAtomIndicesB;

};

} // MTKpp namespace

#endif // SUPERIMPOSE_H
