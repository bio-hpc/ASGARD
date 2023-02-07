/*!
   \file mmPotential.h
   \brief Base class for MM potentials
   \author Martin Peters

   \todo Use templates in MM Library

   $Date: 2010/03/29 20:28:34 $
   $Revision: 1.11 $

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

#ifndef MMPOTENTIAL_H
#define MMPOTENTIAL_H

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>

#ifdef HAVE_ZLIB
#include "zlib.h"
#endif

#include "Utils/constants.h"

/*!
\namespace MTKpp
\brief MTK++ namespace
*/
namespace MTKpp
{
// ============================================================
// Class : mmPotential()
// ------------------------------------------------------------
/*! 
   \class mmPotential

   \brief Base class for MM potentials

   \author Martin Peters

   \version 0.1

   \date 2005

   \section usefulMath Some useful math for MM potentials

   \f{eqnarray*}
    r_{ij}      &=& r_i - r_j \\
    |r_{ij}|    &=& \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2} \\
    |r_{ij}|    &=& \left((x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)\right)^{1/2} \\
    |r_{ij}|    &=& \sqrt{r_i \cdot r_j} \\
    \hat r_{ij} &=& {r_{ij} \over |r_{ij}|}
   \label{dist:r}
   \f}

   where \f$r_i\f$ is the position of atom \f$i\f$ and \f$x_i\f$ is the \f$x\f$ component of the position vector of atom \f$i\f$.


   \f{eqnarray*}
    {d\theta \over d\cos\theta} & = & \left({d\cos\theta \over d\theta}\right)^{-1} \\
                            & = & -{1 \over \sin\theta}
   \label{dt:dcost}
   \f}

   \f{eqnarray*}
    \nabla_i = \hat x {\partial \over \partial x_i} + \hat y {\partial \over \partial y_i} + \hat z {\partial \over \partial z_i}
   \f}
   where \f$\hat x\f$ is a unit vector parallel to axis of the reference coordinate system.

*/
// ============================================================
class mmPotential
{
public:
    /*!
       \brief mmPotential Constructor
    */
    mmPotential();

    /*!
       \brief mmPotential Destructor
    */
    virtual ~mmPotential();

    /*!
       \brief mmPotential initializer
    */
    void initialize();


    /// START ALLOCATION ///

    /*!
       \brief Set number of atoms
       \param nAtoms number of atoms
       \return success
    */
    int setNumAtoms(int nAtoms);

    /*!
       \brief Set number of distinct atom type
       \param n number of distinct atom type
       \return success
    */
    int setNumTypes(int n);

    /*!
       \brief Set number of residues
       \param n number of residues
       \return success
    */
    int setNumResidues(int n);

    /*!
       \brief Set number of excluded atoms
       \param n number of excluded atoms
       \param n2 number of excluded 1-4 atoms
       \return success
    */
    int setExcludedSize(int n, int n2);

    /// END ALLOCATION ///

    /*!
       \brief Get number of atoms
       \return number of atoms
    */
    int getNumAtoms();

    /*!
       \brief Get residue names
       \return residue names
    */
    char* getResidueNames();

    /*!
       \brief Get residue pointers
       \return residue pointers
    */
    int* getResiduePointers();

    /*!
       \brief Set number of atoms in the largest residue
       \param n number of atoms in the largest residue
       \return success
    */
    int setNumAtomsBigRes(int n);

    /*!
       \brief Get coordinate array
       \return coordinate array pointer
    */
    double* getCoords();

    /*!
       \brief Get atomic charge array pointer
       \return charges array pointer
    */
    double* getCharges();

    /*!
       \brief Get atomic symbols array pointer
       \return symbols array pointer
    */
    char* getSymbols();

    /*!
       \brief Get atomic symbols array pointer
       \return symbols array pointer
    */
    char* getAtomNames();

    /*!
       \brief Get number of unique types
       \return number of unique types
    */
    int getNumTypes();

    /*!
       \brief Get atom integer types
       \return integer types array pointer
    */
    int* getIntTypes();

    /*!
       \brief Get atom character types
       \return character types array pointer
    */
    char* getCharTypes();

    /*!
       \brief Get number of excluded atoms array pointer
       \return number of excluded atoms array pointer
    */
    int* getNumExcluded();

    /*!
       \brief Get excluded atoms array pointer
       \return excluded atoms array pointer
    */
    int* getExcluded();

    /*!
       \brief Get number of excluded 1-4 atoms array pointer
       \return number of excluded atoms array pointer
    */
    int* getNumExcluded14();

    /*!
       \brief Get excluded 1-4 atoms array pointer
       \return excluded atoms array pointer
    */
    int* getExcluded14();

    /*!
       \brief Get atom flags array pointer
       \return atom flags array pointer
    */
    int* getAtomFlags();

    /*!
       \brief Get gradient array pointer
       \return gradient array pointer
    */
    double* getGradients();

    /*!
       \brief Set number of bonds
       \param nBonds number of bonds
       \return success
    */
    int setNumBonds(int nBonds);

    /*!
       \brief Get number of bonds
       \return number of bonds
    */
    int getNumBonds();

    /*!
       \brief Set number of bonds containing Hydrogen
       \param n number of bonds containing Hydrogen
       \return success
    */
    int setNumBondsWithH(int n);

    /*!
       \brief Set number of bonds not containing Hydrogen
       \param n number of bonds not containing Hydrogen
       \return success
    */
    int setNumBondsWithOutH(int n);

    /*!
       \brief Get bond array pointer
       \return bond array pointer
    */
    int* getBonds();

    /*!
       \brief Get bond parameter array pointer
       \return bond parameter array pointer
    */
    double* getBondParams();

    /*!
       \brief Set number of unique bonds
       \param n number of unique bonds
       \return success
    */
    int setNumUniqueBonds(int n);

    /*!
       \brief Set number of angles
       \param nAngles number of angles
       \return success
    */
    int setNumAngles(int nAngles);

    /*!
       \brief Get number of angles
       \return number of angles
    */
    int getNumAngles();

    /*!
       \brief Set number of angles containing Hydrogen
       \param n number of angles containing Hydrogen
       \return success
    */
    int setNumAnglesWithH(int n);

    /*!
       \brief Set number of angles not containing Hydrogen
       \param n number of angles not containing Hydrogen
       \return success
    */
    int setNumAnglesWithOutH(int n);

    /*!
       \brief Get angle array pointer
       \return angle array pointer
    */
    int* getAngles();

    /*!
       \brief Get angle parameter array pointer
       \return angle parameter array pointer
    */
    double* getAngleParams();

    /*!
       \brief Set number of unique angles
       \param n number of unique angles
       \return success
    */
    int setNumUniqueAngles(int n);

    /*!
       \brief Set number of torsions
       \param nTorsions number of torsions
       \return success
    */
    int setNumTorsions(int nTorsions);

    /*!
       \brief Get number of torsions
       \return number of torsions
    */
    int getNumTorsions();

    /*!
       \brief Get torsion array pointer
       \return torsion array pointer
    */
    int* getTorsions();

    /*!
       \brief Get torsion parameter array pointer
       \return torsion parameter array pointer
    */
    double* getTorsionParams();

    /*!
       \brief Set number of impropers
       \param nImpropers number of impropers
       \return success
    */
    int setNumImpropers(int nImpropers);

    /*!
       \brief Get number of impropers
       \return number of impropers
    */
    int getNumImpropers();

    /*!
       \brief Get improper array pointer
       \return improper array pointer
    */
    int* getImpropers();

    /*!
       \brief Get improper parameter array pointer
       \return improper parameter array pointer
    */
    double* getImproperParams();

    /*!
       \brief Set number of dihedrals containing Hydrogen
       \param n number of dihedrals containing Hydrogen
       \return success
    */
    int setNumDihedralsWithH(int n);

    /*!
       \brief Set number of dihedrals not containing Hydrogen
       \param n number of dihedrals not containing Hydrogen
       \return success
    */
    int setNumDihedralsWithOutH(int n);

    /*!
       \brief Set number of unique dihedrals
       \param n number of unique dihedrals
       \return success
    */
    int setNumUniqueDihedrals(int n);

    /*!
       \brief Set number of non-bonded interactions
       \param nNonBonded number of non-bonded interactions
       \return success
    */
    int setNumNonBonded(int nNonBonded);

    /*!
       \brief Set number of non-bonded interactions
       \return number of non-bonded interactions
    */
    int getNumNonBonded();

    /*!
       \brief Get nonbonded array pointer
       \return nonbonded array pointer
    */
    int* getNonBonded();

    /*!
       \brief Get nonbonded parameter array pointer
       \return nonbonded parameter array pointer
    */
    double* getNonBondedParams();

    /*!
       \brief Get nonbonded indexing array pointer
       \return nonbonded indexing array pointer
    */
    int* getNonBondedPtrs();

    /*!
       \brief Get 1-4 nonbonded indexing array pointer
       \return 1-4 nonbonded indexing array pointer
    */
    int* getNonBonded14Ptrs();

    /*!
       \brief Get L-J R6 Parameters
       \return L-J R6 Parameters
    */
    double* getR6Params();

    /*!
       \brief Get L-J R12 Parameters
       \return L-J R12 Parameters
    */
    double* getR12Params();

    /*!
       \brief Get nonbonded indexing array pointer
       \return nonbonded indexing array pointer
    */
    int* getNonBondedParameterIndex();

    /*!
       \brief Set the potential to be used
       \param bBond turn on bond calculation
       \param bAngle turn on angle calculation
       \param bTorsion turn on torsion calculation
       \param bImproper turn on improper calculation
       \param bNonBonded turn on non-bonded (Ele & L-J) calculation
       \param bHBond turn on H-bond calculation
    */
    void setPotential(bool bBond, bool bAngle, bool bTorsion,
                      bool bImproper, bool bNonBonded, bool bHBond);

    /*!
       \brief Get the potential to be used
       \param bBond Bond switch
       \param bAngle Angle switch
       \param bTorsion Torsion switch
       \param bImproper Improper switch
       \param bNonBonded Nonbonded switch
       \param bHBond H-Bond switch
    */
    void getPotential(bool& bBond, bool& bAngle, bool& bTorsion,
                      bool& bImproper, bool& bNonBonded, bool& bHBond);

    /*!
       \brief Turn on gradient calculation
       \param i on/off
    */
    void calcForces(int i);

    /*!
       \brief Calculate energy and/or gradient
    */
    void calcEnergy();

    /*!
       \brief Decompose energy
       \param bBond Bond switch
       \param bAngle Angle switch
       \param bTorsion Torsion switch
       \param bImproper Improper switch
       \param bNonBonded Nonbonded switch
       \param bHBond H-Bond switch
       \param f Output file
    */
    void decomposeEnergy(bool bBond, bool bAngle, bool bTorsion,
                         bool bImproper, bool bNonBonded, bool bHBond,
                         std::string f);

    /*!
       \brief Get output flags array for energy decomposition
       \return output flags array
    */
    int* getOutputFlags();

    /*!
       \brief Print energy to screen (verbose)
    */
    void printEnergy();

    /*!
       \brief Print energy (verbose)
       \param os print stream
    */
    void printEnergy(std::ostream& os);

    /*!
       \brief Print energy terms in single line to the screen
    */
    void printEnergy2();

    /*!
       \brief Print energy terms in single line to the screen
       \param os print stream
    */
    void printEnergy2(std::ostream& os);

    /*!
       \brief Print energy terms in single line to the screen
       \param s A string thats added to the begining of the line
    */
    void printEnergy2(std::string s);

    /*!
       \brief Print energy terms in single line to the screen
       \param s A string thats added to the begining of the line
       \param s2 A string thats added to the end of the line
    */
    void printEnergy2(std::string s, std::string s2);

    /*!
       \brief Print energy terms in single line to the screen
       \param os print stream
       \param s A string thats added to the begining of the line
       \param s2 A string thats added to the end of the line
    */
    void printEnergy2(std::ostream& os, std::string s, std::string s2);

    /*!
       \brief Print coords
    */
    void printCoords();

    /*!
       \brief Calculate bond energy and/or gradient
       \return bond energy
    */
    virtual double calcBondEnergy();

    /*!
       \brief Calculate angle energy and/or gradient
       \return angle energy
    */
    virtual double calcAngleEnergy();

    /*!
       \brief Calculate torsion energy and/or gradient
       \return torsion energy
    */
    virtual double calcTorsionEnergy();

    /*!
       \brief Calculate improper energy and/or gradient
       \return improper energy
    */
    virtual double calcImproperEnergy();

    /*!
       \brief Calculate non-bonded energy and/or gradient
       \return non-bonded energy
    */
    virtual double calcNonBondedEnergy();

    /*!
       \brief Calculate H-bond energy and/or gradient
       \return H-Bond energy
    */
    virtual double calcHBondEnergy();

    /*!
       \brief Get the total energy
       \return total energy
    */
    double getTotalEnergy();

    /*!
       \brief Get the bond energy
       \return bond energy
    */
    double getBondEnergy();

    /*!
       \brief Get the total van der Waals energy
       \return total vdW energy
    */
    double getTotalVDWEnergy();

    /*!
       \brief Get the total electrostatic energy
       \return ELE energy
    */
    double getTotalEleEnergy();

    /*!
       \brief Set gradient matrix to zero
    */
    void resetGMatrix();

    /*!
       \brief Set non-hydrogen atom gradients zero
    */
    void setNonHGradsToZero();

    /*!
       \brief Get the gMax value, gMaxAtom, and gNorm
       \param gMax gMax value
       \param gMaxAtom gMaxAtom index
       \param gNorm gNorm value
    */
    void getGradNorm(double& gMax, int& gMaxAtom, double& gNorm);

protected: // DATA

    //! Number of atoms
    int nAtoms;

    //! Number of distinct atom types
    int nDistinctTypes;

    //! Number of bonds containing Hydrogen
    int nBondsWithH;

    //! Number of bonds not containing hydrogen
    int nBondsWithOutH;

    //! Number of angles containing hydrogen
    int nAnglesWithH;

    //! Number of angles not containing hydrogen
    int nAnglesWithOutH;

    //! Number of dihedrals containing hydrogen
    int nDihedralsWithH;

    //! Number of dihedrals not containing hydrogen
    int nDihedralsWithOutH;

    /*!
        \brief Number of excluded atoms
        This is a sum over all atoms excluded atoms
    */
    int nExcluded;

    /*!
        \brief Number of excluded 1-4 atoms
        This is a sum over all atoms excluded atoms
    */
    int nExcluded14;

    //! Number of residues
    int nResidues;

    //! Number of unique bond types
    int nUniqueBonds;

    //! Number of unique angle types
    int nUniqueAngles;

    //! Number of unique angle types
    int nUniqueDihedrals;

    //! Number of unique angle types
    int nAtomsBigRes;

    //! Array of atom coordinates
    double *xyz;

    //! Atom symbols
    char  *symbols;

    //! Atom names
    char  *names;

    //! Array of atom integer types
    int   *iTypes;

    //! Array of atom integer types
    char   *cTypes;

    //! Array of atom charges
    double *charges;

    //! Array of atom charges
    double *masses;

    //! Array of atom coordinates
    double *gradients;

    //! Array of number of excluded atoms per atom
    int   *numExcluded;

    //! Array of number of excluded 1-4 atoms per atom
    int   *numExcluded14;

    //! Array of atom flags
    int   *atomFlags;

    //! Array of atom flags for energy decomposition
    int   *outFlags;

    //! Array of number of excluded atoms per atom
    int   *excluded;

    //! Array of number of excluded 1-4 atoms per atom
    int   *excluded14;

    //! L-J A Coefficient
    double *r12Params;

    //! L-J B Coefficient
    double *r6Params;

    //! Residue names
    char *resNames;

    //! Residue pointers
    int *resPointers;

    //! Non-Bonded index
    int *nonBondedParameterIndex;

    //! Number of bonds
    int nBonds;

    //! Array of bonds (A-B, atom indices)
    int *bonds;

    //! Bond parameters array
    double *bondParams;

    //! Number of angles
    int nAngles;

    //! Array of angles (A-B-C, atom indices)
    int *angles;

    //! Angle parameters array
    double *angleParams;

    //! Number of torsions
    int nTorsions;

    //! Array of torsions (A-B-C-D, atom indices)
    int *torsions;

    //! Torsion parameters array
    double *torsionParams;

    //! Number of Impropers
    int nImpropers;

    // Array of impropers
    int *impropers;

    //! Improper parameter array
    double *improperParams;

    //! Total number of non bonded interactions
    int nNonBonded;

    //! non-bonded atom indices array
    int *nonBonded;

    //! non-bonded array
    double *nonBondedParams;

    //! non-bonded pointers
    int *nonBondedPtrs;

    //! Total number of non bonded 1-4 interactions
    int nNonBonded14;

    //! non-bonded 1-4 atom indices array
    int *nonBonded14;

    //! non-bonded 1-4 array
    double *nonBonded14Params;

    //! non-bonded 1-4 pointers
    int *nonBonded14Ptrs;

protected: // DATA
    //! Calculate only E
    bool           bEnergy;

    //! Calculate both E and G
    bool           bGradient;

    //! Switch for bond part of the energy function
    bool           bBond;

    //! Switch for angle part of the energy function
    bool           bAngle;

    //! Switch for torsion part of the energy function
    bool           bTorsion;

    //! Switch for improper part of the energy function
    bool           bImproper;

    //! Switch for nonbonded part of the energy function
    bool           bNonBonded;

    //! Switch for H-bond part of the energy function
    bool           bHBond;

    //! Decompose bond energy
    bool           bBondDecompose;

    //! Decompose angle energy
    bool           bAngleDecompose;

    //! Decompose torsion energy
    bool           bTorsionDecompose;

    //! Decompose improper energy
    bool           bImproperDecompose;

    //! Decompose non-bonded energy
    bool           bNonBondedDecompose;

    //! Decompose H-Bonding energy
    bool           bHBondDecompose;

    //! Decomposed energy output file
    std::string    pwdFile;

    //! Output File Stream
    std::ofstream  pwdFileStream;

    //! Decomposed energy output file
    std::string    pwdGZFile;

#ifdef HAVE_ZLIB
    //! Output File Stream
    gzFile         pwdGZFileStream;
#endif
    //! Total energy:
    double         totalEnergy;

    //! bond energy
    double         bondEnergy;

    //! angle energy
    double         angleEnergy;

    //! torsion energy
    double         torsionEnergy;

    //! improper energy
    double         improperEnergy;

    //! van der Waals energy
    double         vdWEnergy;

    //! 1-4 part of the vdW energy
    double         vdW14Energy;

    //! R6 component of the vdW energy
    double         R6;

    //! R12 component of the vdW energy
    double         R12;

    //! electrostatic energy
    double         eleEnergy;

    //! 1-4 part of the ele energy
    double         ele14Energy;

    //! non bonded energy
    double         nonBondedEnergy;

    //! H-Bond energy
    double         hBondEnergy;
};
} // MTKpp namespace

#endif // MMPOTENTIAL_H
