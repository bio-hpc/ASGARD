/*!
   \file hybridize.h
   \brief Determines hybridizations of atoms in a molecule
   \author Martin Peters

   $Date: 2010/03/29 20:43:22 $
   $Revision: 1.10 $

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

#ifndef HYBRIDIZE_H
#define HYBRIDIZE_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
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
namespace MTKpp
{

class collection;
class molecule;
class submolecule;
class atom;
struct Bond;
class vector3d;

// ============================================================
// Class : hybridize()
// ------------------------------------------------------------
/*!
   \class hybridize

   \brief Determines hybridizations of atoms in a molecule

   See the \ref bondTypeDetermination page for more details.

*/

/*!
   \page bondTypeDetermination Bond Type Determination
    There are two algorithms implemented to preceive atom hybridizations, bond order and formal charge

    The first is the Meng algorithm published in 1991. (<em>J. Comp. Chem.</em> <b>1991</b>, 12, 891-898)

    The second and default algorithm used is in close agreement with that published by Labute in 2005.
    (<em>J. Chem. Inf. Comput. Sci.</em> <b>2005</b>, 45, 215-221)

    \section LabuteAlgorithm The Labute Algorithm

     \image html hybridizeStart.png "PPARg Inhibitor"
     \image latex hybridizeStart.eps "PPARg Inhibitor" width=8cm

     Let \f$x_1, ..., x_n \f$ denote the 3D coordinates of n atoms with atomic 
     number \f$ Z_1, ..., Z_n \f$, number of bonded atoms, 
     \f$ Q_i \f$ and \f$r_{ij} = |x_i - x_j| \f$.

     \subsection LabuteStep1 Step 1
      Bonds are perceived by first producing a candidate list and then refining it using geometry.
      Covalent radii from Meng are used in the following: \f$ 0.1 < r_{ij} < R_i + R_j + 0.4 \f$

      For each atom, i, a "dimension", d_i, is assigned based on a principal
      component analysis of the Gram Matrix (the current atom and its bonded
      atoms).

      \f[
         D = \sum_{i=0}^{k}(q_i-\bar{q})(q_i-\bar{q})^{T}
      \f]

      \f[
         \bar{q} = {1 \over k}\sum_{i=0}^{k}q_i
      \f]

      \f$ d_i \f$ is set to k if k < 2 otherwise, \f$ d_i \f$ is the number of
      positive eigenvalues of \f$ D \f$ with square root greater than 0.2.
      \f$ d_i \f$ will be 0 for isolated atoms, 1 for terminal and linear
      atoms with at least 2 bonds, 2 for planar atoms (e.g., sp2 or square
      planar), and 3 otherwise (e.g. tetrahedral or sp3d)

     \image html hybridizeStep1.png "Step 1: D values"
     \image latex hybridizeStep1.eps "Step 1: D values" width=8cm

     \subsection LabuteStep2 Step 2
     An upper bound \f$ B_i \f$ for the number of bonds allowed by an atom is
     determined using \f$ d_i \f$ and \f$ Z_i \f$ as follows:
      - 1. \f$ B_i = 0 \f$ if \f$ d_i = 0 \f$
      - 2. \f$ B_i = 1 \f$ if \f$ Z_i < 3 (H, He) \f$
      - 3. \f$ B_i = 2 \f$ if \f$ d_i = 1, Z_i > 2 \f$ (sp hybridized and linear)
      - 4. \f$ B_i = 3 \f$ if \f$ d_i = 2, Z_i < 11 \f$ (sp2 hybridized for 2nd row elements)
      - 5. \f$ B_i = 4 \f$ if \f$ d_i = 2, Z_i > 10 \f$ or \f$ d_i = 3, Z_i < 11 \f$ (square planar or sp3 hybridized)
      - 6. \f$ B_i = 7 \f$ otherwise

     Only the shortest \f$ B_i \f$ are retained. At this point all atom hybridizations and bond orders are set to zero or undefined.

     \subsection LabuteStep3 Step 3
      This step assigns obvious hybridization based on \f$ d, Z, and Q \f$. Each substep is applied to atom of zero hybridization:
      - 1. \f$ Z_i = 1,2 \f$ atoms are set to sp3
      - 2. \f$ Q_i > 4, Z_i = (Group 5) \f$ and \f$ Q_i = 5, Z_i = {Group 4,5,6,7,8} \f$  atoms are set to sp3d
      - 3. \f$ Q_i > 4, Z_i = (Group 6) \f$ and \f$ Q_i = 6, Z_i = {Group 4,5,6,7,8} \f$  atoms are set to sp3d2
      - 4. \f$ Q_i > 4, Z_i = (Group 7) \f$ and \f$ Q_i = 7, Z_i = {Group 4,5,6,7,8} \f$  atoms are set to sp3d3
      - 5. \f$ Q_i = 4, Z_i > 10, d_i = 2 \f$ atoms are set to sp3d2 (square planar)
      - 6. \f$ Z_i = (Transition Metal) \f$ atoms are set to sp3d2
      - 7. \f$ Q_i > 4, Z_i > 10 \f$ and not {Si, P, S, Se} atoms are set to sp3d2, sp3 otherwise
      - 8. \f$ (Q_i = 4) and (Q_i = 3, d_i = 3)  \f$ atoms are set to sp3 (tetrahedral)
      - 9. \f$ (Q_i > 2, Z_i = (Group 6,7,8)) \f$ atoms are set to sp3
      - 10. \f$ Z_i not (C,N,O,Si,P,S,Se) \f$ atoms are set to sp3
      - 11. All atoms that none of its bonded atoms have zero hybridization are set to sp3

     \image html hybridizeStep3.png "Step 3: Hybridization Values"
     \image latex hybridizeStep3.eps "Step 3: Hybridization Values" width=8cm

     \subsection LabuteStep4 Step 4
      Only atoms with unassigned hybridizations have \f$ d < 3, Z = (C,N,O,Si,P,S,Se),Q < 4 \f$ and at least one bonded neighbor
      with an unassigned hybridization.
      At this stage all bond orders, \f$ b_ij \f$ in which atom \f$ i \f$ or \f$ j \f$ has non-zero hybridization are set to 1.

     \subsection LabuteStep5 Step 5
      A dihedral test is used to identify bonds of order 1.  The smallest out-of-plane dihedral is computed using:
      \f$ \min_{j,k} {|P_{ijkl}| , |\pi-P_{ijkl}|, |-\pi-P_{ijkl}|} \f$
      If this dihedral is greater than 15 degrees then \f$ b_{ij} \f$ is set to 1.

     \sa torsion()

     \image html hybridizeStep5.png "Step 5: Lines in red represent bonds of order 1"
     \image latex hybridizeStep5.eps "Step 5: Lines in red represent bonds of order 1" width=8cm

     \subsection LabuteStep6 Step 6
      The following table of lower bound single bond lengths:

    \f[
      \begin{array}{lc||lc}
      bond  &  dist &  bond &  dist \\
      C-C & 1.54 & C-N & 1.47 \\
      C-O & 1.43 &  C-Si & 1.86 \\
      C-P & 1.85 &  C-S & 1.75 \\
      C-Se & 1.97 & N-N & 1.45 \\
      N-O & 1.43 &  N-Si & 1.75 \\
      N-P & 1.68 &  N-S & 1.76 \\
      N-Se & 1.85 & O-O & 1.47 \\
      O-Si & 1.63 & O-P & 1.57 \\
      O-S & 1.57 &  O-Se & 1.97 \\
      Si-Si & 2.36 & Si-P & 2.26 \\
      Si-S & 2.15 & Si-Se & 2.42 \\
      P-P & 2.26 &  P-S & 2.07 \\
      P-Se & 2.27 & S-S & 2.05 \\
      S-Se & 2.19 & Se-Se & 2.34 \\
     \end{array}
    \f]

      and \f$ |x_i - x_j| > L_{ij} - 0.05 \f$, where \f$ L_{ij} \f$ is the reference bond length, is used to identify single bonds.

     \image html hybridizeStep6.png "Step 6: Lines in red represent bonds of order 1"
     \image latex hybridizeStep6.eps "Step 6: Lines in red represent bonds of order 1" width=8cm

     \subsection LabuteStep7 Step 7
      After steps 5 and 6 the hybridizations of all uncharacterized atoms not involved in a bond of unknown order are set to sp3.

     \image html hybridizeStep7.png "Step 7: Atom Hybridizations"
     \image latex hybridizeStep7.eps "Step 7: Atom Hybridizations" width=8cm

     \subsection LabuteStep8 Step 8
      A molecular graph is formed including only atoms (vertices) that have undefined hybridization and bonds (edges) that
      have unknown order.  This graph is then divided into components or subgraphs.  Each subgraph is analyzed independently and
      bond orders are assigned.

     \image html hybridizeStep8a.png "Step 8 A: Subgraphs of the Molecular Graph"
     \image latex hybridizeStep8a.eps "Step 8 A: Subgraphs of the Molecular Graph" width=8cm

     Assign edge weighs using \f$ w_{ij} = u_i + u_j + 2\delta(r_{ij} < L_{ij} - 0.11) + \delta(r_{ij} < L_{ij} - 0.25) \f$
     and the following atom parameters, \f$ u \f$ (3rd and 4th row elements are mapped to the corresponding 2nd row
     with 0.1 been subtracted, -20.0 for all other atoms):

    \f[
      \begin{array}{lccc}
      atom  &  Q=1 &  Q=2 &  Q=3 \\
      C-O   &  1.3 &  4.0 &  4.0 \\
      C-N   & -6.9 &  4.0 &  4.0 \\
      C     &  0.0 &  4.0 &  4.0 \\
      N-C-O & -2.4 & -0.8 & -7.0 \\
      N-C-N & -1.4 &  1.3 & -3.0 \\
      N     &  1.2 &  1.2 &  0.0 \\
      O-C-O & 4.2 & -8.1 & -20.0 \\
      O-C-O & 4.2 & -8.1 & -20.0 \\
      O     & 0.2 & -6.5 & -20.0 \\
     \end{array}
    \f]

     \image html hybridizeStep8b.png "Step 8 B: Edge Weights "
     \image latex hybridizeStep8b.eps "Step 8 B: Edge Weights " width=8cm

     A Maximum Weighted Matching Algorithm is employed to find the best arrangement of double/triple bonds in each subgraph.

     \image html hybridizeStep8c.png "Step 8 C: Maximum Weighted Matching Result"
     \image latex hybridizeStep8c.eps "Step 8 C: Maximum Weighted Matching Result" width=8cm

     \subsection LabuteStep9 Step 9
      Ionization states and formal charges are perceived from the connectivity and bond order.

      The formal charge of atom i, f_i, is calculated as follows:
      f_i = c_i - o_i + b_i
      where:
        c_i = atom group in periodic table
        o_i = nominal octet (2 for hydrogen, 6 for boron, 8 for carbon and all other sp3 atoms in groups 5,6,7,8)
        b_i = sum of the atom bond orders

      If the system doesn't contain hydrogen atoms, then the following steps are applied to set the ionization state:
      - 1. \f$ Z_i = 1 \f$ atoms are set to 0
      - 2. \f$ Z_i = (transition metals), Q_i > 0 \f$ atoms are set to \f$ f_i \f$
      - 3. \f$ Z_i = 1 \f$ atoms are set to 0
      - 4. \f$ Z_i = 1 \f$ atoms are set to 0

     \subsection LabuteStep10 Step 10

*/
// ============================================================
class hybridize
{
public:
    /*!
       \brief hybridize Constructor
       \param parent molecule pointer
       \param paramFile parameter file
    */
    hybridize(molecule *parent, std::string paramFile);

    //! hybridize Destructor
    virtual ~hybridize();

    /*!
       \brief Hybridize using bond types present in the file such as mol, mol2

       \return success
    */
    int run();

    /*!
       \brief Hybridize using the Meng algorithm

        J. Comp. Chem. 1991, 12, 891-898

       \return success
    */
    int runMeng();

    /*!
       \brief Hybridize using the Labute algorithm

       J. Chem. Inf. Model. 2005, 45, 215-221

       \return success
    */
    int runLabute();

protected: // FUNCTIONS

    //! Read in Labute parameters
    int readParameters();

    //! split string
    void lSplitString(std::string& s, const std::string s2, std::vector<std::string>& v, int i);

    //! Do step 1 of Labute algorithm
    int runLabute1();

    //! Do step 2 of Labute algorithm
    int runLabute2();

    //! Do step 3 of Labute algorithm
    int runLabute3();

    //! Center of mass
    void COM(std::vector<atom*> atoms, vector3d* c);

    //! set all atoms hybridizations
    void setHybridizations();

    //! set all atoms formal charges
    void setFormalCharges();

protected: // DATA

    //---------------//
    // -  POINTERS - //
    //---------------//

    //! molecule pointer
    molecule*                pParent;

    //! Parameter file
    std::string              parameterFile;

    //! smallest out-of-plane dihedral
    double                   smallOutOfPlaneDihedral;

    //! single bond map
    std::map<std::string, double> singleBonds;

    //! single bond map iterator
    typedef std::map<std::string, double>::iterator mapIterator;

    //! log likelihoods ratios
    std::map<std::string, std::vector<double> > logLikelihoodRatios;

    //! single bond map iterator
    typedef std::map<std::string, std::vector<double> >::iterator logLikelihoodRatiosIterator;

    //! Single bond edge weight
    double                   edgeWeightSingle;

    //! Double bond edge weight
    double                   edgeWeightDouble;

    //! Triple bond parameter
    double                   tripleBondParameter;

    //! atom pointer
    atom*                    pAtom;

    //! Bond pointer
    Bond*                    pBond;

    //! atom list
    std::vector<atom*> atomList;

    //! number of atoms in the molecule
    unsigned int nAtoms;

    //! a "dimension" for each atom
    std::vector<int> dim;

    //! upper bound on the number of bonds each atom can have
    std::vector<int> B;

    //! number of connections of each atom
    std::vector<int> Q;

    //! atomic number of each atom
    std::vector<int> atNumbers;

    //! group which each atom is a part of
    std::vector<int> atGroups;

    //! period which each atom is a part of
    std::vector<int> atPeriods;

    //! formal charge of each atom
    std::vector<int> formalCharges;

    //! atom symbols
    std::vector<std::string> atSymbols;

    //! atom electronegativities
    std::vector<double> atENs;

    //! Planarity flag of each atom
    std::vector<int> planarFlags;

    //! hybrizations of each atom
    std::vector<int> hybridizations;

    //! bonds
    std::map<int, Bond*> bonds;

    //! Bond map iterator
    typedef std::map<int, Bond*>::iterator BondMapIterator;

    //! bonding indices
    std::vector<std::vector<int> > bdAtoms;

    //! Error message
    std::string errorMessage;
};

} // MTKpp namespace

#endif // HYBRIDIZE_H

