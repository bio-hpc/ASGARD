/*!
   \file functionalize.h
   \brief Breaks a molecule into its functional group
   \author Martin Peters

   $Date: 2010/03/29 20:43:22 $
   $Revision: 1.8 $

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

#ifndef FUNCTIONALIZE_H
#define FUNCTIONALIZE_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

#include "Utils/constants.h"

namespace MTKpp
{

class collection;
class molecule;
class submolecule;
class atom;
struct Bond;

class vector3d;

class stdLibrary;
class stdGroup;
class stdFrag;
struct stdAtom;
struct stdBond;

class fingerPrint;

// ============================================================
// Struct : funcGroup
// ------------------------------------------------------------
/*!
   \struct funcGroup
   \brief Container for functional group info
   \author Martin Peters
   \date 2006
*/
// ============================================================
struct funcGroup
{
    //! Atoms in function group
    std::map<stdAtom*, atom*>atomMap;

    //! standard fragment pointer
    stdFrag*                 pStdFrag;
};

// ============================================================
// Class : functionalize()
// ------------------------------------------------------------
/*! 
   \class functionalize

   \brief Determines functional groups in molecules using predefined fragments

   \author Martin Peters

   \date 2006

   See the \ref functionalGroupRec page for more details.
*/

/*!
   \page functionalGroupRec Functional Group Recognition

    The algorithm used is in close agreement with that published by Ullman,
    (<em>J. ACM</em> <b>1976</b>, 23, 31-42) and Willett, Wilson, and Reddaway,
    (<em>J. Chem. Inf. Comput. Sci.</em> <b>1991</b>, 31, 225-233).

    To functionalize a molecule involves searching it for chemical substructures.

    Substructures searching is known as the subgraph isomorphism problem of graph theory
    and belongs to the class of NP-complete computational problems.

    Due to the NP-complete nature of substructure searching usually a screen is carried out
    to eliminate subgraphs that cannot be contained in the molecule.

   \section functionalizeMotivations My Motivations to Write this Code:
     - To carry out functional group alignment of drug-like molecules.
     - To develope a tool for de novo design.
     - To optimize fragment positions in drug-protein complexes.
     - An analysis tool for a variety of project, for example, NMR.

    \subsection bruteForceSubgrapgh Brute Force Subgraph Isomorphism
     Generating the adjacency matrices A and B for the fragment and the molecule containing \f$P_A\f$ and \f$P_B\f$ 
     atoms respectively.

     An exhautive search involves generating \f$ P_B!/[P_A!(P_B-P_A)!] \f$ combinations of \f$P_A\f$.

     Ullmann first noticed that using a depth-first backtracking search dramatically increases efficiency.

     Using a labeled graph and a non-binary connection table increases algorithm speed.

   \subsection ullmanAlgorithmSec Ullman Subgraph Isomorghism

    The following example was taken from Molecular Modelling, Principles and Applications 2nd Edition by Andrew R. Leach.

    - Take for example the following fragment:
    \image html frag_functionalize.png
    \image latex frag_functionalize.eps

    - Determine the fragments' adjacency matrix
    \f[F =
     \left(\begin{array}{cccc}
      0 & 1 & 0 & 0 \\
      1 & 0 & 1 & 0 \\
      0 & 1 & 0 & 1 \\
      0 & 0 & 1 & 0
     \end{array}\right)
    \f]

    - Take for example the following molecule:
    \image html mol_functionalize.png
    \image latex mol_functionalize.eps

    - Determine the molecules' adjacency matrix
    \f[M =
     \left(\begin{array}{cccccc}
      0 & 1 & 0 & 0 & 0 & 0 \\
      1 & 0 & 1 & 1 & 0 & 0 \\
      0 & 1 & 0 & 0 & 0 & 0 \\
      0 & 1 & 0 & 0 & 0 & 0 \\
      0 & 0 & 0 & 1 & 0 & 0 \\
      0 & 0 & 0 & 1 & 0 & 0 \\
     \end{array}\right)
    \f]

    - The Ullman algorithm tries to find the matrix A which satisfy \f$A(AM)^T\f$
    \image html super_functionalize.png
    \image latex super_functionalize.eps

    \f[A =
     \left(\begin{array}{cccccc}
      0 & 0 & 1 & 0 & 0 & 0 \\
      0 & 1 & 0 & 0 & 0 & 0 \\
      0 & 0 & 0 & 1 & 0 & 0 \\
      0 & 0 & 0 & 0 & 0 & 1 \\
     \end{array}\right)
    \f]


    \f[
     A(AM)^T = 
     \left(\begin{array}{cccccc}
      0 & 0 & 1 & 0 & 0 & 0 \\
      0 & 1 & 0 & 0 & 0 & 0 \\
      0 & 0 & 0 & 1 & 0 & 0 \\
      0 & 0 & 0 & 0 & 0 & 1 \\
     \end{array}\right)
     \left(\begin{array}{cccc}
      0 & 1 & 0 & 0 \\
      1 & 0 & 1 & 0 \\
      0 & 1 & 0 & 0 \\
      0 & 1 & 0 & 1 \\
      0 & 0 & 1 & 0 \\
      0 & 0 & 1 & 0
     \end{array}\right) =
     \left(\begin{array}{cccc}
      0 & 1 & 0 & 0 \\
      1 & 0 & 1 & 0 \\
      0 & 1 & 0 & 1 \\
      0 & 0 & 1 & 0
     \end{array}\right) = F
    \f]



    This depth-first backtracking algorithm uses a General match matrix, M that contains all the possible 
    equivalences between atoms from A and B.  The elements of this matrix, 
    \f$m_{ij}(1 \leq i \leq Pa; 1 \leq j \leq Pb) \f$ are such that:
    \f[
      \def\tempa{if the $i$th atom of $A$ can be mapped to the $j$th atom of $B$}
      \def\tempb{otherwise}
       m_{ij} = 
      \left\{\begin{array}{lp}
       1 & \tempa \\
       0 & \tempb
      \end{array}\right
    \f]

    "The Ullmann heuristic states that if a fragment atom \f$a_i\f$ has a neighbor \f$x\f$, and a molecule atom
    \f$b_j\f$ can be mapped to \f$a_i\f$ , then there must exist a neighbor of \f$b_j\f$, \f$y\f$, that can be 
    mapped to \f$x\f$":
    \f[
       m_{ij} = \forall \quad x(1 ...P_A)[(a_{ix} = 1) \quad \Rightarrow \quad \exists \quad y(1 ... P_B)(m_{xy}b_{jy} = 1)]
    \f]

    if at any state during the search an atom \f$i\f$ in \f$A\f$ such that \f$m_{ij} = 0\f$ for all atoms in \f$B\f$ 
    then a mismatch is identified.
    \f[
       mismatch = \exists \quad i(1 ...P_A)[(m_{ij} = 0 \quad \forall \quad j(1 ... P_B)]
    \f]

    \code
     int   Pa;        // Number of atoms in the query pattern
     int   Pb;        // Number of atoms in the database structure
     array A[Pa,Pa];  // adjacency matrix of A
     array B[Pb,Pb];  // adjacency matrix of B
     array M[Pa, Pb]; // general match matrix
     bool isomorphism = false;

     void ullman(int d, array M) {
       array M1 = M;
       bool mismatch;

       for (all unique mappings of atom d) {
         choose new unique mapping for query atom d

         update M accordingly

         refine(M, mismatch);

         if (!mismatch) {
           if (d == Pa) {
             isomorphism = true;
             store M;
           }
           else {
             ullman(d+1, M);
           }
         }
         else {
           M = M1;
         }
       }
     }

     void refine(array M, bool mismatch) {
       mismatch = false;
       bool change = false;
       while (!change || mismatch) {
         for (int i = 0; i < Pa; i++) {
           for (int j = 0; j < Pb; j++) {
             if (M[i][j]) {
               assign mij;
               change = change or (!mij);
             }
           }
         }
         assign mismatch;
       }
     }
     int main() {
       ullman(1, M);
     }
    \endcode

   \section functionalizeAlgorithmSec The Functionalize Algorithm

    - Read all fragments
     -# Generate simple fingerprints

    - Read in molecule to functionalize
     -# Determine rings, add hydrogens, etc.
     -# Generate simple fingerprint
     -# Loop over all fragments in library
     -# Compare molecule to simple fingerprint (Screening)
     -# Match fragment to molecule using the Ullman and Willett algorithm of subgraph isomorphism
     -# if match assign fragment code to molecule

*/
// ============================================================
class functionalize
{
    friend class superimpose;
public:
    /*!
       \brief functionalize Constructor
       \param parent molecule pointer
    */
    functionalize(molecule *parent = 0);

    //! functionalize Destructor
    virtual ~functionalize();

    /*!
       \brief Functionalize
       \param stdLib standard library pointer
       \return success
    */
    int run(stdLibrary* stdLib);

    /*!
       \brief Determine Symmetric Fragments
       \param stdLib standard library pointer
       \return success
    */
    int determineSymmetricFrags(stdLibrary* stdLib);

protected: // FUNCTIONS

    /*!
       \brief Initialize the match matrix between the molecule and the fragment
       \param genMatchMatrix Match matrix
       \param molAtoms Number of atoms in the molecule
       \param molAdjMatSize Size of molecule adjacency matrix
       \param molAdjMatrix Molecule adjacency matrix
       \param molAtomSymbols Atom symbols in molecule
       \param molAtomKinds Atom types in molecule
       \param fragAtoms Number of atoms in the fragment
       \param fragAdjMatSize Size of fragment adjacency matrix
       \param fragAdjMatrix Fragment adjacency matrix
       \param fragAtomSymbols Atom symbols in fragment
       \param fragAtomKinds Atom types in fragment
    */
    void initializeMatchMatrix(int genMatchMatrix[],
                               int molAtoms, int molAdjMatSize,
                               int molAdjMatrix[], char molAtomSymbols[],
                               int molAtomKinds[],
                               int fragAtoms, int fragAdjMatSize,
                               int fragAdjMatrix[], char fragAtomSymbols[],
                               int fragAtomKinds[]);

    /*!
       \brief Determine if the fragment is a subgraph of the molecule
       \param pos Current fragment atom being considered
       \param fragAtoms Number of atoms in the fragment
       \param fragAdjMatrix Fragment adjacency matrix
       \param fragAtomSymbols Fragment atom element symbols
       \param molAtoms Number of atoms in the molecule
       \param molAdjMatrix Molecule adjacency matrix
       \param molAtomSymbols Molecule atom element symbols
       \param matchMatrix Match matrix
       \param subGraph storage vector for subgraph information
       \param subGraphs storage vector for all subgraphs
    */
    void ullmann(int pos, int fragAtoms, int fragAdjMatrix[], char fragAtomSymbols[],
                 int molAtoms, int molAdjMatrix[], char molAtomSymbols[], int matchMatrix[],
                 std::vector<int> &subGraph, std::vector<std::vector<int> > &subGraphs);

    /*!
       \brief Determine if the fragment is a subgraph of the molecule
       \param posMolAtom Which molecule atom is being considered
       \param posFragAtom Which fragment atom is being considered
       \param fragAtoms Number of fragment atoms
       \param molAtoms Number of molecule atoms
       \param matchMatrix Match matrix
    */
    void updateMatchMatrix(int posMolAtom, int posFragAtom, int fragAtoms,
                           int molAtoms, int matchMatrix[]);

    /*!
       \param fragAtoms Number of fragment atoms
       \param fragAdjMatrix Fragment adjacency matrix
       \param fragAtomSymbols Fragment atom element symbols
       \param molAtoms Number of molecule atoms
       \param molAdjMatrix Molecule adjacency matrix
       \param molAtomSymbols Molecule atom element symbols
       \param matchMatrix Match matrix
       \param mismatch Whether or not a mismatch occurs or not in the updated match match, i.e a row that contains all zeros
    */
    void refineMatchMatrix(int fragAtoms, int fragAdjMatrix[], char fragAtomSymbols[], int molAtoms,
                      int molAdjMatrix[], char molAtomSymbols[], int matchMatrix[], bool &mismatch);

protected: // DATA

    //---------------//
    // -  POINTERS - //
    //---------------//

    //! molecule pointer
    molecule*                pParent;

    //! atom pointer
    atom*                    pAtom;

    //! Bond pointer
    Bond*                    pBond;

    //! stdGroup pointer
    stdGroup*                pStdGroup;

    //! stdFrag pointer
    stdFrag*                 pStdFrag;

    //! stdAtom pointer
    stdAtom*                 pStdAtom;

    //! stdBond pointer
    stdBond*                 pStdBond;

    //! fingerPrint pointer
    fingerPrint*             pFingerPrint;

    //-----------------------//
    // - VECTOR CONTAINERS - //
    //-----------------------//

    //! stdGroup list
    std::vector<stdGroup*>   stdGroupList;

    //! stdFrag list
    std::vector<stdFrag*>    stdFragList;

    //----------------------//
    // - VECTOR ITERATORS - //
    //----------------------//

    //! stdGroup iterator
    typedef std::vector<stdGroup*>::iterator stdGroupIterator;

    //! stdFrag iterator
    typedef std::vector<stdFrag*>::iterator stdFragIterator;

    //! int iterator
    typedef std::vector<int>::iterator intIterator;
};

} // MTKpp namespace

#endif // FUNCTIONALIZE_H

