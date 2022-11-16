/*!
   \file ring.h
   \brief Determines rings in a molecule
   \author Martin Peters

   $Date: 2010/03/29 20:45:26 $
   $Revision: 1.15 $

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

#ifndef RING_h
#define RING_h

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>

// #include "Utils/diagonalize.h"
#include <Eigen/Dense>
using namespace Eigen;

namespace MTKpp
{

class molecule;
class atom;
class Bond;
class vector3d;

// ============================================================
// Struct : ring
// ------------------------------------------------------------
/*! 
   \struct ring
   \brief Container for ring info
   \author Martin Peters
   \date 2005
*/
// ============================================================
struct ring
{
    //! internal index
    int index;

    //! Atoms in ring
    std::vector<atom*> atoms;

    //! Size of ring
    int size;

    /*!
       Is the ring planar?
       - 0 = No
       - 1 = Yes
    */
    int planar;

    /*!
       Is the ring aromatic?
       - 0 = No
       - 1 = Yes
    */
    int aromatic;

    /*!
       Is the ring a heterocycle?
       - 0 = No
       - 1 = yes
    */
    int hetero;

    //! Number of heteroatoms
    int nHetero;

    //! Number of Nitrogen atoms
    int nNitrogen;

    //! Number of Oxygen atoms
    int nOxygen;

    //! Number of Sulfur atoms
    int nSulfur;

    //! connection points
    std::vector<atom*> connPts;

    //! Centroid
    //ublas::vector<double> centroid;
    Eigen::VectorXd centroid;

    //! plane and normal

    //ublas::matrix<double, ublas::column_major> planeNormal;
    Eigen::MatrixXd planeNormal;

};

// ============================================================
// Class : rings()
// ------------------------------------------------------------
/*!
   \class rings

   \brief Determines rings in a molecule

   \author Martin Peters

   \date 2006

   See the \ref ringPerception page for more details.
*/

/*!
   \page ringPerception Ring Perception
    The algorithm used is in close agreement with that published by Fan, Panaye, Doucet, and Barbu in 1993.
    (<em>J. Chem. Inf. Comput. Sci.</em> <b>1993</b>, 33, 657-662)

    The functions contained in rings determines the smallest set of smallest rings (SSSR) from a molecule graph.

    The aromaticity of the rings is tested using the 4n+2 pi electron rule of Huckel.

   \section graphTheory_sec Graph Theory

     The algorithm relies on some knowledge of graph theory.

      - A RING can be represented as R(n1,n2,...) where ni are the atomic indices.

      - Each atom is a NODE.

      - The first selected NODE is called the ROOT.

      - A PATH walked that starts and ends at the ROOT NODE is called CLOSED.

      - CONNECTIVITY of a node is defined as the number of links to other nodes.

      - A NODE with a single link to another NODE is called TERMINAL.

      - A BLOCK is a group of NODES such that all links between them are involved in one 
      or more rings

    The SSSR of a molecule is represented as S(m1,m2,...) where mi are the ring sizes.

    An OPEN ACYCLIC NODE is an acyclic atom which is not located between two blocks.

    A CLOSED ACYCLIC NODE is an acyclic atom located between two blocks.

   \section algorithm_sec Graphical Illustration of the SSSR Algorithm

    - Take for example the following molecule:
    \image html ring1.png
    \image latex ring1.eps


    - Highlighting all OPEN ACYCLIC NODES
    \image html ring2.png
    \image latex ring2.eps

    - Remove all OPEN ACYCLIC NODES
    \image html ring3.png
    \image latex ring3.eps

    - Highlighting all CLOSED ACYCLIC NODES
    \image html ring4.png
    \image latex ring4.eps

    - Remove all CLOSED ACYCLIC NODES
    \image html ring5.png
    \image latex ring5.eps

    - Separate into blocks
    \image html ring6.png
    \image latex ring6.eps

    - How many rings are there in the current block?
    \image html ring7.png
    \image latex ring7.eps

    - Let NODE 1 be the ROOT NODE, we can produce numerous ring systems including R(1,2,3,15,13,14), 
      R(1,2,3,15,16,10,11,12,13,14), etc. How do we know we have found the smallest ring?
      The closed path found containing the ROOT NODE is recursively searched until it can not be reduced further, 
      in other words R(1,2,3,15,13,14) is found.
    \image html ring8.png
    \image latex ring8.eps

    - Once an <em>irreducible</em> closed path is found all NODES with two 
      links are removed. NODES 2, 1, and 14 can safely be removed.  The 
      algorithm then picks another ROOT NODES and the next ring is found 
      until all are found.
    \image html ring9.png
    \image latex ring9.eps

   \section huckel Huckel Rule of Aromaticity
     The algorithm used is in close agreement with that published by 
     Roos-Kozel, and Jorgensen in 1981.
     (<em>J. Chem. Inf. Comput. Sci.</em> <b>1981</b>, 21, 101-111)

     - Rings are classified as aromatic (AR), antiaromatic (AA), or
     nonaromatic (NA).  A ring system is aromatic if and only if it
     contains 4n+2 (n=0,1,2,3,4,...) pi electrons and is planar.

     - Nonaromatic are screened using the following rules:
      -# No intra-ring double bonds
      -# contains quaternary atom
      -# contains more than one saturated carbon
      -# contains a monoradical
      -# contains a sulfoxide or sulfone

     - The number of electrons atom contribute in simple aromatics (no 
       exocyclic pi bonds):
      -# 0 = cationic carbon and boron
      -# 2 = saturated heteroatoms
      -# 2 = anionic carbon
      -# 1 = radicals not on pi bond
      -# 1 = atoms on intra-ring pi bond

     - Rings thats contain exocyclic pi bonds are
 
    \image html aro.png
    \image latex aro.eps

   \section ringsCentroid_sec Ring Centroid
    \f[
        XYZ^{c} = {1 \over N} \sum_{i}^{N}XYZ_{i}
    \f]

   \section ringsPlaneNormal_sec Ring Normal and Plane

    - Form the coordinate matrix of the ring atoms
    \f[
        XYZ^{ring} = XYZ^{ring} - XYZ^{c}
    \f]

    - Principal Component Analysis of the Ring Coordinate Correlation Matrix
    \f[
        C = X^{T}X
    \f]

    \f[
        CP = \lambda_{i}P
    \f]

    \f[
        P = P + XYZ^{c}
    \f]

    - The first two eigenvectors define the plane of the ring, while the third defines the ring normal.

*/
// ============================================================
class rings
{
public:
    /*!
       \brief rings Constructor
       \param parent molecule pointer
    */
    rings(molecule *parent = 0);

    //! rings Destructor.
    virtual ~rings();

    //! Determine all rings in the molecule
    void determine();

    /*!
       \brief Determine if ring is aromatic
       \param r ring pointer
    */
    void kekulize(ring* r);

    /*!
       \brief Determine the plane and normal of the ring
       \param r ring pointer
    */
    int getPlaneNormal(ring* r);

    /*!
       \brief Determine the centroid of the ring
       \param r ring pointer
    */
    void calcCentroid(ring* r);

protected: // functions

    //! Removes hydrogen atoms from the cyclicAtom vector
    void removeHydrogenAtoms();

    //! Removes hydrogen containing bonds from the cyclicbonds vector
    void removeHydrogenBonds();

    //! Removes all OPEN ACYCLIC NODES (hydrogen and terminal atoms)
    void removeOpenAcyclic();

    //! Removes all OPEN ACYCLIC bonds
    void removeOpenAcyclicBonds();

    //! Find all neighbouring atoms
    void findNeighbors();

    /*! 
       \brief Find all neighbouring atoms
       \param a vector of atoms
       \param b vector of vector of ints
    */
    void findNeighbors(std::vector<atom*> a, std::vector< std::vector<int> > &b);

    /*! 
       \brief Find all bonds between atoms in a from all possible in b
       \param a vector of atoms
       \param b vector of bonds
       \param c vector of bonds
       \param d vector of int
    */
    void getBonds(std::vector<atom*> a, std::vector<Bond*> b, std::vector<Bond*> &c, std::vector<int> &d);

    //! Removes all CLOSED ACYCLIC NODES (atom joining ring systems)
    void removeClosedAcyclic();

    /*!
      \brief Pick Root Atom
      \param a vector of atoms
      \return index of atom with greatest connectivity
    */
    int pickRootAtom(std::vector<atom*> a);

    /*!
       \brief Separate cyclic atoms into blocks
    */
    void separateBlocks();

    /*!
       \brief Find the smallest set of smallest rings (SSSR)
    */
    void findSSSR();

    void decomposeBlock(std::vector<atom*> block);
    void getIrreducibleClosedPath(std::vector<atom*> block, int root, std::vector<atom*>& path);
    int numberRingsInBlock(std::vector<atom*> block);
    void getNewPath(std::vector<atom*> &path, std::vector<atom*> curPath, bool &done);
    void eliminateReducibleAtoms(std::vector<atom*> block, std::vector<atom*> path);
    void checkPath(std::vector<atom*> &path, bool &ok);

    /*!
       \brief Depth-First-Search of the molecular graph
       \param curAtom search NODE
       \param rootAtom ROOT NODE
       \param vertexes vector of atom pointers or vertexes
       \param edges vector of Bond pointers or edges
       \param neighbours NODES linking each other
       \param vertexesColor vertex colors used to traverse the molecular graph
       \param edgesColor edge colors used to traverse the molecular graph

       \param loop delete
       \param t delete

       \param curPath storage of the current closed cycle
       \param paths storage for all the paths found
    */ 
    void dfs_visit(int curAtom, int rootAtom,
//                   std::vector<atom*> vertexes, std::vector<Bond*> edges,
//                   std::vector< std::vector<int> > neighbours,
                   std::vector<atom*> &vertexes, std::vector<Bond*> edges,
                   std::vector< std::vector<int> > &neighbours,
                   std::vector<int> &vertexesColor, std::vector<int> &edgesColor,
                   bool &loop,int &t,
                   std::vector<atom*> &curPath,
                   std::vector< std::vector<atom*> > &paths);

    void dfs_visitNEW(int curAtom, int rootAtom, bool first, bool &loop,
//                      std::vector<atom*> vertexes, std::vector<Bond*> edges,
//                      std::vector< std::vector<int> > neighbors,
                      std::vector<atom*> &vertexes, std::vector<Bond*> edges,
                      std::vector< std::vector<int> > &neighbors,
                      std::vector<int> &vertexesColor, std::vector<int> &edgesColor,
                      std::vector<atom*> &path);

protected: // data
    //! molecule pointer
    molecule*                pParent;

    //! molecule atom list
    std::vector<atom*>       molAtomList;

    //! molecule Bond map
    std::map<int, Bond*>     molBondMap;

    //! Bond map iterator
    typedef std::map<int, Bond*>::iterator BondMapIterator;

    //! Hydrogen atoms
    std::vector<atom*>       hydrogens;

    //! Open and closed acyclic atoms
    std::vector<atom*>       acyclicAtoms;

    //! cyclic atoms
    std::vector<atom*>       cyclicAtoms;

    //! Cyclic Bonds
    std::vector<Bond*>       cyclicBonds;

    //! Neighboring atoms
    std::vector< std::vector<int> >    neighbors;

    //! vector of blocks
    std::vector< std::vector<atom*> >  blocks;

    //! atom iterator
    std::vector<atom*>::iterator atomIterator;

    // Bond iterator
    std::vector<Bond*>::iterator bondIterator;

    //! atom pointer
    atom*                    pAtom;

    //! atom pointer
    atom*                    pAtom2;

    //! Bond pointer
    Bond*                    pBond;

    std::vector<atom*>       finalPath;
};

} // MTKpp namespace

#endif // RING_h
