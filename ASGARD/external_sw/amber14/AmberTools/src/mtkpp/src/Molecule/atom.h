/*!
   \file atom.h
   \brief Container for atom information
   \author Martin Peters

   $Date: 2010/03/29 20:42:27 $
   $Revision: 1.20 $

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

#ifndef ATOM_H
#define ATOM_H

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <map>

namespace MTKpp
{

class collection;
class molecule;
class submolecule;
struct element;
struct stdAtom;
class vector3d;

// ============================================================
// Class : atom()
// ------------------------------------------------------------
/*! 
   \class atom

   \brief Container for atom info

   \author Martin Peters

    \section AtomTypeDefinitions Atom Type (Internal) Definitions
      - 0 = Undefined
      - 1 = Hydrogen
      - 2 = Terminal Heavy Atom
      - 3 = Open Chain Heavy Atom
      - 4 = Closed Chain Heavy Atom
      - 5 = Ring Heavy Atom
      - 6 = Aromatic Ring Heavy Atom
      - 7 = Chain Atom (not terminal, Ring or Hydrogen)
      - 8 = Chain or terminal Atom

    \section AtomHybridizationDefinitions Atom Hybridization Definitions
      - 0 = undefined
      - 1 = s
      - 2 = sp
      - 3 = sp2
      - 4 = sp3
      - 5 = sp3d
      - 6 = sp3d2
      - 7 = sp3d3
      - 8 = other

    \section AtomHydrophilicDefinitions Atom Hydrophilic Definitions
      - -1  Undefined
      - 0   Hydrophobic
      - 1   N,O
      - 2   S in SH
      - 3   S with double bond
      - 4   <= 1 bond away from charged atom
      - 5   <= 1 bond away from OH, NH, or NH2 with no delocalized electrons
      - 6   <= 1 bond away from SH with no delocalized electrons
      - 7   <= 1 bond away from O with double bond
      - 8   <= 1 bond away from S with valence > 2
      - 9   > one neighboring O or N with no delocalized electrons
        -# A two bond away from O with double bond
        -# B one bond away from S with valence > 2
      - 10  two or more instance of any of the previous two(A,B)

*/
// ============================================================

class atom
{
public:

    /*!
       \brief copy Constructor
       \param pAt atom pointer
       \param parent submolecule pointer
    */
    atom(atom* pAt, submolecule *parent = 0);

    /*!
       \brief atom Constructor
       \param parent submolecule pointer
    */
    atom(submolecule *parent = 0);

    //! atom Destructor.
    virtual ~atom();

    ///////////////////////
    // - SET FUNCTIONS - //
    ///////////////////////

    /*!
       \brief Set element of atom
       \param ele element pointer
    */
    void                     setElement      (element* ele);

    /*!
       \brief Set name of atom
       \param name atom name
    */
    void                     setName         (const std::string name);

    /*!
       \brief Set coordinates of atom
       \param x X coordinate 
       \param y Y coordinate
       \param z Z coordinate
    */
    void                     setCoords       (const double& x,const double& y,const double& z);

    /*!
       \brief Set index of atom
       \param n index
    */
    void                     setIndex        (const int& n); // internal molecule indexing

    /*!
       \brief Set collection index of atom
       \param n index
    */
    void                     setColIndex     (const int& n); // internal collection indexing

    /*!
       \brief Set file id of atom
       \param n file id
    */
    void                     setFileID       (const int& n); // input file indexing

    //void                    setX            (const double&);
    //void                    setY            (const double&);
    //void                    setZ            (const double&);

    /*!
       \brief Set atomic number of atom
       \param n atomic number
    */
    void                     setAtomicNum    (const int& n);

    /*!
       \brief Set nuclear charge of atom
       \param n nuclear charge
    */
    void                     setZcharge      (const double& n);

    /*!
       \brief Set formal charge of atom
       \param n formal charge
    */
    void                     setFormalCharge (const int& n);

    /*!
       \brief Set valene of atom
       \param n atom valence
    */
    void                     setValence(const int& n);

    /*!
       \brief Set hybridization of atom
       \param i atom hybridization
    */
    void                     setHybridization(const int& i);

    /*!
       \brief Set type of atom
       \param n atom type
    */
    void                     setType(const int& n);

    /*!
       \brief Set Meng type
       \param n Meng atom type
    */
    void                     setMengType(const std::string& n);

    /*!
       \brief Set occupancy
       \param o occupancy
    */
    void                     setOccupancy(const double& o);

    /*!
       \brief Set temperature factor
       \param b temperature factor
    */
    void                     setTempFactor(const double& b);

    /*!
       \brief Set setHydrophilicity
       \param h hydrophilicity
    */
    void                     setHydrophilicity(const int& h);

    /*!
       \brief Get Meng type
       \return Meng atom type
    */
    std::string              getMengType();

    /*!
       \brief Add atom property
       \param name name of property
       \param value value of property
    */
     void                     addProperty(const std::string& name, double value);

    /*!
       \brief Add atom property
       \param name name of property
       \param value value of property
    */
     void                     addProperty(const std::string& name, int value);

    //-------------------//
    // - MAP ITERATORS - //
    //-------------------//
    //! property map iterator
    typedef std::map<std::string, double>::iterator   PropertyMapIterator;

    //! property map iterator
    typedef std::map<std::string, int>::iterator      intPropertyMapIterator;


    //--------------------//
    // - MAP CONTAINERS - //
    //--------------------//

    //! Atom properties map
    std::map<std::string, double>                     itsPropertiesMap;

    //! Atom properties map
    std::map<std::string, int>                        itsIntPropertiesMap;

    ///////////////////////
    // - GET FUNCTIONS - //
    ///////////////////////

    /*!
       \brief Get element of atom
       \return atom element
    */
    element*                 getElement();

    /*!
       \brief Get element symbol of atom
       \return atom element
    */
    std::string              getElementSymbol();

    /*!
       \brief Get name of atom
       \return atom name
    */
    std::string              getName();

    /*!
       \brief Get coordinates of atom
       \return atom coordinates 
    */
    vector3d*                getCoords();

    /*!
       \brief Get index of atom
       \return atom index 
    */
    int                      getIndex();

    /*!
       \brief Get collection index of atom
       \return atom index 
    */
    int                      getColIndex();

    /*!
       \brief  Get file id of atom
       \return atom file id 
    */
    int                      getFileID();

    /*!
       \brief Get x coordinate of atom
       \return atom x coordinate 
    */
    double                   getX();

    /*!
       \brief Get y coordinate of atom
       \return atom y coordinate 
    */
    double                   getY();

    /*!
       \brief Get z coordinate of atom
       \return atom z coordinate 
    */
    double                   getZ();

    /*!
       \brief Get parent of atom
       \return submolecule with which the atom is apart of
    */
    submolecule*             getParent();

    /*!
       \brief Get atomic number of atom
       \return atom atomic number
    */
    int                      getAtomicNum();

    /*!
       \brief Get atomic number of atom
       \return atom atomic number
    */
    double                   getAtomicMass();

    /*!
       \brief Get nuclear charge of atom
       \return atom nuclear charge
    */
    double                   getZcharge();

    /*!
       \brief Get formal charge of atom
       \return atom formal charge
    */
    int                      getFormalCharge();

    /*!
       \brief Get valene of atom
       \return atom valence
    */
    int                      getValence();

    /*!
       \brief Get hybridization of atom
       \return atom hybridization
    */
    int                      getHybridization();

    /*!
       \brief Get type of atom
       \return atom type
    */
    int                      getType();

    /*!
       \brief Get number of bonded oxygen atoms
       \return number of bonded oxygen atoms
    */
    int                      numBondedOxygens();

public: // - functions

    /*!
       \brief get atom property
       \param name name of property
       \return value of property
    */
    double                   getProperty(const std::string& name);

    /*!
       \brief get atom properties
       \return map of properties
    */
    std::map<std::string, double> getPropertyMap();

    /*!
       \brief get atom property
       \param name name of property
       \return value of property
    */
    int                      getIntProperty(const std::string& name);

    /*!
       \brief get atom properties
       \return map of properties
    */
    std::map<std::string, int> getIntPropertyMap();

    /*!
       \brief has atom a certain property
       \param name name of property
       \return got it or not
    */
    bool                     hasProperty(const std::string& name);

    /*!
       \brief Get occupancy
       \return occupancy
    */
    double                   getOccupancy();

    /*!
       \brief Get occupancy
       \return occupancy
    */
    double                   getTempFactor();

    /*!
       \brief Get setHydrophilicity
       \return hydrophilicity
    */
    int                      getHydrophilicity();

public: // - Molecular Mechanics functions

    /*!
       \brief Set standard atom
       \param a standard atom
    */
    void                     setStdAtom(stdAtom* a);

    /*!
       \brief Set standard atom
       \return standard atom
    */
    stdAtom*                 getStdAtom();

    /*!
       \brief Add bonded atom
       \param a bonded atom
    */
    void                     addBondedAtom(atom* a);

    /*!
       \brief Add 1-3 bonded atom
       \param a 1-3 bonded atom
    */
    bool                     addBonded13Atom(atom* a);

    /*!
       \brief Add 1-4 bonded atom
       \param a 1-4 bonded atom
    */
    void                     addBonded14Atom(atom* a);

    /*! 
       \brief Vector of bonded atoms
    */
    std::vector<atom*>       bondedAtoms;

    /*! 
       \brief Vector of 1-3 bonded atoms
    */
    std::vector<atom*>       bonded13Atoms;

    /*!
       Vector of 1-4 bonded atoms
    */
    std::vector<atom*>       bonded14Atoms;

    /*!
       \brief Get Bonded Atoms
       \return vector of bonded atom
    */
    std::vector<atom*>       getBondedAtoms();

    /*!
       \brief Get Bonded Heavy Atoms
       \return vector of bonded atom
    */
    std::vector<atom*>       getBondedHeavyAtoms();

    /*!
       \brief Get number of Bonds
       \return number of bonded atom
    */
    int                      getNumBonds();

    /*!
       \brief Get number of 1-3 Bonds
       \return number of bonded atom
    */
    int                      getNum13Bonds();

    /*!
       \brief Get number of 1-4 Bonds
       \return number of bonded atom
    */
    int                      getNum14Bonds();

    /*!
       \brief Check if a bond exist with anyother atom
       \param a bonded atom
    */
    bool                     hasBondedAtom(atom* a);

    /*!
       \brief Check if a 1-3 bond exist with anyother atom
       \param a 1-3 bonded atom
    */
    bool                     has13BondedAtom(atom* a);

    /*!
       \brief Get 1-3 Bonded Atoms
       \return vector of bonded atom
    */
    std::vector<atom*>       get13BondedAtoms();

    /*!
       \brief Check if a 1-4 bond exist with anyother atom
       \param a 1-4 bonded atom
    */
    bool                     has14BondedAtom(atom* a);

    /*!
       \brief Get 1-4 Bonded Atoms
       \return vector of bonded atom
    */
    std::vector<atom*>       get14BondedAtoms();

protected:

    //! atoms' submolecule
    submolecule*             pParent;

    //! atoms' element
    element*                 pElement;

    //! element string
    std::string              itsElement;

    //! name
    std::string              itsName;

    //! coordinates
    vector3d*                pCoords;

    //! collection index
    int                      itsCollectionIndex;

    //! index
    int                      itsIndex;

    //! file id
    int                      itsFileId;

    //! atomic number
    int                      itsAtomicNum;

    //! nuclear charge
    double                   itsZcharge;

    //! formal charge
    int                      itsFormalCharge;

    //! valence
    int                      itsValence;

    /*!
       atom hybridization definitions
       - 0 = undefined
       - 1 = s
       - 2 = sp
       - 3 = sp2
       - 4 = sp3
       - 5 = sp3d
       - 6 = sp3d2
       - 7 = sp3d3
       - 8 = other
    */
    int                      itsHybridization;
/*

used to be:
       - 0 = undefined
       - 1 = sp3
       - 2 = sp2
       - 3 = sp
       - 4 = s
       - 5 = other

http://chemistry.boisestate.edu/rbanks/inorganic/bonding%20and%20hybridization/bonding_hybridization.htm
http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/idatm.html
http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/vdwrad.html

linear	 				sp	 	BeF2, CO2
trigonal planar			sp2	 	BF3, CO32-
tetrahedral	 			sp3	 	CH4, SO42-
trigonal bipyramidal	sp3d	PF5
octahedral		 		sp3d2	F6

One thing that students should be careful about as far as hybridization is concerned is that 
elements in row 2 (carbon, nitrogen, oxygen, and fluorine) must NEVER exceed their octets. 
That is, they should never have more than four bonds. Elements in row 3 and lower do on occasion 
exceed their octets, and they can do this by bringing in d-orbitals when they hybridize. Thus, 
phosphorus (but never nitrogen) can use one of its 3d orbitals to make five sp3d hybrids (leading 
to five substituents and trigonal bipyramidal molecular geometry), and sulfur (but never oxygen) 
can use two 3d orbitals to make six sp3d2 hybrids (leading to six bonds and octahedral geometry).
*/

    /*!
        atom Type Definitions
        - 0 = Undefined
        - 1 = Hydrogen
        - 2 = Terminal Heavy Atom
        - 3 = Open Chain Heavy Atom
        - 4 = Closed Chain Heavy Atom
        - 5 = Ring Heavy Atom
        - 6 = Aromatic Ring Heavy Atom
        - 7 = Chain Atom (not terminal, Ring or Hydrogen)
        - 8 = Chain or terminal Atom
    */
    int                      itsType;

    /*
       MMFF94 Atom Types
         Name   Num   CN   FC    Description
        - CR      1   [4]        Alkyl Carbon
        - C=C     2   [3]        Vinylic C
        - CSP2    2   [3]        Generic sp2 C
        - C=0     3   [3]        Generic Carbonyl C [3]
        - C=N     3   [3]        Imine-type C   [3]
        - CGD     3   [3]        Guanidine C    [3]
        - C=OR    3   [3]        Ketone or Aldehyde Carbonyl C [3]
        - C=ON    3   [3]        Amdide Carbonyl C [3]
        - COO     3   [3]        Carboxylic acid or ester carbonyl carbon [3]
        - COON    3   [3]        Carbamate carbonyl C [3]
        - COOO    3   [3]        Carbonic acid or ester carbonyl C [3]
        - C=OS    3   [3]        Thioester carbonyl C double bonded to O [3]
        - C=S     3   [3]        Thioester carbonyl C double bonded to S [3]
        - C=SN    3   [3]        Thioamide C, double bonded to S [3]
        - CSO2    3   [3]        C in > C=SO2 [3]
        - CS=O    3   [3]        Sulfinyl C in > C=S=O [3]
        - CSS     3   [3]        Thiocarboxylic acid or ester C [3]
        - C=P     3   [3]        C doubly bonded to P [3]
        - CSP     4   [2]        Acetylenic C [2]
        - =C=     4   [2]        Allenic C [2]
        - HC      5   [1]        H attached to C [1]
        - HSI     5   [1]        H attached to Si [1]
        - -O-     6   [2]        Generic divalent O [2]
        - OR      6   [2]        Ether O [2]
        - OC=O    6   [2]        Carboxylic acid or ester O [2]
        - OC=C    6   [2]        Enolic or phenolic O [2]
        - OC=N    6   [2]        O in -O-C=N- moiety [2]
        - OC=S    6   [2]        divalent O in thioacid or ester [2]
        - ONO2    6   [2]        divalent nitrate (ether) O [2]
        - ON=O    6   [2]        divalent nitrate (ether) O [2]
        - OSO3    6              divalent O in sulfate
        - OSO2    6              divalent O in sulfite
        - OSO     6   [2]        one of two divalent O attached to S
        - OS=O    6   [2]        divalent O in R(RO)S=O
        - -OS     6   [2]        other divalent oxygen attached to S
        - OPO3    6   [2]        divalent O in phosphate group
        - OPO2    6   [2]        divalent O in phosphite group
        - OPO     6              divalent O, one of the 2 Os attached to P
        - -OP     6              Other divalent O attached to P
        - O=C     7   [1]        Generic carbonyl O
        - O=CN    7   [1]        Carbonyl O in amides
        - O=CR    7   [1]        Carbonyl O in aldehydes and ketones
        - O=CO    7   [1]        Carbonyl O in acids and esters
        - O=N     7   [1]        Nitroso O
        - O=S     7   [1]        Doubly bonded sulfoxide O
        - O=S=    7   [1]        O=S on S doubly bonded to e.g., C
        - NR      8   [3]        Amine N
        - N=C     9   [2]        Imine N
        - N=N     9   [2]        Azo group N
        - NC=O   10   [3]        Amide N
        - NC=S   10   [3]        Thioamide N
        - NN=C   10   [3]        N in N-N=C moiety with deloc. lp
        - NN=N   10   [3]        N in N-N=N moiety with deloc. lp
        - F      11   [1]        Fluorine
        - Cl     12   [1]        Chlorine
        - Br     13   [1]        Bromine
        - I      14   [1]        Iodine
        - S      15   [2]        Thiol, sulfide, or disulfide S
        - S=C    16              S doubly bonded to C
        - S=O    17   [3]        Sulfoxide S
        - >S=N   17              tricoordinate S doubly bonded to N
        - SO2    18   [4]        Sulfone S
        - SO2N   18   [4]        Sulfonamide S
        - SO3    18   [4]        Sulfonate group S
        - SO4    18   [4]        Sulfate group S
        - =SO2   18              Sulfone S, doubly bonded to C
        - SNO    18              Sulfur in nitrogen analog of a sulfone
        - SI     19   [4]        Silicon
        - CR4R   20   [4]        Aliphatic C in 4-mem ring
        - HOR    21   [1]        Hydroxyl H in alcohols
        - HO     21   [2]        Generic hydroxyl hydrogen
        - CR3R   22   [4]        Aliphatic carbon in 3-mem ring
        - HNR    23   [1]        Generic H on sp3 N, e.g. in amines
        - HPYL   23   [1]        H on N in pyrrole
        - H3N    23   [1]        H in ammonia
        - HNOX   23   [1]        H on N in a N-oxide
        - HOCO   24   [1]        Hydroxyl H in carboxylic acids
        - HOP    24   [1]        Hydroxyl H in H-O-P moiety
        - PO4    25   [4]        Phosphate group P
        - PO3    25   [4]        P with 3 attached Os
        - PO2    25   [4]        P with 2 attached Os
        - PO     25   [4]        Phosphine oxide P
        - PTET   25   [4]        General tetracoordinate phoshorus
        - P      26   [3]        Phosphorus in phosphines
        - HN=C   27   [1]        H on imine N
        - HN=N   27   [1]        H on azo N
        - HNCO   28   [1]        H on amide N
        - HNCS   28   [1]        H on thioamide N
        - HNCC   28   [1]        H on emamine N
        - HNCN   28   [1]        H in H-N-C=N moiety
        - HNNC   28   [1]        H in H-N-N=C moiety
        - HNNN   28   [1]        H in H-N-N=N moiety
        - HNSO   28   [1]        H on NSO, NSO2, or NSO3 N
        - HNC%   28   [1]        H on N triply bonded to C [1]
        - HSP2   28   [1]        Generic H on sp3 N [1]
        - HOCC   29   [1]        Enolic or phenolic hydroxyl H [1]
        - HOCN   29   [1]        Hydroxyl H in HO-C=N moiety [1]
        - CR4E   30   [3]        Olefinic C in 4-mem ring [3]
        - HOH    31              Hydroxyl H in water
        - O2CM   32   [1] {-1/2} O in carboxylate group 
        - ONX    32   [1]        O in N-Oxides
        - O=N    32   [1]        O in nitroso group
        - O2N    32   [1]        O in nitro group
        - O2NO   32   [1]        Nitro-group O in nitrate
        - O3N    32   [1] {-1/3} Nitrate anion O
        - O-S    32   [1]        Single terminal O on tetrahedral S
        - O2S    32   [1]  {var} One of 2 terminal O's on S
        - O3S    32   [1]  {var} One of 3 terminal O's on S
        - O4S    32   [1] {-1/2} Terminal O in sulfate anion
        - OSMS   32   [1] {-1/2} Terminal O in thiosulfinate anion
        - OP     32   [1]        Oxygen in phosphine oxide
        - O2P    32   [1]  {var} One of 2 terminal O's on P
        - O3P    32   [1]  {var} One of 3 terminal O's on P
        - O4P    32   [1]  {var} One of 4 terminal O's on P
        - O4Cl   32   [1] {-1/4} O in perchlorate anion
        - HOS    33   [1]        H on O attached to S
        - NR+    34   [4]    {1} Quaternary N
        - OM     35   [1]   {-1} Oxide O on sp3 C
        - OM2    35   [1]   {-1} Oxide O in sp2 C
        - HNR+   36   [1]        H on quaternary N
        - HIM+   36   [1]        H on imidazolium N
        - HPD+   36   [1]        H on pyridinium N
        - HNN+   36   [1]        H on amidinium N
        - HNC+   36   [1]        H on protonated imine N
        - HGD+   36   [1]        H on guanidinium N
        - CB     37   [3]        Aromatic C e.g. in benzene
        - NPYD   38   [2]        Aromatic N with sigma lp
        - NPYL   39   [2]        Aromatic 5-ring N with pi lp
        - NC=C   40   [3]        Enamine or aniline N, deloc. lp [3]
        - NC=N   40   [3]        N in N-C=N [3]
        - ////NC=N     40 N in N-C=P
        - NC%C   40   [3]        N attached to C-C triple bond
        - CO2M   41   [3]        C in carboxylate anion
        - CS2M   41   [3]        C in thiocarboxylate anion
        - NSP    42   [1]        Triply bonded N
        - NSO2   43   [3]        Sulfonamide N
        - NSO3   43   [3]        Sulfonamide N
        - NC%N   43   [3]        N attached to cyano group
        - STHI   44   [2]        Aromatic 5-ring sulfur with pi lone pair
        - NO2    45   [3]        N in nitro group
        - NO3    45   [3]        N in nitrate group
        - N=O    46   [2]        Nitrogen in nitroso group
        - NAZT   47   [1]        Terminal N in azido or diazo group
        - NSO    48   [2]        Divalent N replacing monovalent O in SO2 group
        - O+     49   [3]    {1} Oxonium O
        - HO+    50   [1]        H on oxonium O
        - O=+    51   [2]    {1} Oxenium O
        - HO=+   52   [1]        H on oxenium O
        - =N=    53   [2]        Central N in C=N=N or N=N=N
        - N+=C   54   [3]    {1} Iminium N
        - N+=N   54   [3]    {1} +ve N doubly bonded to N
        - NCN+   55   [3]  {1/2} Either N in N+=C-N
        - NGD+   56   [3]  {1/3} Guanidinium N
        - CGD+   57   [3]        Guanidinium C 
        - CNN+   57   [3]        C in +N=C-N: resonance structure
        - NPD+   58   [3]    {1} Aromatic N in pyridinium
        - OFUR   59   [2]        Aromatic 5-ring O with pi lone pair
        - C%-    60   [1]        Isonitrile C
        - NR%    61   [2]        Isonitrile N
        - NM     62   [2]   {-1} Anionic divalent N
        - C5A    63   [2]        Aromatic 5-ring C, alpha to N, O, or S
        - C5B    64   [2]        Aromatic 5-ring C, beta to N, O, or S
        - N5A    65   [2]        Aromatic 5-ring N, alpha to N, O, or S
        - N5B    66   [2]        Aromatic 5-ring N, beta to N, O, or S
        - N2OX   67   [3]        sp2-hybridized N-oxide N
        - N3OX   68   [4]        sp3-hybridized N-oxide N
        - NPOX   69   [3]        Pyridinium N-oxide N
        - OH2    70   [2]        O in water
        - HS     71   [1]        H attached to S
        - HS=N   71   [1]        H attached to >S= Sulfur double bonded to N
        - HP     71   [1]        H attached to P
        - S-P    72   [1]        Terminal S bonded to P
        - SM     72   [1]   {-1} Anionic terminal S
        - SSMO   72   [1] {-1/2} Terminal S in thiosulfinate group
        - SO2M   73   [3]        S in anionic sulfinate group
        - SSOM   73   [3]        Tricoordinate S in anionic thiosulfinate group
        - =S=O   74              Sulfinyl S e.g. in C=S=O
        - -P=C   75   [3]        P doubly bonded to C
        - N5M    76   [2]  {var} N in 5-ring aromatic anion
        - CLO4   77   [4]        Perchlorate anion chlorine
        - C5     78   [3]        General C in 5-mem heteroaromatic ring
        - N5     79   [2]        General N in 5-mem heteroaromatic ring
        - CIM+   80              Aromatic C between N's in imidazolium
        - NIM+   81   [3]  {1/2} Aromatic N in imidazolium
        - N5A+   81   [3]    {1} +ve N in 5-ring alpha position
        - N5B+   81   [3]    {1} +ve N in 5-ring beta position
        - N5+    81   [3]    {1} +ve N in other 5-ring position
        - N5AX   82   [3]        N-oxide N in 5-ring alpha position
        - N5BX   82   [3]        N-oxide N in 5-ring beta position
        - N5OX   82   [3]        N-oxide N in other 5-ring position
        - FE+2   87   [0]    {2} Dipositive iron cation
        - FE+3   88   [0]    {3} Tripositive iron cation
        - F-     89   [0]   {-1} Floride anion
        - CL-    90   [0]   {-1} Chloride anion
        - BR-    91   [0]   {-1} Bromide anion
        - LI+    92   [0]    {1} Li cation
        - NA+    93   [0]    {1} Na cation
        - K+     94   [0]    {1} Dipositive K
        - ZN+2   95   [0]    {2} Dipositive Zn
        - CA+2   96   [0]    {2} Dipositive Ca
        - CU+1   97   [0]    {1} Monopositive Cu
        - CU+2   98   [0]    {2} Dipositive Cu
        - MG+2   99   [0]    {2} Dipositive Mg

       GAFF Atom Types
         Name Hybridization  Description
        - c   sp2            Sp2 C carbonyl group
        - c1  sp             Sp C
        - c2  sp2            Sp2 C
        - c3  sp3            Sp3 C
        - ca  sp2            Sp2 C in pure aromatic systems
        - cp  sp2            Head Sp2 C that connect two rings in biphenyl sys
        - cq  sp2            Head Sp2 C that connect two rings in biphenyl sys. identical to cp
        - cc  sp2            Sp2 carbons in non-pure aromatic systems
        - cd  sp2            Sp2 carbons in non-pure aromatic systems, identical to cc
        - ce  sp2            Inner Sp2 carbons in conjugated systems
        - cf  sp2            Inner Sp2 carbons in conjugated systems, identical to ce
        - cg  sp             Inner Sp carbons in conjugated systems
        - ch  sp             Inner Sp carbons in conjugated systems, identical to cg
        - cx  sp3            Sp3 carbons in triangle systems
        - cy  sp3            Sp3 carbons in square systems
        - cu  sp2            Sp2 carbons in triangle systems
        - cv  sp2            Sp2 carbons in square systems
        - n   sp2            Sp2 nitrogen in amide groups
        - n1  sp             Sp N
        - n2  sp2            aliphatic Sp2 N with two connected atoms
        - n3  sp3            Sp3 N with three connected atoms
        - n4  sp3            Sp3 N with four connected atoms
        - na  sp2            Sp2 N with three connected atoms
        - nb  sp2            Sp2 N in pure aromatic systems
        - nc  sp2            Sp2 N in non-pure aromatic systems
        - nd  sp2            Sp2 N in non-pure aromatic systems, identical to nc
        - ne  sp2            Inner Sp2 N in conjugated systems
        - nf  sp2            Inner Sp2 N in conjugated systems, identical to ne
        - nh  sp3            Amine N connected one or more aromatic rings
        - no  sp2            Nitro N
        - o   sp2            Oxygen with one connected atom
        - oh  sp3            Oxygen in hydroxyl group
        - os  sp3            Ether and ester oxygen
        - ow  sp3            Oxygen in water
        - p2  sp3d           Phosphate with two connected atoms
        - p3  sp3d           Phosphate with three connected atoms, such as PH3
        - p4  sp3d           Phosphate with three connected atoms, such as O=P(CH3)2
        - p5  sp3d           Phosphate with four connected atoms, such as O=P(OH)3
        - pb  sp2            Sp2 P in pure aromatic systems
        - pc  sp2            Sp2 P in non-pure aromatic systems
        - pd  sp2            Sp2 P in non-pure aromatic systems, identical to pc
        - pe  sp2            Inner Sp2 P in conjugated systems
        - pf  sp2            Inner Sp2 P in conjugated systems, identical to pe
        - px  sp2            Special p4 in conjugated systems
        - py  sp2            Special p5 in conjugated systems
        - s   sp2            S with one connected atom
        - s2  sp2            S with two connected atom, involved at least one double bond
        - s4  sp2            S with three connected atoms
        - s6  sp3            S with four connected atoms
        - sh  sp3            Sp3 S connected with hydrogen
        - ss  sp3            Sp3 S in thio-ester and thio-ether
        - sx  sp2            Special s4 in conjugated systems
        - sy  sp2            Special s6 in conjugated systems
        - h1                 H bonded to aliphatic carbon with 1 electrwd. group
        - h2                 H bonded to aliphatic carbon with 2 electrwd. group
        - h3                 H bonded to aliphatic carbon with 3 electrwd. group
        - h4                 H bonded to non-sp3 carbon with 1 electrwd. group
        - h5                 H bonded to non-sp3 carbon with 2 electrwd. group
        - ha                 H bonded to aromatic carbon
        - hc                 H bonded to aliphatic carbon without electrwd. group
        - hn                 H bonded to nitrogen atoms
        - ho                 Hydroxyl group
        - hp                 H bonded to phosphate
        - hs                 Hydrogen bonded to sulphur
        - hw                 Hydrogen in water
        - hx                 H bonded to C next to positively charged group
        - f                  Fluorine
        - cl                 Chlorine
        - br                 Bromine
        - i                  Iodine

       Meng/Chimera Atom Types
         Type  Description
        - C3   sp3-hybridized Carbon
        - C2   sp2-hybridized Carbon
        - Car  aromatic Carbon **
        - Cac  carboxylate Carbon
        - C1   sp-hybridized Carbon
        - N3+  sp3-hybridized Nitrogen, formal +ve charge
        - N3   sp3-hybridized Nitrogen, neutral
        - N2   sp2-hybridized, double bonded or aromatic N **
        - Npl  sp2-hybridized but not double bonded N (planar) **
        - Ng+  Guanidinium Nitrogen, partial +ve charge
        - Ntr  Nitro Nitrogen
        - Nox  N-oxide Nitrogen
        - N1   sp-hybridized Nitrogen
        - O3   sp3-hybridized Oxygen
        - O2   sp2-hybridized Oxygen
        - Oar  aromatic O **
        - O3-  resonance-equivalent terminal oxygen on tetrahedral center (phosphate, sulfate, etc) **
        - O2-  resonance-quivalent terminal oxygen on planar center (carboxylate, nitro, nitrate) **
        - S3+  sp3-hybridized Sulfur, formal +ve charge
        - S3   sp3-hybridized Sulfur, neutral
        - S2   sp2-hybridized Sulfur
        - Sac  sulfate, sulfonate, or sulfamate Sulfur
        - Son  sulfone sulfur (>SO2) **
        - Sxd  sulfoxide sulfur (>SO) **
        - S    other Sulfur
        - Bac  borate Boron
        - Box  other oxidized Boron
        - B    other Boron
        - P3+  sp3-hybridized Phosphorous, formal +ve charge
        - Pac  phosphate, phosphonate, or phosphamate Phosphorous
        - Pox  P-oxide Phosphorous
        - P    other Phosphorous
        - HC   Hydrogen bonded to Carbon
        - H    other Hydrogen
    */
    std::string              itsMengType;

    /*!
        hydrophilic type Definitions
        - -1  Undefined
        - 0   Hydrophobic
        - 1   N,O
        - 2   S in SH
        - 3   S with double bond
        - 4   <= 1 bond away from charged atom
        - 5   <= 1 bond away from OH, NH, or NH2 with no delocalized electrons
        - 6   <= 1 bond away from SH with no delocalized electrons
        - 7   <= 1 bond away from O with double bond
        - 8   <= 1 bond away from S with valence > 2
        - 9   > one neighboring O or N with no delocalized electrons
          -# A two bond away from O with double bond
          -# B one bond away from S with valence > 2
        - 10  two or more instance of any of the previous two(A,B)
    */
    int                      itsHydrophilicity;

    /*!
        occupancy
    */
    double                   itsOccupancy;

    /*!
        hydrophilic type Definitions
    */
    double                   itsTempFactor;

protected: // - Molecular Mechanics data
    //! standard atom pointer
    stdAtom*                 pStdAtom;

};

} // MTKpp namespace

#endif // ATOM_H


