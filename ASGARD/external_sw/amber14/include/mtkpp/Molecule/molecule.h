/*!
   \file molecule.h
   \brief Container for submolecules, bonds, angles, torsions, and impropers
   \author Martin Peters

   Container for submolecules, bonds, angles, torsions, and impropers

   $Date: 2010/04/29 18:59:17 $
   $Revision: 1.28 $

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

#ifndef MOLECULE_H
#define MOLECULE_H

#include "Utils/constants.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <sstream>

namespace MTKpp
{
// ============================================================
// - classes
class collection;
class submolecule;
class atom;
class vector3d;
class rings;
class conformers;
class protonate;
class fingerPrint;
class functionalize;
class hydrophobize;
class pharmacophore;
class stdLibrary;
class stdFrag;
struct stdAtom;
class parameters;

// ------------------------------------------------------------
// - structs
struct moleculeFlags;
struct Bond;
struct Angle;
struct Torsion;
struct Improper;
struct ring;
struct funcGroup;
struct hydrophobe;
struct conformer;
struct bondParam;
struct angleParam;
struct torsionParam;
struct improperParam;
struct hBondParam;
struct LJ612Param;
struct LJ612SE;
struct clique;

// ============================================================
// Class : molecule()
// ------------------------------------------------------------
/*!
   \class molecule
   \brief Container for submolecules, bonds, angles, torsions, and impropers
*/
// ============================================================

class molecule
{
public:

    /*!
       \brief Molecule Constructor
       \param parent collection pointer
    */
    molecule(collection *parent = 0);

    //! Molecule Destructor
    virtual ~molecule();

    //-----------------------//
    // - GENERAL FUNCTIONS - //
    //-----------------------//

    /*!
       \brief Get collection which molecule is apart of
       \return collection pointer
    */
    collection*              getParent();

    /*!
       \brief Set name of molecule
       \param name molecule name
    */
    void                     setName(const std::string& name);

    /*!
       \brief Get name of molecule
       \return molecule name
    */
    std::string              getName();

    /*!
       \brief Set chain of molecule
       \param chain molecule chain
    */
    void                     setChain(const std::string& chain);

    /*!
       \brief Get chain of molecule
       \return chain molecule chain
    */
    std::string              getChain();

    /*!
       \brief Set kind of molecule
       \param n molecule kind
    */
    void                     setKind(const int& n);

    /*!
       \brief Get kind of molecule
       \return molecule kind
    */
    int                      getKind();

    /*!
       \brief Set id of molecule
       \param id molecule id
    */
    void                     setMolId(const int& id);

    /*!
       \brief Get id of molecule
       \return molecule id
    */
    int                      getMolId();

    /*!
       \brief Set number of atoms in molecule
       \param natoms number of atoms
    */
    void                     setNumAtoms(const int& natoms);

    /*!
       \brief Get number of atoms in molecule
       \return number of atoms in molecule
    */
    int                      getNumAtoms();

    /*!
       \brief Get number of atoms in molecule
       \return number of atoms in molecule
    */
    int                      getNumHeavyAtoms();

    /*!
       \brief Set number of bonds in molecule
       \param nbonds number of bonds
    */
    void                     setNumBonds(const int& nbonds);

    /*!
       \brief Get number of bonds in molecule
       \return number of bonds in molecule
    */
    int                      getNumBonds();

    /*!
       \brief Set atom index
       \param n atom index
    */
    void                     setAtomIndex(const int& n);

    /*!
       \brief Get atom index
       \return atom index
    */
    int                      getAtomIndex();

    /*!
       \brief Set maximum file id
       \param n maxFileID
    */
    void                     setMaxFileID(const int& n);

    /*!
       \brief Get maximum file id
       \return max file id
    */
    int                      getMaxFileID();

    /*!
       \brief Get number of electrons in molecule
       \return number of electrons in molecule
    */
    int                      getNumElectrons();

    /*!
       \brief Move the center of mass
       \param center New center of mass
    */
    void                     moveCenterOfMass(vector3d* center);

    /*!
       \brief Calculate the center of mass
       \param center Center of mass
    */
    void                     centerOfMass(vector3d* center);

    /*!
       \brief Calculate the center of mass
       \param center Center of mass
    */
    void                     centerOfMass(double center[3]);

    /*!
       \brief Set atomic coordinates
       \param coords vector of vector3d's
    */
    void                     setCoordinates(std::vector< vector3d > &coords);

    /*!
       \brief Get all atomic coordinates
       \param coords double array of coordinates
    */
    void                     getCoordinates(double coords[][3]);

    /*!
       \brief Get all atomic coordinates
       \param coords array of coordinates [x1, y1, z1, x2, y2, z2, ..., xN, yN, zN]
    */
    int                      getCoordinates(double coords[]);

    /*!
       \brief Get vector of atomic coordinates
       \param coords vector of vector3d's
    */
    void                     getCoordinates(std::vector< vector3d > &coords);

    /*!
       \brief Get heavy atom indices
       \return int array of indices
    */
    int*                     getHeavyAtomIndices();

    //
    //   \brief Get all atomic symbols
    //   \param symbols char array of symbols
    //
    //void                     getAtomSymbols(char symbols[][2]);

    /*!
       \brief Get all atomic symbols
       \return char array of symbols
    */
    char*                    getAtomSymbols();

    /*!
       \brief Get residue 1-Letter Symbols
       \param codes1L char array of 1-L symbols
    */
    void                     getRes1LSymbols(char codes1L[]);

    /*!
       \brief Get residue 3-Letter Symbols
       \param codes3L char array of 3-L symbols
    */
    void                     getRes3LSymbols(char codes3L[]);

    /*!
       \brief Get all atomic symbols
       \return char array of symbols
    */
    char*                    getHeavyAtomSymbols();

    /*!
       \brief Get all atom types
       \return char array of types
    */
    char*                    getAtomTypes();

    /*!
       \brief Get heavy atom types
       \return char array of types
    */
    char*                    getHeavyAtomTypes();

    /*!
       \brief Get atom kinds
       \return atom kinds
    */
    int*                     getAtomKinds();

    /*!
       \brief Get all atomic charges
       \param charges double array of charges
    */
    void                     getAtomCharges(double charges[]);

    /*!
       \brief Get total formal charge
       \return total charge
    */
    int                      getFormalCharge();

    /*!
       \brief Get total charge from stdFragments
       \return total charge
    */
    double                   getStdCharge();

    /*!
       \brief Set Total charge
       \param c molecule charge
    */
    void                     setTotalCharge(int c);

    /*!
       \brief Get Total charge
       \return molecule charge
    */
    int                      getTotalCharge();

    /*!
       \brief Get Molecular Weight
       \return molecular weight
    */
    double                      getMolecularWeight();

    //---------------------------//
    // - SUBMOLECULE FUNCTIONS - //
    //---------------------------//

    /*!
       \brief Add submolecule to molecule
       \return submolecule pointer
    */
    virtual submolecule*     addSubMolecule();

    /*!
       \brief Get submolecule from list of submolecules based on id
       \param id submolecule id
       \return submolecule pointer
       \sa getSubMoleculeList()
    */
    submolecule*             getSubMolecule(int id);

    /*!
       \brief Get submolecule from list of submolecules based on id
       \param id submolecule id
       \param smolIndex return submolecule based on index (default)
       \param fileId return submolecule based on file id
       \return submolecule pointer
    */
    submolecule*             getSubMolecule(int id, bool smolIndex, bool fileId);

    /*!
       \brief Get submolecule from list of submolecules based on name
       \param name submolecule name
       \return submolecule pointer
    */
    submolecule*             getSubMolecule(std::string name);

    /*!
       \brief Get list of submolecules
       \return vector of submolecule pointers
    */
    virtual std::vector<submolecule*> getSubMoleculeList();

    /*!
       \brief Get list of submolecules by name
       \param name submolecule name
       \return vector of submolecule pointers
    */
    std::vector<submolecule*> getSubMoleculeList(std::string name);

    /*!
       \brief Get list of submolecules by name and num
       \param name submolecule name
       \param resId submolecule number
       \return vector of submolecule pointers
    */
    std::vector<submolecule*> getSubMoleculeList(std::string name, int resId);

    /*!
       \brief Get list of submolecules by name
       \return vector of submolecule pointers
    */
    std::vector<submolecule*> getSubMoleculeList(int resId);

    /*!
       \brief Get number of submolecules
       \return number of submolecules in the molecule
    */
    int                      getNumSubMolecules();

    //--------------------//
    // - ATOM FUNCTIONS - //
    //--------------------//

    /*!
       \brief Get list of atoms
       \return vector of atom pointers
    */
    std::vector<atom*>       getAtomList();

    /*!
       \brief Get list of atoms
       \return vector of atom pointers
    */
    std::vector<atom*>       getHeavyAtomList();

    /*!
       \brief Get atom
       \param number integer value
       \param atomIndex return atom based on index (default)
       \param fileId return atom based on file id
       \param atomColIndex atom's collection index
       \return atom pointer
    */
    atom*                    getAtom(int number, bool atomIndex = true,
                                                 bool fileId = false,
                                                 bool atomColIndex = false);

    /*!
       \brief Get atom
       \param pStdAtom standard atom pointer
       \return atom pointer
    */
    atom*                    getAtom(stdAtom* pStdAtom);

    /*!
       \brief Get atom
       \param name atom name
       \return atom pointer
    */
    atom*                    getAtom(std::string name);

    /*!
       \brief Delete atom
       \param a atom pointer
    */
    void                     delAtom(atom* a);

    //----------//
    // - BOND - //
    //----------//

    /*!
       \brief Add Bond to molecule
       \param at1 atom pointer
       \param at2 atom pointer
       \param type integer value designating single, double, triple etc., bond type
       \param stereo integer value designating bond stereo
       \param topology integer value assigning ring membership or not
       \param size distance between at1 and at2
       \return Bond pointer
    */
    virtual Bond*            addBond(atom* at1, atom* at2,
                                     const int& type = 0,
                                     const int& stereo = 0,
                                     const int& topology = 0,
                                     const double& size = 0.0);

    /*!
       \brief Get Bond
       \param at1 atom pointer
       \param at2 atom pointer
       \return Bond pointer
    */
    Bond*                    getBond(atom* at1, atom* at2);

    /*!
       \brief Get bond list
       \param bonds array of atom indices
    */
    int                      getBonds(int bonds[]);

    /*!
       \brief Get bond parameters
       \param bondParams array of bond parameters
    */
    int                      getBondParams(double bondParams[]);

    /*!
       \brief Does molecule contain a certain bond
       \param at1 atom pointer
       \param at2 atom pointer
       \return boolean
    */
    bool                     hasBond(atom* at1, atom* at2);

    /*!
       \brief Get number of bonds
       \return number of bonds
    */
    int                      numBonds();

    /*!
       \brief Deletes bond from the molecules' list of bonds
       \param at1 atom pointer
       \param at2 atom pointer
    */
    void                     delBond(atom* at1,atom* at2);

    /*!
       \brief Add 1-3 Bond to molecule
       \param at1 atom pointer
       \param at2 atom pointer
    */
    void                     add13Bond(atom* at1, atom* at2);

    /*!
       \brief Does molecule contain a certain 13 bond
       \param at1 atom pointer
       \param at2 atom pointer
       \return boolean
    */
    bool                     has13Bond(atom* at1, atom* at2);

    /*!
       \brief Add 1-4 Bond to molecule
       \param at1 atom pointer
       \param at2 atom pointer
    */
    void                     add14Bond(atom* at1, atom* at2);

    /*!
       \brief Does molecule contain a certain 14 bond
       \param at1 atom pointer
       \param at2 atom pointer
       \return boolean
    */
    bool                     has14Bond(atom* at1, atom* at2);

    /*!
       \brief Get Bond map
       \return bond map <int Bond*>
    */
    std::map<int, Bond*>     getBondMap();

    //-----------//
    // - ANGLE - //
    //-----------//

    /*!
       \brief Add Angle to molecule
       \param at1 atom pointer
       \param at2 atom pointer
       \param at3 atom pointer
       \param ang Angle between at1, at2, and at3
       \return Angle pointer
    */
    Angle*                   addAngle(atom* at1,atom* at2,atom* at3,
                                      const double& ang);

    /*!
       \brief Does molecule contain a certain angle
       \param at1 atom pointer
       \param at2 atom pointer
       \param at3 atom pointer
       \return boolean
    */
    bool                     hasAngle(atom* at1,atom* at2,atom* at3);

    /*!
       \brief Get angle
       \param at1 atom pointer
       \param at2 atom pointer
       \param at3 atom pointer
       \return Angle pointer
    */
    Angle*                   getAngle(atom* at1,atom* at2,atom* at3);

    /*!
       \brief Deletes angle from the molecules' list of angles
       \param at1 atom pointer
       \param at2 atom pointer
       \param at3 atom pointer
    */
    void                     delAngle(atom* at1,atom* at2,atom* at3);

    /*!
       \brief Get Angle map
       \return angle map <ULONG_KIND Angle*>
    */
    std::map<ULONG_KIND, Angle*>    getAngleMap();

    /*!
       \brief Get number of angles
       \return number of angles
    */
    int                      numAngles();

    /*!
       \brief Get angle list
       \param angles array of atom indices
    */
    int                      getAngles(int angles[]);

    /*!
       \brief Get angle parameter list
       \param angleParams array of angle parameters
    */
    int                      getAngleParams(double angleParams[]);

    //-------------//
    // - TORSION - //
    //-------------//

    /*!
       \brief Add Torsion to molecule
       \param at1 atom pointer
       \param at2 atom pointer
       \param at3 atom pointer
       \param at4 atom pointer
       \param tor Torsion between at1, at2, at3, and at4
       \return Torsion pointer
    */
    Torsion*                 addTorsion(atom* at1,atom* at2,
                                        atom* at3,atom* at4,
                                        const double& tor);

    /*!
       \brief Deletes torsion from the molecules' list of torsions
       \param at1 atom pointer
       \param at2 atom pointer
       \param at3 atom pointer
       \param at4 atom pointer
    */
    void                     delTorsion(atom* at1,atom* at2,
                                        atom* at3,atom* at4);

    /*!
       \brief Does molecule contain a certain torsion
       \param at1 atom pointer
       \param at2 atom pointer
       \param at3 atom pointer
       \param at4 atom pointer
       \return boolean
    */
    bool                     hasTorsion(atom* at1,atom* at2,
                                        atom* at3,atom* at4);

    /*!
       \brief Get torsion
       \param at1 atom pointer
       \param at2 atom pointer
       \param at3 atom pointer
       \param at4 atom pointer
       \return Torsion pointer
    */
    Torsion*                 getTorsion(atom* at1,atom* at2,
                                        atom* at3,atom* at4);

    /*!
       \brief Get torsion
       \param at1 atom pointer
       \param at2 atom pointer
       \param at3 atom pointer
       \param at4 atom pointer
       \param t torsion size
    */
    void                     setTorsion(atom* at1,atom* at2,
                                        atom* at3,atom* at4,
                                        double t);

    /*!
       \brief Get Torsion map
       \return torsion map <int Torsion*>
    */
    std::map<ULONG_KIND, Torsion*>  getTorsionMap();

    /*!
       \brief Get number of torsions
       \return number of torsions
    */
    int                      numTorsions();

    /*!
       \brief Get torsion list
       \param torsions array of atom indices
       \return success
    */
    int                      getTorsions(int torsions[]);

    /*!
       \brief Get number of torsion parameters
       \return number of torsions
    */
    int                      getNumMMTorsions();

    /*!
       \brief Get torsion parameters
       \param torsions array of atom indices
       \param torsionParams array of torsion parameters
       \return success
    */
    int                      getMMTorsions(int torsions[], double torsionParams[]);


    //--------------//
    // - IMPROPER - //
    //--------------//

    /*!
       \brief Add Improper to molecule
       \param at1 atom pointer
       \param at2 atom pointer
       \param at3 atom pointer
       \param at4 atom pointer
       \param imp Improper between at1, at2, at3, and at4
       \return Improper pointer
    */
    Improper*                addImproper(atom* at1,atom* at2,atom* at3,atom* at4,const double& imp);

    /*!
       \brief Deletes improper from the molecules' list of impropers
       \param at1 atom pointer
       \param at2 atom pointer
       \param at3 atom pointer
       \param at4 atom pointer
    */
    void                     delImproper(atom* at1,atom* at2,atom* at3,atom* at4);

    /*!
       \brief Does molecule contain a certain improper
       \param at1 atom pointer
       \param at2 atom pointer
       \param at3 atom pointer
       \param at4 atom pointer
       \return boolean
    */
    bool                     hasImproper(atom* at1,atom* at2,atom* at3,atom* at4);

    /*!
       \brief Get Improper map
       \return improper map <int Improper*>
    */
    std::map<int, Improper*> getImproperMap();

    /*!
       \brief Get number of impropers
       \return number of impropers
    */
    int                      numImpropers();

    /*!
       \brief Get improper list
       \param impropers array of atom indices
    */
    int                      getImpropers(int impropers[]);

    /*!
       \brief Get number of improper parameters
       \return number of MM impropers
    */
    int                      getNumMMImpropers();

    /*!
       \brief Get impropers and parameters
       \param impropers arrary of impropers
       \param improperParams array of improper parameters
       \return ??
    */
    int                      getMMImpropers(int impropers[], double improperParams[]);

    //-----------------------------//
    // - NON BONDED INTERACTIONS - //
    //-----------------------------//

    /*!
       \brief Get number of MM non bonded pairs
    */
    int                      getNumMMnonBondedPairs(double cutoff);

    /*!
       \brief Get MM non bonded pairs
    */
    int                      getMMnonBondedPairs(int nonBonded[], int nonBondedPtrs[], double nonBondedParams[],
                                                 int nonBonded14Ptrs[], double cutoff);

    //----------------//
    // - MM CHARGES - //
    //----------------//

    /*!
       \brief Get charges list
       \param charges array of atom charges
    */
    int                      getMMCharges(double charges[]);

    //-----------------//
    // - MM RESIDUES - //
    //-----------------//

    /*!
       \brief Get residue list
       \param residues array of residue indices
    */
    int                      getMMResidues(int residues[]);

    //-----------//
    // - RINGS - //
    //-----------//

    //! Determine rings
    void                     determineRings();

    //! Determines if rings is aromatic or not
    void                     kekulizeRings();

    //! Determines if rings is aromatic or not
    std::string              getRingInfo();

    /*!
       \brief Add a new ring
       \param r vector of atoms in the ring
    */
    void                     addRing(std::vector<atom*> r);

    /*!
       \brief Get rings
       \return vector of rings
    */
    std::vector<ring*>       getRings();

    /*!
       \brief Get last added ring
       \return ring pointer
    */
    ring*                    getLastAddedRing();

    //----------------//
    // - CONFORMERS - //
    //----------------//

    /*!
       \brief Get rotatable bonds
       \todo
          Definition of Rotatable Bonds
          A bond is considered as being rotatable, if:
           - it is not in a ring, or if it is in a ring larger than 8 --- done
           - it is a single or a triple bond --- only single bond
           - none of the two atoms is terminal  --- done
           - it is not adjacent to a triple bond  ----
           - it is not an amide, amidine, guanidine bond ---
       \return rotatable bond vector
    */
    std::vector<Bond*>       getRotatableBonds();

    /*!
       \brief Get rotatable dihedrals
       \param tors rotatable dihedrals vector
    */
    void                     setRotatableTorsions(std::vector<Torsion*> tors);

    /*!
       \brief Get rotatable dihedrals
       \return rotatable dihedrals vector
    */
    std::vector<Torsion*>    getRotatableTorsions();

    /*!
       \brief Set up conformers
       \return vector of torsions
    */
    std::vector<Torsion*>    setupConformers();

    /*!
       \brief Set up conformers
       \param frozenAtoms Atoms which cannot move
       \return vector of torsions
    */
    std::vector<Torsion*>    setupConformers(std::vector<atom*> frozenAtoms);

    /*!
       \brief Set up conformers
       \param rotBonds user defined rotatable bonds
       \return vector of torsions
    */
    std::vector<Torsion*>    setupConformers(std::vector<Bond*> rotBonds);

    /*!
       \brief Set up conformers
       \param rotTorsions user defined rotatable torsions
       \return vector of torsions
    */
    std::vector<Torsion*>    setupConformers(std::vector<Torsion*> rotTorsions);

    /*!
       \brief Add conformer to the list of conformers
       \param name name of conformer
       \param t vector of torsions in radians
    */
    void                     addConformer(std::string name, std::vector<double> t);

    /*!
       \brief Updates the conformers torsional angles
       \param name name of conformer
       \param t vector of torsions in radians
    */
    void                     updateConformer(std::string name, std::vector<double> t);

    /*!
       \brief Delete conformer from the list of conformers
       \param name name of conformer
    */
    void                     delConformer(std::string name);

    /*!
       \brief Get conformer based on its name
       \param name conformer name
       \return conformer pointer
    */
    conformer*               getConformer(std::string name);

    /*!
       \brief Get conformers pointer
       \return conformers pointer
    */
    conformers*              getConformers();

    /*!
       \brief Get conformer names
       \return conformer names
    */
    std::vector<std::string> getConformerNames();

    //------------------//
    // - FINGERPRINTS - //
    //------------------//

    /*!
       \brief generate simple fingerprint
    */
    void                     generateSimpleFP();

    /*!
       \brief get simple fingerprint
       \return simple fingerprint
    */
    std::vector<unsigned int>getSimpleFP();

    /*
       \brief Generate fragment fingerprint
    */
    //void                     generateFragmentFP();

    /*!
       \brief Get fragment fingerprint
       \return fragment fingerprint
    */
    std::string              getFragmentFP();

    //------------------------//
    // - HYDROPHOBIC GROUPS - //
    //------------------------//

    /*!
       \brief Determine hydrophobic groups in the molecule
       \return success
    */
    int determineHydrophobicGroups();

    /*!
       \brief Add hydrophobic groups to the molecule
       \param ats atom in hydrophobe
       \return success
    */
    int addHydrophobicGroup(std::vector<atom*> ats);

    /*!
       \brief Get hydrophobic groups
       \return vector of hydrophobes
    */
    std::vector<hydrophobe*> getHydrophobicGroups();

    //-----------------------//
    // - FUNCTIONAL GROUPS - //
    //-----------------------//

    /*!
       \brief Determine functional groups in the molecule
       \param pStdLib stdLib pointer
       \return success
    */
    int                      determineFunctionalGroups(stdLibrary* pStdLib);

    /*!
       \brief Add functional group to the molecule
       \param fg funcGroup pointer
    */
    void                     addFunctionalGroup(funcGroup* fg);

    /*!
       \brief Add functional group to the molecule
       \param fragAtoms Number of atoms in fragment
       \param pStdFrag stdFrag pointer
       \param molAtoms Number of atoms in molecule
       \param subGraph subGraph from functionalize
    */
    void                     addFunctionalGroup(int fragAtoms, stdFrag* pStdFrag,
                                                int molAtoms, std::vector<int> subGraph);

    /*!
       \brief Delete functional group from the molecule
       \param fG functional group pointer
       \param i number of fG
    */
    void                     delFunctionalGroup(funcGroup* fG, int i);

    /*!
       \brief Get functional groups
       \return functional groups
    */
    std::vector<funcGroup*>  getFunctionalGroups();

    /*!
       \brief Get functional group to which the atom is apart of
       \param pAt atom pointer
       \return functional group
    */
    funcGroup*               getFunctionalGroup(atom* pAt);

    /*!
       \brief Determine adjacency matrix
       \return success
    */
    int                      generateAdjMatrix();

    /*!
       \brief Determine heavy atom adjacency matrix
       \return success
    */
    int                      generateHeavyAdjMatrix();

    /*!
       \brief Get adjacency matrix
       \return adjacency matrix
    */
    int*                     getAdjMatrix();

    /*!
       \brief Get heavy atom adjacency matrix
       \return heavy atom adjacency matrix
    */
    int*                     getHeavyAdjMatrix();

    /*!
       \brief Get adjacency matrix Size
       \return adjacency matrix size
    */
    int                      getAdjMatrixSize();

    /*!
       \brief Get heavy adjacency matrix Size
       \return adjacency matrix size
    */
    int                      getHeavyAdjMatrixSize();

    /*!
       \brief Determine molecular feature distance matrix
       \return success
    */
    int                      generateFeatureDistMatrix();

    /*!
       \brief Determine molecular feature distance matrix
       \return success
    */
    int                      generateFeatureDistMatrix(std::vector< vector3d >
                                                       coords);

    /*!
       \brief Determine molecular feature distance matrix
       \return success
    */
    double*                  getFeatureDistMatrix();

    /*!
       \brief Get feature groups
       \return feature groups
    */
    std::vector<std::string> getFeatureGroups();

    /*!
       \brief Get feature labels
       \return feature labels
    */
    std::vector<std::string> getFeatureLabels();

    /*!
       \brief Get feature coordinates
       \return feature coordinates
    */
    std::vector<vector3d>    getFeatureCoordinates();

    /*!
       \brief Determine the Maximum Common Pharmacophore (MCP) between two molecules
       \param pMol molecule pointer
       \param mcp atom indices of the maximum common pharmacophore
       \param mcpCoords The maximum common pharmacophore coordinates
       \param distMax Feature distance parameter
       \return success
    */
    int                      findPharmacophore(molecule* pMol,
                             std::vector<std::vector<unsigned int> > &mcp,
                             std::vector<vector3d> &mcpCoords,
                             double distMax = 0.5);

    /*!
       \brief Determine all Pharmacophores between two molecules
       \param pMol molecule pointer
       \param cliqueList List of cliques
       \param distMax Feature distance parameter
       \return success
    */
    int                      findPharmacophore(molecule* pMol,
                             std::vector<clique*> &cliqueList,
                             double distMax = 0.5);


    //----------------//
    // - PROPERTIES - //
    //----------------//

    /*!
       \brief Add molecular property
       \param name name of property
       \param value value of property
    */
    void                     addProperty(const std::string& name,
                                         const std::string& value);

    /*!
       \brief Get molecular property
       \param name name of property
       \return value of property
    */
    std::string              getProperty(const std::string& name);

    /*!
       \brief Delete molecular property
       \param name name of property
    */
    void                     delProperty(const std::string& name);

    /*!
       \brief get molecular properties
       \return map of molecular properties
    */
    std::map<std::string, std::string>  getProperties();

    /*!
       \brief has molecule a certain property
       \param name name of property
       \return got it or not
    */
    bool                     hasProperty(const std::string& name);

    /*
      \brief Set protonation state of molecule

      deprotonate carboxylic, phosphonic, and sulfonic acids
      protonates quaternary nitrogens, amidine and guanidine moieties
    */
    void                     setInSolution();

    //-----------------//
    // -  HYDROGENS  - //
    //-----------------//

    /*!
       \brief Add hydrogens
    */
    void                     addHydrogens();

    /*!
       \brief Remove Hydrogens
    */
     void                    removeHydrogens();

    //----------//
    // - MISC - //
    //----------//

    /*!
       \brief Determine atom valences
    */
    void                     determineValences();

    /*!
       \brief Determine atom hybridizations
       \param method 1 - Meng, 2 - Labute
       \param f parameter file name
    */
    void                     determineHybridizations(int method, const std::string f);

    //-----------------//
    // -    FLAGS    - //
    //-----------------//

    //! Input file type, mol, mol2, pdb, etc.
    std::string inFileType;

    //! Are the bonds assigned?
    bool bBondsAssigned;

    //! Are the disulfide bonds assigned?
    bool bDisulfideBondsAssigned;

    //! Are the angles assigned?
    bool bAnglesAssigned;

    //! Are the torsions assigned?
    bool bTorsionsAssigned;

    //! Are the impropers assigned?
    bool bImpropersAssigned;

    //! Are the bond types: Single, double, triple, assigned?
    bool bBondTypes1Assigned;

    //! Are the Aromatic bonds assigned?
    bool bBondTypes2Assigned;

    //! Are the MM atom types assigned
    bool bMMAtomTypesAssigned;

    //! Are the atomic valences assigned?
    bool bAtomValencesAssigned;

    //! Are the atomic hybridizations assigned?
    bool bAtomHybridizationsAssigned;

    //! Were the rings assigned?
    bool bRingsAssigned;

    //! Were the rings kekulized?
    bool bRingsKekulized;

    //! Is the molecule considered in solution?
    bool bInSolution;

    //! Were hydrogens added?
    bool bHydrogensAdded;

    //! Error Flag
    bool bError;

protected:

    //-------------------//
    // - MOLECULE DATA - //
    //-------------------//

    //! molecule name
    std::string              itsName;

    //! molecule chain
    std::string              itsChain;

    /*!
       molecule kind definitions
       - 0 = undefined
       - 1 = biomolecule, e.g protein, DNA, etc
       - 2 = small molecule, e.g. drug
       - 3 = other
    */
    int                      itsKind;

    //! molecule id
    int                      itsMolId;

    //! number of atoms
    int                      itsNumAtoms;

    //! number of bonds
    int                      itsNumBonds;

    //! submolecule index
    int                      itsSMolIndex;

    //! atom index
    int                      itsAtomIndex;

    //! maximum atom file id
    int                      itsMaxFileID;

    //! Total Charge
    int                      totalCharge;

    //! center of mass
    vector3d*                pCenterMass;

    //-----------------------//
    // - BACKBONE POINTERS - //
    //-----------------------//

    //! molecules' collection
    collection*              pParent;

    //! submolecule pointer
    submolecule*             pSubMolecule;

    //! atom pointer
    atom*                    pAtom;

    //! Bond pointer
    Bond*                    pBond;

    //! parameters pointer
    parameters*              pParameters;

    //! Bond parameter pointer
    bondParam*               pBondParam;

    //! Angle parameter pointer
    angleParam*              pAngleParam;

    //! Torsion parameter pointer
    torsionParam*            pTorsionParam;

    //! H Bond Param
    hBondParam*              pHBondParam;

    //! LJ 6-12 Param
    LJ612Param*              pLJ612Param;

    //! LJ 6-12 Sigma and Epsilon
    LJ612SE*                 pLJ612SE;

    //! Torsion parameter list
    std::vector<torsionParam*> torsionParamList;

    //! Torsion parameter iterator
    typedef std::vector<torsionParam*>::iterator torsionParamIterator;

    //! Improper parameter list
    std::vector<improperParam*> improperParamList;

    //! Improper parameter iterator
    typedef std::vector<improperParam*>::iterator improperParamIterator;

    //! Improper parameter pointer
    improperParam*           pImproperParam;

    //! Angle pointer
    Angle*                   pAngle;

    //! Torsion pointer
    Torsion*                 pTorsion;

    //! Improper pointer
    Improper*                pImproper;

    //! rings pointer
    rings*                   pRings;

    //! ring pointer
    ring*                    pRing;

    //! conformers pointer
    conformers*              pConformers;

    //! conformer pointer
    conformer*               pConformer;

    //! protonate pointer
    protonate*               pProtonate;

    //! fingerPrint pointer
    fingerPrint*             pFingerPrint;

    //! functionalize pointer
    functionalize*           pFunctionalize;

    //! hydrophobize pointer
    hydrophobize*            pHydrophobize;

    //! hydrophobe pointer
    hydrophobe*              pHydrophobe;

    //! pharmacophore pointer
    pharmacophore*           pPharmacophore;

    //----------------------//
    // - VECTOR ITERATORS - //
    //----------------------//

    //! submolecule iterator
    typedef std::vector<submolecule*>::iterator  sMolIterator;

    //! atom iterator
    typedef std::vector<atom*>::iterator         AtomIterator;

    //! Bond iterator
    typedef std::vector<Bond*>::iterator         BondIterator;

    //! conformer iterator
    typedef std::vector<conformer*>::iterator    ConformerIterator;

    //! ring iterator
    typedef std::vector<ring*>::iterator         RingIterator;

    //! funcGroup iterator
    typedef std::vector<funcGroup*>::iterator    FuncGroupIterator;

    //! hydrophobe iterator
    typedef std::vector<hydrophobe*>::iterator   HydrophobeIterator;

    //-------------------//
    // - MAP ITERATORS - //
    //-------------------//

    //! Bond map iterator
    typedef std::map<int, Bond*>::iterator       BondMapIterator;

    //! Angle map iterator
    typedef std::map<ULONG_KIND, Angle*>::iterator      AngleMapIterator;

    //! Torsion map iterator
    typedef std::map<ULONG_KIND, Torsion*>::iterator    TorsionMapIterator;

    //! Improper map iterator
    typedef std::map<int, Improper*>::iterator   ImproperMapIterator;

    //! property map iterator
    typedef std::map<std::string, std::string>::iterator PropertyMapIterator;

    //-----------------------//
    // - VECTOR CONTAINERS - //
    //-----------------------//

    //! submolecule vector
    std::vector<submolecule*>                    itsSubMoleculeList;

    //! atom vector
    std::vector<atom*>                           itsAtomList;

    //! heavy atom vector
    std::vector<atom*>                           itsHeavyAtomList;

    //! 13 Bond vector
    std::vector<int>                             itsBond13List;

    //! 14 Bond vector
    std::vector<int>                             itsBond14List;

    //! rotatable bonds
    std::vector<Bond*>                           itsRotatableBonds;

    //! rotatable torsions
    std::vector<Torsion*>                        itsRotatableTorsions;

    //! Ring vector
    std::vector<ring*>                           itsRings;

    //! funcGroup list
    std::vector<funcGroup*>                      itsFuncGroups;

    //! hydrophobe list
    std::vector<hydrophobe*>                     itsHydrophobes;

    //! Conformer vector
    std::vector<conformer*>                      itsConformers;

    //! Simple fingerprint
    std::vector<unsigned int>                    itsSimpleFP;

    //! Fragment fingerprint
    std::string                                  itsFragmentFP;

    //! 1D fingerprint
    //std::vector<unsigned int>                    its1DFP;

    //! Heavy Atom indices
    int                                          *heavyAtomIndices;

    //! Adjacency Matrix
    int                                          *adjMatrix;

    //! Heavy Atom Adjacency Matrix
    int                                          *heavyAdjMatrix;

    //! Adjacency Matrix Size
    int                                          adjMatrixSize;

    //! Heavy Adjacency Matrix Size
    int                                          heavyAdjMatrixSize;

    //! Feature Distance Matrix
    double                                       *featureDistMatrix;

    //! Adjacency Matrix Size
    int                                          featureDistMatrixSize;

    //! Feature Groups
    std::vector<std::string>                     featureGroups;

    //! Feature Labels
    std::vector<std::string>                     featureLabels;

    //! Feature coordinates
    std::vector<vector3d>                        featureCoordinates;

    //! Atom symbols
    char                                         *atomSymbols;

    //! Atom symbols
    char                                         *heavyAtomSymbols;

    //! Atom types
    char                                         *atomTypes;

    //! Heavy atom types
    char                                         *heavyAtomTypes;

    //! Atom kinds
    int                                          *atomKinds;

    //--------------------//
    // - MAP CONTAINERS - //
    //--------------------//

    //! Bond map
    std::map<int, Bond*>                         itsBondMap;

    //! Angle map
    std::map<ULONG_KIND, Angle*>                 itsAngleMap;

    //! Torsion map
    std::map<ULONG_KIND, Torsion*>               itsTorsionMap;

    //! Improper map
    std::map<int, Improper*>                     itsImproperMap;

    //! Molecule properties map
    std::map<std::string, std::string>           itsPropertiesMap;
};

} // MTKpp namespace

#endif // MOLECULE_H

