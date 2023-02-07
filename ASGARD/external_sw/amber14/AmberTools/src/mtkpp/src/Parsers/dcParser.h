/*! 
   \file dcParser.h
   \brief Parses Divcon files
   \author Martin Peters

   Reads divcon output files and write divcon input files

   $Date: 2010/03/29 20:39:34 $
   $Revision: 1.17 $

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

#ifndef DCPARSER_H
#define DCPARSER_H

#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

#include "baseParser.h"

namespace MTKpp
{

class collection;
class molecule;
class submolecule;
class atom;
class bond;
class vector3d;

// ============================================================
// Class : dcParser()
// ------------------------------------------------------------
/*! 
   \class dcParser
   \brief Reads divcon output files and write divcon input files
   \author Martin Peters
   \date 2005

    \section divconOptions_sec DivCon Input Options
       - Hamiltonians 
        -# MNDO
        -# AM1
        -# PM3
        -# PM3 ZnB : biological Zinc paramters to be used in pm3
        -# MNDO/d
        -# PDDG-PM3

       - Geometry optimizations
        -# steep: steepest descent
        -# conjgrad: conjgrate gradient
        -# bfgs: 
        -# lbfgs: limited memory BFGS
        -# maxopt: do a maximum of X cycles of geometry optimization
        -# ts: transition state geometry calculation
        -# tsprfo: baker method
        -# tsqna: quasi-newton method
        -# tsnr: pure newton method
        -# nomodefollow: mode following not activated in ts search

       - Charge Models Available
        -# mulliken: use mulliken for scrf/pme calculations
        -# cm1: use mulliken for scrf/pme calculations
        -# cm2: use mulliken for scrf/pme calculations

       - Self Consistent Reaction Field:
        -# scrf: requests a scrf calculation with Delphi
        -# scrf-inside: only fill cavities with a continuum
        -# drad: used with scrf-inside
        -# indi: internal dielectric set to X
        -# exdi: external dielectric set to X
        -# scale: grid space set to X grid per angstrom
        -# perfill: molecule occupies X % of the grid
        -# ion: ion strength set to X
        -# probrab: probe radius set to X
        -# water: use water params for scrf calculation (Gnp)
        -# octanol: use octanol params for scrf calculation (Gnp)
        -# divpb: use divpb instead of delphi for scrf calculation
        -# soldmx: use guess dmx for scrf calculation

       - Useful Keywords
        -# freq: carry out frequency analysis
        -# thermo: report thermodynamic quantities from freq calculations
        -# direct: do not store 2-e integrals in scf calculations
        -# residue: store residue pointers (needed for D&C)
        -# pdump: write a restart density matrix every X scf cycles
        -# dump: write a restart file every X opt cycles
        -# shift: initial dynamic level shift set to X
        -# trajectory: dump coordinates to a trajectory file at restart points
        -# cutbond: cutoff bonding beyond X Angstroms
        -# 1SCF: perform only one scf iteration
        -# guess: build initial density matrix from one or more dmx files
        -# filenames: allow unique file names
        -# vdw: calculate van der Waals interactions

       - Gradient  --- options not supported
        -# gradient: output final gradient
        -# central: use central difference in gradient calculation
        -# inter: include only intermolecular contributions
        -# recipintra: consider intramolecular contributions

       - Convergence criteria  --- options not supported
        -# descf: user defined scf energy criterion [descf == ecrit]
        -# dpscf: user defined scf density criterion [dpscf == dcrit]
        -# etest: energy criterion
        -# gtest: gradient criterion
        -# ecrit: energy criterion
        -# dcrit: density matrix criterion
        -# force-it: force geometry opt criteria to be satisfied

       - Subsetting   --- cluster & no-overlap supported
        -# cluster: cluster bases subsetting
        -# dbuff1:  buffer 1 size
        -# dbuff2:  buffer 2 size
        -# no-overlap: off diagonal terms of density matrix are set to zero
        -# atgrid: atom wise grid based subsetting
        -# resgrid: residue wise grid based subsetting
        -# mixgrid: residue wise, grid based subsetting for core 
                    and atom wise for grid based subsetting for buffers
        -# combsub: combination subsetting
        -# standard: no divide-and-conquer

        - Parameterization   ---- options not supported
         -# powell: use powell algorithm to optimize parameters
         -# ga: use genetic algorithm to optimize parameters
         -# files: number of files to loop
         -# testpoint: test best parameter set from ga results
         -# external: read parameters from a file

       - Monte Carlo:   ---- options not supported
        -# mc, mcecrit, mcdcrit, mctemp, mcpres, inter, recipintra,frozenmc 
        -# seed, smarta, dangle, dbox, dr, nmove, nsubmove
        -# upd, ad, ratio, nwrite, npt, nvt,

       - Particle Mesh Ewald:  ---- options not supported
        -# k1pme, k2pme, k3pme, pme,betapme-hbox, betapme, nspline

       - Output  ---- options not supported
        -# prtsub: prints subsystem atom list
        -# prtvec: prints final eigenvectors
        -# prtcoords: prints atomic coordinates
        -# prtpar: prints semi-empirical parameters
        -# prtvdwaals: output van der Waals radii
        -# snapgeom: print out every optimization step into separate files

       - Misc.  ---- options not supported
        -# fullscf: performs full diagonalizations in SCF
        -# addmm: adds MM correction to peptide rotational barrier [default]
        -# nomm: does not use addmm
        -# maxit: maximum number of scf iterations is 100
        -# double: do a double scf step for every scf iteration
        -# lsearch: linear search in lbfgs
        -# diagterm: initial diagonal term in lbfgs opt 
        -# intgls: stewart/talman to be used [talman]
        -# tempk: set temperature to X
        -# tmax: request a maximum CPU time of X seconds
        -# rmin: minimum distance between atoms is X 
        -# cutrepul: cut off for repulsion integrals
        -# chkres: check inter-residue distances for each residue
        -# cartesian: requires x,y,z data
        -# internal: requires a z-matrix
        -# dos: density of states
        -# dipole: calculates dipole moment
        -# ip: calculate ionization potential
        -# homolumo: calculates homo-lumo gap
        -# screen: output info to screen
        -# zmake: generates a z-matrix
        -# qmailgn: generate output for QM-QSAR
        -# rotate: rotational angle for barrier is X
        -# push: push cluster groups apart
*/
// ============================================================
class dcParser : public baseParser
{
public:

    /*!
       \brief dcParser Constructor
    */
    dcParser();

    //! dcParser Destructor
    virtual ~dcParser();

    /*!
       \brief Read DivCon formatted file
       \param i Input file
       \param m collection pointer
    */
    void Read(const std::string &i, molecule* m);

    /*!
       \brief Read DivCon formatted file
       \param i Input file
       \param c collection pointer
       \param m molecule id
    */
    //void Read(const std::string &i, collection* c, const int &m);

    /*!
       \brief Parse heat of formation from the divcon output file
       \param i Input file
       \return Heat of Formation
    */
    double ReadHOF(const std::string &i);

    /*!
       \brief Parse electronic energy from the divcon output file
       \param i Input file
       \return Electronc energy
    */
    double ReadElecEnergy(const std::string &i);

    /*!
       - Hamiltonians
        -# MNDO
        -# AM1
        -# PM3
        -# PM3 ZnB : biological Zinc paramters to be used in pm3
        -# MNDO/d
        -# PDDG-PM3
    */
    void setHamiltonian(std::string ham) {
      hamiltonian = ham;
    }

    /*!
       - Geometry optimizations
        -# steep: steepest descent
        -# conjgrad: conjgrate gradient
        -# bfgs: 
        -# lbfgs: limited memory BFGS
    */
    void setOptimizer(std::string opt) {
      bOpt = true;
      optimizer = opt;
    }

    /*!
       - geometry optimizations
        -# maxopt: do a maximum of X cycles of geometry optimization
        -# snapgeom: print out every optimization step into separate files
        -# lsearch: linear search in lbfgs
        -# diagterm: initial diagonal term in lbfgs opt
        -# force-it: force geometry opt criteria to be satisfied
    */
    void setMaxOpt(int i) {
      bMaxOpt = true;
      maxopt = i;
    }

    //void setSnapGeom(int i) snapGeom = i;
    // void setLSearch();
    // void setDiagTerm();
    // void setForceIt();

    /*!
       - geometry optimizations
        -# ts: transition state geometry calculation
        -# tsprfo: baker method
        -# tsqna: quasi-newton method
        -# tsnr: pure newton method
        -# nomodefollow: mode following not activated in ts search
    */
    //void setTransitionStateOptimization(std::string ts) tsOpt = ts;

    /*!
       - Gradient
        -# gradient: output final gradient
        -# central: use central difference in gradient calculation
        -# inter: include only intermolecular contributions
        -# recipintra: consider intramolecular contributions
    void setGradient();
    void setCentral();
    void setInter();
    void setRecipintra();
    */

    /*!
       - Convergence criteria
        -# descf: user defined scf energy criterion [descf == ecrit]
        -# dpscf: user defined scf density criterion [dpscf == dcrit]
        -# etest: energy criterion
        -# gtest: gradient criterion
        -# ecrit: energy criterion
        -# dcrit: density matrix criterion
    
    void setDESCF(double d);
    void setDPSCF(double d);
    void setETEST(double d);
    void setGTEST(double d);
    void setECRIT(double d);
    void setDCRIT(double d);
    */

    /*!
       - Charge Models Available
        -# mulliken: use mulliken for scrf/pme calculations
        -# cm1: use mulliken for scrf/pme calculations
        -# cm2: use mulliken for scrf/pme calculations
    */
    void setChargeModel(std::string chg) {
      bChargeModel = true;
      chargeModel = chg;
    }

    /*!
       freq: carry out frequency analysis
    */
    void setFrequency() {
      bFreq = true;
    }

    /*!
       thermo: report thermodynamic quantities from freq calculations
    */
    void setThermo() {
      bThermo = true;
    }

     /*!
       - Subsetting
        -# cluster: cluster bases subsetting
        -# dbuff1:  buffer 1 size
        -# dbuff2:  buffer 2 size
        -# no-overlap: off diagonal terms of density matrix are set to zero
        -# atgrid: atom wise grid based subsetting
        -# resgrid: residue wise grid based subsetting
        -# mixgrid: residue wise, grid based subsetting for core 
                    and atom wise for grid based subsetting for buffers
        -# combsub: combination subsetting
        -# standard: no divide-and-conquer

     void setAtGrid();
     void setResGrid();
     void setMixGrid();
     void setCombSub();
     */

    void setCluster(int i) {
      if (i) {
        bCluster = true;
        bStandard = false;
      }
      else {
        bCluster = false;
      }
    }

    void setDBuff1(double d) {
      setCluster(1);
      dBuff1 = d;
    }

    void setDBuff2(double d) {
      setCluster(1);
      dBuff2 = d;
    }

    void setNoOverlap() {
      bNoOverlap = true;
    }

    void setStandard() {
      bStandard = true;
      bCluster = false;
    }

    /*!
       - Self Consistent Reaction Field:
        -# scrf: requests a scrf calculation with Delphi
        -# scrf-inside: only fill cavities with a continuum
        -# drad: used with scrf-inside
        -# indi: internal dielectric set to X
        -# exdi: external dielectric set to X
        -# scale: grid space set to X grid per angstrom
        -# perfill: molecule occupies X % of the grid
        -# ion: ion strength set to X
        -# probrab: probe radius set to X
        -# water: use water params for scrf calculation (Gnp)
        -# octanol: use octanol params for scrf calculation (Gnp)
        -# divpb: use divpb instead of delphi for scrf calculation
        -# soldmx: use guess dmx for scrf calculation
    */
    void doScrf() {
      bScrf = true;
    }

    void setScrfScale(double i) {
      scrfScale = i;
    }

    /*!
      direct: do not store 2-e integrals in scf calculations
    */
    void setDirect(int i) {
      if (i) {
        bDirect = true;
      }
      else {
        bDirect = false;
      }
    }

    /*!
      residue: store residue pointers (needed for D&C)
    */
    void setResidue(int i) {
      if (i) {
        bResidue = true;
      }
      else {
        bResidue = false;
      }
    }

    /*!
      pdump: write a restart density matrix every X scf cycles
    */
    void setPDump(int i) {
      bPDump = true;
      pDump = i;
    }

    /*!
      dump: write a restart file every X opt cycles
    */
    void setDump(int i) {
      bDump = true;
      dump = i;
    }

    /*!
      shift: initial dynamic level shift set to X
    */
    void setShift(double d) {
      if (d > -1) {
        bShift = true;
        shift = d;
      }
      else {
        bShift = false;
      }
    }

    /*!
      trajectory: dump coordinates to a trajectory file at restart points
    */
    //void setTrajectory(int i);

    /*!
      cutbond: cutoff bonding beyond X Angstroms
    */
    void setCutBond(double d) {
      bCutBond = true;
      cutbond = d;
    }

    /*!
      1SCF: perform only one scf iteration
    */
    //void set1Scf();

    /*!
      guess: build initial density matrix from one or more dmx files
    */
    void setGuess() {
      bGuess = true;
    }

    /*!
      vdw: calculate van der Waals interactions
    */
    //void setVdw();

    /*!
       - nmr
        -# nmr: nmr shielding calculation
    */
    void setNMR() {
      bNMR = true;
    }

    /*!
       - nmr
        -# calnuc = 1 proton calculation
        -# calnuc = 2 carbon calculation
    */
    void setCalNuc(unsigned int i) {
      calnum = i;
      this->setNMR();
    }

    /*!
       - PairWise Energy Decomposition
        -# pwd: perform a pwd
    */
    void setPWD() {
      bPwd = true;
    }

    /*!
       - PairWise Energy Decomposition
        -# pwd_atom: perform a pwd and output atom info only
    */
    void setPWDatom() {
      bPwdAtom = true;
    }

    /*!
       - PairWise Energy Decomposition
        -# pwd_residue: perform a pwd and output residue only
    */
    void setPWDresidue() {
      bPwdResidue = true;
    }

    /*!
       - Output
        -# prtsub: prints subsystem atom list
    */
    void setPrtSub() {
      bPrtSub = true;
    }

    /*!
       - Output
        -# prtvec: prints final eigenvectors
    */
    void setPrtVec() {
      bPrtVec = true;
    }

    /*!
       - Output
        -# prtcoords: prints atomic coordinates
    */
    void setPrtCoords() {
      bPrtCoords = true;
    }

    /*!
       - Output
        -# prtpar: prints semi-empirical parameters
    */
    void setPrtPar() {
      bPrtPar = true;
    }

    /*!
       - Output
        -# prtvdwaals: output van der Waals radii
    */
    void setPrtVdw() {
      bPrtVdw = true;
    }

    /*!
      addmm: adds MM correction to peptide rotational barrier [default]
    */
    void setAddMM() {
      bAddMM = true;
    }

    /*!
      nomm: does not use addmm
    */
    void setNoMM() {
      bNoMM = true;
    }

    /*!
       - Coordinate System
        -# cartesian: requires x,y,z data
        -# internal: requires a z-matrix
    */
    void setCoord(std::string coordSystem) {
      coord = coordSystem;
    }

    /*!
      fullscf: performs full diagonalizations in SCF
    */
    void setFullSCF() {
      bFullSCF = true;
    }

    /*!
      maxit: maximum number of scf iterations is 100
    */
    void setMaxIt(int i) {
      bMaxIt = true;
      maxIt = i;
    }

    /*!
      double: do a double scf step for every scf iteration
    */
    void setDouble(int d) {
      bDouble = true;
      iDouble = d;
    }

    /*!
      intgls: stewart/talman to be used [talman]
    */
    void setIntegrals(std::string s) {
      bIntegrals = true;
      integrals = s;
    }

    /*!
       tempk: set temperature to X
    */
    void setTempK(double k) {
      bTempK = true;
      tempK = k;
    }

    /*!
       tmax: request a maximum CPU time of X seconds
    */
    void setMaxTime(int s) {
      bMaxTime = true;
      maxTime = s;
    }

    /*!
       rmin: minimum distance between atoms is X 
    */
    void setMinR(double d) {
      bMinR = true;
      minR = d;
    }

    /*!
       cutrepul: cut off for repulsion integrals
    */
    void setCutRepul(double d) {
      bCutRepul = true;
      cutRepul = d;
    }

    /*!
       chkres: check inter-residue distances for each residue
    */
    void setChkRes() {
      bChkRes = true;
    }

    /*!
       dos: density of states
    */
    void setDos() {
      bDOS = true;
    }

    /*!
       dipole: calculates dipole moment
    */
    void setDipole() {
      bDipole = true;
    }

    /*!
       ip: calculate ionization potential
    */
    void setIP() {
      bIP = true;
    }

    /*!
       homolumo: calculates homo-lumo gap
    */
    void setHomoLumo() {
      bHomoLumo = true;
    }

    /*!
       screen: output info to screen
    */
    void setScreen() {
      bScreen = true;
    }

    /*!
       zmake: generates a z-matrix
    */
    void setZMake() {
      bZmake = true;
    }

    /*!
       qmailgn: generate output for QM-QSAR
    */
    void setQMAlign() {
      bQMAlign = true;
      this->setIntegrals("stewart");
    }

    /*!
       imult: sets multiplicity
    */
    void setImult(int i) {
      imult = i;
    }

    /*!
       uhf: sets uhf
    */
    void setUhf(int i) {
      if (i) {
        bUhf = true;
      }
      else {
        bUhf = false;
      }
    }

    /*!
       External: sets multiplicity
    */
    void setExternal(int i) {
      if (i) {
        bExternal = true;
      }
      else {
        bExternal = false;
      }
    }

    /*!
       rotate: rotational angle for barrier is X
    */
    //void setRotate();

    /*!
       push: push cluster groups apart
    */
    //void setPush();

     /*!
       \brief Write DivCon input file
       \param i Input file
       \param c collection pointer
       \param success success boolean
    */
    void Write(const std::string &i, collection* c, bool &success);

     /*!
       \brief Write DivCon input file
       \param i Input file
       \param m molecule pointer
       \param success success boolean
    */
    void Write(const std::string &i, molecule* m, bool &success);

     /*!
       \brief Write DivCon input file
       \param i Input file
       \param m molecule pointer
       \param coordinates Vector of coordinates which overrides the molecules coordinates
       \param success success boolean
    */
    void Write(const std::string &i, molecule* m,
         std::vector< vector3d > &coordinates, bool &success);

protected: // functions

    /*!
       \brief Read coordinate section of a divcon output
       \param idc input file stream
       \param m molecule pointer
    */
    void readCoordinates(std::ifstream &idc, molecule* m);

    /*!
       \brief Read qmqsar section of a divcon output
       \param idc input file stream
       \param m molecule pointer
    */
    void readQMQSAR(std::ifstream &idc, molecule* m);

    /*!
       \brief Read nmr section of a divcon output
       \param idc input file stream
       \param m molecule pointer
    */
    void readNMR(std::ifstream &idc, molecule* m);

    /*!
       \brief Do some error checking before writing the file
    */
    bool check();

    /*!
       \brief Open Input File
    */
    void openInputFile(const std::string &i);

    /*!
       \brief Write Head of a DivCon Input File
    */
    void writeHead();

    /*!
       \brief Write Coordinate section of a DivCon Input File
    */
    void writeCoords();

    /*!
       \brief Write Coordinate section of a DivCon Input File
    */
    void writeCoords(std::vector< vector3d > &coordinates);

    /*!
       \brief Write Tail of a DivCon Input File
    */
    void writeTail();

protected: // data
 
    //! collection pointer
    collection*        pCollection;

    //! molecule pointer
    molecule*          pMolecule;

    //! submolecule pointer
    submolecule*       pSubMolecule;

    //! atom pointer
    atom*              pAtom;

    //! molecule iterator
    typedef std::vector<molecule*>::iterator molIterator;

    //! submolecule iterator
    typedef std::vector<submolecule*>::iterator sMolIterator;

    //! atom iterator
    typedef std::vector<atom*>::iterator AtomIterator;

    //! string iterator
    typedef std::vector<std::string>::iterator strIterator;

     //! divcon input file name
    std::string        inputFileName;

     //! divcon output file name
    std::string        outputFileName;

    //! Input File Stream
    std::ofstream                 inputFileStream;

    //! cartesion/internal
    std::string coord;

    //! standard (no DnC)
    bool bStandard;

    //! DnC calculation
    bool bCluster;

    //! size of first buffer in DnC
    double dBuff1;

    //! size of second buffer in DnC
    double dBuff2;

    //! store residue pointers
    bool bResidue;

    //! calculate integrals on the fly
    bool bDirect;

    //! Which Hamiltonian to use
    std::string hamiltonian;

    //! Turn on a charge model
    bool bChargeModel;

    //! Which Charge model to use
    std::string chargeModel;

    //! Use a cutoff for repulsion integrals
    bool bCutBond;

    //! size of cutoff
    double cutbond;

    //! use shift (helps convergence)
    bool bShift;

    //! size of shift
    double shift;

    //! Use a guess file
    bool bGuess;

    //! output density matrix file
    bool bPDump;

    //! 
    int pDump;

    //! output restart coordinate file
    bool bDump;

    //! output restart coordinate file every n step
    int dump;

    //! double scf
    bool bDouble;

    //! 
    int iDouble;

    //! maximum number of optimization steps
    bool bMaxOpt;

    //! maximum number of optimization steps
    unsigned int maxopt;

    //! optimize structure
    bool bOpt;

    //! which optimizer
    std::string optimizer;

    //! run solvent calculation
    bool bScrf;

    //! run solvent calculation
    double scrfScale;

    //! solvent is water
    bool bWater;

    //! solvent is octanol
    bool bOctanol;

    //! 
    bool bNoOverlap;

    //! Output screen info
    bool bScreen;

    //! 
    bool bVdw;

    //! calculate dipoles
    bool bDipole;

    //! calculate frequencies
    bool bFreq;

    //! calculate thermochemical properties
    bool bThermo;

    //! calculate NMR shieldings
    bool bNMR;

    //! which nmr nuclei to calculate
    int calnum;

    //! do energy decomposition
    bool bPwd;

    //! do energy decomposition at atomic level
    bool bPwdAtom;

    //! do energy decomposition at residue level
    bool bPwdResidue;

    //! print subsystem info
    bool bPrtSub;

    //! print eigenvector info
    bool bPrtVec;

    //! print coordinates
    bool bPrtCoords;

    //! print parameter info
    bool bPrtPar;

    //! print van der Waals info
    bool bPrtVdw;

    //! add MM correction [default]
    bool bAddMM;

    //! no MM correction
    bool bNoMM;

    //! 
    bool bFullSCF;

    //! 
    bool bMaxIt;

    //! 
    int  maxIt;

    //! 
    bool bCutRepul;

    //! 
    double cutRepul;

    //! 
    bool bIntegrals;

    //! 
    std::string integrals;

    //! 
    bool bTempK;

    //! 
    double tempK;

    //! 
    bool bMaxTime;

    //! 
    int maxTime;

    //! 
    bool bMinR;

    //! 
    double minR;

    //! 
    bool bQMAlign;

    //! 
    bool bChkRes;

    //! 
    bool bDOS;

    //! 
    bool bIP;

    //! 
    bool bHomoLumo;

    //! 
    bool bZmake;

    //! 
    int imult;

    //! 
    bool bUhf;

    //! 
    bool bExternal;
};

} // MTKpp namespace

#endif // DCPARSER_H

