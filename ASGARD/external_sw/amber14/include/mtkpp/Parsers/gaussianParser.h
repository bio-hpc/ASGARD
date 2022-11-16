/*!
   \file gaussianParser.h
   \brief Parses gaussian files
   \author Martin Peters

   Reads and writes guassian files

   $Date: 2010/07/22 20:41:34 $
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

#ifndef GAUSSIANPARSER_H
#define GAUSSIANPARSER_H

#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include "baseParser.h"

#include "Statistics/sheet.h"
#include "Statistics/table.h"

namespace MTKpp
{

class collection;
class molecule;
class submolecule;
class atom;
class element;
class connections;
struct Bond;
struct Angle;
struct Torsion;
struct Improper;
class vector3d;
class zmatParser;

// ============================================================
// Class : gaussianParser()
// ------------------------------------------------------------
/*! 
   \class gaussianParser
   \brief Reads and writes gaussian format files
   \author Martin Peters
   \version 0.1
   \date 2006
*/
// ============================================================
class gaussianParser : public baseParser

{
public:

    /*!
       \brief gaussianParser Constructor
    */
    gaussianParser();

    //! gaussianParser Destructor
    ~gaussianParser();

    /*!
       \brief Read gaussian formatted file
       \param i Input file
       \param pSheet sheet pointer
    */
    void           Read(const std::string &i, sheet* pSheet = 0);

    /*!
       \brief Write gaussian formatted input file
       \param o Output file
       \param m molecule pointer
       \todo Write function
    */
    void           Write(const std::string &o, molecule* m);

    /*!
       \brief Write gaussian formatted input file
       \param o Output file
       \param m molecule pointer
       \param coordinates Vector of coordinates which overrides the molecules coordinates
       \todo Write function
    */
    void           Write(const std::string &o, molecule* m, std::vector< vector3d > &coordinates);

    /*!
       \brief Write gaussian formatted input file
       \param o Output file
       \param c collection pointer
       \param molID The ID of the molecule in the collection
       \todo Write function
    */
    void           Write(const std::string &o, collection* c, const int &molID);

    /*!
       \brief Generate Z-Matrix
       \param m molecule pointer
    */
    int            generateZMatrix(molecule* m);

    /*!
       \brief Read zmatrix mapping file
       \param f mapping file name
    */
    void           readMappingFile(std::string f);

    /*!
       \brief Write zmatrix mapping file
       \param m mapping file name
    */
    void           writeMappingFile(std::string m);

    /*!
       \brief Add command
       \param c command name
    */
    void           addCommandOption(const std::string &c);

    /*!
       \brief Add command and value
       \param c command name
       \param o command value
    */
    void           addCommandOption(const std::string &c, const std::string &o);

    /*!
       \brief Add command and vector of values
       \param c command name
       \param o command value
    */
    void           addCommandOption(const std::string &c, std::vector<std::string> &o);

    /*!
       \brief Remove command and value
       \param c command name
    */
    std::vector<std::string> getCommandOption(const std::string &c);

    /*!
       \brief Remove command and value
       \param c command name
    */
    void           removeCommandOption(const std::string &c);

    /*!
       \brief Add iop
       \param i iop command
    */
    void           addIop(const std::string &i);

    /*!
       \brief clear iop
    */
    void           clearIop();

    /*!
       \brief Set check point file
       \param c chkpt file
    */
    void           setChkPt(const std::string &c);

    /*!
       \brief Set memory requirement
       \param m memory
    */
    void           setMem(const std::string &m);

    /*!
       \brief Set number of processors to use
       \param n number of processors
    */
    void           setNProc(const std::string &n);

    /*!
       \brief Set level of theory
       \param t theory
    */
    void           setTheory(const std::string &t);

    /*!
       \brief Set basis set
       \param b basis set
    */
    void           setBasisSet(const std::string &b);

    /*!
       \brief Set basis set file
       \param b basis set file
    */
    void           setBasisSetFile(const std::string &b);

    /*!
       \brief Set pseudopotential file
       \param b pseudopotential file
    */
    void           setPseudoPotentialFile(const std::string &b);

    /*!
       \brief Set modredundant file
       \param b modredundant file
    */
    void           setModRedundantFile(const std::string &b);

    /*!
       \brief Set charge
       \param i charge
    */
    void           setCharge(const int &i);

    /*!
       \brief Set multiplicity
       \param i multiplicity
    */
    void           setMultiplicity(const int &i);

    /*!
       \brief Set verbosity
       \param i verbosity
    */
    void           setVerbosity(const std::string &i);

    /*!
       \brief Write internal coordinates
       \param i on/off switch
    */
    void           setInternal(int i);

    /*!
       \brief Write Cartesian coordinates
       \param i on/off switch
    */
    void           setCartesian(int i);

    /*!
       \brief Do not write coordinates
    */
    void           setNoCoords();

    /*!
       \brief Write molecule name
       \param i switch 
    */
    void           setWriteMoleculeName(int i);

    /*!
       \brief Write charge and multiplicity
       \param i switch 
    */
    void           setWriteChargeAndMult(int i);

    /*!
       \brief Set Merz-Kollman radii
       \param e element symbol
       \param d value
    */
    void           setMKRadii(std::string e, double d);

    /*!
       \brief Read Formatted Check Point File
       \param fchkPtFile formatted check point file
       \param pSheet sheet pointer
    */
    void           readFormattedChkPtFile(std::string fchkPtFile, sheet* pSheet);

    /*!
       \brief Get Force Constant
       \param pSheet sheet pointer
       \param i atom i index
       \param j atom j index
       \param r bond distance
       \param fc force constant
       \return success
    */
    int            getForceConstant(sheet* pSheet, int i, int j, double& r, double& fc);

    /*!
       \brief Get Force Constant
       \param pSheet sheet pointer
       \param i atom i index
       \param j atom j index
       \param k atom k index
       \param a angle
       \param fc force constant
       \return success
    */
    int            getForceConstant(sheet* pSheet, int i, int j, int k, double& a, double& fc);

    /*!
       \brief Get Force Constant
       \param i atom i index
       \param j atom j index
       \param r bond distance
       \param fc force constant
       \return success
    */
    int            getForceConstantZMAT(int i, int j, double& r, double& fc);

    /*!
       \brief Get Force Constant
       \param i atom i index
       \param j atom j index
       \param k atom k index
       \param a angle
       \param fc force constant
       \return success
    */
    int            getForceConstantZMAT(int i, int j, int k, double& a, double& fc);

    /*!
       \brief Get Frequencies
    */
    std::vector<double> getFrequencies();

protected:
    //! chk pt file
    std::string    itsChkPtFile;

    //! Write check point file
    bool           bChkPt;

    //! memory requirement
    std::string    itsMem;

    //! Turn on mem
    bool           bMem;

    //! number of processors
    std::string    itsNProc;

    //! Turn on number of processors
    bool           bNProc;

    //! Use internal or cartesian coordinates
    bool           bWriteInternalCoords;

    //! Use internal or cartesian coordinates
    bool           bWriteCartCoords;

    //! Write molecule name
    bool           bWriteMoleculeName;

    //! Write charge and multiplicity
    bool           bWriteChargeAndMult;

    //! Level of theory
    std::string    itsTheory;

    //! Basis Set
    std::string    itsBasisSet;

    //! Basis Set File
    std::string    itsBasisSetFile;

    //! PseudoPotential File
    std::string    itsPseudoPotentialFile;

    //! modredundant File
    std::string    itsModRedundantFile;

    //! Charge
    int            itsCharge;

    //! Multiplicity
    int            itsMultiplicity;

    //! Logging level (verbosity) - T, N or P?
    std::string    itsVerbosity;

    //! command options map
    std::map<std::string, std::vector<std::string> > itsCommandOptions;

    //! Merz-Kollman radii
    std::map<std::string, double> itsMKRadii;

    //! command options map iterator
    typedef std::map<std::string, std::vector<std::string> >::iterator mapIterator;

    //! MK radii iterator
    typedef std::map<std::string, double>::iterator dMapIterator;

    //!
    std::vector<std::string> iops;

    // Z-MATRIX
    //! zmatParser pointer
    zmatParser* pZmatParser;

    //! z-matrix was generated
    bool bZMatrixGenerated;

    //! zmatrix
    std::vector<std::vector<std::string> >  zmatrix;

    //! zmatrix bonds
    std::vector<std::vector<int> >          zmatrixBonds;

    //! zmatrix angles
    std::vector<std::vector<int> >          zmatrixAngles;

    //! zmatrix torsions
    std::vector<std::vector<int> >          zmatrixTorsions;

    //! 
    std::map<std::string, double>           zmatData;

    //! 
    std::vector<std::string>                zmatDataAngles;

    //! 
    std::vector<std::string>                zmatDataTorsions;

    //! 
    std::vector<std::string>                zmatDataNames;

    //! 
    std::map<std::string, int>              zmatDataType;

    //! Force constant map
    std::map<std::string, double>           forceConstants;

    //! Frequencies
    std::vector<double>                     frequencies;

    //! mapping between molecule and zmatrix
    std::map<int, int>                      atomMap;

    //! int, int iterator
    typedef std::map<int, int>::iterator    iMapIterator;

    //! Redundant coordinates found in output file
    bool                                    bModRedundant;

    bool                                    g09corhigh;

    //! redundant bonds
    std::map<int, std::string>              modRedBonds;

    //! redundant angles
    std::map<ULONG_KIND, std::string>       modRedAngles;

    //! redundant dihedrals
    std::map<int, std::string>              modRedDihedrals;

    //! labels
    std::vector<std::string>                atomLabels;

    //! coordinates
    std::vector<vector3d*>                  molCoords;

    //! Eigenvectors
    //ublas::matrix<double, ublas::column_major> eigenvectors;
    Eigen::Matrix<double, Dynamic, Dynamic> eigenvectors;

    //! Atom labels
    table<std::string>*                     pGaussAtoms;

    //! Coordinates
    table<double>*                          pGaussCoords;

    //! Eigenvalue table
    table<double>*                          pGaussEValues;

    //! Eigenvector table
    table<double>*                          pGaussEVectors;
};

} // MTKpp namespace

#endif // GAUSSIANPARSER_H
