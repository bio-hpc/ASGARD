/*!
   \file seqAlign.h
   \brief Sequence Alignment
   \author Martin Peters

   Literature/WWW
   - http://en.wikipedia.org/wiki/Dynamic_programming
   - http://en.wikipedia.org/wiki/Sequence_alignment
   - http://en.wikipedia.org/wiki/Needleman-Wunsch
   - http://en.wikipedia.org/wiki/Smith-Waterman
   - http://en.wikipedia.org/wiki/Gap_penalty

   $Date: 2010/03/29 20:45:26 $
   $Revision: 1.5 $

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

#ifndef SEQALIGN_H
#define SEQALIGN_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

#include "Utils/constants.h"

namespace MTKpp
{

class molecule;
class submolecule;
class atom;

// ============================================================
// Class : seqAlign()
// ------------------------------------------------------------
/*!
   \class seqAlign
   \brief Sequence Alignment
*/
// ============================================================

class seqAlign
{
public:

    //! seqAlign Constructor
    seqAlign();

    //! seqAlign Destructor.
    virtual ~seqAlign();

    /*!
       \brief Align Query onto Template
       \return success
    */
    int run();

    /*!
       \brief Set the template molecule
       \param pMolecule molecule pointer
       \return success
    */
    int setTemplate(molecule* pMolecule);

    /*!
       \brief Set the query molecule
       \param pMolecule molecule pointer
       \return success
    */
    int setQuery(molecule* pMolecule);

    /*!
       \brief Set the PAM matrix to be used
       \param pam PAM matrix file

       PAM: Point Accepted Mutation
    */
    void setPAMFile(std::string pam);

    /*!
       \brief Set the PAM matrix size
       \param s size
    */
    void setPAMSize(int s);

    /*!
       \brief Set a PAM type
       \param i row index
       \param s pam type
    */
    void setPAMType(int i, std::string s);

    /*!
       \brief Set the value of the PAM matrix
       \param i row index
       \param j column index
       \param v value
    */
    void setPAMValue(int i, int j, double v);

    /*!
       \brief Print PAM matrix
    */
    void printPAM();

    /*!
       \brief Set the Gap Penalty Type
       \param i type
    */
    void setAlgorithmType(int i);

    /*!
       \brief Set the Gap Penalty Type
       \param i type
    */
    void setGapPenaltyType(int i);

    /*!
       \brief Set the Gap-Open Parameter
       \param d gap-open parameter
    */
    void setGapOpen(double d);

    /*!
       \brief Set the Gap-Entend Parameter
       \param d gap-extend parameter
    */
    void setGapExtend(double d);

    /*!
       \brief Get the correspondence map
       \return correspondence map
    */
    int* getCorrMap();

    /*!
       \brief Print the results
    */
    void printResults();

    /*!
       \brief Get the PAM Score given two characters
       \param i first residue
       \param j second residue
       \return pam score
    */
    double getSimScore(std::string i, std::string j);

protected: // functions

    /*!
       \brief Initial Match Matrix
       \return success
    */
    int initMatchMatrix();

    /*!
       \brief Trace back the Match Matrix
       \return success
    */
    int traceBack();

protected: // data

    //! Template molecule
    molecule* tMol;

    //! Number of residues in template
    int tNRes;

    //! Template character sequence
    char* tCSeq;

    //! Template int sequence
    int* tISeq;

    //! Query molecule
    molecule* qMol;

    //! Number of residues in template
    int qNRes;

    //! Query character sequence
    char* qCSeq;

    //! Template int sequence
    int* qISeq;

    /*!
       \brief Algorith type
       - 1 Needleman-Wunsch
       - 2 Smith-Waterman
    */
    int algorithmType;

    /*!
       \brief gap penalty type
        g = gap size
        d = gapOpen
        e = gapExtend

       - 1 Constant gap penalty (gapPenalty = d)
       - 2 Linear gap penalty   (gapPenalty = -g*d)
       - 3 Affine gap penalty   (gapPenalty = -d - (g-1) * e, where e < d)
    */
    int gapPenaltyType;

    //! gap-open parameter
    double gapOpen;

    //! gap-extend parameter
    double gapExtend;

    //! maximum gap size
    int maxGap;

    //! current gap size
    int gapSize;

    //! PAM Matrix File
    std::string pamFile;

    //! PAM Size
    int pamSize;

    //! PAM Size
    char* pamTypes;

    //! PAM Matrix
    double* pam;

    //! Match Matrix
    double* matchMatrix;

    //! Index to max value in matchMatrix for current cell
    int* maxIndices;

    //! Row index of largest match matrix value
    int bestI;

    //! Column index of largest match matrix value
    int bestJ;

    //! Best Match Matrix value
    double bestMatchMatrixValue;

    //! Alignment of the template
    std::string alignmentA;

    //! Alignment of the query
    std::string alignmentB;

    //! Residues of the query mapped onto the residue of the template
    int* corrMap;
};

} // MTKpp namespace

#endif // COLLECTION_H
