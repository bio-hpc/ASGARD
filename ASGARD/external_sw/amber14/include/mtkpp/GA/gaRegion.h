/*! 
   \file gaRegion.h
   \brief Container for gaPopulations
   \author Martin Peters

   $Date: 2010/03/29 20:24:52 $
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

#ifndef GAREGION_H
#define GAREGION_H

#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "Utils/constants.h"

namespace MTKpp
{

class gaWorld;
class gaPopulation;
class gaIndividual; //temp
class gaGene;
class gaOperators;
class gaOutput;

// ============================================================
// Class : gaRegion()
// ------------------------------------------------------------
/*! 
   \class gaRegion
   \brief Container for gaPopulations
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================

class gaRegion
{
    friend class gaPopulation;
    friend class gaGene;
    friend class gaOperators;
    friend class gaOutput;
public:

    /*!
       \brief gaRegion Constructor
       \param parent gaWorld pointer
    */
    gaRegion(gaWorld *parent = 0);

    //! gaRegion Destructor
    virtual ~gaRegion();

    /*!
       \brief Setup gaRegion
       \param maxInds Maximum number of gaIndividuals in the gaRegion
       \param seed Random number generator seed
       \param maxGens Maximum number of generations to be carried out
       \param elitism Elitism to be employed or not
       \param selectionPressure Parameter used in the selection of gaIndividuals
       \param curPop Current gaPopulation index
       \param popKeep Number of gaPopulations to be stored
       \param nChild Number of children to be produced per pair of gaIndividuals
       \param chreDiff ????
       \param pKeep Percentage of the current gaPopulation to be retained for the next generation
       \param pCrossover Percentage of the next generation to be created by crossover
       \param pMutate Percentage of the next generation to be created by mutation
       \param pAverage Percentage of the next generation to be created by averaging
       \param selection Selection method to be used
       \param crossover Crossover method to be used
       \param mutate Mutation method to be used 
       \param average Averaging method to be used
       \param funcDir Either minimize or maximize the fitness function
       \param maxParameters Maximum parameters of the genetic information
       \param minParameters Minimum parameters of the genetic information
       \param stepSize Step size for the genetic information
    */
    void setup(const int& maxInds, const int& seed, const int& maxGens, 
               const int& elitism, const double& selectionPressure, 
               const int& curPop, const int& popKeep, const int& nChild, const double& chreDiff,
               const double& pKeep, const double& pCrossover, 
               const double& pMutate, const double& pAverage,
               const std::string& selection, const std::string& crossover,
               const std::string& mutate, const std::string& average,
               const std::string& funcDir,
               std::vector<double> maxParameters, std::vector<double> minParameters,
               std::vector<double> stepSize);

    /*!
       \brief Add gaPopulation to gaRegion
       \return gaPopulation pointer
    */
    gaPopulation* addPopulation();

    /*!
       \brief Delete gaPopulation from gaRegion
       \param pop gaPopulation pointer
    */
    void delPopulation(gaPopulation* pop);

    /*!
       \brief Get gaPopulation from gaRegion
       \param number gaPopulation number, either id or index
       \param id bool return by id
       \param index bool return by index
       \return gaPopulation pointer
    */
    gaPopulation* getPopulation(const int& number, bool id, bool index);

    /*!
       \brief Get getPopulations in the gaRegion
       \return vector of gaPopulation pointers
    */
    std::vector<gaPopulation*> getPopulations();

    /*!
       \brief Get number of gaPopulations in gaRegion
       \return number of gaPopulations in gaRegion
    */
    int getNumPopulations();

    /*!
       \brief Initialize gaRegion
    */
    void initialize();

    /*!
       \brief Form the next generation in the gaRegion
    */
    void nextGeneration();

    //------------//
    // - get/set -//
    //------------//

    /*!
       \brief Set gaRegion id
       \param id gaRegion id
    */
    void setId(int id);

    /*!
       \brief Get gaRegion id
       \return gaRegion id
    */
    int getId();

    /*!
       \brief Set gaRegion name
       \param name gaRegion name
    */
    void setName(std::string name);

    /*!
       \brief Get gaRegion name
       \return gaRegion name
    */
    std::string getName();

    /*!
       \brief Get gaWorld which gaRegion is a member of
       \return gaWorld pointer
    */
    gaWorld* getParent();

    /*!
       \brief Set Current energy value
       \param d current fitness value
    */
    void setFitness(double d);

    /*!
       \brief Get Convergence data
       \return convergence data
    */
     std::vector<double> getFitness();

    /*!
       \brief Set min parameters
       \param m min parameter
    */
    void setMinParameters(std::vector<double> m);

    /*!
       \brief Get min parameter for certain gaGene bit
       \param chr gaChromosome index in gaIndividual
       \param gen gaGene index in gaChromosome
       \param pos bit index in gaGene
       \return minimum value
    */
    double getMinParameter(const int& chr, const int& gen, const int& pos);

    /*!
       \brief Set max parameters
       \param m max parameter
    */
    void setMaxParameters(std::vector<double> m);

    /*!
       \brief Get max parameter for certain gaGene bit
       \param chr gaChromosome index in gaIndividual
       \param gen gaGene index in gaChromosome
       \param pos bit index in gaGene
       \return maximum value
    */
    double getMaxParameter(const int& chr, const int& gen, const int& pos);

    /*!
       \brief Set step sizes
       \param s step sizes
    */
    void setStepSizes(std::vector<double> s);

    /*!
       \brief Get step size for certain gaGene bit
       \param chr gaChromosome index in gaIndividual
       \param gen gaGene index in gaChromosome
       \param pos bit index in gaGene
       \return minimum value
    */
    double getStepSize(const int& chr, const int& gen, const int& pos);

    /*!
       \brief Get max number of generations
       \return maximum number of generations
    */
    int getMaxGens();

    ////// WARNING: THE NEXT TWO FUNCTIONS ONLY WORK
    //////          FOR 1 CHR PER INDS WITH, N GENES IN THE CHR
    /*!
       \brief Set the possible gene values (works for 1 chr per inds with, N genes in the chr)
       \param gVs gene values
    */
    void setGeneValues(std::vector<std::vector<double> > gVs);

    /*!
       \brief Get possible values for gene (works for 1 chr per inds with, N genes in the chr)
       \param gen gaGene index in gaChromosome
       \return possible gene value
    */
    std::vector<double> getGeneValues(const int& gen);

protected:
    //! gaPopulation iterator
    typedef std::vector<gaPopulation*>::iterator populationIterator;

    //! gaPopulation vector
    std::vector<gaPopulation*>    itsPopulationList;

    //! gaWorld pointer
    gaWorld*                      pParent;

    //! gaPopulation pointer
    gaPopulation*                 pGaPopulation;

    //! gaPopulation pointer
    gaPopulation*                 pPreviousGaPop;

    /*!
       \brief gaOperator pointer
       \sa gaOperator for more details
    */
    gaOperators*                  pGaOperators;

    //-----------------------------//
    // - Options from Input file - //
    //-----------------------------//

    //! Maximum number of gaIndividuals in gaRegion
    int                           maxInds;

    //! Random number generation seed
    int                           seed;

    //! Maximum number of generations to be carried out
    int                           maxGens;

    //! Elitism is carried out or not
    int                           elitism;

    /*!
       \brief Parameter used in the selection of gaIndividuals

        Should be between 0 and 1. The lower the number means that you are
        more likely to pick a fitter gaIndividual.

       \sa gaSelection gaGaussian
    */
    double                        selectionPressure;

    //! Current gaPopulation index
    int                           curPop;

    //! Number of children allowed to be produced per pair of gaIndividual
    int                           nChild;

    //! Chromosome or Genes must differ by at least this value
    double                        chreDiff;

    //! Number of gaPopulation to be stored
    int                           popKeep;

    /*!
       \brief Percentage of the current gaPopulation to be retained for the next generation
       \sa gaOperator
    */
    double                        pKeep;

    /*!
       \brief Percentage of the next generation to be created by crossover
       \sa gaOperator
    */
    double                        pCrossover;

    /*!
       \brief Percentage of the next generation to be created by mutation
       \sa gaOperator
    */
    double                        pMutate;

    /*!
       \brief Percentage of the next generation to be created by averaging
       \sa gaOperator
    */
    double                        pAverage;

    /*!
       \brief Selection method to be used
       \sa gaSelection
    */
    std::string                   selection;

    /*!
       \brief Crossover method to be used
       \sa gaCrossover
    */
    std::string                   crossover;

    /*!
       \brief Mutation method to be used
       \sa gaMutate
    */
    std::string                   mutate;

    /*!
       \brief Averaging method to be used
       \sa gaAverage
    */
    std::string                   average;

    /*!
       \brief Minimize or maximize the fitness function
    */
    std::string                   funcDir;

    //! Maximum parameters
    std::vector<double>           maxParameters;

    //! Minimum parameters
    std::vector<double>           minParameters;

    //! Step sizes
    std::vector<double>           stepSize;

    //! possible gene values
    std::vector<std::vector<double> > geneValues;

    //-----------------------------//
    // - Options from world file - //
    //-----------------------------//

    //! gaRegion name
    std::string                   itsName;

    //! gaRegion id
    int                           itsId;

    //! Convergence tracker of gaRegion
    std::vector<double>           popConv;

};

} // MTKpp namespace

#endif // GAREGION_H
