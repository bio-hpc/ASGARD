/*!
   \file amber.h
   \brief AMBER MM potential
   \author Martin Peters

   $Date: 2010/03/29 20:28:34 $
   $Revision: 1.7 $

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

#ifndef AMBER_H
#define AMBER_H

#include "mmPotential.h"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace MTKpp
{

class amberBond;
class amberAngle;
class amberTorsion;
class amberImproper;
class amberNonBonded;

// ============================================================
// Class : amber()
// ------------------------------------------------------------
/*! 
   \class amber

   \brief AMBER MM potentials

   \author Martin Peters

   \version 0.1

   \date 2005

   \section amberEF AMBER Energy Function

   \f{eqnarray*}
    E_{\rm total} & = & \sum_{\rm bonds} K_r \left(r - r_{eq}\right)^2 +
                        \sum_{\rm angles} K_\theta \left(\theta - \theta_{eq}\right)^2 + \nonumber \\
                  &   &  \sum_{\rm dihedrals} {V_n \over 2} \left[1 + \cos\left(n\phi - \gamma\right)\right] +
                        \sum_{i<j} \left[{A_{ij} \over r_{ij}^{12}} - {B_{ij} \over r_{ij}^6} + 
                                         {q_iq_j \over \varepsilon r_{ij}} \right] \nonumber
   \label{AMBER:FF}
   \f}

    The AMBER energy function contains bond, angle, dihedral, and non-bonded terms.
    The bond and angle terms are represented by diagonal harmonic expressions.
    The van der Waals term is a 6-12 potential, and the electrostatic is expressed as a 
    Coulombic interaction with atom centered point charges.
    \f$K_r (kcal/mol \AA^2)\f$, \f$K_\theta (kcal/(mol radian^{2}))\f$, \f$V_n\f$ is the magnitude of the torsion,
    \f$n\f$ is the periodicity, \f$\gamma\f$ is the phase difference, \f$A_{ij} = \varepsilon r_{ij}^{*12}\f$ and 
    \f$B_{ij} = 2\varepsilon r_{ij}^{*6}\f$, \f$r_{ij} = r_{i}^{*} + r_{j}^{*}\f$ in \f$\AA\f$, \f$r\f$ is the van der Waals 
    radius for atom \f$i\f$, and \f$\varepsilon_{ij} = \sqrt{\varepsilon_i * \varepsilon_j}\f$

    The uncompressed energy function:

   \f{eqnarray*}
    E_{\rm total} & = & \sum_{\rm bonds} K_r (r_{ij} - r_{eq})^2 +
                        \sum_{\rm angles} K_\theta (\theta_{ijk} - \theta_{eq})^2 +
                        \sum_{\rm dihedrals} \sum_{\rm n} V_n[1 + {\rm cos}(n\phi_{ijkl} - \gamma)] + \nonumber \\
                  &   & \sum_{i<j} \varepsilon_{ij} \left[ \left({r_{ij}^{*} \over r_{ij}}\right)^{12} -
                                  2\left({r_{ij}^{*} \over r_{ij}}\right)^6 \right] +
                       {1 \over VDW } \sum_{i<j}^{1-4} \varepsilon_{ij} \left[ \left({r_{ij}^{*} \over r_{ij}} \right)^{12} -
                                  2\left({r_{ij}^{*} \over r_{ij}}\right)^6 \right] + \nonumber \\
                  &   & \sum_{i<j} {q_iq_j \over \varepsilon r_{ij}} + {1 \over ELE} \sum_{i<j}^{1-4} {q_iq_j \over 
                        \varepsilon r_{ij}} \nonumber
   \label{AMBER:FF2}
   \f}
    where \f$VDW == ELE == 1.2\f$.
*/
// ============================================================
class amber : public mmPotential
{
    friend class amberNonBonded;
public:
    /*!
       \brief amber Constructor
    */
    amber();

    /*!
       \brief amber Destructor
    */
    virtual ~amber();

    /*!
       \brief Initialize amber
    */
    void amberInitialize();

    /*!
       \brief Calculate bond energy
    */
    virtual double calcBondEnergy();

    /*!
       \brief Calculate angle energy
    */
    virtual double calcAngleEnergy();

    /*!
       \brief Calculate torsion energy
    */
    virtual double calcTorsionEnergy();

    /*!
       \brief Calculate improper energy
    */
    virtual double calcImproperEnergy();

    /*!
       \brief Calculate Non bonded energy
    */
    virtual double calcNonBondedEnergy();

    /*!
       \brief Calculate H-Bond Energy
    */
    virtual double calcHBondEnergy();

    /*!
       \brief Set R^6 L-J energy
    */
    void setR6(double);

    /*!
       \brief Set R^12 L-J energy
    */
    void setR12(double);

    /*!
       \brief Set 1-4 van der Waals energy
    */
    void setVDW14(double);

    /*!
       \brief Set van der Waals energy
    */
    void setVDW(double);

    /*!
       \brief Set 1-4 Electrostatic energy
    */
    void setEle14(double);

    /*!
       \brief Set Electrostatic energy
    */
    void setEle(double);

    /*
       \brief print pairwise decomposition data
    */
    virtual void printPWD(int i, int j, double r6, double r12, double vdw_energy, double electrostatic, double r);

protected: // DATA
    //! amberBond pointer
    amberBond*          pAmberBond;

    //! amberAngle pointer
    amberAngle*         pAmberAngle;

    //! amberTorsion pointer
    amberTorsion*       pAmberTorsion;

    //!
    //amberImproper*      pAmberImproper;

    //! amberNonBonded pointer
    amberNonBonded*     pAmberNonBonded;

};

} // MTKpp namespace
#endif // AMBER_H

