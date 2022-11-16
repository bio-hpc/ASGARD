/*! 
   \file amberAngle.h
   \brief AMBER angle energy and gradient
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

#ifndef AMBERANGLE_H
#define AMBERANGLE_H

namespace MTKpp
{

class amber;

// ============================================================
// Class : amberAngle()
// ------------------------------------------------------------
/*! 
   \class amberAngle

   \brief AMBER angle energy and gradient

   \author Martin Peters

   \version 0.1

   \date 2006

   \section AMBERAngleEnG AMBER Angle Energy And Gradient

   \f{eqnarray*}
    E_{\rm angle} & = & K_\theta (\theta - \theta_{eq})^2
   \label{AMBER:angle}
   \f}
    where \f$K_\theta\f$ is the angle force constant, \f$\theta\f$ is the angle (\f$0 \leq \theta \leq \pi\f$), 
    \f$\theta_{eq}\f$ is the standard angle.

   \f{eqnarray*}
    \theta       & = & \arccos\left({{r_{ij} \cdot r_{kj}} \over {|r_{ij}||r_{kj}|}}\right) \\
    \cos{\theta} & = & {{r_{ij} \cdot r_{kj}} \over {|r_{ij}||r_{kj}|}}
   \label{cos:theta}
   \f}

   \f{eqnarray*}
    \nabla_iE & = & {\partial E \over \partial\theta} \cdot 
                    {\partial\theta \over \partial\cos\theta} \cdot 
                    \nabla_i \cos \theta \\
              & = &  2K_\theta(\theta - \theta_{eq}) \cdot -\left({1 \over \sin\theta}\right) \cdot \nabla_i\cos\theta \\
   \label{dAngle:E}
   \f}

   How to determine \f$\nabla_i\cos\theta\f$?

   Considering that \f$\cos\theta\f$ is a function of \f$r_{ij}\f$ and \f$r_{kj}\f$ both of which are functions of \f$r_i\f$,
   \f$r_j\f$, and \f$r_k\f$.  Therefore you need to use the chain rule:

   \f{eqnarray*}
    {\partial\cos\theta \over \partial x_i} & = & \hat x \nabla_i\cos\theta \\
          & = &  {\partial\cos\theta \over \partial x_{ij}} {\partial x_{ij} \over \partial x_i} +
                 {\partial\cos\theta \over \partial y_{ij}} {\partial y_{ij} \over \partial x_i} +
                 {\partial\cos\theta \over \partial z_{ij}} {\partial z_{ij} \over \partial x_i} + \nonumber \\
          &   &  {\partial\cos\theta \over \partial x_{kj}} {\partial x_{kj} \over \partial x_i} +
                 {\partial\cos\theta \over \partial y_{kj}} {\partial y_{kj} \over \partial x_i} +
                 {\partial\cos\theta \over \partial z_{kj}} {\partial z_{kj} \over \partial x_i} \\
          & = &  {\partial\cos\theta \over \partial x_{ij}}
   \label{dcost:dxi}
   \f}

   There is a total of 9 such expression similar to the previous equation which lead to the following

   \f{eqnarray*}
   {\nabla_i\cos\theta} & = & \hat x {\partial\cos\theta \over \partial x_{ij}} +
                           \hat y {\partial\cos\theta \over \partial y_{ij}} +
                           \hat z {\partial\cos\theta \over \partial z_{ij}} \\
   {\nabla_j\cos\theta} & = & \hat x \left({\partial\cos\theta \over \partial x_{ij}} +
                                        {\partial\cos\theta \over \partial x_{kj}}\right) + 
                           \hat y \left({\partial\cos\theta \over \partial y_{ij}} +
                                        {\partial\cos\theta \over \partial y_{kj}}\right) + \nonumber \\
                        &   & \hat z \left({\partial\cos\theta \over \partial z_{ij}} +
                                        {\partial\cos\theta \over \partial z_{kj}}\right) \\
   {\nabla_k\cos\theta} & = & \hat x {\partial\cos\theta \over \partial x_{kj}} +
                           \hat y {\partial\cos\theta \over \partial y_{kj}} +
                           \hat z {\partial\cos\theta \over \partial z_{kj}}
   \label{nine:angle}
   \f}

   Finally, how do you determine \f$d\cos\theta/dx_{ij}\f$?

   \f{eqnarray*}
   \cos{\theta}                & = & {{r_{ij} \cdot r_{kj}} \over {|r_{ij}||r_{kj}|}}
   \f}

   \f{eqnarray*}
    {\partial\cos\theta \over \partial x_{ij}} & = & { |r_{ij}||r_{kj}| {d\left(r_{ij}\cdot {r_{kj}} \right) \over dx_{ij}} - \left(r_{ij}\cdot{r_{kj}}\right) {d\left(|r_{ij}||r_{kj}|  \right) \over dx_{ij}} \over {|r_{ij}|^2|r_{kj}|^2} } \nonumber \\
                            & = & {|r_{ij}||r_{kj}|r_{kj} - \left(r_{ij} \cdot r_{kj}\right){r_{ij} \over |r_{ij}|}|r_{kj}| \over {|r_{ij}|^2|r_{kj}|^2} } \nonumber \\
                            & = & {r_{kj} \over |r_{ij}||r_{kj}|} - {\left( r_{ij}\cdot r_{kj}\right) {r_{ij} \over |r_{ij}|} \over |r_{ij}|^2|r_{kj}|} \nonumber \\
                            & = & {1 \over |r_{ij}|} \left[{r_{kj} \over |r_{kj}|} - {\left(r_{ij}\cdot r_{kj}\right) \over |r_{ij}||r_{kj}|} {r_{ij} \over |r_{ij}|} \right] \nonumber \\
                            & = & {1 \over |r_{ij}|} \left[\hat r_{kj} - \cos\theta\hat r_{ij} \right] \\
    {\partial\cos\theta \over \partial x_{kj}} & = & {1 \over |r_{kj}|} \left[\hat r_{ij} - \cos\theta\hat r_{kj} \right]
   \label{dcos:dxs}
   \f}
*/
// ============================================================
class amberAngle
{
public:
    /*!
       \brief amberAngle Constructor
    */
    amberAngle();

    /*!
       \brief amberAngle Constructor
       \param pAmber parent potential
    */
    amberAngle(amber* pAmber);

    /*!
       \brief amberAngle Destructor
    */
    virtual ~amberAngle();

    /*!
       \brief Calculate angle Energy
       \return angle energy
    */
    double calculateE();

    /*!
       \brief Calculate angle Energy and Gradient
       \return angle energy
    */
    double calculateG();

protected: // FUNCTIONS
    /*!
       \brief Determine the angle between atoms 1, 2, and 3
       \param x1 atom 1 x coordinate
       \param y1 atom 1 y coordinate
       \param z1 atom 1 z coordinate
       \param x2 atom 2 x coordinate
       \param y2 atom 2 y coordinate
       \param z2 atom 2 z coordinate
       \param x3 atom 3 x coordinate
       \param y3 atom 3 y coordinate
       \param z3 atom 3 z coordinate
       \param bError error boolean
       \return angle between atoms 1, 2, and 3
    */
    double angle(double x1, double y1, double z1,
                 double x2, double y2, double z2,
                 double x3, double y3, double z3,
                 bool &bError);

protected: // DATA
    //! parent potential
    amber*         pAmber;

    // bond energy
    double         energy;

    //! angle
    double         deltaAngle;

    //!
    bool           bError;
};

} // MTKpp namespace

#endif // AMBERANGLE_H
