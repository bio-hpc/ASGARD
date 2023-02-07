/*! 
   \file vector3d.h
   \brief A 3-dimensional vector with many common
          functions needed for the manipulation
          of 3d vectors.  Everything is inline so
          as to speed things up.

   \author Andrew Wollacott
   \author Martin Peters

   Andrew Wollacott (primary author)

   Martin Peters
    - Added rotation functions:
    -# set
    -# rotateX
    -# rotateY
    -# rotateZ
    -# formRotMat
    -# rotateXYZ
    -# buildCoord

   $Date: 2010/03/29 20:33:22 $
   $Revision: 1.9 $

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


#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <math.h>
#include <cmath>
#include <iostream>

#include "constants.h"

namespace MTKpp
{

// ============================================================
// Class : vector3d()
// ------------------------------------------------------------
/*! 
   \class vector3d
   \brief A 3-dimensional vector class
   \author Andrew Wollacott
   \author Martin Peters
   \version 0.1
   \date 2005
*/
// ============================================================
class vector3d
{
public:

   /*!
       \brief vector3d Constructor

       Initializes coordinates to 0,0,0
   */
   inline vector3d() 
          : X(0.0), Y(0.0), Z(0.0) {};

   /*!
       \brief vector3d Constructor
       \param x x coordinate
       \param y y coordinate
       \param z z coordinate
   */
   inline vector3d(double x, double y, double z)
          : X(x), Y(y), Z(z) {};

   /*!
       \brief vector3d Constructor
       \param xyz Sets all three coordinates to this value
   */
   inline vector3d(double xyz)
          : X(xyz), Y(xyz), Z(xyz) {};

   /*!
       \brief Set x, y, and z coordinates
       \param x x value
       \param y y value
       \param z z value
   */
   inline void set(double x, double y, double z) {
          this->X = x;
          this->Y = y;
          this->Z = z;}

   /*!
       \brief Set x coordinate to this value
       \param x x value
   */
   inline void setX(double x) {
          this->X = x;}

   /*!
       \brief Set y coordinate to this value
       \param y y value
   */
   inline void setY(double y) {
          this->Y = y;}

   /*!
       \brief Set z coordinate to this value
       \param z z value
   */
   inline void setZ(double z) {
          this->Z = z;}

   /*!
       \brief Get x coordinate
       \return x value
   */
   inline double getX() {
          return this->X;}

   /*!
       \brief Get y coordinate
       \return y value
   */
   inline double getY() {
          return this->Y;}

   /*!
       \brief Get z coordinate
       \return z value
   */
   inline double getZ() {
          return this->Z;}

   /*!
       \brief Access coordinate as array
       \param i Coordinate to be returned
       \return Coordinate i
   */
   inline double &operator[](int i) {
          return  i == 0 ? X
                 :i == 1 ? Y
                 :i == 2 ? Z
                 : X; }

   /*!
       \brief Set vector equal to constant
       \param value Value to be set
       \return vector3d
   */
   inline vector3d& operator=(const double &value) {
            X = value;
            Y = value;
            Z = value;
            return *this; }

   /*!
       \brief Set negative of the vector
       \param lhs negative of this vector
       \return vector3d
   */
   inline friend vector3d operator-(const vector3d &lhs) {
            return vector3d(-lhs.X, -lhs.Y, -lhs.Z); }

   /*!
       \brief Add another vector3d object to this
       \param rhs vector3d object
   */
   inline void operator+=(const vector3d &rhs) {
            X += rhs.X;
            Y += rhs.Y;
            Z += rhs.Z; }

   /*!
       \brief Subtract another vector3d object to this
       \param rhs vector3d object
   */
   inline void operator-=(const vector3d &rhs) {
            X -= rhs.X;
            Y -= rhs.Y;
            Z -= rhs.Z; }

   /*!
       \brief Multiply vector3d by a constant
       \param value Value thats multiplied
   */
   inline void operator*=(const double &value) {
            X *= value;
            Y *= value;
            Z *= value; }

   /*!
       \brief Divide vector3d by a constant
       \param value Value thats divided
   */
   inline void operator/=(const double &value) {
            X /= value;
            Y /= value;
            Z /= value; }

   /*!
       \brief Checks equality of two vectors
       \param lhs Left hand side vector3d object
       \param rhs Right hand side vector3d object
   */
   inline friend bool operator==(const vector3d &lhs, const vector3d &rhs) {
            return (lhs.X == rhs.X) && (lhs.Y == rhs.Y) && (lhs.Z == rhs.Z); }

   /*!
       \brief Checks equality of two vectors
       \param lhs Left hand side vector3d object
       \param rhs Right hand side vector3d object
   */
   inline friend bool operator!=(const vector3d &lhs, const vector3d &rhs) {
            return (lhs.X != rhs.X) || (lhs.Y != rhs.Y) || (lhs.Z != rhs.Z); }

   /*!
       \brief Add two vectors together
       \param lhs Left hand side vector3d object
       \param rhs Right hand side vector3d object
   */
   inline friend vector3d operator+(const vector3d &lhs, const vector3d &rhs) {
            return vector3d( (lhs.X + rhs.X), (lhs.Y + rhs.Y), (lhs.Z + rhs.Z) ); }

   /*!
       \brief Subtract two vectors
       \param lhs Left hand side vector3d object
       \param rhs Right hand side vector3d object
   */
   inline friend vector3d operator-(const vector3d &lhs, const vector3d &rhs) {
            return vector3d( (lhs.X - rhs.X), (lhs.Y - rhs.Y), (lhs.Z - rhs.Z) ); }

   /*!
       \brief Calculate the dot-product of two vectors
       \param lhs Left hand side vector3d object
       \param rhs Right hand side vector3d object
       \return dot-product
   */
   inline friend double operator*(const vector3d &lhs, const vector3d &rhs) {
            return lhs.X*rhs.X + lhs.Y*rhs.Y + lhs.Z*rhs.Z; }

   /*!
       \brief Add a vector to a scalar
       \param lhs Left hand side vector3d object
       \param scalar Scalar to be added
       \return vector3d object
   */
   inline friend vector3d operator+(const vector3d &lhs, const double &scalar) {
            return vector3d( lhs.X + scalar, lhs.Y + scalar, lhs.Z + scalar); }

   /*!
       \brief Add a scalar to a vector
       \param scalar Scalar to be added
       \param rhs Right hand side vector3d object
       \return vector3d object
   */
   inline friend vector3d operator+(const double &scalar, const vector3d &rhs) {
            return vector3d( rhs.X + scalar, rhs.Y + scalar, rhs.Z + scalar); }

   /*!
       \brief Subtract a scalar from a vector
       \param lhs Left hand side vector3d object
       \param scalar Scalar to be subtracted
       \return vector3d object
   */
   inline friend vector3d operator-(const vector3d &lhs, const double &scalar) {
            return vector3d( lhs.X - scalar, lhs.Y - scalar, lhs.Z - scalar); }

   /*!
       \brief Subtract a vector from a scalar
       \param scalar Scalar to be subtracted
       \param rhs Right hand side vector3d object
       \return vector3d object
   */
   inline friend vector3d operator-(const double &scalar, const vector3d &rhs) {
            return vector3d( scalar - rhs.X, scalar - rhs.Y, scalar - rhs.Z); }

   /*!
       \brief Return the product of a vector with a scalar
       \param lhs Left hand side vector3d object
       \param scalar Scalar to be multiplied
       \return vector3d object
   */
   inline friend vector3d operator*(const vector3d &lhs, const double &scalar) {
            return vector3d( scalar*lhs.X, scalar*lhs.Y, scalar*lhs.Z ); }

   /*!
       \brief Return the product of a scalar with a vector
       \param scalar Scalar to be multiplied
       \param rhs Right hand side vector3d object
       \return vector3d object
   */
   inline friend vector3d operator*(const double &scalar, const vector3d &rhs) {
            return vector3d( scalar*rhs.X, scalar*rhs.Y, scalar*rhs.Z) ; }

   /*!
       \brief Divide a vector with a scalar
       \param lhs Left hand side vector3d object
       \param scalar Scalar to be multiplied
       \return vector3d object
   */
   inline friend vector3d operator/(const vector3d &lhs, const double &scalar) {
            return vector3d( lhs.X/scalar, lhs.Y/scalar, lhs.Z/scalar ); }

   /*!
       \brief Output
       \param os stream
       \param v vector3d object
   */
   inline friend std::ostream& operator<< (std::ostream& os, vector3d& v) {
      os << v.X << " " << v.Y << " " << v.Z;
      return os; }

   /*!
       \brief Return the magnitude of a vector
       \return length of vector3d
   */
   inline double length() {
            return sqrt(X*X + Y*Y + Z*Z); }

   /*!
       \brief Returns the square of the magnitude
       \return length of vector3d squared
   */
   inline double length_squared() {
            return (X*X + Y*Y + Z*Z); }

   /*!
       \brief Returns the unit vector
       \return vector3d object
   */
   inline vector3d unit() {
             return vector3d(X, Y, Z)/length(); }

   /*!
       \brief Return the cross product of two vectors
       \return vector3d object
   */
   inline friend vector3d cross(const vector3d &lhs, const vector3d &rhs) {
            return vector3d( (lhs.Y*rhs.Z - rhs.Y*lhs.Z),
                             (lhs.Z*rhs.X - rhs.Z*lhs.X),
                             (lhs.X*rhs.Y - rhs.X*lhs.Y) ); }

   /*!
       \brief Return the "distance" between two vectors
       \return distance
   */
   inline double dist(const vector3d &v2) {
            return sqrt( (X-v2.X)*(X-v2.X) + (Y-v2.Y)*(Y-v2.Y) + (Z-v2.Z)*(Z-v2.Z) ); }

   /*!
       \brief Return the square of the "distance" between two vectors
       \return distance squared
   */
   inline double dist2(const vector3d &v2) {
            return ( (X-v2.X)*(X-v2.X) + (Y-v2.Y)*(Y-v2.Y) + (Z-v2.Z)*(Z-v2.Z) ); }

   /*!
       \brief Return the angle between three vectors
       \param a vector3d a
       \param b vector3d b
       \param c vector3d c
       \return angle in radians
   */
   inline friend double angle(const vector3d &a, const vector3d &b,
                              const vector3d &c) {
            vector3d ab = a - b;
            vector3d cb = c - b;

            ab = ab.unit();
            cb = cb.unit();

            return acos(ab*cb); }

   /*!
       \brief Return the torsion between four vectors
       \param a vector3d a
       \param b vector3d b
       \param c vector3d c
       \param d vector3d d
       \return angle in radians

       \f[
          \cos{\phi}  & = & {{t \cdot u} \over {|t||u|}}
       \f]
       where:
       \f{eqnarray*}
          r_{ab} & = & r_{a} - r_{b} \\
          r_{cb} & = & r_{c} - r_b \\
          r_{cd} & = & r_c - r_d \\
          t      & = & r_{ab} \times r_{cb} \\
          u      & = & r_{cb} \times r_{cd} \nonumber
       \f}

   */
   inline friend double torsion(const vector3d& a, const vector3d& b,
                                const vector3d& c, const vector3d& d) {
            vector3d ab = a - b;
            vector3d cb = c - b;
            vector3d cd = c - d;

            vector3d p,q;
            p = cross(ab, cb);
            q = cross(cb, cd);

            p = p.unit();
            q = q.unit();

            double pXq = p*q;
            if (std::abs(pXq) > 1) {
              pXq = 1;
            }
            else if (pXq < -1) {
              pXq = -1;
            }

            // acos calculates the arc cosine
            // Return value in range [0, PI].
            double ang = acos(pXq);
            double s   = cb*(cross(p, q));

            if (s < 0) ang = -ang;

            return ang;   }

   /*!
       \brief Return the torsion between four vectors
       \param a vector3d a
       \param b vector3d b
       \param c vector3d c
       \param d vector3d d
       \return angle in radians

       \f[
          \phi & = & { \arctan2(|b2|b1 \cdot [b2 \times b3], [b1 \times b1] \cdot [b2 \times b3] )}
       \f]
       where:
       \f{eqnarray*}

          b1 & = & r_a - r_b \\
          b2 & = & r_c - r_b \\
          b3 & = & r_c - r_d \nonumber
       \f}
   */
   inline friend double torsion2(const vector3d& a, const vector3d& b,
                                 const vector3d& c, const vector3d& d) {
            vector3d b1 = a - b;
            vector3d b2 = c - b;
            vector3d b3 = c - d;

            double b2Length = b2.length();
            vector3d v1 = b2Length * b1;

            vector3d b2Xb3 = cross(b2, b3);
            vector3d b1Xb2 = cross(b1, b2);

            double lhs = v1 * b2Xb3;
            double rhs = b1Xb2 * b2Xb3;

            // atan2 calculates the arc tangent of lhs/rhs in radians
            // To compute the value, the function uses the sign of both 
            //   arguments to determine the quadrant.
            // Return value in range [-PI, PI].
            double ang = atan2(lhs, rhs);

            return ang;   }

   /*!
       \brief Rotate about the X-axis
       \param ang Angle to be rotated by
   */
   inline void rotateX(const double &ang) {
            double cosAng = cos(ang);
            double sinAng = sin(ang);
            double Yt     = Y;
            double Zt     = Z;

            Y = Yt*cosAng - Zt*sinAng;
            Z = Yt*sinAng + Zt*cosAng;  }

   /*!
       \brief Rotate about the X-axis, with origin o
       \param ang Angle to be rotated by
       \param o origin
   */
   inline void rotateX(const double &ang, const vector3d &o) {
            double dY = Y - o.Y;
            double dZ = Z - o.Z;
            double cosAng = cos(ang);
            double sinAng = sin(ang);
            double Yt     = Y;
            double Zt     = Z;           

            Yt = dY*cosAng - dZ*sinAng;
            Zt = dY*sinAng + dZ*cosAng;

            Y = Yt + o.Y;
            Z = Zt + o.Z;  }

   /*!
       \brief Rotate about the Y-axis
       \param ang Angle to be rotated by
   */
   inline void rotateY(const double &ang) {
            double cosAng = cos(ang);
            double sinAng = sin(ang);
            double Zt     = Z;
            double Xt     = X;
            
            X = Zt*sinAng + Xt*cosAng;
            Z = Zt*cosAng - Xt*sinAng;  }

   /*!
       \brief Rotate about the Y-axis, with origin o 
       \param ang Angle to be rotated by
       \param o origin
   */
   inline void rotateY(const double &ang, const vector3d &o) {
            double dX = X - o.X;
            double dZ = Z - o.Z;
            double cosAng = cos(ang);
            double sinAng = sin(ang);
            double Zt     = Z;
            double Xt     = X;

            Xt = dZ*sinAng + dX*cosAng;
            Zt = dZ*cosAng - dX*sinAng;
            X = Xt + o.X;
            Z = Zt + o.Z;  }

   /*!
       \brief Rotate about the Z-axis
       \param ang Angle to be rotated by
   */
   inline void rotateZ(const double &ang) {
            double cosAng = cos(ang);
            double sinAng = sin(ang);
            double Xt     = X;
            double Yt     = Y;
            
            X = Xt*cosAng - Yt*sinAng;
            Y = Xt*sinAng + Yt*cosAng;  }

   /*!
       \brief Rotate about the Z-axis, with origin o
       \param ang Angle to be rotated by
       \param o origin
   */
   inline void rotateZ(const double &ang, const vector3d &o) {
            double dX = X - o.X;
            double dY = Y - o.Y;
            double cosAng = cos(ang);
            double sinAng = sin(ang);
            double Xt     = X;
            double Yt     = Y;

            Xt = dX*cosAng - dY*sinAng;
            Yt = dX*sinAng + dY*cosAng;

            X = Xt + o.X;
            Y = Yt + o.Y;  }

   /*!
       \brief Form rotation matrix R from a torsion and rotational angle in radians

       This is axis angle rotation

       \param a Coordinate 1
       \param b Coordinate 2
       \param c Coordinate 3
       \param d Coordinate 4
       \param angle Angle of rotation
       \param rotMat Rotation matrix which gets returned
       \return successful or not
       \code

         1) Measure torsion, t, a-b-c-d

         2) Make sure t lies between 0 and 2*PI

         3) Calculate rotAngle = angle - t

         4) Calculate cos(rotAngle) and sin(rotAngle), and Normalize the vector b - c

         5) Calculate rotation matrix, R:

               +-                               -+
               | t*x*x+c    t*x*y+s*z  t*x*z-s*y |
         R =   | t*x*y-s*z  t*y*y+c    t*y*z+s*x |
               | t*x*z+s*y  t*y*z-s*x  t*z*z+c   |
               +-                               -+

         where:
           c = cos(angle)
           s = sin(angle)
           t = 1 - c
           x,y,z = normalized b-c x,y,z coordinates

       \endcode
   */
   inline friend int formRotMat(const vector3d& a, const vector3d& b,
                                const vector3d& c, const vector3d& d,
                                double angle, double rotMat[9]) {

            // Measure Torsion, t
            double torsionAngle = torsion(a,b,c,d);
            if (torsionAngle < 0) torsionAngle = torsionAngle + 2*PI;

            // Calculate Rotational Angle
            double rotAngle = angle - torsionAngle;

            double cosAng = 0;
            double sinAng = 0;
            double t = 0;
            vector3d bc;

            if (std::abs(rotAngle) > 0.0001) {
              cosAng = cos(rotAngle);
              sinAng = sin(rotAngle);
              t = 1.0 - cosAng;

              // Normalize the vector b - c
              bc = b - c;
              bc = bc.unit();
            }
            else {
              return 0;
            }

            // Form R
            rotMat[0] = t*bc.X*bc.X + cosAng;
            rotMat[1] = t*bc.X*bc.Y + sinAng*bc.Z;
            rotMat[2] = t*bc.X*bc.Z - sinAng*bc.Y;

            rotMat[3] = t*bc.X*bc.Y - sinAng*bc.Z;
            rotMat[4] = t*bc.Y*bc.Y + cosAng;
            rotMat[5] = t*bc.Y*bc.Z + sinAng*bc.X;

            rotMat[6] = t*bc.X*bc.Z + sinAng*bc.Y;
            rotMat[7] = t*bc.Y*bc.Z - sinAng*bc.X;
            rotMat[8] = t*bc.Z*bc.Z + cosAng; 
            return 1;}

   /*!
       \brief Rotate xyz according to the rotation matrix R with origin o
       \param R rotation matrix
       \param o origin
   */
   inline void rotateXYZ(const double R[9], const vector3d& o) {
            this->X -= o.X;
            this->Y -= o.Y;
            this->Z -= o.Z;
            double newX = this->X*R[0] + this->Y*R[1] + this->Z*R[2];
            double newY = this->X*R[3] + this->Y*R[4] + this->Z*R[5];
            double newZ = this->X*R[6] + this->Y*R[7] + this->Z*R[8];
            this->X = newX + o.X;
            this->Y = newY + o.Y;
            this->Z = newZ + o.Z;}

   /*!
      \brief Generate coordinate a from b, c, d, dist, angle, and torsion
      \param a coordinate to be generated
      \param b bonded coordinate
      \param c angle coordinate
      \param d torsion coordinate
      \param dist distance of a-b (angstrom)
      \param angle angle a-b-c (radians)
      \param torsion torsion a-b-c-d (radians)
      \code
          b -- c
         /      \
        a        d
      \endcode
   */
   inline friend void buildCoord(vector3d& a, const vector3d& b,
                                 const vector3d& c, const vector3d& d,
                                 double dist, double angle, double torsion) {
            // Normalize the vector b - c 
            vector3d bc = b - c;
            bc = bc.unit();

            // Normalize the vector d - c 
            vector3d dc = d - c;
            dc = dc.unit();

            vector3d r = cross(bc, dc);
            r = r.unit();

            vector3d s = cross(r, bc);
            s = s.unit();

            a = b - dist * cos(angle) * bc + 
                    dist * sin(angle) * sin(torsion) * r +
                    dist * sin(angle) * cos(torsion) * s;}

protected:

   //! X coordinate
   double X;

   //! Y coordinate
   double Y;

   //! Z coordinate
   double Z;
};

} // MTKpp namespace

#endif
