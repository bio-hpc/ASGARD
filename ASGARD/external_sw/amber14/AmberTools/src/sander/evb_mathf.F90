! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  EVB math functions                                                     |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

   module evb_math

   implicit none

   contains

!  .............................................................................
!  :  Outer product of two vectors                                             :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   function outer ( u, v )

   intrinsic :: size
   intrinsic :: matmul, reshape 
   _REAL_ , intent(in), dimension(:) :: u, v
   _REAL_ , dimension(size(u),size(v)) :: outer

   outer = matmul( reshape( u, (/ size(u), 1 /) ) &
                 , reshape( v, (/ 1, size(v) /) ) )

   end function outer

!  .............................................................................
!  :  Cross product of two vectors                                             :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   function cross ( u, v )

   _REAL_ , intent(in) :: u(3), v(3)
   _REAL_ :: cross(3)

   cross(1) = u(2) * v(3) - u(3) * v(2)
   cross(2) = u(3) * v(1) - u(1) * v(3)
   cross(3) = u(1) * v(2) - u(2) * v(1)

   end function cross

!  .............................................................................
!  :  Vector from atom j to i                                                  :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   function vec ( v, i, j )

   integer, intent(in) :: i, j
   _REAL_ , intent(in), dimension(:,:) :: v
   _REAL_ :: vec(3)

   vec(:) = v(:,i) - v(:,j)

   end function vec

!  .............................................................................
!  :  Distance between atoms i & j                                             :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   function dist ( v, i, j )

   integer, intent(in) :: i, j
   _REAL_ , intent(in), dimension(:,:) :: v
   intrinsic :: sqrt, dot_product
   _REAL_ :: dist

!dEVB   dist = sqrt( ( v(1,i) - v(1,j) )**2 + ( v(2,i) - v(2,j) )**2 &
!dEVB                                       + ( v(3,i) - v(3,j) )**2 )

      dist = sqrt( dot_product( v(:,i) - v(:,j), v(:,i) - v(:,j) ) )

   end function dist

!  .............................................................................
!  :  Unit vector from atom j to i                                             :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   function uvec ( v, i, j )

   integer, intent(in) :: i, j
   _REAL_ , intent(in), dimension(:,:) :: v
   intrinsic :: sqrt, dot_product
   _REAL_ :: uvec(3)

   uvec(:) = v(:,i) - v(:,j)
   uvec(:) = uvec(:) / sqrt( dot_product( uvec(:), uvec(:) ) )

   end function uvec

!  .............................................................................
!  :  Vector rotation of i about axis j-k                                      :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   function vrot ( v, i, j, k )

   integer, intent(in) :: i, j, k
   _REAL_ , intent(in), dimension(:,:) :: v
   intrinsic :: dot_product
   _REAL_ :: vrot(3)

   vrot(:) = cross ( uvec(v,i,j), uvec(v,k,j) ) / ( dist(v,i,j) &
           * ( 1.0d0 - dot_product( uvec(v,i,j), uvec(v,j,k) )**2 ) )

   end function vrot

!  .............................................................................
!  :  Angle formed by atoms i, j, and k                                        :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   function angle ( v, i, j, k )

   integer, intent(in) :: i, j, k
   _REAL_ , intent(in), dimension(:,:) :: v
   intrinsic :: acos, dot_product
   _REAL_ :: angle

   angle = acos( dot_product( uvec(v,i,j), uvec(v,k,j) ) )

   end function angle

!  .............................................................................
!  :  Unit vector perpendicular to the plane of atoms i, j, and k              :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   function perp ( v, i, j, k )

   integer, intent(in) :: i, j, k
   _REAL_ , intent(in), dimension(:,:) :: v
   intrinsic :: sqrt, dot_product
   _REAL_ :: perp(3)

   perp(:) = cross ( vec(v,i,j), vec(v,k,j) )
   perp(:) = perp(:) / sqrt( dot_product( perp(:), perp(:) ) )

   end function perp

!  .............................................................................
!  :  Unit vector from j to i orthogonalized to j, k                           :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   function orthovec ( v, i, j, k )

   integer, intent(in) :: i, j, k
   _REAL_ , intent(in), dimension(:,:) :: v
   intrinsic :: sqrt, dot_product
   _REAL_ :: orthovec(3)

   orthovec(:) = ( uvec(v,i,j) - uvec(v,j,k) * dot_product( uvec(v,i,j) &
                 , uvec(v,j,k) ) ) / sqrt( 1.0d0 - dot_product( uvec(v,i,j) &
                 , uvec(v,j,k) )**2 )

   end function orthovec

!  .............................................................................
!  :  Dihedral angle formed by atoms i, j, k and l                             :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   function dihedral ( v, i, j, k, l )

   integer, intent(in) :: i, j, k, l
   _REAL_ , intent(in), dimension(:,:) :: v
   intrinsic :: acos, atan2, dot_product, abs
   _REAL_ :: dihedral, pi, perp_ijk(3), perp_jkl(3) 

   pi = acos( -1.0d0 )

   perp_ijk(:) = perp(v,i,j,k)
   perp_jkl(:) = perp(v,j,k,l)

   dihedral = atan2( dot_product( cross( perp_ijk(:),  uvec(v,j,k) ) &
            , perp_jkl(:) ), dot_product( perp_ijk(:), perp_jkl(:) ) )

   if( abs( abs(dihedral) - pi ) <= 0.0001d0 ) dihedral = pi 

   end function dihedral

!  .............................................................................
!  :  1-D indexing for a symmetric matrix stored in lower triangle format      :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   function lt_ndx ( i, j )

   integer, intent(in) :: i, j
   integer :: lt_ndx

   if( i > j ) then
      lt_ndx = i * ( i - 1 ) / 2 + j
   else
      lt_ndx = j * ( j - 1 ) / 2 + i
   endif

   end function lt_ndx

!  .............................................................................
!  :  Difference of internal coordinates, correcting for dihedral periodicity  :
!  :  (assumes internal are ordered as bonds, angle & dihedrals)               :
!  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

   function dqint ( q1, q2, nbond, nangle, ndihed, ncoord, ncart )

   integer , intent(in) :: nbond, nangle, ndihed, ncoord, ncart
   _REAL_  , intent(in) :: q1(ncoord), q2(ncoord)
   intrinsic :: acos, abs, sign
   _REAL_  :: dqint(ncoord), pi 
   integer :: n, nn

   pi = acos( -1.0d0 )

   nn = 0

   do n = 1, nbond
      nn = nn + 1
      dqint(nn) = q1(nn) - q2(nn)
   enddo

   do n = 1, nangle
      nn = nn + 1
      dqint(nn) = q1(nn) - q2(nn)
   enddo

   do n = 1, ndihed
      nn = nn + 1
      dqint(nn) = q1(nn) - q2(nn)
      if( abs( dqint(nn) ) >= pi ) &
         dqint(nn) = dqint(nn) - 2.0d0 * sign( pi, dqint(nn) )
   enddo

   do n = 1, ncart
      nn = nn + 1
      dqint(nn) = q1(nn) - q2(nn)
   enddo

  
   end function dqint

   end module evb_math

