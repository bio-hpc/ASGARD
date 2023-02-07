! A C-program for MT19937, with initialization improved 2002/1/26.
! Coded by Takuji Nishimura and Makoto Matsumoto.

! Code converted to Fortran 95 by Josi Rui Faustino de Sousa
! Date: 2002-02-01

! Before using, initialize the state by using init_genrand(seed)
! or init_by_array(init_key, key_length).

! This library is free software.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

! Copyright (C) 1997, 2002 Makoto Matsumoto and Takuji Nishimura.
! Any feedback is very welcome.
! email: matumoto@math.keio.ac.jp

!This is a slightly modified version of the code obtained here:
!http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/FORTRAN/mt19937ar.f90

! The license for this code is taken from this web site:
!  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/elicense.html
!  and is reproduced here:
!
!  Commercial Use of Mersenne Twister
!  2001/4/6
!
! Until 2001/4/6, MT had been distributed under GNU Public License, but after
! 2001/4/6, we decided to let MT be used for any purpose, including commercial
! use. 2002-versions mt19937ar.c, mt19937ar-cok.c are considered to be usable
! freely.

#include "ncsu-utils.h"
#include "ncsu-config.h"

#define ASSUME_GFORTRAN yes

module mt19937

#if defined(MPI) && defined(BINTRAJ)

  implicit none

  intrinsic :: bit_size

  private

  integer,  parameter  :: intg = selected_int_kind( 9 )
  integer,  parameter  :: long = selected_int_kind( 18 )
  integer,  parameter  :: flot = selected_real_kind( 6, 37 )
  integer,  parameter  :: dobl = selected_real_kind( 15, 307 )

  integer,  public, parameter :: wi = intg
  integer,  public, parameter :: wl = long
  integer,  public, parameter :: wr = dobl

  ! Period parameters
  integer( kind = wi ), parameter :: n = 624_wi
  integer( kind = wi ), parameter :: m = 397_wi
  integer( kind = wi ), parameter :: hbs = bit_size( n ) / 2_wi
  integer( kind = wi ), parameter :: qbs = hbs / 2_wi
  integer( kind = wi ), parameter :: tbs = 3_wi * qbs

  type, public :: mt19937_t
    private
    integer(kind = wi) :: mt(n) ! the array for the state vector
    logical(kind = wi) :: mtinit = .false._wi ! means mt[N] is not initialized
    integer(kind = wi) :: mti = n + 1_wi ! mti==N+1 means mt[N] is not initialized
  end type mt19937_t

  public :: init_by_seed
  public :: init_by_array

  public :: random_int32
  public :: random_int31

  public :: random_real1
  public :: random_real2
  public :: random_real3

  public :: random_res53

# ifndef NCSU_NO_NETCDF
  private :: check_ncrc
  public :: mt19937_save
  public :: mt19937_load
# endif /* NCSU_NO_NETCDF */

# ifdef MPI
  public :: mt19937_bcast
# endif /* MPI */

  contains

  elemental function uiadd( a, b ) result( c )

    implicit none

    intrinsic :: ibits, ior, ishft

    integer( kind = wi ), intent( in )  :: a, b

    integer( kind = wi )  :: c

    integer( kind = wi )  :: a1, a2, b1, b2, s1, s2

    a1 = ibits( a, 0, hbs )
    a2 = ibits( a, hbs, hbs )
    b1 = ibits( b, 0, hbs )
    b2 = ibits( b, hbs, hbs )
    s1 = a1 + b1
    s2 = a2 + b2 + ibits( s1, hbs, hbs )
    c  = ior( ishft( s2, hbs ), ibits( s1, 0, hbs ) )

  end function uiadd

  elemental function uisub( a, b ) result( c )

    implicit none

    intrinsic :: ibits, ior, ishft

    integer( kind = wi ), intent( in )  :: a, b

    integer( kind = wi )  :: c

    integer( kind = wi )  :: a1, a2, b1, b2, s1, s2

    a1 = ibits( a, 0, hbs )
    a2 = ibits( a, hbs, hbs )
    b1 = ibits( b, 0, hbs )
    b2 = ibits( b, hbs, hbs )
    s1 = a1 - b1
    s2 = a2 - b2 + ibits( s1, hbs, hbs )
    c  = ior( ishft( s2, hbs ), ibits( s1, 0, hbs ) )

  end function uisub

  elemental function uimlt( a, b ) result( c )

    implicit none

    intrinsic :: ibits, ior, ishft

    integer( kind = wi ), intent( in )  :: a, b

    integer( kind = wi )  :: c

    integer( kind = wi )  :: a0, a1, a2, a3
    integer( kind = wi )  :: b0, b1, b2, b3
    integer( kind = wi )  :: p0, p1, p2, p3

    a0 = ibits( a, 0, qbs )
    a1 = ibits( a, qbs, qbs )
    a2 = ibits( a, hbs, qbs )
    a3 = ibits( a, tbs, qbs )
    b0 = ibits( b, 0, qbs )
    b1 = ibits( b, qbs, qbs )
    b2 = ibits( b, hbs, qbs )
    b3 = ibits( b, tbs, qbs )
    p0 = a0 * b0
    p1 = a1 * b0 + a0 * b1 + ibits( p0, qbs, tbs )
    p2 = a2 * b0 + a1 * b1 + a0 * b2 + ibits( p1, qbs, tbs )
    p3 = a3 * b0 + a2 * b1 + a1 * b2 + a0 * b3 + ibits( p2, qbs, tbs )
    c  = ior( ishft( p1, qbs ), ibits( p0, 0, qbs ) )
    c  = ior( ishft( p2, hbs ), ibits( c, 0, hbs ) )
    c  = ior( ishft( p3, tbs ), ibits( c, 0, tbs ) )

  end function uimlt

  ! initializes mt[N] with a seed
  subroutine init_by_seed(self, s)

    implicit none

    intrinsic :: iand, ishft, ieor, ibits

    type(mt19937_t), intent(inout) :: self
    integer( kind = wi ), intent( in )  :: s

    integer( kind = wi )  :: i, mult_a

#ifndef ASSUME_GFORTRAN
    data mult_a /z'6C078965'/ ! gfortran does not like this
#else
    mult_a = ieor(ishft(z'6C07', 16), z'8965') ! but this is okay
#endif /* ASSUME_GFORTRAN */

    self%mtinit = .true._wi
    self%mt(1) = ibits( s, 0, 32 )
    do i = 2, n, 1
      self%mt(i) = ieor( self%mt(i-1), ishft(self%mt(i-1), -30 ) )
      self%mt(i) = uimlt( self%mt(i), mult_a)
      self%mt(i) = uiadd( self%mt(i), uisub( i, 1_wi ) )
      ! See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.
      ! In the previous versions, MSBs of the seed affect
      ! only MSBs of the array mt[].
      ! 2002/01/09 modified by Makoto Matsumoto
      self%mt(i) = ibits( self%mt(i), 0, 32 )
      ! for >32 bit machines
    end do
    self%mti = n + 1_wi

  end subroutine init_by_seed

  ! initialize by an array with array-length
  ! init_key is the array for initializing keys
  ! key_length is its length
  subroutine init_by_array(self, init_key)

    implicit none

    intrinsic :: iand, ishft, ieor

    type(mt19937_t), intent(inout) :: self
    integer( kind = wi ), intent( in )  :: init_key(:)

    integer( kind = wi )  :: i, j, k, tp, key_length
    integer( kind = wi )  :: seed_d, mult_a, mult_b, msb1_d

    data seed_d /z'12BD6AA'/
    data mult_a /z'19660D'/

#ifndef ASSUME_GFORTRAN
    data mult_b /z'5D588B65'/
    data msb1_d /z'80000000'/
#else
    mult_b = ieor(ishft(z'5D58', 16), z'8B65')
    msb1_d = 1
    msb1_d = ishft(msb1_d, 31)
#endif /* ASSUME_GFORTRAN */

    key_length = size( init_key, dim = 1 )
    call init_by_seed(self, seed_d)
    i = 2_wi
    j = 1_wi
    do k = max( n, key_length ), 1, -1
      tp = ieor( self%mt(i-1), ishft( self%mt(i-1), -30 ) )
      tp = uimlt( tp, mult_a )
      self%mt(i) = ieor( self%mt(i), tp )
      self%mt(i) = uiadd( self%mt(i), uiadd( init_key(j), uisub( j, 1_wi ) ) ) ! non linear
      self%mt(i) = ibits( self%mt(i), 0, 32 ) ! for WORDSIZE > 32 machines
      i = i + 1_wi
      j = j + 1_wi
      if ( i > n ) then
        self%mt(1) = self%mt(n)
        i = 2_wi
      end if
      if ( j > key_length) j = 1_wi
    end do
    do k = n-1, 1, -1
      tp = ieor( self%mt(i-1), ishft( self%mt(i-1), -30 ) )
      tp = uimlt( tp, mult_b )
      self%mt(i) = ieor( self%mt(i), tp )
      self%mt(i) = uisub( self%mt(i), uisub( i, 1_wi ) ) ! non linear
      self%mt(i) = ibits( self%mt(i), 0, 32 ) ! for WORDSIZE > 32 machines
      i = i + 1_wi
      if ( i > n ) then
        self%mt(1) = self%mt(n)
        i = 2_wi
      end if
    end do
    self%mt(1) = msb1_d ! MSB is 1; assuring non-zero initial array
  end subroutine init_by_array

  ! generates a random number on [0,0xffffffff]-interval
  function random_int32(self) result( y )

    implicit none

    type(mt19937_t), intent(inout) :: self

    intrinsic :: iand, ishft, ior, ieor, btest, ibset, mvbits

    integer( kind = wi )  :: y

    integer( kind = wi )  :: kk
    integer( kind = wi )  :: seed_d
    data seed_d   /z'5489'/

#ifndef ASSUME_GFORTRAN
    integer( kind = wi)   :: matrix_a, matrix_b, temper_a, temper_b

    data matrix_a /z'9908B0DF'/
    data matrix_b /z'0'/
    data temper_a /z'9D2C5680'/
    data temper_b /z'EFC60000'/
#else
#   define matrix_a z'9908B0DF'
#   define matrix_b z'0'
#   define temper_a z'9D2C5680'
#   define temper_b z'EFC60000'
#endif /* ASSUME_GFORTRAN */

    if ( self%mti > n ) then ! generate N words at one time
      if ( .not. self%mtinit ) call init_by_seed(self, seed_d) ! if init_genrand() has not been called, a default initial seed is used
      do kk = 1, n-m, 1
        y = ibits( self%mt(kk+1), 0, 31 )
        call mvbits( self%mt(kk), 31, 1, y, 31 )
        if ( btest( y, 0 ) ) then
          self%mt(kk) = int(ieor( ieor( self%mt(kk+m), ishft( y, -1 ) ), matrix_a ))
        else
          self%mt(kk) = int(ieor( ieor( self%mt(kk+m), ishft( y, -1 ) ), matrix_b ))
        end if
      end do
      do kk = n-m+1, n-1, 1
        y = ibits( self%mt(kk+1), 0, 31 )
        call mvbits( self%mt(kk), 31, 1, y, 31 )
        if ( btest( y, 0 ) ) then
          self%mt(kk) = int(ieor( ieor( self%mt(kk+m-n), ishft( y, -1 ) ), matrix_a ))
        else
          self%mt(kk) = int(ieor( ieor( self%mt(kk+m-n), ishft( y, -1 ) ), matrix_b ))
        end if
      end do
      y = ibits( self%mt(1), 0, 31 )
      call mvbits( self%mt(n), 31, 1, y, 31 )
      if ( btest( y, 0 ) ) then
        self%mt(kk) = int(ieor( ieor( self%mt(m), ishft( y, -1 ) ), matrix_a ))
      else
        self%mt(kk) = int(ieor( ieor( self%mt(m), ishft( y, -1 ) ), matrix_b ))
      end if
      self%mti = 1_wi
    end if
    y = self%mt(self%mti)
    self%mti = self%mti + 1_wi
    ! Tempering
    y = ieor( y, ishft( y, -11) )
    y = int(ieor( y, iand( ishft( y, 7 ), temper_a ) ))
    y = int(ieor( y, iand( ishft( y, 15 ), temper_b ) ))
    y = ieor( y, ishft( y, -18 ) )

  end function random_int32

  ! generates a random number on [0,0x7fffffff]-interval
  function random_int31(self) result( i )

    implicit none

    type(mt19937_t), intent(inout) :: self

    intrinsic :: ishft

    integer( kind = wi )  :: i

    i = ishft(random_int32(self), -1)

  end function random_int31

  ! generates a random number on [0,1]-real-interval
  function random_real1(self) result( r )

    implicit none

    type(mt19937_t), intent(inout) :: self

    real( kind = wr )  :: r

    integer( kind = wi )  :: a, a1, a0

    a = random_int32(self)
    a0 = ibits( a, 0, hbs )
    a1 = ibits( a, hbs, hbs )
    r = real( a0, kind = wr ) / 4294967295.0_wr
    r = real( a1, kind = wr ) * ( 65536.0_wr / 4294967295.0_wr ) + r
    ! divided by 2^32-1

  end function random_real1

  ! generates a random number on [0,1)-real-interval
  function random_real2(self) result( r )

    implicit none

    intrinsic :: ibits

    type(mt19937_t), intent(inout) :: self

    real( kind = wr )  :: r

    integer( kind = wi )  :: a, a1, a0

    a = random_int32(self)
    a0 = ibits( a, 0, hbs )
    a1 = ibits( a, hbs, hbs )
    r = real( a0, kind = wr ) / 4294967296.0_wr
    r = real( a1, kind = wr ) / 65536.0_wr + r
    ! divided by 2^32

  end function random_real2

  ! generates a random number on (0,1)-real-interval
  function random_real3(self) result( r )

    implicit none

    type(mt19937_t), intent(inout) :: self

    real( kind = wr )  :: r

    integer( kind = wi )  :: a, a1, a0

    a = random_int32(self)
    a0 = ibits( a, 0, hbs )
    a1 = ibits( a, hbs, hbs )
    r = ( real( a0, kind = wr ) + 0.5_wr ) / 4294967296.0_wr
    r = real( a1, kind = wr ) / 65536.0_wr + r
    ! divided by 2^32

  end function random_real3

  ! generates a random number on [0,1) with 53-bit resolution
  function random_res53(self)  result( r )

    implicit none

    intrinsic :: ishft

    type(mt19937_t), intent(inout) :: self

    real( kind = wr )  :: r

    integer( kind = wi )  :: a, a0, a1
    integer( kind = wi )  :: b, b0, b1

    a = ishft(random_int32(self), -5 )
    a0 = ibits( a, 0, hbs )
    a1 = ibits( a, hbs, hbs )
    b = ishft(random_int32(self), -6 )
    b0 = ibits( b, 0, hbs )
    b1 = ibits( b, hbs, hbs )
    r = real( a1, kind = wr ) / 2048.0_wr
    r = real( a0, kind = wr ) / 134217728.0_wr + r
    r = real( b1, kind = wr ) / 137438953472.0_wr + r
    r = real( b0, kind = wr ) / 9007199254740992.0_wr + r

  end function random_res53
  ! These real versions are due to Isaku Wada, 2002/01/09 added

#ifndef NCSU_NO_NETCDF
subroutine check_ncrc(rc, sbrtn, filename)

   use netcdf
   use ncsu_constants, only : ERR_UNIT
   use ncsu_sander_proxy, only : terminate

   implicit none

   integer, intent(in) :: rc
   character(*), intent(in) :: sbrtn ! this is FORTRAN, after all :-)
   character(*), intent(in) :: filename

   character(len = 80) :: errmsg

   if (rc.ne.nf90_noerr) then
      errmsg = nf90_strerror(rc)
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a,a,a/)') NCSU_ERROR, &
         trim(sbrtn), '(filename=''', trim(filename), ''') : ', trim(errmsg)
      call terminate()
   end if

end subroutine check_ncrc

subroutine mt19937_save(self, filename)

   use netcdf

   implicit none

   type(mt19937_t), intent(in) :: self
   character(*), intent(in) :: filename

   character(*), parameter :: sbrtn = 'mt19937%save'

   integer :: rc, setid, imtinit

   rc = nf90_create(filename, cmode = nf90_clobber, ncid = setid)
   call check_ncrc(rc, sbrtn, filename)

   rc = nf90_put_att(setid, nf90_global, 'mt', self%mt)
   call check_ncrc(rc, sbrtn, filename)

   if (self%mtinit) then
      imtinit = 1
   else
      imtinit = 0
   end if

   rc = nf90_put_att(setid, nf90_global, 'mtinit', imtinit)
   call check_ncrc(rc, sbrtn, filename)

   rc = nf90_put_att(setid, nf90_global, 'mti', self%mti)
   call check_ncrc(rc, sbrtn, filename)

   rc = nf90_close(setid)
   call check_ncrc(rc, sbrtn, filename)

end subroutine mt19937_save

subroutine mt19937_load(self, filename)

   NCSU_USE_AFAILED

   use netcdf

   implicit none

   type(mt19937_t), intent(inout) :: self
   character(*), intent(in) :: filename

   character(*), parameter :: sbrtn = 'mt19937%load'

   integer :: rc, setid, imtinit

   rc = nf90_open(filename, nf90_nowrite, setid)
   call check_ncrc(rc, sbrtn, filename)

   rc = nf90_get_att(setid, nf90_global, 'mt', self%mt)
   call check_ncrc(rc, sbrtn, filename)

   rc = nf90_get_att(setid, nf90_global, 'mtinit', imtinit)
   call check_ncrc(rc, sbrtn, filename)

   select case(imtinit)
      case (0)
         self%mtinit = .false.
      case (1)
         self%mtinit = .true.
      case default
         ncsu_assert_not_reached()
         continue
   end select

   rc = nf90_get_att(setid, nf90_global, 'mti', self%mti)
   call check_ncrc(rc, sbrtn, filename)

   rc = nf90_close(setid)
   call check_ncrc(rc, sbrtn, filename)

end subroutine mt19937_load
#endif /* NCSU_NO_NETCDF */

#ifdef MPI
subroutine mt19937_bcast(self, comm, root)

   use ncsu_utils

   implicit none

   type(mt19937_t), intent(inout) :: self
   integer, intent(in) :: comm, root

#  include "ncsu-mpi.h"

   integer :: error

   ncsu_assert(comm.ne.MPI_COMM_NULL)

   call mpi_bcast(self%mt, size(self%mt), MPI_INTEGER, root, comm, error)
   ncsu_assert(error.eq.0)

   call mpi_bcast(self%mtinit, 1, MPI_INTEGER, root, comm, error)
   ncsu_assert(error.eq.0)

   call mpi_bcast(self%mti, 1, MPI_INTEGER, root, comm, error)
   ncsu_assert(error.eq.0)

end subroutine mt19937_bcast
#endif /* MPI */

#endif /* defined(MPI) && defined(BINTRAJ) */
end module mt19937

