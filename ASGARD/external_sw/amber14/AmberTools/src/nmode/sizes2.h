
!  --- sizes for normal mode analysis program:

!    MAXATOM is maximum number of atoms
!    MAXINT  is maximum number of internal coordinates on which
!            projections are to be made
!    MAXVEC  is the maximum number of normal mode eigenvectors that
!            will be used
!    MEMDRV  is a general workspace; program will complain if this
!            is too small and tell you how large to make it.

parameter (memdrv=10000000)
parameter (maxatom=8000)
parameter (maxint=500)
parameter (maxvec=500)

!     ----- SET THE LIMITS OF SOME ARRAY BOUNDS -----

parameter (maxdih = 5500)
parameter (maxdia = 3500)
parameter (maxinb = 9000)
parameter (maxbon = 2600)
parameter (maxbnh = 2000)
parameter (maxang = 3000)
parameter (maxanh = 3000)
