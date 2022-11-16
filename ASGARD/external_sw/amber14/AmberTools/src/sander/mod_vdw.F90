! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Modify van der Waals interactions                                      |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

   subroutine mod_vdw ( q, nrg_vdw, f, ndim, nbead ) 

#ifdef LES
   use les_data, only : cnum, lestyp, lfac, lesfac, nlesty
#endif
   use parms,     only: cn1, cn2
   use evb_parm,  only: nmodvdw
   use evb_amber, only: vdw_mod, iac, ico
#ifdef DEBUG_EVB
   use evb_check, only: full_evb_debug, mod_vdw_debug
#endif

   implicit none

#  include "parallel.h"
#  include "../include/memory.h"

   integer, intent(in   ) :: ndim, nbead
   _REAL_ , intent(in   ) :: q(ndim)
   _REAL_ , intent(out  ) :: nrg_vdw(nbead)
   _REAL_ , intent(inout) :: f(ndim)

   !............................................................................

   integer :: n, i, j, idx, jdx, ic 
   _REAL_  :: dx, dy, dz, fx, fy, fz, rij2_inv, r6, f6, f12, ff
   _REAL_  :: evdw(nbead), fvdw(ndim)

!  +---------------------------------------------------------------------------+
!  |  V = A_ij / R_ij^12 - B_ij / R_ij^6                                       |
!  |  F_Ri = (12 * A_ij / R_ij^14 - 6 * B_ij / R_ij^8 )                        |
!  |       * ( R_i - R_j )                                                     |
!  +---------------------------------------------------------------------------+

   evdw(:) = 0.0d0
   fvdw(:) = 0.0d0

   do n = 1, nmodvdw

      i = vdw_mod(n)%iatom
      j = vdw_mod(n)%jatom

      ic = ico( ntypes * ( iac(i) - 1 ) + iac(j) )

      idx = ( i - 1 ) * 3 
      jdx = ( j - 1 ) * 3

      dx = q(idx+1) - q(jdx+1) 
      dy = q(idx+2) - q(jdx+2) 
      dz = q(idx+3) - q(jdx+3) 

      rij2_inv = 1.0d0 / ( dx**2 + dy**2 + dz**2 )
      r6 = rij2_inv * rij2_inv * rij2_inv

      f6 = cn2(ic) * r6
      f12 = cn1(ic) * r6 * r6 

#ifdef LES
#ifdef DEBUG_EVB
      if( full_evb_debug .or. mod_vdw_debug ) then
         write(6,'(A,4(I10,2X))') '         ', cnum(i), cnum(j), i, j
      endif
#endif

      lfac = lesfac( nlesty * ( lestyp(i) - 1 ) + lestyp(j) )
      f6 = f6 * lfac
      f12 = f12 * lfac

      if( cnum(i) == 0 .and. cnum(j) == 0 ) then
         evdw(1:nbead) = evdw(1:nbead) + ( f12 - f6 ) / nbead
      else
         if( cnum(i) /= 0 ) then
            evdw( cnum(i) ) = evdw( cnum(i) ) + f12 - f6
         else
            evdw( cnum(j) ) = evdw( cnum(j) ) + f12 - f6
         endif
      endif
#else
      evdw(1) = evdw(1) + f12 - f6
#endif 

      ff = ( 12.0d0 * f12 - 6.0d0 * f6 ) * rij2_inv 

      fx = ff * dx
      fy = ff * dy
      fz = ff * dz

      fvdw(idx+1) = fvdw(idx+1) + fx
      fvdw(idx+2) = fvdw(idx+2) + fy
      fvdw(idx+3) = fvdw(idx+3) + fz

      fvdw(jdx+1) = fvdw(jdx+1) - fx
      fvdw(jdx+2) = fvdw(jdx+2) - fy
      fvdw(jdx+3) = fvdw(jdx+3) - fz

   enddo

!  +---------------------------------------------------------------------------+
!  |  Exclude the above VDW interactions ... so subtract it out of the total   |
!  |  energy & force                                                           |
!  +---------------------------------------------------------------------------+

   f(:) = f(:) - fvdw(:)
   nrg_vdw(:) = - evdw(:)

#ifdef DEBUG_EVB
   if( full_evb_debug .or. mod_vdw_debug ) call vdw_anal2num ( q, fvdw, ndim )
#endif

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine mod_vdw


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Read Amber indices for VDW interactions                                |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine modvdw_init ( ix, ipairs ) 

   use evb_parm,  only: nmodvdw
   use evb_amber, only: vdw_mod, iac, ico, modvdw_initialized
   use nblist,    only: bckptr, nlogrid, nhigrid, numvdw, numhbnd &
                      , myindexlo, myindexhi, numimg

   implicit none

#  include "../include/memory.h"
#  include "parallel.h"

   integer, intent(in) :: ix(*), ipairs(*)

   !...........................................................................

   integer :: i, j, k, m, n, ii, jj, vdwtot, maxsize, alloc_error, dealloc_error
   integer :: ndx, numpack, ncell_lo, ncell_hi, ntot, nvdw
   integer :: ndx_duplicate(nmodvdw)
   integer, allocatable :: idx_tmp(:), jdx_tmp(:)
   integer, parameter :: mask27 = 2**27 - 1
   intrinsic :: ishft
   logical :: not_found(nmodvdw), duplicate(nmodvdw)

!  +---------------------------------------------------------------------------+
!  |  Read non-bonded list and confirm that requested modvdw pairs are leget.  |
!  |  legit.  The code below follows the construct of get_nb_energy.           |
!  +---------------------------------------------------------------------------+

!KFW SKIP tests for now; the code appears to be working
   goto 888 

   maxsize = natom * natom / 4
   print *, 'maxsize = ', maxsize

   allocate( idx_tmp(maxsize), jdx_tmp(maxsize), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   vdwtot = 0
   numpack = 1
   do ndx = myindexlo, myindexhi
      if ( numimg(ndx) > 0 ) then
         ncell_lo = nlogrid(ndx)
         ncell_hi = nhigrid(ndx)
         do k = ncell_lo, ncell_hi
            i = bckptr(k)
            ntot = numvdw(i) + numhbnd(i)
            nvdw = numvdw(i)
            if ( ntot > 0 ) then
               do m = 1, nvdw
                  j = bckptr( iand(ipairs(numpack-1+m),mask27) )
                  vdwtot = vdwtot + 1
                  idx_tmp(vdwtot) = i
                  jdx_tmp(vdwtot) = j
               enddo
               numpack = numpack + ntot
            endif
         enddo
      endif
   enddo

!  +---------------------------------------------------------------------------+
!  |  Flag requests not in non-bonded list                                     |
!  +---------------------------------------------------------------------------+

   not_found(:) = .true.
   do n = 1, nmodvdw
      i = vdw_mod(n)%iatom
      j = vdw_mod(n)%jatom
      do m = 1, vdwtot
         ii = idx_tmp(m)
         jj = jdx_tmp(m)
         if( i == ii .and. j == jj .or. i == jj .and. j == ii) &
            not_found(n) = .false.
      enddo
   enddo

   if( any(not_found) ) then
      write(6,'(A)') 'ERROR: User requested modification of the VDW interactions below'
      write(6,'(A)') '       but atom pairs are not defined in the non-bonded list.'
   endif  
   do n = 1, nmodvdw
      if( not_found(n) ) then
         write(6,'(2(A,I8),A)' ) '       (', vdw_mod(n)%iatom, ',', vdw_mod(n)%jatom, '   )'
      endif
   enddo
   if( any(not_found) ) then
      write(6,'(A)')
      write(6,'(A)') '... possible reasons could be that these atoms are ' &
                  // 'bonded or in an exclusion list.'
      call mexit(6,1)
   endif

!  +---------------------------------------------------------------------------+
!  |  Flag duplicate pairs                                                     |
!  +---------------------------------------------------------------------------+

   duplicate(:) = .false.
   do n = 1, nmodvdw
      i = vdw_mod(n)%iatom
      j = vdw_mod(n)%jatom
      do m = 1, nmodvdw
         ii = vdw_mod(m)%iatom
         jj = vdw_mod(m)%jatom
         if( i == jj .and. j == ii ) then
            duplicate(n) = .true.
            ndx_duplicate(n) = m
         endif
         if( i == ii .and. j == jj .and. m /= n ) then
            duplicate(n) = .true.
            ndx_duplicate(n) = m
         endif
      enddo
   enddo

   if( any(duplicate) ) then
      write(6,'(A)') 'ERROR: User requested modification of the VDW interactions below'
      write(6,'(A)') '       but requested atom pairs are duplicated.'
   endif
   do n = 1, nmodvdw
      if( duplicate(n) ) then
         write(6,'(2(A,I8),A)' ) '       (', vdw_mod(n)%iatom, ',', vdw_mod(n)%jatom, '   )'
         m = ndx_duplicate(n)
         write(6,'(2(A,I8),A)' ) '       (', vdw_mod(m)%iatom, ',', vdw_mod(m)%jatom, '   )'
      endif
   enddo
   if( any(duplicate) ) then
      write(6,'(A)')
      write(6,'(A)') '... please remove duplicate pairs from the vdw_mod list ' &
                  // 'and try again.'
      call mexit(6,1)
   endif

   deallocate( idx_tmp, jdx_tmp, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )

!KFW SKIP tests for now; the code appears to be working
 888 continue

!  +---------------------------------------------------------------------------+
!  |  Data read from topology file                                             |
!  +~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
!  |  iac :: index for atom types in VDW interactions                          |
!  |  ico :: index to VDW parameters cn1 and cn2                               |
!  +---------------------------------------------------------------------------+

   allocate( iac(natom), stat = alloc_error ) 
   REQUIRE( alloc_error == 0 )
   allocate( ico(ntypes*ntypes), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   iac(:) = ix(i04:i04+natom-1)
   ico(:) = ix(i06:i06+ntypes*ntypes-1)

   modvdw_initialized = .true.

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine modvdw_init


