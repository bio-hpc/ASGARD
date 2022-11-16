! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Main interface to EVB                                                  |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

#ifdef LES
   subroutine evb_ntrfc ( x, f, ener, ix, ipairs, vel0_nrg_sum )
#else
   subroutine evb_ntrfc ( x, f, ener, mass, ix, ipairs)
#endif

   use evb_amber, only: morsify_initialized, modvdw_initialized 
   use evb_parm,  only: xch_type, nmorse, nmodvdw 
   use evb_data,  only: evb_Hmat, evb_frc
   use constants, only: A_TO_BOHRS, AU_TO_KCAL
   use les_data
   use schlegel,  only: ncoord, vmm, dvmm, ddvmm
   use state

#if defined(LES)
   use pimd_vars, only: natomCL, nbead, nbead_inv, nrg_all
   use evb_pimd,  only: evb_begin, evb_end, lpimd_dcrypt, natomPCL, natomPQM, bead_dcrypt &
                      , atomCL_dcrypt, atomQM_dcrypt &
                      , pie, pif, nrg_bead, fpimd &
                      , evb_mat_bead, evb_vec0_bead, vel0_bead
   use evb_data, only : evb_vel0 
   use evb_parm, only : nevb
#else
   use evb_amber, only : xq, xnrg, xf
#endif /* LES */

   implicit none

#  include "parallel.h"
   include 'mpif.h'
#  include "extra.h"
#  include "../include/memory.h"

   integer, intent(in   ) :: ix(*), ipairs(*)
   _REAL_ , intent(in   ) :: x(3,natom)
#ifdef LES
   _REAL_ , intent(inout) :: vel0_nrg_sum
#else
   _REAL_ , intent(in   ) :: mass(natom)
#endif
   _REAL_ , intent(inout) :: f(3,natom)
   type(state_rec)        :: ener
  
   !  ..........................................................................

   _REAL_  :: v, dv(ncoord), ddv(ncoord,ncoord), qint(ncoord)
 
#ifdef LES
   integer :: mm, nslice, npimd_set 

   _REAL_  :: f0_bead(3,natomCL,nbead)
   _REAL_  :: evb_vec0_lpimd(nevb,lpimd_size), evb_vec0_tmp(nevb,nbead) &
            , evb_mat_lpimd(nevb,nevb,lpimd_size), evb_mat_tmp(nevb,nevb,nbead) &
            , vel0_lpimd(lpimd_size), vel0_tmp(nbead) &
            , f0_lpimd(3,natomCL,lpimd_size)

   _REAL_  :: nrg_diabatic(nevb) 
   _REAL_  :: nrg_morsify(nbead), nrg_vdw(nbead)
   _REAL_  :: piff(3,natomCL,nevb), piqq(3,natomCL)
   _REAL_  :: pimd_vec0sq(nevb)

   _REAL_  :: vmm_bead(nbead), dvmm_bead(ncoord,nbead) &
            , ddvmm_bead(ncoord,ncoord,nbead)

   _REAL_  :: vmm_evb(nbead,nevb), dvmm_evb(ncoord,nbead,nevb) &
            , ddvmm_evb(ncoord,ncoord,nbead,nevb)
   integer :: m, n, nn, nnn


#else
   integer :: nbead
   _REAL_  :: nrg_morsify(1), nrg_vdw(1)
#endif /* LES */

   integer :: ierr

#ifndef LES
   nbead = 1
#endif

!  +---------------------------------------------------------------------------+
!  |  (0) Ultimately, this should be done at the parmtop level.                |
!  |  (1) Change traditional harmonic bond interaction to Morse.               |
!  |  (2) Modify certain VDW interactions involved in the topology change;     |
!  |      the simpliest case is to exclude them.                               |
!  |  ener(2) :: van der Waals                                                 |
!  |  ener(5) :: bond energy                                                   |
!  +:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::+
!  |  Note: morsify & mod_vdw operates on all the coordinates (even for LES)   |
!  |  & only the masters perform this force field change.  For EVB/(LES)-PIMD, |
!  |  the morsify and mod_vdw parameters are replicated for the other beads    |
!  |  based on the classical EVB input specifications in evb_pimd_init.        |
!  +---------------------------------------------------------------------------+

   if( master ) then

      nrg_morsify(:) = 0.0d0
      if( nmorse > 0 ) then
         if( .not. morsify_initialized ) call morsify_init ( ix )
         call morsify ( x, nrg_morsify, f, natom*3, nbead )
         ener%pot%bond = ener%pot%bond + sum( nrg_morsify(:) )
         ener%pot%tot  = ener%pot%tot  + sum( nrg_morsify(:) )
      endif

      nrg_vdw(:) = 0.0d0
      if( nmodvdw > 0 ) then
         if( .not. modvdw_initialized ) call modvdw_init ( ix, ipairs )
         call mod_vdw ( x, nrg_vdw, f, natom*3, nbead )
         ener%pot%vdw = ener%pot%vdw + sum( nrg_vdw(:) )
         ener%pot%tot = ener%pot%tot + sum( nrg_vdw(:) )
      endif

   endif 

#if defined(LES)

!  +---------------------------------------------------------------------------+
!  |  Each master has a diabatic state energy and corresponding forces.        |
!  |  MPI_ALLGATHER ( commmaster ) these into arrays of size (:,nevb)          |
!  +---------------------------------------------------------------------------+

   if( master ) then
      nrg_bead(:) = nrg_all(:) + nrg_morsify(:) + nrg_vdw(:)
      if( trim( adjustl( xch_type ) ) /= "dist_gauss" ) then
         do m = 1, natomPCL
            mm = atomCL_dcrypt(m)
            f(:,mm) = f(:,mm) * nbead_inv
         enddo
         call mpi_allgather ( nrg_bead, nbead, MPI_DOUBLE_PRECISION, pie, nbead &
                            , MPI_DOUBLE_PRECISION, commmaster, ierr )
         call mpi_allgather ( f, natom*3, MPI_DOUBLE_PRECISION, pif &
                            , natom*3, MPI_DOUBLE_PRECISION, commmaster, ierr )
      endif
!  +---------------------------------------------------------------------------+
!  |  Compute diabatic states for each bead                                    |
!  +---------------------------------------------------------------------------+

      do n = 1, nbead

         if( trim( adjustl( xch_type ) ) == "dist_gauss" ) then
            do m = 1, natomCL
               mm = bead_dcrypt(m,n)
               piqq(:,m  ) = x(:,mm  )
            enddo
         endif

  if( trim( adjustl( xch_type ) ) == "dist_gauss" ) then

         call cart2internal ( piqq * A_TO_BOHRS, qint )
         call schlegel_vmm ( qint, v, dv, ddv )

           vmm_bead(    n) =   v
          dvmm_bead(  :,n) =  dv(  :)
         ddvmm_bead(:,:,n) = ddv(:,:)

  endif

      enddo

  if( trim( adjustl( xch_type ) ) == "dist_gauss" ) then

      call mpi_allgather ( vmm_bead, nbead, MPI_DOUBLE_PRECISION, vmm_evb &
                         , nbead, MPI_DOUBLE_PRECISION, commmaster, ierr )
      call mpi_allgather ( dvmm_bead, ncoord*nbead, MPI_DOUBLE_PRECISION, dvmm_evb &
                         , ncoord*nbead, MPI_DOUBLE_PRECISION, commmaster, ierr )
      call mpi_allgather ( ddvmm_bead, ncoord*ncoord*nbead, MPI_DOUBLE_PRECISION, ddvmm_evb &
                         , ncoord*ncoord*nbead, MPI_DOUBLE_PRECISION, commmaster, ierr )

  endif

   endif
!  +---------------------------------------------------------------------------+
!  |  EVB/PIMD requires all PEs to have the diabatic state energies & forces   |
!  |  ... so MPI_BCAST.  Note the high communication cost of pif ... instead   |
!  |  we keep the diabatic state forces on the masters and formulate each      |
!  |  diabatic state contribution to the Hellman-Feynman forces and the reduce |
!  +---------------------------------------------------------------------------+

   call mpi_bcast ( pie, nbead*nevb, MPI_DOUBLE_PRECISION, 0, commsander, ierr )

  if( trim( adjustl( xch_type ) ) == "dist_gauss" ) then

   call mpi_bcast ( vmm_evb, nbead*nevb, MPI_DOUBLE_PRECISION, 0, commsander, ierr )
   call mpi_bcast ( dvmm_evb, ncoord*nbead*nevb, MPI_DOUBLE_PRECISION, 0, commsander, ierr )
   call mpi_bcast ( ddvmm_evb, ncoord*ncoord*nbead*nevb, MPI_DOUBLE_PRECISION, 0, commsander, ierr )

  endif
  
!  +---------------------------------------------------------------------------+
!  |  Scatter the EVB calculation for each bead onto PEs defined based on      |
!  |  arrays evb_begin & evb_end                                               |
!  +---------------------------------------------------------------------------+

   if( comm_lpimd /= MPI_COMM_NULL ) then

      vel0_tmp(:) = 0.0d0
      evb_vec0_tmp(:,:) = 0.0d0
      evb_mat_tmp(:,:,:) = 0.0d0

      npimd_set = 0

   fpimd (:,:) = 0.0d0
   pimd_vec0sq(:) = 0.0d0

      do n = evb_begin(lpimd_rank+1), evb_end(lpimd_rank+1)

  if( trim( adjustl( xch_type ) ) == "dist_gauss" ) then

           vmm(    :) =   vmm_evb(    n,:)
          dvmm(  :,:) =  dvmm_evb(  :,n,:)
         ddvmm(:,:,:) = ddvmm_evb(:,:,n,:)

  endif
         npimd_set = npimd_set + 1

         if( trim( adjustl( xch_type ) ) /= "constant" ) then
            do m = 1, natomCL
               mm = bead_dcrypt(m,n)
               piqq(:,m  ) = x(:,mm  )

               if( cnum(mm) == 0 ) then
                  piff(:,m,:) = pif(:,mm,:) * nbead_inv
               else
                  piff(:,m,:) = pif(:,mm,:)
               endif

            enddo

         endif

         call evb_matrix ( pie(n,:), piqq, natomCL*3 )

!        if( trim( adjustl( xch_type ) ) /= "constant" ) &
#ifdef LES
            call evb_force ( piff, piqq, natomCL*3 )
#else
            call evb_force ( piff, piqq, mass, natomCL*3 )
#endif

!  +---------------------------------------------------------------------------+
!  |  MPI_GATHER onto root of comm_lpimd                                       |
!  +---------------------------------------------------------------------------+

         call mpi_gather ( evb_Hmat%evb_vec0, nevb, MPI_DOUBLE_PRECISION &
                         , evb_vec0_lpimd, nevb, MPI_DOUBLE_PRECISION, 0 &
                         , comm_lpimd, ierr )
         call mpi_gather ( evb_Hmat%evb_mat, nevb*nevb, MPI_DOUBLE_PRECISION &
                         , evb_mat_lpimd, nevb*nevb, MPI_DOUBLE_PRECISION, 0 &
                         , comm_lpimd, ierr )
         call mpi_gather ( evb_vel0%evb_nrg, 1, MPI_DOUBLE_PRECISION &
                         , vel0_lpimd, 1, MPI_DOUBLE_PRECISION, 0 &
                         , comm_lpimd, ierr )

         if( trim( adjustl( xch_type ) ) /= "constant" ) &
            call mpi_gather ( evb_vel0%evb_f, natomCL*3, MPI_DOUBLE_PRECISION &
                            , f0_lpimd, natomCL*3, MPI_DOUBLE_PRECISION, 0 &
                            , comm_lpimd, ierr )

!  +---------------------------------------------------------------------------+
!  |  Decrypt bead contributions on root(comm_lpimd)                           |
!  +---------------------------------------------------------------------------+

         do nn = 1, lpimd_size
            nslice = lpimd_dcrypt(nn,npimd_set)

                vel0_tmp(    nslice) = vel0_lpimd(nn)
            evb_vec0_tmp(:  ,nslice) = evb_vec0_lpimd(:,nn)
             evb_mat_tmp(:,:,nslice) = evb_mat_lpimd(:,:,nn)
            if( trim( adjustl( xch_type ) ) /= "constant" ) &
               f0_bead(:,:,nslice) = f0_lpimd(:,:,nn)

         enddo
      enddo
   endif

!  +---------------------------------------------------------------------------+
!  |  Since root(comm_lpimd) corresponds to root(commmaster), MPI_ALLREDUCE    |
!  |  will synchronize all PIMD bead data on the masters                       |
!  +---------------------------------------------------------------------------+

   if( master ) then

      call mpi_allreduce ( vel0_tmp, vel0_bead, nbead, MPI_DOUBLE_PRECISION &
                         , MPI_SUM, commmaster, ierr )
      call mpi_allreduce ( evb_vec0_tmp, evb_vec0_bead, nevb*nbead &
                         , MPI_DOUBLE_PRECISION, MPI_SUM, commmaster, ierr )
      call mpi_allreduce ( evb_mat_tmp, evb_mat_bead, nevb*nevb*nbead &
                         , MPI_DOUBLE_PRECISION, MPI_SUM, commmaster, ierr )

   endif

!  +---------------------------------------------------------------------------+
!  |  Now do the LES-PIMD forces                                               |
!  +---------------------------------------------------------------------------+
!  .............................................................................
!  :  For coordinate-dependent EVB exchange type                               :
!  `````````````````````````````````````````````````````````````````````````````

   if( master ) then

      fpimd(:,:) = 0.0d0
      pimd_vec0sq(:) = 0.0d0

      if( trim( adjustl( xch_type ) ) /= "constant" ) then

!        do n = 1, nbead
         do n = evb_begin(lpimd_rank+1), evb_end(lpimd_rank+1)
            do m = 1, natomCL
               mm = bead_dcrypt(m,n)
               fpimd(:,mm) = fpimd(:,mm) + f0_bead(:,m,n)
            enddo
         enddo

!  .............................................................................
!  :  for "constant" EVB exchange type                                         :
!  `````````````````````````````````````````````````````````````````````````````

      else

         do n = 1, nbead
            pimd_vec0sq(:) = pimd_vec0sq(:) + evb_vec0_bead(:,n)**2
         enddo

!  +---------------------------------------------------------------------------+
!  |  Form Hellmann-Feynman forces                                             |
!  |                                                                           |
!  |  F(:) = - sum(k,l) C_k1 * C_l1 * grad( H_kl )                             |
!  +---------------------------------------------------------------------------+

!  +,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,+
!  |  Loop over the classical atoms                                            |
!  +'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''+

         nn = masterrank + 1
         do m = 1, natomPCL
            mm = atomCL_dcrypt(m)
            fpimd(:,mm) = fpimd(:,mm) + pimd_vec0sq(nn) * pif(:,mm,nn)
         enddo

!  +,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,+
!  |  Loop over the quantized atoms                                            |
!  +'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''+

         nn = masterrank + 1
         do m = 1, natomPQM
            mm  = atomQM_dcrypt(m,1)
            nnn = atomQM_dcrypt(m,2)
            fpimd(:,mm) = fpimd(:,mm) + evb_vec0_bead(nn,nnn)**2 * pif(:,mm,nn)
         enddo

      endif

      call mpi_allreduce ( fpimd, f, natom*3, MPI_DOUBLE_PRECISION, MPI_SUM &
                         , commmaster, ierr )

   endif

!  +---------------------------------------------------------------------------+
!  |  Broadcast PIMD forces. Since the forces on the masters are synchronized, |
!  |  we broadcast from commsander root (master) to the respective PEs         |
!  +---------------------------------------------------------------------------+

   if( master ) then
      if( trim( adjustl( xch_type ) ) == "dist_gauss" ) then
         do nn = 1, nevb
            nrg_diabatic(nn) = sum( evb_mat_bead(nn,nn,:) )
         enddo
         ener%pot     = null_potential_energy_rec
         ener%pot%tot = nrg_diabatic(masterrank+1)
      endif
      vel0_nrg_sum = sum( vel0_bead(:) )
      evb_vel0%evb_nrg = vel0_nrg_sum
      evb_frc%evb_nrg = vel0_nrg_sum

   endif

!  call mpi_bcast ( f, natom*3, MPI_DOUBLE_PRECISION, 0, commsander, ierr )
!  call mpi_bcast ( x, natom*3, MPI_DOUBLE_PRECISION, 0, commsander, ierr )

   call mpi_bcast ( f, natom*3, MPI_DOUBLE_PRECISION, 0, commworld, ierr )
   call mpi_bcast ( x, natom*3, MPI_DOUBLE_PRECISION, 0, commworld, ierr )

#else /* LES */

   if( master ) then 

      xq(:) = reshape( x(:,:), (/  3 * natom /) )

      if( trim( adjustl( xch_type ) ) /= "dist_gauss" ) then
         call mpi_allgather ( ener%pot%tot, 1, MPI_DOUBLE_PRECISION, xnrg, 1 &
                            , MPI_DOUBLE_PRECISION, commmaster, ierr )
         call mpi_allgather ( f, natom*3, MPI_DOUBLE_PRECISION, xf, natom*3 &
                            , MPI_DOUBLE_PRECISION, commmaster, ierr )
         call mpi_barrier ( commmaster, ierr )
         call evb_matrix ( xnrg, xq, natom*3 )
      else
!  +---------------------------------------------------------------------------+
!  |  Compute diabatic states                                                  |
!  +---------------------------------------------------------------------------+

         call cart2internal ( xq * A_TO_BOHRS, qint )
         call schlegel_vmm ( qint, v, dv, ddv )

         call mpi_allgather ( v, 1, MPI_DOUBLE_PRECISION, vmm, 1 &
                            , MPI_DOUBLE_PRECISION, commmaster, ierr )
         call mpi_allgather ( dv, ncoord, MPI_DOUBLE_PRECISION, dvmm, ncoord &
                            , MPI_DOUBLE_PRECISION, commmaster, ierr )
         call mpi_allgather ( ddv, ncoord*ncoord, MPI_DOUBLE_PRECISION, ddvmm &
                            , ncoord*ncoord, MPI_DOUBLE_PRECISION, commmaster, ierr )

         call mpi_bcast ( v, 1, MPI_DOUBLE_PRECISION, 0, commsander, ierr )
         call mpi_bcast ( dv, ncoord, MPI_DOUBLE_PRECISION, 0, commsander, ierr )
         call mpi_bcast ( ddv, ncoord*ncoord, MPI_DOUBLE_PRECISION, 0, commsander, ierr )

      endif

#ifdef LES
      call evb_force ( xf, xq, natom*3 )
#else
      call evb_force ( xf, xq, mass, natom*3 )
#endif
      f(:,:) = reshape( evb_frc%evb_f(:), (/ 3, natom /) )

      if( trim( adjustl( xch_type ) ) == "dist_gauss" ) then
         ener%pot     = null_potential_energy_rec
         ener%pot%tot = evb_Hmat%evb_mat(masterrank+1,masterrank+1)
      endif
   endif 

   call mpi_bcast( f, natom*3, MPI_DOUBLE_PRECISION, 0, commworld, ierr )
   call mpi_bcast( x, natom*3, MPI_DOUBLE_PRECISION, 0, commworld, ierr )

#endif /* LES */

end subroutine evb_ntrfc


!  +---------------------------------------------------------------------------+
!  |  Allocate storage space for xq, xf, xnrg for EVB                          |
!  +---------------------------------------------------------------------------+

   subroutine evb_amber_alloc

   use evb_parm,  only: nevb 
   use evb_amber, only: xnrg, xq, xf 
#ifdef LES
   use pimd_vars, only: natomCL
#endif

   implicit none

#  include "../include/memory.h"

   integer :: alloc_error

   allocate( xnrg(nevb), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
#ifdef LES
   allocate( xq(natomCL*3), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( xf(natomCL*3,nevb), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
#else
   allocate( xq(natom*3), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
   allocate( xf(natom*3,nevb), stat = alloc_error )
   REQUIRE( alloc_error == 0 )
#endif /*LES*/

   end subroutine evb_amber_alloc

!  +---------------------------------------------------------------------------+
!  |  Deallocate storage space for xnrg, xq, xf for EVB                        |
!  +---------------------------------------------------------------------------+

   subroutine evb_amber_dealloc

   use evb_amber, only: xnrg, xq, xf 

   implicit none

   integer :: dealloc_error

   deallocate( xnrg, xq, xf, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )

   end subroutine evb_amber_dealloc
