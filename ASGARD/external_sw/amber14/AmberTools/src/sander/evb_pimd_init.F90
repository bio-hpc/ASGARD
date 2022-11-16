! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Initialize EVB-PIMD variables                                          |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

#include "../include/assert.fh"
#include "../include/dprec.fh"

   module evb_pimd

   use les_data
   use pimd_vars, only: ipimd
   use evb_parm,  only: nevb, nbias, evb_dyn, dbonds_RC, bond_RC
   use miller,    only: ndiv, div_ndx

   implicit none

#  include "../include/memory.h"
   save 

   integer :: natomCL, natomPCL, natomPQM, nbead
   integer, allocatable :: bead_dcrypt(:,:) &
                         , atomCL_dcrypt(:), atomQM_dcrypt(:,:)
   integer, allocatable :: PE_slice(:,:)
   integer, allocatable :: master_worldrank(:) 

   integer, allocatable :: evb_begin(:), evb_end(:), lpimd_dcrypt(:,:)

   _REAL_ :: nbead_inv

   _REAL_ , allocatable :: pie(:,:), piq(:,:), pif(:,:,:) &
                         , nrg_bead(:), fpimd(:,:) &
                         , evb_mat_bead(:,:,:) &
                         , evb_vec0_bead(:,:), vel0_bead(:)

   logical :: evb_pimd_ready = .false.    !! already allocated (ready == .true.)

   integer :: jobs_per_node, npimd_nodes

   integer :: nslice_per_group, nslice_per_node, lpimd_group

   contains

! %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
! %|  Initialize coupling to PIMD                                            |%
! %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%

   subroutine evb_pimd_init 

   use evb_parm,  only: C_evb, nmorse, nmodvdw
   use evb_amber, only: morse, k_harm, r0_harm, vdw_mod
   use evb_nml,   only: morsify, modvdw
   use evb_check, only: full_evb_debug
#ifdef LES
   use evb_parm,  only: inc_bond_RC, inc_dbonds_RC
   use miller,    only: i_qi, gradRC, gradRC_norm
#endif

   implicit none 
#ifdef MPI
#  include "parallel.h"
   include 'mpif.h'
#endif
#  include "../include/memory.h"
#  include "extra.h"

   !  .........................................................................

   integer :: n, m, mm, nn, alloc_error, dealloc_error
   integer :: i, j
#ifdef MPI
   integer :: ierr
#endif

!  +---------------------------------------------------------------------------+
!  |  Determine the worldrank of the masters                                   |
!  +---------------------------------------------------------------------------+

#ifdef MPI
   call mpi_bcast ( mastersize, 1, MPI_INTEGER, 0, commworld, ierr )

   allocate( master_worldrank(mastersize), stat = alloc_error )

   REQUIRE( alloc_error == 0 )

   if( master ) then

      call mpi_gather ( worldrank, 1, MPI_INTEGER, master_worldrank &
                      , 1, MPI_INTEGER, 0, commmaster, ierr )

   endif
   call mpi_bcast ( master_worldrank, mastersize, MPI_INTEGER, 0, commworld, ierr )
#endif

!  +---------------------------------------------------------------------------+
!  |  Allocate arrays for EVB-PIMD                                             |
!  :...........................................................................:
!  |  natomCL  :: number of atoms in each bead (includes cnum == 1 )           |
!  |  natomPCL :: number of pure classical atoms (excludes cnum == 1 )         |
!  |  natomPQM :: number of quantized atoms                                    |
!  +---------------------------------------------------------------------------+

   nbead = ncopy
   natomCL = 0

   do n = 1, natom
      if( cnum(n) == 0 .or. cnum(n) == 1 ) natomCL = natomCL + 1
   enddo

   natomPCL = 0
   do n = 1, natom
      if( cnum(n) == 0 ) natomPCL = natomPCL + 1
   enddo

   natomPQM = natom - natomPCL

   if( .not. evb_pimd_ready ) then
      call evb_pimd_alloc
      evb_pimd_ready = .true.
   endif
   if( full_evb_debug ) &
      write(6,'(A)') '| done allocating for evb_pimd'

!  +---------------------------------------------------------------------------+
!  |  Scale C_evb for PIMD; scaling of the exchange term is done in evb_matrix |
!  +---------------------------------------------------------------------------+

   nbead_inv = 1.0d0 / dble( nbead )

   write(6,'(/)')
   write(6,'(A)') '| Initializing EVB-PIMD: scaling the diabatic ' &
               // 'energy shifts'
   write(6, '( A, (8F10.5, 2X) )' ) '| OLD C_evb = ', C_evb(:)
   C_evb(:) = C_evb(:) * nbead_inv 
   write(6, '( A, (8F10.5, 2X) )' ) '| NEW C_evb = ', C_evb(:)

!  +---------------------------------------------------------------------------+
!  |  Set index array for distributing the PIMD jobs over CPUs                 |
!  +---------------------------------------------------------------------------+
!ifdef MPI
!  if( nbead > worldsize ) then
!     if( mod( nbead, worldsize ) /= 0 ) then 
!        write(6,'(A)') 'Error: the number of beads is not divisible ' &
!                     //'by the number of processors'
!        write(6,'(A,I8)') 'nbead = ', nbead
!        write(6,'(A,I8)') 'no. cpu  = ', worldsize
!        call mexit(6,1)
!     endif 
!     jobs_per_node = nbead / worldsize
!     nsize = worldsize 
!  else
!     jobs_per_node = 1
!     nsize = nbead
!  endif
!endif

!  allocate( evb_begin(nsize), evb_end(nsize), stat = alloc_error )
!  REQUIRE( alloc_error == 0 )

!  do n = 1, nsize
!     evb_begin(n) = ( n - 1 ) * jobs_per_node + 1
!     evb_end  (n) = n * jobs_per_node
!  enddo

!  allocate( PE_slice(nsize,jobs_per_node), stat = alloc_error )
!  REQUIRE( alloc_error == 0 )

!  do m = 1, nsize
!     do n = evb_begin(m), evb_end(m)
!        PE_slice(m,n) = n
!     enddo
!  enddo

!  do n = 1, jobs_per_node
!     do m = 1, nsize
!        PE_slice(m,n) = evb_begin(m) + n - 1
!     enddo
!  enddo

   write(6, '( A, I8 )' ) '| nbead         = ', nbead
   write(6, '( A, I8 )' ) '| natom         = ', natom
   write(6, '( A, I8 )' ) '| natomCL       = ', natomCL
   write(6, '( A, I8 )' ) '| natomPCL      = ', natomPCL
   write(6, '( A, I8 )' ) '| natomPQM      = ', natomPQM
#ifdef MPI
   write(6, '( A, I8 )' ) '| worldsize     = ', worldsize
   write(6, '( A, I8 )' ) '| jobs_per_node = ', jobs_per_node
#endif
   write(6,'(A)') '|'
   write(6, '( A, 32(I4, 2X) )' ) '| evb_begin = ', evb_begin(:)
   write(6, '( A, 32(I4, 2X) )' ) '| evb_end   = ', evb_end(:) 
   write(6, '( A, 32(I4, 2X) )' ) '| nslice/node = ', nslice_per_node

#ifdef MPI
   do n = 1, nslice_per_node
      write(6,'( A,I3,A,34(I4, 2X) )') '| lpimd_dcrypt(1:lpimd_size,', n, ') = ' &
             , lpimd_dcrypt(:,n), lpimd_rank, worldrank
   enddo
#endif

   do n = 1, jobs_per_node
      do m = 1, npimd_nodes
         write(6,'(3(A,I4))') 'PE_slice(', m, ',', n, ') = ', PE_slice(m,n)
      enddo
   enddo

!  +---------------------------------------------------------------------------+
!  |  Map the LES atoms to the corresponding PIMD slices                       |
!  +---------------------------------------------------------------------------+

   do n = 1, nbead
      mm = 0
      do m = 1, natom
         if( cnum(m) == 0 .or. cnum(m) == n ) then
            mm = mm + 1
            bead_dcrypt(mm,n) = m
         endif
      enddo
      if( mm /= natomCL) then
         write(6,'(A)') 'ERROR:: max index for bead_dcrypt(:,n) /= natomCL'
      endif
   enddo

   mm = 0
   nn = 0
   do m = 1, natom
      if( cnum(m) == 0 ) then
         mm = mm + 1
         atomCL_dcrypt(mm) = m
      else
         nn = nn + 1
         atomQM_dcrypt(nn,1) = m
         atomQM_dcrypt(nn,2) = cnum(m)
      endif
   enddo

!  do n = 1, nbead 
!     do mm = 1, natomCL
!        write(6,*) bead_dcrypt(mm,n), mm, n
!     enddo
!  enddo

!  +---------------------------------------------------------------------------+
!  |  Re-map RC indices for PIMD and QI                                        |
!  +---------------------------------------------------------------------------+

#ifdef LES

   div_ndx(1) = ncopy
   div_ndx(2) = ncopy / 2
   select case( trim( adjustl( evb_dyn ) ) )
      case( "qi_bond_pmf" )
         do n = 1, nbias
            bond_RC(n)%iatom = bead_dcrypt( bond_RC(n)%iatom, div_ndx(n) )
            bond_RC(n)%jatom = bead_dcrypt( bond_RC(n)%jatom, div_ndx(n) )
         enddo
         i_qi = 1
         allocate( gradRC(3*natomCL,ndiv), stat = alloc_error )
         allocate( gradRC_norm(ndiv), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
      case( "qi_bond_dyn" )
         do n = 1, nbias
            bond_RC(n)%iatom = bead_dcrypt( bond_RC(n)%iatom, div_ndx(n) )
            bond_RC(n)%jatom = bead_dcrypt( bond_RC(n)%jatom, div_ndx(n) )
         enddo
         i_qi = 2
         allocate( gradRC(3*natomCL,ndiv), stat = alloc_error )
         allocate( gradRC_norm(ndiv), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
      case( "qi_dbonds_pmf" )
         do n = 1, nbias
            dbonds_RC(n)%iatom = bead_dcrypt( dbonds_RC(n)%iatom, div_ndx(n) )
            dbonds_RC(n)%jatom = bead_dcrypt( dbonds_RC(n)%jatom, div_ndx(n) )
            dbonds_RC(n)%katom = bead_dcrypt( dbonds_RC(n)%katom, div_ndx(n) )
         enddo
         i_qi = 1
         allocate( gradRC(3*natomCL,ndiv), stat = alloc_error )
         allocate( gradRC_norm(ndiv), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
      case( "qi_dbonds_dyn" )
         do n = 1, nbias
            dbonds_RC(n)%iatom = bead_dcrypt( dbonds_RC(n)%iatom, div_ndx(n) )
            dbonds_RC(n)%jatom = bead_dcrypt( dbonds_RC(n)%jatom, div_ndx(n) )
            dbonds_RC(n)%katom = bead_dcrypt( dbonds_RC(n)%katom, div_ndx(n) )
         enddo
         i_qi = 2
         allocate( gradRC(3*natomCL,ndiv), stat = alloc_error )
         allocate( gradRC_norm(ndiv), stat = alloc_error )
         REQUIRE( alloc_error == 0 )
      case( "bond_umb" )
!3132009 if( ipimd == 2 ) then
         if( ipimd == 2 .or. ipimd == 3 ) then
!           write(6,*) '<<< iatom, jatom = ', bond_RC(1)%iatom, bond_RC(1)%jatom
            do n = 1, nbias
               bond_RC(n)%iatom = bead_dcrypt( bond_RC(n)%iatom, 1 )
               bond_RC(n)%jatom = bead_dcrypt( bond_RC(n)%jatom, 1 )
            enddo
!           write(6,*) '>>> iatom, jatom = ', bond_RC(1)%iatom, bond_RC(1)%jatom
         endif

      case( "dbonds_umb" )
!3132009 if( ipimd == 2 ) then
         if( ipimd == 2 .or. ipimd == 3 ) then
!           write(6,*) '<<< iatom, jatom, katom = ', dbonds_RC(1)%iatom &
!                    , dbonds_RC(1)%jatom, dbonds_RC(1)%katom
            do n = 1, nbias
               dbonds_RC(n)%iatom = bead_dcrypt( dbonds_RC(n)%iatom, 1 )
               dbonds_RC(n)%jatom = bead_dcrypt( dbonds_RC(n)%jatom, 1 )
               dbonds_RC(n)%katom = bead_dcrypt( dbonds_RC(n)%katom, 1 )
            enddo
!           write(6,*) '>>> iatom, jatom, katom = ', dbonds_RC(1)%iatom &
!                    , dbonds_RC(1)%jatom, dbonds_RC(1)%katom
         endif

      case( "groundstate" )

         if( inc_dbonds_RC ) then
!3132009 if( ipimd == 2 ) then
         if( ipimd == 2 .or. ipimd == 3 ) then
            do n = 1, nbias
               dbonds_RC(n)%iatom = bead_dcrypt( dbonds_RC(n)%iatom, 1 )
               dbonds_RC(n)%jatom = bead_dcrypt( dbonds_RC(n)%jatom, 1 )
               dbonds_RC(n)%katom = bead_dcrypt( dbonds_RC(n)%katom, 1 )
            enddo

         endif
         endif

         if( inc_bond_RC ) then
!3132009 if( ipimd == 2 ) then
         if( ipimd == 2 .or. ipimd == 3 ) then
            do n = 1, nbias
               bond_RC(n)%iatom = bead_dcrypt( bond_RC(n)%iatom, 1 )
               bond_RC(n)%jatom = bead_dcrypt( bond_RC(n)%jatom, 1 )
            enddo
         endif
         endif

   end select

#endif /* LES */

!  +---------------------------------------------------------------------------+
!  |  Replicate the morsification based on the first PIMD slice                |
!  +---------------------------------------------------------------------------+

!goto 888

   if( nmorse > 0 ) then
      deallocate( morse, k_harm, r0_harm, stat = dealloc_error )
      REQUIRE( dealloc_error == 0 )

      allocate( morse(nmorse*nbead), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
      allocate( k_harm(nmorse*nbead), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
      allocate( r0_harm(nmorse*nbead), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
   endif

   mm = 0
   do n = 1, nmorse
      i = morsify(n)%iatom
      j = morsify(n)%jatom
      write(6,'(A,2(2X,I8))') 'morsify root = ', i, j
      do m = 1, nbead
         mm = mm + 1
         morse(mm)%iatom = bead_dcrypt(i,m)
         morse(mm)%jatom = bead_dcrypt(j,m)
         write(6,'(A,2(2X,I8))') 'bead_dcrypt = ', bead_dcrypt(i,m) &
                                                 , bead_dcrypt(j,m)
         morse(mm)%D  = morsify(n)%D / dble(ncopy)
         morse(mm)%a  = morsify(n)%a
         morse(mm)%r0 = morsify(n)%r0
      enddo
   enddo

   nmorse = nmorse * nbead

   write(6,'(A)') 'MORSIFIED BONDS'
   do mm = 1, nmorse 
      write(6,'(3(2X,I8),3(2X,F14.8))') mm, morse(mm)%iatom, morse(mm)%jatom &
                                      , morse(mm)%D, morse(mm)%a, morse(mm)%r0
   enddo

!  +---------------------------------------------------------------------------+
!  |  Replicate the VDW exclusion based on the first PIMD slice                |
!  +---------------------------------------------------------------------------+

   if( nmodvdw > 0 ) then
      deallocate( vdw_mod, stat = dealloc_error )
      REQUIRE( dealloc_error == 0 )

      allocate( vdw_mod(nmodvdw*nbead), stat = alloc_error )
      REQUIRE( alloc_error == 0 )
   endif

   mm = 0
   do n = 1, nmodvdw
      i = modvdw(n)%iatom
      j = modvdw(n)%jatom
      write(6,'(A,2(2X,I8))') 'mod_vdw root = ', i, j
      do m = 1, nbead
         mm = mm + 1
         vdw_mod(mm)%iatom = bead_dcrypt(i,m)
         vdw_mod(mm)%jatom = bead_dcrypt(j,m)
         write(6,'(A,2(2X,I8))') 'bead_dcrypt = ', bead_dcrypt(i,m), bead_dcrypt(j,m)
      enddo
   enddo

   nmodvdw = nmodvdw * nbead

   write(6,'(A)') 'VDW EXCLUSIONS'
   do mm = 1, nmodvdw
      write(6,'(3(2X,I8))') mm, vdw_mod(mm)%iatom, vdw_mod(mm)%jatom
   enddo

!888 continue

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
 
   end subroutine evb_pimd_init


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Allocate storage space for EVB-PIMD                                    |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine evb_pimd_alloc

   use evb_parm, only: nevb

   implicit none

#  include "../include/memory.h"

   !  ..........................................................................

   integer :: alloc_error

  
   allocate( pie(nbead,nevb), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( piq(3,natom), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( pif(3,natom,nevb), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( nrg_bead(nbead), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( fpimd(3,natom), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( evb_mat_bead(nevb,nevb,nbead), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( evb_vec0_bead(nevb,nbead), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( vel0_bead(nbead), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( bead_dcrypt(natomCL,nbead), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( atomCL_dcrypt(natomPCL), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

   allocate( atomQM_dcrypt(natomPQM,2), stat = alloc_error )
   REQUIRE( alloc_error == 0 )

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_pimd_alloc


! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  Deallocate storage space for EVB-PIMD                                  |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine evb_pimd_dealloc

   use evb_check, only: full_evb_debug

   implicit none

   !  ..........................................................................

   integer :: dealloc_error

   deallocate( pie, piq, pif, nrg_bead, fpimd, evb_mat_bead &
             , bead_dcrypt, stat = dealloc_error )
   REQUIRE( dealloc_error == 0 )

   if( full_evb_debug ) &
      write(6,'(A)') '| DONE deallocating for evb_pimd'

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine evb_pimd_dealloc

   end module evb_pimd


