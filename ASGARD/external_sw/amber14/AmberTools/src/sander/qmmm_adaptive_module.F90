#include "../include/dprec.fh"
module qmmm_adaptive_module
! ----------------------------------------------------------------
! Adaptive QM/MM implementation, for the variable qm treatment of
! solvent molecules
!
! Reference:
! Bulo, Sikkema, Ensing, Visscher, JCTC, 2009, 5, 2212
!
! Similar to the variable solvent implementation by Ross Walker,
! with as main difference that no force jumps occur throughout the
! simulations (forces are continuous).
!
! This is achieved by introducing a transition region within which
! solvent molecules gradually aquire/loose QM character. A sigmoid
! switching function determines the QM character within the
! transition region.
!
! The AMBER implementation uses multisander
! Each group calculates the forces for its partitioning
! The master thread performs the force mixing and sends the resulting
! adaptive QM/MM force to the slave groups
! Propagation is then done by each group
!
! A corrected adaptive QM/MM energy can be calculated along the
! trajectory. This may require re-calculating QM/MM energies
! for QM/MM partitionings that (dis)appear during the course
! of the simulation.
!
! Implementation by
! Andreas Goetz (SDSC)
! Rosa Bulo (VU University Amsterdam)
! Kyoyeon Park (UCSD)
! 
! Date: January 2011 to October 2011
!
! The module has the following public routines.
!
!  adaptive_qmmm()               :  Called on each thread.
!                                   Main driver routine
!
! Everything else is private to this module.
! Check the comments in the code for details.
!
!  Following routines operate on the data type defined in this module:
!
!  new()              : Allocate arrays
!  delete()           : deallocate arrays (allocation should be done within this module)
!                       THIS GOT LOST SOMEWHERE... so arrays are currently not deallocated
!
! 
! NOTE:
!
! This module uses and modifies data in qmmm_vsolv_type (qmmm_vsolv_module)
!
!
! NOTES ABOUT MPI RANKS AND COMMUNICATORS:
!
! mastersize : The number of masters (sander copies that are running)
!              runs from 0 to (number of copies - 1)
!
! commaster  : The communicator for the master threads
!
! In our implementation, the master with ID 0 will do all the
! force mixing. The code for this master to receive, for example,
! an integer from all other masters and store it in an array awg(:)
! according to master ID, would be:
!
! if ( masterrank == 0 ) then
!         ! Master of first group needs to receive data                               
!         do i = 1, mastersize-1
!            call mpi_recv(awg(i), 1, mpi_integer, &
!                          i, 0, commmaster, ist, ierr)
!         end do
! else if ( masterrank > 0 ) then
!         ! Masters of other groups need to send data to master of group 0            
!         call mpi_send(my_int, 1, mpi_integer, &
!                       0, 0, commmaster, ist, ierr)
! end if
!
! Master thread with ID 0 runs the permanent QM system
! Master thread with ID 1 runs the permanent QM + 1st adaptive QM solvent system
! ...
! Master thread with ID N runs the permanent QM + Nth adaptive QM solvent system
!
! ----------------------------------------------------------------

  implicit none

  private

  public :: adaptive_qmmm
  public :: qmmm_adaptive_type
  public :: adaptive_reset

  type qmmm_adaptive_type

     logical :: debug

     integer :: verbosity

     ! flag indicating whether the energy correction term Wbk shall
     ! be calculated. 
     ! calc_wbk=0 : do not calculate Wbk (default)
     ! calc_wbk=1 : use two-timestep approximation
     ! calc_wbk=2 : use published equation, 3-point central difference
     integer :: calc_wbk

     ! bookkeeping term 
     _REAL_ :: wbk

!     ! A specific QM center is used to determine the distance of the
!     ! solvent molecules from the QM region. 
!     ! This is required for the distance to be a smooth function of
!     ! the coordinates.
!     ! Eventually, we may also use the center of weight but this is not
!     ! implemented yet.
!     integer :: qm_center_id

     ! The number of partitions in the adaptive QM/MM scheme
     ! QM region: fixed QM region + qmmm_vsolv_type%nearest_qm_solvent
     ! T region : transition region 
     ! 1 = QM region only (not adaptive, no T region)
     ! 2 = QM region + 1 solvent molecule in T region
     ! 3 = QM region + 2 solvent moleculues in T region
     ! etc
     ! This has to be identical to the number of master threads
     integer :: n_partition

     ! RA : outer radius of active region
     ! RT : outer radius of transition region
     _REAL_ :: RA, RT

     ! print_qm_coords : whether qm_coords should be printed
     logical :: print_qm_coords

     ! Additional information (which atoms are fixed, and the QM center)
     ! are stored in the qmmm_vsolv_type data type, see qmmm_vsolv_module
     ! There are n_partition weights:
     !   weights(1) : QM region
     !   weights(2) : QM region + 1 solvent molecule in T region
     ! We store (n_partition + 1) solvent IDs and distances, that is, ts_solvent_pointers()
     ! and ts_solvent_distances() contains not only information about the solvent molecules
     ! in the T region, but also the closest solvent molecules in the QM region and
     ! the MM region:
     ! ts_solvent_pointers(1) = outermost solvent in QM region
     ! ts_solvent_pointers(2) = first solvent in T region
     ! ...
     ! ts_solvent_pointers(n_partition+1) = first sovlent in MM region
     integer, dimension(:), pointer :: ts_solvent_pointers => null()
     _REAL_ , dimension(:), pointer :: ts_solvent_distances => null()
     _REAL_ , dimension(:), pointer :: weights => null()

     ! Store data from previous steps to correct the potential energy
     _REAL_ , dimension(:,:), pointer :: weights_prev => null()
     integer, dimension(:,:), pointer :: ts_solvent_pointers_prev => null()
     integer, dimension(:,:), pointer :: nearest_solvent_pointers_prev => null()
     _REAL_ , dimension(:), pointer :: energies => null()
     _REAL_ , dimension(:), pointer :: energies_prev => null()
     _REAL_ , dimension(:), pointer :: energies_recalculated => null()

     integer, dimension(:), pointer :: match_partition => null()
     integer, dimension(:), pointer :: calc_partition => null()

     _REAL_, dimension(:,:), pointer :: coords_prev => null()

  end type qmmm_adaptive_type

  logical, save :: first_call = .true.

#ifdef MPI
  interface new
     module procedure new_qmmm_adaptive_type
  end interface

  interface delete
     module procedure delete_qmmm_adaptive_type
  end interface

  interface broadcast
     module procedure broadcast_qmmm_adaptive_type
  end interface

#endif

contains

#ifdef MPI
  ! All work related to adaptive QM/MM is done only by the 
  ! group master threads (masterrank >= 0)
  ! Only for an additional call to force(), which may be required for
  ! the energy correction term Wbk, also the slave threads have to do work.
  ! Only the "main master" (masterrank == 0) does the final force mixing and
  ! other global work such as determining which partitionings that may have (dis)appeared
  ! wrt previous time steps as required for the calculation of the energy correction
  subroutine adaptive_qmmm (nstep, natom, unimaged_coords,&
                            force, poten, ntpr, ntwx, &
                            xx, ix, ih, ipairs, qsetup, do_list_update, &
                            corrected_energy, aqmmm_flag)

     use qmmm_module, only: qmmm_vsolv
     implicit none

#include "parallel.h"
     include 'mpif.h'

     integer, intent(in)    :: nstep
     integer, intent(in)    :: natom
     _REAL_,  intent(in)    :: unimaged_coords(3,natom)
     _REAL_,  intent(inout) :: force(3,natom)
     _REAL_,  intent(in)    :: poten
     integer, intent(in)    :: ntpr
     integer, intent(in)    :: ntwx

     ! required for call to force() in do_recalculate()
     _REAL_,           intent(in)    :: xx(*)
     integer,          intent(in)    :: ix(*)
     character(len=4), intent(in)    :: ih(*)
     integer,          intent(in)    :: ipairs(*)
     logical,          intent(inout) :: qsetup, do_list_update
     _REAL_,           intent(out)   :: corrected_energy
     integer,          intent(out)   :: aqmmm_flag

     ! Local variables
     integer :: atom_id1, atom_id2, ierr
     integer :: recalc_type
     type(qmmm_adaptive_type), save :: self
     _REAL_ :: coords(3,natom)   ! FIXME could use a pointer here, but pointee would need target attribute
     _REAL_ :: distance1, distance2
     _REAL_ :: poten_recalculated
     _REAL_, allocatable :: force_backup(:,:)


     ! only main master keeps all adaptive data
     ! only general adqmm namelist settings are broadcast
     if (first_call) then
        first_call = .false.
        call setup_adaptive(self, natom)
        if ( ( self%debug .or. (self%verbosity > 1) ) .and. ( masterrank >=0 ) ) then
           write(6,'(a)') 'adQMMM setup done'
        end if
     end if
     
     if (self%debug) then
        write(6,'(a)') '>>>>> entered adaptive_qmmm (qmmm_adaptive_module)'
        call flush(6)
     end if

     aqmmm_flag = self%calc_wbk

     ! all masters of each group:
     ! determine outermost atom IDs
     if ( masterrank >= 0 ) then
        call get_qm_and_mm_boundaries(self, atom_id1, distance1, &
                                     atom_id2, distance2)  
        if ( self%debug .or. (self%verbosity > 1) ) then
           write(6,'(a,i3)')    'I am the master of group ', masterrank
           write(6,'(a,i6)')    'Outermost atom_id  = ', atom_id1
           write(6,'(a,f12.6)') 'Outermost distance = ', dsqrt(distance1)
        end if
     end if

     if ( masterrank > 0 ) then
        call send_atid_dist(self, atom_id1, distance1, atom_id2, distance2)
     else if ( masterrank == 0 ) then
        ! Store the local id's, distances 
        self%ts_solvent_pointers(1)  = atom_id1
        self%ts_solvent_distances(1) = distance1
        ! Receive data from all groups
        call receive_atid_dist(self)
     end if

     if (masterrank == 0) then
        self%ts_solvent_distances(:) = dsqrt(self%ts_solvent_distances(:))
        call check_partition(self)
        call calculate_weights(self,nstep)
     end if

     if (masterrank > 0) then
        ! All groups send their forces to the main master
        call send_forces(self, natom, force)
     else if (masterrank == 0) then
        ! The main master receives forces from all groups and mixes the forces
        call receive_and_mix_forces(self, natom, force)
     end if

     if (self%calc_wbk > 0) then

        ! All groups:
        if (masterrank > 0) then

           ! Send potential energy to the main master
           call send_energy(self, poten)
           ! receive flag whether energy needs to be recalculated
           ! (either with different partitioning or different geometry)
           if (nstep+1 > self%calc_wbk) then
              call receive_recalculate(self, qmmm_vsolv%recalculate, recalc_type)
           end if
           
        ! Main master:
        else if (masterrank == 0) then

           ! Receive potential energies from all groups
           self%energies(1) = poten
           call receive_energies(self,0)
           ! Determine whether the energy has to be recalculated for any partitioning
           ! and send result to group masters
           if (nstep+1 > self%calc_wbk) then
              call check_match_partition(self, qmmm_vsolv%nearest_solvent_pointers, & 
                                         qmmm_vsolv%nearest_qm_solvent)
              call check_and_send_recalculate(self, qmmm_vsolv%recalculate, recalc_type)
           end if

        end if
    
        ! Calculate missing energy partitionings
        ! First broadcast recalculate info to all slaves within each group
        call mpi_bcast(qmmm_vsolv%recalculate, 1, MPI_LOGICAL, 0, commsander, ierr)
        call mpi_bcast(recalc_type, 1, MPI_INTEGER, 0, commsander, ierr)

        if (qmmm_vsolv%recalculate) then
           ! Recalculate energy for the previous partitioning
           ! Note: This has to be done by *all* threads in a group for which 
           !       recalculate is .true.
           if (self%calc_wbk == 1) then
              coords(1:3,1:natom) = unimaged_coords(1:3,1:natom)
           else 
              coords(1:3,1:natom) = self%coords_prev(1:3,1:natom)
           endif

           ! make backup of our mixed forces (currently held by master)
           if (masterrank == 0) then
              allocate(force_backup(3,natom), stat=ierr)
              if (ierr /= 0) then
                 call sander_bomb('adaptive_qmmm (qmmm_adaptive_module)', &
                      'Allocation error for force_backup(:,:)', &
                      'Will quit now')
              end if
              force_backup = force
           end if

           call do_recalculate(poten_recalculated, masterrank, recalc_type, &
                qmmm_vsolv%nearest_qm_solvent, self%nearest_solvent_pointers_prev, &
                qmmm_vsolv%nearest_solvent_pointers, &
                xx, ix, ih, ipairs, coords, &
                qsetup, do_list_update, &
                self%debug)

           ! restore force backup (will be distributed to other groups below)
           if (masterrank == 0) then
              force = force_backup
              deallocate(force_backup, stat=ierr)
              if (ierr /= 0) then
                 call sander_bomb('adaptive_qmmm (qmmm_adaptive_module)', &
                      'Deallocation error for force_backup(:,:)', &
                      'Will quit now')
              end if
           end if

        end if

        ! reset qmmm_vsolv%recalculate to .false. on all threads
        qmmm_vsolv%recalculate = .false.
    
        ! Group masters:
        if ( masterrank > 0 ) then

           ! Send recalculated potential energies to the main master
           call send_energy(self, poten_recalculated)

        ! Main master:
        else if (masterrank == 0) then

           ! Receive recalculated potential energies from all groups  
           self%energies_recalculated(:) = 0.0d0
           self%energies_recalculated(1) = poten_recalculated
           call receive_energies(self,1)

           if (nstep+1 > self%calc_wbk) then
           ! Calculate missing contributions to wbk
              call calc_corrected_energy(self, corrected_energy)
           end if

        end if

        ! backup coordinates from this step
        if ( self%calc_wbk > 1 ) then
           self%coords_prev = unimaged_coords
        end if

     endif ! self%calc_wbk

     ! Main master:
     if (masterrank == 0) then

        ! Send mixed forces
        call send_mixed_forces(self, natom, force)

     ! All threads:
     ! (each thread works only on the subset of forces of the atoms it owns
     !  but it's easier to simply broadcast everything)
     else

        ! Receive mixed forces
        call receive_mixed_forces(self, natom,force)

     end if

     ! All masters:
     if (masterrank >= 0) then

        ! Store data of relevance for adQMMM
        ! This is presently done on all master threads, although only
        ! the main master has/uses all information
        call store_present_data(self, qmmm_vsolv%nearest_solvent_pointers)

     end if

     ! print coordinates for the QM region
     if ( (masterrank >= 0) .and. self%print_qm_coords ) then
        if ( mod(nstep,ntwx) == 0 )then
           call print_qm_coords(nstep, masterrank, unimaged_coords, natom)
        end if
     end if 

     if ( masterrank == 0 ) then
        if ( self%verbosity > 0 ) then
           call print_adqmmm_info(self, nstep, ntpr)
        end if
     end if
 
     if (self%debug) then
        write(6,'(a)') '<<<<< leaving adaptive_qmmm (qmmm_adaptive_module)'
        call flush(6)
     end if

  end subroutine adaptive_qmmm


  ! Calculate the adaptive QM/MM energy
  ! E(ad) = sum_i sigma_i E_i - Wbk
  ! where Wbk is the book-keeping correction along the trajectory
  ! Wbk may need re-calculation of the energy for partitionings
  ! that (dis)appear along the trajctory.
  subroutine calc_corrected_energy(self, corrected_energy)

    implicit none

    type(qmmm_adaptive_type), intent(inout) :: self
    _REAL_, intent(out) :: corrected_energy

    ! Local variables
    _REAL_ :: mixed_energy
    integer :: i

    if (self%debug) then
       write(6,'(a)') '>>>>> entered calc_corrected_energy (qmmm_adaptive_module)'
       call flush(6)
    end if

    mixed_energy = 0.0d0

    do i=1,self%n_partition
       if (self%calc_wbk == 1) then
          mixed_energy = mixed_energy + self%weights(i)*self%energies(i)
       else 
          mixed_energy = mixed_energy + self%weights_prev(i,1)*self%energies_prev(i)
       end if
    end do

    call calculate_wbk(self)

    corrected_energy = mixed_energy - self%wbk

    if (self%debug) then
       write(6,'(a,f20.5)')'mixed energy    ',mixed_energy
       write(6,'(a,f20.5)')'corrected energy',corrected_energy
       write(6,'(a)') '<<<<< leaving calc_corrected_energy (qmmm_adaptive_module)'
       call flush(6)
    end if

  end subroutine calc_corrected_energy


  ! Calculate Wbk, the book-keeping correction along the trajectory
  ! Wbk needs re-calculated energies for partitionings
  ! that may have (dis)appeared along the trajctory.
  subroutine calculate_wbk(self)

    implicit none

    type(qmmm_adaptive_type), intent(inout) :: self

    ! Local variable
    _REAL_ :: wbk_tmp
    integer :: i

    if (self%debug) then
       write(6,'(a)') '>>>>> entered calculate_wbk (qmmm_adaptive_module)'
       call flush(6)
    end if

    wbk_tmp = 0.0d0

    if (self%calc_wbk == 1) then
       do i=1,self%n_partition
          if (self%match_partition(i) == 1) then
             wbk_tmp = wbk_tmp + (self%weights(i) - self%weights_prev(i,1)) * self%energies(i)
          else
             wbk_tmp = wbk_tmp + self%weights(i) * self%energies(i)
             wbk_tmp = wbk_tmp - self%weights_prev(i,1) * self%energies_recalculated(i)
          end if
       end do
    else if (self%calc_wbk == 2) then
       do i=1,self%n_partition
          if (self%match_partition(i) == 1) then
             wbk_tmp = wbk_tmp + (self%weights(i) - self%weights_prev(i,2)) &
                                  * self%energies_prev(i) / 2.0d0
          else if (self%match_partition(i) == 2) then
             wbk_tmp = wbk_tmp - self%weights_prev(i,2) * self%energies_prev(i) / 2.0d0
             wbk_tmp = wbk_tmp + self%weights(i) * self%energies_recalculated(i) / 2.0d0
          else if (self%match_partition(i) == 3) then
             wbk_tmp = wbk_tmp + self%weights(i) * self%energies_prev(i) / 2.0d0
             wbk_tmp = wbk_tmp - self%weights_prev(i,2) * self%energies_recalculated(i) / 2.0d0
          else if (self%match_partition(i) == 4) then
             wbk_tmp = wbk_tmp + (self%weights(i) - self%weights_prev(i,2)) &
                                  * self%energies_recalculated(i) / 2.0d0
          end if
       end do
    else
       call sander_bomb('calculate_wbk (qmmm_adaptive_module)', &
            'Wrong option for calc_wbk',&
            'bye bye!')
    end if

    self%wbk = self%wbk + wbk_tmp

    if (self%debug) then
       write(6,'(a,f10.3)')'wbk term : ', self%wbk
       write(6,'(a)') '<<<<< leaving calculate_wbk (qmmm_adaptive_module)'
       call flush(6)
    end if


  end subroutine calculate_wbk


  ! Recalculate the energy with the previous partitioning
  ! Output is only the potential energy
  ! (forces are not needed)
  subroutine do_recalculate(poten_recalculated, masterrank, recalc_type, &
       nearest_qm_solvent, nearest_solvent_pointers_prev, nearest_solvent_pointers, &
       xx, ix, ih, ipairs, coords, &
       qsetup, do_list_update, &
       debug)

    use state, only: state_rec

    implicit none
#include "../include/memory.h"

    _REAL_, intent(out) :: poten_recalculated
    integer, intent(in) :: masterrank
    integer, intent(in) :: recalc_type
    integer, intent(in) :: nearest_qm_solvent
    integer, intent(in) :: nearest_solvent_pointers_prev(:,:)
    integer, intent(inout) :: nearest_solvent_pointers(:)

    _REAL_, intent(in) :: xx(*)
    integer, intent(in) :: ix(*)
    character(len=4), intent(in) :: ih(*)
    integer, intent(in) :: ipairs(*)
     _REAL_, intent(in) :: coords(3,natom) ! natom defined in "../include/memory.h"
    logical, intent(inout) :: qsetup, do_list_update

    logical, intent(in) :: debug

    ! LOCAL
    type(state_rec) :: ener_scratch
    integer :: nearest_solvent_pointers_backup(nearest_qm_solvent+1)
    
    if (masterrank > -1) then

       if (debug) then
          write(6,'(a)') '>>>>> entered subroutine do_recalculate (qmmm_adaptive_module)'
          write(6,*) 'I AM MASTER WITH ID ', masterrank
          write(6,*) 'nearest_qm_solvent = ', nearest_qm_solvent
          call flush(6)
       end if

       ! The master of this group needs to make a backup of the present
       ! partitioning, use the partitioning from last step and later restore
       ! the present partitioning
       ! NOTE: The partitioning will be used in qmmm_vsolv_update() which is
       !       called from force() to update the QM/MM region and partitioning
       !       information also on the slaves of this group
       ! BACKUP
       nearest_solvent_pointers_backup(:) = nearest_solvent_pointers(:)
       if (recalc_type == 0) then
          ! for calc_wbk = 1
          ! use old solvent pointers from previous step
          nearest_solvent_pointers(:) = nearest_solvent_pointers_prev(:,1)
       else if (recalc_type == 2 .or. recalc_type == 4) then
          ! for calc_wbk = 2
          ! use solvent pointers from [t+1] step
          ! no need to change the solvent pointers
       else if (recalc_type == 3) then
          ! for calc_wbk = 2
          ! use solvent pointers from [t-1] step
          nearest_solvent_pointers(:) = nearest_solvent_pointers_prev(:,2)
       end if

    end if

    ! All threads of this group will recalculate the forces now
    ! The QM/MM partitioning will be set up in subroutine
    ! qmmm_vsolv_update() using the solvent pointers from the previous step
    ! x = unimaged_coords
    ! xx(lforce) = force array
    !NOTE: cannot pass an array force(3,natom) since xx(lforce) is actually larger
    !      (in parallel code dies with heap errors)
    !      (see also locmem.F90 where lforce is computed)
    call force(xx, ix, ih, ipairs, coords, &
         xx(lforce), ener_scratch, ener_scratch%vir,  &
         xx(l96), xx(l97), xx(l98), xx(l99), qsetup, &
         do_list_update, 0)

    poten_recalculated = ener_scratch%pot%tot

    if (masterrank > -1) then
       ! RESTORE
       nearest_solvent_pointers(:) = nearest_solvent_pointers_backup(:)

       if (debug) then
          write(6,'(a,i3,a,f20.10)') 'potential energy for master ID ', masterrank, &
               ' = ', poten_recalculated
          write(6,'(a)') '<<<<< leaving subroutine do_recalculate'
          call flush(6)
       end if

    end if

  end subroutine do_recalculate


  subroutine receive_recalculate(self, recalculate, recalc_type)

    implicit none
#include "parallel.h"
    include 'mpif.h'

    type(qmmm_adaptive_type), intent(in) :: self
    integer, intent(out) :: recalc_type
    logical, intent(out) :: recalculate

    integer :: ist(MPI_STATUS_SIZE), ierr

    if (self%debug) then
       write(6,'(a)') '>>>>> entered receive_recalculate (qmmm_adaptive_module)'
       call flush(6)
    end if

    call mpi_recv(recalculate, 1, MPI_LOGICAL, 0, 0, commmaster, ist, ierr)
    call mpi_recv(recalc_type, 1, MPI_INTEGER, 0, 0, commmaster, ist, ierr)

    if (self%debug) then
       write(6,'(a,l)') 'received recalculate = ', recalculate
       write(6,'(a)') '<<<<< leaving receive_recalculate (qmmm_adaptive_module)'
       call flush(6)
    end if

  end subroutine receive_recalculate


  ! Check whether the energy for a (dis)appeared partitioning needs
  ! to be recalculated
  ! Send results to group masters
  subroutine check_and_send_recalculate(self, recalculate, recalc_type)

    implicit none

#include "parallel.h"
    include 'mpif.h'

    type(qmmm_adaptive_type), intent(in) :: self
    integer, intent(out) :: recalc_type
    logical, intent(out) :: recalculate

    integer :: i
    integer :: ist(MPI_STATUS_SIZE), ierr
    integer :: recalc_type_others
    logical :: recalculate_others

    if (self%debug) then
       write(6,'(a)') '>>>>> entered check_and_send_recalculate (qmmm_adaptive_module)'
       call flush(6)
    end if

    if (self%match_partition(1) /= 1) then
       recalculate = .true.
    end if
    recalc_type=self%match_partition(1)

    do i = 2, self%n_partition
       recalculate_others = .false.
       if (self%match_partition(i) /= 1) then
          recalculate_others = .true.
       end if
       recalc_type_others=self%match_partition(i)
       call mpi_send(recalculate_others, 1, MPI_LOGICAL, i-1, 0, commmaster, ist, ierr)
       call mpi_send(recalc_type_others, 1, MPI_INTEGER, i-1, 0, commmaster, ist, ierr)
       if (self%debug) then
          write(6,'(a,l,a,i3)') 'sending recalculate = ', recalculate_others, ' to master ID = ', i-1
          call flush(6)
       end if
    end do

    if (self%debug) then
       write(6,'(a)') '<<<<< leaving check_and_send_recalculate (qmmm_adaptive_module)'
       call flush(6)
    end if

  end subroutine check_and_send_recalculate



  ! Store data of present MD step
  subroutine store_present_data(self, nearest_solvent_pointers)

    implicit none
    type(qmmm_adaptive_type), intent(inout) :: self
    integer, intent(in) :: nearest_solvent_pointers(:)
    
    if (self%debug) then
       write(6,'(a)') '>>>>> entered store_present_data (qmmm_adaptive_module)'
       call flush(6)
    end if

    ! Store weights, solvent pointers and energies
    ! for energy correction in following steps 
    self%nearest_solvent_pointers_prev(:,2) = self%nearest_solvent_pointers_prev(:,1)
    self%nearest_solvent_pointers_prev(:,1) = nearest_solvent_pointers(:)
    self%weights_prev(:,2)               = self%weights_prev(:,1)
    self%weights_prev(:,1)               = self%weights(:)
    self%ts_solvent_pointers_prev(:,2)   = self%ts_solvent_pointers_prev(:,1)
    self%ts_solvent_pointers_prev(:,1)   = self%ts_solvent_pointers(:)
    self%energies_prev(:)                = self%energies
    
    if (self%debug) then
       write(6,'(a)') '<<<<< leaving store_present_data (qmmm_adaptive_module)'
       call flush(6)
    end if

  end subroutine store_present_data
    


  ! The adaptive_data arrays are allocated here. This is called only
  ! once in the adaptive md run (the initilization phase)
  subroutine setup_adaptive (self, natom)

    use qmmm_module, only: qmmm_vsolv
    implicit none

    integer, intent(in) :: natom

#include "parallel.h"

    type(qmmm_adaptive_type), intent(out) :: self

    if ( masterrank == 0 ) then
       call get_adqmmm_namelist(self)
    end if

    call broadcast(self)

    if ( masterrank >= 0 ) then

       call print_adqmmm_namelist(self)
       if (self%debug) then
          write(6,'(a,i3)') '  >>> Running Adaptive QMMM on Group Master ID ',masterrank,' <<<'
          write(6,'(a)') '>>>>> entered in setup_adaptive (qmmm_adaptive_module)'
          call flush(6)
       end if

    end if
    
    ! Allocate
    call new(self, natom)

    if ( masterrank >= 0 ) then

       if (self%n_partition /= mastersize) then
          call sander_bomb('setup_adaptive (qmmm_adaptive_module)', &
               'The number of groups differs from n_partition.', &
               'bye bye!')
       end if

       if ( self%debug ) then
          write (6,'(a,i6)')'self%n_partition = ', self%n_partition
          write (6,'(a,i6)')'qmmm_vsolv%nearest_qm_solvent = ', qmmm_vsolv%nearest_qm_solvent
          write (6,'(a)') '<<<<< leaving setup_adaptive (qmmm_adaptive_module)'
          call flush(6)
       end if

    end if

  end subroutine setup_adaptive


  ! Read ad_qmmm namelist values from file mdin,
  ! use default values if none are present.
  subroutine get_adqmmm_namelist(self)

    use qmmm_module, only : qmmm_nml
    use file_io_dat

    implicit none

    ! Passed in
    type(qmmm_adaptive_type), intent(out) :: self

    ! Local
    _REAL_  :: RA, RT
    integer :: debug, verbosity, calc_wbk, n_partition, print_qm_coords
    integer :: ierr
    namelist /adqmmm/ debug, verbosity, calc_wbk, n_partition, RA, RT, print_qm_coords

    ! Set default values for adqmmm namelist values
    debug           = 0
    verbosity       = 0
    calc_wbk        = 0
    n_partition     = 1
    print_qm_coords = 0

    RA = -1.d0
    RT = -1.d0
    
    ! Read namelist
    rewind 5
    read(5, nml=adqmmm, iostat=ierr)

    if ( ierr > 0 ) then
       call sander_bomb('get_adqmmm_namelist (qmmm_adaptive_module)', &
            '&adqmmm namelist read error', &
            'Please check your input.')
    else if ( ierr < 0 ) then
       write(6,'(a)') '&adqmmm namelist read success'
    end if

    if (qmmm_nml%vsolv == 3) then
       if ((RA < 0.d0) .or. (RT < 0.d0)) then
          call sander_bomb('get_adqmmm_namelist (qmmm_adaptive_module)', &
            '&adqmmm namelist read error', &
            'RA and RT are not defined!')
       end if
    end if

    ! debug
    self%debug = .false.
    if (debug > 0) then
       self%debug = .true.
    end if

    ! print QM coords
    self%print_qm_coords = .false.
    if (print_qm_coords > 0) then
       self%print_qm_coords = .true.
    end if

    ! verbosity
    self%verbosity = verbosity

    ! calc_wbk
    self%calc_wbk = calc_wbk

    ! n_partition
    self%n_partition = n_partition

    ! RA and RT
    self%RA = RA
    self%RT = RT

  end subroutine get_adqmmm_namelist


  ! print namelist values
  subroutine print_adqmmm_namelist(self)

    type(qmmm_adaptive_type), intent(in) :: self

!    write(6,'(/80("-")/"  3.2 ADAPTIVE SOLVENT QM/MM CALCULATION INFO",/80("-"))')

    write (6, '(/a)') &
         '| Citations for ADAPTIVE SOLVENT QM/MM run:'

    write (6, '(4(/a))') &
         '| A. W. G"otz, K. Park, R. E. Bulo, F. Paesani, R. C. Walker', &
         '| "Efficient adaptive QM/MM implementation: Application to', &
         '|  ion binding by peptides in solution"', &
         '| in preparation.'

    write (6, '(4(/a))') &
         '| R. E. Bulo, B. Ensing, J. Sikkema, L. Visscher', &
         '| "Toward a practical method for adaptive QM/MM simulations"', &
         '| J. Chem. Theory Comput. 9 (2009) 2212-2221.', &
         '| DOI: 10.1021/ct900148e'


    write(6,'(/,a)') 'QMMM ADQMMM options: (check also QMMM VSOLV options above)'

    write(6,'(3x,a,l5)')   'debug           = ', self%debug
    write(6,'(3x,a,i5)')   'verbosity       = ', self%verbosity
    write(6,'(3x,a,l5)')   'print_qm_coords = ', self%print_qm_coords
    write(6,'(3x,a,l5)')   'calc_wbk        = ', self%print_qm_coords
    write(6,'(3x,a,i5)')   'n_partition     = ', self%n_partition
    write(6,'(3x,a,f5.2)') 'RA              = ', self%RA
    write(6,'(3x,a,f5.2)') 'RT              = ', self%RT

  end subroutine print_adqmmm_namelist


  ! Allocate and initialize arrays
  subroutine new_qmmm_adaptive_type(self, natom)

    use constants, only: zero
    use qmmm_module, only: qmmm_vsolv
    implicit none

    type(qmmm_adaptive_type), intent(inout) :: self
    integer, intent(in) :: natom

#include "parallel.h"

    integer :: ierr

    if (self%debug) then
       write(6,'(a)') '>>>>> entered new_qmmm_adaptive_type (qmmm_adaptive_module)'
       call flush(6)
    end if

    self%wbk = 0.0d0

    if ( self%calc_wbk > 1 ) then
       allocate ( self%coords_prev(3, natom), stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('new_qmmm_adaptive_type (qmmm_adaptive_module)', &
               'Allocation error for coords_prev(:,:)', &
               'Will quit now')
       end if
    end if

    ! slaves do not need storage for data below
    if ( masterrank < 0 ) return

    if ( self%n_partition > 0 ) then

       allocate ( self%weights(self%n_partition), stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('new_qmmm_adaptive_type (qmmm_adaptive_module)', &
               'Allocation error for weights(:)', &
               'Will quit now')
       end if
       self%weights = zero

       allocate ( self%weights_prev(self%n_partition,2), stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('new_qmmm_adaptive_type (qmmm_adaptive_module)', &
               'Allocation error for weights_prev(:)', &
               'Will quit now')
       end if
       self%weights_prev = zero

       allocate ( self%ts_solvent_pointers_prev(self%n_partition + 1, 2), &
                  stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('new_qmmm_adaptive_type (qmmm_adaptive_module)', &
               'Allocation error for ts_solvent_pointers_prev(:)', &
               'Will quit now')
       end if
       self%ts_solvent_pointers_prev = 0

       allocate ( self%nearest_solvent_pointers_prev(qmmm_vsolv%nearest_qm_solvent+1, 2), &
                  stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('new_qmmm_adaptive_type (qmmm_adaptive_module)', &
               'Allocation error for ts_solvent_pointers_prev(:)', &
               'Will quit now')
       end if
       self%nearest_solvent_pointers_prev = 0

       allocate ( self%energies(self%n_partition), &
                  stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('new_qmmm_adaptive_type (qmmm_adaptive_module)', &
               'Allocation error for energies(:)', &
               'Will quit now')
       end if
       self%energies = zero

       allocate ( self%energies_prev(self%n_partition), &
                  stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('new_qmmm_adaptive_type (qmmm_adaptive_module)', &
               'Allocation error for energies_prev(:)', &
               'Will quit now')
       end if
       self%energies_prev = zero

       allocate ( self%energies_recalculated(self%n_partition), &
                  stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('new_qmmm_adaptive_type (qmmm_adaptive_module)', &
               'Allocation error for energies_recalculated(:)', &
               'Will quit now')
       end if
       self%energies_recalculated = zero

       allocate ( self%ts_solvent_pointers(self%n_partition + 1), stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('new_qmmm_adaptive_type (qmmm_adaptive_module)', &
               'Allocation error for ts_solvent_pointers(:)', &
               'Will quit now')
       end if
       self%ts_solvent_pointers = 0

       allocate ( self%ts_solvent_distances(self%n_partition + 1), stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('new_qmmm_adaptive_type (qmmm_adaptive_module)', &
               'Allocation error for ts_solvent_distances(:)', &
               'Will quit now')
       end if
       self%ts_solvent_distances = zero

       allocate ( self%match_partition(self%n_partition), stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('new_qmmm_adaptive_type (qmmm_adaptive_module)', &
               'Allocation error for match_partition(:)', &
               'Will quit now')
       end if
       self%match_partition = 0

       allocate ( self%calc_partition(self%n_partition), stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('new_qmmm_adaptive_type (qmmm_adaptive_module)', &
               'Allocation error for calc_partition(:)', &
               'Will quit now')
       end if
       self%calc_partition = 0

    endif

    if (self%debug) then
       write(6,'(a)') '<<<<< leaving new_qmmm_adaptive_type (qmmm_adaptive_module)'
       call flush(6)
    end if

  end subroutine new_qmmm_adaptive_type


  ! Dellocate arrays
  subroutine delete_qmmm_adaptive_type(self)

    implicit none

    type(qmmm_adaptive_type), intent(inout) :: self

    integer :: ierr

    if ( associated(self%coords_prev) ) then
       deallocate ( self%coords_prev, stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('delete_adaptive_type (qmmm_adaptive_module)', &
               'Deallocation error for coords_prev(:,:)', &
               'Will quit now')
       end if
    end if

    if ( associated(self%weights) ) then
       deallocate ( self%weights, stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('delete_adaptive_type (qmmm_adaptive_module)', &
               'Deallocation error for weights(:)', &
               'Will quit now')
       end if
    end if

    if ( associated(self%ts_solvent_pointers) ) then
       deallocate ( self%ts_solvent_pointers, stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('delete_adaptive_type (qmmm_adaptive_module)', &
               'Deallocation error for ts_solvent_pointers(:)', &
               'Will quit now')
       end if
    end if
    
    if ( associated(self%ts_solvent_distances) ) then
       deallocate ( self%ts_solvent_distances, stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('delete_adaptive_type (qmmm_adaptive_module)', &
               'Deallocation error for ts_solvent_distances(:)', &
               'Will quit now')
       end if
    end if

    if ( associated(self%weights_prev) ) then
       deallocate ( self%weights_prev, stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('delete_adaptive_type (qmmm_adaptive_module)', &
               'Deallocation error for weights_prev(:)', &
               'Will quit now')
       end if
    end if

    if ( associated(self%ts_solvent_pointers_prev) ) then
       deallocate ( self%ts_solvent_pointers_prev, stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('delete_adaptive_type (qmmm_adaptive_module)', &
               'Deallocation error for ts_solvent_pointers_prev(:)', &
               'Will quit now')
       end if
    end if

    if ( associated(self%nearest_solvent_pointers_prev) ) then
       deallocate ( self%nearest_solvent_pointers_prev, stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('delete_adaptive_type (qmmm_adaptive_module)', &
               'Deallocation error for nearest_solvent_pointers_prev(:,:)', &
               'Will quit now')
       end if
    end if

    if ( associated(self%energies) ) then
       deallocate ( self%energies, stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('delete_adaptive_type (qmmm_adaptive_module)', &
               'Deallocation error for energies(:)', &
               'Will quit now')
       end if
    end if

    if ( associated(self%energies_prev) ) then
       deallocate ( self%energies_prev, stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('delete_adaptive_type (qmmm_adaptive_module)', &
               'Deallocation error for energies_prev(:)', &
               'Will quit now')
       end if
    end if

    if ( associated(self%energies_recalculated) ) then
       deallocate ( self%energies_recalculated, stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('delete_adaptive_type (qmmm_adaptive_module)', &
               'Deallocation error for energies_recalculated(:)', &
               'Will quit now')
       end if
    end if

    if ( associated(self%match_partition) ) then
       deallocate ( self%match_partition, stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('delete_adaptive_type (qmmm_adaptive_module)', &
               'Deallocation error for match_partition(:)', &
               'Will quit now')
       end if
    end if

    if ( associated(self%calc_partition) ) then
       deallocate ( self%calc_partition, stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('delete_adaptive_type (qmmm_adaptive_module)', &
               'Deallocation error for calc_partition(:)', &
               'Will quit now')
       end if
    end if


  end subroutine delete_qmmm_adaptive_type


  ! This subroutine gets called by the master of each group
  ! The main master has to use it later, to determine the coefficients for
  ! the adaptive QM/MM force mixing scheme.
  ! The slaves have to communicate it to the main master
  ! This subroutine determines the distance and id of the solvent QM atom with
  ! the largest distance from the QM center.
  subroutine get_qm_and_mm_boundaries (self, largest_qm_id, largest_qm_dist, &
                                      smallest_mm_id, smallest_mm_dist)

    use qmmm_module, only: qmmm_vsolv
    implicit none

    type(qmmm_adaptive_type), intent(in) :: self
    integer, intent(out) :: largest_qm_id
    integer, intent(out) :: smallest_mm_id
    _REAL_ , intent(out) :: largest_qm_dist
    _REAL_ , intent(out) :: smallest_mm_dist

    if (self%debug) then
       write(6,'(a)') '>>>>> entered get_qm_and_mm_boundaries (qmmm_adaptive_module)'
       call flush(6)
    end if

    if (qmmm_vsolv%qm_center_atom_id == 0) then
       call sander_bomb('get_qm_and_mm_boundaries (qmmm_adaptive_module)', &
               'qmmm_vsolv%qm_center_atom_id is zero', &
               'Please choose non-zero value for adaptive QM/MM!!!')
    end if


    ! Calculate the distance for the outermost QM solvent

    if (qmmm_vsolv%nearest_qm_solvent > 0) then
       largest_qm_id = qmmm_vsolv%nearest_solvent_pointers(qmmm_vsolv%nearest_qm_solvent)
       largest_qm_dist = qmmm_vsolv%nearest_solvent_distances(qmmm_vsolv%nearest_qm_solvent)
    else
       largest_qm_id = qmmm_vsolv%fixed_iqmatoms(1)
       largest_qm_dist = 0.d0
    end if

    smallest_mm_id = qmmm_vsolv%nearest_solvent_pointers(qmmm_vsolv%nearest_qm_solvent+1)
    smallest_mm_dist = qmmm_vsolv%nearest_solvent_distances(qmmm_vsolv%nearest_qm_solvent+1)

    if (self%debug) then
       write(6,'(a)') '>>>>> leaving get_qm_and_mm_boundaries (qmmm_adaptive_module)'
       call flush(6)
    end if

  end subroutine get_qm_and_mm_boundaries



  ! called by masters of all groups but the main master
  subroutine send_atid_dist (self, id1, distance1, id2, distance2)

    implicit none

#include "parallel.h"
    include 'mpif.h'

    type(qmmm_adaptive_type), intent(in) :: self
    integer, intent(in) :: id1, id2
    _REAL_ , intent(in) :: distance1, distance2

    integer :: ist(MPI_STATUS_SIZE), ierr

    if (self%debug) then
       write(6,'(a)') '>>>>> entered send_atid_dist (qmmm_adaptive_module)'
       call flush(6)
    end if

    ! Send adaptive QMMM data to the main master
    ! which has ID = 0 in the commmaster communicator
    call mpi_send(id1, 1, MPI_INTEGER, &
                  0, 0, commmaster, ist, ierr)
    call mpi_send(id2, 1, MPI_INTEGER, &
                  0, 0, commmaster, ist, ierr)
    call mpi_send(distance1, 1, MPI_DOUBLE_PRECISION, &
                  0, 0, commmaster, ist, ierr)
    call mpi_send(distance2, 1, MPI_DOUBLE_PRECISION, &
                  0, 0, commmaster, ist, ierr)

    if (self%debug) then
       write(6,'(a)') '<<<<< leaving send_atid_dist (qmmm_adaptive_module)'
       call flush(6)
    end if

  end subroutine send_atid_dist


  subroutine receive_atid_dist (self)

    implicit none

#include "parallel.h"
    include 'mpif.h'

    type(qmmm_adaptive_type), intent(inout) :: self

    integer :: i, atom_id1, atom_id2, slave_id
    _REAL_  :: dist1, dist2
    integer :: ist(MPI_STATUS_SIZE), ierr
    
    if (self%debug) then
       write(6,'(a)') '>>>>> entered master_receive_atid_dist (qmmm_adaptive_module)'
       call flush(6)
    end if

    if ( self%debug .or. (self%verbosity > 1) ) then
       write(6,'(a)') 'group id, outermost atom id, outermost distance:'
    end if

    ! Loop over all the slave master threads and receive the
    ! distances and atom ids of the outermost QM atoms
    do i = 2, self%n_partition 
       slave_id = i - 1
       call mpi_recv(atom_id1, 1, MPI_INTEGER, &
                     slave_id, 0, commmaster, ist, ierr)
       call mpi_recv(atom_id2, 1, MPI_INTEGER, &
                     slave_id, 0, commmaster, ist, ierr)
       call mpi_recv(dist1, 1, MPI_DOUBLE_PRECISION, &
                     slave_id, 0, commmaster, ist, ierr)
       call mpi_recv(dist2, 1, MPI_DOUBLE_PRECISION, &
                     slave_id, 0, commmaster, ist, ierr)
       self%ts_solvent_pointers(i)  = atom_id1
       self%ts_solvent_distances(i) = dist1
       if ( self%debug .or. (self%verbosity > 1) ) then
          write(6,'(i5,i6,2f12.6)') slave_id, atom_id1, dist1, dsqrt(dist1)
       end if
       if (i == self%n_partition) then
          self%ts_solvent_pointers(self%n_partition+1) = atom_id2
          self%ts_solvent_distances(self%n_partition+1) = dist2
       end if
    end do


    if (self%debug) then
       write(6,'(a)') '<<<<< leaving master_receive_atid_dist (qmmm_adaptive_module)'
       call flush(6)
    end if

  end subroutine receive_atid_dist


  ! This routine estabishes the outer border (RT) of the T-region.
  ! RT is defined by the 1st solvent molecule that lies beyond the
  ! nearest_qm_solvent+n_partition solvent molecule from the qm_center atom
  subroutine get_outer_ts_border(self)

    use qmmm_module, only: qmmm_vsolv
    implicit none

    type(qmmm_adaptive_type), intent(inout) :: self

    if (self%debug) then
       write(6,'(a)') '>>>>> entered get_outer_ts_border (qmmm_adaptive_module)'
       call flush(6)
    end if
   
    self%ts_solvent_pointers(self%n_partition+1)  =  &
                       qmmm_vsolv%nearest_solvent_pointers(qmmm_vsolv%nearest_qm_solvent+1)

    self%ts_solvent_distances(self%n_partition+1) =  &
                       qmmm_vsolv%nearest_solvent_distances(qmmm_vsolv%nearest_qm_solvent+1)

    if (self%debug) then
       write(6,'(a)') '<<<<< leaving get_outer_ts_border (qmmm_adaptive_module)'
       call flush(6)
    end if

  end subroutine get_outer_ts_border

  
  subroutine check_partition(self)

    use qmmm_module, only : qmmm_nml

    implicit none

    type(qmmm_adaptive_type), intent(inout) :: self

    integer :: i
    _REAL_, pointer :: tsd(:)
 
    tsd => self%ts_solvent_distances
    self%calc_partition(:) = 0

    if (qmmm_nml%vsolv == 2) then

       self%calc_partition(:) = 1

    else if (qmmm_nml%vsolv == 3) then

       if (self%RA < tsd(1)) then
          call sander_bomb('adaptive_qmmm (qmmm_adaptive_module)', &
          'Radius of QM region is smaller than the first partition. ', &
          'Reduce the number of solvent molecules in the active region or increase RA.')
       end if
       if (self%RT > tsd(self%n_partition + 1)) then
          call sander_bomb('adaptive_qmmm (qmmm_adaptive_module)', &
          'Radius of T region is larger than the biggest partition. ', &
          'Increase the number of partitions or reduce RT.')
       end if

       do i = 1, self%n_partition-1
          if (tsd(i+1) >= self%RA .and. tsd(i+1) <= self%RT) then
             self%calc_partition(i) = 1
          endif
       end do

       if ( tsd(self%n_partition) <= self%RT ) then
           self%calc_partition(i) = 1
       end if 

    end if


  end subroutine check_partition



  ! This subroutine computes the weights for all partitionings,
  ! based on the distances (which are actually squared!!).
   subroutine calculate_weights (self,nstep)

    use qmmm_module, only : qmmm_nml

    implicit none

    type(qmmm_adaptive_type), intent(inout) :: self

    integer :: i
    integer, intent(in) :: nstep
    _REAL_, pointer :: tsd(:)
    _REAL_  :: lambda(self%n_partition + 1)

    if (self%debug) then
       write(6,'(a)') '>>>>> entered calculate_weights (qmmm_adaptive_module)'
       call flush(6)
    end if

    tsd => self%ts_solvent_distances

    if (qmmm_nml%vsolv == 2 ) then
       self%RA = tsd(1)
       self%RT = tsd(self%n_partition + 1)
    end if

    if ( self%debug .or. (self%verbosity > 1) ) then
       write(6,'(a,2f10.3)') 'RA and RT : ', self%RA, self%RT
       write(6,'(a,i6)') 'T-region distances, Lambda values at step:',nstep
    end if

    ! Determine lambda values
    do i = 1, self%n_partition + 1

       call calculate_lambda(tsd(i), self%RA, self%RT, lambda(i))

       if ( self%debug .or. (self%verbosity > 1) ) then
          write(6,'(2i5,f12.6,f12.6)') i, self%ts_solvent_pointers(i), &
                                       tsd(i), lambda(i)
       end if

    end do

    if ( self%debug .or. (self%verbosity > 1) ) then
       write(6,'(a,i6)') 'Sigma values at step:', nstep
    end if

    ! Determine sigma values
    do i = 1, self%n_partition

       self%weights(i) = lambda(i+1) - lambda(i)

       if ( self%debug .or. (self%verbosity > 1) ) then
          write(6,'(2i5,2f12.6)') i, self%ts_solvent_pointers(i),self%weights(i)
       end if
    end do

    if (self%debug) then
       write(6,'(a)') '<<<<< leaving calculate_weights (qmmm_adaptive_module)'
       call flush(6)
    end if

  end subroutine calculate_weights


  ! This subroutine computes the MM character (lambda) of
  ! a solvent molecule in the T-region
  subroutine calculate_lambda (dist, RA, RT, lambda)

    use constants, only : two, three
    implicit none

    !Passed in
    _REAL_, intent(in)  :: dist, RA, RT
    _REAL_, intent(out) :: lambda

    if (dist < RA) then
       lambda = 0.d0
    else if ((dist >= RA) .and. (dist <= RT)) then
       lambda = ( (dist - RA)**two ) * ( three*RT - RA - two*dist ) / (RT - RA)**three
    else if (dist > RT) then
       lambda = 1.d0
    end if


  end subroutine calculate_lambda


  ! Send the forces to the main master, to be mixed (only called by group masters)
  subroutine send_forces (self, natom, force)

    implicit none

#include "parallel.h"
    include 'mpif.h'

    type(qmmm_adaptive_type), intent(in) :: self
    integer                 , intent(in) :: natom
    _REAL_                  , intent(in) :: force(3*natom)

    ! Local
    integer :: ist(MPI_STATUS_SIZE), ierr

    if (self%debug) then
       write(6,'(a)') '>>>>> entered send_forces (qmmm_adaptive_module)'
       call flush(6)
    end if

    call mpi_send(force, 3*natom, MPI_DOUBLE_PRECISION, &
                  0, 0, commmaster, ist, ierr)
    if (self%debug) then
       write(6,'(a)') '<<<<< leaving send_forces (qmmm_adaptive_module)'
       call flush(6)
    end if

  end subroutine send_forces

  ! Send the potential energies to the main master, to be mixed (only called by group masters)
  subroutine send_energy (self, energy)

    implicit none

#include "parallel.h"
    include 'mpif.h'

    type(qmmm_adaptive_type), intent(in) :: self
    _REAL_, intent(in) :: energy

    ! Local
    integer :: ist(MPI_STATUS_SIZE), ierr

    if (self%debug) then
       write(6,'(a)') '>>>>> entered send_energy (qmmm_adaptive_module)'
       call flush(6)
    end if

    call mpi_send(energy, 1, MPI_DOUBLE_PRECISION, &
                  0, 0, commmaster, ist, ierr)
    if (self%debug) then
       write(6,'(a)') '<<<<< leaving send_energy (qmmm_adaptive_module)'
       call flush(6)
    end if

  end subroutine send_energy


  ! This subroutine combines all the forces, using the previously
  ! computed weights.
  ! The forces from the master nodes are passed in (mixed_force), and then
  ! mixed with the forces from the other groups, and sent back out.
  subroutine receive_and_mix_forces (self, natom, mixed_force)

    use constants, only: zero

    implicit none

#include "parallel.h"
    include 'mpif.h'

    type(qmmm_adaptive_type), intent(inout) :: self
    integer                 , intent(in) :: natom
    _REAL_, intent(inout) :: mixed_force(3,natom)

    integer :: i_partition
    _REAL_  :: sigma
    _REAL_  :: forces(3,natom)
    integer :: ist(MPI_STATUS_SIZE), ierr

    if (self%debug) then
       write(6,'(a)') '>>>>> entered receive_and_mix_forces (qmmm_adaptive_module)'
       call flush(6)
    end if

    sigma  = self%weights(1)
    mixed_force(:,:) = sigma * mixed_force(:,:)

    do i_partition = 2,self%n_partition

       call mpi_recv(forces, 3*natom, MPI_DOUBLE_PRECISION, &
                     i_partition-1, 0, commmaster, ist, ierr)
       sigma = self%weights(i_partition)

       mixed_force(:,:) = mixed_force(:,:) + (sigma * forces(:,:))

    end do

    if (self%debug) then
       write(6,'(a)') '>>>>> leaving receive_and_mix_forces (qmmm_adaptive_module)'
       call flush(6)
    end if

  end subroutine receive_and_mix_forces


  ! Check which partitionings between the present and the previous time steps match
  !
  ! For self%wbk = 1:
  ! (one-sided differentation)
  ! If a QM/MM partitioning from the previous time step has disappeared,
  ! its corresponding energy needs to be calculated at the present time step
  !
  ! For self%wbk = 2:
  ! (two-sided central-difference)
  ! We can distinguish 6 cases out of which we only need to consider 4
  ! (O means a partitioning that is present, X means it is not present)
  ! (The present time step is t+1, but we need to calculate the correction Wbk
  !  at step t)
  !
  !    (t-1) (t) (t+1)
  ! 1)   O    O    O
  ! 2)   O    O    X
  ! 3)   X    O    O
  ! 4)   O    X    O
  ! 5)   O    X    X  <- equivalent to 3)
  ! 6)   X    X    O  <- equivalent to 2)
  !
  ! 1): nothing needs to be done
  ! 2): We have the energy from step (t)
  !     For 6) we need the energy of the present partitioning (t+1)
  !     with coordinates from step (t)
  ! 3): We have the energy from step (t)
  !     For 5) we need the energy of the previous partitioning (t-1)
  !     with the coordinates from step (t)
  ! 4): We need the energy of the present partitioning (t+1)
  !     with coordinates from step (t)
  ! 
  subroutine check_match_partition (self, nearest_solvent_pointers, nearest_qm_solvent)

    use constants, only: zero, half

    implicit none

    type(qmmm_adaptive_type), intent(inout) :: self
    integer, intent(in) :: nearest_solvent_pointers(:)
    integer, intent(in) :: nearest_qm_solvent

    integer :: i_partition
    integer :: i, j
    integer :: match
    integer :: nqs, nts, nads  ! solvent moleculs full QM, in ts, and all adaptive

    ! local for convenience, may be stored somewhere else later on
    integer, dimension(1:nearest_qm_solvent+self%n_partition-1) :: solv_p_1
    integer, dimension(1:nearest_qm_solvent+self%n_partition-1) :: solv_p_2
    integer, dimension(1:nearest_qm_solvent+self%n_partition-1) :: solv_p_3


    if (self%debug) then
       write(6,'(a)') '>>>>> entered check_match_partition (qmmm_adaptive_module)'
       call flush(6)
    end if

    nqs = nearest_qm_solvent
    nts = self%n_partition -1
    nads= nqs + nts
    solv_p_1(1:nqs)     = nearest_solvent_pointers(1:nqs)
    solv_p_1(nqs+1:nads)= self%ts_solvent_pointers(2:self%n_partition)

    solv_p_2(1:nqs)     = self%nearest_solvent_pointers_prev(1:nqs,1)
    solv_p_2(nqs+1:nads)= self%ts_solvent_pointers_prev(2:self%n_partition,1)

    solv_p_3(1:nqs)     = self%nearest_solvent_pointers_prev(1:nqs,2)
    solv_p_3(nqs+1:nads)= self%ts_solvent_pointers_prev(2:self%n_partition,2)

    if (self%debug) then
       write(6,'(a,i3)') ' nqs  = ', nqs
       write(6,'(a,i3)') ' nts  = ', nts
       write(6,'(a,i3)') ' nads = ', nads
       write(6,'(a,25i3)') ' solv_p_1  = ', solv_p_1(:)
       write(6,'(a,25i3)') ' solv_p_2  = ', solv_p_2(:)
       if (self%calc_wbk == 2) then
          write(6,'(a,25i3)') ' solv_p_3  = ', solv_p_3(:)
       end if
    end if

    self%match_partition(1:self%n_partition) = 0

    ! Loop over all QM/MM partitionings
    do i_partition = 1,self%n_partition

       ! identify which partitionings are identical in the 
       ! present and previous time steps and
       ! make a list of partitionings that have disappeared
       nts = i_partition - 1
       nads = nqs + nts
       match = 0
       do i = 1, nads
          do j = 1, nads
             if (solv_p_1(i) == solv_p_2(j)) then
                match = match + 1
             end if
          end do
       end do
       if (match == nads) then
          ! This partitioning existed also in the previous time step
          self%match_partition(i_partition) = 1
       end if

    end do

    if (self%calc_wbk == 2) then

       ! Loop over all QM/MM partitionings
       do i_partition = 1,self%n_partition
    
          ! identify which partitionings are identical in the 
          ! present and previous time steps and
          ! make a list of partitionings that have disappeared
          nts = i_partition - 1
          nads = nqs + nts
          match = 0
          do i = 1, nads
             do j = 1, nads
                if (solv_p_2(i) == solv_p_3(j)) then
                   match = match + 1
                end if
             end do
          end do
          if (match == nads) then
             if (self%match_partition(i_partition) == 1) then
                ! for the case (t-1) (t) (t+1)
                !                O    O    O
                self%match_partition(i_partition) = 1
             else
                ! for the case (t-1) (t) (t+1)
                !                O    O    X
                self%match_partition(i_partition) = 2
             end if
          else
             if (self%match_partition(i_partition) == 1) then
                ! for the case (t-1) (t) (t+1)
                !                X    O    O
                self%match_partition(i_partition) = 3
             else
                ! for the case (t-1) (t) (t+1)
                !                O    X    O
                self%match_partition(i_partition) = 4
             endif
          end if
       
       end do

    end if

    if (self%debug) then
       write(6,'(a,25i3)') 'match_partition = ', (self%match_partition(i), i=1,self%n_partition)
       write(6,'(a)') '>>>>> leaving check_match_partition (qmmm_adaptive_module)'
       call flush(6)
    end if

  end subroutine check_match_partition


  ! This subroutine receives energies from all other master threads
  subroutine receive_energies (self,iflag)

    implicit none

#include "parallel.h"
    include 'mpif.h'

    type(qmmm_adaptive_type), intent(inout) :: self
    integer, intent(in) :: iflag

    integer :: i_partition
    integer :: ist(MPI_STATUS_SIZE), ierr
    _REAL_ :: energies

    if (self%debug) then
       write(6,'(a)') '>>>>> entered receive_energies (qmmm_adaptive_module)'
       call flush(6)
    end if

    do i_partition = 2,self%n_partition

       call mpi_recv(energies, 1, MPI_DOUBLE_PRECISION, &
                     i_partition-1, 0, commmaster, ist, ierr)

       ! iflag = 0 : send potential energies of each partition for the current step
       ! iflag = 1 : send recalculated potential energies for Wbk term
       if (iflag == 0) then
          self%energies(i_partition) = energies
       else if (iflag == 1) then
          self%energies_recalculated(i_partition) = energies
       end if

    end do

    if (self%debug) then
       write(6,'(a)') '>>>>> leaving receive_energies (qmmm_adaptive_module)'
       call flush(6)
    end if
    
  end subroutine receive_energies


  ! Send the mixed forces to the slaves.
  subroutine send_mixed_forces (self, natom, mixed_forces)

    implicit none

#include "parallel.h"
    include 'mpif.h'

    type(qmmm_adaptive_type), intent(in) :: self
    integer                 , intent(in) :: natom
    _REAL_                  , intent(in) :: mixed_forces(3,natom)

    integer :: i
    integer :: ist(MPI_STATUS_SIZE), ierr

    if (self%debug) then
       write(6,'(a)') '>>>>> entered send_mixed_forces (qmmm_adaptive_module)'
       write(6,'(a)') ' mixed forces to be sent:'
       do i = 1, natom
          write(6,'(i5,3f12.6)') i, mixed_forces(1:3,i)
       end do
       call flush(6)
    end if

    do i = 2, worldsize
       call mpi_send(mixed_forces, 3*natom, MPI_DOUBLE_PRECISION, &
                     i-1, 0, commworld, ist, ierr)
    end do

    if (self%debug) then
       write(6,'(a)') '<<<<< leaving send_mixed_forces (qmmm_adaptive_module)'
       call flush(6)
    end if

  end subroutine send_mixed_forces


  ! Recieve the mixed forces from the master, and accept
  ! them as the new forces
  subroutine receive_mixed_forces (self, natom, mixed_forces)

    implicit none

#include "parallel.h"
    include 'mpif.h'

    type(qmmm_adaptive_type), intent(in) :: self
    integer                 , intent(in)  :: natom
    _REAL_                  , intent(out) :: mixed_forces(3,natom)

    integer :: i
    integer :: ist(MPI_STATUS_SIZE), ierr

    if (self%debug) then
       write(6,'(a)') '>>>>> entered receive_mixed_forces (qmmm_adaptive_module)'
       call flush(6)
    end if

    call mpi_recv(mixed_forces, 3*natom, MPI_DOUBLE_PRECISION, &
                  0, 0, commworld, ist, ierr)

    if (self%debug) then
       write(6,'(a)') ' received following mixed forces:'
       do i = 1, natom
          write(6,'(i5,3f12.6)') i, mixed_forces(1:3,i)
       end do
       write(6,'(a)') '<<<<< leaving receive_mixed_forces (qmmm_adaptive_module)'
       call flush(6)
    end if

  end subroutine receive_mixed_forces


  ! Broadcast qmmm_adaptive data
  ! ATTENTION: GROUP MASTERS RECEIVE ONLY CONTENT OF &adqmmm namelist
  ! Data that changes during the course of the simulations such as
  ! distances of particles in the QM region from the QM center
  ! are not broadcast since only the "main master" needs these
  ! for force mixing
  subroutine broadcast_qmmm_adaptive_type(self)

    implicit none

#include "parallel.h"
    include 'mpif.h'

    type(qmmm_adaptive_type), intent(inout) :: self

    integer :: ierr

    ! write(6,'(a)') '>>>>> entered broadcast_qmmm_adaptive_type (qmmm_adaptive_module)'

    call mpi_bcast(self%debug,        1, mpi_logical, 0, commworld, ierr)
    call mpi_bcast(self%verbosity,    1, mpi_integer, 0, commworld, ierr)
    call mpi_bcast(self%n_partition,  1, mpi_integer, 0, commworld, ierr)
    call mpi_bcast(self%calc_wbk,     1, mpi_integer, 0, commworld, ierr)

    call mpi_bcast(self%print_qm_coords, 1, mpi_logical, 0, commworld, ierr)

    ! write(6,'(a)') '<<<<< leaving broadcast_qmmm_adaptive_type (qmmm_adaptive_module)'

  end subroutine broadcast_qmmm_adaptive_type


  ! print qm coordinates in xyz format
  subroutine print_qm_coords(nstep, masterrank, coords, natom)

    use qmmm_module, only : qmmm_struct
    use molecule, only: mol_info

    implicit none

    integer, intent(in) :: nstep, masterrank
    integer, intent(in) :: natom
    _REAL_, intent(in)  :: coords(3, natom)

    ! Local variable
    _REAL_ :: wrapped_coords(3,natom)
    character(len=2) :: atomtype(106)
    integer :: i, iat, iatyp
    character(len=17) :: filename

#include "box.h"  /* Provides: ifbox  : box type */

    data atomtype / &
         'H ', 'He', & 
         'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
         'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', &
         'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', &
         'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', &
         'Br', 'Kr',                                     &
         'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', &
         'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', &
         'I ', 'Xe',                                     &
         'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', &
         'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
         'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', &
         'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', &
         'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U ', 'Np', 'Pu', &
         'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', &
         'Lr', 'Rf', 'Db', 'Sg' /      

    write(filename, '(a,i3.3,a)') 'QM_coords.', masterrank+1, '.xyz'

    if (nstep == 0) then
       call system('rm -f '//filename)
       return
    end if

    open(77,file=filename,position='append')

    write(77,'(i5)')   qmmm_struct%nquant_nlink
    write(77,'(a,i9)') 'step = ',nstep


    wrapped_coords = coords
    call wrap_molecules(mol_info%nres,mol_info%natom_res,wrapped_coords)
    if(ifbox == 2) call wrap_to(mol_info%nres,mol_info%natom_res,wrapped_coords,box)

    do i=1,qmmm_struct%nquant_nlink

       iat   = qmmm_struct%iqmatoms(i)
       iatyp = qmmm_struct%iqm_atomic_numbers(i)
       if ( iatyp < 107 ) then
          write(77,'(a,3f15.3)') atomtype(iatyp), wrapped_coords(1:3,iat)
       else
          write(77,'(a)') 'unknown atom type !! '
          write(77,'(a)') 'modify print_qm_coords (qmmm_adaptive_module)'
       end if

    end do

    close(77)

  end subroutine print_qm_coords


  subroutine print_adqmmm_info(self, nstep, ntpr)

    implicit none

    type(qmmm_adaptive_type), intent(in) :: self
    integer, intent(in) :: nstep
    integer, intent(in) :: ntpr

    integer :: i, idum1
    logical :: switch_partition

    if (self%debug) then
       write(6,'(a)') '>>>>> entered print_adqmmm_info (qmmm_adaptive_module)'
       call flush(6)
    end if

    if ( nstep == 0 ) then
       call system('rm -f adqmmm_*.dat')
    end if

    switch_partition = .false.

    if ( self%calc_wbk > 0 ) then
       do i=1, self%n_partition
          if ( nstep > self%calc_wbk ) then
             if ( self%match_partition(i) /= 1 ) then
                switch_partition = .true.
             end if
          end if
       end do
    end if

    if ( switch_partition ) then
       open(79, file='adqmmm_switch.dat', position='append')
       write(79,'(i9,5x,20i2)')   nstep, self%match_partition(1:self%n_partition)
       close(79)
    end if

    if (self%calc_wbk == 2) then
       idum1 = nstep + 1
    else 
       idum1 = nstep
    end if

    if ( mod(idum1, ntpr) == 0 ) then
       open(80, file='adqmmm_weights.dat', position='append')
       write(80,'(i9,5x,20f8.3)')nstep, self%weights(1:self%n_partition)
       close(80)
    end if

    if ( mod(idum1, ntpr) == 0 ) then
       open(81, file='adqmmm_res_distances.dat', position='append')
       write(81,'(i9,5x,20f8.3)')nstep, &
                       self%ts_solvent_distances(1:self%n_partition+1)
       close(81)
    end if

    if (self%debug) then
       write(6,'(a)') '>>>>> entered print_adqmmm_info (qmmm_adaptive_module)'
       call flush(6)
    end if

  end subroutine print_adqmmm_info



#else

  ! DUMMY ROUTINE FOR SERIAL VERSION
  subroutine adaptive_qmmm (natom, unimaged_coords, force)

      implicit none

      integer :: natom
      _REAL_ :: unimaged_coords(3,natom)
      _REAL_ :: force(3,natom)

      call sander_bomb('adaptive_qmmm (qmmm_adaptive_module)', &
           'Adaptive QM/MM needs parallel version (multisander)', &
           'Will quit now')

   end subroutine adaptive_qmmm

#endif

   subroutine adaptive_reset

      implicit none

      first_call = .true.

   end subroutine adaptive_reset

end module qmmm_adaptive_module
