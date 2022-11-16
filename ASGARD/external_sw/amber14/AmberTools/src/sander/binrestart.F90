! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"

!! Module for generating binary restart-type output in NetCDF format
!! Developed by Dan Roe
!! 2010-01-10

module binrestart
   private

   integer, save :: coordVID, velocityVID, cellAngleVID, cellLengthVID
   integer, save :: timeVID, TempVID
   integer, save :: remd_indices_var_id, remd_groups_var_id

   public write_nc_restart, &
          read_nc_restart_box, &
          read_nc_restart, &
          read_nc_restart_extents
contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write Netcdf restart file.
!-------------------------------------------------------------------
!     --- WRITE_NC_RESTART ---
!-------------------------------------------------------------------
!     Write Netcdf Restart file with given filename and title. 
!     owrite indicates overwrite status (N is no overwrite), natom 
!     is the # atoms, ntb>0 indicates presence of box coords, first
!     indicates the file should be created and set-up, Coords and 
!     Velo are the coordinates and velocities, temp0 is the current
!     temperature and Time is the current simulation time. If hasV
!     is false, no velocities will be written (for e.g. during min)
subroutine write_nc_restart(filename,title,owrite,natom,ntb,first,Coords,Velo,&
                            Time,hasVin&
#ifdef MPI
                            , temp0, rem, remd_dimension, remd_types, group_num &
                            , replica_indexes, stagid &
#endif
                           )
   use AmberNetcdf_mod
#ifdef BINTRAJ
   use netcdf
   use nblist, only: a,b,c,alpha,beta,gamma
#endif
#ifdef MPI
   use sgld, only: trxsgld
#endif

   implicit none
   ! Input variables
   character(len=*), intent(in)      :: filename
   character(len=*), intent(in)      :: title
   character, intent(in)             :: owrite
   integer, intent(in)               :: natom, ntb
   logical, intent(in)               :: first
   _REAL_, dimension(*), intent(in)  :: Coords, Velo
   _REAL_, intent(in)                :: Time
   logical, intent(in)               :: hasVin
#  ifdef MPI
   _REAL_, intent(in)                :: temp0
   integer, intent(in)               :: rem, remd_dimension
   integer, intent(in), dimension(:) :: remd_types, group_num, replica_indexes
   integer, intent(in)               :: stagid
#  endif
#ifdef BINTRAJ
   ! Local vars
   integer :: ncid, natom3
   logical :: hasTemperature = .false.
   integer :: frcVID ! dummy variable
#  ifdef MPI
   hasTemperature = (rem.ne.0 .and. rem.ne.3)
#  endif
   if (first) then
      ! If first call, create the file and set up all dimensions and vars
      ! owrite status code: 'N', 'O', 'R', 'U' = new, old, replace, unknown
      ! sander flag -O='R', -A='U', default='N'
      if (NC_create(filename, owrite, .true., natom, .true.,&
                    hasVin, (ntb.gt.0), hasTemperature, .true., .false., &
                    title, ncid, timeVID, coordVID, velocityVID, frcVID, &
                    cellLengthVID, cellAngleVID, TempVID)) call mexit(6,1)
#     ifdef MPI
      ! MREMD indices
      if (rem .eq. -1) then
         ! For now, use indexes only for multid remd
         if (NC_defineRemdIndices(ncid, remd_dimension, remd_indices_var_id,&
                                 remd_types, .true.,&
                                 remd_groupsVID=remd_groups_var_id ))&
           call mexit(6,1)
      endif
#     endif
   else
      ! If not the first call, just reopen the existing file
      if (NC_openWrite(filename, ncid)) then
        write (6,'(a)') 'write_nc_restart(): Could not open restart'
        call mexit(6,1)
      endif
   endif

   natom3 = natom * 3
   ! Write time
   call checkNCerror(nf90_put_var(ncid, timeVID, Time), 'write time')
   ! Write coords
   call checkNCerror(nf90_put_var(ncid,coordVID, Coords(1:natom3), &
                  start = (/ 1, 1 /), count = (/ 3, natom /) ), 'write atom coords')
   ! Write velocities TODO: Should this check velocityVID instead?
   if (hasVin) then
      call checkNCerror(nf90_put_var(ncid,velocityVID, Velo(1:natom3), &
                     start = (/ 1, 1 /), count = (/ 3, natom /) ), 'write velocities')
   endif
   ! Write box information
   if (ntb > 0) then
      call checkNCerror(nf90_put_var(ncid,cellLengthVID, &
              (/ a, b, c /), start = (/ 1 /), count = (/ 3 /) ), 'write cell lengths')
      call checkNCerror(nf90_put_var(ncid,cellAngleVID, &
              (/ alpha,beta,gamma /), start = (/ 1 /), count = (/ 3 /) ), &
           'write cell angles')
   endif
#  ifdef MPI
   ! Write replica temperature, indices
   if (rem.ne.0) then
      ! multi-D remd: Store indices of this replica in each dimension
      if (rem .eq. -1) then
         call checkNCerror(nf90_put_var(ncid, remd_indices_var_id, &
                           replica_indexes(:), &
                           start = (/ 1 /), count = (/ remd_dimension /)), &
                           'write replica index for each dimension')
         call checkNCerror(nf90_put_var(ncid, remd_groups_var_id, group_num(:), &
                           start = (/ 1 /), count = (/ remd_dimension /)), &
                           'write replica group for each dimension')
      endif
      if (trxsgld) then
         call checkNCerror(nf90_put_var(ncid, TempVID, REAL(stagid)), &
                           'write SGLD replica index')
      else
         call checkNCerror(nf90_put_var(ncid, TempVID, temp0), 'write temp0')
      endif
   endif
#  endif
   ! Close restart file       
   call NC_close(ncid)
#else
   call NC_NoNetcdfError(6) 
   call mexit(6,1)
#endif
end subroutine write_nc_restart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read box information from a Netcdf restart file.
!-------------------------------------------------------------------
!     --- READ_NC_RESTART_BOX ---
!-------------------------------------------------------------------
!     Read box information from the Netcdf Restart file with 
!     specified filename.
!     The box read is called from load_ewald_info() in ew_setup.f
!     and is separate from the coord/velocity read since the box
!     information is needed to set up certain ewald parameters.
subroutine read_nc_restart_box(filename,a,b,c,alpha,beta,gamma)
   use AmberNetcdf_mod 
   implicit none

   character(len=*), intent(in) :: filename
   _REAL_, intent(out) :: a,b,c,alpha,beta,gamma
   ! local
   integer ncid
#ifdef BINTRAJ
   if (NC_openRead(filename, ncid)) call mexit(6,1)
   if (NC_readRestartBox(ncid,a,b,c,alpha,beta,gamma)) call mexit(6,1)
   call NC_close(ncid)
   write(6,'(a)') '| NetCDF restart box info found'
#else
   call NC_NoNetcdfError(6)
   call mexit(6,1)
#endif
end subroutine read_nc_restart_box

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read coord/velocity information from a Netcdf restart file.
!-------------------------------------------------------------------
!     --- READ_NC_RESTART ---
!-------------------------------------------------------------------
!     Read coordinates and velocities from the Netcdf Restart file 
!     with specified filename. This is called from getcor.f. Title
!     will be read in and set. 
!     ntx specifies whether coords and velo or just coords will be 
!     read. ntx=1 means read coords only, and ntx=5 means read 
!     coords and velocities.
!     parmatoms is the expected number of atoms in the restart. If -1, we assume
!        it is not set (if the prmtop hasn't been read yet, such as for the API).
!        In this case, parmatoms is set to the number of atoms in the inpcrd
!        file before returning.
!     Coords and Velo are the coordinates and 
!     velocities, temp0 is the temperature (if present) and Time
!     is the time.
!     NOTE: Box info is not read here; it is obtained using 
!     read_nc_restart_box in load_ewald_info.
subroutine read_nc_restart(filename,title,ntx,parmatoms,Coords,Velo,temp0,Time)
   use AmberNetcdf_mod
#ifdef BINTRAJ
   use netcdf
#endif
   use constants, only : NO_INPUT_VALUE_FLOAT
   implicit none

   character(len=*), intent(in)      :: filename
   character(len=80), intent(out)    :: title
   integer, intent(in)               :: ntx
   integer, intent(in out)           :: parmatoms
   _REAL_, dimension(*), intent(out) :: Coords, Velo
   _REAL_, intent(out)               :: temp0, Time
#ifdef BINTRAJ
   integer :: ncid, ncatom, ncatom3

   ! ---=== Open file
   if (NC_openRead(filename, ncid)) call mexit(6,1)
   ! Setup restart: title, coordVID, velocityVID, time, TempVID
   if (NC_setupRestart(ncid, title, ncatom, coordVID, velocityVID, &
                       TempVID, Time)) call mexit(6,1)
   ! Check that number of atoms matches
   if (ncatom /= parmatoms .and. parmatoms /= -1) then
      write(6,'(2x,a)') "FATAL: NATOM mismatch in restart and topology files."
      call mexit(6,1)
   else if (parmatoms == -1) then
      parmatoms = ncatom
   endif
   ncatom3 = ncatom * 3
   ! ---=== Get Coords
   if (NC_error(nf90_get_var(ncid, coordVID, Coords(1:ncatom3), &
                             start = (/ 1, 1 /), count = (/ 3, ncatom /)),&
                'reading restart coordinates')) call mexit(6,1)
   ! ---=== Get velocities
   !        ntx=1 No Velocity Read 
   !        ntx=5 Read Velocity
   if (ntx == 5) then
      if (velocityVID .eq. -1) then
        write(6,'(2x,a)') "FATAL: ntx=5 specified but no velocities in INPCRD"
        call mexit(6,1)
      endif
      if (NC_error(nf90_get_var(ncid, velocityVID, Velo(1:ncatom3), &
                                start = (/ 1, 1 /), count = (/ 3, ncatom /)),&
                   'reading restart velocities')) call mexit(6,1)
   endif
   ! ---=== Replica Temperature
   if (TempVID .ne. -1) then
      if (NC_error(nf90_get_var(ncid, TempVID, temp0),&
                   "read_nc_restart(): Getting restart temperature")) call mexit(6,1)
   else
     temp0=NO_INPUT_VALUE_FLOAT
   endif 

   ! NOTE: TO BE ADDED
   !labelDID;
   !int cell_spatialDID, cell_angularDID;
   !int spatialVID, cell_spatialVID, cell_angularVID;
  
   ! ---=== Close file
   call NC_close(ncid)
#else
   call NC_NoNetcdfError(6)
   call mexit(6,1)
#endif
end subroutine read_nc_restart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read coord/velocity information from a Netcdf restart file.
!-------------------------------------------------------------------
!     --- READ_NC_RESTART ---
!-------------------------------------------------------------------
!     Read coordinates and velocities from the Netcdf Restart file 
!     with specified filename. This is called from getcor.f. Title
!     will be read in and set. 
!     ntx specifies whether coords and velo or just coords will be 
!     read. ntx=1 means read coords only, and ntx=5 means read 
!     coords and velocities.
!     parmatoms is the expected number of atoms in the restart.
!     Coords and Velo are the coordinates and 
!     velocities, temp0 is the temperature (if present) and Time
!     is the time.
!     NOTE: Box info is not read here; it is obtained using 
!     read_nc_restart_box in load_ewald_info.
subroutine read_nc_restart_extents(filename, extents)
   use AmberNetcdf_mod
#ifdef BINTRAJ
   use netcdf
#endif
   use constants, only : NO_INPUT_VALUE_FLOAT
   implicit none

   character(len=*), intent(in)        :: filename
   _REAL_, dimension(3,2), intent(out) :: extents
#ifdef BINTRAJ
   _REAL_, dimension(:), allocatable :: Coords
   _REAL_  :: Time
   integer :: ncid, ncatom, ncatom3, ierror, j, i
   character(len=256) :: title

   title = ''
   ! ---=== Open file
   if (NC_openRead(filename, ncid)) call mexit(6,1)
   ! Setup restart: title, coordVID, velocityVID, time, TempVID
   if (NC_setupRestart(ncid, title, ncatom, coordVID, velocityVID, &
                       TempVID, Time)) call mexit(6,1)
   ! Check that number of atoms matches
   ncatom3 = ncatom * 3
   ! Allocate the coordinate array
   allocate(Coords(ncatom3), stat=ierror)
   REQUIRE( ierror == 0 )
   ! ---=== Get Coords
   if (NC_error(nf90_get_var(ncid, coordVID, Coords(1:ncatom3), &
                             start = (/ 1, 1 /), count = (/ 3, ncatom /)),&
                'reading restart coordinates')) call mexit(6,1)
   ! ---=== Close file
   call NC_close(ncid)
   ! Now go through and find the max/min of x, y, and z
   extents(1, 1) = Coords(1) ! Min X
   extents(2, 1) = Coords(2) ! Min Y
   extents(3, 1) = Coords(3) ! Min Z
   extents(1, 2) = Coords(1) ! Max X
   extents(2, 2) = Coords(2) ! Max Y
   extents(3, 2) = Coords(3) ! Max Z
   j = 4
   do i = 2, ncatom
      if ( Coords(j) < extents(1, 1) ) extents(1, 1) = Coords(j)
      if ( Coords(j) > extents(1, 2) ) extents(1, 2) = Coords(j)
      j = j + 1
      if ( Coords(j) < extents(2, 1) ) extents(2, 1) = Coords(j)
      if ( Coords(j) > extents(2, 2) ) extents(2, 2) = Coords(j)
      j = j + 1
      if ( Coords(j) < extents(3, 1) ) extents(3, 1) = Coords(j)
      if ( Coords(j) > extents(3, 2) ) extents(3, 2) = Coords(j)
      j = j + 1
   end do
   ! Deallocate our work array
   deallocate(Coords, stat=ierror)
   REQUIRE( ierror == 0 )
#else
   call NC_NoNetcdfError(6)
   call mexit(6,1)
#endif
end subroutine read_nc_restart_extents

end module binrestart
