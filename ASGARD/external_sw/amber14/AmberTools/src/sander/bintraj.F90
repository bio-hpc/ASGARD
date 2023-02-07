! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"

!! Module for generating binary trajectory-type output in NetCDF format
!! Developed by John Mongan <jmongan@mccammon.ucsd.edu>, November 2005

module bintraj
   private
   integer :: crd_ncid, vel_ncid, frc_ncid
   integer :: FrameDimID ! Needed for setup_remd_indices
   integer :: CoordVarID, crd_TimeVarID, crd_frame
   integer :: TempVarID, remd_indices_var_id, Cell_lengthVarID, Cell_angleVarID
   integer :: VelocVarID, vel_TimeVarID, vel_frame
   integer :: FrcVarID, frc_TimeVarID, frc_frame

   public open_binary_files, close_binary_files, write_binary_traj, &
          end_binary_frame, setup_remd_indices, check_atom_mismatch

contains

!> MODULE BINTRAJ FUNCTION CHECK_ATOM_MISMATCH
!> @brief Check if # atoms in parm matches # atoms in netcdf
subroutine check_atom_mismatch(parmAtoms, ncatom)
  implicit none
  integer, intent(in) :: parmAtoms, ncatom
  if (ncatom .ne. parmAtoms) then
    write(6,'(2(a,i8))') 'Error: NetCDF traj file has ', ncatom, ' atoms, but&
                         & current topology has ', parmAtoms
    call mexit(6,1)
  endif
end subroutine check_atom_mismatch 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Open the coordinate, velocity file(s)
subroutine open_binary_files
   use AmberNetcdf_mod
#ifdef BINTRAJ
   use file_io_dat
#  ifdef MPI
      use remd, only: rem
#  endif 

   implicit none

#  include "../include/md.h"
#  include "../include/memory.h"
#  include "box.h"
   ! Local vars
   character(80) :: vel_title
   logical       :: crd_file_exists = .false.
   logical       :: vel_file_exists = .false.
   logical       :: frc_file_exists = .false.
   logical       :: hasTemperature = .false.
   integer       :: ierr, atomCnt, ncatom

   crd_ncid = -1
   vel_ncid = -1
   frc_ncid = -1  !AWG added 1/17/14
   if (ntwv == 0 .and. ntwx == 0 .and. ntwf == 0) return

   atomCnt = natom
   if (ntwprt > 0) atomCnt = ntwprt
#  ifdef MPI
   ! Is this a run requiring replica temperature?
   hasTemperature = (rem.ne.0 .and. rem/=3)
#  endif
   ! If append specified, determine if files exist.
   if (facc == 'A') then
     crd_file_exists = NetcdfFileExists( mdcrd )
     vel_file_exists = NetcdfFileExists( mdvel )
     frc_file_exists = NetcdfFileExists( mdfrc )
   endif
   ! Set up output trajectory
   if (ntwx.gt.0) then
      ! If write or append but file does not exist, create trajectory.
      if (facc == 'W' .or. .not.crd_file_exists) then
         ! Create NetCDF trajectory file
         if ( NC_create( mdcrd, owrite, .false., atomCnt, &
                         .true., ntwv.lt.0, ntb.gt.0, hasTemperature, .true., &
                         ntwf.lt.0, title, crd_ncid, crd_TimeVarID, CoordVarID, &
                         VelocVarID, FrcVarID, Cell_lengthVarID, Cell_angleVarID, &
                         TempVarID, frameDID=FrameDimID ) ) then
            write(6,'(a)') "Creation of NetCDF trajectory file failed."
            call mexit(6,1)
         endif
         crd_frame = 1
      else
         ! Not write (i.e. append) and file exists. Set up trajectory.
         ! NOTE: Does this need to be openRead?
         if (NC_openWrite( mdcrd, crd_ncid )) call mexit(6,1)
         if (NC_setupMdcrd(crd_ncid, title, crd_frame, ncatom, &
                           CoordVarID, VelocVarID, crd_TimeVarID, &
                           Cell_lengthVarID, Cell_angleVarID, TempVarID)) call mexit(6,1)
         ! Check for natoms mismatch
         call check_atom_mismatch(atomCnt, ncatom)
         ! Check for box mismatch
         if ( (ntb.gt.0).neqv.(Cell_lengthVarID.ne.-1) ) then
            write(6,'(a)') 'Error: Cannot append to netcdf traj, box info mismatch.'
            call mexit(6,1)
         endif
         ! Check for temperature mismatch 
         if ( (TempVarID.ne.-1).neqv.hasTemperature ) then
            write(6,'(a)') 'Error: Cannot append to netcdf traj, T info mismatch.'
            call mexit(6,1)
         endif
         ! NOTE: Do not set up replica indices yet since REMD setup
         !       has not yet occured.
      endif
   endif      
   ! Set up output velocities
   if (ntwv.gt.0) then
      if (facc == 'W' .or. .not.vel_file_exists) then
         ! Create velocity traj - ierr is a dummy var
         if (NC_create( mdvel, owrite, .false., atomCnt, &
                        .false., .true., .false., .false., .true., &
                        .false., title, vel_ncid, vel_TimeVarID, ierr, &
                        VelocVarID, ierr, ierr, ierr, ierr, ierr ) ) then
            write(6, '(a)') "Creation of NetCDF velocity file failed."
            call mexit(6,1)
         endif
         vel_frame = 1
       else
         ! Not write (i.e. append) and file exists. Set up velocities
         if (NC_openWrite( mdvel, vel_ncid )) call mexit(6,1)
         if (NC_setupMdcrd(vel_ncid, vel_title, vel_frame, ncatom, &
                           ierr, VelocVarID, vel_TimeVarID, &
                           ierr, ierr, ierr)) call mexit(6,1)
         ! Check for natoms mismatch
         call check_atom_mismatch(atomCnt, ncatom)
         ! Check that there is velocity info
         if (VelocVarID .eq. -1) then
            write(6, '(a)') "Error: Cannot append netcdf velocities, no velocity info present."
            call mexit(6,1)
         endif
      endif 
   else if (ntwv.lt.0) then
      if (ntwx.eq.0 .or. crd_ncid.eq.-1) then
         write(6, '(a)') "ntwx = 0 and ntwv < 0 not allowed."
         call mexit(6,1)
      endif
      ! Combined mdcrd/mdvel
      vel_ncid = crd_ncid
      vel_frame = crd_frame
   endif

   ! Set up output forces
   if (ntwf > 0) then
      if (facc == 'W' .or. .not. frc_file_exists) then
         ! Create force trajectory -- ierr is a dummy var
         if (NC_create( mdfrc, owrite, .false., atomCnt, &
                        .false., .false., .false., .false., .true., &
                        .true., title, frc_ncid, frc_TimeVarID, ierr, &
                        ierr, FrcVarID, ierr, ierr, ierr, ierr ) ) then
            write(6, '(a)') "Creation of NetCDF force dump file failed."
            call mexit(6, 1)
         end if
         frc_frame = 1
      else
         write(6, *) 'Error: Cannot append netcdf forces.'
         call mexit(6, 1)
      end if
   else if (ntwf < 0) then
      ! Combined mdcrd/mdfrc
      frc_ncid = crd_ncid
      frc_frame = crd_frame
   end if
#else
   call NC_NoNetcdfError(6) 
   call mexit(6,1)
#endif
end subroutine open_binary_files
!----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Set up dimension information for multi-D REMD. This is in a subroutine 
!+ separate from open_binary_files since REMD setup occurs after the call
!+ to open_binary_files.
subroutine setup_remd_indices
  use AmberNetcdf_mod
#ifdef MPI
   use remd,        only: rem, remd_dimension, remd_types
   use file_io_dat, only: ioutfm
#endif
   implicit none
#ifdef BINTRAJ
#  ifdef MPI
   ! Only setup remd indices if ioutfm and multid remd
   if ( (ioutfm.ne.0) .and. (rem.eq.-1) ) then
      if ( NC_defineRemdIndices(crd_ncid, remd_dimension, remd_indices_var_id, &
                               remd_types, .false., frameDID=FrameDimID) ) call mexit(6,1)
   endif
#  endif
#else
   call NC_NoNetcdfError(6)
   call mexit(6,1)
#endif
   return
end subroutine setup_remd_indices

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Close the binary coordinate, velocity file(s).
subroutine close_binary_files
   use AmberNetcdf_mod
#ifdef BINTRAJ
   use file_io_dat

   implicit none

   if (ntwx > 0) call NC_close(crd_ncid)
   if (ntwv > 0) call NC_close(vel_ncid)
   if (ntwf > 0) call NC_close(frc_ncid)

#else
   call NC_NoNetcdfError(6) 
   call mexit(6,1)
#endif

end subroutine close_binary_files
!----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Emit coordinates or velocities, r(istart:n), to  netcdf.
#ifdef BINTRAJ 
subroutine write_binary_traj(r,istart,n,unit) ! FIXME: Split into crd,vel,box
#ifdef MPI
   use constantph, only : target_ph
#endif
   use file_io_dat
   use netcdf
   use AmberNetcdf_mod
   use nblist, only: a,b,c,alpha,beta,gamma
#  ifdef MPI
      use remd, only: rem, my_remd_data, replica_indexes, remd_dimension, stagid
      use sgld, only: trxsgld
#  endif 
   
   implicit none
#  include "box.h"
   integer, intent(in) :: istart,n,unit
   _REAL_, intent(in) ::  r(n)

   select case (unit)
      case (MDCRD_UNIT)
         if (n == 3 .and. ntb > 0) then       ! Assume this is box (fails on one atom systems)
            call checkNCerror(nf90_put_var(crd_ncid, Cell_lengthVarID, &
                                      (/ a,b,c/), start = (/ 1, crd_frame /), &
                                      count = (/ 3, 1 /)), &
                              'write cell lengths')
            call checkNCerror(nf90_put_var(crd_ncid, Cell_angleVarID, &
                                      (/ alpha,beta,gamma /), &
                                      start = (/ 1, crd_frame /), &
                                      count = (/ 3, 1 /)), &
                              'write cell angles')
         else
            call checkNCerror(nf90_put_var(crd_ncid,CoordVarID, r(istart:n), &
                                      start = (/ 1, 1, crd_frame /), &
                                      count = (/ 3, (n-istart+1)/3, 1 /)), &
                              'write atom coords')
#  ifdef MPI
            ! If this is a replica run write temp0
            if (rem.ne.0) then
               ! multi-D remd: Store indices of this replica in each dimension
               if (rem .eq. -1) then
                  call checkNCerror(nf90_put_var(crd_ncid, remd_indices_var_id, &
                                    replica_indexes(:), &
                                    start = (/ 1, crd_frame /), &
                                    count = (/ remd_dimension, 1 /)), &
                                    'write replica index for each dimension')
               endif
               if (rem==4) then
                  call checkNCerror(nf90_put_var(crd_ncid,TempVarID,target_ph, &
                                                 start = (/ crd_frame /) ), &
                                    'write replica pH')
               else if (trxsgld) then
                  call checkNCerror(nf90_put_var(crd_ncid,TempVarID, &
                                                 REAL(stagid), &
                                                 start = (/ crd_frame /) ), &
                                    'write SGLD replica index')
               else
                  call checkNCerror(nf90_put_var(crd_ncid,TempVarID, &
                                                 my_remd_data%mytargettemp, &
                                                 start = (/ crd_frame /) ), &
                                    'write replica mytargettemp')
               endif
            endif
#  endif
         end if
      case (MDVEL_UNIT)
         call checkNCerror(nf90_put_var(vel_ncid,VelocVarID, r(istart:n), &
                                        start = (/ 1, 1, vel_frame /), &
                                        count = (/ 3, (n-istart+1)/3, 1 /)), &
                           'write velocities')
      case (MDFRC_UNIT)
         call checkNCerror(nf90_put_var(frc_ncid,FrcVarID, r(istart:n), &
                                        start = (/ 1, 1, frc_frame /), &
                                        count = (/ 3, (n-istart+1)/3, 1 /)), &
                           'write forces')
      case default
         write (6,*) 'Error: unhandled unit ',unit,' selected for output in bintraj'
   end select

end subroutine write_binary_traj

#else
subroutine write_binary_traj(r,istart,n,unit)
   use AmberNetcdf_mod
   integer, intent(in) :: istart,n,unit
   _REAL_, intent(in) ::  r(n)
   call NC_NoNetcdfError(6)
   call mexit(6,1)
end subroutine write_binary_traj
#endif

!-----------------------------------------------------------------------
   
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write scalar data and increment frame counter
subroutine end_binary_frame(unit)
   use AmberNetcdf_mod
#ifdef BINTRAJ
   use file_io_dat
   use netcdf
   implicit none
#  include "../include/md.h"
   integer, intent(in) :: unit
 select case (unit)
      case (MDCRD_UNIT)
         call checkNCerror(nf90_put_var(crd_ncid, crd_TimeVarID, &
               (/ t /), start = (/ crd_frame /), count = (/ 1 /)), 'write time')
         
         call checkNCerror(nf90_sync(crd_ncid))

         crd_frame = crd_frame + 1
         ! Sync frame counters if combined trajectory
         if (ntwv.lt.0) vel_frame = vel_frame + 1
         if (ntwf.lt.0) frc_frame = frc_frame + 1
      case (MDVEL_UNIT)
         call checkNCerror(nf90_put_var(vel_ncid, vel_TimeVarID, &
               (/ t /), start = (/ vel_frame /), count = (/ 1 /)), 'write time')
         
         call checkNCerror(nf90_sync(vel_ncid))
         
         vel_frame = vel_frame + 1
      case (MDFRC_UNIT)
         call checkNCerror(nf90_put_var(frc_ncid, frc_TimeVarID, &
               (/ t /), start = (/ frc_frame /), count = (/ 1 /)), 'write time')
         call checkNCerror(nf90_sync(frc_ncid))

         frc_frame = frc_frame + 1
      case default
         write (6,*) 'Error: unhandled unit ',unit,' selected for end frame in bintraj'
   end select

#else
   integer, intent(in) :: unit
   call NC_NoNetcdfError(6) 
   call mexit(6,1)
#endif

end subroutine end_binary_frame

end module bintraj
