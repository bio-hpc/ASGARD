! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"

module trajenemod

! MODULE: TRAJENEMOD
! ================== TRAJECTORY ENERGY POST-PROCESSING ===================
! Daniel R. Roe, 2009
! Based on original implementation by Carlos Simmerling

implicit none

contains

!*********************************************************************
!               SUBROUTINE TRAJENE
!*********************************************************************
! carlos add trajene routine for processing trajectory energies

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine trajene here]
subroutine trajene(x,ix,ih,ipairs,ene,ok,qsetup)
                                                       !\/-----\/--For Debug
   use nblist,only: fill_tranvec,a,b,c,alpha,beta,gamma!,nbflag,cutoffnb
   use constants,only: one,half
   use state
   use file_io_dat, only : INPTRAJ_UNIT, inptraj, ioutfm, MDCRD_UNIT, &
                           ntwx, title
   use AmberNetcdf_mod
   use memory_module, only : lcrd, natom, lforce, lvel, iibh, ijbh, l50, &
                             lwinv, ibellygp, l95
#ifdef BINTRAJ
   use netcdf
   use bintraj,only: check_atom_mismatch, end_binary_frame
#endif

   implicit none
   ! INPUT VARIABLES
   integer ipairs(*),ix(*)
   _REAL_ x(*)
   type(state_rec) :: ene
   character(len=4) ih(*)
   logical, intent(out) :: ok
   logical, intent(in)  :: qsetup

   ! INTERNAL VARIABLES
   integer member,j,xstop
   integer newnfft1,newnfft2,newnfft3
   _REAL_ carrms,oldbox(3)
   logical loutfm, is_netcdf

#ifdef BINTRAJ
   ! For netcdf files
   integer ncid,ncframe,ncatom,coordVID,velocityVID,cellLengthVID,cellAngleVID
   integer err
#endif

#  include "tgtmd.h"
! DAN ROE: Is extra.h still needed?
#  include "extra.h"
! needed for ntb, box
#  include "box.h"
! Needed for debug info 
!#  include "ew_cntrl.h"
!#  include "ew_erfc_spline.h"
! needed for nfft
#  include "ew_pme_recip.h"
!------------------------------------------------
   loutfm = ioutfm <= 0
   is_netcdf=.false.

   ! If PBC, save original nfft sizes. PME allocates memory for nfft based 
   ! on the original box size, so if the box becomes much larger than the 
   ! original box size the new nfft values would cause a memory error.
   if (ntb>0) then
      write (6,'(a,3i6)') "TRAJENE: Original NFFTs: ",nfft1,nfft2,nfft3
   endif

   ! This is for indexing the x array to get coordinates
   xstop=lcrd+(natom*3)-1;

   ! Open the inptraj file for coordinate reading
   if (NC_checkTraj(inptraj)) then
#ifdef BINTRAJ
      ! Open NETCDF format trajectory file
      if (NC_openRead(inptraj, ncid)) then
         write(6,'(a,a)') "TRAJENE: Could not open netcdf file: ",inptraj
         return
      endif
      ! Get netcdf file information. err is dummy var for time and temp
      if (NC_setupMdcrd( ncid, title, ncframe, ncatom, coordVID, velocityVID,&
                         err, cellLengthVID, cellAngleVID, err)) then
        write(6,'(a)') "TRAJENE: Could not set up netcdf file for reading."
        return
      endif
      write (6,'(a,i5)') "TRAJENE: Frames in trajectory= ",ncframe
      write (6,'(a,i5)') "TRAJENE: Atoms in trajectory= ",ncatom
      ! Check for atom mismatch
      call check_atom_mismatch( natom, ncatom )
      ! Check for box mismatch
      if ( (ntb.gt.0).and.(cellLengthVID.eq.-1) ) then
         write(6,'(a)') 'Error: ntb > 0 and no box coordinates in NetCDF traj.'
         call mexit(6,1)
      endif
      is_netcdf=.true.
#else
      ! No NETCDF support
      call NC_NoNetcdfError(6)
      return
#endif
   else 
      ! Standard Amber Trajectory Format
      call amopen(INPTRAJ_UNIT,inptraj,'O','F','R')
      read(INPTRAJ_UNIT,'(a80)') title
   endif
   write (6,'(a80)') title

   ! Now loop over trajectory file, exiting only on error or end of file
   member=1
   do while ( .true. )

      !       --- read next coordinate set from trajectory
      if (.not.is_netcdf) then
         read(INPTRAJ_UNIT,'(10f8.3)',end=1000,err=1010) (x(j),j=lcrd,xstop)
#ifdef BINTRAJ
      else
         if (member>ncframe) goto 1000 
           
         if (NC_error(nf90_get_var(ncid,coordVID,x(lcrd:xstop), &
                                   start = (/ 1, 1, member /), &
                                   count = (/ 3, (xstop-lcrd+1)/3, 1 /)),&
                      'reading netcdf coordinates')) goto 1010
#endif
      endif

      ! DAN ROE:
      ! If box coords are present in trajectory (ntb>0), read them in and
      !  update the unit cell and grid size.
      if (ntb>0) then
         ! 1- Save old box coords.
         oldbox(1) = box(1)
         oldbox(2) = box(2)
         oldbox(3) = box(3)
         !write(6,*) "DEBUG: OLDBOX: ",oldbox(1),oldbox(2),oldbox(3)
         !write(6,*) "DEBUG: OLDabc: ",a,b,c

         ! 2- Read in current box coords.
         if (.not.is_netcdf) then
            read(INPTRAJ_UNIT,'(3f8.3)',end=1000,err=1020) box(1), box(2), box(3)
#ifdef BINTRAJ
         else
            if (NC_error(nf90_get_var(ncid,cellLengthVID,box(1:3), &
                                      start = (/ 1, member /), &
                                      count = (/ 3, 1 /)),&
                         'reading netcdf box coordinates')) goto 1020 
#endif
         endif
         !write(6,*) "DEBUG: BOX: ",box(1),box(2),box(3)

         ! 3- If the box size has changed, some parameters have to be recalc.
         if (box(1)/=oldbox(1).or.box(2)/=oldbox(2).or.box(3)/=oldbox(3)) then

            !write (6,'(a)') "TRAJENE: Updating Box parameters"
            ! 3.a - Update PME cell parameters
            call fill_ucell(box(1),box(2),box(3),alpha,beta,gamma)
            !call fill_tranvec()

            ! 3c- Recompute grid sizes
            ! If the grid sizes change the non-bonded energy will probably be
            ! significantly different than it was in the simulation.
            !
            !   RCFFT needs an even dimension for x direction
            call compute_nfft((a + one)*half   ,newnfft1)
            newnfft1=newnfft1*2
            call compute_nfft(b,newnfft2)
            call compute_nfft(c,newnfft3)
            !write (6,'(a,3i6)') "TRAJENE: New NFFTs: ",newnfft1,newnfft2,newnfft3
            if ( (newnfft1/=nfft1) .or. (newnfft2/=nfft2) .or. (newnfft3/=nfft3) ) then
               write(6,'(a)') 'TRAJENE WARNING: nfft size would change with this box!'
               write (6,'(a,3i6)') "TRAJENE: New NFFTs: ",newnfft1,newnfft2,newnfft3
            endif
         endif ! box sizes have changed

         ! DAN ROE: More Debug for PME
         !write(6,'(/a)') 'Ewald parameters:'
         !write(6,'(5x,4(a,i8))') 'verbose =',verbose, &
         !      ', ew_type =',ew_type,', nbflag  =',nbflag, &
         !      ', use_pme =',use_pme
         !write(6,'(5x,4(a,i8))') 'vdwmeth =',vdwmeth, &
         !      ', eedmeth =',eedmeth,', netfrc  =',netfrc
         !write(6, 9002) a, b, c
         !write(6, 9003) alpha, beta, gamma
         !write(6, 9004) nfft1, nfft2, nfft3
         !write(6, 9006) cutoffnb, dsum_tol
         !write(6, 9007) ew_coeff
         !write(6, 9005) order
         !9002 format (5x,'Box X =',f9.3,3x,'Box Y =',f9.3,3x,'Box Z =',f9.3)
         !9003 format (5x,'Alpha =',f9.3,3x,'Beta  =',f9.3,3x,'Gamma =',f9.3)
         !9004 format (5x,'NFFT1 =',i5  ,7x,'NFFT2 =',i5  ,7x,'NFFT3 =',i5)
         !9005 format (5x,'Interpolation order =',i5)
         !9006 format (5x,'Cutoff=',f9.3,3x,'Tol   =',e9.3)
         !9007 format (5x,'Ewald Coefficient =',f9.5)
      endif ! ntb>0

      write (6,'(a,i6)') 'minimizing coord set #',member
      member=member+1

      call runmin(x,ix,ih,ipairs,x(lcrd),x(lforce),x(lvel), &
            ix(iibh),ix(ijbh),x(l50),x(lwinv),ix(ibellygp), &
            x(l95),ene,carrms,qsetup)

      write (6,364) ene%pot%tot,carrms
      364 format ('minimization completed, ENE=',1x,e14.8, &
            1x,'RMS=',1x,e12.6)

      if (master .and. itgtmd == 1) then
         write (6,'(a,f8.3)') "Final RMSD from reference: ",rmsdvalue
      end if

      ! write the frame to the mdcrd file if the user has set ntwx=1
      ! don't worry about imaging etc since we write the same way it came in
      ! NOTE: Eventually make ntwx work properly, i.e. frames can be skipped

      if (master .and. ntwx >= 1) then
         call corpac(x(lcrd),1,natom*3,MDCRD_UNIT,loutfm)
         if (ntb > 0)  call corpac(box,1,3,MDCRD_UNIT,loutfm)
#ifdef BINTRAJ
         if (is_netcdf) call end_binary_frame(MDCRD_UNIT)
#endif
      !elseif (master) then
      !   write (6,*) "Not writing coordinates to mdcrd due to NTWX value"
      endif

      !       ---loop for next coordinate set

   end do

   !     ---end of trajectory file

   1000 write (6,'(a)') "TRAJENE: Trajectory file ended"
   ok=.true.
   goto 1500

   1010 write (6,'(a)') "TRAJENE: Error reading trajectory coordinates"
   ok=.false.
   goto 1500

   1020 write (6,'(a)') "TRAJENE: Error reading box coordinates"
   ok=.false.
   goto 1500

   1500 write (6,'(a)') "TRAJENE: Trajene complete."
#ifdef BINTRAJ
   if (is_netcdf) call NC_close(ncid)
#endif
   return
end subroutine trajene

end module trajenemod
