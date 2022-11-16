! <compile=optimized>
#include "copyright.h"
#include "../include/assert.fh"
#include "../include/dprec.fh"

!*******************************************************************************
!
! Module: scaledMD_mod
!
! Description: 
!
! Module for controlling accelerated molecular dynamics calculations
! Specifically for Scaled MD.
!
! Written by Romelia Salomon-Ferrer, 5/2013
!              
!*******************************************************************************

module scaledMD_mod

  use file_io_dat, only : MAX_FN_LEN, scaledMDLOG_UNIT, scaledMDlog

  implicit none


! Variables
!

!ScaledMD
  integer, save             :: scaledMD
  _REAL_, save              :: scaledMD_lambda
  _REAL_, save              :: scaledMD_energy, scaledMD_weight, scaledMD_unscaled_energy
  integer, save             :: scaledMD_ntwx
 
contains

!*******************************************************************************
!
! Subroutine: scaledMD_setup
!
! Description: Sets up the scaledMD run. Opens the log file, sets up exchange tables
!              etc.
!
!*******************************************************************************

subroutine scaledMD_setup(ntwx)
  

  implicit none
#  include "parallel.h"

! Formal arguments:
  integer               :: ntwx



!  scaledMD_lambda = 1.0
!  scaledMD_energy = 0.0
!  scaledMD_weight = 0.0
!  scaledMD_unscaled_energy = 0.0
  scaledMD_ntwx = ntwx

   ! Open and write out header to scaledMDlog file. Only overall master
   !  deals with the scaledMDlog.
   if (worldrank == 0) then
      call amopen(scaledMDLOG_UNIT,scaledMDlog,'U','F','W')
      write(scaledMDLOG_UNIT,'(a)') "# Accelerated Molecular Dynamics log file"
      write(scaledMDLOG_UNIT,'(a)') "# ntwx,total_nstep, &
            &scaledMD_unscaled_energy,scaledMD_lambda,scaledMD_energy,scaledMD_weight"
   endif

  return
end subroutine scaledMD_setup

!*******************************************************************************
!
! Subroutine:  calculate_scaledMD_total_weights
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine scaledMD_scale_frc(atm_cnt,tot_potenergy,frc)

  implicit none
#  include "parallel.h"

! Formal arguments:
  integer             :: atm_cnt
  _REAL_              :: tot_potenergy
  _REAL_              :: frc(*)

  integer             :: i

!scaledMD
! Scaled energy
  scaledMD_energy = tot_potenergy * scaledMD_lambda
  scaledMD_weight = -(tot_potenergy * (1.0-scaledMD_lambda))
  scaledMD_unscaled_energy = tot_potenergy 

  do i = 1, atm_cnt
    frc(i*3 - 2) = frc(i*3 - 2) * scaledMD_lambda
    frc(i*3 - 1) = frc(i*3 - 1) * scaledMD_lambda
    frc(i*3) = frc(i*3) * scaledMD_lambda
  enddo 
  return

end subroutine scaledMD_scale_frc


!*******************************************************************************
!
! Subroutine:  write_scaledMD_log
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine write_scaledMD_log(ntwx,total_nstep)


  implicit none

! Formal arguments:

  integer               :: ntwx
  integer               :: total_nstep

! Local variables:

  write(scaledMDLOG_UNIT,'(2x,2i10,4f22.12)') ntwx,total_nstep, &
    scaledMD_unscaled_energy,scaledMD_lambda,scaledMD_energy,scaledMD_weight
  
  return

end subroutine write_scaledMD_log

!*********************************************************************
!               SUBROUTINE scaledMD_CLEANUP
!*********************************************************************
!  Close scaledMD files
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine scaledMD_cleanup()

   implicit none
#  include "parallel.h"

   ! Close scaledMDlog file
   if (worldrank==0) then
      close(scaledMDLOG_UNIT)
   endif


end subroutine scaledMD_cleanup

end module scaledMD_mod
