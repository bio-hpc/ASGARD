! <compile=optimized>
#include "copyright.h"
#include "../include/assert.fh"
#include "../include/dprec.fh"

!*******************************************************************************
!
! Module: amd_mod
!
! Description: 
!
! Module for controlling accelerated molecular dynamics calculations
!
! Written by Romelia Salomon-Ferrer, 2/2012
!              
!*******************************************************************************

module amd_mod

  use file_io_dat, only : MAX_FN_LEN, AMDLOG_UNIT, amdlog

  implicit none


! Variables
!

!AMD
  logical                   :: amd_save_all_weights
  integer, save             :: iamd,iamdlag,w_amd
  _REAL_, save              :: w_sign
  _REAL_, save              :: EthreshD,alphaD,EthreshP,alphaP
  _REAL_, save              :: EthreshD_w,alphaD_w,EthreshP_w,alphaP_w
  _REAL_, save              :: totalenergy,fwgt,totdih,dihsum
  _REAL_, save              :: amd_dih_noH, amd_dih_H,fwgtd
  integer, save             :: num_amd_recs, num_amd_lag, amd_ntwx
  _REAL_, allocatable, save :: amd_weights_and_energy(:)
  _REAL_, save              :: E_dih_boost ! Dihedral boost energy in kcal/mol
  _REAL_, save              :: E_total_boost ! Total Potential boost energy in kcal/mol

contains

!*******************************************************************************
!
! Subroutine: amd_setup
!
! Description: Sets up the AMD run. Opens the log file, sets up exchange tables
!              etc.
!
!*******************************************************************************

subroutine amd_setup(ntwx)
  

  implicit none
#  include "parallel.h"

! Formal arguments:
  integer               :: ntwx

! Local variables:
  integer :: alloc_failed


  allocate(amd_weights_and_energy(6), stat = alloc_failed)
  REQUIRE(alloc_failed==0)

  fwgt = 1.0
  fwgtd = 1.0
  E_total_boost = 0.0
  E_dih_boost = 0.0
  num_amd_recs = 1
  num_amd_lag = 0
  amd_ntwx = ntwx

  w_sign = 1.0d0
  if(w_amd .ne. 0) w_sign = -1.0d0

   ! Open and write out header to amdlog file. Only overall master
   !  deals with the amdlog.
   if (worldrank == 0) then
      call amopen(AMDLOG_UNIT,amdlog,'U','F','W')
      write(AMDLOG_UNIT,'(a)') "# Accelerated Molecular Dynamics log file"
      write(AMDLOG_UNIT,'(a)') "# All energy terms stored in units of kcal/mol"
      write(AMDLOG_UNIT,'(a)') "# ntwx,total_nstep,Unboosted-Potential-Energy, &
            &Unboosted-Dihedral-Energy,Total-Force-Weight,Dihedral-Force-Weight,Boost-Energy-Potential,Boost-Energy-Dihedral"
   endif

  return
end subroutine amd_setup

!*******************************************************************************
!
! Subroutine:  calculate_amd_dih_weights
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine calculate_amd_dih_weights(totdih_ene)

  use constants, only: KB
  implicit none

! Formal arguments:
  _REAL_              :: totdih_ene
! Local variables
  _REAL_              :: EV
  _REAL_              :: windowed_factor, windowed_derivative
  _REAL_              :: windowed_exp

! Calculate the boosting weight for amd

  totdih = 0.0d0
  E_dih_boost = 0.0d0
  fwgtd   = 1.0d0
  totdih = totdih_ene

  EV = (EthreshD - totdih_ene) * w_sign 
  if( EV .gt. 0.0d0 )then
!  if(totdih.le.EthreshD)then
    if(num_amd_lag .eq. 0)then
    !EV = EthreshD - totdih_ene
    E_dih_boost = (EV**2)/ &
             (alphaD + EV) ! Dih boost E in kcal/mol
    fwgtd = (alphaD**2)/((alphaD + EV)**2)
      if(w_amd .ne.0 )then
        windowed_exp = exp((EthreshD_w - totdih_ene) * w_sign/alphaD_w)
        windowed_factor = 1.0d0/(1.0d0 + windowed_exp)
        windowed_derivative =  E_dih_boost * windowed_factor * windowed_factor * windowed_exp / alphaD_w
        E_dih_boost = E_dih_boost * windowed_factor
        fwgtd = fwgtd * windowed_factor + windowed_derivative
      endif
    end if
  end if

  return

end subroutine calculate_amd_dih_weights

!*******************************************************************************
!
! Subroutine:  calculate_amd_total_weights
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine calculate_amd_total_weights(atm_cnt,tot_potenergy,totdih_ene,E_boost,frc,temp0)

  use constants, only: KB
  implicit none
#  include "parallel.h"

! Formal arguments:
  integer             :: atm_cnt
  _REAL_              :: tot_potenergy
  _REAL_              :: totdih_ene
  _REAL_, intent(out) :: E_boost
  _REAL_              :: frc(*)
  _REAL_              :: temp0
  _REAL_              :: EV, ONE_KB
  _REAL_              :: windowed_factor, windowed_derivative
  _REAL_              :: windowed_exp

  integer             :: i

!AMD DUAL BOOST CALC START
   if(iamd.gt.0)then
     totalenergy = 0.0d0
     E_total_boost = 0.0d0
     fwgt   = 1.0d0
     totalenergy = tot_potenergy + E_dih_boost
     EV = (EthreshP - totalenergy ) * w_sign
     if (((iamd == 1).or.(iamd == 3)) .and. ( EV .gt. 0.0d0 )) then
    !if (((iamd == 1).or.(iamd == 3)) .and. (totalenergy.le.EthreshP)) then
       if (num_amd_lag .eq. 0) then
       !EV = EthreshP - totalenergy 
 ! Potential Energy boost in Kcal/mol 
       E_total_boost = (EV**2)/ &
                 (alphaP + EV)
       fwgt = (alphaP**2)/((alphaP + EV)**2)
        if(w_amd .ne.0 )then
           windowed_exp = exp((EthreshP_w - totalenergy) * w_sign/alphaP_w)
           windowed_factor = 1.0d0/(1.0d0 + windowed_exp)
           windowed_derivative =  E_total_boost * windowed_factor * windowed_factor * windowed_exp / alphaP_w
           E_total_boost = E_total_boost * windowed_factor
           fwgt = fwgt * windowed_factor + windowed_derivative
        endif
       do i = 1, atm_cnt
         frc(i*3 - 2) = frc(i*3 - 2)*fwgt
         frc(i*3 - 1) = frc(i*3 - 1)*fwgt
         frc(i*3) = frc(i*3)*fwgt
       enddo 
       end if
     end if
     if (worldrank ==0 .and. (num_amd_recs.eq.amd_ntwx)) then
       ONE_KB = 1.0d0 / (temp0 * KB)
       amd_weights_and_energy(1) = tot_potenergy 
       amd_weights_and_energy(2) = totdih_ene
       amd_weights_and_energy(3) = fwgt
       amd_weights_and_energy(4) = fwgtd
       amd_weights_and_energy(5) = E_total_boost
       amd_weights_and_energy(6) = E_dih_boost
!       amd_weights_and_energy(5) = E_total_boost * ONE_KB
!       amd_weights_and_energy(6) = E_dih_boost * ONE_KB
     endif
     E_boost = E_total_boost + E_dih_boost
     num_amd_recs = num_amd_recs +1
     if(num_amd_lag .eq. iamdlag)then
       num_amd_lag = 0
     else
       num_amd_lag = num_amd_lag +1
     endif 
   end if

  return

end subroutine calculate_amd_total_weights


!*******************************************************************************
!
! Subroutine:  write_amd_weights
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine write_amd_weights(ntwx,total_nstep)


  implicit none

! Formal arguments:

  integer               :: ntwx
  integer               :: total_nstep

! Local variables:

  integer               :: istep


  istep = total_nstep - ntwx
  write(AMDLOG_UNIT,'(2x,2i10,6f22.12)') ntwx,total_nstep, &
    amd_weights_and_energy(1), amd_weights_and_energy(2), &
    amd_weights_and_energy(3), amd_weights_and_energy(4), &
    amd_weights_and_energy(5), amd_weights_and_energy(6)
  
  num_amd_recs = 1

  return

end subroutine write_amd_weights

!*********************************************************************
!               SUBROUTINE AMD_CLEANUP
!*********************************************************************
!  Close AMD files
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine amd_cleanup()

   implicit none
#  include "parallel.h"

   integer ier

   ! Close amdlog file
   if (worldrank==0) then
      close(AMDLOG_UNIT)
   endif

! dealocate AMD vectors
   deallocate(amd_weights_and_energy, stat=ier)
   REQUIRE(ier==0)


end subroutine amd_cleanup

end module amd_mod
