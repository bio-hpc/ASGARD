      subroutine setupSEcorrections(symbol2,ierror)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!     Author   :: Antoine MARION 
!!     Date     :: 2013-12-10
!!     Function :: Setup the se_correction_module to add
!!                 correction to the SemiEmpirical Hamiltonians
!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      use se_corrections_module
      implicit none
#include "sebomd.dim"
#include "se_constants.h"
#include "se_elemnt.h"
#include "se_fulsys.h"
#include "se_logic.h"
#include "se_param.h"
#include "se_xyzcrd.h"
      integer ierror
      character*4 symbol2(natoms)

      integer ctype
      logical mm_pep, let_pm3
      double precision k_pep

!     mm_pep  = .false.
      let_pm3 = .false. ! exit if no intermolecular PIFx parameters
!      let_pm3 = .true.
!     ctype = 0
!     k_pep = 9.8526d0
      call se_corrections_info_from_sander(mm_pep, ctype, k_pep)

      call SEcorr_init(natoms,ctype,mm_pep,k_pep,let_pm3,
     &             iatnum,zchg,symbol2,
     &             agaus,bgaus,cgaus,pbc,
     &             xyz,hartree2eV,Bohr2Ang,ev2kcal,ierror)

      end subroutine setupSEcorrections

