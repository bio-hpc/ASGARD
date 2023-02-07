      subroutine se_setup(ierror)
!
!     does all the preliminary setup to prepare the program for the
!     first energy evaluation.
!
      implicit none
#include "sebomd.dim"
#include "se_electr1.h"
#include "se_elemnt.h"
#include "se_endecomp.h"
#include "se_esave.h"
#include "se_fulsys.h"
#include "se_parproc.h"
#include "se_param.h"
#include "se_logic.h"

      integer ierror

      integer norbs, netchg, nzchg
      integer i, iai, iat

      call se_debugmsg0("in setup")

      ierror = 0
!
!     read and process all information supplied in files.
!
      call se_rdall(ierror)
      if(ierror.ne.0) return
      norbs = 0
      netchg = netcharge
      nzchg = 0
      do i=1,natoms
        iai = iatnum(i)
        if(iai.ne.0) then
          nzchg = nzchg + int(zchg(iai))
          norbs = norbs + natorb(iai)
        endif
      enddo
      nelec = nzchg - netchg
      nelecef(1) = nelec

      if(mod(nelec,2).ne.0)then
        ierror = 1
        write(iout,'(/" ODD NUMBER OF ELECTRONS IS NOT ALLOWED FOR",
     .        " A CLOSED-SHELL CALCULATION")')
      endif
      if(nelec.gt.2*norbs)then
        ierror = 1
        write(iout,'(/" NUMBER OF ELECTRONS IS MORE THAN TWICE THE",
     .        " NUMBER OF ATOMIC ORBITALS")')
      endif
      if(ierror.ne.0) return
!
!     assign semiempirical parameters based on the chosen hamiltonian.
!     these are the intramolecular contributions, needed to build
!     up the 1 electron part of the fock matrix, so getpar needs
!     to be called in any case 
!
      call se_getpar
      tempk = 1000.0d0
!
!     check if some parameters are missing (based on Uss)
!
      do i = 1, natoms
        iat = iatnum(i)
        if (ucore(0,iat).eq.0.0d0) then
          write(iout,'("")')
          write(iout,'(" ERROR: NO PARAMETER FOR ATOM ",i6,
     .               " (ATOMIC NUMBER =",i3,")")') i, iat
          call mexit(iout,1)
        end if
      end do
      call se_debugmsg0("out setup")
      return
      end
