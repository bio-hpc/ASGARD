! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"

!Interface module to provide option for doing QM/MM with divcon.
!Written by Ross Walker (TSRI, 2005) from initial implementation by SAH

#if defined(MPI) || defined(NO_SANDER_DIVCON)
subroutine qm_div()
!DUMMY DIVCON ROUTINE
  call sander_bomb('qm_div','IDC>0 Error','Sander was compiled without support for Divcon.')
  return
end subroutine qm_div

#else

subroutine qm_div(x, ix, f, escf, atom_type)

  use qmmm_module, only : qmmm_div, qmmm_struct, qmmm_nml, allocate_qmmm_pair_list

  implicit none

#include "../include/memory.h"

!Passed in
  _REAL_ x(*), f(*), escf
  integer ix(*)
  character(len=4), intent(in) :: atom_type(*)

!Local
  integer ier
  integer i, j, k, m, ir, jr
  integer mm_no,lnk_no, qm_no
  _REAL_ :: forcemod(3)
  integer jjp,jlnk,iyes,iii,ii,ntotqmres,jlnkindx,jjpindx

  integer, dimension(2*nres) :: npqmres !nres is too much, but at least
  integer, dimension(nres) :: lnkres  !its dynamic-ish now
  integer, dimension(nres) :: nqmres

  _REAL_ , dimension(2,3) :: bxbnd
  _REAL_ , dimension(3*qmmm_struct%nlink) :: qmlink


  if (qmmm_struct%qm_mm_first_call) then
        !qm_xcrd only actually needs to be 4,qm_mm_pairs long...
        allocate (qmmm_struct%qm_xcrd(4,natom), stat=ier )
        !REQUIRE(ier == 0) !Dealocated in deallocate qmmm
        call allocate_qmmm_pair_list( natom )
        allocate ( qmmm_struct%qm_coords(3,qmmm_struct%nquant+qmmm_struct%nlink), stat=ier )
                   !Stores the REAL and link atom qm coordinates
        !REQUIRE(ier == 0)
        allocate(qmmm_struct%dxyzcl(3,natom),stat=ier)
        allocate(qmmm_struct%dxyzqm(3,qmmm_struct%nquant_nlink),stat=ier)
        ntotqmres=1
        allocate(qmmm_struct%jqatms(qmmm_struct%nquant_nlink),stat=ier)
  endif
  if(qmmm_nml%idc == 2)then
!      qmmm_struct%ipres = ix(i02)
      k = 1
      jjp = 1
      jjpindx = 1
      jlnkindx = 1
      jlnk = 1
      npqmres(1) = 1
      do m=1,nres
         lnkres(k) = 0
         iyes = 0
         ii = ix(m) !assign ii to first number in residue i
         iii = ix(m+1)-1 !assign iii to last number in residue i
         do j=ii,iii
            if(jjpindx > qmmm_struct%nquant)jjpindx=qmmm_struct%nquant
            if(j == qmmm_struct%iqmatoms(jjpindx)) then !current atm in QM region
               iyes = 1
               jjp = jjp + 1   !atom(not link atm) counter
               jjpindx = jjpindx + 1
            end if
            if(qmmm_struct%nlink .gt. 0)then
               if(jlnk.gt.qmmm_struct%nlink) jlnkindx=qmmm_struct%nlink
               if(j == qmmm_struct%link_pairs(1,jlnkindx)) then !current MM atm is touching QM region
                  lnkres(k) = lnkres(k) + 1 !sah 4-6-05 increment lnkres so residue can't have more than one link atom
                  jlnk = jlnk + 1           !link atom needed counter
                  jlnkindx = jlnkindx + 1
              if(lnkres(k) .gt. 1) then
!                  WRITE(*,*)"atom ",j,"residue ",m
!                  call sander_bomb('resnba<resnba.f>',&
!                      'there is a residue(',m,') with two link atoms',&
!                      'adjust QM atoms and rerun')
               end if !sah 4-6-05 added bomb so that one residue couldn't have two link atoms
               end if
            end if
         end do  !end loop through residue
         if(iyes == 1) then
            nqmres(k) = i  !tells res# in pdb file
            k = k + 1
            npqmres(k) = jjp+jlnk-1 !tells first atm of res in divcon w/link atm
         end if
         !bw if we have run through all qm atoms, just exit the loop.
         if(jjp >= qmmm_struct%nquant) exit
      end do
      ntotqmres = k - 1
  endif  !if(idc==2)
           call qm_fill_qm_xcrd( x, natom, qmmm_struct%scaled_mm_charges)
           call position_link_atoms(x)

!Print link atom coords to mdout file
           if(qmmm_struct%qm_mm_first_call)then
              call print_link_atom_info(qmmm_struct%qm_coords,atom_type)
           endif

      do ir=1, qmmm_struct%nlink
         do jr=1,3
         qmlink((ir-1)*3+jr)=                                   &
          qmmm_struct%qm_coords(jr,qmmm_struct%nquant+ir)
         enddo
      enddo

      if(qmmm_struct%qm_mm_first_call .and. qmmm_nml%writepdb)then
         call qm_write_pdb('qmmm_region.pdb')
      endif

      qmmm_struct%jqatms(1:qmmm_struct%nquant) = qmmm_struct%iqmatoms(1:qmmm_struct%nquant)
      qmmm_div%ntotatm = natom

      call qm_mm_div(x,qmlink, &
             f,escf,qmmm_nml%idc,npqmres,lnkres,     &
             ntotqmres,qmmm_struct%nlink)
!*********distribute link atom forces*************
      do k=1,qmmm_struct%nlink
        mm_no = 3*qmmm_struct%link_pairs(1,k)-2  !location of atom in qm pair list
        lnk_no = qmmm_struct%link_pairs(2,k) !Nquant number of QM atom bound to link atom
        qm_no = 3*qmmm_struct%iqmatoms(lnk_no)-2

        !Note this routine uses the flink in the form -flink. 
        call distribute_lnk_f(forcemod,qmmm_struct%dxyzqm(1,qmmm_struct%nquant+k),x(mm_no), &
                              x(qm_no),qmmm_nml%lnk_dis)

        !NOTE: forces are reversed in QM calc with respect to amber force array
        !so we subtract forcemod from MM atom and add it to QM atom.
        m = (qmmm_struct%link_pairs(1,k)-1)*3 !Natom number of MM link pair.
        !MM atom's new force = FMM(x,y,z) - FORCEMOD(x,y,z)
        f(m+1) = f(m+1) - forcemod(1)
        f(m+2) = f(m+2) - forcemod(2)
        f(m+3) = f(m+3) - forcemod(3)

        m = (qmmm_struct%iqmatoms(lnk_no)-1)*3
        !QM atom's new force = FQM(x,y,z) - Flink(x,y,z) + FORCEMOD(x,y,z)
        !Note QM forces should be subtracted from sander F array to leave total force.
        f(m+1) = f(m+1) - qmmm_struct%dxyzqm(1,qmmm_struct%nquant+k) + forcemod(1)
        f(m+2) = f(m+2) - qmmm_struct%dxyzqm(2,qmmm_struct%nquant+k) + forcemod(2)
        f(m+3) = f(m+3) - qmmm_struct%dxyzqm(3,qmmm_struct%nquant+k) + forcemod(3)

      end do
!**********end of force distribution***************


  return

end subroutine qm_div

#endif

