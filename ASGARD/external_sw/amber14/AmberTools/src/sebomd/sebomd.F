!------------------------------------------------------------------------
      subroutine sebomd_energy(
     &  iflag, iflagch, ntb, nat, virial, edc,
     &  boxsander, xyzsander, gradsander,
     &  nresid, presid, chgsander, symbol2,
     &  nmol,pmolid,
     &  pdmx)

!     iflag =  -1  error occured
!     .         0  initialize sebomd variables, do not subset
!     .            and do not calculate energy (occurs only once)
!     .         1  subset and calculate energy
!     .         2  do no subset (use old subsetting scheme) and calculate
!     .            energy
!
!     iflagch = 0 do not write charges to sebomd.ch
!     .         1 write charges to sebomd.ch
!
!     ntb = 0 do not use periodic boundaries conditions
!     .   > 0 use PBC = .true.
!
!     nat = the number of atoms
!
!     virial = the virial (in kcal/mol)
!
!     edc = total energy (in kcal/mol)
!
!     boxsander = box dimensions (3 dim.)
!
!     xyzsander = the coordinates, format xyzsander(3*nat)
!
!     gradsander = the gradient, format gradsander(3*nat) (in kcal/A)
!
!     nresid = number of residues
!
!     presid = residue pointers
!
!     chgsander = starting charges used when PME is on
!
!     symbol2 = atom types of the system
!
!     nmol = number of molecules
!
!     pmolid = array: number of atoms per molecule (size=nmol)
!
!     pdmx = 0 do not print density matrix
!          = 1 print density matrix to sebomd.dmx file
!

      use se_inter, only : se_inter_init
      implicit none
#include "sebomd.dim"
#include "se_logic.h"
#include "se_xyzcrd.h"
#include "se_bpairs.h"
#include "se_esave.h"
#include "se_pbcbox.h"
#include "se_subpar.h"
#include "se_parproc.h"
#include "se_fulsys.h"
#include "se_resdef.h"
#include "se_glbmat1.h"
#include "se_glbmat2.h"
#include "se_glbmat4.h"
#include "se_constants.h"
#include "se_opt.h"
#ifdef MPI
#include "mpif.h"
#endif
      integer iflag, iflagch, ntb, nat
      integer nresid
      double precision edc
      double precision xyzsander(3*nat), boxsander(3), gradsander(3*nat)
      double precision gradpme(3*nat), chgsander(nat), virial(4)
      logical first
      integer presid(maxres), pdmx
      integer nmol, pmolid(nmol)
      character*4 symbol2(nat)
      data first /.true./
      save first

      integer i, ierror, j, ii

      ierror = 0

      call se_init(nat, boxsander, xyzsander, gradsander, chgsander,
     &             maxatm, natoms, xyz, chgpme, dbox, dhalf, boxvol)

      if (iflag.eq.0) then
!       initialize stuff
        pole = .true.
        call se_info_from_sander(netcharge,
     &    ipolyn, fullscf, screen, wrtscr, prtsub,
     &    pme, pmeqm, mullewald, cmewald, chewald,
     &    stand, clust, resclust, smallclust,
     &    debug,
     &    mndo, am1, pm3, rm1, am1d, pm3pddg,
     &    nncore,  dcbuff1, dcbuff2,
     &    diag_routine, dpmaxcrt)


#ifdef MPI
        call se_mpi_vars(nproc, myid, commsebomd)
        if ( myid .ne. 0) then
          open(unit=6,file='/dev/null')
        endif

        if (nproc.gt.maxproc) then
          write(iout,'(" ERROR: NUMBER OF PROCESSORS IS GREATER",
     &           " THAN ",
     &           "MAXPROC",
     &           /"        INCREASE MAXPROC IN sebomd.dim ",
     &           "AND RECOMPILE")')
          ierror = 1
          goto 1000
        endif
#endif

        master = myid.eq.0

        call se_debugmsg0("in sebomd")

        ierror = 0
        call se_setunit

        if (ierror.ne.0) goto 1000
        autosub = clust.or.resclust

!       set the following variables by default
        gradient = .true.

!       periodic boundary conditions
        if (ntb.ne.0) then
          pbc = .true.
        else
          pbc = .false.
        endif

!       type of atoms
        call se_rdelem(symbol2, ierror)

!       information about residues

        if (smallclust) then
          nres = 1
          irpnt(nres) = 1
          do i = 2, natoms
!           if (iatnum(i) /= 1) then
            if ((iatnum(i) == 8).or.    ! start of a water molecule
     .          ((iatnum(i) /= 1).and.(irpnt(nres)+1 < i))) then
              nres = nres + 1
              irpnt(nres) = i
              if(nres .gt. maxres) then
                ierror = 1
                WRITE(IOUT,
     .             '(/" MAXIMUM ALLOWED NUMBER OF RESIDUES REACHED",
     .           /" -- INCREASE MAXRES PARAMETER IN sebomd.dim",
     .            " AND RECOMPILE")')
                ierror = 1
                goto 1000
              endif
            endif
          end do
        else
          nres = nresid
          if(nres .gt. maxres) then
            ierror = 1
            WRITE(IOUT,
     .         '(/" MAXIMUM ALLOWED NUMBER OF RESIDUES REACHED",
     .       /" -- INCREASE MAXRES PARAMETER IN sebomd.dim",
     .        " AND RECOMPILE")')
            ierror = 1
            goto 1000
          endif
          do i=1,nres
            irpnt(i) = presid(i)
          enddo
        endif
        irpnt(nres+1) = nat + 1

!       information about molecules
        call se_inter_init(natoms,nmol,pmolid,symbol2)

        if (ierror.ne.0) goto 1000

!       If problem with convergence, use a new density matrix

        if (.not.first) guess = .false.
        call se_initp(0,ierror)
      endif   ! iflag = 0

      call se_debugmsg0("in sebomd (after iflag = 0 section")

!     if (pbc) then
!       call se_backinbox(natoms, nres, xyz, dbox, irpnt)
!     endif

      call se_edriver(iflag, nat,
     $        gradsander, ierror, gradpme, symbol2)
      do i=1,4
        virial(i) = vir(i)
      enddo
      if (pmeqm) then
        do i=1,3
          virial(i) = virial(i) - virlrqm*eV2kcal/3.0d0
        enddo
      endif
      edc = eheat
      if (ierror.ne.0) then
        write(iout,'(/" ERROR IN CALCULATING SEBOMD ENERGY")')
        goto 1000
      endif

#ifdef MPI
      if (myid .eq. 0) then
#endif

        if (iflagch.eq.1.and.iflag.ge.1) then
          call se_write_charges(natoms, edc, iatnum, xyz, atchg,
     .                              atchg2, atchg3, dipole_atm)
        endif
#ifdef MPI
      endif
#endif

      if (pdmx .eq. 1) call se_wrtdmx(ierror)

      if ((iflag.eq.1) .and. ipolyn) then
        call se_updatedmxinter()
      endif ! iflag = 1 and ipolyn

 1000 continue
      first = .false.
      if(ierror.ne.0) then
        write(iout,'(/" SEBOMD ERROR")')
        iflag = -1
        call se_debugmsg0("out sebomd (with ierror != 0")
        return
      endif


      call se_debugmsg0("out sebomd")
      call se_debugmsg1f("sebomd: energy=",edc)
      return
      end
