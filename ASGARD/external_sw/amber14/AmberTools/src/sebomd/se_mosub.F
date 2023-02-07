      subroutine se_mosub(goodf,iter,pseudo,ierror)
!
!     determines the molecular orbitals for all subsystems.
!
!     goodf =  flag for a good estimate of the fermi energy.  when
!              the density matrix has just been initialized to
!              diagonal form, then we don't have a good estimate
!              of efermi.  in this case, two rounds of diagonalizations
!              will be necessary to determine the global density
!              matrix.
!
!     iter   = scf iteration number: 1,2,...
!
!     pseudo = logical flag do carry out pseudo-diagonalization
!
!     returned error codes:
!
!     ierror = 0 --> okay
!     ierror = 1 --> (1) qr iteration in diag didn't converge, or
!                    (2) fermi energy in doferm didn't converge
!
!     note: this routine uses the /work/ common block for temporary
!           storage.  anything in ff, evec, or ww below will be destroyed.
!
!
!
      implicit none
#include "sebomd.dim"
#include "se_bpairs.h"
#include "se_glbmat1.h"
#include "se_glbmat2.h"
#include "se_logic.h"
#include "se_parproc.h"
#include "se_subsys.h"
#include "se_fulsys.h"
#include "se_electr.h"
#include "se_electr1.h"
#include "se_frozen.h"
#include "se_work.h"
#include "se_pass1.h"
#include "se_bal_mpi.h"
#include "se_elemnt.h"
#include "se_opt.h"
      integer iter, ierror
      logical goodf,pseudo

#ifdef MPI 
#include "mpif.h"
      integer ier, ipmax, num_pdiag, num_pdiat
#endif
      integer npass, ipass, iflag
      integer iorbmx, i, k, iatm, norbsi, ijglb
      integer iiorb, ijloc0, ijloc, if, jf
      integer npairs, natms, ii, jj, j, jatm, norbsj
      integer nk, nocc, k1, ijk, iimax
      integer norbsk, iatm0, iatm1, ijaddr
      double precision arg, deltae, tolera, betah, efermix
      logical first
      data first /.true./
      save first
!
!     bolzmann constant in ev/kelvin:
!
      double precision boltz
      data boltz /8.617335408d-5/
      save boltz

      betah = 1.0d0/(boltz*tempk)
      nfull = 0
      nfroz = 0
      iflag = 0  ! todo: check and remove need for iflag in this subroutine

      iimax = iimat(natoms+1)-1

#ifdef MPI 
! JV Init eval in case we changed subsetting and they
! JV have garbage from other subsystems - OK for serial version
      iorbmx = iorbpt(nsub+1)-1
      do i=1,iorbmx
        eval(i) = 0.0d0
        evecsq(i) = 0.0d0
      enddo
      if ((iter.eq.1.) .and. (iflag.eq.0)) then
        do i=1,iorbmx
          fermi(i) = 0.0d0
        enddo
      endif

!C JV Balance by number of orbitals per subsytem as initial guess
      call se_balance_bynorbs
!CJV Init timers for next load balancing (by completion time)
!     do i = 1,nsub
!       time_subs(i) = 0.0
!     enddo
#endif

!     when no reliable estimate of efermi exists, then two diagonalization
!     passes are necessary.  on the first pass, the eigenvectors are
!     used to get the quantities evecsq, which are in turn used to
!     determine the fermi energy.  then the eigenvectors are determined
!     again, and the global density matrix is constructed.  this approach
!     avoids the necessity of storing all the subsystem eigenvectors
!     simultaneously.  under normal circumstances, multiple passes are
!     only necessary on the first call at program startup.  on later
!     scf iterations and for later geometries, the fermi energy from the
!     previous call may be used.
!
      if(goodf.or.stand.or.nsub.eq.1)then
        npass = 1
      else
        npass = 2
      endif

      do ipass=1,npass
!
!       have diag compute eigenvalues only on the first pass.
!
        doeval = ipass.eq.1
!
!       assemble fock matrix for each subsystem by accessing the correct
!       blocks of the global fock matrix.
!

        tolera = 1.0d-16
 
#ifdef MPI 
        do ijk= 1,my_numsubs
           if (iter .eq. 1) call se_etimer(watch(15))
           k = my_subs(ijk)
#else
        do k=1,nsub
#endif
!
!         diagonal blocks first.
!
          norbsk = iorbpt(k+1)-iorbpt(k)
          iiorb = 0
          ijloc0 = 0
          do i=iatom1(k),iatom1(k+1)-1
            iatm = iatoms(i)
            norbsi = natorb(iatnum(iatm))
            ijglb = iimat(iatm)
            do if=1,norbsi
              ijloc = ijloc0 + if + iiorb
              do jf=1,if
                ff(ijloc) = fdiag(ijglb)+fpmediag(ijglb)
                ijglb = ijglb + 1
                ijloc = ijloc + norbsk
              enddo
            enddo
            iiorb = iiorb + norbsi
            ijloc0 = iiorb*norbsk
          enddo
!
!         now off-diagonal blocks.
!
          npairs = ip1(natoms+1)-1
          natms = iatom1(k+1) - iatom1(k)
          if(natms.gt.1)then
            iatm0 = iatom1(k) - 1
            iatm1 = iatoms(iatom1(k))
            ii = natorb(iatnum(iatm1))
            do i=2,natms
              iatm = iatoms(iatm0+i)
              norbsi = natorb(iatnum(iatm))
              jj = 0
              ijloc0 = 0
              do j=1,i-1
                jatm = iatoms(iatm0+j)
!
!               get position in pairlist of (iatm,jatm) pair.
!
                call se_ijfind(npairs,iatm,jatm,ijaddr)
                norbsj = natorb(iatnum(jatm))
!
!               copy global diatomic block fdiat to ff, but only if iatm
!               and jatm are close enough to bond.  otherwise, set the
!               block of ff to zero.
!
                if(ijaddr.ne.0)then
                  ijglb = ijmat(ijaddr)
                  do if=1,norbsi
                    ijloc = ijloc0 + ii + if
                    do jf=1,norbsj
                      ff(ijloc) = fdiat(ijglb)+fpmediat(ijglb)
                      ijloc = ijloc + norbsk
                      ijglb = ijglb + 1
                    enddo
                  enddo
                else
                  do if=1,norbsi
                    ijloc = ijloc0 + ii + if
                    do jf=1,norbsj
                      ff(ijloc) = 0.0d0
                      ijloc = ijloc + norbsk
                    enddo
                  enddo
                endif
                jj = jj + norbsj
                ijloc0 = jj*norbsk
              enddo
              ii = ii + norbsi
            enddo
          endif
          nk = norbsk
          nocc = nelec/2
          k1 = iorbpt(k)

          nfull = nfull + 1
!     .   either a pseudo diagonalization or a full diagonalization
          if(pseudo)then
             call se_diagp(nk,nocc,ff,eval(k1),ww,vv,evec)
          else
             if (diag_routine.eq.1) then
              call se_diag(nk,ff,nk,tolera,ww,eval(k1),
     &            idegen,evec,ierror)
              if(ierror.ne.0)then
                 write(iout,'(/" ERROR IN DIAG -- NO CONVERGENCE IN",
     .               " QR ITERATION")')
                 return
              endif
             elseif (diag_routine.eq.2.) then
               do ii = 1,nk*nk
                 evec(ii) = ff(ii)
               end do
               call dsyev('V', 'L', nk, evec, nk, eval(k1),
     &                    ww,nk*nk, ierror)
               if(ierror.ne.0)then
                  write(iout,'(/" ERROR IN DSYEV")')
                  return
               endif
             else
               write(6,'("Error in mosub.F")')
               write(6,'("diag_routine should be <= 1")')
               write(6,'("Exiting...")')
               call mexit(6,1)
             endif
          endif

!
!         compute global density matrix for standard calculation,
!         or, if multiple subsystems and final pass, add subsystem
!         contribution to global density matrix.  if, for some
!         strange reason, there is only one subsystem, and the user
!         has not specified 'standard', then the density matrix
!         cannot be determined until after the call to doferm.
!
          if(stand)then
             call se_denful(norbsk,evec)
          elseif(ipass.eq.npass.and.nsub.gt.1)then
!
!           if automatic subsetting is being used, then the pointers
!           in the array fermi can change from one geometry to the
!           next.  to be safe, recompute the fermi occupation numbers
!           on the first scf iteration of all geometries after the
!           initial one.
!

                efermix = efermi(1)
             if(autosub.and..not.first.and.iter.eq.1)then
              do i=iorbpt(k),iorbpt(k+1)-1
                deltae = eval(i) - efermix
                arg = betah*deltae
                if(arg.gt.50.0d0)then
                  fermi(i) = 0.0d0
                elseif(arg.lt.-50.0d0)then
                  fermi(i) = 1.0d0
                else
                  fermi(i) = 1.0d0/(1.0d0 + exp(arg))
                endif
              enddo
            endif

            call se_densub(k,norbsk,evec)
          endif
!
!         for the first pass of a d&c calculation, fill in this
!         subsystem's contribution to evecsq.  when the whole
!         evecsq array has been assigned, the fermi energy may be
!         determined.
!
          if(ipass.eq.1.and..not.stand)then
            call se_esqr(k,norbsk,evec)

! jv moved call to doferm outside 500 loop
!
!           if this is the last subsystem, then we can determine
!           the fermi energy.
!
          endif
!
!         add to global density matrix for the strange situation where
!         'standard' has not been defined, but there is only one
!         subsystem.
!
          if(.not.stand.and.nsub.eq.1.and.ipass.eq.npass)then
            call se_doferm(ierror)
            call se_densub(k,norbsk,evec)
          endif

#ifdef MPI 
! JV Time susbystems for balancing only on first iteration
          if (iter .eq. 1) then
             call se_etimer(watch(16))
             time_subs(k) = watch(16)-watch(15)
          endif
#endif
        enddo

#ifdef MPI 
!-RDC C JV Gather fermi and eval so we can rebalance and still keep all info
!        if ( iter .eq. 1 ) then

        call mpi_allreduce(fermi,tmp_mpi,iorbmx,
     +           mpi_double_precision,mpi_sum,commsebomd,ier)
        do ijk = 1,iorbmx
          fermi(ijk) = tmp_mpi(ijk)
        enddo

        call mpi_allreduce(eval,tmp_mpi,iorbmx,
     +           mpi_double_precision,mpi_sum,commsebomd,ier)
        do ijk = 1,iorbmx
          eval(ijk) = tmp_mpi(ijk)
        enddo

        call mpi_allreduce(evecsq,tmp_mpi,iorbmx,
     +           mpi_double_precision,mpi_sum,commsebomd,ier)
        do ijk = 1,iorbmx
          evecsq(ijk) = tmp_mpi(ijk)
        enddo
!
! jv load balance based on time once we've done a single pass
!           call mpi_allreduce(time_subs,tmp_mpi,nsub,
!    +           mpi_double_precision,mpi_sum,commsebomd,ier)
!           do i = 1,nsub
!              time_subs(i) = tmp_mpi(i)
!           enddo
!            call se_balance_bytime
!         endif
#endif

        if(ipass.eq.1.and..not.stand .and. (iflag .eq. 0))then
! jv init all fermi values to zero in case we rearranged subsystems
! jv ok to leave this for serial version
          iorbmx = iorbpt(nsub+1)-1
          do ijk = 1,iorbmx
            fermi(ijk) = 0.0
          enddo
          call se_doferm(ierror)
          if(ierror.ne.0)then
            write(iout,'(/" ERROR IN DOFERM -- NO CONVERGENCE",
     .                       " IN FERMI ENERGY DETERMINATION")')
            stop
          endif
        endif
      enddo

#ifdef MPI 
! JV  Sum PDIAG and PDIAT from all PE's
      num_pdiag = iimat(natoms+1)-1
      ipmax = ip1(natoms+1)
      num_pdiat = ijmat(ipmax)-1
      call mpi_allreduce(pdiag,tmp_mpi,num_pdiag,
     +   mpi_double_precision,mpi_sum,commsebomd,ier)

      do ijk = 1,num_pdiag
        pdiag(ijk) = tmp_mpi(ijk)
      enddo
      call mpi_allreduce(pdiat,tmp_mpi,num_pdiat,
     +   mpi_double_precision,mpi_sum,commsebomd,ier)
      do ijk = 1,num_pdiat
        pdiat(ijk) = tmp_mpi(ijk)
      enddo
#endif
!     the first call to mosub is over:
      first = .false.
      return
      end
