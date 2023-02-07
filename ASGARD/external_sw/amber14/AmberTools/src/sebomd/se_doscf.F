      subroutine se_doscf(eelect1,ecore1,etot1,ierror)
!
!     driver routine to carry out the iterative self-consistent field
!     calculation.  energies are returned in ev.
!
!     eelect1 = electronic energy
!     ecore1  = core-core repulsions
!     etot1   = eelect1 + ecore1
!
      implicit none
#include "sebomd.dim"
#include "se_glbmat1.h"
#include "se_glbmat2.h"
#include "se_logic.h"
#include "se_bpairs.h"
#include "se_electr.h"
#include "se_electr1.h"
#include "se_fulsys.h"
#include "se_crit.h"
#include "se_scfstuff.h"
#include "se_timings.h"
#include "se_parproc.h"
#include "se_convrg.h"
#include "se_constants.h"
#include "se_opt.h"
      double precision eelect1, ecore1, etot1
      integer ierror

      logical first,pseudo,doshft,goodf
      data first /.true./
      save first,eeold,ushift,descf,dpscf,etold

      logical frst_call

      double precision ushift, descf, dpscf, ashift
      double precision egap, dii, dij, pold, pnew
      double precision denaux, deprev, dpprev, dpmax
      double precision etold, eeold, eaux2, dpold
      double precision ecrt2, dcrt2, eeaux, dee, daux2, eshift
      integer iimax, ipmax, ijmax
      integer maxit, noscil
      integer k, iter, init, idiag, nocc
      integer ii, ij

      frst_call = first

      if(first)then
        eelect1 = 99999.0d0
        etot1 = 999999.0d0
        etold = 999999.0d0
        ushift = 0.0d0
        descf = -1.0d0
        dpscf = -1.0d0
      endif
      egap = 0.0d0
      ashift = ushift
      doshft = ushift.gt.1.0d-6

      ierror = 0
      iimax = iimat(natoms+1)-1
      ipmax = ip1(natoms+1)
      ijmax = ijmat(ipmax)-1
!
!     energy (ev) and density matrix convergence criteria:
!
      if(descf.gt.0.0d0) eecrt = descf
      if(dpscf.gt.0.0d0) dencrt = dpscf
!
!     override above criteria if density matrix has very good
!     convergence.
!
      eeaux = 10.0d0*eecrt
      denaux = dencrt/10.0d0
      if(descf.gt.0.0d0) eeaux = descf
      if(dpscf.gt.0.0d0) denaux = dpscf

      deprev = -1.0d0
      dpprev = 1.0d0
      dpmax = 1.0d0

      pseudo = .false.
!
!     do scf iterations.
!
      maxit = maxitscf

      noscil = 0
      iter = 0

      ndiag = 0
      do k=1,maxit

!       idiag will be a flag to keep track of when density matrix has
!       been initialized to diagonal form.  if idiag.eq.1 then mosub
!       will be told that there is no good estimate of the fermi energy.
!       also, if fshift is called, then no shift will be applied because
!       p doesn't even come close to commuting with f.
!
        if(first.and.k.eq.1)then
          idiag = 1
        else
          idiag = 0
        endif
!
!       reset density matrix and invoke shift (if it's not already in
!       effect) if we're undergoing large-scale oscillations after 10
!       iterations.
!
        if(iter.gt.10)then

          if(etot1.gt.(etold+1.0d0).and.abs(dpmax).gt.0.1d0)then
            if (screen) then
              write(iout,'(" resetting density matrix, ",
     &                   "invoking shift")')
            endif

            iter = 0
            init = 0
            call se_initp(init, ierror)
            if (ierror.ne.0) return
            idiag = 1
            eelect1 = 99999.0d0
            etot1 = 999999.0d0
            doshft = .true.
            ashift = max(2.0*ashift,ushift,10.0d0)
            deprev = 1.0d0
            dpprev = 1.0d0
            noscil = 0
          endif
        endif
             
        iter = iter + 1
!
!       compute global fock matrix, and, if it's the first scf iteration,
!       compute core-core repulstions.
!
        if (pmeqm) then
          call se_fockpme(k,ecore1)
        else
          call se_fock(k,ecore1)
        endif
!
!       get electronic energy using the current fock matrix and orbitals.
!
        eeold = eelect1
        etold = etot1
        call se_escf(eelect1)
        etot1 = eelect1 + ecore1
        call se_debugmsg1f("doscf eelect1 = ", eelect1*ev2kcal)
        dee = eelect1 - eeold
!       dee = etot1 - etold
! 05/2008  convergence criteria relative to the total energy 
!       dee = (etot1 - etold) / etot1
! 06/2011 search for machine precision
        call se_deeprec(etot1,eecrt)
        eeaux = 10.0d0*eecrt

!
!       check for convergence.
!
!       if(abs(dee).lt.eecrt.and.abs(dpmax).lt.dencrt)then
!         ecrt2 = eecrt*2.0d0
!         dcrt2 = dencrt*2.0d0
!         if(abs(deprev).lt.ecrt2.and.abs(dpprev).lt.dcrt2)then
!           goto 1000
!         endif
!       elseif(abs(dee).lt.eeaux.and.abs(dpmax).lt.denaux)then
!         eaux2 = eeaux*2.0d0
!         daux2 = denaux*2.0d0
!         if(abs(deprev).lt.eaux2.and.abs(dpprev).lt.daux2)then
!           goto 1000
!         endif
!       endif
        if (screen) then
          write(iout,'("SEBOMD: iter=",2i3," E(SCF)=",f20.10," dee=",
     .      f20.10," dpmax=",f20.10," efermi=",f20.10)')
     .      k,iter,etot1, dee,dpmax,efermi(1)
        endif
        if (abs(dpmax).le.dpmaxcrt) then
          goto 1000
        endif
        deprev = dee
        dpprev = dpmax
!
!       see if a pseudo-diagonalization can be done instead of
!       a full diagonalization.
!
        if(stand.and..not.fullscf.and.iter.gt.1)then
!         pseudo = abs(dee/eelect1).lt.0.05d0.and.abs(dpmax).lt.0.1d0
          pseudo = abs(dee/etot1).lt.0.05d0.and.abs(dpmax).lt.0.1d0
        endif
!
!       do a full diagonalization every 5 or 10 iterations, depending
!       on whether a dynamic level shift is being used.
!
        if(pseudo.and.iter.gt.10)then
          if(mod(iter,10).eq.0)then
            pseudo = .false.
          elseif(doshft.and.mod(iter,5).eq.0)then
            pseudo = .false.
          endif
        endif
!
!       keep track of energy oscillations.
!
        if(dee*deprev.lt.0.0d0.and.iter.gt.5) noscil = noscil + 1
!
!       fock matrix mixing scheme to accelerate convergence:
!
!
!     . dynamic or static level shift:
!
        if(doshft)then
          if(iter.eq.1.or..not.stand)then
            egap = 0.0d0
          else
!
!     .     use previous homo-lumo gap if we're doing a
!     .     pseudo-diagonalization.
!
            if(.not.pseudo)then
              nocc = nelec/2
              egap = eval(nocc+1) - eval(nocc)
            endif
          endif
          if(egap.lt.ashift)then
            eshift = 0.5d0*(egap - ashift)
!
!     .     fshift will not actually do a shift if the density matrix
!     .     has just been initialized to diagonal form (idiag.eq.1).
!
            call se_fshift(idiag,eshift)
          endif
        endif
!
!       mosub will build new density matrix, so save old one to check
!       for convergence.
!
        do ii=1,iimax
          piiold(ii) = pdiag(ii)
        enddo
        do ij=1,ijmax
          pijold(ij) = pdiat(ij)
        enddo
!
!       determine mo's and density matrix.
!
!       a good estimate of fermi energy exists if idiag.ne.1
!
        goodf = idiag.ne.1
        call se_mosub(goodf,iter,pseudo,ierror)

        ndiag = ndiag + 1
        if(ierror.ne.0) return
!
!       check for convergence in density matrix and energy.
!

        dpold = dpmax
        dpmax = 0.0d0
        do ii=1,iimax
          dii = pdiag(ii) - piiold(ii)
          if(abs(dii).gt.abs(dpmax)) dpmax = dii
        enddo
        do ij=1,ijmax
          dij = pdiat(ij) - pijold(ij)
          if(abs(dij).gt.abs(dpmax)) dpmax = dij
        enddo
!
!     . density matrix mixing scheme to accelerate convergence:
!

!       call se_pmix(iter)
!
!     . use a 50-50 density matrix mixture on first iteration to
!     . avoid oscillations.
!
        if(doshft.and.iter.eq.1)then

          pold = 0.5d0
          pnew = 0.5d0
          do ii=1,iimax
            pdiag(ii) = pnew*pdiag(ii) + pold*piiold(ii)
          enddo
          do ij=1,ijmax
            pdiat(ij) = pnew*pdiat(ij) + pold*pijold(ij)
          enddo
        endif
      enddo
!     ierror = 1
      write(iout,'(/" no convergence in scf calculation")')
 1000 first = .false.
!     if (wrtedmx) then
!        call se_wrtdmx(ierror)
!        if (ierror.ne.0) return
!     endif
      return
      end
