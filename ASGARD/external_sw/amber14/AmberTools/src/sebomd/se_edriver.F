      subroutine se_edriver(iflag, nat,
     $     gradsander, ierror,gradpme, symbol2)

      use se_corrections_module, only : Apply_se_corrections
      implicit none
#include "sebomd.dim"
#include "se_logic.h"
#include "se_esave.h"
#include "se_xyzcrd.h"
#include "se_fulsys.h"
#include "se_parproc.h"
#include "se_timings.h"
#include "se_frozen.h"
#include "se_electr.h"
#include "se_electr1.h"
#include "se_constants.h"
#include "se_gradnt.h"
#ifdef MPI
#include "mpif.h"
      integer ier
      integer tmp_mpi_data
#endif

      integer iflag, nat
      integer ierror
      double precision gradsander(*),gradpme(*)
      double precision gradpme1(maxatm)
      double precision gradpme2(maxatm)
      double precision grad1(3,maxatm)
      double precision grade(3*maxatm)

!     locals:
      logical newsub
      integer i, ii
      double precision edir, eself, erec
!  Antoine MARION :: 2013-12-10 Start
      double precision gradSEcorrections(3,3*maxatm)
      double precision vir_SE(3)
      character*4 symbol2(nat)
      double precision E_SEcorrections, E_SEpep
      logical SE_corrections
      double precision vir_se_corrections

      call se_debugmsg0('in edriver')

      SE_corrections = .true.
!  Antoine MARION :: 2013-12-10 End


!     this part is not parallelized, run it by the master node

      ierror = 0

      if (iflag.eq.1) then
        newsub = .true.
      else
        newsub = .false.
      endif


      if (iflag.eq.0) then
        call se_setup(ierror)

        if (SE_corrections) then
          call setupSEcorrections(symbol2,ierror)
        endif
          
        if(ierror.ne.0) return
      endif

      if (nat.ne.natoms) then
        write(iout,'(" ERROR: THE NUMBER OF ATOMS BETWEN SANDER AND",
     &        " SEBOMD IS NOT IDENTICAL")')
        ierror = 1
        return
      endif

      nscf = 0
      call se_setbox
      if (iflag.eq.1) then
        if (pme) then
          call se_pme_setup
          call se_pme_calcb
          call se_pme_calctheta
        endif
      endif
   
      if (iflag.eq.0) return

      if (newsub) then
        call se_gensub(ierror)
        if(ierror.ne.0) return
      endif
      call se_energy(newsub,eheat,ierror)
!  Antoine MARION :: 2013-12-10 Start
      if (SE_corrections) then
        call Apply_se_corrections(XYZ,E_SEcorrections,E_SEpep,
     &                         gradSEcorrections,vir_SE)
        ecore = ecore + E_SEcorrections
        etot  = etot + E_SEcorrections
        eheat = eheat + (E_SEcorrections + E_SEpep)*ev2kcal
      endif
!  Antoine MARION :: 2013-12-10 End
      if(ierror.ne.0) return
      call se_atmchg
      if (pme) then
        if (gradient) then
          call se_gcart(grad1)
          ii=-3
          do i=1,natoms
            ii=ii+3
            grade(ii+1)=grad1(1,i)
            grade(ii+2)=grad1(2,i)
            grade(ii+3)=grad1(3,i)
          enddo
        endif

        do i=1,3*natoms
          gradpme1(i)=0.0d0
          gradpme2(i)=0.0d0
        enddo
        call se_pme_direct2(edir,eself,ecoul,gradpme1)
        call se_pme_recip(erec,gradpme2)

        do i=1,3*natoms               
          gradpme(i)=gradpme1(i)+gradpme2(i)
          grade(i)=grade(i)+gradpme(i)
        enddo   
           
!     .   add PME energies to total energy
!       elr = edir + eself + erec - ecoul
        elr = edir + erec + eself - ecoul
        etot = etot + elr
        
!     .   add PME contribution to heat of formation
!         eheat = eheat + 23.061*elr
        eheat = eheat + eV2kcal*elr
        do i=1,3
          vir(i)=vir(i)+eV2kcal*(-elr-ecoul)/3.0d0
        enddo     
!       do i=1,3*natoms
!         gradpme(i)=23.061*gradpme(i)
!       enddo
        
      endif

      if(screen) then
        write(iout,'(" SEBOMD FINAL QUANTITIES:",
     .           /" -----------------")')
#ifdef MPI
        call mpi_allreduce(nfull,tmp_mpi_data,1,MPI_INTEGER,MPI_SUM,
     .                      commsebomd,ier)
        nfull = tmp_mpi_data
#endif
        write(iout,
     $        '(" NUMBER OF FULLY DIAGONALIZED SUBSYSTEMS: ",I6,
     &        /" NUMBER OF FROZEN SUBSYSTEMS:             ",I6,
     &        /" NUMBER OF SCF CALCULATIONS = ",i6,
     &        /" NUMBER OF DIAGONALIZATIONS = " ,i6)')
     &        nfull, nfroz, nscf, ndiag
        write(iout,'(/" ELECTRONIC ENERGY     = ",f25.8," EV",
     .        /" CORE-CORE REPULSIONS  = ",f25.8," EV")')
     .        eelect,ecore
        if (pme) then
          write(iout,'("")')
          write(iout,'(" PME DIRECT ENERGY     = ",f25.8," EV")') edir
          write(iout,'(" PME RECIPROCAL ENERGY = ",f25.8," EV")') erec
          write(iout,'(" PME SELF ENERGY       = ",f25.8," EV")') eself
          write(iout,'(" CLAS. COULOMB ENERGY  = ",f25.8," EV")') ecoul
          write(iout,'(" TOTAL LONG RANGE      = ",f25.8," EV")') elr
          write(iout,'("")')
          write(iout,'(" PME DIRECT ENERGY     = ",f25.8," kcal/mol")')
     .      edir*ev2kcal
          write(iout,'(" PME RECIPROCAL ENERGY = ",f25.8," kcal/mol")')
     .      erec*ev2kcal
          write(iout,'(" PME SELF ENERGY       = ",f25.8," kcal/mol")')
     .      eself*ev2kcal
          write(iout,'(" CLAS. COULOMB ENERGY  = ",f25.8," kcal/mol")')
     .      ecoul*ev2kcal
          write(iout,'(" TOTAL LONG RANGE      = ",f25.8," kcal/mol")')
     .      elr*ev2kcal
          write(iout,'("")')
!    &        edir,erec,eself,ecoul, elr
!    &        pmedir,pmerec,pmeself,ecoul, elr
        endif
        write(iout,'(" TOTAL ENERGY          = ",f25.8," EV",
     .        /" HEAT OF FORMATION     = ",f25.8," KCAL/MOL",
     .        /" VIRIAL                = ",f25.8)')
     .        etot,eheat,vir(4)
        write(iout,'(" FERMI ENERGY          = ",f25.8," EV")')
     &           efermi(1)
        IF(ionpot) write (iout,
     $        '(" IONIZATION POTENTIAL  = ",f25.8," EV")')
     $        -eval(nelec/2)
        write(iout,'(" -----------------")')
      endif

      if(GRADIENT)then
        if (.not.pme) then
          call se_gcart(grad1)
          ii=-3
          do i=1,natoms
            ii=ii+3
            grade(ii+1)=grad1(1,i)
            grade(ii+2)=grad1(2,i)
            grade(ii+3)=grad1(3,i)
          enddo
        endif
!  Antoine MARION :: 2013-12-10 Start
        if (SE_corrections) then
          ii=-3
          do i=1,natoms
            ii=ii+3
            grade(ii+1)=grade(ii+1)+gradSEcorrections(1,i)
            grade(ii+2)=grade(ii+2)+gradSEcorrections(2,i)
            grade(ii+3)=grade(ii+3)+gradSEcorrections(3,i)
          enddo
          vir_se_corrections = vir_se(1)+vir_se(2)+vir_se(3)
        endif
!  Antoine MARION :: 2013-12-10 End
        vir(4)=vir(1)+vir(2)+vir(3)+vir_se_corrections
        if (ierror.ne.0) return
      endif
     
      do i=1,3*natoms
        gradsander(i) = grade(i)
      enddo

      if (pole) then
        call se_dipole
      endif
      call se_debugmsg0('out edriver')
      return
      end
