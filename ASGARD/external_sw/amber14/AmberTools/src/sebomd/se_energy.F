      subroutine se_energy(newsub,eheat1,ierror)
!
!     driver routine to determine heat of formation eheat1.
!
      implicit none
#include "sebomd.dim"
#include "se_parproc.h"
#include "se_pmeqm.h"
#include "se_logic.h"
#include "se_gradnt.h"
#include "se_esave.h"
#include "se_fulsys.h"
#include "se_subsys.h"
#include "se_frozen.h"
#include "se_param.h"
#include "se_bal_mpi.h"
#include "se_constants.h"

      logical newsub
      double precision eheat1
      integer ierror
#ifdef MPI 
#include "mpif.h"
      integer ijk, ier
#endif

      logical first
      data first /.true./
      save first,hfatms
      double precision dedir(3),derec(3),deelecsr(3)
      integer init, i, iai, ij, j, ji
      double precision hfatms, edir, eself, eclmb, erec, eelecsr
      integer nstepd
      data nstepd /0/
      save nstepd
 
      nstepd = nstepd+1
!
!     density matrix will get initialized to diagonal form on first call.
!
      if(first)then
        init = 0
      else
        init = 1
      endif
!
!
!     update subsystems if automatic subsetting is being used.
!
      if (newsub) then
        do i=1,nsub
          notfroz(i) = .true.
        enddo
!
!       pairlist should be updated for d&c calculation with automatic
!       subsetting, or for any calculation in which the cutbond keyword
!       has been specified, or for the very first call to energy.
!
        if(autosub.or.first)then
!
!         do pre-pairlist setup.
!
          call se_sproc1(ierror)
          if(ierror.ne.0) return
!
!         generate pairlist.
!
          call se_bpair(ierror)
          if(ierror.ne.0)then
            write(iout,'(/" maximum allowed number of bonded atom",
     .              " pairs exceeded at atom ",i5,
     .              /" increase mbpair parameter in sebomd.dim",
     .              " and recompile")') ierror
            return
          endif
!
!       assign pointers to global matrices (some depend on pairlist).
!
          call se_glbpnt(ierror)
          if(ierror.ne.0) return

!         the density matrix should be re-initialized whenever the pairlist
!         is updated because its structure and size depends upon the
!         pairlist.
!
          call se_initp(init,ierror)
          if (ierror.ne.0) return
!
!
!         create arrays that allow quick location of pairlist addesses
!         of specified atom pairs.
!
          call se_ijmake
!
!         get numbers of subsystems shared by each atom pair.
!
          call se_sproc2
        endif
      endif

      if ((nstepd.gt.3).and.ipolyn) then
        call se_dodmxinter ! density matrix interpolation
      endif

!
!     on first call, compute heat of formation of separated atoms.  also
!     set gradient norm to zero so that default convergence criteria
!     will be used in doscf.  return with ierror=-1 if user has specified
!     the 'testrun' keyword.
!

      if(first)then
        first = .false.
        hfatms = 0.0d0
        do i=1,natoms
          iai = iatnum(i)
          if(iai.ne.0) then
            hfatms = hfatms + hfatm(iai) - ev2kcal*eeatm(iai)
          endif
        enddo
      endif

!
!     get scf electronic energy.
!
      if (pmeqm) then
!       compute ewald pair matrix
        do ij=1,natoms*natoms
          phipme(ij) = 0.0d0
          dxphipme(ij) = 0.0d0
          dyphipme(ij) = 0.0d0
          dzphipme(ij) =  0.0d0
          phisr(ij) =  0.0d0
          dxphisr(ij) =  0.0d0
          dyphisr(ij) =  0.0d0
          dzphisr(ij) = 0.0d0
        enddo
        call se_pme_setup
#ifdef MPI 
        do i=myid+1,natoms,nproc
#else
        do i=1,natoms
#endif
          call se_pme_qm_direct(edir,eself,eclmb,dedir,i,i)
          call se_pme_qm_recip(erec,derec,i,i)
          ij=(i-1)*natoms+i
          phipme(ij) = erec + edir
          dxphipme(ij) = derec(1)+dedir(1)
          dyphipme(ij) = derec(2)+dedir(2)
          dzphipme(ij) = derec(3)+dedir(3)
     
          do j=1,i-1
            call se_pme_qm_direct(edir,eself,eclmb,dedir,i,j)
            call se_pme_qm_recip(erec,derec,i,j)
            call se_elecsr(eelecsr,deelecsr,i,j)
            ij=(i-1)*natoms+j
            ji=(j-1)*natoms +i
            phipme(ij) = erec + edir
            phisr(ij)=eelecsr
            phipme(ji)= phipme(ij)
            phisr(ji)=phisr(ij)
            dxphipme(ij) = derec(1)+dedir(1)
            dyphipme(ij) = derec(2)+dedir(2)
            dzphipme(ij) = derec(3)+dedir(3)
            dxphipme(ji) = -dxphipme(ij)
            dyphipme(ji) = -dyphipme(ij)
            dzphipme(ji) = -dzphipme(ij)
           
            dxphisr(ij) = deelecsr(1)
            dyphisr(ij) = deelecsr(2)
            dzphisr(ij) = deelecsr(3)
            dxphisr(ji) = -dxphisr(ij)
            dyphisr(ji) = -dyphisr(ij)
            dzphisr(ji) = -dzphisr(ij)
!           write(59,*) edir , erec,phipme(ij)
          enddo
        enddo
#ifdef MPI 

        call mpi_allreduce(phipme,tmp_mpi,natoms*natoms,
     +       mpi_double_precision,mpi_sum,commsebomd,ier)
        do ijk = 1,natoms*natoms
          phipme(ijk) = tmp_mpi(ijk)
        enddo

        call mpi_allreduce(dxphipme,tmp_mpi,natoms*natoms,
     +       mpi_double_precision,mpi_sum,commsebomd,ier)
        do ijk = 1,natoms*natoms
          dxphipme(ijk) = tmp_mpi(ijk)
        enddo

        call mpi_allreduce(dyphipme,tmp_mpi,natoms*natoms,
     +       mpi_double_precision,mpi_sum,commsebomd,ier)
        do ijk = 1,natoms*natoms
          dyphipme(ijk) = tmp_mpi(ijk)
        enddo

        call mpi_allreduce(dzphipme,tmp_mpi,natoms*natoms,
     +       mpi_double_precision,mpi_sum,commsebomd,ier)
        do ijk = 1,natoms*natoms
          dzphipme(ijk) = tmp_mpi(ijk)
        enddo

        call mpi_allreduce(phisr,tmp_mpi,natoms*natoms,
     +       mpi_double_precision,mpi_sum,commsebomd,ier)
        do ijk = 1,natoms*natoms
          phisr(ijk) = tmp_mpi(ijk)
        enddo

        call mpi_allreduce(dxphisr,tmp_mpi,natoms*natoms,
     +       mpi_double_precision,mpi_sum,commsebomd,ier)
        do ijk = 1,natoms*natoms
          dxphisr(ijk) = tmp_mpi(ijk)
        enddo

        call mpi_allreduce(dyphisr,tmp_mpi,natoms*natoms,
     +       mpi_double_precision,mpi_sum,commsebomd,ier)
        do ijk = 1,natoms*natoms
          dyphisr(ijk) = tmp_mpi(ijk)
        enddo

        call mpi_allreduce(dzphisr,tmp_mpi,natoms*natoms,
     +       mpi_double_precision,mpi_sum,commsebomd,ier)
        do ijk = 1,natoms*natoms
          dzphisr(ijk) = tmp_mpi(ijk)
        enddo
#endif

!        do i=1,natoms
!         do j=1,i-1
!           call se_elecsr(eelecsr,deelecsr,i,j)
!           ij=(i-1)*natoms+j
!           ji=(j-1)*natoms +i
  
!          phisr(ij)=eelecsr
!          phisr(ji)=phisr(ij)

!          dxphisr(ij) = deelecsr(1)
!          dyphisr(ij) = deelecsr(2)
!          dzphisr(ij) = deelecsr(3)
!          dxphisr(ji) = -dxphisr(ij)
!          dyphisr(ji) = -dyphisr(ij)
!          dzphisr(ji) = -dzphisr(ij)
!         enddo
!        enddo

      endif

      call se_doscf(eelect,ecore,etot,ierror)

      if (ierror.ne.0) return
!
!     heat of formation in kcal/mol:
!
      eheat1 =  ev2kcal*etot + hfatms
      nscf = nscf + 1
      return
      end
