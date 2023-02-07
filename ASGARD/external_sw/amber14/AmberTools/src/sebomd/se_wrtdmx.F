      subroutine se_wrtdmx(ierror)
!     By Arjan van der Vaart

      implicit none
#include "sebomd.dim"
#include "se_bpairs.h"
#include "se_elemnt.h"
#include "se_fopen.h"
#include "se_fulsys.h"
#include "se_glbmat1.h"
#include "se_glbmat2.h"
#include "se_parproc.h"
      integer ierror

      integer i, iai, norbsi, ij, ii, k
      integer ie, ipmax

!     text file instead of binary, so easy portable between
!     different architectures

      if (myid.eq.0) then
         call se_opnfil(1,ierror)
         if (ierror.ne.0) then
            write(iout,'(/" ERROR: COULD NOT OPEN ",A20)') fname(1)
            return
         endif
      endif

!     number of atoms
      write(idmx,'(i8)') natoms

!     atomkinds
      write(idmx,'(10i8)') (iatnum(i),i=1,natoms)

!     pairlist
      write(idmx,'(10i8)') (ip1(i),i=1,natoms+1)
      ie = ip1(natoms+1)-1
      write(idmx,'(10i8)') (ipair(i),i=1,ie)
      write(idmx,'(10i8)') (iimat(i),i=1,natoms)
      ipmax = ip1(natoms+1)
      write(idmx,'(10i8)') (ijmat(i),i=1,ipmax)

!     density matrix
      do i=1,natoms
        iai = iatnum(i)
        norbsi = natorb(iai)
        if(norbsi.ne.0) then
          ij = iimat(i)
          do ii=1,norbsi
            write(idmx,'(4F25.20)') (pdiag(k),k=ij,ij+ii-1)
            ij = ij + ii
          enddo
        endif
      enddo

      write(idmx,'(4F25.20)') (pdiat(i),i=1,ijmat(ipmax)-1)
      close (idmx)

      end
