      subroutine se_setbox

      implicit none
#include "sebomd.dim"
#include "se_elemnt.h"
#include "se_fulsys.h"
#include "se_pbcbox.h"
#include "se_resdef.h"
#include "se_xyzcrd.h"

!     local:
      integer i1, iat, i, j, i2, ires
      double precision xyztot(3)
      double precision totm, dd, d

!     in the following "local" refers to residue i
!
!     xyztot -> local vector sum of all coordinates
!     xyzmin -> global minimum
!     xyzmax -> global maximum
!     gcmin  -> global minimum coordinate geometric center residues
!     gcmax  -> global maximum coordinate geometric center residues


!     because center of mass is needed for virial calc, we
!     calculate COM instead of geom. center

      i1 = irpnt(1)
      iat = iatnum(i1)
      totm = emass(iat)
      do j=1,3
        xyztot(j) = emass(iat)*xyz(j,i1)
        xyzmin(j) = xyz(j,i1)
        xyzmax(j) = xyz(j,i1)
      enddo
      i2 = irpnt(2)-1
      do i=i1+1,i2
        iat = iatnum(i)
        totm = totm + emass(iat)
        do j=1,3
          xyztot(j) = xyztot(j) + emass(iat)*xyz(j,i)
          xyzmin(j) = min(xyzmin(j), xyz(j,i))
          xyzmax(j) = max(xyzmax(j), xyz(j,i))
        end do
      end do
!        nat = irpnt(2) - i1
      do j=1,3
        gc(j,1) = xyztot(j)/totm
        gcmin(j) = gc(j,1)
        gcmax(j) = gc(j,1)
      end do
      do i=i1,i2
        do j=1,3
          gcres(j,i) = gc(j,1)
        end do
      end do

!     do the same for all the other residues

      do ires=2,nres
        i1 = irpnt(ires)
        iat = iatnum(i1)
        totm = emass(iat)
        do j=1,3
          xyztot(j) = emass(iat)*xyz(j,i1)
          xyzmin(j) = min(xyzmin(j),xyz(j,i1))
          xyzmax(j) = max(xyzmax(j),xyz(j,i1))
        end do
        i2 = irpnt(ires+1)-1
        do i=i1+1,i2
          iat = iatnum(i)
          totm = totm + emass(iat)
          do j=1,3
            xyztot(j) = xyztot(j) + emass(iat)*xyz(j,i)
            xyzmin(j) = min(xyzmin(j),xyz(j,i))
            xyzmax(j) = max(xyzmax(j),xyz(j,i))
          enddo
        enddo
!       nat = irpnt(ires+1) - i1
        do j=1,3
          gc(j,ires) = xyztot(j)/totm
          gcmin(j) = min(gcmin(j), gc(j,ires))
          gcmax(j) = max(gcmax(j), gc(j,ires))
        enddo
        do i=i1,i2
          do j=1,3
            gcres(j,i) = gc(j,ires)
          enddo
        enddo
      enddo

!     check the PBC box, calculate boxmin, boxmax
      do j=1,3
        dxyz(j) = xyzmax(j) - xyzmin(j)
        dd = gcmax(j)-gcmin(j)
        d = 0.5*(dbox(j) - dd)
        boxmin(j) = gcmin(j) - d
        boxmax(j) = gcmax(j) + d
      end do

      end


