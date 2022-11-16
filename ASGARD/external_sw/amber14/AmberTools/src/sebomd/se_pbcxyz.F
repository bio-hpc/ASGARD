      subroutine se_pbcxyz(i,j,xj,yj,zj)

!     returns minimum image coordinates for atom j based on its
!     location relative to atom i (atomwise)
!
!     note that in case of push, the values of gcres correspond to the
!     original (non-pushed) configuration. this way, the correct geometry
!     will be taken (since else dxij would always be larger than dhalf)

      implicit none
#include "sebomd.dim"
#include "se_pbcbox.h"
#include "se_xyzcrd.h"
      integer i, j
      double precision xj, yj, zj
      double precision dxij, dyij, dzij
      double precision gcdxij, gcdyij, gcdzij
      double precision dlayer
      double precision ax, ay, az, amin

      double precision se_switch
      external se_switch

      xj = xyz(1,j)
      yj = xyz(2,j)
      zj = xyz(3,j)

!     center-of-residue coordinate differences:

      dxij = xyz(1,j) - xyz(1,i)
      dyij = xyz(2,j) - xyz(2,i)
      dzij = xyz(3,j) - xyz(3,i)

      gcdxij = gcres(1,j) - gcres(1,i)
      gcdyij = gcres(2,j) - gcres(2,i)
      gcdzij = gcres(3,j) - gcres(3,i)
!     write(0,'("gcres: ",6f12.6)') gcres(1,i),gcres(2,i),gcres(3,i),
!    . gcres(1,j),gcres(2,j),gcres(3,j)

!     use minimum image for j if the residue centers are separated
!     by more than half the box dimension.

!     if(abs(gcdxij).gt.dhalf(1)) xj = xj + sign(dbox(1),dxij)
!     if(abs(gcdyij).gt.dhalf(2)) yj = yj + sign(dbox(2),dyij)
!     if(abs(gcdzij).gt.dhalf(3)) zj = zj + sign(dbox(3),dzij)

      dxij = dxij - dbox(1)*anint(gcdxij/dbox(1))
      dyij = dyij - dbox(2)*anint(gcdyij/dbox(2))
      dzij = dzij - dbox(3)*anint(gcdzij/dbox(3))
      xj = xyz(1,i) + dxij
      yj = xyz(2,i) + dyij
      zj = xyz(3,i) + dzij

      return
!---------------------------------------------------------------------------

!     apply switching (don't work yet)
      dlayer = 2.0d0
      dxij = xj - xyz(1,i)
      dyij = yj - xyz(2,i)
      dzij = zj - xyz(3,i)

      write(0,'("pbcxyz old: ",9f12.6)') xj, yj, zj, dxij, dyij, dzij,
     . xyz(1,i), xyz(2,i), xyz(3,i)
!     ax = 0.5d0*erfc((abs(dxij)-(dhalf(1)-dlayer))*0.5d0)
!     ay = 0.5d0*erfc((abs(dyij)-(dhalf(2)-dlayer))*0.5d0)
!     az = 0.5d0*erfc((abs(dzij)-(dhalf(3)-dlayer))*0.5d0)
      ax = se_switch(dxij,dhalf(1)-0.0d0*dlayer,dlayer)
      ay = se_switch(dyij,dhalf(2)-0.0d0*dlayer,dlayer)
      az = se_switch(dzij,dhalf(3)-0.0d0*dlayer,dlayer)
      amin = min(min(ax,ay),az)
      xj = xyz(1,i) + dxij/amin
      yj = xyz(2,i) + dyij/amin
      zj = xyz(3,i) + dzij/amin
      write(0,'("pbcxyz new: ",6f12.6)') xj, yj, zj, ax, ay, az

      return
      end

