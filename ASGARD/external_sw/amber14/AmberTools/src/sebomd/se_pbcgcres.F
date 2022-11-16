      subroutine se_pbcgcres(i,j,xj,yj,zj)

!     returns minimum image coordinates for atom j based on its
!     location relative to atom i (atomwise)

!     note that in case of push, the values of gcres correspond to the
!     original (non-pushed) configuration. this way, the correct geometry
!     will be taken (since else dxij would always be larger than dhalf)

      implicit none
#include "sebomd.dim"
#include "se_pbcbox.h"

      integer i, j
      double precision xj, yj, zj
      double precision dxij, dyij, dzij

      xj = gcres(1,j)
      yj = gcres(2,j)
      zj = gcres(3,j)

!     center-of-residue coordinate differences:

      dxij = gcres(1,i) - gcres(1,j)
      dyij = gcres(2,i) - gcres(2,j)
      dzij = gcres(3,i) - gcres(3,j)

!     use minimum image for j if the residue centers are separated
!     by more than half the box dimension.

      if(abs(dxij).gt.dhalf(1)) xj = xj + sign(dbox(1),dxij)
      if(abs(dyij).gt.dhalf(2)) yj = yj + sign(dbox(2),dyij)
      if(abs(dzij).gt.dhalf(3)) zj = zj + sign(dbox(3),dzij)

      return
      end

