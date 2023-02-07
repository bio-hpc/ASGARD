      subroutine se_getpar
!
!     assigns semiempirical parameters based on which hamiltonian the
!     user has selected.
!
      implicit none
#include "sebomd.dim"
#include "se_elemnt.h"
#include "se_logic.h"
#include "se_ovraux.h"
#include "se_param.h"
#include "se_param1.h"
#include "se_param2.h"
#include "se_param3.h"
#include "se_param4.h"
#include "se_param5.h"
#include "se_param6.h"
#include "se_param7.h"

      integer iat, l, igaus, na, la
      double precision sqrhlf, power

      if(mndo)then
        do iat=1,83
          do l=0,2
            ucore(l,iat) = ucore1(l,iat)
            expnt(l,iat) = expnt1(l,iat)
            al(l,iat) = al1(l,iat)
            if(l.ne.0) dl(l,iat) = dl1(l,iat)
            beta(l,iat) = beta1(l,iat)
          enddo
          gss(iat) = gss1(iat)
          gpp(iat) = gpp1(iat)
          gsp(iat) = gsp1(iat)
          gp2(iat) = gp21(iat)
          hsp(iat) = hsp1(iat)
          acore(iat) = acore1(iat)
          eeatm(iat) = eeatm1(iat)
          hfatm(iat) = hfatm1(iat)
        enddo
      elseif(am1)then
        do iat=1,83
          do l=0,2
            ucore(l,iat) = ucore2(l,iat)
            expnt(l,iat) = expnt2(l,iat)
            al(l,iat) = al2(l,iat)
            if(l.ne.0) dl(l,iat) = dl2(l,iat)
            beta(l,iat) = beta2(l,iat)
          enddo
          gss(iat) = gss2(iat)
          gpp(iat) = gpp2(iat)
          gsp(iat) = gsp2(iat)
          gp2(iat) = gp22(iat)
          hsp(iat) = hsp2(iat)
          acore(iat) = acore2(iat)
          do igaus=1,4
            agaus(igaus,iat) = agaus2(igaus,iat)
            bgaus(igaus,iat) = bgaus2(igaus,iat)
            cgaus(igaus,iat) = cgaus2(igaus,iat)
          enddo
          eeatm(iat) = eeatm2(iat)
          hfatm(iat) = hfatm2(iat)
        enddo
      elseif(rm1)then
        do iat=1,83
          do l=0,2
            ucore(l,iat) = ucore5(l,iat)
            expnt(l,iat) = expnt5(l,iat)
            al(l,iat) = al5(l,iat)
            if(l.ne.0) dl(l,iat) = dl5(l,iat)
            beta(l,iat) = beta5(l,iat)
          enddo
          gss(iat) = gss5(iat)
          gpp(iat) = gpp5(iat)
          gsp(iat) = gsp5(iat)
          gp2(iat) = gp25(iat)
          hsp(iat) = hsp5(iat)
          acore(iat) = acore5(iat)
          do igaus=1,4
            agaus(igaus,iat) = agaus5(igaus,iat)
            bgaus(igaus,iat) = bgaus5(igaus,iat)
            cgaus(igaus,iat) = cgaus5(igaus,iat)
          enddo
          eeatm(iat) = eeatm5(iat)
          hfatm(iat) = hfatm5(iat)
        enddo
      elseif(am1d)then
        do iat=1,83
          do l=0,2
            ucore(l,iat) = ucore6(l,iat)
            expnt(l,iat) = expnt6(l,iat)
            al(l,iat) = al6(l,iat)
            if(l.ne.0) dl(l,iat) = dl6(l,iat)
            beta(l,iat) = beta6(l,iat)
          enddo
          gss(iat) = gss6(iat)
          gpp(iat) = gpp6(iat)
          gsp(iat) = gsp6(iat)
          gp2(iat) = gp26(iat)
          hsp(iat) = hsp6(iat)
          acore(iat) = acore6(iat)
          do igaus=1,4
            agaus(igaus,iat) = agaus6(igaus,iat)
            bgaus(igaus,iat) = bgaus6(igaus,iat)
            cgaus(igaus,iat) = cgaus6(igaus,iat)
          enddo
          eeatm(iat) = eeatm6(iat)
          hfatm(iat) = hfatm6(iat)
        enddo
      elseif(pm3pddg)then
        do iat=1,83
          do l=0,2
            ucore(l,iat) = ucore7(l,iat)
            expnt(l,iat) = expnt7(l,iat)
            al(l,iat) = al7(l,iat)
            if(l.ne.0) dl(l,iat) = dl7(l,iat)
            beta(l,iat) = beta7(l,iat)
          end do
          gss(iat) = gss7(iat)
          gpp(iat) = gpp7(iat)
          gsp(iat) = gsp7(iat)
          gp2(iat) = gp27(iat)
          hsp(iat) = hsp7(iat)
          acore(iat) = acore7(iat)
          do igaus=1,4
            agaus(igaus,iat) = agaus7(igaus,iat)
            bgaus(igaus,iat) = bgaus7(igaus,iat)
            cgaus(igaus,iat) = cgaus7(igaus,iat)
          end do
          eeatm(iat) = eeatm7(iat)
          hfatm(iat) = hfatm7(iat)
        end do
      elseif(pm3)then
        do iat=1,83
          do l=0,2
            ucore(l,iat) = ucore3(l,iat)
            expnt(l,iat) = expnt3(l,iat)
            al(l,iat) = al3(l,iat)
            if(l.ne.0) dl(l,iat) = dl3(l,iat)
            beta(l,iat) = beta3(l,iat)
          enddo
          gss(iat) = gss3(iat)
          gpp(iat) = gpp3(iat)
          gsp(iat) = gsp3(iat)
          gp2(iat) = gp23(iat)
          hsp(iat) = hsp3(iat)
          acore(iat) = acore3(iat)
          do igaus=1,2
            agaus(igaus,iat) = agaus3(igaus,iat)
            bgaus(igaus,iat) = bgaus3(igaus,iat)
            cgaus(igaus,iat) = cgaus3(igaus,iat)
          enddo
          eeatm(iat) = eeatm3(iat)
          hfatm(iat) = hfatm3(iat)
        enddo
      endif
!     compute auxiliary quantities for overlap integrls.
      sqrhlf = dsqrt(0.5d0)
      do iat=1,83
        na = nquant(iat)
        power = na + 0.5d0
        do la=0,min(2,na)
          exhalf(la,iat) = sqrhlf*expnt(la,iat)**power
        enddo
      enddo
      return
      end
