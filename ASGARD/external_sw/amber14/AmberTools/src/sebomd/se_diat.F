      subroutine se_diat(iai,norbsi,npqi,xisi,xipi,xidi,ai0,ai1,ai2,di1,
     .                di2,bsi,bpi,bdi,aci,agi,bgi,cgi,xi,yi,zi,znuci,
     .                iaj,norbsj,npqj,xisj,xipj,xidj,aj0,aj1,aj2,dj1,
     .                dj2,bsj,bpj,bdj,acj,agj,bgj,cgj,xj,yj,zj,znucj,
     .                rij,ilevel,ennij,hij,repij,eenuc)
!
!     ilevel codes:  0 ==> do ennij, hij, repij, eenuc
!                    1 ==> do ennij,      repij, eenuc
!                    2 ==> do ennij, hij, repij
!                    3 ==> do ennij,      repij
!                    4 ==> do             repij
!
      implicit none
#include "sebomd.dim"
#include "se_logic.h"
#include "se_param7.h"
      integer iai, norbsi, npqi
      integer iaj, norbsj, npqj
      double precision xisi, xipi, xidi, ai0, ai1, ai2
      double precision xisj, xipj, xidj, aj0, aj1, aj2
      double precision di1, di2, bsi, bpi, bdi
      double precision dj1, dj2, bsj, bpj, bdj
      double precision aci, acj
      double precision xi, yi, zi, znuci
      double precision xj, yj, zj, znucj
      double precision agi(4),bgi(4),cgi(4)
      double precision agj(4),bgj(4),cgj(4)
      double precision rij, ennij
      double precision hij(9,9),repij(10,10),eenuc(2,10)
      integer ilevel
!
!     local arrays:
!
      integer i, j, k, nimax, njmax, iamax, iamin
      double precision betai, betaj, argi, argj, repss, eij
      double precision term1, term2, term3
      double precision gterm, hterm, xterm
      double precision zij, expi, expj, dnapnb
      integer igaus, jgaus, ni, nj
      logical orbsi,orbsj
      double precision sij(9,9)
      integer irepul(10)
      data irepul /1,2,5,3,6,8,4,7,9,10/
      save irepul
!
!
      do i=1,10
       do j=1,10
        repij(j,i) = 0.0d0
       enddo
      enddo
      do i=1,9
       do j=1,9
        hij(j,i)  = 0.0d0
        sij(j,i)  = 0.0d0
       enddo
      enddo
        ennij = 0.0d0

      orbsi = norbsi.gt.0
      orbsj = norbsj.gt.0
!
!     get diatomic overlaps and assign off-diagonal block of
!     1-electron matrix (non-sparkle atoms only).
!
      if(orbsi.and.orbsj.and.(ilevel.eq.0.or.ilevel.eq.2))then
         call se_overlp(iai,norbsi,npqi,xisi,xipi,xidi,xi,yi,zi,
     .              iaj,norbsj,npqj,xisj,xipj,xidj,xj,yj,zj,
     .              rij,sij)
        betai = bsi
        do i=1,norbsi
          if(i.eq.2)then
            betai = bpi
          elseif(i.eq.5)then
            betai = bdi
          endif
          betaj = bsj
          do j=1,norbsj
            if(j.eq.2)then
              betaj = bpj
            elseif(j.eq.5)then
              betaj = bdj
            endif
            hij(i,j) = 0.5d0*(betai+betaj)*sij(i,j)
          enddo
        enddo
      endif
!
!     get electron-electron repulsions.  treat sparkles as monopoles.
!
      nimax = max(norbsi,1)
      njmax = max(norbsj,1)
      call se_repul(iai,nimax,ai0,ai1,ai2,di1,di2,xi,yi,zi,
     .           iaj,njmax,aj0,aj1,aj2,dj1,dj2,xj,yj,zj,
     .           rij,repij)
!
!     we are done if only repij has been requested.
!
      if(ilevel.eq.4) return
!
!     assemble electron-nuclear attraction integrals.
!
!     electrons on i attracted to core of j:
!
      if(ilevel.eq.0.or.ilevel.eq.1)then
        ni = (norbsi*(norbsi+1))/2
        if(ni.gt.0)then
          do k=1,ni
            eenuc(1,k) = -znucj*repij(irepul(k),1)
          enddo
        endif
!
        nj = (norbsj*(norbsj+1))/2
!
!       electrons on j attracted to core of i:
!
        if(nj.gt.0)then
          do k=1,nj
            eenuc(2,k) = -znuci*repij(1,irepul(k))
          enddo
        endif
      endif
!
!     core-core repulsion:
!
      iamax = max(iai,iaj)
      iamin = min(iai,iaj)
      zij = znuci*znucj
      repss = repij(1,1)
      eij = zij*repss
      argi = aci*rij
      if(argi.gt.25.0d0)then
        expi = 0.0d0
      else
        expi = exp(-argi)
      endif
      argj = acj*rij
      if(argj.gt.25.0d0)then
        expj = 0.0d0
      else
        expj = exp(-argj)
      endif
      if(iamin.eq.1.and.(iamax.eq.7.or.iamax.eq.8))then
!
!       n-h or o-h interactions:
!
        if(iai.eq.1)then
          hterm = expi
          xterm = expj
        else
          hterm = expj
          xterm = expi
        endif
        ennij = eij*(1.0d0 + rij*xterm + hterm)

      else
!
!       regular interaction:
!
        ennij = eij*(1.0d0 + expi + expj)

      endif
!
!     add gaussians for am1 or pm3 (non-sparkle atoms only).
!
      if(am1)then
        gterm = 0.0d0
        if(orbsi)then
          do igaus=1,4
            argi = bgi(igaus)*(rij-cgi(igaus))**2
            if(argi.lt.25.0d0)then
              gterm = gterm + agi(igaus)*exp(-argi)
            endif
          enddo
        endif
        if(orbsj)then
          do jgaus=1,4
            argj = bgj(jgaus)*(rij-cgj(jgaus))**2
            if(argj.lt.25.0d0)then
              gterm = gterm + agj(jgaus)*exp(-argj)
            endif
          enddo
        endif
        ennij = ennij + gterm*zij/rij
      elseif(rm1)then
        gterm = 0.0d0
        if(orbsi)then
          do igaus=1,4
            argi = bgi(igaus)*(rij-cgi(igaus))**2
            if(argi.lt.25.0d0)then
              gterm = gterm + agi(igaus)*exp(-argi)
            endif
          enddo
        endif
        if(orbsj)then
          do jgaus=1,4
            argj = bgj(jgaus)*(rij-cgj(jgaus))**2
            if(argj.lt.25.0d0)then
              gterm = gterm + agj(jgaus)*exp(-argj)
            endif
          enddo
        endif
        ennij = ennij + gterm*zij/rij
      elseif(am1d)then
        gterm = 0.0d0
        if(orbsi)then
          do igaus=1,4
            argi = bgi(igaus)*(rij-cgi(igaus))**2
            if(argi.lt.25.0d0)then
              gterm = gterm + agi(igaus)*exp(-argi)
            endif
          enddo
        endif
        if(orbsj)then
          do jgaus=1,4
            argj = bgj(jgaus)*(rij-cgj(jgaus))**2
            if(argj.lt.25.0d0)then
              gterm = gterm + agj(jgaus)*exp(-argj)
            endif
          enddo
        endif
        ennij = ennij + gterm*zij/rij
      elseif(pm3pddg)then
        gterm = 0.0d0
        if(orbsi)then
          do igaus=1,4
            argi = bgi(igaus)*(rij-cgi(igaus))**2
            if(argi.lt.25.0d0)then
              gterm = gterm + agi(igaus)*exp(-argi)
            endif
          end do
        endif
        if(orbsj)then
          do jgaus=1,4
            argj = bgj(jgaus)*(rij-cgj(jgaus))**2
            if(argj.lt.25.0d0)then
              gterm = gterm + agj(jgaus)*exp(-argj)
            endif
          end do
        endif
        ennij = ennij + gterm*zij/rij
        dnapnb = 1.0d0/(znuci+znucj)
        do igaus = 1,2
          do jgaus = 1,2
             term1 = znuci*pa7(igaus,iai) + znucj*pa7(jgaus,iaj)
             term2 = rij - da7(igaus,iai) - da7(jgaus,iaj)
             term3 = dnapnb*term1*exp(-10.0d0*term2*term2)
             ennij = ennij + term3
          end do
        end do
      elseif(pm3)then
        gterm = 0.0d0
        if(orbsi)then
          do igaus=1,2
            argi = bgi(igaus)*(rij-cgi(igaus))**2
            if(argi.lt.25.0d0)then
              gterm = gterm + agi(igaus)*exp(-argi)
            endif
          enddo
        endif
        if(orbsj)then
          do jgaus=1,2
            argj = bgj(jgaus)*(rij-cgj(jgaus))**2
            if(argj.lt.25.0d0)then
              gterm = gterm + agj(jgaus)*exp(-argj)
            endif
          enddo
        endif
        ennij = ennij + gterm*zij/rij
      endif
      return
      end
