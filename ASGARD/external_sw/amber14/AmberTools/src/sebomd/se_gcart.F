      subroutine se_gcart(gxyz)
!
!     subroutine to compute cartesian energy gradient using variational
!     finite-difference derivatives.  returns gradient in kcal/angstrom.
!     see subroutine fock for more extensive comments on energy terms.
!
      implicit none
#include "sebomd.dim"
#include "se_bpairs.h"
#include "se_param.h"
#include "se_constants.h"
#include "se_parproc.h"
#include "se_logic.h"
#include "se_pmeqm.h"
#include "se_xyzcrd.h"
#include "se_esave.h"
#include "se_bal_mpi.h"
#include "se_fulsys.h"
#include "se_subsys.h"
#include "se_glbmat1.h"
#include "se_glbmat2.h"
#include "se_elemnt.h"
#ifdef MPI 
#include "mpif.h"
      integer ijk, ier
#endif
      double precision gxyz(3,maxatm)
!
!     local variables:
!
      integer i, j, ijfull, ijbond, ig, jg, ijphi
      integer ii, ij, jj, iii, jjj, ip, jp
      integer ilevel, ijp
      integer norbsi, nrepi, npqi
      integer norbsj, nrepj, npqj
      integer igaus, icart
      integer jjg, ir, iai
      integer ijg, jr, iaj
      double precision xi, yi, zi, znuci
      double precision xj, yj, zj, znucj
      double precision xisj, xipj, xidj
      double precision ai0, ai1, ai2
      double precision aj0, aj1, aj2
      double precision di1, di2
      double precision dj1, dj2
      double precision dii, dij, djj, dgij
      double precision dennij, dxennij, dyennij, dzennij
      double precision xisi, xipi, xidi
      double precision bsi, bpi, bdi
      double precision bsj, bpj, bdj
      double precision aci
      double precision acj, qiqj
      double precision delect, dgrad1, ddij, fcoul
      double precision dxesrqm, dyesrqm, dzesrqm
      double precision dxepmeqm, dyepmeqm, dzepmeqm
      double precision rij, ennij
      logical orbsi,orbsj,orbsij,share1,share2
      double precision agi(4),bgi(4),cgi(4),agj(4),bgj(4),cgj(4)
      double precision hij(9,9),
     .          repij(10,10),eenuc(2,10),dhii(9,9),dhjj(9,9),dhij(9,9),
     .          dfii(9,9),dfjj(9,9),dfij(9,9),drepij(10,10),pii(9,9),
     .          pjj(9,9),pij(9,9)
      double precision dxrepij(10,10),dyrepij(10,10),dzrepij(10,10)
      double precision dxhij(9,9),dyhij(9,9),dzhij(9,9)
      integer irepul(4,4)
      data irepul /1,2,3,4,
     .             2,5,6,7,
     .             3,6,8,9,
     .             4,7,9,10/
      save irepul

      call se_debugmsg0("in gcart")
 
      do i=1,natoms
        do j=1,3
          gxyz(j,i) = 0.0d0
#ifdef MPI 
          gtmp_mpi(j,i) = 0.0d0
#endif
        enddo
      enddo

      do i=1,4
        vir(i)=0.0d0
#ifdef MPI 
        virtmp_mpi(i) = 0.0d0
#endif
      enddo


      if(natoms.eq.1) return
!
!     loop over all possible pairs of atoms.
!
      ijfull = 1
      ijbond = 1
!
!     ijfull is a pointer to position in list of all possible atom pairs.
!     ijbond is a pointer to postiion in bonded atom pairlist ipair.
!
#ifdef MPI 
      do i=myid+2,natoms,nproc
        ijbond = ip1(i)
#else
      do i=2,natoms
#endif
        iai = iatnum(i)
        if(iai.ne.0) then
!
!         unperturbed atom i coordinates:
!
          xi = xyz(1,i)
          yi = xyz(2,i)
          zi = xyz(3,i)

          norbsi = natorb(iai)
          orbsi = norbsi.gt.0
          nrepi = max((norbsi*(norbsi+1))/2,1)
          npqi = nquant(iai)

          xisi = expnt(0,iai)
          xipi = expnt(1,iai)
          xidi = expnt(2,iai)
          ai0 = al(0,iai)
          ai1 = al(1,iai)
          ai2 = al(2,iai)
          di1 = dl(1,iai)
          di2 = dl(2,iai)
          bsi = beta(0,iai)
          bpi = beta(1,iai)
          bdi = beta(2,iai)
          aci = acore(iai)
          do igaus=1,4
            agi(igaus) = agaus(igaus,iai)
            bgi(igaus) = bgaus(igaus,iai)
            cgi(igaus) = cgaus(igaus,iai)
          enddo
          znuci = zchg(iai)
!
!         temporarily store diagonal block of density matrix for atom i.
!
          if(orbsi)then
            ij = iimat(i)
            do ii=1,norbsi
              do jj=1,ii
                pii(ii,jj) = pdiag(ij)
                pii(jj,ii) = pdiag(ij)
                ij = ij + 1
              enddo
            enddo
          endif
!
!         get intersection with all other atoms j < i.
!
          do j=1,i-1
            iaj = iatnum(j)
            if(iaj.ne.0) then
              norbsj = natorb(iaj)
              orbsj = norbsj.gt.0
              orbsij = orbsi.and.orbsj
              if(ipair(ijbond).eq.j)then
                share1 = nshare(1,ijbond).ne.0
                share2 = .true.
              else
                share1 = .false.
                share2 = .false.
              endif
              nrepj = max((norbsj*(norbsj+1))/2,1)
              npqj = nquant(iaj)
              xisj = expnt(0,iaj)
              xipj = expnt(1,iaj)
              xidj = expnt(2,iaj)
              aj0 = al(0,iaj)
              aj1 = al(1,iaj)
              aj2 = al(2,iaj)
              dj1 = dl(1,iaj)
              dj2 = dl(2,iaj)
              bsj = beta(0,iaj)
              bpj = beta(1,iaj)
              bdj = beta(2,iaj)
              acj = acore(iaj)
              do igaus=1,4
                agj(igaus) = agaus(igaus,iaj)
                bgj(igaus) = bgaus(igaus,iaj)
                cgj(igaus) = cgaus(igaus,iaj)
              enddo
              if(pbc)then
                call se_pbcxyz(i,j,xj,yj,zj)
              else
                xj = xyz(1,j)
                yj = xyz(2,j)
                zj = xyz(3,j)
              endif
              znucj = zchg(iaj)
!            
!             temporarily store diagonal block of density matrix for atom j.
!            
              if(orbsj)then
                ij = iimat(j)
                do ii=1,norbsj
                  do jj=1,ii
                    pjj(ii,jj) = pdiag(ij)
                    pjj(jj,ii) = pdiag(ij)
                    ij = ij + 1
                  enddo
                enddo
              endif
!
!            temporarily store i-j off-diagonal block of density matrix.
!
              if(orbsij.and.share1)then
                ij = ijmat(ijbond)
                do ii=1,norbsi
                  do jj=1,norbsj
                    pij(ii,jj) = pdiat(ij)
                    ij = ij + 1
                  enddo
                enddo
              endif
!          
!             compute core-core, two-electron, and possibly diatomic
!             one-electron integrals for i-j intersection with perturbed
!             coordinate for atom i.  core-electron may be skipped
!             because they are reformed from two-electron integrals.
!          
              rij = dsqrt((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)
              ilevel = 2
!          
!             skip hij if pij is zero.
!          
              if(.not.share1) ilevel = 3
!          
              call se_gdiat(iai,norbsi,npqi,xisi,xipi,xidi,ai0,ai1,ai2,
     .               di1,di2,bsi,bpi,bdi,aci,agi,bgi,cgi,xi,yi,zi,znuci,
     .               iaj,norbsj,npqj,xisj,xipj,xidj,aj0,aj1,aj2,dj1,
     .               dj2,bsj,bpj,bdj,acj,agj,bgj,cgj,xj,yj,zj,znucj,
     .               rij,ilevel,ennij,hij,repij,dxrepij,dyrepij,dzrepij,
     .                                     dxhij,dyhij,dzhij,
     .                                     dxennij,dyennij,dzennij,
     .               eenuc)
              do icart=1,3
                if (icart.eq.1) then
                  do iii=1,10
                    do jjj=1,10
                      drepij(jjj,iii)= dxrepij(jjj,iii)
                    enddo
                  enddo
                  do iii=1,9
                    do jjj=1,9
                      dhij(jjj,iii)= dxhij(jjj,iii)
                    enddo
                  enddo
                  dennij= dxennij/bohr2ang
                elseif (icart.eq.2) then
                  do iii=1,10
                    do jjj=1,10
                      drepij(jjj,iii)= dyrepij(jjj,iii)
                    enddo
                  enddo
                  do iii=1,9
                    do jjj=1,9
                      dhij(jjj,iii)= dyhij(jjj,iii)
                    enddo
                  enddo
                  dennij=dyennij/bohr2ang
                elseif (icart.eq.3) then
                  do iii=1,10
                    do jjj=1,10
                      drepij(jjj,iii)= dzrepij(jjj,iii)
                    enddo
                  enddo
                  do iii=1,9
                    do jjj=1,9
                      dhij(jjj,iii)= dzhij(jjj,iii)
                    enddo
                  enddo
                  dennij= dzennij/bohr2ang
                endif
                if(orbsij.and.share1)then
                  do ii=1,norbsi
                    do jj=1,norbsj
!                     dhij(ii,jj) = (hij(ii,jj) - dhij(ii,jj))*0.5d0
                      dfij(ii,jj) = dhij(ii,jj)
                    enddo
                  enddo
                endif
!
!           since electron-nuclear attractions are not stored in any
!           sort of global array, the changes in the diagonal blocks
!           of the 1-electron matrix must be computed directly from
!           the changes in the electron repulsions.  copy changes to
!           fock matrix.
!
                if(orbsi)then
                  do ii=1,norbsi
                    do jj=1,ii
                      dhii(ii,jj) = -znucj*drepij(irepul(ii,jj),1)
                      dfii(ii,jj) = dhii(ii,jj)
                    enddo
                  enddo
                endif

                if(orbsj)then
                  do ii=1,norbsj
                    do jj=1,ii
                      dhjj(ii,jj) = -znuci*drepij(1,irepul(ii,jj))
                      dfjj(ii,jj) = dhjj(ii,jj)
                    enddo
                  enddo
                endif
!
!           add on 2-center 2-electron terms to the changes in the
!           off-diagonal block of the fock matrix
!
                if(orbsij.and.share1)then
                  do ig=1,norbsi
                    do jg=1,norbsj
                      dgij = 0.0d0
                      do ip=1,norbsi
                        ir = irepul(ig,ip)
                        do jp=1,norbsj
                          jr = irepul(jg,jp)
                          dgij = dgij + pij(ip,jp)*drepij(ir,jr)
                        enddo
                      enddo
                      dfij(ig,jg) = dfij(ig,jg) - 0.5d0*dgij
                    enddo
                  enddo
                endif
!
!           add on 2-center 2-electron terms to changes in diagonal
!           blocks of fock matrix.
!
                if(orbsij)then
                  do ig=1,norbsi
                    do jg=1,ig
                      jjg = irepul(jg,jg)
                      ijg = irepul(ig,jg)
                      dgij = 0.0d0
                      do ip=1,norbsj
                        do jp=1,norbsj
                          ijp = irepul(ip,jp)
                          dgij = dgij + drepij(ijg,ijp)*pjj(ip,jp)
                        enddo
                      enddo
                      dfii(ig,jg) = dfii(ig,jg) + dgij
                    enddo
                  enddo
                  do ig=1,norbsj
                    do jg=1,ig
                      jjg = irepul(jg,jg)
                      ijg = irepul(ig,jg)
                      dgij = 0.0d0
                      do ip=1,norbsi
                        do jp=1,norbsi
                          ijp = irepul(ip,jp)
                          dgij = dgij + drepij(ijp,ijg)*pii(ip,jp)
                        enddo
                      enddo
                      dfjj(ig,jg) = dfjj(ig,jg) + dgij
                    enddo
                  enddo
                endif
!
!           we now have all changes in the 1-electron and fock matrices
!           due to the i-j interaction.  now compute the corresponding
!           change in the total energy.
!
                dij = 0.0d0
                dii = 0.0d0
                djj = 0.0d0
                if(orbsi)then
                  do ii=1,norbsi
                    do jj=1,ii
                      ddij = (dhii(ii,jj) + dfii(ii,jj))*pii(ii,jj)
                      dij = dij + ddij
                    enddo
                    dii = dii + ddij
                  enddo
                endif
                if(orbsj)then
                  do ii=1,norbsj
                    do jj=1,ii
                      ddij = (dhjj(ii,jj) + dfjj(ii,jj))*pjj(ii,jj)
                      dij = dij + ddij
                    enddo
                    djj = djj + ddij
                  enddo
                endif
                delect = dij - 0.5d0*(dii + djj)
                if(orbsij.and.share1)then
                  do ii=1,norbsi
                    do jj=1,norbsj
                      delect = delect + (dhij(ii,jj)
     &                                +  dfij(ii,jj))*pij(ii,jj)
                    enddo
                  enddo
                endif
!
!               update gradient entries for atoms i and j:
!
                dgrad1 = delect + dennij
     
                gxyz(icart,i) = gxyz(icart,i) + dgrad1
                gxyz(icart,j) = gxyz(icart,j) - dgrad1
                call se_pbcxyz(i,j,xj,yj,zj)
                if (pme) then
                  qiqj=chgpme(i)*chgpme(j)*elec2eV!14.39964521!*332.0806806d0
                  if (icart.eq.1) then
                    fcoul=qiqj*(xyz(1,i)-xj)/rij**3
                    vir(1) = vir(1) + ((dgrad1+fcoul)*ev2kcal)
     &                   *(xyz(1,i)-xj)
                  elseif (icart.eq.2) then
                    fcoul=qiqj*(xyz(2,i)-yj)/rij**3
                    vir(2) = vir(2) + ((dgrad1+fcoul)*ev2kcal)
     &                   *(xyz(2,i)-yj)
                  elseif (icart.eq.3) then
                    fcoul=qiqj*(xyz(3,i)-zj)/rij**3
                    vir(3) = vir(3) + ((dgrad1+fcoul)*ev2kcal)
     &                   *(xyz(3,i)-zj)
                  endif
                else
                  if (icart.eq.1) then
                    vir(1)=vir(1)+dgrad1*(xyz(1,i)-xj)*ev2kcal
                  elseif (icart.eq.2) then
                    vir(2)=vir(2)+dgrad1*(xyz(2,i)-yj)*ev2kcal
                  elseif (icart.eq.3) then
                    vir(3)=vir(3)+dgrad1*(xyz(3,i)-zj)*ev2kcal
                  endif
                endif
              enddo
!
!         update pointers for full and bond-based pairlist.
!
              ijfull = ijfull + 1
              if(share2) ijbond = ijbond + 1
            endif
          enddo
        endif
      enddo
#ifdef MPI 
      call mpi_allreduce(gxyz,gtmp_mpi,natoms*3,
     +       mpi_double_precision,mpi_sum,commsebomd,ier)
      do ijk = 1,natoms
        do i = 1, 3
          gxyz(i,ijk) = gtmp_mpi(i,ijk)
        end do
      end do
      call mpi_allreduce(vir,virtmp_mpi,4,
     +       mpi_double_precision,mpi_sum,commsebomd,ier)
      do i = 1, 4
        vir(i)=virtmp_mpi(i)
      enddo
#endif
      if (pmeqm) then
        ijphi=0
        do i=1,natoms
          xi = xyz(1,i)
          yi = xyz(2,i)
          zi = xyz(3,i)
          iai = iatnum(i)
          do  j=1,natoms
            ijphi=ijphi+1
            if (chewald) then
              dxepmeqm =  dxphipme(ijphi)*chgpme(j)*chgpme(i)
              dyepmeqm =  dyphipme(ijphi)*chgpme(j)*chgpme(i)
              dzepmeqm =  dzphipme(ijphi)*chgpme(j)*chgpme(i)
            
              dxesrqm =  dxphisr(ijphi)*chgpme(j)*chgpme(i)
              dyesrqm =  dyphisr(ijphi)*chgpme(j)*chgpme(i)
              dzesrqm =  dzphisr(ijphi)*chgpme(j)*chgpme(i)
            elseif (mullewald) then
              dxepmeqm =  dxphipme(ijphi)*atchg(j)*atchg(i)
              dyepmeqm =  dyphipme(ijphi)*atchg(j)*atchg(i)
              dzepmeqm =  dzphipme(ijphi)*atchg(j)*atchg(i)
             
              dxesrqm =  dxphisr(ijphi)*atchg(j)*atchg(i)
              dyesrqm =  dyphisr(ijphi)*atchg(j)*atchg(i)
              dzesrqm =  dzphisr(ijphi)*atchg(j)*atchg(i)
            elseif (cmewald) then
              dxepmeqm =  dxphipme(ijphi)*atchg2(j)*atchg2(i)
              dyepmeqm =  dyphipme(ijphi)*atchg2(j)*atchg2(i)
              dzepmeqm =  dzphipme(ijphi)*atchg2(j)*atchg2(i)
             
              dxesrqm =  dxphisr(ijphi)*atchg2(j)*atchg2(i)
              dyesrqm =  dyphisr(ijphi)*atchg2(j)*atchg2(i)
              dzesrqm =  dzphisr(ijphi)*atchg2(j)*atchg2(i)
            else 
              write(iout,*)'definition of the set of charges (gcartana)'
              stop
            endif
            gxyz(1,i) = gxyz(1,i) + 0.5d0*dxepmeqm 
            gxyz(2,i) = gxyz(2,i) + 0.5d0*dyepmeqm  
            gxyz(3,i) = gxyz(3,i) + 0.5d0*dzepmeqm

            gxyz(1,j) = gxyz(1,j) - 0.5d0*dxepmeqm  
            gxyz(2,j) = gxyz(2,j) - 0.5d0*dyepmeqm  
            gxyz(3,j) = gxyz(3,j) - 0.5d0*dzepmeqm

            call se_pbcxyz(i,j,xj,yj,zj)

            vir(1)=vir(1)+0.5d0*dxesrqm*(xyz(1,i)-xj)*ev2kcal
            vir(2)=vir(2)+0.5d0*dyesrqm*(xyz(2,i)-yj)*ev2kcal
            vir(3)=vir(3)+0.5d0*dzesrqm*(xyz(3,i)-zj)*ev2kcal 
          enddo
        enddo
      endif
!
!     convert energy differences to derivatives in kcal/angstrom.
!
      do i=1,natoms
        do j=1,3
          gxyz(j,i) = gxyz(j,i)*ev2kcal
        enddo
      enddo

      call se_debugmsg0("out gcart")
      return
      end
