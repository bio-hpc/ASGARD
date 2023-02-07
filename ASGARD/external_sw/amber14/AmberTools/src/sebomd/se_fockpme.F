      subroutine se_fockpme(iter,ecore1)
!
!     bug fix: 7/26/99 -- s. dixon.  when sparkles appear after the first
!     atom, then the ijbond pointer can misidentify entries in the bonded
!     pairlist.  this was fixed by making sure that atoms i and j are both
!     valid candidates for the pairlist before using the ijbond pointer.
!
      implicit none
#include "sebomd.dim"
#include "se_pmeqm.h"
#include "se_glbmat1.h"
#include "se_glbmat2.h"
#include "se_parproc.h"
#include "se_param.h"
#include "se_xyzcrd.h"
#include "se_logic.h"
#include "se_bal_mpi.h"
#include "se_bpairs.h"
#include "se_fulsys.h"
#include "se_subsys.h"
#include "se_elemnt.h"
#include "se_esave.h"

      integer iter
      double precision ecore1
#ifdef MPI 
#include "mpif.h"
      integer ijk, ier
#endif
!
!     local:
!
      integer i, iai, norbsi, ii, nrepi, npqi, ig, ir, iorb, iorbi
      integer j, iaj, norbsj, jj, nrepj, npqj, jg, jr, iorbj
      integer ijfull, jjg, ijg, ijp, k, iii, kk, jjj
      double precision gssi, gppi, gspi, hspi
      double precision hgplus, ggplus, ggminu
      integer ipxpx, ipypx, ipzpx
      integer        ipypy, ipzpy
      integer               ipzpz
      double precision ppys, ppxs, ppzs
      double precision ppxpx, ppypx, ppzpx
      double precision        ppypy, ppzpy
      double precision               ppzpz
      double precision pxy, pxz, pyz, pxyz
      double precision hterm, gterm, ghalf, gij
      logical orbsi,orbsj,orbsij,iter1,share1,share2
      double precision agi(4),bgi(4),cgi(4),agj(4),bgj(4),cgj(4)
      double precision hij(9,9),
     .          repij(10,10),eenuc(2,10),pii(9,9),pjj(9,9),pij(9,9)
      integer ilevel
      integer iimax, ijmax, ipmax, ij, ijbond, ii0, ijaddr
      integer kjaddr, ijphi, npairs, kiaddr
      integer li, igaus, ip, jp
      integer iss, ipxs, ipys, ipzs
      double precision pss, di1, di2, pssgh, rij, ennij
      double precision xi, yi, zi, znuci, xisj, xipj, xidj
      double precision xj, yj, zj, znucj, xisi, xipi, xidi
      double precision dj1, dj2, gp2i
      double precision bsi, bpi, bdi, ai0, ai1, ai2, aci
      double precision bsj, bpj, bdj, aj0, aj1, aj2, acj
      double precision sumqphi(natoms), sumf, fpmediatn
      double precision dhki
      double precision dhkj

      integer irepul(4,4),lquant(9)
      data irepul /1,2,3,4,
     .             2,5,6,7,
     .             3,6,8,9,
     .             4,7,9,10/
      data lquant /0,1,1,1,2,2,2,2,2/
      save irepul,lquant

!
!
!     initialize fock matrix.
!
      iimax = iimat(natoms+1)-1
      ipmax = ip1(natoms+1)
      ijmax = ijmat(ipmax)-1
      do ii=1,iimax
        fdiag(ii) = 0.0d0
        fpmediag(ii) = 0.0d0
      enddo
      do ij=1,ijmax
        fdiat(ij) = 0.0d0
        fpmediat(ij) = 0.0d0
      enddo
!
!     initialize 1-electron matrix if this is the first scf iteration.
!
      iter1 = iter.eq.1
      if(iter1)then
        ecore1 = 0.0

! jv first init everything to zero, then go back and do diagonal of hdiag
! can leave this in serial version, but done this way for parallel version
        do ii=1,iimax
          hdiag(ii) = 0.0d0
        enddo
        do ij=1,ijmax
          hdiat(ij) = 0.0d0
        enddo

! jv  must set whole arrays to zero before doing any of this in parallel
! because after first cycle arrays will have garabage all over the place
! want everything zero, but diagonal done only by pieces
! can't have all pe's sum filled in diagonal again
#ifdef MPI 
        do i=myid+1,natoms,nproc
#else
        do i=1,natoms
#endif
          iai = iatnum(i)
          norbsi = natorb(iai)
!
!         skip if it's not a real atom.
!
          if(norbsi.eq.0) cycle
          ij = iimat(i)
          do ii=1,norbsi
            li = lquant(ii)
            do jj=1,ii
              if(jj.eq.ii) hdiag(ij) = ucore(li,iai)
              ij = ij + 1
            enddo
          enddo
        enddo
      endif
!
!     fock matrix contributions from 1-center, 2-electron terms.
!
! jv would it be possible to split this up by atoms
! jv then just do a global sum on fdiag, like densub ?
#ifdef MPI 
      do i=myid+1,natoms,nproc
#else
      do i=1,natoms
#endif
        iai = iatnum(i)
        norbsi = natorb(iai)
        if(norbsi.eq.0) cycle
        gssi = gss(iai)
        iss = iimat(i)
        pss = pdiag(iss)
        fdiag(iss) = fdiag(iss) + 0.5d0*pss*gssi
        if(norbsi.eq.1) cycle
        gppi = gpp(iai)
        gspi = gsp(iai)
        gp2i = gp2(iai)
        hspi = hsp(iai)
        hgplus = hspi + gspi
        ggplus = gppi + gp2i
        ggminu = gppi - gp2i
!
!       indexing scheme for lower triangle of diagonal block:
!
!       iimat(i)                                        iss
!       iimat(i)+1 iimat(i)+2                       --\ ipxs ipxpx
!       iimat(i)+3 iimat(i)+4 iimat(i)+5            --/ ipys ipypx ipypy
!       iimat(i)+6 iimat(i)+7 iimat(i)+8 iimat(i)+9     ipzs ipzpx ipzpy ipzpz
!
        ipxs = iss + 1
        ppxs = pdiag(ipxs)
        ipxpx = ipxs + 1
        ppxpx = pdiag(ipxpx)
        ipys = ipxpx + 1
        ppys = pdiag(ipys)
        ipypx = ipys + 1
        ppypx = pdiag(ipypx)
        ipypy = ipypx + 1
        ppypy = pdiag(ipypy)
        ipzs = ipypy + 1
        ppzs = pdiag(ipzs)
        ipzpx = ipzs + 1
        ppzpx = pdiag(ipzpx)
        ipzpy = ipzpx + 1
        ppzpy = pdiag(ipzpy)
        ipzpz = ipzpy + 1
        ppzpz = pdiag(ipzpz)
!
        pxy = pdiag(ipxpx) + pdiag(ipypy)
        pxz = pdiag(ipxpx) + pdiag(ipzpz)
        pyz = pdiag(ipypy) + pdiag(ipzpz)
        pxyz = pxy + pdiag(ipzpz)
        fdiag(iss) = fdiag(iss) + pxyz*(gspi - 0.5d0*hspi)
!
        hterm = 2.0d0*hspi - 0.5d0*hgplus
        fdiag(ipxs) = fdiag(ipxs) + ppxs*hterm
        fdiag(ipys) = fdiag(ipys) + ppys*hterm
        fdiag(ipzs) = fdiag(ipzs) + ppzs*hterm
!
        pssgh = pss*(gspi - 0.5d0*hspi)
        ghalf = gppi*0.5d0
        gterm = gp2i - 0.25d0*ggminu
        fdiag(ipxpx) = fdiag(ipxpx) + pssgh + ppxpx*ghalf + pyz*gterm
        fdiag(ipypy) = fdiag(ipypy) + pssgh + ppypy*ghalf + pxz*gterm
        fdiag(ipzpz) = fdiag(ipzpz) + pssgh + ppzpz*ghalf + pxy*gterm
!
        gterm = ggminu - 0.25d0*ggplus
        fdiag(ipypx) = fdiag(ipypx) + ppypx*gterm
        fdiag(ipzpx) = fdiag(ipzpx) + ppzpx*gterm
        fdiag(ipzpy) = fdiag(ipzpy) + ppzpy*gterm
      enddo

!
!     off-diagonal blocks and diatomic contributions to diagonal blocks:
!
      if(natoms.eq.1) return
!
      ijfull = 1
      ijbond = 1
!
!     ijfull is a pointer to position in list of all possible atom pairs.
!     ijbond is a pointer to postiion in bonded atom pairlist ipair.
!
! jv think about doing just atoms of the subsystem here - but maybe can't
! if they require interaction with all other atoms
#ifdef MPI 
      do i=myid+2,natoms,nproc
! jv make ijbond start at the beginning of the pairlist for this atom
        ijbond = ip1(i)
#else
      do i=2,natoms
#endif
        ii0 = iimat(i)
        iai = iatnum(i)
        if(iai.eq.0) cycle
        norbsi = natorb(iai)
        orbsi = norbsi.gt.0
        nrepi = max((norbsi*(norbsi+1))/2,1)

        ai0 = al(0,iai)
        ai1 = al(1,iai)
        ai2 = al(2,iai)
        di1 = dl(1,iai)
        di2 = dl(2,iai)
        bsi = beta(0,iai)
        bpi = beta(1,iai)
        bdi = beta(2,iai)
        if(iter1)then
          npqi = nquant(iai)
          xisi = expnt(0,iai)
          xipi = expnt(1,iai)
          xidi = expnt(2,iai)
          aci = acore(iai)
          do igaus=1,4
            agi(igaus) = agaus(igaus,iai)
            bgi(igaus) = bgaus(igaus,iai)
            cgi(igaus) = cgaus(igaus,iai)
          enddo
        endif

        znuci = zchg(iai)
!
!       temporarily store diagonal block of density matrix for atom i.
!
        if(orbsi)then
          ij = iimat(i)
          do ip=1,norbsi
            do jp=1,ip
              pii(ip,jp) = pdiag(ij)
              pii(jp,ip) = pdiag(ij)
              ij = ij + 1
            enddo
          enddo
        endif
        xi = xyz(1,i)
        yi = xyz(2,i)
        zi = xyz(3,i)
!
!       get i-j intersection in fock matrix.  the diatomic block fdiat
!       will only be filled if i and j are a "bonded" pair stored
!       in ipair.  however i and j always affect the diagonal block
!       fdiag of the other atom.
!
        do j=1,i-1
          iaj = iatnum(j)
          if(iaj.eq.0) cycle

          norbsj = natorb(iaj)
          orbsj = norbsj.gt.0
          orbsij = orbsi.and.orbsj
!
!         see if (i,j) is the next pair in the bonded pairlist and
!         assign "share" status.
!
! **********************************************************************
!         bug fix - an additional condition was added:
!
!              (previous)                       (fixed)
!         if(ipair(ijbond).eq.j) --> if(orbsij.and.ipair(ijbond).eq.j)
!
!         when atom i is a sparkle, then (i,j) is not a valid candidate
!         for the bonded pairlist, so the ipair(ijbond).eq.j test is
!         not valid.  the orbsij condition was added to make sure a
!         coincidental match doesn't throw off the share1/share2 status
!         and the fock indexing scheme.
! **********************************************************************
!
          if(orbsij.and.ipair(ijbond).eq.j)then
            share1 = nshare(1,ijbond).ne.0
            share2 = .true.
          else
            share1 = .false.
            share2 = .false.
          endif
!
!         share1 indicates whether or not atoms i and j ever appear in the
!         same subsystem as non-buffer atoms.  if not, then the i-j block
!         of pdiat will be zero.  this information may be used to skip over
!         the 2-electron portion of fdiat.
!
!         share2 indicates whether the atoms ever appear with any status
!         in the same subsystems.  if not, then the i-j block of hdiat
!         will be ignored and there is no need to compute overlap integrals.
!
!         note that both share1 and share2 are affected if the cutbond
!         option is used.  in that case, atoms i and j must be in the
!         same subsystem and separated by no more than the cutoff to
!         be considered.
!
          nrepj = max((norbsj*(norbsj+1))/2,1)
          aj0 = al(0,iaj)
          aj1 = al(1,iaj)
          aj2 = al(2,iaj)
          dj1 = dl(1,iaj)
          dj2 = dl(2,iaj)
          bsj = beta(0,iaj)
          bpj = beta(1,iaj)
          bdj = beta(2,iaj)
          if(iter1)then
            npqj = nquant(iaj)
            xisj = expnt(0,iaj)
            xipj = expnt(1,iaj)
            xidj = expnt(2,iaj)
            acj = acore(iaj)
            do igaus=1,4
              agj(igaus) = agaus(igaus,iaj)
              bgj(igaus) = bgaus(igaus,iaj)
              cgj(igaus) = cgaus(igaus,iaj)
            enddo
          endif
          if(pbc)then
            call se_pbcxyz(i,j,xj,yj,zj)
          else
            xj = xyz(1,j)
            yj = xyz(2,j)
            zj = xyz(3,j)
          endif
          znucj = zchg(iaj)
!
!         temporarily store diagonal block of density matrix for atom j.
!
          if(orbsj)then
            ij = iimat(j)
            do ip=1,norbsj
              do jp=1,ip
                pjj(ip,jp) = pdiag(ij)
                pjj(jp,ip) = pdiag(ij)
                ij = ij + 1
              enddo
            enddo
          endif
!
!         temporarily store off-diagonal block of density matrix for i-j
!         intersection.
!
          if(orbsij.and.share1)then
            ij = ijmat(ijbond)
            do ip=1,norbsi
              do jp=1,norbsj
                pij(ip,jp) = pdiat(ij)
                ij = ij + 1
              enddo
            enddo
          endif
!
!         integral computation
!         --------------------
!
!         first scf iteration:    compute 2-electron integrals repij,
!                                 1-electron integrals eenuc, and
!                                 core-core repulsions ennij.
!                                 compute diatomic 1-electron integrals
!                                 hij only if atoms i and j share a
!                                 subsystem (any status).
!
!
!         subsequent iterations:  if direct calculation, compute only repij.
!                                 if not direct, then compute no integrals.
!
          rij = dsqrt((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)
          ilevel = 0
          if(.not.iter1)then
!
!           direct calculation, and not first iteration.  get only
!           two-electron integrals.
!
            ilevel = 4
          elseif(.not.share2)then
!
!           first iteration, but not bonded.  need only core-electron
!           integrals eenuc and core-core repulsion ennij.
!
            ilevel = 1
          endif

          call se_diat(iai,norbsi,npqi,xisi,xipi,xidi,ai0,ai1,ai2,di1,
     .              di2,bsi,bpi,bdi,aci,agi,bgi,cgi,xi,yi,zi,znuci,
     .              iaj,norbsj,npqj,xisj,xipj,xidj,aj0,aj1,aj2,dj1,
     .              dj2,bsj,bpj,bdj,acj,agj,bgj,cgj,xj,yj,zj,znucj,
     .              rij,ilevel,ennij,hij,repij,eenuc)

!
!         on first scf iteration, store 1-electron matrix and
!         core-core repulsions.
!
          if(iter1)then
!
!           1-electron (diagonal, i):
!
            if(orbsi)then
              ij = iimat(i)
              do ii=1,nrepi
                hdiag(ij) = hdiag(ij) + eenuc(1,ii)
                ij = ij + 1
              enddo
            endif
!
!           1-electron (diagonal, j):
!
            if(orbsj)then
              ij = iimat(j)
              do jj=1,nrepj
                hdiag(ij) = hdiag(ij) + eenuc(2,jj)
                ij = ij + 1
              enddo
            endif
!
!           1-electron (off-diagonal, ij):
!
            if(orbsij.and.share2)then
              ij = ijmat(ijbond)
              do ii=1,norbsi
                do jj=1,norbsj
                  hdiat(ij) = hij(ii,jj)
                  ij = ij + 1
                enddo
              enddo
            endif
!
!           core-core repulsions:
!
            ecore1 = ecore1 + ennij
          endif
!
!         determine off-diagonal block of fock matrix for i-j intersection.
!         block will be zero, however, if atoms i and j do not share a
!         subsystem as non-buffer atoms.
!         
          if(orbsij.and.share1)then
            ij = ijmat(ijbond)
            do ig=1,norbsi
              do jg=1,norbsj
                gij = 0.0d0
                do ip=1,norbsi
                  ir = irepul(ig,ip)
                  do jp=1,norbsj
                    jr = irepul(jg,jp)
                    gij = gij + pij(ip,jp)*repij(ir,jr)
                  enddo
                enddo
                fdiat(ij) = fdiat(ij) - 0.5d0*gij
                ij = ij + 1
              enddo
            enddo
          endif
!         
!         add 2-center 2-electron terms arising from atom j to diagonal
!         block of atom i.
!         
          ij = ii0
          do ig=1,norbsi
            do jg=1,ig
              jjg = irepul(jg,jg)
              ijg = irepul(ig,jg)
              gij = 0.0d0
              do ip=1,norbsj
                do jp=1,norbsj
                  ijp = irepul(ip,jp)
                  gij = gij + repij(ijg,ijp)*pjj(ip,jp)
                enddo
              enddo
              fdiag(ij) = fdiag(ij) + gij
              ij = ij + 1
            enddo
          enddo
          
          
!         
!         add 2-center 2-electron terms arising from atom i to diagonal
!         of atom j.
!         
          ij = iimat(j)
          do ig=1,norbsj
            do jg=1,ig
              jjg = irepul(jg,jg)
              ijg = irepul(ig,jg)
              gij = 0.0d0
              do ip=1,norbsi
                do jp=1,norbsi
                  ijp = irepul(ip,jp)
                  gij = gij + repij(ijp,ijg)*pii(ip,jp)
                enddo
              enddo
              fdiag(ij) = fdiag(ij) + gij
              ij = ij + 1
            enddo
          enddo
!         
!         update pointers for full and bond-based pairlist.
!         
          
          ijfull = ijfull + 1
          if(share2) ijbond = ijbond + 1
        enddo
      enddo


#ifdef MPI 
! jv sum fdiag and fdiat globally
      call mpi_allreduce(fdiag,tmp_mpi,iimax,
     +                   mpi_double_precision,mpi_sum,commsebomd,ier)
      do ijk = 1,iimax
        fdiag(ijk) = tmp_mpi(ijk)
      enddo
      call mpi_allreduce(fdiat,tmp_mpi,ijmax,
     +                   mpi_double_precision,mpi_sum,commsebomd,ier)
      do ijk = 1,ijmax
       fdiat(ijk) = tmp_mpi(ijk)
      enddo
! jv do this only if iter1 also
      if(iter1) then
        call mpi_allreduce(hdiag,tmp_mpi,iimax,
     +                     mpi_double_precision,mpi_sum,commsebomd,ier)
        do ijk = 1,iimax
          hdiag(ijk) = tmp_mpi(ijk)
        enddo
        call mpi_allreduce(hdiat,tmp_mpi,ijmax,
     +                     mpi_double_precision,mpi_sum,commsebomd,ier)
        do ijk = 1,ijmax
          hdiat(ijk) = tmp_mpi(ijk)
        enddo
! jv sum ecore from all pe's - only on first iter (first subsys)
        call mpi_allreduce(ecore1,tmp_mpi(1),1,
     +                     mpi_double_precision,mpi_sum,commsebomd,ier)
        ecore1 = tmp_mpi(1)
      endif
#endif

      if (chewald) then
        ijphi=0
        do i=1,natoms
          iai = iatnum(i)
          norbsi = natorb(iai)
          sumf=0.0d0
          do j=1,natoms
            ijphi=ijphi+1
            sumf = sumf + phipme(ijphi)*chgpme(j)
          enddo
          ij = iimat(i)-1
          do iorb=1,norbsi
            ij=ij+iorb
            fpmediag(ij) =  - sumf
          enddo
        enddo
      elseif (mullewald) then
!       get the mulliken charges     
        call se_pme_calcchg
        ijphi=0
        do i=1,natoms
          iai = iatnum(i)
          norbsi = natorb(iai) 
          sumf=0.0d0
          do j=1,natoms
             ijphi=ijphi+1
             sumf = sumf + phipme(ijphi)*atchg(j)
          enddo
          ij = iimat(i)-1
          do iorb=1,norbsi
             ij=ij+iorb
             fpmediag(ij) =  - sumf
!            fdiag(ij) = fdiag(ij) - sumf 
          enddo
        enddo
      elseif (cmewald) then
!       get the cm1 charges      
        call se_pme_calccm1
        do i=1,natoms
          sumqphi(i)=0.0d0
          ijphi=(i-1)*natoms+1
          do j=1,natoms
            sumqphi(i)=sumqphi(i)+phipme(ijphi)*atchg2(j)
            ijphi=ijphi+1
          enddo
        enddo
        ijphi=0
        ijbond = 1
        npairs = ip1(natoms+1)-1
#ifdef MPI 
        do  i=myid+1,natoms,nproc
#else
        do  i=1,natoms
#endif
          iai = iatnum(i)
          norbsi = natorb(iai) 
          ij = iimat(i)-1
          do iorb=1,norbsi
            ij=ij+iorb
            fpmediag(ij) =  fpmediag(ij)*sumqphi(i)
!           fdiag(ij) = fdiag(ij) - sumf 
          enddo
          if (i.ge.2) then
            do j=1,i-1
              ii = max(i,j)
              jj = min(i,j)
              call se_ijfind(npairs,ii,jj,ijaddr)
              if (ijaddr.ne.0) then
                iaj = iatnum(j)
                norbsj = natorb(iaj) 
                ij = ijmat(ijaddr)
                fpmediatn = 0.0d0
                call se_pme_getdhk(iai,iaj,dhkj)
                call se_pme_getdhk(iaj,iai,dhki)
                do k=1,natoms
                  if ((k.ne.i).and.(k.ne.j)) then
                    iii=max(i,k)
                    kk=min(i,k)
                    kiaddr=(iii*(iii-3))/2 + kk + 1
                    jjj=max(j,k)
                    kk=min(j,k)
                    kjaddr=(jjj*(jjj-3))/2 + kk + 1
                    fpmediatn = fpmediatn - 
     &              sumqphi(k)*(bkk1(kiaddr)*dhkj+bkk1(kjaddr)*dhki)
                  endif
                enddo                       
                call se_ijfind(npairs,i,j,ijaddr)
                ij=ijmat(ijaddr)
                do iorbi=1,norbsi  
                  do iorbj=1,norbsj  
                    fpmediat(ij) = fpmediat1(ij)*sumqphi(i)
     &                   + fpmediat2(ij)*sumqphi(j)
     &                   + fpmediatn*pdiat(ij)
                    ij=ij+1
                  enddo
                enddo
              endif
              ijbond=ijbond+1
            enddo
          endif
        enddo
! gm (beg)
#ifdef MPI 
        call mpi_allreduce(fpmediag,tmp_mpi,iimax,
     +                     mpi_double_precision,mpi_sum,commsebomd,ier)
        do ijk = 1,iimax
          fpmediag(ijk) = tmp_mpi(ijk)
        enddo
        call mpi_allreduce(fpmediat,tmp_mpi,ijmax,
     +                     mpi_double_precision,mpi_sum,commsebomd,ier)
        do ijk = 1,ijmax
          fpmediat(ijk) = tmp_mpi(ijk)
        enddo
#endif
! gm (end)
      else 
        write(iout,*)'definition of the set of charges (fockpme)'
        stop
      endif                     !(set of charges)
!
!     add in 1-electron matrix.
!
      do ii=1,iimax
        fdiag(ii) = fdiag(ii) + hdiag(ii)
      enddo
      do ij=1,ijmax
        fdiat(ij) = fdiat(ij) + hdiat(ij)
      enddo
      return
      end
