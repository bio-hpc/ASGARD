      subroutine addspace_atm1st( ipick, curlesid )
      implicit none
c     
#include "SIZE.h"
#include "TOP.h"
#include "TTOP.h"
#include "MISC.h"
      integer :: size,ncopy,oatom,tempint,chunk,iexcl

      integer, allocatable :: copy(:)

      real*8 velfac, randnum

      logical size_error

      allocate( copy(natom) )

      size=0
      do i=1,natom
         if( ipick(i) > 1 ) then
            size=size+1
            copy(size)=i
         end if
      end do

      ncopy=0
      tnatom=natom
      do i=1,natom
         if(ipick(i) > 1 ) then
            if(ncopy==0) ncopy=ipick(i)
            if(ncopy/=ipick(i)) then
               print *, 'inconsistent copy number'
               stop
            end if
            tnatom=tnatom+1
            revpoint(i)=tnatom
            chrg(i)=chrg(i) / dble(ncopy)
            iacfac(i)=iacfac(i) * ncopy
            nlev(i)=nlev(i)+1
            if(.not.allmas) amass(i)=amass(i)/dble(ncopy)
         else
            revpoint(i)=0
         end if
      end do

      chunk=ncopy/npack
      tnatom=natom
      do i = 2, ncopy
         do j=1, size
            oatom = copy(j)
            tnatom = tnatom + 1
            poiatom(tnatom)=oatom
            if(rcrdok.or.rcvdok.or.rcvbok.or.rcbdok) then
               if( npack.gt.1 ) then
                  tempint = ( ( i - 1 ) / chunk ) * natom
                  x(tnatom)=x(tempint+oatom)
                  y(tnatom)=y(tempint+oatom)
                  z(tnatom)=z(tempint+oatom)
               else
                  x(tnatom)=x(oatom)
                  y(tnatom)=y(oatom)
                  z(tnatom)=z(oatom)
               endif
            endif
            
            if( rcvdok.or.rcvbok) then
               velfac=1.0d0
               if(.not.allmas.and.ipick(i).gt.1) then
                  velfac=velfac*sqrt(float(ipick(i)))
               end if
               if(.not.nomodv.and.ipick(i).gt.1) then
                  call amrand(randnum)
                  velfac=velfac*(randnum*0.2+0.9)
               end if
               vx(tnatom)=vx(oatom)*velfac
               vy(tnatom)=vy(oatom)*velfac
               vz(tnatom)=vz(oatom)*velfac
            end if

            nlev(tnatom)=nlev(oatom)
            lesid(tnatom,nlev(tnatom))=curlesid

            chrg(tnatom)=chrg(oatom)
            amass(tnatom)=amass(oatom)
            rborn(tnatom)=rborn(oatom)
            fs(tnatom)=fs(oatom)
            igraph(tnatom)=igraph(oatom)
            isymbl(tnatom)=isymbl(oatom)
            itree(tnatom)=itree(oatom)
            iac(tnatom)=iac(oatom)
            iacfac(tnatom)=iacfac(oatom)
         end do
      end do
c     
c     residue
c     
      do i=2,ncopy
         nres=nres+1
         ipres(nres)=natom + (i-2)*size + 1
         labres(nres)='CPY'
      end do
c     
c     
c     bond with hydrogen
c     
      tnbonh = nbonh
      do j=1,ncopy
         do i=1,nbonh
            if( ipick(ibh(i))>1.or.ipick(jbh(i))>1 ) then
               if(j==1) then
                  bfach(i)=bfach(i)*ncopy
               else
                  tnbonh = tnbonh + 1
                  if( ipick(ibh(i)) == 1 ) then
                     ibh(tnbonh) = ibh(i)
                  else
                     ibh(tnbonh) = revpoint(ibh(i)) + (j-2) * size
                  end if

                  if( ipick(jbh(i)) == 1 ) then
                     jbh(tnbonh) = jbh(i)
                  else
                     jbh(tnbonh) = revpoint(jbh(i)) + (j-2) * size
                  end if
                  icbh(tnbonh)=icbh(i)
                  bfach(tnbonh)=bfach(i)
               end if
            end if
         end do
      end do
      nbonh = tnbonh
c     
c     bond without hydrogen
c     
      tnbona = nbona
      do j=1,ncopy
         do i=1,nbona
            if( ipick(ib(i))>1.or.ipick(jb(i))>1 ) then
               if(j==1) then
                  bfac(i)=bfac(i)*ncopy
               else
                  tnbona = tnbona + 1
                  if( ipick(ib(i))==1 ) then
                     ib(tnbona) = ib(i)
                  else
                     ib(tnbona) = revpoint(ib(i)) + (j-2) * size
                  end if

                  if( ipick(jb(i))==1 ) then
                     jb(tnbona) = jb(i)
                  else
                     jb(tnbona) = revpoint(jb(i)) + (j-2) * size
                  end if
                  icb(tnbona)=icb(i)
                  bfac(tnbona)=bfac(i)
               end if
            end if
         end do
      end do
      nbona=tnbona
      mbona=tnbona
c     
c     angle with hydrogen
c     
      tntheth=ntheth
      do j=1,ncopy
         do i=1,ntheth
            if(ipick(ith(i))>1.or.ipick(jth(i))>1.or. 
     c         ipick(kth(i))>1) then
               if(j==1) then
                  afach(i)=afach(i)*ncopy
               else
                  tntheth=tntheth+1
                  if(ipick(ith(i))==1) then
                     ith(tntheth)=ith(i)
                  else
                     ith(tntheth)=revpoint(ith(i)) + (j-2)*size
                  end if

                  if(ipick(jth(i))==1) then
                     jth(tntheth)=jth(i)
                  else
                     jth(tntheth)=revpoint(jth(i)) + (j-2)*size
                  end if

                  if(ipick(kth(i))==1) then
                     kth(tntheth)=kth(i)
                  else
                     kth(tntheth)=revpoint(kth(i)) + (j-2)*size
                  end if
                  icth(tntheth)=icth(i)
                  afach(tntheth)=afach(i)
               end if
            end if
         end do
      end do
      ntheth=tntheth

      numbb=0
      tnbper=0
      numba=0
      tngper=0
      numbd=0
      tndper=0
c     
c     
c     angle without hydrogen
c     
c     
      tntheta=ntheta
      do j=1,ncopy
         do i=1,ntheta
            if(ipick(it(i))>1.or.ipick(jt(i))>1.or.
     c         ipick(kt(i))>1) then
               if(j==1) then
                  afac(i)=afac(i)*ncopy
               else
                  tntheta=tntheta+1
                  if(ipick(it(i))==1) then
                     it(tntheta)=it(i)
                  else
                     it(tntheta)=revpoint(it(i)) + (j-2)*size
                  end if

                  if(ipick(jt(i))==1) then
                     jt(tntheta)=jt(i)
                  else
                     jt(tntheta)=revpoint(jt(i)) + (j-2)*size
                  end if

                  if(ipick(kt(i))==1) then
                     kt(tntheta)=kt(i)
                  else
                     kt(tntheta)=revpoint(kt(i)) + (j-2)*size
                  end if
                  ict(tntheta)=ict(i)
                  afac(tntheta)=afac(i)
               end if
            end if
         end do
      end do

      ntheta = tntheta
      mtheta = tntheta
c     
c     torsion with hydrogen
c     
      tnphih=nphih

      do j=1,ncopy
         do i=1,nphih
            if( ipick(iph(i))>1 .or. ipick(jph(i))>1 .or.
     &          ipick(kph(i))>1 .or. ipick(lph(i))>1 ) then
               if(j==1) then
                  tfach(i)=tfach(i)*ncopy
               else
                  tnphih=tnphih+1
                  if(ipick(iph(i))==1) then
                     iph(tnphih)=iph(i)
                  else
                     iph(tnphih)=revpoint(iph(i)) + (j-2)*size
                  end if

                  if(ipick(jph(i))==1) then
                     jph(tnphih)=jph(i)
                  else
                     jph(tnphih)=revpoint(jph(i)) + (j-2)*size
                  end if

                  if(ipick(kph(i))==1) then
                     kph(tnphih)=kph(i)
                  else
                     kph(tnphih)=revpoint(kph(i)) + (j-2)*size
                  end if
               
                  if(ipick(lph(i))==1) then
                     lph(tnphih)=lph(i)
                  else
                     lph(tnphih)=revpoint(lph(i)) + (j-2)*size
                  end if

                  inegh(tnphih)=inegh(i)
                  jnegh(tnphih)=jnegh(i)
                  knegh(tnphih)=knegh(i)
                  lnegh(tnphih)=lnegh(i)

                  if(ipick(iph(i)).eq.1.and.ipick(lph(i)).eq.1) then
                     knegh(tnphih)=-1
                  end if

                  icph(tnphih)=icph(i)
                  tfach(tnphih)=tfach(i)
               end if
            end if
         end do
      end do
      nphih=tnphih
c     
c     
c     torsion without hydrogen
c     
      tnphia=nphia
      do j=1,ncopy
         do i=1,nphia
            if( ipick(ip(i))>1 .or. ipick(jp(i))>1 .or.
     &          ipick(kp(i))>1 .or. ipick(lp(i))>1 ) then
               if(j==1) then
                  tfac(i)=tfac(i)*ncopy
               else
                  tnphia=tnphia+1
                  if(ipick(ip(i))==1) then
                     ip(tnphia)=ip(i)
                  else
                     ip(tnphia)=revpoint(ip(i)) + (j-2)*size
                  end if

                  if(ipick(jp(i))==1) then
                     jp(tnphia)=jp(i)
                  else
                     jp(tnphia)=revpoint(jp(i)) + (j-2)*size
                  end if

                  if(ipick(kp(i))==1) then
                     kp(tnphia)=kp(i)
                  else
                     kp(tnphia)=revpoint(kp(i)) + (j-2)*size
                  end if
               
                  if(ipick(lp(i))==1) then
                     lp(tnphia)=lp(i)
                  else
                     lp(tnphia)=revpoint(lp(i)) + (j-2)*size
                  end if

                  ineg(tnphia)=ineg(i)
                  jneg(tnphia)=jneg(i)
                  kneg(tnphia)=kneg(i)
                  lneg(tnphia)=lneg(i)

                  if( ipick(ip(i)).eq.1.and.ipick(lp(i)).eq.1) then
                     knegh(tnphia)=-1
                  end if
                  icp(tnphia)=icp(i)
                  tfac(tnphia)=tfac(i)
               end if
            end if
         end do
      end do
      nphia=tnphia
      mphia=tnphia
c     
c     
c     nonbond exclusion list
c     
c     
      tnexcl = next

      iexcl=0
      do i=1,natom
         tnatom = revpoint(i)
         if(ipick(i)>1) then
            numex(tnatom)=0
         end if

         do k=iexcl + 1, iexcl + numex(i)
            if(natex(k)==0) then
               numex(tnatom)=numex(tnatom)+1
               tnexcl=tnexcl+1
               natex(tnexcl)=0
            else if(ipick(i)>1.and.ipick(natex(k))>1) then
               numex(tnatom)=numex(tnatom)+1
               tnexcl = tnexcl + 1
               natex(tnexcl)= revpoint(natex(k))
            end if
         end do
         iexcl = iexcl+numex(i)
      end do

      do j=3,ncopy
         tempint=natom+(j-2)*size
         do i=1,size
            numex(tempint+i)=numex(natom+i)
         end do

         tempint=next+(j-2)*(tnexcl-next)
         do i=1,(tnexcl-next)
            if(natex(next+i)==0) then
               natex(tempint+i)=0
            else
               natex(tempint+i)=natex(next+i)+(j-2)*size
            end if
         end do
      end do
      
      natom = natom + size * (ncopy-1)
      next = next + (ncopy-1)*(tnexcl-next)

      deallocate( copy )
      call checksz( size_error )
      if (size_error) then
         stop
      endif


      end subroutine
