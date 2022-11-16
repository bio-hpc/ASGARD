      integer maxid,maxtyp,maxlev
      parameter (maxid=1000,maxtyp=1000,maxlev=10)
c
      logical rprmok,wprmok,rcrdok,wcrdok,rcvdok,rcvbok,rcbdok,wnmrok
      logical save,match,match2,inlist,lespert,pimd,atm1st
c itimass = flag for TI_MASS section in the topology file
      logical find,allmas,bigmas,lreal(maxnatom),
     .        treal(maxnatom),itimass
          integer dscale,ascale,bscale, npack
      logical nomodv
c
      character name*10,title*80,crdline*80 
c
      real*8 fac,lesfac(maxtyp*maxtyp), newm
      real*8 x(maxnatom),y(maxnatom),z(maxnatom)
      real*8 tx(maxnatom),ty(maxnatom),tz(maxnatom),omass(maxnatom)
      real*8 vx(maxnatom),vy(maxnatom),vz(maxnatom),crdbox(3)
      real*8 tvx(maxnatom),tvy(maxnatom),tvz(maxnatom)
      real*8 box1,box2,getd,ochrg(maxnatom),ocgper(maxnatom)
c
      integer stdi,stderr,urprm,uwprm,urcrd,uwcrd,uw2prm,uw3prm,uwnmr
      integer i,j,k,idx,nlestyp,numlev(maxtyp),lestyp(maxnatom)
      integer lesid(maxnatom,maxlev),typid(maxtyp,maxlev)
      integer nlev(maxnatom),curlesid,ncopies(maxid),numcop
      integer bfach(maxbnd),bfac(maxbnd),afach(maxang),afac(maxang)
      integer bfacp(maxbnd),afacp(maxang),tbfacp(maxbnd),tafacp(maxang)
      integer tfacp(maxdih),ttfacp(maxdih),iacpfac(maxnatom)
      integer tiacpfac(maxnatom),imass(maxnatom),ipert(maxnatom)
      integer tfach(maxdih),tfac(maxdih),ncop
      integer tbon1(maxbndt),tbon2(maxbndt),tang1(maxangt),
     .        tang2(maxangt)
      integer tphi1(maxdiht),tphi2(maxdiht),nrem(maxnatom)
      integer revpoint(maxnatom),ipick(maxnatom),poiatom(maxnatom)
      integer tlesid(maxnatom,maxlev),tnlev(maxnatom),iacfac(maxnatom)
      integer tbfach(maxbnd),tbfac(maxbnd),tafach(maxang),tafac(maxang)
      integer ttfach(maxdih),ttfac(maxdih),newiac(maxnatom)
      integer tiacfac(maxnatom),nexcl,tnexcl,newiacfac(maxnatom)
      integer j1,k2,i2,i3,k5,i4,namel,npick,realtors,scaltyp
      integer geti,of,origpt(maxnatom),torigpt(maxnatom)
      integer totcop(maxnatom),numbb,numba,numbd,dumdih,dumiac
c	
      common /misc1i/ stdi,stderr,iacfac,tiacfac,origpt,torigpt
      common /misc2i/ lesid,nlev,bfach,bfac,afach,afac,tfach,tfac
      common /misc3i/ nlestyp,lestyp,ncopies,totcop,numcop,nrem
      common /misc4i/ bfacp,tbfacp,afacp,tafacp,tfacp,ttfacp,iacpfac
      common /misc5i/ tiacpfac,numbb,numba,numbd,dumdih,dumiac
      common /misc6i/ newiac,newiacfac,imass,realtors,scaltyp,ipert
      common /misc7i/ dscale,ascale,bscale,atm1st
c
      common /misc1r/ x,y,z,vx,vy,vz,lesfac,ochrg,ocgper,omass,crdbox
      common /misc2r/ box1,box2,newm
c
      common /misc1l/ rcrdok,allmas,bigmas,lespert,lreal,treal,itimass
      common /misc2l/ rcvdok,rcvbok, rcbdok, nomodv
c
      common /misc1c/ crdline, pimd, npack
