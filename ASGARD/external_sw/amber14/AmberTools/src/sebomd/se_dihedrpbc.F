C------------------------------------------------------------------------
!     Modified by Antoine MARION :: 2013-11-08
!          --  add PBC
C------------------------------------------------------------------------
      subroutine se_dihedrpbc(pbc,xyz,i,ia,ib,ic,dih)
c
c     determines the i-ia-ib-ic dihedral angle in radians.  the
c     angle dih is positive if ic is located clockwise from i when
c     viewing from ia through ib.
c
      implicit none
      logical pbc
      double precision xyz(3,*)
      integer i, ia, ib, ic
      double precision dih

      double precision ix,iy,iz
      double precision iax,iay,iaz
      double precision ibx,iby,ibz
      double precision icx,icy,icz

      double precision aiix, aiiy, aiiz
      double precision bcx, bcy, bcz
      double precision abx, aby, abz
      double precision dot1, dot2, dot3, dot4
      double precision proj1, proj2, absqr
      double precision aibcx, aibcy, aibcz
      double precision direc
      double precision aileng, bcleng, cosine

      ix = xyz(1,i)
      iy = xyz(2,i)
      iz = xyz(3,i)


      if (pbc) then
       call se_pbcxyz(i,ia,iax,iay,iaz)
       call se_pbcxyz(i,ib,ibx,iby,ibz)
       call se_pbcxyz(i,ic,icx,icy,icz)
      else
       iax = xyz(1,ia)
       iay = xyz(2,ia)
       iaz = xyz(3,ia)

       ibx = xyz(1,ib)
       iby = xyz(2,ib)
       ibz = xyz(3,ib)

       icx = xyz(1,ic)
       icy = xyz(2,ic)
       icz = xyz(3,ic)
      endif

c
c     shift ia-i and ib-ic bond vectors to a common origin.
c

      aiix = ix - iax
      aiiy = iy - iay
      aiiz = iz - iaz
      bcx  = icx - ibx
      bcy  = icy - iby
      bcz  = icz - ibz

!      aiix = xyz(1,i) - xyz(1,ia)
!      aiiy = xyz(2,i) - xyz(2,ia)
!      aiiz = xyz(3,i) - xyz(3,ia)
!      bcx = xyz(1,ic) - xyz(1,ib)
!      bcy = xyz(2,ic) - xyz(2,ib)
!      bcz = xyz(3,ic) - xyz(3,ib)
c
c     form the ia-ib bond axis vector.
c

      abx = ibx - iax
      aby = iby - iay
      abz = ibz - iaz

!      abx = xyz(1,ib) - xyz(1,ia)
!      aby = xyz(2,ib) - xyz(2,ia)
!      abz = xyz(3,ib) - xyz(3,ia)
c
c     remove from (aiix,aiiy,aiiz) and (bcx,bcy,bcz) any projection along
c     the (abx,aby,abz) axis.
c
      dot1 = aiix*abx + aiiy*aby + aiiz*abz
      absqr = abx**2 + aby**2 + abz**2
      proj1 = dot1/absqr
      aiix = aiix - proj1*abx
      aiiy = aiiy - proj1*aby
      aiiz = aiiz - proj1*abz
      dot2 = bcx*abx + bcy*aby + bcz*abz
      proj2 = dot2/absqr
      bcx = bcx - proj2*abx
      bcy = bcy - proj2*aby
      bcz = bcz - proj2*abz
c
c     compute the cross-product (aiix,aiiy,aiiz) x (bcx,bcy,bcz).  store
c     it in the vector (aibcx,aibcy,aibcz).
c
      aibcx = aiiy*bcz - aiiz*bcy
      aibcy = aiiz*bcx - aiix*bcz
      aibcz = aiix*bcy - aiiy*bcx
c
c     if (aibcx,aibcy,aibcz) points in the same direction as
c     (abx,aby,abz) then ic is located clockwise from i when
c     viewed from ia toward ib.  thus, in moving along the path
c     i-ia-ib-ic, a clockwise or positive angle is observed.
c     to determine whether these vectors point in the same or
c     opposite directions, compute their dot product.
c
      dot3 = aibcx*abx + aibcy*aby + aibcz*abz
      direc = sign(1.0d0,dot3)
c
c     compute the dihedral angle dih.
c
      dot4 = aiix*bcx + aiiy*bcy + aiiz*bcz
      aileng = sqrt(aiix**2 + aiiy**2 + aiiz**2)
      bcleng = sqrt(bcx**2 + bcy**2 + bcz**2)
      if(abs(aileng*bcleng).lt.1.0d-5)then
        cosine = 1.0d0
      else
        cosine = dot4/(aileng*bcleng)
      endif
      cosine = max(-1.0d0,cosine)
      cosine = min(1.0d0,cosine)
      dih = acos(cosine)*direc
      return
      end
