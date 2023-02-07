

!     -----------------------------------------------------------------
!     All of the ewald code and supporting routines were written and
!     contributed by Tom Darden from the National Institute of
!     Environmental Health Sciences division of the NIH.
!     Originally written with a modified version of AMBER 3A, the code
!     was updated during the summer of 1994 to be compatible with
!     AMBER 4.1.

!     This section contains code necessary to perform 3D FFTs where
!     libraries are not available.  It is based on piecing together a
!     series of 1D FFTs and is probably not super efficient.  The 1D
!     FFT code is a double precision version of fftpack from netlib,
!     written by Paul N. Swartztrauber at NCAR Boulder Colorado.

!     See www.netlib.org/fftpack

!     -----------------------------------------------------------------

!     Revisions: tec3 (added comments, modified source to fit it with
!     default Makefile, Compile, MACHINE source handling scheme, added
!     CPP selectable precision, added CPP control of SGI multiprocessing,
!     changed output formats to rid write(6,*), streamlined source,
!     merged routines, etc...)

!     NOTE: Some of the routine names within here are rather long so
!     may choke lame compilers.  This code was developed on SGI
!     computers so should surely work on these; it has also been
!     limitedly tested on Cray and HP machines.

!     The following C preprocessor code (only visible in the
!     untransformed source) is a hack to allow single precision
!     versions of an originally all double precision code.
!     It is recommended however that users run with double
!     precision... (by default sander is double precision)


!     The following routines are defined, in alphabetical order:
!     CFFTB
!     CFFTB1
!     CFFTF
!     CFFTF1
!     CFFTI
!     CFFTI1
!     PASSB
!     PASSB2
!     PASSB3
!     PASSB4
!     PASSB5
!     PASSF
!     PASSF2
!     PASSF3
!     PASSF4
!     PASSF5
!     PUBZ3D
!     PUBZ3DI


!     --- CFFTB ---

! Backward complex Transform


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cfftb here]
subroutine cfftb (n,c,wsave)
   implicit double precision (a-h,o-z)
   dimension       c(*)       ,wsave(*)
   if (n == 1) return
   iw1 = n+n+1
   iw2 = iw1+n+n
   call cfftb1 (n,c,wsave,wsave(iw1),wsave(iw2))
   return
end subroutine cfftb 

!     --- CFFTB1 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cfftb1 here]
subroutine cfftb1 (n,c,ch,wa,ifac)
   implicit double precision (a-h,o-z)
   dimension       ch(*)      ,c(*)       ,wa(*)      ,ifac(*)
   nf = ifac(2)
   na = 0
   l1 = 1
   iw = 1
   do 116 k1=1,nf
      ip = ifac(k1+2)
      l2 = ip*l1
      ido = n/l2
      idot = ido+ido
      idl1 = idot*l1
      if (ip /= 4) goto 103
      ix2 = iw+idot
      ix3 = ix2+idot
      if (na /= 0) goto 101
      call passb4 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
      goto 102
      101 call passb4 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
      102 na = 1-na
      goto 115
      103 if (ip /= 2) goto 106
      if (na /= 0) goto 104
      call passb2 (idot,l1,c,ch,wa(iw))
      goto 105
      104 call passb2 (idot,l1,ch,c,wa(iw))
      105 na = 1-na
      goto 115
      106 if (ip /= 3) goto 109
      ix2 = iw+idot
      if (na /= 0) goto 107
      call passb3 (idot,l1,c,ch,wa(iw),wa(ix2))
      goto 108
      107 call passb3 (idot,l1,ch,c,wa(iw),wa(ix2))
      108 na = 1-na
      goto 115
      109 if (ip /= 5) goto 112
      ix2 = iw+idot
      ix3 = ix2+idot
      ix4 = ix3+idot
      if (na /= 0) goto 110
      call passb5 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
      goto 111
      110 call passb5 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
      111 na = 1-na
      goto 115
      112 if (na /= 0) goto 113
      call passb (nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
      goto 114
      113 call passb (nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
      114 if (nac /= 0) na = 1-na
      115 l1 = l2
      iw = iw+(ip-1)*idot
   116 continue
   if (na == 0) return
   n2 = n+n
   do 117 i=1,n2
      c(i) = ch(i)
   117 continue
   return
end subroutine cfftb1 

!     --- CFFTF ---

! Forward complex Transform


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cfftf here]
subroutine cfftf (n,c,wsave)
   implicit double precision (a-h,o-z)
   dimension       c(*)       ,wsave(*)
   if (n == 1) return
   iw1 = n+n+1
   iw2 = iw1+n+n
   call cfftf1 (n,c,wsave,wsave(iw1),wsave(iw2))
   return
end subroutine cfftf 

!     --- CFFTF1 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cfftf1 here]
subroutine cfftf1 (n,c,ch,wa,ifac)
   implicit double precision (a-h,o-z)
   dimension       ch(*)      ,c(*)       ,wa(*)      ,ifac(*)
   nf = ifac(2)
   na = 0
   l1 = 1
   iw = 1
   do 116 k1=1,nf
      ip = ifac(k1+2)
      l2 = ip*l1
      ido = n/l2
      idot = ido+ido
      idl1 = idot*l1
      if (ip /= 4) goto 103
      ix2 = iw+idot
      ix3 = ix2+idot
      if (na /= 0) goto 101
      call passf4 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
      goto 102
      101 call passf4 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
      102 na = 1-na
      goto 115
      103 if (ip /= 2) goto 106
      if (na /= 0) goto 104
      call passf2 (idot,l1,c,ch,wa(iw))
      goto 105
      104 call passf2 (idot,l1,ch,c,wa(iw))
      105 na = 1-na
      goto 115
      106 if (ip /= 3) goto 109
      ix2 = iw+idot
      if (na /= 0) goto 107
      call passf3 (idot,l1,c,ch,wa(iw),wa(ix2))
      goto 108
      107 call passf3 (idot,l1,ch,c,wa(iw),wa(ix2))
      108 na = 1-na
      goto 115
      109 if (ip /= 5) goto 112
      ix2 = iw+idot
      ix3 = ix2+idot
      ix4 = ix3+idot
      if (na /= 0) goto 110
      call passf5 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
      goto 111
      110 call passf5 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
      111 na = 1-na
      goto 115
      112 if (na /= 0) goto 113
      call passf (nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
      goto 114
      113 call passf (nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
      114 if (nac /= 0) na = 1-na
      115 l1 = l2
      iw = iw+(ip-1)*idot
   116 continue
   if (na == 0) return
   n2 = n+n
   do 117 i=1,n2
      c(i) = ch(i)
   117 continue
   return
end subroutine cfftf1 

!     --- CFFTI ---

! Initialization for CFFTF and CFFTB

!******************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cffti here]
subroutine cffti(n,wsave)
   
   !******************************************************************
   
   ! subroutine cffti initializes the array wsave which is used in
   ! both cfftf and cfftb. the prime factorization of n together with
   ! a tabulation of the trigonometric functions are computed and
   ! stored in wsave.
   
   ! input parameter
   
   ! n       the length of the sequence to be transformed
   
   ! output parameter
   
   ! wsave   a work array which must be dimensioned at least 4*n+15
   !         the same work array can be used for both cfftf and cfftb
   !         as long as n remains unchanged. different wsave arrays
   !         are required for different values of n. the contents of
   !         wsave must not be changed between calls of cfftf or cfftb.

   implicit none
   integer      n
   double precision       wsave(*)

   integer      iw1
   integer      iw2

   if (n == 1) return
   iw1 = n+n+1
   iw2 = iw1+n+n
   call cffti1 (n,wsave(iw1),wsave(iw2))
   return
end subroutine cffti 

!     --- CFFTI1 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cffti1 here]
subroutine cffti1 (n,wa,ifac)
   implicit double precision (a-h,o-z)
   dimension       wa(*)      ,ifac(*)    ,ntryh(4)
   data ntryh(1),ntryh(2),ntryh(3),ntryh(4)/3,4,2,5/
   nl = n
   nf = 0
   j = 0
   101 j = j+1
   if (j-4) 102,102,103
   102 ntry = ntryh(j)
   goto 104
   103 ntry = ntry+2
   104 nq = nl/ntry
   nr = nl-ntry*nq
   if (nr) 101,105,101
   105 nf = nf+1
   ifac(nf+2) = ntry
   nl = nq
   if (ntry /= 2) goto 107
   if (nf == 1) goto 107
   do 106 i=2,nf
      ib = nf-i+2
      ifac(ib+2) = ifac(ib+1)
   106 continue
   ifac(3) = 2
   107 if (nl /= 1) goto 104
   ifac(1) = n
   ifac(2) = nf
   tpi = 6.28318530717959d0
   argh = tpi/dble(n)
   i = 2
   l1 = 1
   do 110 k1=1,nf
      ip = ifac(k1+2)
      ld = 0
      l2 = l1*ip
      ido = n/l2
      idot = ido+ido+2
      ipm = ip-1
      do 109 j=1,ipm
         i1 = i
         wa(i-1) = 1.d0
         wa(i) = 0.d0
         ld = ld+l1
         fi = 0.d0
         argld = dble(ld)*argh
         do 108 ii=4,idot,2
            i = i+2
            fi = fi+1.d0
            arg = fi*argld
            wa(i-1) = cos(arg)
            wa(i) = sin(arg)
         108 continue
         if (ip <= 5) goto 109
         wa(i1-1) = wa(i-1)
         wa(i1) = wa(i)
      109 continue
      l1 = l2
   110 continue
   return
end subroutine cffti1 

!     --- PASSB ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passb here]
subroutine passb (nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
   implicit double precision (a-h,o-z)
   dimension       ch(ido,l1,ip)          ,cc(ido,ip,l1)          , &
         c1(ido,l1,ip)          ,wa(*)      ,c2(idl1,ip), &
         ch2(idl1,ip)
   idot = ido/2
   nt = ip*idl1
   ipp2 = ip+2
   ipph = (ip+1)/2
   idp = ip*ido
   
   if (ido < l1) goto 106
   do 103 j=2,ipph
      jc = ipp2-j
      do 102 k=1,l1
         do 101 i=1,ido
            ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
            ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
         101 continue
      102 continue
   103 continue
   do 105 k=1,l1
      do 104 i=1,ido
         ch(i,k,1) = cc(i,1,k)
      104 continue
   105 continue
   goto 112
   106 do 109 j=2,ipph
      jc = ipp2-j
      do 108 i=1,ido
         do 107 k=1,l1
            ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
            ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
         107 continue
      108 continue
   109 continue
   do 111 i=1,ido
      do 110 k=1,l1
         ch(i,k,1) = cc(i,1,k)
      110 continue
   111 continue
   112 idl = 2-ido
   inc = 0
   do 116 l=2,ipph
      lc = ipp2-l
      idl = idl+ido
      do 113 ik=1,idl1
         c2(ik,l) = ch2(ik,1)+wa(idl-1)*ch2(ik,2)
         c2(ik,lc) = wa(idl)*ch2(ik,ip)
      113 continue
      idlj = idl
      inc = inc+ido
      do 115 j=3,ipph
         jc = ipp2-j
         idlj = idlj+inc
         if (idlj > idp) idlj = idlj-idp
         war = wa(idlj-1)
         wai = wa(idlj)
         do 114 ik=1,idl1
            c2(ik,l) = c2(ik,l)+war*ch2(ik,j)
            c2(ik,lc) = c2(ik,lc)+wai*ch2(ik,jc)
         114 continue
      115 continue
   116 continue
   do 118 j=2,ipph
      do 117 ik=1,idl1
         ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
      117 continue
   118 continue
   do 120 j=2,ipph
      jc = ipp2-j
      do 119 ik=2,idl1,2
         ch2(ik-1,j) = c2(ik-1,j)-c2(ik,jc)
         ch2(ik-1,jc) = c2(ik-1,j)+c2(ik,jc)
         ch2(ik,j) = c2(ik,j)+c2(ik-1,jc)
         ch2(ik,jc) = c2(ik,j)-c2(ik-1,jc)
      119 continue
   120 continue
   nac = 1
   if (ido == 2) return
   nac = 0
   do 121 ik=1,idl1
      c2(ik,1) = ch2(ik,1)
   121 continue
   do 123 j=2,ip
      do 122 k=1,l1
         c1(1,k,j) = ch(1,k,j)
         c1(2,k,j) = ch(2,k,j)
      122 continue
   123 continue
   if (idot > l1) goto 127
   idij = 0
   do 126 j=2,ip
      idij = idij+2
      do 125 i=4,ido,2
         idij = idij+2
         do 124 k=1,l1
            c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
            c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
         124 continue
      125 continue
   126 continue
   return
   127 idj = 2-ido
   do 130 j=2,ip
      idj = idj+ido
      do 129 k=1,l1
         idij = idj
         do 128 i=4,ido,2
            idij = idij+2
            c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
            c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
         128 continue
      129 continue
   130 continue
   return
end subroutine passb 

!     ---- PASSB2 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passb2 here]
subroutine passb2 (ido,l1,cc,ch,wa1)
   implicit double precision (a-h,o-z)
   dimension       cc(ido,2,l1)           ,ch(ido,l1,2)           , &
         wa1(*)
   if (ido > 2) goto 102
   do 101 k=1,l1
      ch(1,k,1) = cc(1,1,k)+cc(1,2,k)
      ch(1,k,2) = cc(1,1,k)-cc(1,2,k)
      ch(2,k,1) = cc(2,1,k)+cc(2,2,k)
      ch(2,k,2) = cc(2,1,k)-cc(2,2,k)
   101 continue
   return
   102 do 104 k=1,l1
      do 103 i=2,ido,2
         ch(i-1,k,1) = cc(i-1,1,k)+cc(i-1,2,k)
         tr2 = cc(i-1,1,k)-cc(i-1,2,k)
         ch(i,k,1) = cc(i,1,k)+cc(i,2,k)
         ti2 = cc(i,1,k)-cc(i,2,k)
         ch(i,k,2) = wa1(i-1)*ti2+wa1(i)*tr2
         ch(i-1,k,2) = wa1(i-1)*tr2-wa1(i)*ti2
      103 continue
   104 continue
   return
end subroutine passb2 

!     --- PASSB3 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passb3 here]
subroutine passb3 (ido,l1,cc,ch,wa1,wa2)
   implicit double precision (a-h,o-z)
   dimension       cc(ido,3,l1)           ,ch(ido,l1,3)           , &
         wa1(*)     ,wa2(*)
   data taur,taui /-.5d0,.866025403784439d0/
   if (ido /= 2) goto 102
   do 101 k=1,l1
      tr2 = cc(1,2,k)+cc(1,3,k)
      cr2 = cc(1,1,k)+taur*tr2
      ch(1,k,1) = cc(1,1,k)+tr2
      ti2 = cc(2,2,k)+cc(2,3,k)
      ci2 = cc(2,1,k)+taur*ti2
      ch(2,k,1) = cc(2,1,k)+ti2
      cr3 = taui*(cc(1,2,k)-cc(1,3,k))
      ci3 = taui*(cc(2,2,k)-cc(2,3,k))
      ch(1,k,2) = cr2-ci3
      ch(1,k,3) = cr2+ci3
      ch(2,k,2) = ci2+cr3
      ch(2,k,3) = ci2-cr3
   101 continue
   return
   102 do 104 k=1,l1
      do 103 i=2,ido,2
         tr2 = cc(i-1,2,k)+cc(i-1,3,k)
         cr2 = cc(i-1,1,k)+taur*tr2
         ch(i-1,k,1) = cc(i-1,1,k)+tr2
         ti2 = cc(i,2,k)+cc(i,3,k)
         ci2 = cc(i,1,k)+taur*ti2
         ch(i,k,1) = cc(i,1,k)+ti2
         cr3 = taui*(cc(i-1,2,k)-cc(i-1,3,k))
         ci3 = taui*(cc(i,2,k)-cc(i,3,k))
         dr2 = cr2-ci3
         dr3 = cr2+ci3
         di2 = ci2+cr3
         di3 = ci2-cr3
         ch(i,k,2) = wa1(i-1)*di2+wa1(i)*dr2
         ch(i-1,k,2) = wa1(i-1)*dr2-wa1(i)*di2
         ch(i,k,3) = wa2(i-1)*di3+wa2(i)*dr3
         ch(i-1,k,3) = wa2(i-1)*dr3-wa2(i)*di3
      103 continue
   104 continue
   return
end subroutine passb3 

!     --- PASSB4 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passb4 here]
subroutine passb4 (ido,l1,cc,ch,wa1,wa2,wa3)
   implicit double precision (a-h,o-z)
   dimension       cc(ido,4,l1)           ,ch(ido,l1,4)           , &
         wa1(*)     ,wa2(*)     ,wa3(*)
   if (ido /= 2) goto 102
   do 101 k=1,l1
      ti1 = cc(2,1,k)-cc(2,3,k)
      ti2 = cc(2,1,k)+cc(2,3,k)
      tr4 = cc(2,4,k)-cc(2,2,k)
      ti3 = cc(2,2,k)+cc(2,4,k)
      tr1 = cc(1,1,k)-cc(1,3,k)
      tr2 = cc(1,1,k)+cc(1,3,k)
      ti4 = cc(1,2,k)-cc(1,4,k)
      tr3 = cc(1,2,k)+cc(1,4,k)
      ch(1,k,1) = tr2+tr3
      ch(1,k,3) = tr2-tr3
      ch(2,k,1) = ti2+ti3
      ch(2,k,3) = ti2-ti3
      ch(1,k,2) = tr1+tr4
      ch(1,k,4) = tr1-tr4
      ch(2,k,2) = ti1+ti4
      ch(2,k,4) = ti1-ti4
   101 continue
   return
   102 do 104 k=1,l1
      do 103 i=2,ido,2
         ti1 = cc(i,1,k)-cc(i,3,k)
         ti2 = cc(i,1,k)+cc(i,3,k)
         ti3 = cc(i,2,k)+cc(i,4,k)
         tr4 = cc(i,4,k)-cc(i,2,k)
         tr1 = cc(i-1,1,k)-cc(i-1,3,k)
         tr2 = cc(i-1,1,k)+cc(i-1,3,k)
         ti4 = cc(i-1,2,k)-cc(i-1,4,k)
         tr3 = cc(i-1,2,k)+cc(i-1,4,k)
         ch(i-1,k,1) = tr2+tr3
         cr3 = tr2-tr3
         ch(i,k,1) = ti2+ti3
         ci3 = ti2-ti3
         cr2 = tr1+tr4
         cr4 = tr1-tr4
         ci2 = ti1+ti4
         ci4 = ti1-ti4
         ch(i-1,k,2) = wa1(i-1)*cr2-wa1(i)*ci2
         ch(i,k,2) = wa1(i-1)*ci2+wa1(i)*cr2
         ch(i-1,k,3) = wa2(i-1)*cr3-wa2(i)*ci3
         ch(i,k,3) = wa2(i-1)*ci3+wa2(i)*cr3
         ch(i-1,k,4) = wa3(i-1)*cr4-wa3(i)*ci4
         ch(i,k,4) = wa3(i-1)*ci4+wa3(i)*cr4
      103 continue
   104 continue
   return
end subroutine passb4 

!     --- PASSB5 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passb5 here]
subroutine passb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
   implicit double precision (a-h,o-z)
   dimension       cc(ido,5,l1)           ,ch(ido,l1,5)           , &
         wa1(*)     ,wa2(*)     ,wa3(*)     ,wa4(*)
   data tr11,ti11,tr12,ti12 /.309016994374947d0, &
         .951056516295154d0, &
         -.809016994374947d0,.587785252292473d0/
   if (ido /= 2) goto 102
   do 101 k=1,l1
      ti5 = cc(2,2,k)-cc(2,5,k)
      ti2 = cc(2,2,k)+cc(2,5,k)
      ti4 = cc(2,3,k)-cc(2,4,k)
      ti3 = cc(2,3,k)+cc(2,4,k)
      tr5 = cc(1,2,k)-cc(1,5,k)
      tr2 = cc(1,2,k)+cc(1,5,k)
      tr4 = cc(1,3,k)-cc(1,4,k)
      tr3 = cc(1,3,k)+cc(1,4,k)
      ch(1,k,1) = cc(1,1,k)+tr2+tr3
      ch(2,k,1) = cc(2,1,k)+ti2+ti3
      cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
      ci2 = cc(2,1,k)+tr11*ti2+tr12*ti3
      cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
      ci3 = cc(2,1,k)+tr12*ti2+tr11*ti3
      cr5 = ti11*tr5+ti12*tr4
      ci5 = ti11*ti5+ti12*ti4
      cr4 = ti12*tr5-ti11*tr4
      ci4 = ti12*ti5-ti11*ti4
      ch(1,k,2) = cr2-ci5
      ch(1,k,5) = cr2+ci5
      ch(2,k,2) = ci2+cr5
      ch(2,k,3) = ci3+cr4
      ch(1,k,3) = cr3-ci4
      ch(1,k,4) = cr3+ci4
      ch(2,k,4) = ci3-cr4
      ch(2,k,5) = ci2-cr5
   101 continue
   return
   102 do 104 k=1,l1
      do 103 i=2,ido,2
         ti5 = cc(i,2,k)-cc(i,5,k)
         ti2 = cc(i,2,k)+cc(i,5,k)
         ti4 = cc(i,3,k)-cc(i,4,k)
         ti3 = cc(i,3,k)+cc(i,4,k)
         tr5 = cc(i-1,2,k)-cc(i-1,5,k)
         tr2 = cc(i-1,2,k)+cc(i-1,5,k)
         tr4 = cc(i-1,3,k)-cc(i-1,4,k)
         tr3 = cc(i-1,3,k)+cc(i-1,4,k)
         ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
         ch(i,k,1) = cc(i,1,k)+ti2+ti3
         cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
         ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
         cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
         ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         dr3 = cr3-ci4
         dr4 = cr3+ci4
         di3 = ci3+cr4
         di4 = ci3-cr4
         dr5 = cr2+ci5
         dr2 = cr2-ci5
         di5 = ci2-cr5
         di2 = ci2+cr5
         ch(i-1,k,2) = wa1(i-1)*dr2-wa1(i)*di2
         ch(i,k,2) = wa1(i-1)*di2+wa1(i)*dr2
         ch(i-1,k,3) = wa2(i-1)*dr3-wa2(i)*di3
         ch(i,k,3) = wa2(i-1)*di3+wa2(i)*dr3
         ch(i-1,k,4) = wa3(i-1)*dr4-wa3(i)*di4
         ch(i,k,4) = wa3(i-1)*di4+wa3(i)*dr4
         ch(i-1,k,5) = wa4(i-1)*dr5-wa4(i)*di5
         ch(i,k,5) = wa4(i-1)*di5+wa4(i)*dr5
      103 continue
   104 continue
   return
end subroutine passb5 

!     --- PASSF ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passf here]
subroutine passf (nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
   implicit double precision (a-h,o-z)
   dimension       ch(ido,l1,ip)          ,cc(ido,ip,l1)          , &
         c1(ido,l1,ip)          ,wa(*)      ,c2(idl1,ip), &
         ch2(idl1,ip)
   idot = ido/2
   nt = ip*idl1
   ipp2 = ip+2
   ipph = (ip+1)/2
   idp = ip*ido
   
   if (ido < l1) goto 106
   do 103 j=2,ipph
      jc = ipp2-j
      do 102 k=1,l1
         do 101 i=1,ido
            ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
            ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
         101 continue
      102 continue
   103 continue
   do 105 k=1,l1
      do 104 i=1,ido
         ch(i,k,1) = cc(i,1,k)
      104 continue
   105 continue
   goto 112
   106 do 109 j=2,ipph
      jc = ipp2-j
      do 108 i=1,ido
         do 107 k=1,l1
            ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
            ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
         107 continue
      108 continue
   109 continue
   do 111 i=1,ido
      do 110 k=1,l1
         ch(i,k,1) = cc(i,1,k)
      110 continue
   111 continue
   112 idl = 2-ido
   inc = 0
   do 116 l=2,ipph
      lc = ipp2-l
      idl = idl+ido
      do 113 ik=1,idl1
         c2(ik,l) = ch2(ik,1)+wa(idl-1)*ch2(ik,2)
         c2(ik,lc) = -wa(idl)*ch2(ik,ip)
      113 continue
      idlj = idl
      inc = inc+ido
      do 115 j=3,ipph
         jc = ipp2-j
         idlj = idlj+inc
         if (idlj > idp) idlj = idlj-idp
         war = wa(idlj-1)
         wai = wa(idlj)
         do 114 ik=1,idl1
            c2(ik,l) = c2(ik,l)+war*ch2(ik,j)
            c2(ik,lc) = c2(ik,lc)-wai*ch2(ik,jc)
         114 continue
      115 continue
   116 continue
   do 118 j=2,ipph
      do 117 ik=1,idl1
         ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
      117 continue
   118 continue
   do 120 j=2,ipph
      jc = ipp2-j
      do 119 ik=2,idl1,2
         ch2(ik-1,j) = c2(ik-1,j)-c2(ik,jc)
         ch2(ik-1,jc) = c2(ik-1,j)+c2(ik,jc)
         ch2(ik,j) = c2(ik,j)+c2(ik-1,jc)
         ch2(ik,jc) = c2(ik,j)-c2(ik-1,jc)
      119 continue
   120 continue
   nac = 1
   if (ido == 2) return
   nac = 0
   do 121 ik=1,idl1
      c2(ik,1) = ch2(ik,1)
   121 continue
   do 123 j=2,ip
      do 122 k=1,l1
         c1(1,k,j) = ch(1,k,j)
         c1(2,k,j) = ch(2,k,j)
      122 continue
   123 continue
   if (idot > l1) goto 127
   idij = 0
   do 126 j=2,ip
      idij = idij+2
      do 125 i=4,ido,2
         idij = idij+2
         do 124 k=1,l1
            c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)+wa(idij)*ch(i,k,j)
            c1(i,k,j) = wa(idij-1)*ch(i,k,j)-wa(idij)*ch(i-1,k,j)
         124 continue
      125 continue
   126 continue
   return
   127 idj = 2-ido
   do 130 j=2,ip
      idj = idj+ido
      do 129 k=1,l1
         idij = idj
         do 128 i=4,ido,2
            idij = idij+2
            c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)+wa(idij)*ch(i,k,j)
            c1(i,k,j) = wa(idij-1)*ch(i,k,j)-wa(idij)*ch(i-1,k,j)
         128 continue
      129 continue
   130 continue
   return
end subroutine passf 

!     --- PASSF2 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passf2 here]
subroutine passf2 (ido,l1,cc,ch,wa1)
   implicit double precision (a-h,o-z)
   dimension       cc(ido,2,l1)           ,ch(ido,l1,2)           , &
         wa1(*)
   if (ido > 2) goto 102
   do 101 k=1,l1
      ch(1,k,1) = cc(1,1,k)+cc(1,2,k)
      ch(1,k,2) = cc(1,1,k)-cc(1,2,k)
      ch(2,k,1) = cc(2,1,k)+cc(2,2,k)
      ch(2,k,2) = cc(2,1,k)-cc(2,2,k)
   101 continue
   return
   102 do 104 k=1,l1
      do 103 i=2,ido,2
         ch(i-1,k,1) = cc(i-1,1,k)+cc(i-1,2,k)
         tr2 = cc(i-1,1,k)-cc(i-1,2,k)
         ch(i,k,1) = cc(i,1,k)+cc(i,2,k)
         ti2 = cc(i,1,k)-cc(i,2,k)
         ch(i,k,2) = wa1(i-1)*ti2-wa1(i)*tr2
         ch(i-1,k,2) = wa1(i-1)*tr2+wa1(i)*ti2
      103 continue
   104 continue
   return
end subroutine passf2 

!     --- PASSF3 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passf3 here]
subroutine passf3 (ido,l1,cc,ch,wa1,wa2)
   implicit double precision (a-h,o-z)
   dimension       cc(ido,3,l1)           ,ch(ido,l1,3)           , &
         wa1(*)     ,wa2(*)
   data taur,taui /-.5d0,-.866025403784439d0/
   if (ido /= 2) goto 102
   do 101 k=1,l1
      tr2 = cc(1,2,k)+cc(1,3,k)
      cr2 = cc(1,1,k)+taur*tr2
      ch(1,k,1) = cc(1,1,k)+tr2
      ti2 = cc(2,2,k)+cc(2,3,k)
      ci2 = cc(2,1,k)+taur*ti2
      ch(2,k,1) = cc(2,1,k)+ti2
      cr3 = taui*(cc(1,2,k)-cc(1,3,k))
      ci3 = taui*(cc(2,2,k)-cc(2,3,k))
      ch(1,k,2) = cr2-ci3
      ch(1,k,3) = cr2+ci3
      ch(2,k,2) = ci2+cr3
      ch(2,k,3) = ci2-cr3
   101 continue
   return
   102 do 104 k=1,l1
      do 103 i=2,ido,2
         tr2 = cc(i-1,2,k)+cc(i-1,3,k)
         cr2 = cc(i-1,1,k)+taur*tr2
         ch(i-1,k,1) = cc(i-1,1,k)+tr2
         ti2 = cc(i,2,k)+cc(i,3,k)
         ci2 = cc(i,1,k)+taur*ti2
         ch(i,k,1) = cc(i,1,k)+ti2
         cr3 = taui*(cc(i-1,2,k)-cc(i-1,3,k))
         ci3 = taui*(cc(i,2,k)-cc(i,3,k))
         dr2 = cr2-ci3
         dr3 = cr2+ci3
         di2 = ci2+cr3
         di3 = ci2-cr3
         ch(i,k,2) = wa1(i-1)*di2-wa1(i)*dr2
         ch(i-1,k,2) = wa1(i-1)*dr2+wa1(i)*di2
         ch(i,k,3) = wa2(i-1)*di3-wa2(i)*dr3
         ch(i-1,k,3) = wa2(i-1)*dr3+wa2(i)*di3
      103 continue
   104 continue
   return
end subroutine passf3 

!     --- PASSF4 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passf4 here]
subroutine passf4 (ido,l1,cc,ch,wa1,wa2,wa3)
   implicit double precision (a-h,o-z)
   dimension       cc(ido,4,l1)           ,ch(ido,l1,4)           , &
         wa1(*)     ,wa2(*)     ,wa3(*)
   if (ido /= 2) goto 102
   do 101 k=1,l1
      ti1 = cc(2,1,k)-cc(2,3,k)
      ti2 = cc(2,1,k)+cc(2,3,k)
      tr4 = cc(2,2,k)-cc(2,4,k)
      ti3 = cc(2,2,k)+cc(2,4,k)
      tr1 = cc(1,1,k)-cc(1,3,k)
      tr2 = cc(1,1,k)+cc(1,3,k)
      ti4 = cc(1,4,k)-cc(1,2,k)
      tr3 = cc(1,2,k)+cc(1,4,k)
      ch(1,k,1) = tr2+tr3
      ch(1,k,3) = tr2-tr3
      ch(2,k,1) = ti2+ti3
      ch(2,k,3) = ti2-ti3
      ch(1,k,2) = tr1+tr4
      ch(1,k,4) = tr1-tr4
      ch(2,k,2) = ti1+ti4
      ch(2,k,4) = ti1-ti4
   101 continue
   return
   102 do 104 k=1,l1
      do 103 i=2,ido,2
         ti1 = cc(i,1,k)-cc(i,3,k)
         ti2 = cc(i,1,k)+cc(i,3,k)
         ti3 = cc(i,2,k)+cc(i,4,k)
         tr4 = cc(i,2,k)-cc(i,4,k)
         tr1 = cc(i-1,1,k)-cc(i-1,3,k)
         tr2 = cc(i-1,1,k)+cc(i-1,3,k)
         ti4 = cc(i-1,4,k)-cc(i-1,2,k)
         tr3 = cc(i-1,2,k)+cc(i-1,4,k)
         ch(i-1,k,1) = tr2+tr3
         cr3 = tr2-tr3
         ch(i,k,1) = ti2+ti3
         ci3 = ti2-ti3
         cr2 = tr1+tr4
         cr4 = tr1-tr4
         ci2 = ti1+ti4
         ci4 = ti1-ti4
         ch(i-1,k,2) = wa1(i-1)*cr2+wa1(i)*ci2
         ch(i,k,2) = wa1(i-1)*ci2-wa1(i)*cr2
         ch(i-1,k,3) = wa2(i-1)*cr3+wa2(i)*ci3
         ch(i,k,3) = wa2(i-1)*ci3-wa2(i)*cr3
         ch(i-1,k,4) = wa3(i-1)*cr4+wa3(i)*ci4
         ch(i,k,4) = wa3(i-1)*ci4-wa3(i)*cr4
      103 continue
   104 continue
   return
end subroutine passf4 

!     --- PASSF5 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passf5 here]
subroutine passf5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
   implicit double precision (a-h,o-z)
   dimension       cc(ido,5,l1)           ,ch(ido,l1,5)           , &
         wa1(*)     ,wa2(*)     ,wa3(*)     ,wa4(*)
   data tr11,ti11,tr12,ti12 /.309016994374947d0, &
         -.951056516295154d0, &
         -.809016994374947d0,-.587785252292473d0/
   if (ido /= 2) goto 102
   do 101 k=1,l1
      ti5 = cc(2,2,k)-cc(2,5,k)
      ti2 = cc(2,2,k)+cc(2,5,k)
      ti4 = cc(2,3,k)-cc(2,4,k)
      ti3 = cc(2,3,k)+cc(2,4,k)
      tr5 = cc(1,2,k)-cc(1,5,k)
      tr2 = cc(1,2,k)+cc(1,5,k)
      tr4 = cc(1,3,k)-cc(1,4,k)
      tr3 = cc(1,3,k)+cc(1,4,k)
      ch(1,k,1) = cc(1,1,k)+tr2+tr3
      ch(2,k,1) = cc(2,1,k)+ti2+ti3
      cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
      ci2 = cc(2,1,k)+tr11*ti2+tr12*ti3
      cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
      ci3 = cc(2,1,k)+tr12*ti2+tr11*ti3
      cr5 = ti11*tr5+ti12*tr4
      ci5 = ti11*ti5+ti12*ti4
      cr4 = ti12*tr5-ti11*tr4
      ci4 = ti12*ti5-ti11*ti4
      ch(1,k,2) = cr2-ci5
      ch(1,k,5) = cr2+ci5
      ch(2,k,2) = ci2+cr5
      ch(2,k,3) = ci3+cr4
      ch(1,k,3) = cr3-ci4
      ch(1,k,4) = cr3+ci4
      ch(2,k,4) = ci3-cr4
      ch(2,k,5) = ci2-cr5
   101 continue
   return
   102 do 104 k=1,l1
      do 103 i=2,ido,2
         ti5 = cc(i,2,k)-cc(i,5,k)
         ti2 = cc(i,2,k)+cc(i,5,k)
         ti4 = cc(i,3,k)-cc(i,4,k)
         ti3 = cc(i,3,k)+cc(i,4,k)
         tr5 = cc(i-1,2,k)-cc(i-1,5,k)
         tr2 = cc(i-1,2,k)+cc(i-1,5,k)
         tr4 = cc(i-1,3,k)-cc(i-1,4,k)
         tr3 = cc(i-1,3,k)+cc(i-1,4,k)
         ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
         ch(i,k,1) = cc(i,1,k)+ti2+ti3
         cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
         ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
         cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
         ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         dr3 = cr3-ci4
         dr4 = cr3+ci4
         di3 = ci3+cr4
         di4 = ci3-cr4
         dr5 = cr2+ci5
         dr2 = cr2-ci5
         di5 = ci2-cr5
         di2 = ci2+cr5
         ch(i-1,k,2) = wa1(i-1)*dr2+wa1(i)*di2
         ch(i,k,2) = wa1(i-1)*di2-wa1(i)*dr2
         ch(i-1,k,3) = wa2(i-1)*dr3+wa2(i)*di3
         ch(i,k,3) = wa2(i-1)*di3-wa2(i)*dr3
         ch(i-1,k,4) = wa3(i-1)*dr4+wa3(i)*di4
         ch(i,k,4) = wa3(i-1)*di4-wa3(i)*dr4
         ch(i-1,k,5) = wa4(i-1)*dr5+wa4(i)*di5
         ch(i,k,5) = wa4(i-1)*di5-wa4(i)*dr5
      103 continue
   104 continue
   return
end subroutine passf5 

!     --- PUBZ3D ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine pubz3d here]
subroutine pubz3d(isign,n1,n2,n3,w,ld1,ld2,table,ntable, &
      work,nwork)
   implicit none

   integer n1,n2,n3,ld1,ld2,isign,ntable,nwork
   double complex w(ld1,ld2,n3)
   double complex work( nwork)
   double precision table(ntable,3)

   integer i,j,k
   !     ...note: ntable should be 4*max(n1,n2,n3) +15
   !              nwork should be max(n1,n2,n3)

   !     ...transform along X  first
   
   do 100 k = 1, n3
      do 90 j = 1, n2
         do 70 i = 1,n1
            work(i) = w(i,j,k)
         70 continue
         if ( isign == -1) call cfftf(n1,work,table(1,1))
         if ( isign == 1) call cfftb(n1,work,table(1,1))
         do 80 i = 1,n1
            w(i,j,k) = work(i)
         80 continue
      90 continue
   100 continue
   
   !     ...transform along Y then
   
   do 200 k = 1,n3
      do 190 i = 1,n1
         do 170 j = 1,n2
            work(j) = w(i,j,k)
         170 continue
         if ( isign == -1) call cfftf(n2,work,table(1,2))
         if ( isign == 1) call cfftb(n2,work,table(1,2))
         do 180 j = 1,n2
            w(i,j,k) = work(j)
         180 continue
      190 continue
   200 continue
   
   !     ...transform along Z finally
   
   do 300 i = 1, n1
      do 290 j = 1, n2
         do 270 k = 1,n3
            work(k) = w(i,j,k)
         270 continue
         if ( isign == -1) call cfftf(n3,work,table(1,3))
         if ( isign == 1) call cfftb(n3,work,table(1,3))
         do 280 k = 1,n3
            w(i,j,k) = work(k)
         280 continue
      290 continue
   300 continue

   return
end subroutine pubz3d 

!     --- PUBZ3DI ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine pubz3di here]
subroutine pubz3di(n1,n2,n3,table,ntable)
   implicit none
   integer n1,n2,n3,ntable
   double precision table(ntable,3)
   !     NOTE: ntable should be 4*max(n1,n2,n3) +15

   call cffti(n1,table(1,1))
   call cffti(n2,table(1,2))
   call cffti(n3,table(1,3))

   return
end subroutine pubz3di 
