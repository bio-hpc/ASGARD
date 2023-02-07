
!---------------------------------------------------------------------


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine corffft here]
subroutine corffft(ntrun, nvect, vect, ipair1, ipair2, &
      cftmp, p2cftmp, rcftmp, &
      cf, p2cf, rcf, ndata, data, table, &
      rave, r3iave, r6iave, avcrd)
   
   !     ----- If ntrun = 9, calculates ired correlation functions
   !             (see JACS 2002, 124, 4522, eqs. A17 and A18,
   !              via Fast Fourier Transforms).
   !           cftmp is used to store a(m,l,t):
   !             a(0, ..., m-1; -L,...,L; 0, ..., T-1).
   !           cf contains the values of the correlation
   !             functions for modes ibeg to iend as:
   !             cf(interval * nof_modes + modenumber).
   !     ----- If ntrun = 10, calculates vector correlation functions
   !             (see ...).
   !           cftmp is used to store Y(p,l,t)/r**3:
   !             Y(0, ..., npairs-1; -L,...,L; 0, ..., T-1).
   !           p2cftmp is used to store Y(p,l,t).
   !           rcftmp is used to store 1/r**3(p,t).
   !           cf contains the values of the correlation
   !             functions for pairs 1 to npairs as:
   !             cf(interval * nof_pairs + pair).
   !           Same p2cf and rcf.
   !     ----- data is used for calculation of FFT
   
   !     ----- HG 11/04/2002
   
   implicit double precision (a-h,o-z)
#  include "sizes2.h"
#  include "anal.h"
#  include "files.h"
#  include "infoa.h"
   integer ntrun, nvect, ipair1, ipair2, ndata
   real crd, box
   double precision vect, cftmp, p2cftmp, rcftmp
   double precision cf, p2cf, rcf, data
   double precision rave, r3iave, r6iave, avcrd
   character(len=40) fmt, fmt2
   character(len=3)  chno, chno2
   logical hasbox
   parameter(pi=3.141592654)
   dimension ipair1(*), ipair2(*)
   dimension vect(npairs,*), cftmp(*), p2cftmp(*), rcftmp(*)
   dimension cf(*), p2cf(*), rcf(*), data(*), table(*)
   dimension rave(*), r3iave(*), r6iave(*), avcrd(*)
   dimension crd(3*maxatom), box(3)
   common/tcor/ tmax, tintvl, intname(maxint)
   common/ired/ iorder, npairs
   
   !     --- Define Function
   
   k(i,j) = 3*(i-1) + j
   
   !     --- Init parameters / arrays
   
   intvals = int(tmax/tintvl) + 1
   nveccalc = iend - ibeg + 1
   mtot = 2 * iorder + 1
   nsnap = int((last - first + 1) / iskip)
   fac = sqrt(dble(mtot) / (4.0d0 * pi))
   !c      write(6,*) intvals, nveccalc, mtot, nsnap, fac
   
   if(ntrun == 9) then
      nch  = iend - ibeg + 2
   else if(ntrun == 10) then
      nch  = npairs + 1
      nch2 = npairs
   end if
   if (nch < 10) then
      write(chno,'(i1)') nch
   else if (nch < 100) then
      write(chno,'(i2)') nch
   else if (nch < 1000) then
      write(chno,'(i3)') nch
   else
      write(0,*) 'nch is too big: ', nch
      call mexit(6,1)
   end if
   fmt = '(' // chno // 'f9.4)'
   if(ntrun == 10) then
      if (nch2 < 10) then
         write(chno2,'(i1)') nch2
      else if (nch2 < 100) then
         write(chno2,'(i2)') nch2
      else if (nch2 < 1000) then
         write(chno2,'(i3)') nch2
      else
         write(0,*) 'nch2 is too big: ', nch2
         call mexit(6,1)
      end if
      fmt2 = '(A9,' // chno2 // 'f9.4)'
   end if
   
   if(ntrun == 9) then
      do i=1,2*nsnap*mtot*nveccalc
         cftmp(i) = 0.0d0
      end do
      do i=1,nveccalc*intvals
         cf(i) = 0.0d0
      end do
   else if(ntrun == 10) then
      do i=1,2*nsnap*mtot*npairs
         cftmp(i) = 0.0d0
         p2cftmp(i) = 0.0d0
      end do
      do i=1,2*nsnap*npairs
         rcftmp(i) = 0.0d0
      end do
      do i=1,npairs*intvals
         cf(i) = 0.0d0
         p2cf(i) = 0.0d0
         rcf(i) = 0.0d0
      end do
      do i=1,npairs
         rave(i) = 0.0d0
         r3iave(i) = 0.0d0
         r6iave(i) = 0.0d0
         avcrd(k(i,1)) = 0.0d0
         avcrd(k(i,2)) = 0.0d0
         avcrd(k(i,3)) = 0.0d0
      end do
   end if
   
   !     --- Open the input file
   !     ---   If iftyp != 2, check if additional arrays are
   !     ---   necessary for readfile routine below.
   iftyp = 2
   hasbox = .false.
   call openinp(iftyp)
   
   !     --- Loop over snapshots
   
   ierr = 0
   ieof = 0
   isnap = 0
   do isp=1,9999999
      call readfile(natom,nread,ener,crd,vel,box,iftyp, &
            hasbox,ierr,ieof)
      if(ierr == 1) then
         write(6,*) 'Error during read: natom=',natom,' nread=',nread
         call mexit(6,1)
      end if
      if(ieof == 1) goto 30
      if(isp > last) goto 30
      if(isp >= first .and. mod(isp,iskip) == 0) then
         isnap = isnap + 1
         if(mod(isnap,50) == 0) write(6,100) isnap
         if(ntrun == 9) then
            indsnap = nveccalc*mtot*(isnap-1)
         else if(ntrun == 10) then
            indsnap  = npairs*mtot*(isnap-1)
            indsnap2 = npairs*(isnap-1)
         end if
         !c          write(6,*) 'Indsnap: ', isnap-1, indsnap
         
         !         --- Loop over all vectors
         
         do npair=1,npairs
            indpair = npair
            iat1 = ipair1(npair)
            iat2 = ipair2(npair)
            xcrd = crd(k(iat2,1)) - crd(k(iat1,1))
            ycrd = crd(k(iat2,2)) - crd(k(iat1,2))
            zcrd = crd(k(iat2,3)) - crd(k(iat1,3))
            r = sqrt(xcrd*xcrd + ycrd*ycrd + zcrd*zcrd)
            if(ntrun == 10) then
               avcrd(k(npair,1)) = avcrd(k(npair,1)) + xcrd
               avcrd(k(npair,2)) = avcrd(k(npair,2)) + ycrd
               avcrd(k(npair,3)) = avcrd(k(npair,3)) + zcrd
               rave(npair)  = rave(npair)  + r
               r3 = r*r*r
               r3i = 1.0d0 / r3
               r3iave(npair) = r3iave(npair) + r3i
               r6iave(npair) = r6iave(npair) + r3i*r3i
            end if
            
            !           --- Loop over m=0, ..., +L
            !               (Caveat: m stands for the mode in Bruschweiler paper)
            
            do m=0,iorder
               
               !             --- Calc spherical harmonics
               !                 (for iorder <= 2, hard-coded spherical harmonics in
               !                  terms of cartesian coordinates (e.g. Merzbacher,
               !                  Quantum Mechanics, p. 186) are used)
               
               if(ntrun == 9) then
                  indplus = nveccalc * (iorder + m)
               else if(ntrun == 10) then
                  indplus = npairs * (iorder + m)
               end if
               call spherharm(iorder, m, xcrd, ycrd, zcrd, &
                     r, dplusreal, dplusimg)
               
               if(m > 0) then
                  if(ntrun == 9) then
                     indminus = nveccalc * (iorder - m)
                  else if(ntrun == 10) then
                     indminus = npairs * (iorder - m)
                  end if
                  call spherharm(iorder, -m, xcrd, ycrd, zcrd, &
                        r, dminusreal, dminusimg)
               end if
               
               if(ntrun == 9) then
                  
                  !               --- Loop over all eigenvectors
                  
                  do imode=1,nveccalc
                     indmode = imode
                     q = vect(npair,ibeg + imode - 1)

                     indtot = 2 * (indsnap + indplus + indmode)
                     cftmp(indtot - 1) = cftmp(indtot - 1) + q*dplusreal
                     cftmp(indtot    ) = cftmp(indtot    ) + q*dplusimg
                     if(m > 0) then
                        indtot = 2 * (indsnap + indminus + indmode)
                        cftmp(indtot - 1) = cftmp(indtot - 1) + q*dminusreal
                        cftmp(indtot    ) = cftmp(indtot    ) + q*dminusimg
                     end if
                  end do
               else if(ntrun == 10) then
                  indtot = 2 * (indsnap + indplus + indpair)
                  cftmp(indtot - 1) = cftmp(indtot - 1) + r3i*dplusreal
                  cftmp(indtot    ) = cftmp(indtot    ) + r3i*dplusimg
                  p2cftmp(indtot - 1) = p2cftmp(indtot - 1) + dplusreal
                  p2cftmp(indtot    ) = p2cftmp(indtot    ) + dplusimg
                  if(m > 0) then
                     indtot = 2 * (indsnap + indminus + indpair)
                     cftmp(indtot - 1) = cftmp(indtot - 1) + r3i*dminusreal
                     cftmp(indtot    ) = cftmp(indtot    ) + r3i*dminusimg
                     p2cftmp(indtot - 1) = p2cftmp(indtot - 1) + dminusreal
                     p2cftmp(indtot    ) = p2cftmp(indtot    ) + dminusimg
                  else if(m == 0) then
                     indtot2 = 2 * (indsnap2 + indpair)
                     rcftmp(indtot2 - 1) = rcftmp(indtot2 - 1) + r3i
                     rcftmp(indtot2    ) = 0.0d0
                  end if
               end if  ! (ntrun == 9)
            end do  !  m=0,iorder
         end do  !  npair=1,npairs
      end if  ! (isp >= first .and. mod(isp,iskip) == 0)
   end do  !  isp=1,9999999
   
   !     --- Post-process
   
   30 continue
   
   !     --- Init FFT
   
   call cffti(ndata/2, table)
   
   if(ntrun == 9) then
      
      !       --- Loop over all eigenvectors
      
      do imode=1,nveccalc
         indmode = imode
         indmode2 = intvals * (imode - 1)
         
         !         --- Loop over all m=-L, ..., L
         
         do m=0,mtot-1
            indm = nveccalc * m
            
            !           --- Loop over all snapshots and fill up with zero's
            
            do isp=1,isnap
               indsnap = nveccalc*mtot*(isp-1)
               indtot = 2 * (indsnap + indm + indmode)
               data(2 * isp - 1) = cftmp(indtot - 1)
               data(2 * isp    ) = cftmp(indtot    )
            end do
            do i=2*isnap+1,ndata
               data(i) = 0.0d0
            end do
            
            !           --- Calc correlation function (= C(m,l,t) in Bruschweiler paper)
            
            call corffft2(ndata, data, table)
            
            !           --- Sum into cf (= C(m,t) in Bruschweiler paper)
            
            do nint=1,intvals
               cf(indmode2 + nint) = cf(indmode2 + nint) + data(2*nint-1)
            end do
         end do
      end do  !  imode=1,nveccalc
   else if(ntrun == 10) then
      
      !     --- Loop over all three correlation functions
      
      do i=1,3
         
         !         --- Loop over all pairs
         
         do npair=1,npairs
            indpair  = npair
            indpair2 = intvals * (npair - 1)
            
            !           --- Loop over all m=-L, ..., L
            
            do m=0,mtot-1
               indm = npairs * m
               
               !             --- Loop over all snapshots and fill up with zero's
               
               do isp=1,isnap
                  if(i == 1 .or. i == 2) then
                     indsnap = npairs*mtot*(isp-1)
                     indtot = 2 * (indsnap + indm + indpair)
                     if(i == 1) then
                        data(2 * isp - 1) = cftmp(indtot - 1)
                        data(2 * isp    ) = cftmp(indtot    )
                     else if(i == 2) then
                        data(2 * isp - 1) = p2cftmp(indtot - 1)
                        data(2 * isp    ) = p2cftmp(indtot    )
                     end if
                  else if(i == 3 .and. m == 0) then
                     indsnap2 = npairs*(isp-1)
                     indtot2 = 2 * (indsnap2 + indpair)
                     data(2 * isp - 1) = rcftmp(indtot2 - 1)
                     data(2 * isp    ) = rcftmp(indtot2    )
                  end if
               end do
               do j=2*isnap+1,ndata
                  data(j) = 0.0d0
               end do
               
               !             --- Calc correlation function
               
               call corffft2(ndata, data, table)
               
               !             --- Sum into cf
               
               do nint=1,intvals
                  if(i == 1) then
                     cf(indpair2 + nint) = cf(indpair2 + nint) + &
                           data(2*nint-1)
                  else if(i == 2) then
                     p2cf(indpair2 + nint) = p2cf(indpair2 + nint) + &
                           data(2*nint-1)
                  else if(i == 3 .and. m == 0) then
                     rcf(indpair2 + nint) = rcf(indpair2 + nint) + &
                           data(2*nint-1)
                  end if
               end do
            end do  !  m=0,mtot-1
         end do  !  npair=1,npairs
      end do  !  i=1,3
   end if  ! (ntrun == 9)
   
   !     --- Normalize and output
   
   if(ntrun == 9) then
      write(6,*)
      write(6,*) '***** IRED corr func   *****'
      do nint=1,intvals
         dnorm = 1.0d0 / dble(isnap - nint + 1)
         do imode=1,nveccalc
            ind = intvals * (imode-1) + nint
            cf(ind) = cf(ind) * dnorm
         end do
         write(6,fmt) dble(nint - 1)*tintvl, &
               (cf(intvals * (imode-1) + nint),imode=1,nveccalc)
      end do
   else if(ntrun == 10) then
      do nint=1,intvals
         dnorm = 1.0d0 / dble(isnap - nint + 1)
         do npair=1,npairs
            ind = intvals * (npair-1) + nint
            cf(ind) = cf(ind) * dnorm
            p2cf(ind) = p2cf(ind) * dnorm
            rcf(ind) = rcf(ind) * dnorm
         end do
      end do
      
      dnorm = 1.0d0 / dble(isnap)
      do npair=1,npairs
         rave(npair)   = rave(npair)   * dnorm
         r3iave(npair) = r3iave(npair) * dnorm
         r6iave(npair) = r6iave(npair) * dnorm
         avcrd(k(npair,1)) = avcrd(k(npair,1)) * dnorm
         avcrd(k(npair,2)) = avcrd(k(npair,2)) * dnorm
         avcrd(k(npair,3)) = avcrd(k(npair,3)) * dnorm
      end do
      
      write(6,*)
      write(6,*) '***** Vector length    *****'
      write(6,fmt2) 'RAVE', (rave(npair),npair=1,npairs)
      write(6,fmt2) 'RRIG', &
            (sqrt(avcrd(k(npair,1))*avcrd(k(npair,1)) + &
            avcrd(k(npair,2))*avcrd(k(npair,2)) + &
            avcrd(k(npair,3))*avcrd(k(npair,3))), &
            npair=1,npairs)
      write(6,fmt2) '1/R3AVE', (r3iave(npair),npair=1,npairs)
      write(6,fmt2) '1/R6AVE', (r6iave(npair),npair=1,npairs)
      write(6,*)
      write(6,*) '***** <P2 / (r^3*r^3)> *****'
      do nint=1,intvals
         write(6,fmt) dble(nint - 1)*tintvl, &
               (cf(intvals * (npair-1) + nint) / &
               cf(intvals * (npair-1) + 1),npair=1,npairs)
      end do
      write(6,*)
      write(6,*) '***** <P2>             *****'
      do nint=1,intvals
         write(6,fmt) dble(nint - 1)*tintvl, &
               (p2cf(intvals * (npair-1) + nint) / &
               p2cf(intvals * (npair-1) + 1),npair=1,npairs)
      end do
      write(6,*)
      write(6,*) '***** <1 / (r^3*r^3)>  *****'
      do nint=1,intvals
         write(6,fmt) dble(nint - 1)*tintvl, &
               (rcf(intvals * (npair-1) + nint) / &
               rcf(intvals * (npair-1) + 1),npair=1,npairs)
      end do
   end if  ! (ntrun == 9)
   
   100 format(' Processing snapshot: ',i8)
   
   return
end subroutine corffft 

!---------------------------------------------------------------------


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine spherharm here]
subroutine spherharm(l,m,x,y,z,r,dreal,dimg)
   
   !     Calc spherical harmonics of order l=0,1,2
   !     and -l<=m<=l with cartesian coordinates as input
   !     (see e.g. Merzbacher, Quantum Mechanics, p. 186)
   
   implicit double precision(a-h, o-z)
   
   parameter(sh00=0.28209479)
   parameter(sh10=0.48860251)
   parameter(sh11=0.34549415)
   parameter(sh20=0.31539157)
   parameter(sh21=0.77254840)
   parameter(sh22=0.38627420)
   dreal = 0.0d0
   dimg = 0.0d0
   ri = 1.0d0 / r
   
   if(l == 0 .and. m == 0) then
      dreal = sh00
   else if(l == 1) then
      if(m == 0) then
         dreal = sh10 * z * ri
      else
         dreal = -m * sh11 * x * ri
         dimg  = -    sh11 * y * ri
      end if
   else if(l == 2) then
      if(m == 0) then
         dreal = sh20 * (2.0*z*z - x*x - y*y) * ri * ri
      else if(abs(m) == 1) then
         dreal = -m * sh21 * x * z * ri * ri
         dimg  = -    sh21 * y * z * ri * ri
      else
         dreal = sh22 * (x*x - y*y) * ri * ri
         dimg  = m * sh22 * x * y * ri * ri
      end if
   end if
   
   return
end subroutine spherharm 

!---------------------------------------------------------------------

