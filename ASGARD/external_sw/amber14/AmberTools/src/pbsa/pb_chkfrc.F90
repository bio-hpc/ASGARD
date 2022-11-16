#include "copyright.h"
#include "../include/dprec.fh"

subroutine chkfrc(xx,ix,ih,ipairs,x,fg,ene,r_stack,i_stack, &
              atm1,atm2,ntrns,rotopt,nrotx,nsph,mrot,mrotx,msph,delta)

   use memory_module
   use poisson_boltzmann, only: savh, nfocus
   use random, only: amrset_gen, amrand_gen, rand_gen_state

   implicit none

#  include "pb_constants.h"
#  include "../include/md.h"
#  include "pb_md.h"

   ! Passed
   integer atm1, atm2, ntrns, rotopt, nrotx, nsph, mrot, mrotx, msph
   integer ix(*), ipairs(*), i_stack(*)
   character(len=4) ih(*)
   _REAL_ delta
   _REAL_ xx(*), x(*), fg(*)
   _REAL_ ene(51), r_stack(*)

   ! Local
   integer, parameter :: MAXSPH = 1000
   integer i, j, k, l, m, n, s
   integer ncrcl, mcrcl, nrot, nfopt, n_force_calls
   _REAL_ theta(MAXSPH), psi(MAXSPH)
   _REAL_ dsplc(3), rtn_c(3), rtn_s(3)
   _REAL_ crd(3*natom), xcrd(3*natom), ycrd(3*natom), zcrd(3*natom)
   _REAL_ f(3*natom), nf(3*natom)
   _REAL_ tf(3*natom), tfsq(3*natom)
   _REAL_ tnf(3*natom), tnfsq(3*natom)
   _REAL_ h, sxangle, xangle, randval, vir(4), enep, enem, rms
   _REAL_ aene(51)
   _REAL_ prob
   logical do_list_update
   type (rand_gen_state) :: rand_gen
   _REAL_ crn(3,20), crn1(3), crn2(3)

   if ( nfocus /= 1 ) then
      print *, "nfocus must be 1"
      call mexit(6,0)
   end if

   ! initialization
   prob = 1.5d0
   do_list_update = .true.
   h = savh(nfocus)
   vir = ZERO
   nf = ZERO
   tf = ZERO; tfsq = ZERO
   tnf = ZERO; tnfsq = ZERO
   n_force_calls = 0
   rms = ZERO
   crd(1:3*natom) = x(1:3*natom)
   call amrset_gen(rand_gen,ig)
   if ( rotopt == 1 ) then
      call crcl(ncrcl,nsph,theta,psi)   
      sxangle = TWOPI/nrotx
      nrot = nrotx*ncrcl
   else
      call crcl(mcrcl,msph,theta,psi)   
      sxangle = TWOPI/mrotx
      nrot = mrot
   end if
   xangle = ZERO
!  s = 0
!  k = 0

!  do l = 1, mcrcl
!     write(400,'(2f20.10)') theta(l), psi(l)
!  end do

   ! let the main axis coincide/parallel with z axis
!  m = 3*(atm1-1)
!  dsplc = -crd(m+1:m+3)
!  call translate(natom,crd,dsplc)
!  m = 3*(atm2-1)
!  n = 3*(atm1-1)
!  rtn_c = crd(m+1:m+3) - crd(n+1:n+3)
!  call sphe(rtn_c,rtn_s)
   rtn_s = ZERO
!  call rotate(natom,crd,3,-rtn_s(3))
!  call rotate(natom,crd,2,-rtn_s(2))
!  open (unit=96,file='rot.dat')
   do l = 1, nrot
      if ( rotopt == 1 ) then
         s = floor(real(l-1)/real(ncrcl))+1
         k = l-(s-1)*ncrcl
         ! rotate wrt the main axis (z axis)
         xangle = (s-1)*sxangle
         xcrd = crd
         call rotate(natom,xcrd,3,xangle)
         !rotate the main axis
         ycrd = xcrd
         call rotate(natom,ycrd,2,theta(k))
         call rotate(natom,ycrd,3,psi(k))
      else
         call amrand_gen(rand_gen,randval)
         do while ( randval == ZERO ) 
            call amrand_gen(rand_gen,randval)
         end do
         ! rotate wrt the main axis (z axis)
         randval = randval*100.d0
         s = floor(randval)
         xangle = s*sxangle
         xcrd = crd
         call rotate(natom,xcrd,3,xangle)
         !rotate the main axis
         k = ceiling((randval-s)*1000.d0)
         ycrd = xcrd
         call rotate(natom,ycrd,2,theta(k))
         call rotate(natom,ycrd,3,psi(k))
         write(96,*) s, k
      end if

      !translate the molecule
      do j = 1, ntrns
         do i = 1, 3
            call amrand_gen(rand_gen,randval)
            dsplc(i) = (randval-HALF)*h
         end do
!        dsplc = ZERO
         zcrd = ycrd
         call translate(natom,zcrd,dsplc)
         x(1:3*natom) = zcrd(1:3*natom)
         fg(1:3*natom) = ZERO
            
!        write(300,*) '--------------------------------'
!        write(300,*) 'rotation',s,k,'translation',j
!        write(300,*) natom
!        do n = 1, natom-1, 2
!           m = 3*(n-1)
!           write(300,'(6f12.7)') zcrd(m+1),zcrd(m+2),zcrd(m+3),zcrd(m+4),zcrd(m+5),zcrd(m+6)
!        end do
!        if ( natom == n ) then
!           m = 3*(n-1)
!           write(300,'(3f12.7)') zcrd(m+1),zcrd(m+2),zcrd(m+3)
!        end if
!        cycle
      
!        do n = 1, natom
!           m = 3*(n-1)
!           write(301,'(3f20.15)') zcrd(m+1),zcrd(m+2),zcrd(m+3)
!        end do

         ! analytical force
         ene = ZERO
         pbgrid = .true.
         call force(xx,ix,ih,ipairs,x,fg,ene,vir,r_stack,i_stack, &
                    xx(l96),xx(l97),xx(l98),do_list_update)
         aene = ene

!        do m = 1, natom
!           n=(m-1)*3
!           write(555,'(3f20.10)') fg(n+1),fg(n+2),fg(n+3)
!        end do

         ! numerical force
         nfopt = 0
         if ( nfopt == 1 ) then
            ! only for cg test
            do n = 1, 3
!              print *, l, k, j, n
               x(1:3*natom) = zcrd(1:3*natom)
               do m = 20, 35
                  i = (m-1)*3+n
                  x(i) = x(i) + delta
               end do
               ene = ZERO
               pbgrid = .true.
               call force(xx,ix,ih,ipairs,x,f,ene(23),vir,r_stack,i_stack, &
                          xx(l96),xx(l97),xx(l98),do_list_update)
               enep = ene(23)
               x(1:3*natom) = zcrd(1:3*natom)
               do m = 20, 35
                  i = (m-1)*3+n
                  x(i) = x(i) - delta
               end do
               ene = ZERO
               pbgrid = .true.
               call force(xx,ix,ih,ipairs,x,f,ene(23),vir,r_stack,i_stack, &
                          xx(l96),xx(l97),xx(l98),do_list_update)
               enem = ene(23)
               nf(n) = -(enep-enem)/(TWO*delta) 
            end do
         elseif ( nfopt == 2 ) then 
            do m = 1, natom
               do n = 1, 3
                  i = (m-1)*3+n
                  x(1:3*natom) = zcrd(1:3*natom)
                  x(i) = x(i) + delta
                  ene = ZERO
                  pbgrid = .true.
                  call force(xx,ix,ih,ipairs,x,f,ene(23),vir,r_stack,i_stack, &
                             xx(l96),xx(l97),xx(l98),do_list_update)
                  enep = ene(23)
                  x(1:3*natom) = zcrd(1:3*natom)
                  x(i) = x(i) - delta
                  ene = ZERO
                  pbgrid = .true.
                  call force(xx,ix,ih,ipairs,x,f,ene(23),vir,r_stack,i_stack, &
                             xx(l96),xx(l97),xx(l98),do_list_update)
                  enem = ene(23)
                  nf(i) = -(enep-enem)/(TWO*delta) 
               end do
            end do
         else 
         end if

         ! rotate back force
         call rotate(natom,fg,3,-psi(k))
         call rotate(natom,fg,2,-theta(k))
         call rotate(natom,fg,3,-xangle)
!        call rotate(natom,fg,2,rtn_s(2))
!        call rotate(natom,fg,3,rtn_s(3))
         if ( nfopt /= 0 ) then
            call rotate(natom,nf,3,-psi(k))
            call rotate(natom,nf,2,-theta(k))
            call rotate(natom,nf,3,-xangle)
!           call rotate(natom,nf,2,rtn_s(2))
!           call rotate(natom,nf,3,rtn_s(3))
         end if

         ! print out force
!        write(6,'(3f20.10)') dsplc(1),dsplc(2),dsplc(3)
!        write(6,*) zcrd(1:3*natom)
         call writefile(natom,s,k,j,fg,nf, &
                        tf,tfsq,tnf,tnfsq, &
                        nfopt)

         ! print out energy
         n_force_calls = n_force_calls + 1
         call report_min_progress( n_force_calls, rms, fg, aene, ih(m04) )  
      end do
   end do

   call writesummary(natom,nrot,ntrns, &
                     tf,tfsq,tnf,tnfsq, &
                     nfopt)

end subroutine chkfrc

subroutine sphe(car,sph)

   implicit none

#  include "pb_constants.h"

   _REAL_ car(3), sph(3)

   sph(1) = sqrt(car(1)**2+car(2)**2+car(3)**2)
   sph(2) = acos(car(3)/sph(1))
   if ( car(1) == 0 ) then
      sph(3) = HALF*PI
   else
      sph(3) = atan(abs(car(2)/car(1)))
   end if

   if ( car(1) >= 0 ) then
      if ( car(2) >= 0 ) then
         sph(3) = sph(3)
      else
         sph(3) = TWOPI - sph(3)
      end if
   else
      if ( car(2) >= 0 ) then
         sph(3) = PI - sph(3)
      else
         sph(3) = PI + sph(3)
      end if
   end if

end subroutine sphe

subroutine translate(n,car,dsplc)

   implicit none

   integer n
   _REAL_ car(3,n), dsplc(3)

   car(1,1:n) = car(1,1:n) + dsplc(1)
   car(2,1:n) = car(2,1:n) + dsplc(2)
   car(3,1:n) = car(3,1:n) + dsplc(3)

end subroutine translate

subroutine rotate(n,car,axis,angle)

   implicit none

#  include "pb_constants.h"

   integer i, n, axis
   _REAL_ car(3,n), angle
   _REAL_ R(3,3), tmp(3)

   ! Warning eliminator
   R = ZERO

    if ( axis == 1 ) then
       R(1,1:3) = (/ONE,ZERO,ZERO/)
       R(2,1:3) = (/ZERO,cos(angle),-sin(angle)/)
       R(3,1:3) = (/ZERO,sin(angle),cos(angle)/)
    else if ( axis == 2 ) then
       R(1,1:3) = (/cos(angle),ZERO,sin(angle)/)
       R(2,1:3) = (/ZERO,ONE,ZERO/)
       R(3,1:3) = (/-sin(angle),ZERO,cos(angle)/)
    else if ( axis == 3 ) then
       R(1,1:3) = (/cos(angle),-sin(angle),ZERO/)
       R(2,1:3) = (/sin(angle),cos(angle),ZERO/)
       R(3,1:3) = (/ZERO,ZERO,ONE/)
    else
       print *, "rotation axis error"
       call mexit(6,0)
    end if

    do i = 1, n
       tmp = car(1:3,i)
       car(1,i) = R(1,1)*tmp(1)+R(1,2)*tmp(2)+R(1,3)*tmp(3)
       car(2,i) = R(2,1)*tmp(1)+R(2,2)*tmp(2)+R(2,3)*tmp(3)
       car(3,i) = R(3,1)*tmp(1)+R(3,2)*tmp(2)+R(3,3)*tmp(3)
    end do

end subroutine rotate

subroutine crcl(n,maxsph,theta,psi)

   implicit none

#  include "pb_constants.h"

   ! Passed variables

   integer n, maxsph
   _REAL_ theta(maxsph), psi(maxsph)

   ! Local variables

   integer ntheta,npsi,npsimax,i,nt,np
   _REAL_ thtstp,tht,stheta,psistp

   ! begin code

   ntheta = int(sqrt(PI*maxsph/FOUR))
   npsimax = int(TWO*ntheta)
   thtstp = PI/ntheta

   i = 1
   do nt = 1, ntheta
      tht = thtstp*nt
      stheta = sin(tht)
      npsi = nint(stheta*npsimax)
      if (npsi == 0) cycle
      psistp = TWOPI/npsi
      do np = 1, npsi
         psi(i) = np*psistp
         theta(i) = tht
         i = i + 1
      enddo
   enddo
   n = i - 1

end subroutine crcl

subroutine writefile(natom,l,k,j,fg,nf, &
                     tf,tfsq,tnf,tnfsq, &
                     nfopt)

   implicit none

   integer l, k, j
   integer m, natom, nfopt
   _REAL_ fg(3,natom), nf(3,natom)
   _REAL_ tf(3,natom), tfsq(3,natom)
   _REAL_ tnf(3,natom), tnfsq(3,natom)
   _REAL_ tf1(3), tnf1(3)

   tf1(1) = sum(fg(1,1:natom))
   tf1(2) = sum(fg(2,1:natom))
   tf1(3) = sum(fg(3,1:natom))
   tnf1(1) = sum(nf(1,1:natom))
   tnf1(2) = sum(nf(2,1:natom))
   tnf1(3) = sum(nf(3,1:natom))

   tf = tf + fg
   tfsq = tfsq+fg*fg
   tnf = tnf + nf 
   tnfsq = tnfsq+nf*nf
   
   write(6,'(a)') '--------------------------------'
   write(6,'(a,i12,i12,a,i12)') 'rotation',l,k,' translation',j
   write(6,'(a)') 'analytical force'
   do m = 1, natom
      write(6,'(3f20.10)') fg(1,m),fg(2,m),fg(3,m)
   end do
   write(6,'(a)') 'total analytical force'
   write(6,'(3f20.10)') tf1(1),tf1(2),tf1(3)

   if ( nfopt /= 0 ) then
      write(6,'(a)') 'numerical force'
      do m = 1, natom
         write(6,'(3f20.10)') nf(1,m),nf(2,m),nf(3,m)
      end do
      write(6,'(a)') 'total numerical force'
      write(6,'(3f20.10)') tnf1(1),tnf1(2),tnf1(3)
   end if

end subroutine writefile

subroutine writesummary(natom,nrot,ntrns, &
                        tf,tfsq,tnf,tnfsq, &
                        nfopt)
   implicit none

   integer natom,nrot,ntrns
   _REAL_ tf(3,natom), tfsq(3,natom)
   _REAL_ tnf(3,natom), tnfsq(3,natom)
   integer nfopt

   integer m, nort
   _REAL_ wt

   nort = nrot*ntrns
   write(6,'(a)') "Summary of finite-difference force"
   write(6,'(a,i12)') "Number of orientations:  ",nort

   wt = 1.d0/dble(nort)
   tf = tf*wt
   tfsq = sqrt(tfsq*wt-tf*tf)
   write(6,'(a)') "Mean of analytical force"
   do m = 1, natom
      write(6,'(3f20.10)') tf(1,m),tf(2,m),tf(3,m)
   end do
   write(6,'(a)') "Standard deviation of analytical force"
   do m = 1, natom
      write(6,'(3f20.10)') tfsq(1,m),tfsq(2,m),tfsq(3,m)
   end do
   
   if ( nfopt /= 0 ) then
      tnf = tnf*wt
      tnfsq = sqrt(tnfsq*wt-tnf*tnf)
      write(6,'(a)') "Mean of numerical force"
      do m = 1, natom
         write(6,'(3f20.10)') tnf(1,m),tnf(2,m),tnf(3,m)
      end do
      write(6,'(a)') "Standard deviation of numerical force"
      do m = 1, natom
         write(6,'(3f20.10)') tnfsq(1,m),tnfsq(2,m),tnfsq(3,m)
      end do
   end if

end subroutine writesummary

subroutine get_arccrn(natom,crd,r,idx,prob,crn1,crn2)

   implicit none

#  include "pb_constants.h"

   ! Passed
   integer natom, idx(3)
   _REAL_ crd(3,natom) 
   _REAL_ prob, r(natom)
   _REAL_ crn1(3), crn2(3)

   ! Local
   integer m, n, mp, np, j
   _REAL_ a1(3), a2(3), a3(3)
   _REAL_ p1(3), p2(3), q1(3), q2(3), p(3), d(3)
   _REAL_ u(3,3), v(3,3), w(3), b(3), t(3)
   _REAL_ TOL, wmax, thresh, h
 
   a1 = crd(1:3,idx(1)); a2 = crd(1:3,idx(2)); a3 = crd(1:3,idx(3))
   p1 = a2 - a1; p2 = a3 - a1
   q1 = (a1+a2)*HALF; q2 = (a1+a3)*HALF
   d(1) = r(idx(1)) + prob
   d(2) = r(idx(2)) + prob
   d(3) = r(idx(3)) + prob
   p(1) = p1(2)*p2(3) - p2(2)*p1(3)
   p(2) = p2(1)*p1(3) - p1(1)*p2(3)
   p(3) = p1(1)*p2(2) - p2(1)*p1(2)
   b(1) = dot_product(p1,q1) + (d(1)+d(2))*(d(1)-d(2))*HALF
   b(2) = dot_product(p2,q2) + (d(1)+d(3))*(d(1)-d(3))*HALF
   b(3) = dot_product(p,a1)
   u(1,1:3) = p1
   u(2,1:3) = p2
   u(3,1:3) = p
!  print *, a1(1), a1(2), a1(3)
!  print *, a2(1), a2(2), a2(3)
!  print *, a3(1), a3(2), a3(3)
!  print *, p1(1), p1(2), p1(3)
!  print *, d(1), d(2), d(3)
!  print *, "u1",u(1,1), u(2,1), u(3,1)
!  print *, "u2",u(1,2), u(2,2), u(3,2)
!  print *, "u3",u(1,3), u(2,3), u(3,3)
!  print *, "b",b(1), b(2), b(3)
  
   mp = 3; np = 3
   m = 3; n = 3
   TOL = 1.d-5
   call svdcmp(u,m,n,mp,np,w,v)
   wmax = ZERO
   do j = 1, n
     if ( w(j) > wmax ) wmax = w(j)
   end do
   thresh = TOL * wmax
   do j = 1, n
     if ( w(j) < thresh ) w(j) = ZERO
   end do
   call svbksb(u,w,v,m,n,mp,np,b,t)
!  print *, "t",t(1), t(2), t(3)

   p = p/sqrt(p(1)**2+p(2)**2+p(3)**2)
   h = sqrt(d(1)**2-dot_product(t-a1,t-a1))
!  print *,d(1),h
!  print *,t(1)-a1(1),t(2)-a1(2),t(3)-a1(3)
!  h = sqrt(d(2)**2-dot_product(t-a2,t-a2))
!  print *,h
!  h = sqrt(d(3)**2-dot_product(t-a3,t-a3))
!  print *,h
   crn1 = t + h*p
   crn2 = t - h*p
!  print *, "c1",crn1(1), crn1(2), crn1(3)
!  print *, "c2",crn2(1), crn2(2), crn2(3)

end subroutine get_arccrn
