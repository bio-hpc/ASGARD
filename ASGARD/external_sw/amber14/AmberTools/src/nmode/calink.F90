
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of program calink here]
program calink
   
   !     Use the "Bahar" model to compute normal modes from CA-positions
   
   !     Usage:  calink -nmode # [-cut <cutoff>] [-novec] [-v vecs]
   
   implicit double precision (a-h,o-z)
   integer nargs,iarg,ftyp,natom
   character(len=80) arg,xyz,vecs,line
   logical novec,hasbox
   real crdr,eig4
   
#  include "calink.h"
   
   dimension crdr(natm3)
   integer iparam(11),ipntr(11)
   logical select(natm3)
   common /q/ x(natm3),workl(lworkl),dr(3), &
         workd(3*natm3),eigval(natm3),vout(natm3,ncvmax),resid(natm3), &
         ad(3,3,natmax),an(3,3,ijmax),ia(natmax+1),ja(ijmax)
   
   call wallclock(time0)
   xyz='xyz'
   vecs='vecs'
   ftyp = 5
   novec = .false.
   last = 999999
   nmode = 999999
   hasbox = .false.
   cut = 7.d0
   
   !     read command line arguments
   
#ifdef HITACHI_GETARG
   iarg = 1
#else
   iarg = 0
#endif
   nargs=iargc()
   if (nargs == iarg) goto 20
   10 continue
   iarg=iarg+1
   if (iarg > nargs) goto 20
   call getarg(iarg,arg)
   if (arg == '-novec') then
      novec = .true.
   else if (arg == '-v') then
      iarg=iarg+1
      call getarg(iarg,vecs)
   else if (arg == '-nmode') then
      iarg=iarg+1
      call getarg(iarg,arg)
      read(arg,*) nmode
   else if (arg == '-cut') then
      iarg=iarg+1
      call getarg(iarg,arg)
      read(arg,*) cut
   else
      if (arg(1:1) == '<' .or. arg(1:1) == '>' .or. &
            arg(1:1) == '|') goto 20
      write(0,'(/,5x,a,a)') 'unknown flag: ',arg
      call usage
   end if
   goto 10

   20 continue
   cut2 = cut*cut
   
   !  --- check number of modes
   
   if (nmode > nmodemax) then
      write(0,*) 'Too many modes specified:',nmode,nmodemax
      call usage
   end if
   
   ! --- read coords in PDB format:
   
   natom = 0
   do iline=1,999999
      read(5,'(a80)',end=25 ) line
      if( line(1:4) == 'ATOM' ) then
         natom = natom + 1
         if (natom > natmax) then
            write(0,*) 'Too many atoms specified:',natom,natmax
            call usage
         end if
         !         read( line,'(30x,3f8.3)' ) x(3*natom-2),x(3*natom-1),
         !    .                               x(3*natom)
         read( line,'(30x,3f13.8)' ) x(3*natom-2),x(3*natom-1), &
               x(3*natom)
      end if
   end do
   25 n = 3*natom
   
   !  --- construct the Hessian in a sparse matrix format:
   
   
   !       zero out all diagonal (i,i) subblocks
   
   do i=1,natom
      do k=1,3
         do l=1,3
            ad(k,l,i) = 0.d0
         end do
      end do
   end do
   
   !       loop over all pairs of atoms:
   
   ij = 1
   do i=1,natom-1
      xi = x(3*i-2)
      yi = x(3*i-1)
      zi = x(3*i  )
      ia(i) = ij
      
      do j=i+1,natom
         dx = x(3*j-2) -xi
         dy = x(3*j-1) -yi
         dz = x(3*j  ) -zi
         r2 = dx*dx + dy*dy + dz*dz
         if( r2 < cut2 ) then
            bi = 1.d0/sqrt(r2)
            dr(1) = dx*bi
            dr(2) = dy*bi
            dr(3) = dz*bi
            
            !             update the (i,i), (j,j) and (i,j) sub-blocks
            
            do k=1,3
               do l=1,3
                  ad(k,l,i) = ad(k,l,i) + dr(k)*dr(l)
                  ad(k,l,j) = ad(k,l,j) + dr(k)*dr(l)
               end do
            end do
            do k=1,3
               do l=1,3
                  an(k,l,ij) = -dr(k)*dr(l)
               end do
            end do
            ja(ij) = j
            ij = ij + 1
            if( ij > ijmax ) then
               write(6,*) 'ij is greater than IJMAX!'
               call usage
            end if
         end if
         
      end do  !  j=i+1,natom
   end do  !  i=1,natom-1
   ia(natom) = ij
   ia(natom+1) = ij
   
   call wallclock(time1)
   write( 6,'(a,f8.2)') '| Time to construct A: ',time1-time0
   time0 = time1
   write(6,'(a,i8,a,f8.2)') 'number of pairs: ', ij, &
         ' for cutoff ',cut
#if 0
   do i=1,natom
      write(6,*) 'diagonal for i=',i
      write(6,*) ad(1,1,i),ad(1,2,i),ad(1,3,i)
      write(6,*) ad(2,1,i),ad(2,2,i),ad(2,3,i)
      write(6,*) ad(3,1,i),ad(3,2,i),ad(3,3,i)
   end do
   do i=1,ij-1
      write(6,*) 'off-diagonal:', ja(i)
      write(6,*) an(1,1,i),an(1,2,i),an(1,3,i)
      write(6,*) an(2,1,i),an(2,2,i),an(2,3,i)
      write(6,*) an(3,1,i),an(3,2,i),an(3,3,i)
   end do
#endif
   
   !  --- diagonalize the matrix:
   
   tol = dlamch('EPS')
   inv_iter = 0
   ido = 0
   nev = nmode
   
   !     best value of ncv is hard to guess; the original idea,
   !     (ncv = 2*nev), seems not to work very well in some instances.
   !     Here we try the larger of 100 and 2*nev+10:
   
   ncv = max( 100, 2*nev+10 )
   info = 0
   iparam(1) = 1
   iparam(3) = 300
   iparam(4) = 1
   iparam(7) = 1
   !     tol = 0.0
   tol = dlamch('EPS')
   write(6,*) 'using tolerance of ',tol
   
   34 continue
   call dsaupd( ido, 'I', n, 'SA', nev, tol, resid, &
         ncv, vout, natm3, iparam, ipntr, workd, &
         workl, lworkl, info )
   if( ido == -1 .or. ido == 1 ) then
      
      !       ----construct y = Ax:
      
      inv_iter = inv_iter + 1
      if( mod(inv_iter,100) == 0 ) &
            write(0,*) 'inverse iteration ',inv_iter
      call matvec(natom,ia,ja,an,ad,workd(ipntr(1)),workd(ipntr(2)))
      
      
      goto 34
   end if
   write(6,'(a,i4,a,i6,a)') 'dsaupd returns: ',info,' after ', &
         inv_iter,' iterations'
   write(6,*) 'number of converged eigenvalues is ',iparam(5)
   call wallclock(time1)
   write( 6,'(a,f8.2)') '| Time for dsaupd : ',time1-time0
   time0 = time1
   if( info < 0 ) call mexit(6,1)
   
   call dseupd( .true., 'A', select, eigval, vout, natm3, sigma, &
         'I', n, 'SA', nev, tol, resid, &
         ncv, vout, natm3, iparam, ipntr, workd, &
         workl, lworkl, info )
   write(6,'(a,i4)') 'dseupd returns: ',info
   call wallclock(time1)
   write( 6,'(a,f8.2)') '| Time for dseupd : ',time1-time0
   time0 = time1
   
   !  --- output the eigenvectors and eigenvalues:
   
#ifdef FORMATTED
   open(unit=61,file=vecs,status='NEW',form='FORMATTED',iostat=ios)
#else
   open(unit=61,file=vecs,status='NEW',form='UNFORMATTED',iostat=ios)
#endif
   if (ios /= 0) then
      close(unit=61)
      write(0,'(/,2x,a,a)') 'Error on OPEN: ',vecs
      call mexit(6,1)
   end if
   
#ifdef FORMATTED
   write(61,'(a,f8.2)') 'eigenvector file, cutoff was ',cut
   write(61,'(i5)') n
   write(61,'(7f11.5)') (x(i),i=1,n)
#else
   write(61) n
   do i=1,n
      crdr(i) = x(i)
   end do
   write(61) (crdr(i),i=1,n)
#endif
   do j = 1,nmode
      !       eigval(j) = sqrt(abs(eigval(j)))
      
      !       ---calculate the collectivity index
      
      m3 = 0
      sum = 0.d0
      do m=1,natom
         u2 = vout(m3+1,j)*vout(m3+1,j) + vout(m3+2,j)*vout(m3+2,j) &
               + vout(m3+3,j)*vout(m3+3,j)
         sum = sum + u2*log(u2)
         m3 = m3 + 3
      end do
      coll = exp(-sum)/natom
      write(6,'(i5,f15.10,f10.5)') j, eigval(j), coll
      
#ifdef FORMATTED
      write (61,'(1x,a)') '****'
      write (61,'(i5,f12.7)') j, eigval(j)
#else
      eig4 = eigval(j)
      write(61) j,eig4
#endif
      do i=1,n
#ifdef FORMATTED
         vout(i,j) = vout(i,j)
#else
         crdr(i) = vout(i,j)
#endif
      end do
#ifdef FORMATTED
      write (61,'(7f11.5)') (vout(i,j),i=1,n)
#else
      write(61) (crdr(i),i=1,n)
#endif
   end do
   close(unit=61)
   50 continue
   
end program calink 
!-------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine usage here]
subroutine usage()
#include "calink.h"
   write(6,'(a,a)') &
         'usage: calink [-nmode #] [-novec | -v <vecs>] [-cut #]'
   write(6,'(a)') 'Compile-in maximum sizes:'
   write(6,'(a,i6)') '      natom: ',natmax
   write(6,'(a,i6)') '      modes: ',nmodemax
   write(6,'(a,i6)') '      ij-s : ',ijmax
   call mexit(6,1)
end subroutine usage 
!-------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine matvec here]
subroutine matvec(n,ia,ja,an,ad,b,c)
   
   implicit double precision (a-h,o-z)
   
   dimension ia(*),ja(*),ad(3,3,n),an(3,3,*),b(*),c(*)
   
   !     translate algorithm 7.18 of Sparse Matrix Algebra for use with 3x3
   !       matrices as elements
   
   i3 = 0
   do i=1,n
      c(i3+1) = ad(1,1,i)*b(i3+1)+ad(1,2,i)*b(i3+2)+ad(1,3,i)*b(i3+3)
      c(i3+2) = ad(2,1,i)*b(i3+1)+ad(2,2,i)*b(i3+2)+ad(2,3,i)*b(i3+3)
      c(i3+3) = ad(3,1,i)*b(i3+1)+ad(3,2,i)*b(i3+2)+ad(3,3,i)*b(i3+3)
      i3 = i3 + 3
   end do
   i3 = 0
   do i=1,n
      iaa = ia(i)
      iab = ia(i+1) - 1
      if (iab < iaa ) goto 30
      u1 = c(i3+1)
      u2 = c(i3+2)
      u3 = c(i3+3)
      z1 = b(i3+1)
      z2 = b(i3+2)
      z3 = b(i3+3)
      do k=iaa,iab
         j = ja(k)
         j3 = 3*(j-1)
         u1 = u1 + an(1,1,k)*b(j3+1) + an(1,2,k)*b(j3+2) &
               + an(1,3,k)*b(j3+3)
         u2 = u2 + an(2,1,k)*b(j3+1) + an(2,2,k)*b(j3+2) &
               + an(2,3,k)*b(j3+3)
         u3 = u3 + an(3,1,k)*b(j3+1) + an(3,2,k)*b(j3+2) &
               + an(3,3,k)*b(j3+3)
         c(j3+1) = c(j3+1) + an(1,1,k)*z1 + an(1,2,k)*z2 + an(1,3,k)*z3
         c(j3+2) = c(j3+2) + an(2,1,k)*z1 + an(2,2,k)*z2 + an(2,3,k)*z3
         c(j3+3) = c(j3+3) + an(3,1,k)*z1 + an(3,2,k)*z2 + an(3,3,k)*z3
      end do
      
      c(i3+1) = u1
      c(i3+2) = u2
      c(i3+3) = u3
      30 i3 = i3 + 3
   end do
   return
end subroutine matvec 
