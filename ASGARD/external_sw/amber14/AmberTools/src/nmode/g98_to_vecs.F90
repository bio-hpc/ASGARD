
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of program g98_to_vecs here]
program g98_to_vecs

   !    convert g98 frequencies to Amber "vecs" format
   
   !    arguments are:  <g98 min file>  <g98 freq file>
   
   implicit none
   integer maxat,nat,i,j,k,l,lmin,lmax,m1,m2,nvect
   parameter( maxat=50 )
   integer iel(3*maxat)
   character(len=80) arg
   character(len=30) spacer
   real x(3*maxat), vect(3*maxat,3*maxat),freq(3*maxat)
   real mass(3*maxat),mu(3*maxat),am
   
   !     ---get the coords from the min. file.
   
   call getarg(1,arg)
   open(unit=10,file=arg,status='OLD',form='FORMATTED')
   
   do i=1,maxat
      read(10,*,end=99) m1,m2,x(3*i-2),x(3*i-1),x(3*i)
   end do
   99 nat = i-1
   nvect = 3*nat - 6
   close(10)
   
   !     --- read in the frequency section from the g98 output
   !         (old style, or high-precision format)
   
   call getarg(2,arg)
   open(unit=11,file=arg,status='OLD',form='FORMATTED')
   
   do lmin=1,nvect,5
      read(11, '(a30)' ) spacer
      read(11, '(a30)' ) spacer
      lmax = min(nvect, lmin+4)
      read(11, '(23x,5f10.4)') (freq(l), l=lmin,lmax)
      read(11, '(23x,5f10.4)') (mu(l), l=lmin,lmax)
      read(11, '(a30)' ) spacer
      read(11, '(a30)' ) spacer
      read(11, '(a30)' ) spacer
      read(11, '(a30)' ) spacer
      read(11, '(a30)' ) spacer
      do k=1,3*nat
         read(11, '(10x,i6,7x,5f10.5)') iel(k),(vect(k,l), l=lmin,lmax)
      end do
   end do
   close(11)
   
   !     --- get atomic masses:
   
   do k=1,3*nat
      if( iel(k) == 1 ) then
         am = 1.008
      else if( iel(k) == 6 ) then
         am = 12.01
      else if( iel(k) == 7 ) then
         am = 14.01
      else if( iel(k) == 8 ) then
         am = 16.00
      else if( iel(k) == 16 ) then
         am = 32.06
      else
         write(6,*) 'unknown element: ',k,iel(k)
         stop
      end if
      mass(k) = am
   end do
   
   !     --- write out the "vecs" file:
   
   write(6,'(a,1x,a)') 'NORMAL COORDINATE FILE', arg(1:40)
   write(6,'(i5)') 3*nat
   write(6,'(7f11.5)') (x(i),i=1,3*nat)
   
   do i=1,nvect
      write(6,'(a)') ' ****'
      write(6,'(i5,f12.5)') i,freq(i)
      write(6,'(7f11.5)') (vect(j,i)/sqrt(mu(i)),j=1,3*nat)
   end do
   
end program g98_to_vecs 
