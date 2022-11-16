
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of program g98_to_delta here]
program g98_to_delta

   !    use g98 frequencies to estimate Raman intensities
   
   !    arguments are:  <s0 pdb file>  <s1 pdb file> <freq file>
   
   implicit none
   integer maxat,nat,i,j,k,l,lmin,lmax,m1,m2,nvect
   parameter( maxat=50 )
   integer iel(3*maxat)
   character(len=80) arg
   character(len=30) spacer
   character(len=1) asym
   real x0(3*maxat), vect(3*maxat,3*maxat),freq(3*maxat)
   real x1(3*maxat), deltar(3*maxat), mass(3*maxat),mu(3*maxat)
   real delta,relax,rintens,relax_tot,sum,am
   
   !     ---get the ground state coords from the min. file.
   
   call getarg(1,arg)
   write(6,*) 'getting ground  state geometry from ',arg
   open(unit=10,file=arg,status='OLD',form='FORMATTED')
   
   do i=1,maxat
      read(10,'(30x,3f8.3)',end=98) x0(3*i-2),x0(3*i-1),x0(3*i)
   end do
   98 nat = i-1
   nvect = 3*nat - 6
   close(10)
   
   !     ---get the excited state coords from the min. file.
   
   call getarg(2,arg)
   write(6,*) 'getting excited state geometry from ',arg
   open(unit=10,file=arg,status='OLD',form='FORMATTED')
   
   do i=1,maxat
      read(10,'(12x,a1,17x3f8.3)',end=99) &
            asym,x1(3*i-2),x1(3*i-1),x1(3*i)
   end do
   99 nat = i-1
   nvect = 3*nat - 6
   close(10)
   
   !     --- get the displacement vector between ground and excited states:
   
   do k=1,3*nat
      deltar(k) = x1(k) - x0(k)
   end do
   !     write(6,'(6f10.5)') (deltar(k), k=1,3*nat)
   
   !     --- read in the frequency section from the g98 output
   !         (old style, or high-precision format)
   !         these are mass-weighted displacements
   
   call getarg(3,arg)
   write(6,*) 'getting ground state modes from ',arg
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
   
   !     ---- check orthogonality:
   
   !     do i=1,nvect
   !       do j=i,nvect
   !         sum = 0.0
   !         do k=1,3*nat
   !           sum = sum + vect(k,i)*vect(k,j)*mass(k)
   !         end do
   !         sum = sum/sqrt(mu(i)*mu(j))
   !         write(6,*) i,j,sum
   !       end do
   !     end do
   
   !     ---- compute dimensionless deltas:
   
   write(6,'(a)') '  mode    freq     DELTA   relaxation  intensity'
   relax_tot = 0.0
   do i=1,nvect
      sum = 0.0
      do k=1,3*nat
         sum = sum + vect(k,i)*deltar(k)*mass(k)
      end do
      delta = 0.171*sqrt(freq(i)/mu(i))*sum
      relax = 0.5*freq(i)*delta*delta
      relax_tot = relax_tot + relax
      rintens = (delta*freq(i))**2
      write(6,'(i5,f10.1,f10.5,2f10.1)') i,freq(i),delta,relax,rintens
   end do
   write(6,'(25x,f10.1)') relax_tot
   
end program g98_to_delta 
