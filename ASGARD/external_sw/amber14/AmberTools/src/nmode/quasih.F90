
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of program quasih here]
program quasih
#define FORMATTED 1
   
   !     read cartesian coords from trajectory files, and compute
   !     quasiharmonic modes
   
   !     usage: quasih -natom # [-f # -m <mass> -x <xyz>} [-novec | -v <vecs>]
   !                 [ -first #] [ -last # ] [ -nmode # ]
   
   !     input file types (for -f argument; default is "5"):
   
   !     1 - .crd2 file
   !     2 - .crd  file
   !     3 - traj  file
   !     4 - .xyz file
   !     5 - binpos file
   !     6 - .crd2 file without velocities
   
   implicit double precision (a-h,o-z)
   integer nargs,iarg,ftyp,natom,first
   character(len=80) arg,mass,xyz,vecs
   logical novec,hasbox
   
   !    maximum number of atoms, snapshots: modes:
   
   parameter (natmax=2500,natm3=3*natmax)
   parameter (nspmax=1500)
   parameter (nmodemax=1000,ncvmax=2*nmodemax, &
         lworkl=ncvmax*(ncvmax+8))
   
   common /q/ at(natm3,nspmax),ener(36),xav(natm3),workl(lworkl), &
         w(nspmax),atmass(natmax),atm(natm3),atmi(natm3), &
         workd(3*natm3),eigval(natm3),vout(natm3,ncvmax),resid(natm3)
   
   !    ---define BIGMEM if you can afford the space for the transpose
   !       of the trajectory, since this might speed up matrix-multiplies
   !       below.  However, on some machines, allocating the extra memory
   !       actually costs in performance....
   
#ifdef BIGMEM
   dimension a(nspmax,natm3)
#endif
   real*4 crdr(natm3),eig4
   integer iparam(11),ipntr(11)
   logical select(natm3)
   
   call wallclock(time0)
   print *,' '
   print *,'------------------------------------------------------'
   print *,'Amber 8  Quasih                      Scripps/UCSF 2004'
   print *,'------------------------------------------------------'
   print *,' '
   mass='mass'
   xyz='xyz'
   vecs='vecs'
   ftyp = 5
   natom = 999999
   novec = .false.
   first = 1
   last = 999999
   nmode = 999999
   hasbox = .false.
   
   !     read command line arguments
   
#ifdef HITACHI_GETARG
   iarg = 1
#else
   iarg = 0
#endif
   nargs=iargc()
   if (nargs == iarg) call usage
   10 continue
   iarg=iarg+1
   if (iarg > nargs) goto 20
   call getarg(iarg,arg)
   if (arg == '-m') then
      iarg=iarg+1
      call getarg(iarg,mass)
   else if (arg == '-novec') then
      novec = .true.
   else if (arg == '-x') then
      iarg=iarg+1
      call getarg(iarg,xyz)
   else if (arg == '-v') then
      iarg=iarg+1
      call getarg(iarg,vecs)
   else if (arg == '-natom') then
      iarg=iarg+1
      call getarg(iarg,arg)
      read(arg,*) natom
   else if (arg == '-nmode') then
      iarg=iarg+1
      call getarg(iarg,arg)
      read(arg,*) nmode
   else if (arg == '-f') then
      iarg=iarg+1
      call getarg(iarg,arg)
      read(arg,*) ftyp
   else if (arg == '-first') then
      iarg=iarg+1
      call getarg(iarg,arg)
      read(arg,*) first
   else if (arg == '-last') then
      iarg=iarg+1
      call getarg(iarg,arg)
      read(arg,*) last
   else
      if (arg(1:1) == '<' .or. arg(1:1) == '>' .or. &
            arg(1:1) == '|') goto 20
      write(0,'(/,5x,a,a)') 'unknown flag: ',arg
      call usage
   end if  ! (arg == '-m')
   goto 10

   20 continue
   
   !  --- check number of atoms
   
   if (natom > natmax) then
      write(0,*) 'Too many atoms specified:',natom,natmax
      call usage
   end if
   n = 3*natom
   print *,'quasih, natom = ',natom
   
   !  --- check number of modes
   
   if (nmode > nmodemax) then
      write(0,*) 'Too many modes specified:',nmode,nmodemax
      call usage
   end if
   print *,'quasih, nmode = ',nmode
   
   !  ---   read the mass file
   
   open(unit=10,file=mass,status='OLD',form='FORMATTED',iostat=ios)
   if (ios /= 0) then
      close(unit=10)
      write(0,'(/,2x,a,a)') 'Error on OPEN: ',mass
      call mexit(6,1)
   end if
   read(10,'(10f8.2)') (atmass(i),i=1,natom)
   close(unit=10)
   k = 0
   do i=1,natom
      atm(k+1) = sqrt(atmass(i))
      atm(k+2) = atm(k+1)
      atm(k+3) = atm(k+1)
      atmi(k+1) = 1./atm(k+1)
      atmi(k+2) = 1./atm(k+1)
      atmi(k+3) = 1./atm(k+1)
      k = k + 3
   end do
   
   !  ---   open the input file
   
   call openinp(ftyp)
   
   !  ---   read input file
   
   ierr=0
   ieof=0
   do i=1,n
      xav(i) = 0.0
   end do
   
   !  --- first pass through: get the average coordinates
   
   isnap = 0
   do isp=1,last
      call readfile(natom,nread,ener,crdr,workl,box,ftyp, &
            hasbox,ierr,ieof)
      if (ierr == 1) then
         write(0,*) 'Error during read: natomal=',natom,' nread=',nread
         call mexit(6,1)
      end if
      if (ieof == 1) goto 30
      if (isp >= first) then
         isnap = isnap + 1
         if( isnap > nspmax ) then
            write(6,*) 'More than ',nspmax,' (allowed) snapshots '
            call mexit(6,1)
         end if
         do i=1,n
            at(i,isnap) = crdr(i)
            xav(i) = xav(i) + at(i,isnap)
         end do
      end if
   end do


   
   30 nsp = isnap
   write(6,'(a,i6,a,i6,a,i6)') 'Using ',nsp,' coordinate sets from', &
         first,' through ',last
   fnorm = 1./float(nsp)
   call dscal(n, fnorm, xav, 1)
   open(unit=10,file=xyz,status='NEW',form='FORMATTED',iostat=ios)
   if (ios /= 0) then
      close(unit=10)
      write(0,'(/,2x,a,a)') 'Error on OPEN: ',xyz
      call mexit(6,1)
   end if
   write(10,*) 'mean coordinates for quasiharmonic analysis'
   write(10,'(i5)') natom
   write(10,'(6f12.7)') (xav(i),i=1,3*natom)
   close (10)
   
   !  --- construct the mass-weighted deviations from the mean structure:
   
   fnorm2 = sqrt(fnorm)
   do isp=1,nsp
      do j =1,n
         at(j,isp) = fnorm2*atm(j)*(at(j,isp) - xav(j))
#ifdef BIGMEM
         a(isp,j) = at(j,isp)
#endif
      end do
   end do
   call wallclock(time1)
   write( 6,'(a,f8.2)') '| Time to construct A: ',time1-time0
   time0 = time1
   
   !  --- diagonalize the matrix:
   
   ido = 0
   nev = min( nmode, n )
   ncv = min( 2*nev, n )
   info = 0
   iparam(1) = 1
   iparam(3) = 300
   iparam(4) = 1
   iparam(7) = 1
   tol = 0.0
   
   34 continue
   call dsaupd( ido, 'I', n, 'LA', nev, tol, resid, &
         ncv, vout, natm3, iparam, ipntr, workd, &
         workl, lworkl, info )
   if( ido == -1 .or. ido == 1 ) then
      
      id = ipntr(2) - 1
      
      !       ----construct w = Ax:
      
      do j=1,nsp
         w(j) = ddot(n,at(1,j),1,workd(ipntr(1)),1)
      end do
      
      !       ---construct y = A'w:
      
      do i=1,n
#ifdef BIGMEM
         workd(id+i) = ddot(nsp,a(1,i),1,w,1)
#else
         workd(id+i) = 0.d0
         do j=1,nsp
            workd(id+i) = workd(id+i) + at(i,j)*w(j)
         end do
#endif
      end do
      
      goto 34
   end if
   write(6,'(a,i4)') 'dsaupd returns: ',info
   call wallclock(time1)
   write( 6,'(a,f8.2)') '| Time for dsaupd : ',time1-time0
   time0 = time1
   if( info < 0 ) call mexit(6,1)
   
   call dseupd( .true., 'A', select, eigval, vout, natm3, sigma, &
         'I', n, 'LA', nev, tol, resid, &
         ncv, vout, natm3, iparam, ipntr, workd, &
         workl, lworkl, info )
   write(6,'(a,i4)') 'dseupd returns: ',info
   call wallclock(time1)
   write( 6,'(a,f8.2)') '| Time for dseupd : ',time1-time0
   time0 = time1
   
   
   !  --- compute the quasiharmonic eigenvalues:
   
   do i=1,nmode
      if (eigval(i) > 0.0) then
         eigval(i) = 108.587*sqrt(0.6/eigval(i))
      else if (eigval(i) < 0.0) then
         eigval(i) = -108.587*sqrt(-0.6/eigval(i))
      else
         write(6,*) 'bad eigenvalue:',i,eigval(i)
      end if
   end do
   
   !  --- output the eigenvectors and eigenvalues:
   
   if (novec) goto 50
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
   print *, 'quasiharmonic eigenvector file'
   write(61,*) 'quasiharmonic eigenvector file'
   write(61,'(i5)') n
   write(61,'(7f11.5)') (xav(i),i=1,n)
#else
   write(61) n
   do i=1,n
      crdr(i) = xav(i)
   end do
   write(61) (crdr(i),i=1,n)
#endif
   j = nmode
   do k = 1,nmode
#ifdef FORMATTED
      write (61,'(1x,a)') '****'
      write (61,'(i5,f12.5)') k, eigval(j)
#else
      eig4 = eigval(j)
      write(61) k,eig4
#endif
      do i=1,n
#ifdef FORMATTED
         vout(i,j) = vout(i,j)*atmi(i)
#else
         crdr(i) = vout(i,j)*atmi(i)
#endif
      end do
#ifdef FORMATTED
      write (61,'(7f11.5)') (vout(i,j),i=1,n)
#else
      write(61) (crdr(i),i=1,n)
#endif
      j = j - 1
   end do
   close(unit=61)
   50 continue
   
   !  --- get thermochemistry
   
   print *,'now do the thermochemistry stuff'
   do i = 1,nmode
      w(i) = eigval(nmode-i+1)
   end do
   temp = 298.15
   pressure = 1.0
   
   call thermo(natom,nmode,1,0,xav,atmass,w,vout(1,1),vout(1,2), &
         vout(1,3),vout(1,4),temp,pressure)
   
   call wallclock(time1)
   write( 6,'(a,f8.2)') '| Time for final processing: ',time1-time0
   
end program quasih 
!-------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine usage here]
subroutine usage()
   
   !    maximum number of atoms, snapshots: modes:
   
   parameter (natmax=2400,natm3=3*natmax)
   parameter (nspmax=50)
   parameter (nmodemax=50,ncvmax=2*nmodemax,lworkl=ncvmax*(ncvmax+8))
   write(6,'(a,a)') &
         'usage: quasih -natom # [-f # -m <mass> -x <xyz> ', &
         ' [-novec | -v <vecs>] [ -first #] [ -last # ] [ -nmode # ]'
   write(6,'(a)')    'Compiled-in maximum sizes:'
   write(6,'(a,i6)') '                    natom: ',natmax
   write(6,'(a,i6)') '                shapshots: ',nspmax
   write(6,'(a,i6)') '                    modes: ',nmodemax
   call mexit(6,1)
end subroutine usage 
