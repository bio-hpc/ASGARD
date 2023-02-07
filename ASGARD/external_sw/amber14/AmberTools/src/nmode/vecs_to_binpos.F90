
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of program vecs_to_binpos here]
program vecs_to_binpos

   !    make a trajectory-player binpos file from a vecs file.
   
   !    At most NVEC separate binpos files are all concatenated on stdout.
   !    Use "split -b#" to split the resulting file into individual
   !    binpos files for each normal mode.
   
   parameter( maxatom3=90000)
   parameter(nvec=100)
   dimension x0(maxatom3), vec(maxatom3), x(maxatom3)
   character(len=80) title,arg
   
   !     --- parse the command line:
   
   scale = 0.1
   indx = iargc()
   if (indx == 0) goto 20
   10 continue
   iarg = iarg + 1
   call getarg(iarg,arg)
   if (arg == '-scale') then
      iarg = iarg + 1
      call getarg(iarg,arg)
      read( arg, '(f5.2)' ) scale
      write( 0, '(a,f8.3)') 'setting scale to ',scale
   else
      if (arg == ' ') goto 20
      write(0,'(/,5x,a,a)') 'unknown flag: ',arg
      call mexit(6, 1)
   end if
   if (iarg < indx) goto 10
   20 continue
   
   read(5,'(a)') title
   write(0,'(a)') title
   read(5,'(i5)') nat3
   if( nat3 > maxatom3 ) then
      write(0,*) 'too many atoms: ',nat3,maxatom3
      call mexit(6,1)
   end if
   write(0,'(i5)') nat3
   natoms = nat3/3
   read(5,'(7f11.5)') (x0(i),i=1,nat3)
   
   !    loop over the first NVEC vectors:
   
   do ivec=1,nvec
      
      read(5,'(a)',end=99) title
      write(0,'(a)') title
      read(5,'(i5,f12.5)') mode, freq
      write(0,'(i5,f12.5)') mode, freq
      read(5,'(7f11.5)') (vec(i),i=1,nat3)
      
      call startbinpos()
      del = -150.0*scale
      do j=1,21
         do i=1,nat3
            x(i) = x0(i) + del*vec(i)
         end do
         call writebinpos(natoms, x)
         del = del + 15.0*scale
      end do
      do j=1,21
         do i=1,nat3
            x(i) = x0(i) + del*vec(i)
         end do
         call writebinpos(natoms, x)
         del = del - 15.0*scale
      end do
      
   end do
   99 continue
   write(0,'(a,i5,a)') 'processed ',ivec-1,' eigenvectors'
   nbytes = 42*( 4 + 4*nat3 ) + 4
   write(0,'(a,i6,a)') 'each binpos file will be ',nbytes, &
         ' bytes in length'
   
end program vecs_to_binpos 
