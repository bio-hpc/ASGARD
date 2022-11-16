program fix_new_inpcrd_vel
   implicit none
   character(len=80) :: filenew,fileold,title
   character(len=80) :: fmt,fmtin,dtype
   integer ier,ionerr,iok,nf,inpcrd_unit,j,n,dim1,num_list
   double precision,allocatable :: crd(:,:),vel(:,:),acc(:,:),old_acc(:,:)
   double precision :: time,a,b,c,alpha,beta,gam
   double precision, parameter ::  &
                     sander_to_tinker_time_convert = 20.455d0
   double precision, parameter ::  &
                     tinker_to_sander_time_convert = 1.d0 / 20.455d0
  double precision :: zero=0.d0
  integer :: values(8)
  character(len=12) :: date,time1,zone
  character(len=8) :: word,word1
  character(len=16) :: word2

   write(6,*)'new_format inpcrd file? (must exist) :'
   read(5,*)fileold
   nf = 9
   open(unit=nf,file=fileold,status='old')

   fmtin = '(a)'
   dtype = 'TITLE'
   ionerr = 1 ! not fatal if missing ---should be there
   call nxtsec(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
   if ( iok == 0 )then
      read(nf,fmt)title
   else
      write(6,*)'Warning : title missing'
      title = ' '
   endif
   fmtin = '(E16.8)'
   dtype = 'ATOMIC_COORDS_SIMULATION_TIME'
   ionerr = 1 ! not fatal if missing ---should be there
   call nxtsec(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
   if ( iok == 0 )then
      read(nf,fmt)time
   else
      write(6,*)'Warning : time missing'
      time = 0.d0
   endif
   ! coordinates
   fmtin = '(I8)'
   dtype = 'ATOMIC_COORDS_NUM_LIST'
   ionerr = 1 ! not fatal if missing
   call nxtsec(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
   if ( iok == 0 )then !this data type found
      read(nf,fmt)num_list
   else
      write(6,*)'ATOMIC_COORDS_NUM_LIST flag missing'
      stop
   endif
   dim1 = 3
   allocate(crd(3,num_list),stat=ier)
   if ( ier /= 0 )then
      write(6,*)'problem allocating crds: dims 3 x ',num_list
      stop
   endif
   ionerr = 0 !fatal if missing
   fmtin = '(5E16.8)'
   dtype = 'ATOMIC_COORDS_LIST'
   call nxtsec(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
   read(nf,fmt)((crd(j,n),j=1,dim1),n=1,num_list)
   ! velocities
   fmtin = '(I8)'
   dtype = 'ATOMIC_VELOCITIES_NUM_LIST'
   ionerr = 1 ! not fatal if missing
   call nxtsec(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
   if ( iok == 0 )then !this data type found
      read(nf,fmt)num_list
   else
      write(6,*)'ATOMIC_VELOCITIES_NUM_LIST flag missing'
      stop
   endif
   dim1 = 3
   allocate(vel(3,num_list),stat=ier)
   if ( ier /= 0 )then
      write(6,*)'problem allocating vel: dims 3 x ',num_list
      stop
   endif
   ionerr = 0 !fatal if missing
   fmtin = '(5E16.8)'
   dtype = 'ATOMIC_VELOCITIES_LIST'
   call nxtsec(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
   read(nf,fmt)((vel(j,n),j=1,dim1),n=1,num_list)
   ! accelerations
   fmtin = '(I8)'
   dtype = 'ATOMIC_ACCELERATIONS_NUM_LIST'
   ionerr = 1 ! not fatal if missing
   call nxtsec(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
   if ( iok == 0 )then !this data type found
      read(nf,fmt)num_list
   else
      write(6,*)'ATOMIC_ACCELERATIONS_NUM_LIST flag missing'
      stop
   endif
   dim1 = 3
   allocate(acc(3,num_list),stat=ier)
   if ( ier /= 0 )then
      write(6,*)'problem allocating acc: dims 3 x ',num_list
      stop
   endif
   ionerr = 0 !fatal if missing
   fmtin = '(5E16.8)'
   dtype = 'ATOMIC_ACCELERATIONS_LIST'
   call nxtsec(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
   read(nf,fmt)((acc(j,n),j=1,dim1),n=1,num_list)
   ! old accelerations
   fmtin = '(I8)'
   dtype = 'OLD_ATOMIC_ACCELERATIONS_NUM_LIST'
   ionerr = 1 ! not fatal if missing
   call nxtsec(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
   if ( iok == 0 )then !this data type found
      read(nf,fmt)num_list
   else
      write(6,*)'OLD_ATOMIC_ACCELERATIONS_NUM_LIST flag missing'
      stop
   endif
   dim1 = 3
   allocate(old_acc(3,num_list),stat=ier)
   if ( ier /= 0 )then
      write(6,*)'problem allocating old_acc: dims 3 x ',num_list
      stop
   endif
   ionerr = 0 !fatal if missing
   fmtin = '(5E16.8)'
   dtype = 'OLD_ATOMIC_ACCELERATIONS_LIST'
   call nxtsec(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
   read(nf,fmt)((old_acc(j,n),j=1,dim1),n=1,num_list)
   ! unit cell info
   fmtin = '(3e16.8)'
   dtype = 'UNIT_CELL_PARAMETERS'
   ionerr = 1 ! not fatal if missing
   call nxtsec(nf,  6,  ionerr,fmtin,  dtype,  fmt,  iok)
   if ( iok == 0 )then !this data type found in prmtop
      read(nf,fmt)a,b,c,alpha,beta,gam
   else 
      a = 0.d0
   endif
   close(nf)
   write(6,*)'new vel units inpcrd file? (will create) :'
   read(5,*)filenew
   open(unit=inpcrd_unit,file=filenew,status='new')

   call date_and_time(date,time1,zone,values)
   write(inpcrd_unit,'(a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a)')  &
           '%VERSION  VERSION_STAMP = V0001.000  DATE = ', &
           values(2),'/',values(3),'/',values(1)-2000, &
           '  ',values(5),':',values(6),':',values(7), ' '
   write(inpcrd_unit,'(a)')'%FLAG TITLE'
  write(inpcrd_unit,'(a)')'%FORMAT(a)'
  write(inpcrd_unit,'(a)')title(1:len_trim(title))
  write(inpcrd_unit,'(a)')'%FLAG ATOMIC_COORDS_SIMULATION_TIME'
  write(inpcrd_unit,'(a)')'%FORMAT(E16.8)'
  write(inpcrd_unit,'(E16.8)')zero
!-----------------------------------------------------------
  write(inpcrd_unit,'(a)')'%FLAG ATOMIC_COORDS_NUM_LIST'
  write(inpcrd_unit,'(a)')'%FORMAT(i8)'
  write(inpcrd_unit,'(I8)')num_list
  write(word,'(I8)')num_list
  word1 = adjustl(word)
  word2 = '(3,'//word1(1:len_trim(word1))//')'
  write(inpcrd_unit,'(a)')'%FLAG ATOMIC_COORDS_LIST'
  write(inpcrd_unit,'(a)')'%COMMENT   dimension = '//word2
  write(inpcrd_unit,'(a)')'%FORMAT(3e20.12)'
!-----------------------------------------------------------
  !if ( dyn_read )then
     write(inpcrd_unit,'(3e20.12)')((crd(j,n),j=1,3),n=1,num_list)
  !else
     !write(inpcrd_unit,'(3e20.12)')((xyz_crd(j,n),j=1,3),n=1,num_list)
  !endif
  ! velocities, accelerations
  !if ( dyn_read )then
     write(inpcrd_unit,'(a)')'%FLAG ATOMIC_VELOCITIES_NUM_LIST'
     write(inpcrd_unit,'(a)')'%FORMAT(i8)'
     write(inpcrd_unit,'(I8)')num_list
     write(inpcrd_unit,'(a)')'%FLAG ATOMIC_VELOCITIES_LIST'
     write(inpcrd_unit,'(a)')'%COMMENT   dimension = '//word2
     write(inpcrd_unit,'(a)')'%FORMAT(3e20.12)'
     write(inpcrd_unit,'(3e20.12)') &
       ((tinker_to_sander_time_convert*vel(j,n),j=1,3),n=1,num_list)
     write(inpcrd_unit,'(a)')'%FLAG ATOMIC_ACCELERATIONS_NUM_LIST'
     write(inpcrd_unit,'(a)')'%FORMAT(i8)'
     write(inpcrd_unit,'(I8)')num_list
     write(inpcrd_unit,'(a)')'%FLAG ATOMIC_ACCELERATIONS_LIST'
     write(inpcrd_unit,'(a)')'%COMMENT   dimension = '//word2
     write(inpcrd_unit,'(a)')'%FORMAT(3e20.12)'
     write(inpcrd_unit,'(3e20.12)') &
       ((tinker_to_sander_time_convert**2*acc(j,n),j=1,3),n=1,num_list)
     write(inpcrd_unit,'(a)')'%FLAG OLD_ATOMIC_ACCELERATIONS_NUM_LIST'
     write(inpcrd_unit,'(a)')'%FORMAT(i8)'
     write(inpcrd_unit,'(I8)')num_list
     write(inpcrd_unit,'(a)')'%FLAG OLD_ATOMIC_ACCELERATIONS_LIST'
     write(inpcrd_unit,'(a)')'%COMMENT   dimension = '//word2
     write(inpcrd_unit,'(a)')'%FORMAT(3e20.12)'
     write(inpcrd_unit,'(3e20.12)') &
       ((tinker_to_sander_time_convert**2*old_acc(j,n),j=1,3), &
                                                             n=1,num_list)
  !endif
!-----------------------------------------------------------
  write(inpcrd_unit,'(a)')'%FLAG UNIT_CELL_PARAMETERS'
  write(inpcrd_unit,'(a,a)')'%COMMENT lengths a,b,c; ', &
                          'then angles alpha,beta,gamma'
  write(inpcrd_unit,'(a)')'%FORMAT(3e20.12)'
  !if ( dyn_read )then
    ! write(inpcrd_unit,'(3e20.12)') &
         !(dyn_cell_length(j),j=1,3),(dyn_cell_angles(j),j=1,3)
  !else
     write(inpcrd_unit,'(3e20.12)')a,b,c,alpha,beta,gam
  !endif
end program fix_new_inpcrd_vel
