program new_to_old_crd
   implicit none
   character(len=80) :: filenew,fileold,title
   character(len=80) :: fmt,fmtin,dtype
   integer ier,ionerr,iok,nf,j,n,dim1,num_list
   double precision,allocatable :: crd(:,:),vel(:,:),acc(:,:),old_acc(:,:)
   double precision :: time,a,b,c,alpha,beta,gam
   double precision, parameter ::  &
                     sander_to_tinker_time_convert = 20.455d0
   double precision, parameter ::  &
                     tinker_to_sander_time_convert = 1.d0 / 20.455d0

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
   write(6,*)'Tinker dyn file? (will create) :'
   read(5,*)filenew
   open(unit=nf,file=filenew,status='new')
   write(nf,*)' Number of Atoms and Title :'
   write(nf,'(i6,2x,a)')num_list,title(1:len_trim(title))
   write(nf,*)' Periodic Box Dimensions :'
   write(nf,'(3d26.16)')a,b,c
   write(nf,'(3d26.16)')alpha,beta,gam
   write(nf,*)' Current Atomic Positions :'
   write(nf,'(3d26.16)')((crd(j,n),j=1,3),n=1,num_list)
   write(nf,*)' Current Atomic Velocities :'
   write(nf,'(3d26.16)') &
       ((sander_to_tinker_time_convert*vel(j,n),j=1,3),n=1,num_list)
   write(nf,*)' Current Atomic Accelerations :'
   write(nf,'(3d26.16)') &
       ((sander_to_tinker_time_convert**2*acc(j,n),j=1,3),n=1,num_list)
   write(nf,*)' Previous Atomic Accelerations :'
   write(nf,'(3d26.16)') &
       ((sander_to_tinker_time_convert**2*old_acc(j,n),j=1,3),n=1,num_list)
end program new_to_old_crd
