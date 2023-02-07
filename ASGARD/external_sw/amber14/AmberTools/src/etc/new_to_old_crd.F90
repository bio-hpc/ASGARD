program new_to_old_crd
   implicit none
   character(len=80) :: filenew,fileold,title
   character(len=80) :: fmt,fmtin,dtype
   integer ier,ionerr,iok,nf,j,n,dim1,num_list
   double precision,allocatable :: crd(:,:)
   double precision :: time,a,b,c,alpha,beta,gam

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
   write(6,*)'old_format inpcrd file? (will create) :'
   read(5,*)filenew
   open(unit=nf,file=filenew,status='new')
   write(nf,'(a)')title
   if ( num_list > 100000 )then
      write(nf,'(i6,e15.7)')num_list,time
   else
      write(nf,'(i5,e15.7)')num_list,time
   endif
   write(nf,'(6f12.7)')((crd(j,n),j=1,3),n=1,num_list)
   if ( a > 0.d0 )write(nf,'(6f12.7)')a,b,c,alpha,beta,gam
end program new_to_old_crd
