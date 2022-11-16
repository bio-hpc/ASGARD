! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"

subroutine mexit(filenum, exitstatus)
   implicit none
   integer filenum
   integer exitstatus

   if (filenum > 0) then  ! close this unit if greater than zero
      close(unit=filenum)
   endif
   ! exit status; error if non-zero
#  if XLF90 || IBM3090 || F2C
   if (exitstatus.ne.0) stop 1; stop 0
#  else
   call exit(exitstatus)
#  endif
end subroutine mexit

subroutine private_getx(natom,pathlen,fpath,x)
   implicit none
   integer natom, pathlen, filenum
   character(len=pathlen) fpath
   _REAL_ x(*)

   integer nr3, ier, i
   character(len=80) title1
   character(len=256) testbuffer
   filenum=9
   nr3=3*natom
   call myopen(filenum,fpath(1:pathlen),'O','F','R')
   read(filenum,'(a80)') title1
   read(filenum,'(a80)') testbuffer
   if( testbuffer(6:6) == ' ' ) then ! this is an old file
      read(testbuffer,'(i5,e15.7)') natom
   else                              ! assume a new file
      read(testbuffer,'(i6,e15.7)') natom
   end if
   read(filenum,'(6f12.7)') (x(i),i=1,nr3)
   close(filenum, iostat=ier)
   return
end subroutine private_getx

subroutine prepb_read(passed_ipb,passed_inp,mycn1,mycn2,mynttyp)
!  use pbtimer_module
   use parms, only : cn1, cn2
   implicit none
#  include "../include/md.h"
#  include "pb_md.h"
#  include "flocntrl.h"
   _REAL_ mycn1(*),mycn2(*)
   integer passed_ipb,passed_inp,mynttyp,ierror
!  character(len=8) initial_date
!  character(len=10) initial_time

   if ( passed_ipb /= 1 .and. passed_ipb /=2 ) then
      print *,"ipb should be either 1 or 2."; call mexit(6,1)
   endif
   ! We need to allocate our CN1 and CN2 arrays here
   allocate(cn1(mynttyp), cn2(mynttyp), stat=ierror)
   cn1(1:mynttyp)=mycn1(1:mynttyp)
   cn2(1:mynttyp)=mycn2(1:mynttyp)
   imin = 0
   igb = 0
   ipb = passed_ipb
   inp = passed_inp
   npbstep = 1
!  call date_and_time( initial_date, initial_time )
!  call pbtimer_init()
!  call pbtimer_start(PBTIME_TOTAL)
   do_pbfd = 1
   do_pbnp = 1
   do_pbdir = 1
   idecomp = 0


   return
end subroutine prepb_read

! The interface of the FORTRAN side
! Author: Mengjuei Hsieh
subroutine mypb_force(natom,nres,ntypes,ipres,iac,ico,exclat,mystep,&
                   cn1,cn2,cg,xx,f,epol,evdw,eelt,esurf,edisp)
   use poisson_boltzmann, only : pb_force, &
                                 xmin, ymin, zmin, xmax, ymax, zmax
   use dispersion_cavity, only : np_force
   use decomp, only : allocate_int_decomp, allocate_real_decomp, &
                      deallocate_int_decomp, deallocate_real_decomp
   implicit none
#  include "../include/md.h"
#  include "pb_md.h"

   integer natom,   nres,  ntypes,npdec,    mystep
   integer ipres(*),iac(*),ico(*),exclat(*)
   _REAL_  cn1(*),  cn2(*),cg(*), xx(*),    f(*)
   _REAL_  evdw,eelt,epol,esurf,edisp
   ! Initialize the cpu timer. Needed for machines where returned cpu times
   f(1:natom*3)=0d0
   if ( mystep /= 0 ) npbstep = mystep
   if ( npbstep < 1 ) then
      ntnbr = 1; ntnba = 1; pbgrid = .true.
   else
      if ( mod(npbstep,nsnbr) == 0 ) ntnbr = 1
      if ( mod(npbstep,nsnba) == 0 ) ntnba = 1
      if ( mod(npbstep,npbgrid) == 0 ) then
         pbgrid = .true.
         xmin = 0d0; ymin = 0d0; zmin = 0d0
         xmax = 0d0; ymax = 0d0; zmax = 0d0
      endif
   endif

   ! To make it compatible with decomposition facility
   if ( idecomp > 0 ) then
      call allocate_int_decomp(natom, nres)
   else
      call allocate_int_decomp(1, 1)
   endif

   if( ipb /= 0 ) then
      call pb_force(natom,nres,ntypes,npdec,ipres,iac,ico,exclat, &
                    cn1,cn2,cg,xx,f,evdw,eelt,epol)
   end if
   if ( pbgrid ) pbgrid = .false.
   if ( pbinit ) pbinit = .false.
   if ( inp /= 0 ) then
      esurf = 0.0d0; edisp = 0.0d0
      call np_force(natom,nres,ntypes,ipres,iac,ico, &
                    cn1,cn2,xx,f,esurf,edisp)
   end if
   ! Yes, I know this sucks, bear with me.
   call deallocate_int_decomp()
   call flush(6)
   return
end subroutine mypb_force

!subroutine clockstop
!   use pbtimer_module
!   implicit none
!   character(len=8) final_date
!   character(len=10) final_time
!   call pbtimer_stop(PBTIME_TOTAL)
!   call date_and_time( final_date, final_time )
!   return
!end subroutine clockstop
