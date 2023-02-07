#include "../include/dprec.fh"
module sebomd_module

  implicit none
  private

  ! subroutines
  public :: sebomd_setup
  public :: sebomd_write_info
  public :: sebomd_check_options
  public :: read_sebomd_namelist
  public :: sebomd_namelist_default
  public :: sebomd_gradient_write
  public :: sebomd_hessian_write
  public :: sebomd_hessian_compute
  public :: sebomd_open_files
  public :: sebomd_close_files
  public :: sebomd_save_forces
#ifdef MPI
  public :: sebomd_bcast_obj
#endif

  type sebomd_structure
    logical :: do_sebomd ! flag for SEBOMD (YES/NO)
    integer :: idcflag
    integer :: iflagch
    integer :: iflagch_old
    integer :: pdmx
    integer :: nhessian
    integer :: diverror
    integer :: ctype
    _REAL_ :: esebomd
    _REAL_, dimension(3) :: dbox
    ! namelist replica {
    character(10) :: hamiltonian
    character(10) :: modif  ! specific modification of the hamiltonian
    integer :: ncore
    _REAL_ :: dbuff1
    _REAL_ :: dbuff2
    character(256) :: charge_out
    _REAL_ :: lambda
    _REAL_ :: peptk   ! peptidic constant
    _REAL_ :: dpmax   ! density matrix convergence criteria
    integer :: method
    integer :: charge
    integer :: longrange
    integer :: fullscf
    integer :: ntwc
    integer :: chtype
    integer :: chewald
    integer :: screen
    integer :: guess
    integer :: pdump
    integer :: ipolyn
    integer :: nresidue
    integer :: ntwh
    integer :: iprec
    integer :: peptcorr  ! boolean flag for peptidic correction
    integer :: debugmsg
    integer :: debugforces
    integer :: diag_routine
    ! } namelist replica
  end type sebomd_structure

  type(sebomd_structure) :: sebomd_obj

  ! objects
  public :: sebomd_obj
contains
!------------------------------------------------------------------------------
  subroutine sebomd_setup()
    implicit none

    sebomd_obj%do_sebomd = .false.
    sebomd_obj%idcflag = 0
    return
  end subroutine sebomd_setup
!------------------------------------------------------------------------------
  subroutine sebomd_namelist_default()
    ! namelist default {
    sebomd_obj%hamiltonian = "PM3"
    sebomd_obj%modif = "none"
    sebomd_obj%method = 0
    sebomd_obj%ncore = 1
    sebomd_obj%dbuff1 = 6.0
    sebomd_obj%dbuff2 = 0.0
    sebomd_obj%charge = 0
    sebomd_obj%longrange = 0
    sebomd_obj%iprec = 4
    sebomd_obj%fullscf = 0
    sebomd_obj%ntwc = 0
    sebomd_obj%chtype = 0
    sebomd_obj%chewald = 0
    sebomd_obj%screen = 0
    sebomd_obj%guess = 0
    sebomd_obj%pdump = 0
    sebomd_obj%ipolyn = 1
    sebomd_obj%nresidue = 0
    sebomd_obj%ntwh = 0
    sebomd_obj%lambda = 1.0d0
    sebomd_obj%peptcorr = 0
    sebomd_obj%peptk = 0.0d0
    sebomd_obj%charge_out = "sebomd.chg"
    sebomd_obj%debugmsg = 0
    sebomd_obj%debugforces = 0
    sebomd_obj%diag_routine = 2
    sebomd_obj%dpmax = 1e-7
    ! } namelist default
    return
  end subroutine  sebomd_namelist_default
!------------------------------------------------------------------------------
  subroutine read_sebomd_namelist()
    implicit none
    integer :: stat

    namelist /sebomd/ hamiltonian, &
                      modif, &
                      ncore, &
                      dbuff1, &
                      dbuff2, &
                      charge_out, &
                      lambda, &
                      peptk, &
                      method, &
                      charge, &
                      longrange, &
                      fullscf, &
                      ntwc, &
                      chtype, &
                      chewald, &
                      screen, &
                      guess, &
                      pdump, &
                      ipolyn, &
                      nresidue, &
                      ntwh, &
                      iprec, &
                      peptcorr, &
                      debugmsg, &
                      debugforces, &
                      diag_routine, &
                      dpmax
   
    character(10) :: hamiltonian
    character(10) :: modif
    integer :: method
    integer :: ncore
    _REAL_ :: dbuff1
    _REAL_ :: dbuff2
    character(256) :: charge_out
    _REAL_ :: lambda
    _REAL_ :: peptk
    _REAL_ :: dpmax
    integer :: charge
    integer :: longrange
    integer :: fullscf
    integer :: ntwc
    integer :: chtype
    integer :: chewald
    integer :: screen
    integer :: guess
    integer :: pdump
    integer :: ipolyn
    integer :: nresidue
    integer :: ntwh
    integer :: peptcorr
    integer :: iprec
    integer :: debugmsg
    integer :: debugforces
    integer :: diag_routine

    hamiltonian =   sebomd_obj%hamiltonian
    modif =   sebomd_obj%modif
    method =   sebomd_obj%method
    ncore =   sebomd_obj%ncore
    dbuff1 =   sebomd_obj%dbuff1
    dbuff2 =   sebomd_obj%dbuff2
    charge =   sebomd_obj%charge
    longrange =   sebomd_obj%longrange
    iprec =   sebomd_obj%iprec
    fullscf =   sebomd_obj%fullscf
    ntwc =   sebomd_obj%ntwc
    chtype =   sebomd_obj%chtype
    chewald =   sebomd_obj%chewald
    screen =   sebomd_obj%screen
    guess =   sebomd_obj%guess
    pdump =   sebomd_obj%pdump
    ipolyn =   sebomd_obj%ipolyn
    nresidue =   sebomd_obj%nresidue
    ntwh =   sebomd_obj%ntwh
    lambda =   sebomd_obj%lambda
    peptcorr =   sebomd_obj%peptcorr
    peptk =   sebomd_obj%peptk
    charge_out = sebomd_obj%charge_out
    debugmsg = sebomd_obj%debugmsg
    debugforces = sebomd_obj%debugforces
    diag_routine = sebomd_obj%diag_routine
    dpmax = sebomd_obj%dpmax


    read(unit=5, nml=sebomd, iostat = stat)
    if (stat /= 0) then
      write(6,'(A)') 'Error reading &sebomd namelist'
      call mexit(6,1)
    end if
 
    sebomd_obj%hamiltonian =   hamiltonian
    sebomd_obj%modif =   modif
    sebomd_obj%method =   method
    sebomd_obj%ncore =   ncore
    sebomd_obj%dbuff1 =   dbuff1
    sebomd_obj%dbuff2 =   dbuff2
    sebomd_obj%charge =   charge
    sebomd_obj%longrange =   longrange
    sebomd_obj%iprec =   iprec
    sebomd_obj%fullscf =   fullscf
    sebomd_obj%ntwc =   ntwc
    sebomd_obj%chtype =   chtype
    sebomd_obj%chewald =   chewald
    sebomd_obj%screen =   screen
    sebomd_obj%guess =   guess
    sebomd_obj%pdump =   pdump
    sebomd_obj%ipolyn =   ipolyn
    sebomd_obj%nresidue =   nresidue
    sebomd_obj%ntwh =   ntwh
    sebomd_obj%lambda =   lambda
    sebomd_obj%peptk =   peptk
    sebomd_obj%peptcorr =   peptcorr
    sebomd_obj%charge_out = charge_out
    sebomd_obj%debugmsg = debugmsg
    sebomd_obj%debugforces = debugforces
    sebomd_obj%diag_routine = diag_routine
    sebomd_obj%dpmax = dpmax
    return
  end subroutine read_sebomd_namelist
!------------------------------------------------------------------------------
  subroutine sebomd_write_info()
    implicit none
   ! assert modif keyword
   if (sebomd_obj%modif(1:len_trim(sebomd_obj%modif)).eq."none") then
     sebomd_obj%ctype = 0
   elseif (sebomd_obj%modif(1:len_trim(sebomd_obj%modif)).eq."PIF2") then
     sebomd_obj%ctype = 1
   elseif (sebomd_obj%modif(1:len_trim(sebomd_obj%modif)).eq."PIF3") then
     sebomd_obj%ctype = 2
   elseif (sebomd_obj%modif(1:len_trim(sebomd_obj%modif)).eq."MAIS1") then
     sebomd_obj%ctype = 3
   elseif (sebomd_obj%modif(1:len_trim(sebomd_obj%modif)).eq."MAIS2") then
     sebomd_obj%ctype = 4
   else
     write(6,*)
     write(6,'(a,a,a)') "ERROR: modif keyword '" &
                       ,sebomd_obj%modif(1:len_trim(sebomd_obj%modif)) &
                       ,"' IS NOT RECOGNIZED"
     call mexit(6,1)
   endif

   ! TODO assert peptcorr and peptk
   if (sebomd_obj%peptcorr.ne.0) then
     ! peptide correction
     ! keep value from user if (s)he added one
     if (sebomd_obj%peptk.eq.0.0d0) then
       ! no input value => let's give one
       if (sebomd_obj%hamiltonian(1:len_trim(sebomd_obj%hamiltonian)).eq."AM1") then
!        original MOPAC Value:
!         sebomd_obj%peptk = 3.3191D0
!        O. Ludwig et al. J. Mol. Model. 1996,2,341-350
         sebomd_obj%peptk = 5.9864D0

       else if (sebomd_obj%hamiltonian(1:len_trim(sebomd_obj%hamiltonian)).eq."PM3") then
!        original MOPAC Value:
!         sebomd_obj%peptk = 7.1853D0
!        O. Ludwig et al. J. Mol. Model. 1996,2,341-350
         sebomd_obj%peptk = 9.8526D0

       else if (sebomd_obj%hamiltonian(1:len_trim(sebomd_obj%hamiltonian)).eq."MNDO") then
         sebomd_obj%peptk = 6.1737D0
       else
         write(6,'("ERROR: no default peptidic correction value for ",a10," hamiltonian")')  &
           sebomd_obj%hamiltonian(1:len_trim(sebomd_obj%hamiltonian))
         write(6,'("Please provide one with the peptk keyword")')
         call mexit(6,1)
       endif
     endif
   endif

!   write(6,'("SEBOMD options:")')
!   write(6,'("Hamiltonian = ",a)') sebomd_obj%hamiltonian
    write(6,9800)
    write(6,'(5x,2(a,a10),2(a,i10))') 'hamiltonian =',adjustr(sebomd_obj%hamiltonian), &
       ', modif = ',adjustr(sebomd_obj%modif), &
       ',  longrange   =',sebomd_obj%longrange, ',  method      =',sebomd_obj%method
    if (sebomd_obj%method > 0) then
       write(6,'(5x,a,i10,a,f10.4,a,f10.4)') 'ncore       =',sebomd_obj%ncore, &
       ',  dbuff1      =', sebomd_obj%dbuff1, ',  dbuff2      =', sebomd_obj%dbuff2
    endif
    write(6,'(5x,a,i10,a,e10.3,a,i10)') 'charge      =',sebomd_obj%charge,',  dpmax       =', &
       sebomd_obj%dpmax, ',  fullscf     =',sebomd_obj%fullscf
    write(6,'(5x,3(a,i10))') 'ipolyn      =',sebomd_obj%ipolyn,',  pdump       =',sebomd_obj%pdump, &
          ',  guess       =',sebomd_obj%guess
    if (sebomd_obj%longrange == 1) then
       write(6,'(5x,4(a,i10))') 'ntwc        =',sebomd_obj%ntwc,',  chtype      =',sebomd_obj%chtype, &
        ',  screen      =', sebomd_obj%screen, &
        ',  ntwh        =', sebomd_obj%ntwh
    elseif (sebomd_obj%longrange == 4) then
       write(6,'(5x,4(a,i10))') 'ntwc        =',sebomd_obj%ntwc,',  chewald     =',sebomd_obj%chewald, &
        ',  screen      =', sebomd_obj%screen, &
        ',  ntwh        =', sebomd_obj%ntwh
    else
       write(6,'(5x,3(a,i10))') 'ntwc        =',sebomd_obj%ntwc, ',  screen      =', sebomd_obj%screen, &
        ',  ntwh        =', sebomd_obj%ntwh
    endif
    write(6,'(5x,a,i10,a,f10.4)') 'peptcorr    =',sebomd_obj%peptcorr, ', peptk =', sebomd_obj%peptk
    write(6,*)
9800 format(/80('-')/,'   SEBOMD  DATA  FOR  THE  RUN',/80('-')/)
    return
  end subroutine sebomd_write_info
!------------------------------------------------------------------------------
#ifdef MPI
  subroutine sebomd_bcast_obj()
    implicit none
#include "parallel.h"
#include "mpif.h"
    integer :: ierr
    call mpi_bcast(sebomd_obj%do_sebomd  ,  1, MPI_LOGICAL, 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%idcflag    ,  1, MPI_INTEGER, 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%iflagch    ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%iflagch_old,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%pdmx       ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%nhessian   ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%diverror   ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%esebomd    ,  1, MPI_DOUBLE_PRECISION , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%dbox       ,  3, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%hamiltonian, 10, MPI_CHARACTER, 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%modif      , 10, MPI_CHARACTER, 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%ncore      , 10, MPI_INTEGER, 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%dbuff1     , 1, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%dbuff2     , 1, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%charge_out ,256, MPI_CHARACTER, 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%lambda     ,  1, MPI_DOUBLE_PRECISION , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%peptk      ,  1, MPI_DOUBLE_PRECISION , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%dpmax      ,  1, MPI_DOUBLE_PRECISION , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%method     ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%charge     ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%longrange  ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%fullscf    ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%ntwc       ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%chtype     ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%chewald    ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%screen     ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%guess      ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%pdump      ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%ipolyn     ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%nresidue   ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%ntwh       ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%peptcorr   ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%ctype      ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%iprec      , 10, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%debugmsg   ,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%debugforces,  1, MPI_INTEGER , 0, commsander, ierr)
    call mpi_bcast(sebomd_obj%diag_routine, 1, MPI_INTEGER , 0, commsander, ierr)
    return
  end subroutine sebomd_bcast_obj
#endif /* MPI */
!------------------------------------------------------------------------------
  subroutine sebomd_check_options()

   implicit none
   
#include "box.h"
#include "parallel.h"

#ifdef MPI
   if (sebomd_obj%method == 0) then 
      write(6,*) 'SEBOMD Parallel version can only be used with method > 0'
      call mexit(6,1)
   endif
#endif
   if (sebomd_obj%method > 3) then
      write(6,*) 'method keyword in &sebomd can only be equal to 0, 1, 2 or 3'
      call mexit(6,1) 
   endif

   if (sebomd_obj%longrange > 4) then
      write(6,*) 'longrange keyword in &sebomd can''t be greater than 3'
      call mexit(6,1)
   endif

   if (sebomd_obj%fullscf > 1) then
      write(6,*) 'fullscf keyword in &sebomd can''t be greater than 1'
      call mexit(6,1)
   endif

   if (sebomd_obj%screen > 2) then
      write(6,*) 'screen keyword in &sebomd can''t be greater than 2'
      call mexit(6,1)
   endif

   call sebomd_check_option_ntp()

   end subroutine sebomd_check_options
!------------------------------------------------------------------------------
   subroutine sebomd_check_option_ntp()
     implicit none
#include "../include/md.h"
     if (ntp > 0 .and. barostat /= 2) then
       write(6,*) 'SEBOMD Error: ntp = ',ntp,' and barostat = ',barostat
       write(6,*) 'only MC Barostat is allowed for NPT simulations'
       call mexit(6,1)
     endif
     return
   end subroutine sebomd_check_option_ntp
!------------------------------------------------------------------------------
  subroutine sebomd_hessian_compute(xx,ix,ih,ipairs,x,ener, &
                             qsetup, do_list_update, nstep)
     use state
     implicit none
     _REAL_ xx(*)
     _REAL_ x(*)
     integer   ix(*)
     integer   ipairs(*)
     character(len=4) ih(*)
     integer nstep
     logical qsetup
     logical do_list_update
     type(state_rec) :: ener
#include "../include/memory.h"
     _REAL_ deltax, oldx, hij
     integer nij
     integer i, j
  
  ! compute hessian for the first 20 atoms if necessary (see locmem for memory allocation)
  ! H(i,j) = d g(i) / dj = (g(xi, xj+dx)-g(xi, xj))/dx    for all i
     deltax = 1e-3
     xx(hessian:hessian+sebomd_obj%nhessian*(sebomd_obj%nhessian+1)/2) = 0.0d0
     do j = 1, sebomd_obj%nhessian
        ! move a point
        oldx = x(j)
        x(j) = x(j) - deltax
        ! compute new forces
        xx(grad3tmp:grad3tmp-1+sebomd_obj%nhessian) = 0.0d0
        call force(xx,ix,ih,ipairs,x,xx(grad3tmp),ener,ener%vir, &
           xx(l96), xx(l97), xx(l98), xx(l99), qsetup, &
           do_list_update, nstep)
        x(j) = x(j) + deltax + deltax
        ! compute new forces
        xx(grad4tmp:grad4tmp-1+sebomd_obj%nhessian) = 0.0d0
        call force(xx,ix,ih,ipairs,x,xx(grad4tmp),ener,ener%vir, &
           xx(l96), xx(l97), xx(l98), xx(l99), qsetup, &
           do_list_update, nstep)
        ! compute hessian vector
        ! H(i.j) = (g(i,j+deltax)-g(i,j))/deltax
        do i = 1, sebomd_obj%nhessian
           hij = (xx(grad3tmp-1+i)-xx(grad4tmp-1+i))/(deltax+deltax)
           nij = max(i,j)*(max(i,j)-1)/2+min(i,j)
           if (i.ne.j) then
           xx(hessian-1+nij) = xx(hessian-1+nij) + hij*0.5    ! double count
           else
           xx(hessian-1+nij) = xx(hessian-1+nij) + hij*1.0    ! no count for diagonal
           endif
  !        write(6,'("hess",3i5,f20.10)') i,j,nij,hij
        end do
        ! back to the original point
        x(j) = oldx
     end do
     call sebomd_hessian_write(xx(hessian), sebomd_obj%nhessian*(sebomd_obj%nhessian+1)/2)
  !  do nij = 1, sebomd_obj%nhessian*(sebomd_obj%nhessian+1)/2
  !     write(6,'("hessian(", i5,") = ", f20.10)') nij, xx(hessian-1+nij)
  !  end do
     return
  end subroutine sebomd_hessian_compute
!------------------------------------------------------------------------------
  subroutine sebomd_hessian_write(hessian, n)
     implicit none
     _REAL_ hessian(*)
     integer n
  
     logical first
     save first
     data first /.true./
  
     integer i
  
     if (first) then
        call amopen(25,'divcon.hess','R','F','W')
        write(25,'(i5)') n
        first = .false.
     else
        call amopen(25,'divcon.hess','O','F','A')
     endif
     write(25,'(6f14.7)') (hessian(i),i=1,n)
     close(25)
     return
  end subroutine sebomd_hessian_write
!------------------------------------------------------------------------------
  subroutine sebomd_gradient_write(forces, n)
     implicit none
     _REAL_ forces(*)
     integer n
  
     logical first
     save first
     data first /.true./
  
     integer i
  
     if (first) then
        call amopen(25,'divcon.grad','R','F','W')
        write(25,'(i5)') n
        first = .false.
     else
        call amopen(25,'divcon.grad','O','F','A')
     endif
     write(25,'(6f14.7)') (-forces(i),i=1,n)
     close(25)
     return
  end subroutine sebomd_gradient_write
!------------------------------------------------------------------------------
  subroutine sebomd_open_files
    use file_io_dat, only : sechgunit
    if (sebomd_obj%ntwc.ne.0) then
      call amopen(sechgunit, sebomd_obj%charge_out,'U','F','W')
    end if
  end subroutine sebomd_open_files
!------------------------------------------------------------------------------
  subroutine sebomd_close_files
    use file_io_dat, only : sechgunit
    if (sebomd_obj%ntwc.ne.0) then
       close(sechgunit)
    end if
  end subroutine
!------------------------------------------------------------------------------
  subroutine sebomd_save_forces(iflag, n, f_mm, f_qm, ftmp1, ftmp2)
    ! routine to handle the saving of forces that must be
    ! apply to the system whatever the lambda term used
    ! What we want is:
    !   f(i) = (1-lambda)*f_MM(i) + lambda*f_QM(i) + f_restraint(i)
    ! Thus
    !   f(i) = (1-lambda)*(f_MM(i) + f_restraint(i))
    !            + lambda*(f_QM(i) + f_restraint(i))
    ! By default, forces.F90 computes f_MM(i) + f_restraint(i) in the same array
    !
    ! Purpose of this subroutine: extract f_restraint(i) then add it to f_QM(i)
    !
    ! parameters:
    !   n: number of atoms
    !   f_mm: the current array of forces (f in forces.F90, it sums up all forces)
    !   f_qm: the sebomd array of forces
    !   ftmp1: temporary array
    !   ftmp2: temporary array
    !   iflag: an integer which specifies which part of the job to do
    !     = 0: initialize
    !     = 1: before any restraint calculation -> store f_mm
    !     = 2: after any restraint calculation -> store f_restraint
    !     = 3: end: update f_qm with f_restraint
    implicit none
    integer :: iflag, n
    double precision, dimension (3*n) :: f_mm, f_qm, ftmp1, ftmp2

    integer :: i

    if (iflag.eq.0) then
      ! step 0: initialization (ftmp2 set to zero)
      do i = 1, 3*n
        ftmp2(i) = 0.0d0
      end do
    else if (iflag.eq.1) then
      ! step 1: we save forces before restraint computation
      do i = 1, 3*n
        ftmp1(i) = f_mm(i)
      end do
    else if (iflag.eq.2) then
      ! step 2: after restraint: we add the restraint forces to ftmp2
      do i = 1, 3*n
        ftmp2(i) = ftmp2(i) + f_mm(i)-ftmp1(i)
      end do
    else if (iflag.eq.3) then
      ! step 3: final step: f_qm is updated with all the restraint forces
      do i = 1, 3*n
        f_qm(i) = f_qm(i) + ftmp2(i)
      end do
    else
      ! should not be there
      continue
    end if
  end subroutine sebomd_save_forces
!------------------------------------------------------------------------------
end module sebomd_module
