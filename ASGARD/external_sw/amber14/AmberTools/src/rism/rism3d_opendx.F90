!<compile=optimized>

!The 3D-RISM-KH software found here is copyright (c) 2010-2012 by 
!Andriy Kovalenko, Tyler Luchko and David A. Case.
!
!This program is free software: you can redistribute it and/or modify it
!under the terms of the GNU General Public License as published by the Free
!Software Foundation, either version 3 of the License, or (at your option)
!any later version.
!
!This program is distributed in the hope that it will be useful, but
!WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
!or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!for more details.
!
!You should have received a copy of the GNU General Public License in the
!../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
!
!Users of the 3D-RISM capability found here are requested to acknowledge
!use of the software in reports and publications.  Such acknowledgement
!should include the following citations:
!
!1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999); 
!ibid. 112:10391-10417 (2000).   
!
!2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by  
!F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.  
!
!3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
!and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010). 

#include "../include/dprec.fh"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Minimimal support for text format OpenDX rectangular grid output.  There is 
!parallel support for writting but this is done through communicating with the
!master node so it is very slow.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module rism3d_opendx
  use safemem
  use rism_report_c
  implicit none
  character(*),parameter,private :: obj1Label = "object 1 class gridpositions counts"
  character(*),parameter,private :: obj2Label = "object 2 class gridconnections counts"
  character(*),parameter,private :: obj3Label = "object 3 class array type double rank 0 items"
  character(*),parameter,private :: originLabel = "origin"
  character(*),parameter,private :: deltaLabel = "delta"

  interface readDXHeader
     module procedure readDXHeader_file, readDXHeader_unit
  end interface readDXHeader

  interface readDX_p
     module procedure readDX_3D_p, readDX_1D_p
  end interface readDX_p

  interface readDX
     module procedure readDX_3D, readDX_1D
  end interface readDX

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Reads in a OpenDX file header.  Expects
!!!a 3D, regularly spaced grid in ASCII format.
!!!IN:
!!!   file   :: file name of unit. Only for error reporting
!!!   origin :: real space coordinate of origin
!!!   delta  :: gridspacing
!!!   npos   :: number of grid points in each dimension
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readDXHeader_file(file,origin,delta,npos)
    use rism_util, only : freeUnit
    implicit none
    character(len=*), intent(in) :: file
    _REAL_, intent(out) :: origin(:), delta(:)
    integer, intent(out) :: npos(:)
    integer :: unit, iostat
    unit = freeUnit()
    open(unit=unit, file=trim(file),status='old', iostat=iostat)
    if(iostat /= 0) then
       call rism_report_error("(a,i4)","readDXHeader: could not open "//trim(file)//":",iostat)
    end if
    call readDXHeader(unit,file,origin,delta,npos)
    close(unit)
  end subroutine readDXHeader_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Reads in a OpenDX file header.  Expects
!!!a 3D, regularly spaced grid in ASCII format.
!!!IN:
!!!   unit   :: open unit
!!!   file   :: file name of unit. Only for error reporting
!!!   origin :: real space coordinate of origin
!!!   delta  :: gridspacing
!!!   npos   :: number of grid points in each dimension
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readDXHeader_unit(unit,file,origin,delta,npos)
    implicit none
    integer, intent(in) :: unit
    character(len=*), intent(in) :: file
    _REAL_, intent(out) :: origin(:), delta(:)
    integer, intent(out) :: npos(:)
    character(len=1024) :: buffer
    !ncon : number of connections between data points in each dimension (should be the same as npos)
    !npt  : number of data points.  Should be product(npos)
    integer :: ncon(3), npt
    integer :: ipt, i,j,k, itemp
    integer :: iostat
    _REAL_ :: temp(3)

    !read header
    read(unit,'(a)') buffer
    if(buffer(1:len(obj1Label)) .ne. obj1Label)then
       call rism_report_error("first line in "//trim(file)//" must be:"&
            //obj1Label//" Nx Ny Nz")
    end if
    read(buffer(len(obj1Label)+1:),*)npos 

    read(unit,'(a)') buffer
    if(buffer(1:len(originLabel)) .ne. originLabel)then
       call rism_report_error("second line in "//trim(file)//" must be:"&
            //originLabel//" Ox Oy Oz")
    end if
    read(buffer(len(originLabel)+1:),*)origin 

    read(unit,'(a)') buffer
    if(buffer(1:len(deltaLabel)) .ne. deltaLabel)then
       call rism_report_error("third line in "//trim(file)//" must be:"&
            //deltaLabel//" dxx dxy dxz")
    end if
    read(buffer(len(deltaLabel)+1:),*)temp
    delta(1)=temp(1) 
    read(unit,'(a)') buffer
    if(buffer(1:len(deltaLabel)) .ne. deltaLabel)then
       call rism_report_error("fourth line in "//trim(file)//" must be:"&
            //deltaLabel//" dyx dyy dyz")
    end if
    read(buffer(len(deltaLabel)+1:),*)temp
    delta(2)=temp(2)
    read(unit,'(a)') buffer
    if(buffer(1:len(deltaLabel)) .ne. deltaLabel)then
       call rism_report_error("fifth line in "//trim(file)//" must be:"&
            //deltaLabel//" dzx dzy dzz")
    end if
    read(buffer(len(deltaLabel)+1:),*)temp
    delta(3)=temp(3) 

    read(unit,'(a)') buffer
    if(buffer(1:len(obj2Label)) .ne. obj2Label)then
       call rism_report_error("sixth line in "//trim(file)//" must be:"&
            //obj2Label//" Nx Ny Nz")
    end if
    read(buffer(len(obj2Label)+1:),*)ncon 

    read(unit,'(a)') buffer
    if(buffer(1:len(obj3Label)) .ne. obj3Label)then
       call rism_report_error("sixth line in "//trim(file)//" must be:"&
            //obj3Label//" Nx Ny Nz")
    end if
    read(buffer(len(obj3Label)+1:),*)npt 
    
    !check that the file is consistent
    if(npt /= product(npos))then
       write(0,'(a,i8,a,i8,a)') "number of data points (",npt,&
            ") is not equal to the number of positions (",product(npos),")"
    end if
  end subroutine readDXHeader_unit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Reads in a OpenDX file, returning data in a pointer, allocating necessary memory.  Expects
!!!a 3D, regularly spaced grid in ASCII format.
!!!IN:
!!!   file :: file name
!!!   origin :: real space coordinate of origin
!!!   delta  :: gridspacing
!!!OUT:
!!!    returns an Nx X Ny X Nz pointer of data.  You have to deallocate it.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function readDX_3D_p(file,origin,delta) result(data)
    use rism_util, only : freeUnit
    implicit none
    character(len=*), intent(in) :: file
    _REAL_, intent(out) :: origin(:), delta(:)
    _REAL_, pointer :: data(:,:,:)
    character(len=1024) :: buffer
    !npos : number of data points in each dimension
    integer :: npos(3)
    integer :: ipt, i,j,k, itemp
    integer :: unit, iostat
    _REAL_ :: temp(3)
    nullify(data)

    unit = freeUnit()
    open(unit=unit, file=file,status='old', iostat=iostat)
    if(iostat /= 0) then
       call rism_report_error("(a,i4)","could not open "//trim(file)//":",iostat)
    end if

    call readDXHeader(unit,file,origin,delta,npos)

    !allocate memory
    data=>safemem_realloc(data,npos(1),npos(2),npos(3),&
         .false.)

    call readDX_data(unit,data,npos,npos)
    close(unit)
  end function readDX_3D_p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Reads in a OpenDX file, writing data into an array.  Expects
!!!a 3D, regularly spaced grid in ASCII format.
!!!IN:
!!!   file   :: file name
!!!   data   :: an Nx X Ny X Nz array of data
!!!   origin :: real space coordinate of origin
!!!   delta  :: gridspacing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readDX_3D(file,data,origin,delta)
    use rism_util, only : freeUnit
    implicit none
    character(len=*), intent(in) :: file
    _REAL_, intent(out) :: origin(:), delta(:)
    _REAL_, intent(out) :: data(:,:,:)
    character(len=1024) :: buffer
    !npos : number of data points in each dimension
    integer :: npos(3)
    integer :: ipt, i,j,k, itemp
    integer :: unit, iostat
    _REAL_ :: temp(3)

    unit = freeUnit()
    open(unit=unit, file=file,status='old', iostat=iostat)
    if(iostat /= 0) then
       call rism_report_error("(a,i4)","could not open "//trim(file)//":",iostat)
    end if

    call readDXHeader(unit,file,origin,delta,npos)

    call readDX_data(unit,data,ubound(data),npos)
    close(unit)
  end subroutine readDX_3D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Reads in a OpenDX file, returning data in a pointer, allocating necessary memory.  Expects
!!!a 3D, regularly spaced grid in ASCII format.
!!!IN:
!!!   file :: file name
!!!   origin :: real space coordinate of origin
!!!   delta  :: gridspacing
!!!   npos   :: number of grid points in each dimension
!!!OUT:
!!!    returns an Nx * Ny * Nz (1D, column major) pointer of data.  You have to deallocate it.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function readDX_1D_p(file,origin,delta,npos) result(data)
    use rism_util, only : freeUnit
    implicit none
    character(len=*), intent(in) :: file
    _REAL_, intent(out) :: origin(:), delta(:)
    _REAL_, pointer :: data(:)
    integer, intent(out) :: npos(:)
    character(len=1024) :: buffer
    integer :: ipt, i,j,k, itemp
    integer :: unit, iostat
    _REAL_ :: temp(3)
    nullify(data)

    unit = freeUnit()
    open(unit=unit, file=file,status='old', iostat=iostat)
    if(iostat /= 0) then
       call rism_report_error("(a,i4)","could not open "//trim(file)//":",iostat)
    end if

    call readDXHeader(unit,file,origin,delta,npos)

    !allocate memory
    data=>safemem_realloc(data,product(npos),&
         .false.)
    call readDX_data(unit,data,npos,npos)
    close(unit)
  end function readDX_1D_p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Reads in a OpenDX file, writing data to an array.  Expects
!!!a 3D, regularly spaced grid in ASCII format.
!!!IN:
!!!   file :: file name
!!!   data   :: array of size product  Nx * Ny * Nz (1D, column major)
!!!   ndata  :: number of elements in each dimension in the data array
!!!   origin :: real space coordinate of origin
!!!   delta  :: gridspacing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readDX_1D(file,data,ndata,origin,delta)
    use rism_util, only : freeUnit
    implicit none
    character(len=*), intent(in) :: file
    integer, intent(in) :: ndata(3)
    _REAL_, intent(out) :: origin(:), delta(:)
    _REAL_, intent(in) :: data(:)
    character(len=1024) :: buffer
    !npos : number of data points in each dimension
    integer :: npos(3)
    integer :: ipt, i,j,k, itemp
    integer :: unit, iostat
    _REAL_ :: temp(3)

    unit = freeUnit()
    open(unit=unit, file=file,status='old', iostat=iostat)
    if(iostat /= 0) then
       call rism_report_error("(a,i4)","could not open "//trim(file)//":",iostat)
    end if

    call readDXHeader(unit,file,origin,delta,npos)

    call readDX_data(unit,data,ndata,npos)
    close(unit)
  end subroutine readDX_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!write to file in open DX format for use in VMD.  When writing in
!!!parallel, each process must call this function with its local
!!!data. Data transfer is handled internally.  We assume decomposition
!!!in the z-axis.
!!!IN:
!!!   file     :: file name to write to
!!!   data     :: data to write in a n(1)*n(2)*n(3) linear array
!!!   rism_box :: box dimensions in angstroms
!!!   n        :: number of points in each dimension of local data
!!!   nz_total :: total number of points in z dimension across all processes (== n(3) if single process)
!!!   ratucm   :: center of mass of the system
!!!   o_rank   :: (optional) mpi process rank
!!!   o_nproc  :: (optional) mpi number of processes
!!!   o_comm   :: (optional) mpi communicator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writeDX (file, data, rism_box, n,nz_total,ratucm, &
       o_rank,o_nproc,o_comm) 
    use rism_util, only : freeUnit, rmExPrec
    implicit none
#if defined(MPI)
    include 'mpif.h' 
#endif /*defined(MPI)*/
    character(len=*), intent(in) :: file
    integer,optional ::  o_rank,o_nproc,o_comm
    integer, intent(in) :: nz_total
    integer,intent(in) ::  n(3)
    _REAL_,intent(in) :: data(n(1),n(2),n(3)), rism_box(3),ratucm(3)
    integer :: rank=0,nproc=1,comm=0
    integer :: i,j,k, irank, err,count,icount
    integer, parameter :: dataperline=3
    integer :: unit, iostat
#ifdef MPI
    _REAL_,pointer :: z_data(:)=>NULL()
    integer, pointer :: nz_offset(:)=>NULL(),nz_local(:)=>NULL()
#endif /*MPI*/

#ifdef RISM_DEBUG

    write(0,*) "writeDX",rism_box,ratucm,rank
    call flush(6)
#endif /*RISM_DEBUG*/
    unit = freeUnit()
    if(present(o_rank)) rank = o_rank
    if(present(o_nproc)) nproc = o_nproc
    if(present(o_comm)) comm = o_comm
#ifdef MPI
    nz_offset => safemem_realloc(nz_offset,nproc,.false.)
    nz_local => safemem_realloc(nz_local,nproc,.false.)
#endif /*MPI*/
    if(rank==0)then
       open(unit=unit,file=file,iostat=iostat,status='REPLACE',form='FORMATTED')
       if(iostat /= 0)then
          call rism_report_error("opening "//trim(file))
       end if
    end if
    if(rank==0)then
       write(unit,"(a,i8,i8,i8)") "object 1 class gridpositions counts",n(1),n(2),nz_total
       write(unit,"(a,3(f15.8))") "origin ",ratucm(1)-rism_box(1)/2,ratucm(2)-rism_box(2)/2,ratucm(3)-rism_box(3)/2
       write(unit,"(a,f15.8,a)") "delta  ",rism_box(1)/n(1)," 0 0"
       write(unit,"(a,f15.8,a)") "delta  0 ",rism_box(2)/n(2)," 0"
       write(unit,"(a,f15.8,a)") "delta  0 0 ",rism_box(3)/nz_total
       write(unit,"(a,i9,i9,i9)") "object 2 class gridconnections counts"&
            ,n(1),n(2),nz_total
       write(unit,"(a,1x,i27,1x,a)") &
            "object 3 class array type double rank 0 items", &
            n(1)*n(2)*nz_total, "data follows"
    end if
    count = 0
#if defined(MPI)
    z_data =>safemem_realloc(z_data,nz_total,.false.)
    call mpi_gather(n(3),1,mpi_integer,nz_local,1,mpi_integer,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_OPENDX: could not gather N")
    nz_offset(1) = 0
    do i = 2, nproc
       nz_offset(i) =sum( nz_local(1:i-1))
    end do
#endif /*defined(MPI)*/
    do i=1,n(1)
       do j=1,n(2)
#if defined(MPI)
          call mpi_gatherv(data(i,j,:),n(3),mpi_double_precision,&
               z_data,nz_local,nz_offset,mpi_double_precision,&
               0,comm,err)
          if(err /=0) call rism_report_error&
               ("RISM3D_OPENDX: could not gather DATA")
#endif /*defined(MPI)*/
          if(rank==0)then
             do k=1,nz_total
#if defined(MPI)
                write(unit,"(E16.5E3)",advance="no") rmExPrec(z_data(k))
#else
                write(unit,"(E16.5E3)",advance="no") rmExPrec(data(i,j,k))
#endif /*defined(MPI)*/
                count = count +1
                if(mod(count,dataperline) == 0)then
                   write(unit,"()")
                   count = 0;
                end if
             end do
          end if
       end do
    end do
    if(rank == 0)then
       if(mod(count,dataperline) /= 0)then
          write(unit,"()")
       end if
       write(unit,"(a)") 'object "Untitled" class field'
       close(unit)
    end if
#ifdef MPI
    if(safemem_dealloc(nz_offset)/=0) call rism_report_error("WRITEDX: nz_offset deallocation failed")
    if(safemem_dealloc(nz_local)/=0) call rism_report_error("WRITEDX: nz_local deallocation failed")
    if(safemem_dealloc(z_data)/=0) call rism_report_error("WRITEDX: z_data deallocation failed")
#endif /*MPI*/
  end subroutine writeDX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!            PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Reads in data from OpenDX file, into pre-allocated memory.  Expects
!!!a 3D, regularly spaced grid in ASCII format.
!!!IN:
!!!   unit   :: open unit
!!!   data   :: array of size product(ndata)
!!!   ndata  :: number of elements in each dimension in the data array
!!!   npos   :: number of elements in each dimension in the file
!!!OUT:
!!!    fills an Nx X Ny X Nz array of data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readDX_data(unit,data,ndata,npos)
    use rism_util, only : freeUnit
    implicit none
    integer,intent(in) :: unit, npos(3), ndata(3)
    _REAL_ :: data(ndata(1),ndata(2),ndata(3))
    character(len=1024) :: buffer
    !npos : number of data points in each dimension
    integer :: ipt, i,j,k, itemp
    integer :: iostat
    _REAL_ :: temp(3)

    !read data
    i=1
    j=1
    k=1
    data=huge(1d0)
    do ipt = 1, product(npos)/3
       read(unit,*) temp
       do itemp=1,3
          data(i,j,k)=temp(itemp)
          !increament array counters
          k = k+1
          j = j+(k-1)/npos(3)
          k = mod(k-1,npos(3))+1
          i = i+(j-1)/npos(2)
          j = mod(j-1,npos(2))+1
       end do
    end do
    !handle the last line special
    if(k /= 1)then
       read(unit,*) data(i,j,k:npos(3))
    end if
    !Check for unassigned values
    if(count(data/=huge(1d0)) /= product(npos))then
       call rism_report_error("OpenDX data incompletely read.  Is your file complete?")
    end if
  end subroutine readDX_data

end module rism3d_opendx
