!<compile=optimized>

!The volslice software found here is copyright (c) 2012 by 
!Tyler Luchko and David A. Case.
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

!program to take an arbitrary 2D slice of data out of a 3D volume
!Initially there are many restrictions.

!Input: only supports regular grids
!       only supports opendx format

!Output: only supports regular grids
!        only supports gnuplot format

!Slice: limited to xy, yz or zx

module vs
  use rism3d_opendx
  implicit none
  character(len=256) :: plane='', input='', output=''
  integer :: index=huge(1)

  !dxOrigin :: origin from the DX file.  Tells us how to translate the solute and
  !            how to interpret volume specifications
  !origin :: origin for integration
  !vec    :: orientation vector for polar averaging
  !radius :: radial extent for averaging
  !zlen   :: z-axis extent for averaging.  Centered on origin.
  !dr     :: radial spacing
  !dz     :: z-axis spacing for polar averaging
  _REAL_ :: dxOrigin(3),delta(3)
  integer :: npos(3)

  _REAL_, pointer :: data(:,:,:)=>NULL()
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Get user options
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getoptions()
    use getopts_c
    implicit none
    integer :: err
    character(len=1024),pointer :: extra(:)=>NULL()
    call getopts_add("plane",plane,max=1,min=1)
    call getopts_add("index",index,max=1,min=1)
    err = getopts_process()
    if(err <0)then
       write(0,'(a,i2)') "ERROR: invalid options:",err
       call usage()
       stop
    end if
    call getopts_get("plane",1,plane)
    if(len_trim(plane) /= 2) then
       write(0,'(a,a,i2)') "ERROR: --plane must be a two character combination of x,y or z", plane, len_trim(plane)
       call usage()
    end if
    if(plane .eq. "yx") plane = "xy"
    if(plane .eq. "zy") plane = "yz"
    if(plane .eq. "xz") plane = "zx"
    call getopts_get("index",1,index)
    if(index == huge(1) .or. index <= 0) then
       write(0,'(a)') "ERROR: a positive index must be specified"
       call usage()
    end if
    extra=> getopts_unread()
    if(size(extra) /= 2)then
       write(0,'(a)') "ERROR: an input and output file must be specified on the command line"
       call usage()
    end if
    input = trim(extra(1))
    output = trim(extra(2))
  end subroutine getoptions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Read volume file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readfile()
    implicit none
    data=> readDX_p(input,dxOrigin,delta)
    npos = ubound(data)
    if(plane .eq. "xy" .and. index > npos(3)) then
       write(0,'(a)') "ERROR: index larger than indices in z dimension"
       call usage()
    elseif(plane .eq. "yz" .and. index > npos(1)) then
       write(0,'(a)') "ERROR: index larger than indices in x dimension"
       call usage()
    elseif(plane .eq. "zx" .and. index > npos(2)) then
       write(0,'(a)') "ERROR: index larger than indices in y dimension"
       call usage()
    end if
  end subroutine readfile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Writes out the volue slice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine writeslice()
    implicit none
    integer :: unit = 99, iostat
    integer :: i,j,k 
    open(unit=unit, file=output,status='replace', iostat=iostat)
    if(plane .eq. "xy")then
       do i=1,npos(1)
          do j=1,npos(2)
             write(unit,'(1p,e16.8e3,1x,e16.8e3,1x,e16.8e3)') dxorigin(1)+(i-1)*delta(1),&
                  dxorigin(2)+(j-1)*delta(2),&
                  data(i,j,index)
          end do
          write(unit,'(a)')
       end do
    end if
    close(unit)
  end subroutine writeslice

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Usage message
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine usage()
    implicit none
    write(0,'(a)') "USAGE:  volslice --plane plane --index index input.dx output.dat"
    write(0,'(a)') "--plane                  xy, yz, or zx"
    write(0,'(a)') "--index                  index of the dimension perpendicular to the plane"
    write(0,'(a)') "input.dx                3D OpenDX file"
    write(0,'(a)') "output.dat              2D Gnuplot file"
    stop
  end subroutine usage
end module vs

program volslice
  use vs
  implicit none
  call getoptions()
  call readfile()
  call writeslice()
!  call cleanup()
end program volslice
