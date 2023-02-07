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
!Reads/writes non-portable binary restart files for 3D-RISM
!calculations. The file format is Fortran unformatted binary and has the form
!
!ngr natv
!tuv(ngr,natv)
!
!Simply put, ngr and natv (integers) are the number of grid points and number of 
!solvent types, giving the size of the density/correlation array tuv.
!This format is not portable and not MPI friendly in anyway and is deprecated.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module rism3d_csv
  use rism_report_c
  implicit none
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Reading unformatted array
!!!IN
!!!  file : file name
!!!  tuv  : correlation function
!!!  ngr  : expected number of data points for each solvent type
!!!  natv : expected number of solvent types 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  rdufma (file, tuv,ngr,natv)
    implicit none
    character(*), intent(in) :: file
    integer,intent(in) ::  ngr,natv
    _REAL_,intent(out) ::  tuv(ngr,natv)
    integer :: ncha=99,iostat,i,j
    integer :: ngro,natvo
    open (ncha,file=file,form='unformatted',status='old',iostat=iostat,position="REWIND")
    if(iostat /= 0) then
       call rism_report_error("on "//file//" open")
    endif
    read (ncha,iostat=iostat)  ngro, natvo
    if(iostat /= 0) then
       call rism_report_error("on "//file//" read")
    endif
    if (ngro /= ngr .OR. natvo /= natv)  then
       close (ncha,iostat=iostat)
       call rism_report_error('(a,i6,a,i3,a,i6,a,i3)','RXRISM:  actual   NGr = ',&
            ngr, ', NatV = ', natv,'; saved   NGr = ',ngro,', NatV = ',natvo)
    endif
    do i=1,natv
       read (ncha,iostat=iostat)  tuv(1:ngr,i)
       if(iostat /= 0) then
          call rism_report_error('(a,i4)','RDUFMA: I/O error:',iostat)
       endif
    end do
    close (ncha,iostat=iostat)
    if(iostat /= 0) then
       call rism_report_error("on "//file//" close")
    endif
  end subroutine rdufma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Reading unformatted array
!!!IN
!!!  file : file name
!!!  tuv  : correlation function
!!!  ngr  : number of data points for each solvent type
!!!  natv : number of solvent types 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  wrufma (file, tuv,ngr,natv)
    implicit none
    character(*), intent(in) :: file
    integer,intent(in) ::   ngr,natv
    _REAL_,intent(in) ::  tuv(ngr,natv)
    integer :: ncha=99,iostat,i
    open (ncha,file=file,form='unformatted',iostat=iostat,status="REPLACE")
    if(iostat /= 0) then
       call rism_report_error("on "//file//" open")
    endif
    write (ncha,iostat=iostat)  ngr, natv
    if(iostat /= 0) then
       call rism_report_error("on "//file//" write")
    endif
    do i=1,natv
       write (ncha,iostat=iostat)  tuv(1:ngr,i)
       if(iostat /= 0) then
          call rism_report_error('(a,i4)','WRUFMA:  I/O error:',iostat)
       endif
    end do
    close (ncha,iostat=iostat)
    if(iostat /= 0) then
       call rism_report_error("on "//file//" close")
    endif
  end subroutine wrufma
  
end module rism3d_csv
