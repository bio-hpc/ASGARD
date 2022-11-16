!<compile=optimized>

!The 3D-RISM-KH software found here is copyright (c) 2011-2012 by Tyler Luchko 
!and David A. Case.
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

!Collection of utility functions for arrays

module array_util
  implicit none
  interface array_index
     module procedure index_c, index_i, index_r8!, index_l
  end interface
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the first element of the array that matches the given value.  If 
!!!'back' then returns the last element.
!!!IN:
!!!   a     : character array
!!!   value : character string
!!!   back  : (optional) search backwards
!!!OUT:
!!!    Index of the element matching the value, -1 if it does not exist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function index_c(a,value,back) result(i)
    implicit none
    character(len=*), intent(in) :: a(:), value
    logical, optional, intent(in) :: back
    integer :: i, start, finish, incr
    start = lbound(a,1)
    finish = ubound(a,1)
    incr = 1
    if(present(back))then
       if(back)then
          start = ubound(a,1)
          finish = lbound(a,1)
          incr = -1
       end if
    end if
    do i=start,finish,incr
       if(trim(a(i)) .eq. trim(value)) return
    end do
    i = -1
  end function index_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the first element of the array that matches the given value.  If 
!!!'back' then returns the last element.
!!!IN:
!!!   a     : integer array
!!!   value : integer value
!!!   back  : (optional) search backwards
!!!OUT:
!!!    Index of the element matching the value, -1 if it does not exist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function index_i(a,value,back) result(i)
    implicit none
    integer, intent(in) :: a(:), value
    logical, optional, intent(in) :: back
    integer :: i, start, finish, incr
    start = lbound(a,1)
    finish = ubound(a,1)
    incr = 1
    if(present(back))then
       if(back)then
          start = ubound(a,1)
          finish = lbound(a,1)
          incr = -1
       end if
    end if
    do i=start,finish,incr
       if(a(i) == value) return
    end do
    i = -1
  end function index_i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the first element of the array that matches the given value.  If 
!!!'back' then returns the last element.
!!!IN:
!!!   a     : real*8 array
!!!   value : real*8 value
!!!   back  : (optional) search backwards
!!!OUT:
!!!    Index of the element matching the value, -1 if it does not exist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function index_r8(a,value,back) result(i)
    implicit none
    real*8, intent(in) :: a(:), value
    logical, optional, intent(in) :: back
    integer :: i, start, finish, incr
    start = lbound(a,1)
    finish = ubound(a,1)
    incr = 1
    if(present(back))then
       if(back)then
          start = ubound(a,1)
          finish = lbound(a,1)
          incr = -1
       end if
    end if
    do i=start,finish,incr
       if(a(i) == value) return
    end do
    i = -1
  end function index_r8
end module array_util
