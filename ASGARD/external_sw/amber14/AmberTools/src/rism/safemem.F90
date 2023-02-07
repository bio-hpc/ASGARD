!<compile=optimized>

!The 3D-RISM-KH software found here is copyright (c) 2010-2012 by
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Functions and subroutines to safely allocate, reallocate (data preserving)
!!!and deallocate pointers. 
!!!
!!!Typical usage example:
!!!
!!!  _REAL_, pointer :: p(:)=>NULL()
!!!  call safemem_realloc(p,N)
!!!
!!!  !somewhere else, not necessarily in the same scope ('p' may be a
!!!  !different pointer handle here), resize the array.
!!!  call safemem_realloc(p,M)
!!!  
!!!  !somewhere else, not necessarily in the same scope ('p' may be a
!!!  !different pointer handle here), deallocate the memory.
!!!  if(safemem_dealloc(p) /= 0) then
!!!    !deal with error
!!!  endif
!!!
!!!  !at any point you can get a summary of the current and maximum
!!!  !allocated memory in bytes
!!!  integer(kind=8) :: memstats(10)
!!!  memstats = memStatus()
!!!  write(unit,'(a)') "Type         Current         Maximum"
!!!  write(unit,'(a,i12,a,f12.5,a)') "Integer  ",memstats(1)," B ",&
!!!       dble(memstats(6))/BYTES_PER_GB," GB"
!!!  write(unit,'(a,i12,a,f12.5,a)') "Real     ",memstats(2)," B ",&
!!!       dble(memstats(7))/BYTES_PER_GB," GB"
!!!  write(unit,'(a,i12,a,f12.5,a)') "Logical  ",memstats(3)," B ",&
!!!       dble(memstats(8))/BYTES_PER_GB," GB"
!!!  write(unit,'(a,i12,a,f12.5,a)') "Character",memstats(4)," B ",&
!!!       dble(memstats(9))/BYTES_PER_GB," GB"
!!!  write(unit,'(a)') "---------------------------------------"
!!!  write(unit,'(a,i12,a,f12.5,a)') "Total    ",memstats(5)," B ",&
!!!       dble(memstats(10))/BYTES_PER_GB," GB"
!!!
!!!Notes on proper usage:
!!!-safemem_realloc preserves data by default.  Since Fortran has no
!!! built in facility for this, it is generally slow.  If you don't
!!! need this feature, use 'preserve=.false.'.
!!!-Uninitiallized pointers should be set to NULL.  Either use
!!! 'p=>NULL()' or 'nullify(p)'.
!!!-As long as you are careful to set unitialized data to NULL, you
!!! can mix these routine calls with your own allocated/deallocate
!!! statements on any pointer.  However, the memory summary will be
!!! invalid.
!!!-If compiled and linked with FFTW, 16-byte aligned memory, useful
!!! for SIMD operations can be allocated by using
!!! 'o_aligned=.true.'. If you allocate aligned memory, you _must_
!!! deallocate it with the aligned option.  It is not possible for
!!! these routines to determine how the memory was allocated.
!!!-The memory summary only tracks heap memory allocated with these
!!! routines.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module safemem
    use rism_report_c
    use FFTW3
    implicit none
    public
    type memTracker
       private
       !int     :: Current number of bytes of integer memory allocated
       !real    :: Current number of bytes of real memory allocated
       !logical :: Current number of bytes of logical memory allocated
       !char    :: Current number of bytes of character memory allocated
       !total   :: Current number of bytes of total memory allocated
       integer(8) :: int=0,real=0,logical=0,char=0,total=0
       !maxint     :: Current number of bytes of integer memory allocated
       !maxreal    :: Current number of bytes of real memory allocated
       !maxlogical :: Current number of bytes of logical memory allocated
       !maxchar    :: Current number of bytes of character memory allocated
       !max        :: Current number of bytes of total memory allocated
       integer(8) :: maxint=0,maxreal=0,maxlogical=0,maxchar=0,max=0
    end type memTracker

    !BYTES_PER_GIGABYTES :: used to convert between bytes and GB
    integer(8),parameter :: BYTES_PER_GB = 1024**3
    !BYTES_PER_MEGABYTES :: used to convert between bytes and MB
    integer(8),parameter :: BYTES_PER_MB = 1024**2
    !BYTES_PER_KILOBYTES :: used to convert between bytes and KB
    integer(8),parameter :: BYTES_PER_KB = 1024**1

    type(memTracker),private,save :: totalMem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! interfaces
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!reallocate memory of a specific pointer during the simulation. This can 
!!!handle up to 5 dimensions and the optiontional last argument can indicate 
!!!whether or not to preserve the contents of the array.  Preservation is the 
!!!default.  Intrinsic function associated is applied to the first argument.
!!!IN:
!!!   p     :: original well defined pointer; the pointer association status
!!!            must be associated or disassociated, ie, nullified, but cannot
!!!            be undefined.
!!!   n1    :: size of the first dimension
!!!   n2    :: size of the second dimension (optional)
!!!   n3    :: size of the third dimension (optional)
!!!   n4    :: size of the fourth dimension (optional)
!!!   n5    :: size of the fifth dimension (optional)
!!!   preserve :: boolean indicating whether or not to preserve the contents of 
!!!           the array (optional)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTERFACE safemem_realloc
       MODULE PROCEDURE safemem_realloc_1d_real
       MODULE PROCEDURE safemem_realloc_2d_real
       MODULE PROCEDURE safemem_realloc_3d_real
       MODULE PROCEDURE safemem_realloc_4d_real
       MODULE PROCEDURE safemem_realloc_5d_real
       MODULE PROCEDURE safemem_realloc_5d_real_interp

       MODULE PROCEDURE safemem_realloc_1d_integer
       MODULE PROCEDURE safemem_realloc_2d_integer

       MODULE PROCEDURE safemem_realloc_1d_character
       MODULE PROCEDURE safemem_realloc_2d_character

       MODULE PROCEDURE safemem_realloc_1d_logical
    END INTERFACE safemem_realloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!checks if the pointer is allocated before deallocating.  Returns any error code  
!!!Intrinsic function associated is applied to the first argument.
!!!IN:
!!!   p :: well defined pointer to be deallocated; the pointer association status
!!!        must be associated or disassociated, ie, nullified, but cannot
!!!        be undefined.
!!!OUT:
!!!   error code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTERFACE safemem_dealloc
       MODULE PROCEDURE safemem_dealloc_pointer_1d_real
       MODULE PROCEDURE safemem_dealloc_pointer_2d_real
       MODULE PROCEDURE safemem_dealloc_pointer_3d_real
       MODULE PROCEDURE safemem_dealloc_pointer_4d_real
       MODULE PROCEDURE safemem_dealloc_pointer_5d_real
       MODULE PROCEDURE safemem_dealloc_pointer_1d_int
       MODULE PROCEDURE safemem_dealloc_pointer_2d_int
       MODULE PROCEDURE safemem_dealloc_pointer_3d_int
       MODULE PROCEDURE safemem_dealloc_pointer_4d_int
       MODULE PROCEDURE safemem_dealloc_pointer_5d_int
       MODULE PROCEDURE safemem_dealloc_pointer_1d_character
       MODULE PROCEDURE safemem_dealloc_pointer_2d_character
       MODULE PROCEDURE safemem_dealloc_pointer_1d_logical
    end INTERFACE safemem_dealloc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns an array of current and maximum memory allocated in bytes so
!!!far.  No distinction is made between different 
!!!kinds of the same type:
!!!(1) : Current integer memory
!!!(2) : Current real memory
!!!(3) : Current logical memory
!!!(4) : Current char memory
!!!(5) : Current logical memory
!!!(6) : Max integer memory
!!!(7) : Max real memory
!!!(8) : Max logical memory
!!!(9) : Max char memory
!!!(10): Max logical memory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    interface memStatus
       module procedure memStatus_obj, memStatus_global
    end interface memStatus


    public safemem_realloc, safemem_dealloc, memStatus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!PRIVATE METHODS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    private memadd_i,memadd_r,memadd_c,memadd_l

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Given a value x,y,z and arrays of coordinates x1, x2, x3, this function locates the index of 
!!!x1,2,3 that is closest to x, y, z
!!!IN:
!!!   x1 :: values in dimension 1 (ascending order)
!!!   x2 :: values in dimension 2 (ascending order)
!!!   x3 :: values in dimension 3 (ascending order)
!!!   n1 :: number of values
!!!   n2 :: number of values
!!!   n3 :: number of values
!!!   x  :: x value
!!!   y  :: y value
!!!   z  :: z value
!!!OUT:
!!!    three 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function nearestIndex(x1,x2,x3,n1,n2,n3,x,y,z)
      implicit none
      integer, intent(in) :: n1, n2, n3
      _REAL_, intent(in) :: x1(n1),x2(n2),x3(n3),x,y,z
      integer :: nearestIndex(3)
      integer :: i
      nearestIndex = 1
      do i = 2, n1
         if(abs(x1(i)-x) < abs(x1(nearestindex(1))-x))then
            nearestIndex(1) = i
         end if
      end do
      do i = 2, n2
         if(abs(x2(i)-y) < abs(x2(nearestindex(2))-y))then
            nearestIndex(2) = i
         end if
      end do
      do i = 2, n3
         if(abs(x3(i)-z) < abs(x3(nearestindex(3))-z))then
            nearestIndex(3) = i
         end if
      end do
      return
    end function nearestIndex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!reallocate memory of a specific pointer during the simulation. This can 
!!!handle 1 dimension and preserves the contents of the array.
!!!IN:
!!!   p     :: original well defined pointer; the pointer association status
!!!            must be associated or disassociated, ie, nullified, but cannot
!!!            be undefined.
!!!   n1    :: size of the first dimension
!!!   o_preserve :: (optional) boolean indicating whether or not to
!!!                 preserve the contents of the array
!!!   o_aligned :: (optional) Use the FFTW library to allocate 16-byte
!!!                aligned memory necessary to SIMD operations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION safemem_realloc_1d_real(p, n1,o_preserve,o_aligned)
      implicit none
      _REAL_, POINTER, DIMENSION(:) :: safemem_realloc_1d_real
      _REAL_, POINTER, DIMENSION(:):: p
      logical ,optional, intent(in) :: o_preserve,o_aligned
      INTEGER, intent(in) :: n1
      logical :: preserve,aligned
      INTEGER :: i,nold1, err

      type(C_PTR) :: cptr=C_NULL_PTR

#ifdef RISM_DEBUG
      write(6,*) "safemem_realloc_1d_real", associated(p),n1,preserve
      call flush(6)
#endif /*RISM_DEBUG*/    
      aligned = .false.
      if(present(o_aligned)) aligned = o_aligned

      preserve = .true.
      if(present(o_preserve)) preserve = o_preserve

      if(.not.preserve)then
         if(safemem_dealloc(p,aligned) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 1D REAL")
      end if
      !get new memory
      if(aligned)then
         err = 0
         cptr = fftw_alloc_real(int(n1, C_SIZE_T))
         if(.not.c_associated(cptr)) then
            err = 1
         else
            call c_f_pointer(cptr,safemem_realloc_1d_real,[n1])
         end if
      else
         ALLOCATE(safemem_realloc_1d_real(n1), STAT=err)
      end if
      IF(err /= 0)&
           call rism_report_error("reallocation_pointer allocate error: 1D REAL")
      call memadd_r(totalmem,product(ubound(safemem_realloc_1d_real))*kind(p))

      !transfer data to the new memory
      IF(ASSOCIATED(p) .and. preserve .and. n1 >0)then
         nold1 = MIN(UBOUND(p,1), n1)
         call dcopy(nold1, p(1), 1, &
              safemem_realloc_1d_real(1),1)
         safemem_realloc_1d_real(1:nold1) = p(1:nold1)
         !free old memory
         if(safemem_dealloc(p,aligned) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 1D REAL")
      end if
    END FUNCTION SAFEMEM_REALLOC_1D_REAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!reallocate memory of a specific pointer during the simulation. This can 
!!!handle 2 dimensions and preserves the contents of the array.
!!!IN:
!!!   p     :: original well defined pointer; the pointer association status
!!!            must be associated or disassociated, ie, nullified, but cannot
!!!            be undefined.
!!!   n1    :: size of the first dimension
!!!   n2    :: size of the second dimension
!!!   o_preserve :: (optional) boolean indicating whether or not to
!!!                 preserve the contents of the array
!!!   o_aligned :: (optional) Use the FFTW library to allocate 16-byte
!!!                aligned memory necessary to SIMD operations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION safemem_realloc_2d_real(p,n1, n2,o_preserve,o_aligned)
      implicit none
      _REAL_, POINTER, DIMENSION(:,:) :: safemem_realloc_2d_real
      _REAL_, POINTER, DIMENSION(:,:) :: p
      logical,optional, intent(in) :: o_preserve,o_aligned
      INTEGER, intent(in) :: n1,n2
      logical :: preserve, aligned
      INTEGER :: nold1, nold2, err
      integer :: i2
      type(C_PTR) :: cptr=C_NULL_PTR

#ifdef RISM_DEBUG
      write(6,*) "safemem_realloc_2d_real", associated(p),n1,n2,preserve
      call flush(6)
#endif /*RISM_DEBUG*/

      aligned = .false.
      if(present(o_aligned)) aligned = o_aligned

      preserve=.true.
      if(present(o_preserve)) preserve = o_preserve

      if(.not.preserve)then
         if(safemem_dealloc(p,aligned) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 2D REAL 1")
      end if
      if(aligned)then
         err = 0
         cptr = fftw_alloc_real(int(n1*n2, C_SIZE_T))
         if(.not.c_associated(cptr)) then
            err = 1
         else
            call c_f_pointer(cptr,safemem_realloc_2d_real,[n1,n2])
         end if
      else
         ALLOCATE(safemem_realloc_2d_real(1:n1,1:n2), STAT=err)
      end if
      IF(err /= 0)&
         call rism_report_error("(a,i4)","reallocation_pointer allocate error: 2D REAL 2",err)
      call memadd_r(totalmem,product(ubound(safemem_realloc_2d_real))*kind(p))
      if(ASSOCIATED(p) .and.preserve .and. n1*n2 > 0)then
         safemem_realloc_2d_real = 0d0
         nold1 = MIN(UBOUND(p,1), n1)
         nold2 = MIN(UBOUND(p,2), n2)
         do i2=1,nold2
            call dcopy(nold1, p(1, i2), 1, &
                 safemem_realloc_2d_real(1, i2),1) 
         end do
         if(safemem_dealloc(p,aligned) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 2D REAL 3")
      end if
    END FUNCTION SAFEMEM_REALLOC_2D_REAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!reallocate memory of a specific pointer during the simulation. This can 
!!!handle 3 dimension and can preserve the contents of the array.
!!!IN:
!!!   p     :: original well defined pointer; the pointer association status
!!!            must be associated or disassociated, ie, nullified, but cannot
!!!            be undefined.
!!!   n1    :: size of the first dimension
!!!   n2    :: size of the second dimension
!!!   n3    :: size of the third dimension
!!!   o_preserve :: (optional) boolean indicating whether or not to
!!!                 preserve the contents of the array
!!!   o_aligned :: (optional) Use the FFTW library to allocate 16-byte
!!!                aligned memory necessary to SIMD operations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION safemem_realloc_3d_real(p, n1, n2, n3, o_preserve,o_aligned)
      implicit none
      _REAL_, POINTER, DIMENSION(:,:,:) :: safemem_realloc_3d_real
      _REAL_, POINTER, DIMENSION(:,:,:):: p
      logical, optional, intent(in) :: o_preserve,o_aligned
      INTEGER, intent(in) :: n1, n2, n3
      logical :: preserve,aligned
      INTEGER :: nold1, nold2, nold3, err
      integer :: i2, i3
      type(C_PTR) :: cptr=C_NULL_PTR

#ifdef RISM_DEBUG
      write(6,*) "safemem_realloc_3d_real", associated(p),n1,n2,n3,preserve
      call flush(6)
#endif /*RISM_DEBUG*/
      aligned = .false.
      if(present(o_aligned)) aligned = o_aligned

      preserve = .true.
      if(present(o_preserve)) preserve = o_preserve

      if(.not.preserve)then
         if(safemem_dealloc(p,aligned) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 3D REAL")
      end if
      if(aligned)then
         err = 0
         cptr = fftw_alloc_real(int(n1*n2*n3, C_SIZE_T))
         if(.not.c_associated(cptr)) then
            err = 1
         else
            call c_f_pointer(cptr,safemem_realloc_3d_real,[n1,n2,n3])
         end if
      else
         ALLOCATE(safemem_realloc_3d_real(1:n1,1:n2,1:n3), STAT=err)
      end if
      IF(err /= 0)&
         call rism_report_error("(a,1x,i4)", "reallocation_pointer allocate error: 3D REAL: ",err)
      call memadd_r(totalmem,product(ubound(safemem_realloc_3d_real))*kind(p))

      if(ASSOCIATED(p) .and. preserve .and. n1*n2*n3 > 0)then 
         safemem_realloc_3d_real = 0d0
         nold1 = MIN(UBOUND(p,1), n1)
         nold2 = MIN(UBOUND(p,2), n2)
         nold3 = MIN(UBOUND(p,3), n3)
         do i3=1,nold3
            do i2=1,nold2
               call dcopy(nold1, p(1, i2, i3), 1, &
                    safemem_realloc_3d_real(1, i2, i3),1)
            end do
         end do
         if(safemem_dealloc(p,aligned) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 3D REAL")
      end if
    END FUNCTION SAFEMEM_REALLOC_3D_REAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!reallocate memory of a specific pointer during the simulation. This can 
!!!handle 4 dimension and can preserve the contents of the array.
!!!IN:
!!!   p     :: original well defined pointer; the pointer association status
!!!            must be associated or disassociated, ie, nullified, but cannot
!!!            be undefined.
!!!   n1    :: size of the first dimension
!!!   n2    :: size of the second dimension
!!!   n3    :: size of the third dimension
!!!   n4    :: size of the fourth dimension
!!!   o_preserve :: (optional) boolean indicating whether or not to
!!!                 preserve the contents of the array
!!!   o_center :: (optional) when preserving the contents of the
!!!               array, center the data in the array instead of
!!!               preserving the value at each index
!!!   o_aligned :: (optional) Use the FFTW library to allocate 16-byte
!!!                aligned memory necessary to SIMD operations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION safemem_realloc_4d_real(p, n1, n2, n3, n4, o_preserve, o_center,o_aligned)
      implicit none
      _REAL_, POINTER, DIMENSION(:,:,:,:) :: safemem_realloc_4d_real
      _REAL_, POINTER, DIMENSION(:,:,:,:):: p
      logical,optional, intent(in) :: o_preserve, o_center,o_aligned
      INTEGER, intent(in) :: n1, n2, n3, n4
      logical :: preserve, center,aligned
      INTEGER :: nold1, nold2, nold3, nold4, err, offset(4)
      integer :: i2,i3,i4
      type(C_PTR) :: cptr=C_NULL_PTR

#ifdef RISM_DEBUG
      write(6,*) "safemem_realloc_4d_real", associated(p),n1,n2,n3,n4,preserve
      call flush(6)
#endif /*RISM_DEBUG*/
      aligned = .false.
      if(present(o_aligned)) aligned = o_aligned

      preserve = .true.
      if(present(o_preserve)) preserve=o_preserve

      center = .false.
      if(present(o_center)) center=o_center

      if(.not.preserve)then
         if(safemem_dealloc(p,aligned) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 4D REAL")
      end if
      if(aligned)then
         err = 0
         cptr = fftw_alloc_real(int(n1*n2*n3*n4, C_SIZE_T))
         if(.not.c_associated(cptr)) then
            err = 1
         else
            call c_f_pointer(cptr,safemem_realloc_4d_real,[n1,n2,n3,n4])
         end if
      else
         ALLOCATE(safemem_realloc_4d_real(n1,n2,n3,n4), STAT=err)
      end if
      IF(err /= 0)&
         call rism_report_error("reallocation_pointer allocate error: 4D REAL")
      call memadd_r(totalmem,product(ubound(safemem_realloc_4d_real))*kind(p))

      if(ASSOCIATED(p) .and. preserve .and. n1*n2*n3*n4 > 0)then
         safemem_realloc_4d_real = 0d0
         nold1 = MIN(UBOUND(p,1), n1)
         nold2 = MIN(UBOUND(p,2), n2)
         nold3 = MIN(UBOUND(p,3), n3)
         nold4 = MIN(UBOUND(p,4), n4)
         if(center)then
            offset(1) = (UBOUND(p,1)- n1)/2
            offset(2) = (UBOUND(p,2)- n2)/2
            offset(3) = (UBOUND(p,3)- n3)/2
            offset(4) = (UBOUND(p,4)- n4)/2
         else
            offset=0
         end if
         do i4=1,nold4
            do i3=1,nold3
               do i2=1,nold2
                  call dcopy(nold1, &
                       p(max(offset(1),0)+1, &
                       max(offset(2),0)+i2, &
                       max(offset(3),0)+i3, &
                       max(offset(4),0)+i4),1, &
                       safemem_realloc_4d_real(&
                       abs(min(offset(1),0))+1,&
                       abs(min(offset(2),0))+i2,&
                       abs(min(offset(3),0))+i3,&
                       abs(min(offset(4),0))+i4),1)
               end do
            end do
         end do
         if(safemem_dealloc(p,aligned) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 4D REAL")
      end if
    END FUNCTION SAFEMEM_REALLOC_4D_REAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!reallocate memory of a specific pointer during the simulation. This can 
!!!handle 5 dimensions. The original contents of the array are 3D polynomial 
!!!interpolated to fill the new array.  That is, only the first three 
!!!dimensions are interpolated
!!!IN:
!!!   p     :: original well defined pointer; the pointer association status
!!!            must be associated or disassociated, ie, nullified, but cannot
!!!            be undefined.
!!!   n1    :: size of the first dimension
!!!   n2    :: size of the second dimension
!!!   n3    :: size of the third dimension
!!!   n4    :: size of the fourth dimension
!!!   n5    :: size of the fifth dimension
!!!   xi1   :: input linear values of the first dimension grid points
!!!   xi2   :: input linear values of the second dimension grid points
!!!   xi2   :: input linear values of the third dimension grid points
!!!   xo1   :: input linear values of the first dimension grid points
!!!   xo2   :: input linear values of the second dimension grid points
!!!   xo2   :: input linear values of the third dimension grid points
!!!   o_aligned :: (optional) Use the FFTW library to allocate 16-byte
!!!                aligned memory necessary to SIMD operations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION safemem_realloc_5d_real_interp(p, n1, n2, n3, n4, n5,&
         xi1,xi2,xi3,xo1,xo2,xo3,o_aligned)
      implicit none
      _REAL_, POINTER, DIMENSION(:,:,:,:,:) :: safemem_realloc_5d_real_interp
      _REAL_, POINTER, DIMENSION(:,:,:,:,:):: p
      _REAL_, DIMENSION(:), intent(in) :: xi1,xi2,xi3,xo1,xo2,xo3
      logical ,optional, intent(in) :: o_aligned
      INTEGER, intent(in) :: n1, n2, n3, n4, n5
      logical :: aligned
      _REAL_ :: dy
      INTEGER :: nold1, nold2,nold3, nold4,nold5,order =3, err
      integer :: i1,i2,i3, i4,i5,index(3), id
      _REAL_ :: x(0:1,0:1,0:1),r,s,t
      type(C_PTR) :: cptr=C_NULL_PTR

#ifdef RISM_DEBUG
      write(6,*) "safemem_realloc_5d_real_interp", associated(p),n1,n2,n3,n4,n5
      call flush(6)
#endif /*RISM_DEBUG*/    
      aligned = .false.
      if(present(o_aligned)) aligned = o_aligned

      nold1 = MIN(Ubound(p,1), n1)
      nold2 = MIN(UBOUND(p,2), n2)
      nold3 = MIN(UBOUND(p,3), n3)
      nold4 = MIN(UBOUND(p,4), n4)
      nold5 = MIN(UBOUND(p,5), n5)
      if(aligned)then
         err = 0
         cptr = fftw_alloc_real(int(n1*n2*n3*n4*n5, C_SIZE_T))
         if(.not.c_associated(cptr)) then
            err = 1
         else
            call c_f_pointer(cptr,safemem_realloc_5d_real_interp,[n1,n2,n3,n4,n5])
         end if
      else
         ALLOCATE(safemem_realloc_5d_real_interp(n1,n2,n3,n4,n5), STAT=err)
      end if
      IF(err /= 0)&
         call rism_report_error("reallocation_pointer allocate error: 5D REAL Interpolate")
      call memadd_r(totalmem,product(ubound(safemem_realloc_5d_real_interp))*kind(p))
      IF(.NOT. ASSOCIATED(p)) RETURN
      nold1 = MIN(Ubound(p,1), n1)
      nold2 = MIN(UBOUND(p,2), n2)
      nold3 = MIN(UBOUND(p,3), n3)
      nold4 = MIN(UBOUND(p,4), n4)
      nold5 = MIN(UBOUND(p,5), n5)

      do i5 = 1 , n5
         do i4 = 1 , n4
            do i3 = 1 , n3
               do i2 = 1 , n2
                  do i1 = 1 , n1
                     !find nearest original point
                     index = nearestindex(xi1,xi2,xi3,ubound(p,1),ubound(p,2),ubound(p,3),xo1(i1),xo2(i2),xo3(i3))
                     if(index(1) == 1 .or. index(2) == 1 .or. index(3) == 1 .or.&
                          index(1) == ubound(p,1) .or. index(2) == ubound(p,2) .or. index(3) == ubound(p,3)) then
                        safemem_realloc_5d_real_interp(i1,i2,i3,i4,i5) = p(index(1),index(2),index(3),i4,i5)
                     else
                        safemem_realloc_5d_real_interp(i1,i2,i3,i4,i5) = p(index(1),index(2),index(3),i4,i5)
                     end if
                  end do
               end do
            end do
         end do
      end do
      if(safemem_dealloc(p,aligned) /= 0)&
           call rism_report_error("reallocation_pointer deallocate error: 5D REAL INTERPOLATE")
    END FUNCTION SAFEMEM_REALLOC_5D_REAL_INTERP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!reallocate memory of a specific pointer during the simulation. This can 
!!!handle 5 dimension and can preserve the contents of the array.
!!!IN:
!!!   p     :: original well defined pointer; the pointer association status
!!!            must be associated or disassociated, ie, nullified, but cannot
!!!            be undefined.
!!!   n1    :: size of the first dimension
!!!   n2    :: size of the second dimension
!!!   n3    :: size of the third dimension
!!!   n4    :: size of the fourth dimension
!!!   n5    :: size of the fifth dimension
!!!   o_preserve :: boolean indicating whether or not to preserve the contents of 
!!!           the array (optional)
!!!   o_center :: when preserving the contents of the array, center the 
!!!             data in the array instead of preserving the value at each
!!!             index 
!!!   o_aligned :: (optional) Use the FFTW library to allocate 16-byte
!!!                aligned memory necessary to SIMD operations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION safemem_realloc_5d_real(p, n1, n2, n3, n4, n5, o_preserve, o_center,o_aligned)
      implicit none
      _REAL_, POINTER, DIMENSION(:,:,:,:,:) :: safemem_realloc_5d_real
      _REAL_, POINTER, DIMENSION(:,:,:,:,:):: p
      logical, optional, intent(in) :: o_preserve, o_center,o_aligned
      INTEGER, intent(in) :: n1, n2, n3, n4, n5
      logical :: preserve,center,aligned
      INTEGER :: nold1, nold2, nold3, nold4, nold5, err, offset(5)
      integer :: i1,i2,i3,i4,i5
      type(C_PTR) :: cptr=C_NULL_PTR

#ifdef RISM_DEBUG
      write(6,*) "safemem_realloc_5d_real", associated(p),n1,n2,n3,n4,n5,preserve
      call flush(6)
#endif /*RISM_DEBUG*/
      aligned = .false.
      if(present(o_aligned)) aligned = o_aligned

      preserve = .true.
      if(present(o_preserve)) preserve = o_preserve

      center=.false.
      if(present(o_center)) center = o_center

      if(.not.preserve)then
         if(safemem_dealloc(p,aligned) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 5D REAL")
      end if
      if(aligned)then
         err = 0
         cptr = fftw_alloc_real(int(n1*n2*n3*n4*n5, C_SIZE_T))
         if(.not.c_associated(cptr)) then
            err = 1
         else
            call c_f_pointer(cptr,safemem_realloc_5d_real,[n1,n2,n3,n4,n5])
         end if
      else
         ALLOCATE(safemem_realloc_5d_real(n1,n2,n3,n4,n5), STAT=err)
      end if
      if(err /= 0)&
         call rism_report_error("reallocation_pointer allocate error: 5D REAL")
      call memadd_r(totalmem,product(ubound(safemem_realloc_5d_real))*kind(p))
      if(ASSOCIATED(p) .and. preserve .and. n1*n2*n3*n4*n5 > 0)then
         safemem_realloc_5d_real = 0d0
         nold1 = MIN(UBOUND(p,1), n1)
         nold2 = MIN(UBOUND(p,2), n2)
         nold3 = MIN(UBOUND(p,3), n3)
         nold4 = MIN(UBOUND(p,4), n4)
         nold5 = MIN(UBOUND(p,5), n5)
         if(center)then
            offset(1) = (UBOUND(p,1)- n1)/2
            offset(2) = (UBOUND(p,2)- n2)/2
            offset(3) = (UBOUND(p,3)- n3)/2
            offset(4) = (UBOUND(p,4)- n4)/2
            offset(5) = (UBOUND(p,5)- n5)/2
         else
            offset = 0
         end if
         do i5=1,nold5
            do i4=1,nold4
               do i3=1,nold3
                  do i2=1,nold2
                     call dcopy(nold1, &
                          p(max(offset(1),0)+1, &
                          max(offset(2),0)+i2, &
                          max(offset(3),0)+i3, &
                          max(offset(4),0)+i4, &
                          max(offset(5),0)+i5),1,&
                          safemem_realloc_5d_real(&
                          abs(min(offset(1),0))+1,&
                          abs(min(offset(2),0))+i2,&
                          abs(min(offset(3),0))+i3,&
                          abs(min(offset(4),0))+i4,&
                          abs(min(offset(5),0))+i5),1)
                  end do
               end do
            end do
         end do
         if(safemem_dealloc(p,aligned) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 5D REAL")
      endif
    END FUNCTION SAFEMEM_REALLOC_5D_REAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!reallocate memory of a specific pointer during the simulation. This can 
!!!handle 1 dimension and preserves the contents of the array.
!!!IN:
!!!   p     :: original well defined pointer; the pointer association status
!!!            must be associated or disassociated, ie, nullified, but cannot
!!!            be undefined.
!!!   n1    :: size of the first dimension
!!!   o_preserve :: boolean indicating whether or not to preserve the contents of 
!!!           the array (optional)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION safemem_realloc_1d_integer(p, n1,o_preserve)
      implicit none
      integer, POINTER, DIMENSION(:) :: safemem_realloc_1d_integer
      integer, POINTER, DIMENSION(:) :: p
      logical, optional, intent(in) :: o_preserve
      INTEGER, intent(in) :: n1
      logical :: preserve
      INTEGER :: nold1, err
      integer :: i1

#ifdef RISM_DEBUG
      write(6,*) "safemem_realloc_5d_integer", associated(p),n1,preserve
      call flush(6)
#endif /*RISM_DEBUG*/
      preserve = .true.
      if(present(o_preserve)) preserve = o_preserve

      if(.not.preserve)then
         if(safemem_dealloc(p) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 1D INTEGER")
      end if
      ALLOCATE(safemem_realloc_1d_integer(1:n1), STAT=err)
      call memadd_i(totalmem,product(ubound(safemem_realloc_1d_integer))*kind(p))
      IF(err /= 0)&
         call rism_report_error("reallocation_pointer allocate error: 1D INTEGER")
      if(associated(p) .and. preserve .and. n1 > 0) then
         nold1 = MIN(UBOUND(p,1), n1)
         !explicit loop to prevent stack overflow with intel compilers
         do i1=1,nold1
            safemem_realloc_1d_integer(i1) = p(i1)
         end do
         if(safemem_dealloc(p) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 1D INTEGER")
      end if
    END FUNCTION SAFEMEM_REALLOC_1D_INTEGER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!reallocate memory of a specific pointer during the simulation. This can 
!!!handle 2 dimension and preserves the contents of the array.
!!!IN:
!!!   p     :: original well defined pointer; the pointer association status
!!!            must be associated or disassociated, ie, nullified, but cannot
!!!            be undefined.
!!!   n1    :: size of the first dimension
!!!   n2    :: size of the second dimension
!!!   o_preserve :: boolean indicating whether or not to preserve the contents of 
!!!           the array (optional)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION safemem_realloc_2d_integer(p, n1, n2,o_preserve)
      implicit none
      integer, POINTER, DIMENSION(:,:) :: safemem_realloc_2d_integer
      integer, POINTER, DIMENSION(:,:) :: p
      logical, optional, intent(in) :: o_preserve
      INTEGER, intent(in) :: n1,n2
      logical :: preserve
      INTEGER :: nold1,nold2, err
      integer :: i1, i2

#ifdef RISM_DEBUG
      write(6,*) "safemem_realloc_2d_integer", associated(p),n1,n2,preserve
      call flush(6)
#endif /*RISM_DEBUG*/
      preserve=.true.
      if(present(o_preserve)) preserve = o_preserve

      if(.not.preserve)then
         if(safemem_dealloc(p) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 2D INTEGER")
      end if
      ALLOCATE(safemem_realloc_2d_integer(1:n1,1:n2), STAT=err)
      call memadd_i(totalmem,product(ubound(safemem_realloc_2d_integer))*kind(p))
      IF(err /= 0)&
         call rism_report_error("reallocation_pointer allocate error: 2D INTEGER")
      if(associated(p) .and. preserve .and. n1*n2 > 0) then
         nold1 = MIN(ubound(p,1), n1)
         nold2 = MIN(ubound(p,2), n2)
         !explicit loop to prevent stack overflow with intel compilers
         do i1=1,nold1
            do i2=1,nold2
               safemem_realloc_2d_integer(i1,i2) = p(i1,i2)
            end do
         end do
         if(safemem_dealloc(p) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 2D INTEGER")
      end if
    END FUNCTION SAFEMEM_REALLOC_2D_INTEGER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!reallocate memory of a specific pointer during the simulation. This can 
!!!handle 1 dimension and preserves the contents of the array.
!!!IN:
!!!   p     :: original well defined pointer; the pointer association status
!!!            must be associated or disassociated, ie, nullified, but cannot
!!!            be undefined.
!!!   length :: character length for each element
!!!   n1    :: size of the first dimension
!!!   o_preserve :: boolean indicating whether or not to preserve the contents of 
!!!           the array (optional)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION safemem_realloc_1d_character(p, length, n1,o_preserve)
      implicit none
      integer, intent(in) :: length
      character(len=length), POINTER, DIMENSION(:) :: safemem_realloc_1d_character
      character(len=length), POINTER, DIMENSION(:) :: p
      logical, optional, intent(in) :: o_preserve
      INTEGER, intent(in) :: n1
      logical :: preserve
      INTEGER :: nold1, err
      integer :: i1

#ifdef RISM_DEBUG
      write(6,*) "safemem_realloc_1d_character", associated(p),n1,o_preserve
      call flush(6)
#endif /*RISM_DEBUG*/    
      preserve = .true.
      if(present(o_preserve)) preserve = o_preserve

      if(.not.preserve)then
         if(safemem_dealloc(p) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 1D CHARACTER")
      end if
      ALLOCATE(safemem_realloc_1d_character(1:n1), STAT=err)
      call memadd_c(totalmem,product(ubound(safemem_realloc_1d_character))*kind(p))
      IF(err /= 0)&
         call rism_report_error("reallocation_pointer allocate error: 1D CHARACTER")
      if(associated(p) .and. preserve .and. n1 > 0) then
         nold1 = MIN(UBOUND(p,1), n1)
         !explicit loop to prevent stack overflow with intel compilers
         do i1=1,nold1
            safemem_realloc_1d_character(i1) = trim(p(i1))
         end do
         if(safemem_dealloc(p) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 1D CHARACTER")
      end if
    END FUNCTION SAFEMEM_REALLOC_1D_CHARACTER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!reallocate memory of a specific pointer during the simulation. This can 
!!!handle 2 dimensions and preserves the contents of the array.
!!!IN:
!!!   p     :: original well defined pointer; the pointer association status
!!!            must be associated or disassociated, ie, nullified, but cannot
!!!            be undefined.
!!!   length :: character length for each element
!!!   n1    :: size of the first dimension
!!!   n2    :: size of the second dimension
!!!   o_preserve :: boolean indicating whether or not to o_preserve the contents of 
!!!           the array (optional)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION safemem_realloc_2d_character(p, length, n1, n2,o_preserve)
      implicit none
      integer, intent(in) :: length
      character(len=length), POINTER, DIMENSION(:,:) :: safemem_realloc_2d_character
      character(len=length), POINTER, DIMENSION(:,:) :: p
      logical, optional, intent(in) :: o_preserve
      INTEGER, intent(in) :: n1, n2
      logical :: preserve
      INTEGER :: nold1, nold2, err
      integer :: i1, i2

#ifdef RISM_DEBUG
      write(6,*) "safemem_realloc_2d_character", associated(p),n1,n2,o_preserve
      call flush(6)
#endif /*RISM_DEBUG*/    
      preserve = .true.
      if(present(o_preserve)) preserve = o_preserve

      if(.not.preserve)then
         if(safemem_dealloc(p) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 2D CHARACTER")
      end if
      ALLOCATE(safemem_realloc_2d_character(1:n1,1:n2), STAT=err)
      call memadd_c(totalmem,product(ubound(safemem_realloc_2d_character))*kind(p))
      IF(err /= 0)&
         call rism_report_error("reallocation_pointer allocate error: 1D CHARACTER")
      if(associated(p) .and. preserve .and. n1*n2 > 0) then
         nold1 = MIN(UBOUND(p,1), n1)
         nold2 = MIN(UBOUND(p,2), n2)
         !explicit loop to prevent stack overflow with intel compilers
         do i1=1,nold1
            do i2=1,nold2
               safemem_realloc_2d_character(i1,i2) = trim(p(i1,i2))
            end do
         end do
         if(safemem_dealloc(p) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 2D CHARACTER")
      end if
    END FUNCTION SAFEMEM_REALLOC_2D_CHARACTER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!reallocate memory of a specific pointer during the simulation. This can 
!!!handle 1 dimension and preserves the contents of the array.
!!!IN:
!!!   p     :: original well defined pointer; the pointer association status
!!!            must be associated or disassociated, ie, nullified, but cannot
!!!            be undefined.
!!!   n1    :: size of the first dimension
!!!   o_preserve :: boolean indicating whether or not to preserve the contents of 
!!!           the array (optional)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    FUNCTION safemem_realloc_1d_logical(p, n1,o_preserve)
      implicit none
      logical, POINTER, DIMENSION(:) :: safemem_realloc_1d_logical
      logical, POINTER, DIMENSION(:) :: p
      logical, optional, intent(in) :: o_preserve
      INTEGER, intent(in) :: n1
      logical :: preserve
      INTEGER :: nold1, err
      integer :: i1

#ifdef RISM_DEBUG
      write(6,*) "safemem_realloc_1d_logical", associated(p),n1,preserve
      call flush(6)
#endif /*RISM_DEBUG*/
      preserve = .true.
      if(present(o_preserve)) preserve = o_preserve

      if(.not.preserve)then
         if(safemem_dealloc(p) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 1D LOGICAL")
      end if
      ALLOCATE(safemem_realloc_1d_logical(1:n1), STAT=err)
      call memadd_l(totalmem,product(ubound(safemem_realloc_1d_logical))*kind(p))
      IF(err /= 0)&
         call rism_report_error("reallocation_pointer allocate error: 1D LOGICAL")
      if(associated(p) .and. preserve .and. n1 > 0) then
         nold1 = MIN(UBOUND(p,1), n1)
         !explicit loop to prevent stack overflow with intel compilers
         do i1=1,nold1
            safemem_realloc_1d_logical(i1) = p(i1)
         end do
         if(safemem_dealloc(p) /= 0)&
              call rism_report_error("reallocation_pointer deallocate error: 1D LOGICAL")
      end if
    END FUNCTION SAFEMEM_REALLOC_1D_LOGICAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!checks if the pointer is allocated before deallocating.  Returns any error code  
!!!IN:
!!!   p :: well defined pointer to be deallocated; the pointer association status
!!!        must be associated or disassociated, ie, nullified, but cannot
!!!        be undefined.
!!!   o_aligned :: (optional) Use the FFTW library to free 16-byte
!!!                aligned memory.  It is not possible to tell how
!!!                memory was originally allocated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function safemem_dealloc_pointer_1d_real(p,o_aligned)
      implicit none
      _REAL_, pointer, dimension(:) :: p
      logical ,optional, intent(in) :: o_aligned
      logical :: aligned
      type(c_ptr) :: cptr
      integer :: safemem_dealloc_pointer_1d_real,temp=0
     
      aligned = .false.
      if(present(o_aligned)) aligned = o_aligned

      !work around for old gfortran bug
      !  if(associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_1d_real)
      if(associated(p))then
         call memadd_r(totalmem,-product(ubound(p))*kind(p))
         if(aligned)then
            cptr = c_loc(p(1))
            call fftw_free(cptr)
            nullify(p)
            temp=0
         else
            deallocate(p,STAT=temp)
         end if
      end if
      safemem_dealloc_pointer_1d_real=temp
    end function safemem_dealloc_pointer_1d_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!checks if the pointer is allocated before deallocating.  Returns any error code  
!!!IN:
!!!   p :: well defined pointer to be deallocated; the pointer association status
!!!        must be associated or disassociated, ie, nullified, but cannot
!!!        be undefined.
!!!   o_aligned :: (optional) Use the FFTW library to free 16-byte
!!!                aligned memory.  It is not possible to tell how
!!!                memory was originally allocated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function safemem_dealloc_pointer_2d_real(p,o_aligned)
      implicit none
      _REAL_, pointer, dimension(:,:) :: p
      logical ,optional, intent(in) :: o_aligned
      logical :: aligned
      type(c_ptr) :: cptr
      integer :: safemem_dealloc_pointer_2d_real,temp=0
     
      aligned = .false.
      if(present(o_aligned)) aligned = o_aligned

      !work around for old gfortran bug
      !  if(associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_2d_real)
      if(associated(p))then
         call memadd_r(totalmem,-product(ubound(p))*kind(p))
         if(aligned)then
            cptr = c_loc(p(1,1))
            call fftw_free(cptr)
            nullify(p)
         else
            deallocate(p,STAT=temp)
         end if
      end if
      safemem_dealloc_pointer_2d_real=temp
    end function safemem_dealloc_pointer_2d_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!checks if the pointer is allocated before deallocating.  Returns any error code  
!!!IN:
!!!   p :: well defined pointer to be deallocated; the pointer association status
!!!        must be associated or disassociated, ie, nullified, but cannot
!!!        be undefined.
!!!   o_aligned :: (optional) Use the FFTW library to free 16-byte
!!!                aligned memory.  It is not possible to tell how
!!!                memory was originally allocated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function safemem_dealloc_pointer_3d_real(p,o_aligned)
      implicit none
      _REAL_, pointer, dimension(:,:,:) :: p
      logical ,optional, intent(in) :: o_aligned
      logical :: aligned
      type(c_ptr) :: cptr
      integer :: safemem_dealloc_pointer_3d_real,temp=0
     
      aligned = .false.
      if(present(o_aligned)) aligned = o_aligned

      !work around for old gfortran bug
      !  if(associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_3d_real)
      if(associated(p))then
         call memadd_r(totalmem,-product(ubound(p))*kind(p))
         if(aligned)then
            cptr = c_loc(p(1,1,1))
            call fftw_free(cptr)
            nullify(p)
         else
            deallocate(p,STAT=temp)
         end if
      end if
      safemem_dealloc_pointer_3d_real=temp
    end function safemem_dealloc_pointer_3d_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!checks if the pointer is allocated before deallocating.  Returns any error code  
!!!IN:
!!!   p :: well defined pointer to be deallocated; the pointer association status
!!!        must be associated or disassociated, ie, nullified, but cannot
!!!        be undefined.
!!!   o_aligned :: (optional) Use the FFTW library to free 16-byte
!!!                aligned memory.  It is not possible to tell how
!!!                memory was originally allocated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function safemem_dealloc_pointer_4d_real(p,o_aligned)
      implicit none
      _REAL_, pointer, dimension(:,:,:,:) :: p
      logical ,optional, intent(in) :: o_aligned
      logical :: aligned
      type(c_ptr) :: cptr
      integer :: safemem_dealloc_pointer_4d_real,temp=0
     
      aligned = .false.
      if(present(o_aligned)) aligned = o_aligned

      !work around for old gfortran bug
      !  if(associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_4d_real)
      if(associated(p))then
         call memadd_r(totalmem,-product(ubound(p))*kind(p))
         if(aligned)then
            cptr = c_loc(p(1,1,1,1))
            call fftw_free(cptr)
            nullify(p)
         else
            deallocate(p,STAT=temp)
         end if
      end if
      safemem_dealloc_pointer_4d_real=temp
    end function safemem_dealloc_pointer_4d_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!checks if the pointer is allocated before deallocating.  Returns any error code  
!!!IN:
!!!   p :: well defined pointer to be deallocated; the pointer association status
!!!        must be associated or disassociated, ie, nullified, but cannot
!!!        be undefined.
!!!   o_aligned :: (optional) Use the FFTW library to free 16-byte
!!!                aligned memory.  It is not possible to tell how
!!!                memory was originally allocated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function safemem_dealloc_pointer_5d_real(p,o_aligned)
      implicit none
      _REAL_, pointer, dimension(:,:,:,:,:) :: p
      logical ,optional, intent(in) :: o_aligned
      logical :: aligned
      type(c_ptr) :: cptr
      integer :: safemem_dealloc_pointer_5d_real,temp=0
     
      aligned = .false.
      if(present(o_aligned)) aligned = o_aligned

       !work around for old gfortran bug
      !  if(associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_5d_real)
      if(associated(p))then
         call memadd_r(totalmem,-product(ubound(p))*kind(p))
         if(aligned)then
            cptr = c_loc(p(1,1,1,1,1))
            call fftw_free(cptr)
            nullify(p)
         else
            deallocate(p,STAT=temp)
         end if
      end if
      safemem_dealloc_pointer_5d_real=temp
    end function safemem_dealloc_pointer_5d_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!checks if the pointer is allocated before deallocating.  Returns any error code  
!!!IN:
!!!   p :: well defined pointer to be deallocated; the pointer association status
!!!        must be associated or disassociated, ie, nullified, but cannot
!!!        be undefined.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function safemem_dealloc_pointer_1d_int(p)
      implicit none
      integer, pointer, dimension(:) :: p
      integer :: safemem_dealloc_pointer_1d_int,temp=0
      !work around for old gfortran bug
      !  if(associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_1d_int)
      if(associated(p))then
         call memadd_i(totalmem,-product(ubound(p))*kind(p))
         deallocate(p,STAT=temp)
      end if
      safemem_dealloc_pointer_1d_int=temp
    end function safemem_dealloc_pointer_1d_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!checks if the pointer is allocated before deallocating.  Returns any error code  
!!!IN:
!!!   p :: well defined pointer to be deallocated; the pointer association status
!!!        must be associated or disassociated, ie, nullified, but cannot
!!!        be undefined.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function safemem_dealloc_pointer_2d_int(p)
      implicit none
      integer, pointer, dimension(:,:) :: p
      integer :: safemem_dealloc_pointer_2d_int,temp=0
      !work around for old gfortran bug
      !  if(associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_2d_int)
      if(associated(p))then
         call memadd_i(totalmem,-product(ubound(p))*kind(p))
         deallocate(p,STAT=temp)
      end if
      safemem_dealloc_pointer_2d_int=temp
    end function safemem_dealloc_pointer_2d_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!checks if the pointer is allocated before deallocating.  Returns any error code  
!!!IN:
!!!   p :: well defined pointer to be deallocated; the pointer association status
!!!        must be associated or disassociated, ie, nullified, but cannot
!!!        be undefined.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function safemem_dealloc_pointer_3d_int(p)
      implicit none
      integer, pointer, dimension(:,:,:) :: p
      integer :: safemem_dealloc_pointer_3d_int,temp=0
      !work around for old gfortran bug
      !  if(associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_3d_int)
      if(associated(p))then
         call memadd_i(totalmem,-product(ubound(p))*kind(p))
         deallocate(p,STAT=temp)
      end if
      safemem_dealloc_pointer_3d_int=temp
    end function safemem_dealloc_pointer_3d_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!checks if the pointer is allocated before deallocating.  Returns any error code  
!!!IN:
!!!   p :: well defined pointer to be deallocated; the pointer association status
!!!        must be associated or disassociated, ie, nullified, but cannot
!!!        be undefined.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function safemem_dealloc_pointer_4d_int(p)
      implicit none
      integer, pointer, dimension(:,:,:,:) :: p
      integer :: safemem_dealloc_pointer_4d_int,temp=0
      !work around for old gfortran bug
      !  if(associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_4d_int)
      if(associated(p))then
         call memadd_i(totalmem,-product(ubound(p))*kind(p))
         deallocate(p,STAT=temp)
      end if
      safemem_dealloc_pointer_4d_int=temp
    end function safemem_dealloc_pointer_4d_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!checks if the pointer is allocated before deallocating.  Returns any error code  
!!!IN:
!!!   p :: well defined pointer to be deallocated; the pointer association status
!!!        must be associated or disassociated, ie, nullified, but cannot
!!!        be undefined.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function safemem_dealloc_pointer_5d_int(p)
      implicit none
      integer, pointer, dimension(:,:,:,:,:) :: p
      integer :: safemem_dealloc_pointer_5d_int,temp=0
      !work around for old gfortran bug
      !  if(associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_5d_int)
      if(associated(p))then
         call memadd_i(totalmem,-product(ubound(p))*kind(p))
         deallocate(p,STAT=temp)
      end if
      safemem_dealloc_pointer_5d_int=temp
    end function safemem_dealloc_pointer_5d_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!checks if the pointer is allocated before deallocating.  Returns any error code  
!!!IN:
!!!   p :: well defined pointer to be deallocated; the pointer association status
!!!        must be associated or disassociated, ie, nullified, but cannot
!!!        be undefined.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function safemem_dealloc_pointer_1d_character(p)
      implicit none
      character(len=*), pointer, dimension(:) :: p
      integer :: safemem_dealloc_pointer_1d_character,temp=0
      !work around for old gfortran bug
      !  if(associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_1d_character)
      if(associated(p))then
         call memadd_c(totalmem,-product(ubound(p))*kind(p))
         deallocate(p,STAT=temp)
      end if
      safemem_dealloc_pointer_1d_character = temp
    end function safemem_dealloc_pointer_1d_character

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!checks if the pointer is allocated before deallocating.  Returns any error code  
!!!IN:
!!!   p :: well defined pointer to be deallocated; the pointer association status
!!!        must be associated or disassociated, ie, nullified, but cannot
!!!        be undefined.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function safemem_dealloc_pointer_2d_character(p)
      implicit none
      character(len=*), pointer, dimension(:,:) :: p
      integer :: safemem_dealloc_pointer_2d_character,temp=0
      !work around for old gfortran bug
      !  if(associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_2d_character)
      if(associated(p))then
         call memadd_c(totalmem,-product(ubound(p))*kind(p))
         deallocate(p,STAT=temp)
      end if
      safemem_dealloc_pointer_2d_character = temp
    end function safemem_dealloc_pointer_2d_character

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!checks if the pointer is allocated before deallocating.  Returns any error code  
!!!IN:
!!!   p :: well defined pointer to be deallocated; the pointer association status
!!!        must be associated or disassociated, ie, nullified, but cannot
!!!        be undefined.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function safemem_dealloc_pointer_1d_logical(p)
      implicit none
      logical, pointer, dimension(:) :: p
      integer :: safemem_dealloc_pointer_1d_logical,temp=0
      !work around for old gfortran bug
      !  if(associated(p)) deallocate(p,STAT=safemem_dealloc_pointer_1d_logical)
      if(associated(p))then
         call memadd_l(totalmem,-product(ubound(p))*kind(p))
         deallocate(p,STAT=temp)
      end if
      safemem_dealloc_pointer_1d_logical = temp
    end function safemem_dealloc_pointer_1d_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns an array of current and maximum memory allocated in bytes so
!!!far for a specific memTracker.  No distinction is made between different 
!!!kinds of the same type:
!!!(1) : Current integer memory
!!!(2) : Current real memory
!!!(3) : Current logical memory
!!!(4) : Current char memory
!!!(5) : Current logical memory
!!!(6) : Max integer memory
!!!(7) : Max real memory
!!!(8) : Max logical memory
!!!(9) : Max char memory
!!!(10): Max logical memory
!!!IN:
!!!   this : memTracker object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function memStatus_obj(this) result(mem)
      implicit none
      type(memTracker), intent(in) :: this
      integer(8) :: mem(10)
      mem = (/this%int, this%real, this%logical, this%char, this%total,&
           this%maxint, this%maxreal, this%maxlogical, this%maxchar, this%max/)
    end function memStatus_obj

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns an array of global current and maximum memory allocated in bytes so
!!!far.  No distinction is made between different kinds of the same type:
!!!(1) : Current integer memory
!!!(2) : Current real memory
!!!(3) : Current logical memory
!!!(4) : Current char memory
!!!(5) : Current logical memory
!!!(6) : Max integer memory
!!!(7) : Max real memory
!!!(8) : Max logical memory
!!!(9) : Max char memory
!!!(10): Max logical memory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function memStatus_global() result(mem)
      implicit none
      integer(8) :: mem(10)
      mem = memStatus_obj(totalMem)
    end function memStatus_global

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Adds the size in bytes of the array and increments the integer counter.  Use
!!!a negative number to subtract.
!!!IN:
!!!   this  : memtracker object
!!!   nbyte : number of bytes 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine memadd_i(this,nbyte)
      implicit none
      type(memTracker),intent(inout) :: this
      integer, intent(in) :: nbyte
      this%int = this%int+nbyte
      this%total = this%total+nbyte
      if(this%int > this%maxint) this%maxint = this%int
      if(this%total > this%max) this%max = this%total
    end subroutine memadd_i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Adds the size in bytes of the array and increments the real counter.  Use
!!!a negative number to subtract.
!!!IN:
!!!   this  : memtracker object
!!!   nbyte : number of bytes 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine memadd_r(this,nbyte)
      implicit none
      type(memTracker),intent(inout) :: this
      integer, intent(in) :: nbyte
      this%real = this%real+nbyte
      this%total = this%total+nbyte
      if(this%real > this%maxreal) this%maxreal = this%real
      if(this%total > this%max) this%max = this%total
    end subroutine memadd_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Adds the size in bytes of the array and increments the character counter.  Use
!!!a negative number to subtract.
!!!IN:
!!!   this  : memtracker object
!!!   nbyte : number of bytes 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine memadd_c(this,nbyte)
      implicit none
      type(memTracker),intent(inout) :: this
      integer, intent(in) :: nbyte
      this%char = this%char+nbyte
      this%total = this%total+nbyte
      if(this%char > this%maxchar) this%maxchar = this%char
      if(this%total > this%max) this%max = this%total
    end subroutine memadd_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Adds the size in bytes of the array and increments the logical counter.  Use
!!!a negative number to subtract.
!!!IN:
!!!   this  : memtracker object
!!!   nbyte : number of bytes 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine memadd_l(this,nbyte)
      implicit none
      type(memTracker),intent(inout) :: this
      integer, intent(in) :: nbyte
      this%logical = this%logical+nbyte
      this%total = this%total+nbyte
      if(this%logical > this%maxlogical) this%maxlogical = this%logical
      if(this%total > this%max) this%max = this%total
    end subroutine memadd_l
  end module safemem
