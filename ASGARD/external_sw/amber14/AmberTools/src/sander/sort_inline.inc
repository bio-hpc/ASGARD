! vi:syntax=fortran
! Fast in-line hybrid QSORT+INSERTION SORT for Fortran. -- Joe Krahn, 2006
!
!------------------------------------------------------------------------
! PURPOSE:
! Sort a one-dimensional array according to a user-defined comparison.
! It is approximately the same as C qsort(), but without function
! calls for every comparison step. Some compilers may already inline
! their own QSORT functions.
!
!------------------------------------------------------------------------
! USAGE:
!
! You must define the following five entities. The first 3 are procedures
! that must be defined as internal routines, because they access ARRAY
! by host association. It is also possible to define them as CPP macros.
! All of them are passed two array indices as arguyments.
!
! 1) LOGICAL FUNCTION SORT_COMPARE(a,b)
!        INTEGER, INTENT(IN) :: a,b
!
!   Return the logical result of (array(a)%key < array(b)%key),
!   assuming an ascending order sort.
!
! 2) SUBROUTINE SORT_SWAP(a,b)
!   Swap array members A and B. (QSORT part)
!
! 3) SUBROUTINE SORT_SHIFT(a,b)
!   Shift array members (A:B-1) to (A+1:B), wrapping member B
!   to the A position. (Insertion sort part)
!
! 4) INTEGER, PARAMETER :: QSORT_THRESHOLD
!
!   The QSORT is used down to the QSORT_THRESHOLD size sorted blocks.
!   Then insertion sort is used for the remainder, because it is faster
!   for small sort ranges.
!
!   Large sorting elements or small cache may make a smaller threshold
!   more useful. You can also set this to a run-time argument/variable
!   with no performance loss.
!
!   A fairly small value, like 8 or 16, seems to be a safe default,
!   getting most of the speed-up, without risking any cache issues on
!   larger data.
!
! 5) INTEGER, PARAMETER :: SORT_DIMENSION
!
!   This defines the sort axis used to determine the array size, lbound
!   and ubound. For one-dimentional arrays, it is just 1.
!
!------------------------------------------------------------------------
! NOTES:
! The procedure uses a optimized combination of QSORT and INSERTION
! sorting. The algorithm is based on code used in GLIBC. 
! A stack is used in place of recursive calls. The stack size must
! be at least as big as the number of bits in the largest array index.
!
! Sorting vectors of a multidimensional allocatable array can be
! VERY slow. In this case, or with large derived types, it is better
! to sort a simple derived type of key/index pairs, then reorder
! tha actual data using the sorted indices.
!
!------------------------------------------------------------------------
! EXAMPLES
!
! Sort an integer array:
!
! subroutine SORT_DATA(array)
!   implicit none
!   integer :: hold, array(:)
!   integer, parameter :: QSORT_THRESHOLD=96, SORT_DIMENSION=1
!   include 'sort_inline.inc'
! contains
!   logical function SORT_COMPARE(a,b)
!     integer, intent(in) :: a,b
!     SORT_COMPARE = ( array(a) < array(b) )
!   end function SORT_COMPARE
!   subroutine SORT_SWAP(a,b)
!      integer, intent(in) :: a,b
!      integer :: hold
!      hold=array(a)
!      array(a)=array(b)
!      array(b)=hold
!   end subroutine SORT_SWAP
!   subroutine SORT_SHIFT(a,b)
!      integer, intent(in) :: a,b
!      integer :: hold
!      hold=array(b)
!      array(a+1:b)=array(a:b-1)
!      array(a)=hold
!   end subroutine SORT_SHIFT
! end subroutine SORT_DATA
!
!
! Sort an index to an integer array.
! First set up the index:  index(:) = (/i,i=1,size(index)/)
!
! subroutine SORT_DATA(array,index)
!   implicit none
!   integer :: hold, array(:), index(:)
!   integer, parameter :: QSORT_THRESHOLD=96, SORT_DIMENSION=1
!   include 'sort_inline.inc'
! contains
!   logical function SORT_COMPARE(a,b)
!     integer, intent(in) :: a,b
!     SORT_COMPARE = ( array(index(a)) < array(index(b)) )
!   end function SORT_COMPARE
!   subroutine SORT_SWAP(a,b)
!      integer, intent(in) :: a,b
!      integer :: hold
!      hold=index(a)
!      index(a)=index(b)
!      index(b)=hold
!   end subroutine SORT_SWAP
!   subroutine SORT_SHIFT(a,b)
!      integer, intent(in) :: a,b
!      integer :: hold
!      hold=index(b)
!      index(a+1:b)=index(a:b-1)
!      index(a)=hold
!   end subroutine SORT_SHIFT
! end subroutine SORT_DATA
!------------------------------------------------------------------------

  integer :: stack_top, right_size,left_size
  integer :: mid, left, right, low, high
  integer :: ARRAY_SIZE, ARRAY_LBOUND, ARRAY_UBOUND

! A stack of 32 can handle the entire extent of a 32-bit index.
! Use 64 for 64-bit indexed arrays, if they might contain
! more than 2^32 elements.
  integer, parameter :: SORT_STACK_SIZE = 32
  type sort_stack; integer :: low, high; end type
  type(sort_stack) :: stack(SORT_STACK_SIZE)

  array_size=size(array,SORT_DIMENSION)
  array_lbound=(lbound(array,SORT_DIMENSION))
  array_ubound=(ubound(array,SORT_DIMENSION))

  if (ARRAY_SIZE > QSORT_THRESHOLD) then
    low = ARRAY_LBOUND
    high = ARRAY_UBOUND
    stack_top = 0
    QSORT_LOOP: do
      mid = (low + high)/2
      if (SORT_COMPARE (mid, low)) then
        call sort_swap(mid,low)
      end if
      if (SORT_COMPARE (high, mid)) then
        call sort_swap(high,mid)
        if (SORT_COMPARE (mid, low)) then
          call sort_swap(mid,low)
        end if
      end if
      left  = low + 1
      right = high - 1

      COLLAPSE_WALLS: do
        do while (SORT_COMPARE (left, mid))
          left=left+1
        end do
        do while (SORT_COMPARE (mid, right))
          right=right-1
        end do
        if (left < right) then
          call sort_swap(left,right)
          if (mid == left) then
            mid = right
          else if (mid == right) then
            mid = left
          end if
          left=left+1
          right=right-1
        else
          if (left == right) then
            left=left+1
            right=right-1
          end if
          exit COLLAPSE_WALLS
        end if
      end do COLLAPSE_WALLS

! Set up indices for the next iteration.
! Determine left and right partition sizes.
! Defer partitions smaller than the QSORT_THRESHOLD.
! If both partitions are significant,
! push the larger one onto the stack.
      right_size = right - low
      left_size = high - left
      if (right_size <= QSORT_THRESHOLD) then
        if (left_size <= QSORT_THRESHOLD) then
          ! Ignore both small partitions: Pop a partition or exit.
          if (stack_top<1) exit QSORT_LOOP
          low=stack(stack_top)%low; high=stack(stack_top)%high
          stack_top=stack_top-1
        else
          ! Ignore small left partition.
          low = left
        end if
      else if (left_size <= QSORT_THRESHOLD) then
        ! Ignore small right partition.
        high = right
      else if (right_size > left_size) then
        ! Push larger left partition indices.
        stack_top=stack_top+1
        stack(stack_top)=sort_stack(low,right)
        low = left
      else
        ! Push larger right partition indices.
        stack_top=stack_top+1
        stack(stack_top)=sort_stack(left,high)
        high = right
      end if
    end do QSORT_LOOP
  end if

! Sort the remaining small partitions using insertion sort,
! which should be faster for partitions smaller than the
! appropriate QSORT_THRESHOLD.

! First, find smallest element in first QSORT_THRESHOLD and place
! it at the array's beginning.  This operation speeds up the
! insertion sort's inner loop.
  low = ARRAY_LBOUND
  high = ARRAY_UBOUND
  left=low
  do right = low+1, MIN(low+QSORT_THRESHOLD,high)
    if (SORT_COMPARE(right,left)) left=right
  end do
  if (left/=low) then
    call sort_swap(left,low)
  end if

! Insertion sort, from left to right.
! (assuming that the left is the lowest numbered index)
  INSERTION_SORT: do right = low+2,high
    left=right-1
    if (SORT_COMPARE(right,left)) then
      do while (SORT_COMPARE (right, left-1))
        left=left-1
      end do
      call sort_shift(left,right)
    end if
  end do INSERTION_SORT

