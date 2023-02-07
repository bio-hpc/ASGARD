! <compile=optimized>
#include "copyright.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ A collection of useful (but difficult-to-place) subroutines.
module sander_lib

   implicit none

contains

! These subroutines should have few-to-no dependencies, or we will quickly
! reach a point where we get cyclic dependencies that kills compilation

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Checks if **s exist in a given inpcrd line and emit helpful messages
subroutine check_inpcrd_overflow(line, periodic)

   implicit none

   ! Passed variables
   character(*), intent(in) :: line     ! line buffer that offends the parser
   logical, intent(in)      :: periodic ! Does our simulation have PBCs?

   ! Local variables
   integer      :: i ! counter

   
   ! My advice to you depends on whether you're running a simulation with
   ! PBCs or in implicit solvent/vacuum
   do i = 1, len_trim(line)
      if (line(i:i) .eq. '*' .and. periodic) then
         write(6, '(a)') '*s in the inpcrd file often indicate an overflow of &
            &the Fortran format used', &
            'to store coordinates in the inpcrd/restart files. &
            &This often happens when', &
            'coordinates are not wrapped into the center cell &
            &(when iwrap = 0) and some', &
            'particles diffuse too far away. Try restarting from your &
            &last good restart', &
            'file and setting iwrap=1 or using a NetCDF restart file &
            &format. See the', &
            'Amber manual for details'
         ! Only print error message once. Bail out here.
         exit
      else if (line(i:i) .eq. '*' .and. .not. periodic) then
         write(6, '(a)') '*s in the inpcrd file often indicate an overflow of &
            &the Fortran format used', &
            'to store coordinates in the inpcrd/restart files. &
            &This often happens when', &
            'particles diffuse very far away from each other. Make &
            &sure you are removing', &
            'center-of-mass translation (nscm /= 0) or check if you &
            &have multiple, mobile', &
            'molecules that have diffused very far away from each other. &
            &This condition is', &
            'highly unusual for non-periodic simulations.'
         ! Only print error message once. Bail out here.
         exit
      end if
   end do

   return

end subroutine check_inpcrd_overflow

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Makes a string completely upper-case
subroutine upper(instring)

  implicit none

! Passed variables

  character(*), intent(in out) :: instring

! Local variables
  
  integer :: i ! counter

  do i = 1, len_trim(instring)

    select case (instring(i:i))
      case('a')
        instring(i:i) = 'A'
      case('b')
        instring(i:i) = 'B'
      case('c')
        instring(i:i) = 'C'
      case('d')
        instring(i:i) = 'D'
      case('e')
        instring(i:i) = 'E'
      case('f')
        instring(i:i) = 'F'
      case('g')
        instring(i:i) = 'G'
      case('h')
        instring(i:i) = 'H'
      case('i')
        instring(i:i) = 'I'
      case('j')
        instring(i:i) = 'J'
      case('k')
        instring(i:i) = 'K'
      case('l')
        instring(i:i) = 'L'
      case('m')
        instring(i:i) = 'M'
      case('n')
        instring(i:i) = 'N'
      case('o')
        instring(i:i) = 'O'
      case('p')
        instring(i:i) = 'P'
      case('q')
        instring(i:i) = 'Q'
      case('r')
        instring(i:i) = 'R'
      case('s')
        instring(i:i) = 'S'
      case('t')
        instring(i:i) = 'T'
      case('u')
        instring(i:i) = 'U'
      case('v')
        instring(i:i) = 'V'
      case('w')
        instring(i:i) = 'W'
      case('x')
        instring(i:i) = 'X'
      case('y')
        instring(i:i) = 'Y'
      case('z')
        instring(i:i) = 'Z'
    end select
  end do

  return

end subroutine upper

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Strips leading whitespace
subroutine strip(instring)

  implicit none

! Passed variables 

  character(*), intent(in out) :: instring

! Local variables

  integer :: i, begin

  begin = 1

  do i = 1, len_trim(instring)

    if (instring(i:i) .gt. ' ') then
      begin = i
      exit
    end if

  end do

  instring = instring(begin:len_trim(instring))

  return

end subroutine strip

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Determines how many words are in a string
subroutine get_num_tokens(string, token_num)

  implicit none

! Passed arguments

  character(*), intent(in) :: string

  integer, intent(out)     :: token_num

! Local variables

  integer :: string_loc  ! our location in the string
  integer :: iend        ! last non-whitespace character location

  string_loc = 1
  iend = len_trim(string)
  token_num = 0

  do while (string_loc .le. iend)

    if ( string(string_loc:string_loc) .le. ' ' ) then
      string_loc = string_loc + 1
    else

      do while ( string(string_loc:string_loc) .gt. ' ' )
        string_loc = string_loc + 1
      end do

      token_num = token_num + 1
    end if
  end do

end subroutine get_num_tokens

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Gets a specific word from a string
subroutine get_token(string, num, token)

  implicit none

! Passed arguments
  
  character(*), intent(in)  :: string  ! The string to parse
  character(*), intent(out) :: token   ! The token to return

  integer, intent(in)       :: num     ! Which token to return

! Local variables

  integer   :: istart
  integer   :: iend
  integer   :: string_loc
  integer   :: token_count

  ! Uncomment the below chunk of code for a "safe" get_num_tokens at the 
  ! expense of calling get_num_tokens() each time a specific token is
  ! pulled from the string. When it's commented out, token will just be
  ! a blank string upon return

! integer   :: num_tokens

! call get_num_tokens(string, num_tokens)

! if (num .gt. num_tokens)
!   write(mdout, *) ' Error in get_token: Looking for more tokens than &
!                     &there are in string'
!   call mexit(6,1)
! end if

  ! Now get the num'th token

  token_count = 0
  istart = 1
  iend = len_trim(string)
  token = ' '

  do while (istart .le. iend)

    if (string(istart:istart) .le. ' ') then

      istart = istart + 1
      
    else

      do string_loc = istart, iend
        if ( string(string_loc:string_loc) .le. ' ' ) exit
      end do

      token_count = token_count + 1

      ! If this is the token we want, store it and return
      if ( token_count .eq. num ) then
        token = string(istart:string_loc-1)
        return
      end if

      istart = string_loc ! Move to the next token

    end if

  end do

end subroutine get_token

end module sander_lib
