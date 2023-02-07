function se_switch(x, half, skin)
  implicit none
  double precision :: se_switch

  double precision :: x
  double precision :: half
  double precision :: skin

  double precision :: xal
  double precision :: xbl
  double precision :: xar
  double precision :: xbr
  double precision :: y

  xal = -half+skin+1e-10
  xbl = -half-1e-10
  xar = half-skin+1e-10
  xbr = half-1e-10

! 3rd order polynomial interpolation between [xbl;xal]
!   and [xar;xbr]
! (first derivatives is zero at xbl, xal, xbr, xar)
!
!   switch
!     ^
!     |   xbl   xal         xar   xbr
!     |    |     |           |     |
!  1  |    |     |-----------|     |
!     |    |    /|           |\    |
!     |    |   / |           | \   |
!     |    |  /  |           |  \  |
!     |    | /   |           |   \ |
!     |    |/    |           |    \|
!  0  |----|     |           |     |-------
!     |    |     |           |     |
!-----+-------------------------------------------------------> distance
!

  if (x <= xbl) then
    se_switch = 0.0d0
  else if (x <= xal) then
    se_switch = ( 2.0d0*x*x*x &
                 -3.0d0*(xal+xbl)*x*x &
                 +6.0d0*xal*xbl*x &
                 +xbl*xbl*xbl &
                 -3.0d0*xal*xbl*xbl )/((xbl-xal)*(xbl-xal)*(xbl-xal))
    y = (x*x-xal*xal)/(xbl*xbl-xal*xal)
    se_switch = 1 + y*y*(2*y-3)
  else if (x <= xar) then
    se_switch = 1.0d0
  else if (x <= xbr) then
    se_switch = ( 2.0d0*x*x*x &
                 -3.0d0*(xar+xbr)*x*x &
                 +6.0d0*xar*xbr*x &
                 +xbr*xbr*xbr &
                 -3.0d0*xar*xbr*xbr )/((xbr-xar)*(xbr-xar)*(xbr-xar))
    y = (x*x-xar*xar)/(xbr*xbr-xar*xar)
    se_switch = 1 + y*y*(2*y-3)
  else
    se_switch = 0.0d0
  endif

! se_switch=max(se_switch,0.10d0)
end function se_switch
