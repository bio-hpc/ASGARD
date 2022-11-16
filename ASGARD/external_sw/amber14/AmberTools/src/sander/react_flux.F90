#include "../include/dprec.fh"

! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
! #|  FLUX_REACT: reactive flux driver                               |#
! #+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

   subroutine react_flux ( x, v, f, mass_inv, temp_fact, dt5, dtx &
                         , natm, nstep, nstlim )

   use wigner, only: rflux
   use wigner, only: init_rflux, back_prop, p_ndx, RC_ndx, nsample &
                   , RC_RS_toler, RC_PS_toler, RC_rmsd_toler, p_space

   implicit none

   integer, intent(in)    :: natm
   _REAL_ , intent(in)    :: mass_inv(natm), temp_fact
   _REAL_ , intent(out)   :: x(natm*3), v(natm*3), f(natm*3)
   _REAL_ , intent(inout) :: dt5, dtx, nstep, nstlim

   !  +---------------------------------------------------------------+

   integer :: n, npt
   integer, save :: RC_ndx_back, RC_ndx_all

   _REAL_  :: pflux 
   _REAL_, save  :: RC_stat(2)
   _REAL_, external :: ddot, dnrm2
   intrinsic :: abs

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  |  Initialize by finding velocity normal to dividing surface    |
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   if( init_rflux ) then

      if( back_prop ) then

         call setvel ( natm, v, mass_inv, temp_fact, 0, 0, 0.0 )

         pflux = ddot( natm*3, v, 1, rflux%normal_vect, 1 ) &
               / dnrm2( natm*3, rflux%normal_vect, 1 )

         write(6,*) 'initial momentum flux = ', pflux 

         if( pflux < 0.0d0 ) then
            pflux = - pflux
            rflux%normal_velv(:) = - v(:)
         endif

         dt5 = - abs( dt5 )
         dtx = - abs( dtx )

         p_ndx  = p_ndx  + 1

      else 

         dt5 = abs( dt5 )
         dtx = abs( dtx )

      endif

      x(:) = rflux%coord_TS(:)
      v(:) = rflux%normal_velv(:)
      f(:) = rflux%force_TS(:)

      init_rflux = .false.
      RC_ndx = 0
      rflux%RC_trj(:) = 0.0d0

   endif

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  |  Check stopping criteria                                      |
!  +:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::+
!  |  [backward propagation] committed to a FE well?  then reorder |
!  |     RC_sampled and set logical for forward propagation        |
!  |                                                               |
!  |  [forward propagation] committed to a FE well?  then analyze  |
!  |     reactive flux trajectory                                  |
!  +---------------------------------------------------------------+

   if( mod(RC_ndx,nsample) == 0 .and. RC_ndx /= 0 ) then

      npt = RC_ndx - nsample

      call RC_armsd( RC_stat, rflux%RC_trj(npt:npt+nsample-1), nsample )

!  :,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:
!  :  Backward propagation                                         :
!  :''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''':

      if( back_prop ) then 

         if( RC_stat(1) < RC_RS_toler .and. &
             RC_stat(2) < RC_rmsd_toler       ) then

            do n = 1, RC_ndx

               rflux%RC_sampled(n) = rflux%RC_trj( RC_ndx - n + 1)

            enddo

            back_prop = .false. 
            init_rflux = .true.

            RC_ndx_back = RC_ndx

!  .................................................................
!  :  Backward propagation ended up in product space               :
!  `````````````````````````````````````````````````````````````````
         else if( RC_stat(1) > RC_PS_toler .and. &
                  RC_stat(2) < RC_rmsd_toler       ) then

            do n = 1, RC_ndx

               rflux%RC_sampled(n) = rflux%RC_trj( RC_ndx + 1 - n)

            enddo

            RC_ndx_all = RC_ndx

            call kappa_keck ( p_ndx, RC_ndx_all )

            init_rflux = .true.
            back_prop  = .true.

         endif

!  :,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:
!  :  Forward propagation                                          :
!  :''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''':

      else 

         if( RC_stat(1) > RC_PS_toler .and. &
             RC_stat(2) < RC_rmsd_toler       ) then

            do n = 1, RC_ndx - 1

               rflux%RC_sampled(RC_ndx_back+n) = rflux%RC_trj(n+1)

            enddo

            RC_ndx_all = RC_ndx + RC_ndx_back - 1

            call kappa_keck ( p_ndx, RC_ndx_all )

            init_rflux = .true.
            back_prop  = .true.

!  .................................................................
!  :  Forward propagation ended up in reactant space               :
!  `````````````````````````````````````````````````````````````````
         else if( RC_stat(1) < RC_RS_toler .and. &
                  RC_stat(2) < RC_rmsd_toler       ) then

            do n = 1, RC_ndx - 1

               rflux%RC_sampled(RC_ndx_back+n) = rflux%RC_trj(n+1)

            enddo

            RC_ndx_all = RC_ndx + RC_ndx_back - 1

            call kappa_keck ( p_ndx, RC_ndx_all )

            init_rflux = .true.
            back_prop  = .true.

         endif 

      endif 

   endif 

   if( init_rflux .and. back_prop ) then

      write(888,*)
      write(888,*) 'p_ndx = ', p_ndx

      do n = 1, RC_ndx_all 
         write(888,*) n, rflux%RC_sampled(n)
      enddo

   endif 

   RC_ndx = RC_ndx + 1

   if( p_ndx > p_space ) then 

      do n = 1, p_space
         write(999,*) n, rflux%xi(n)
      enddo

      nstep = nstlim  ! initiate stop

   endif 

!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x
!  x+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++x

   end subroutine react_flux


!! ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
!! ;  Compute average rmsd of RC                                      ;
!! ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,

   subroutine RC_armsd ( RC_stat, RC, n )

   implicit none

   integer, intent(in) :: n
   _REAL_ , intent(in) :: RC(n)
   _REAL_ :: RC_stat(2)

   !.................................................................

   integer :: nn
   intrinsic :: sum, dble, sqrt
   _REAL_  :: avrg, add
  

   avrg = sum( RC(:) ) / dble(n)

   add = 0.0d0

   do nn = 1, n

      add = add + ( avrg - RC(nn) )**2

   enddo 

   RC_stat(1) = avrg
   RC_stat(2) = sqrt( add / dble(n) )

   end subroutine RC_armsd 

