! <compile=optimized>
!------------------------------------------------------------------------------
! nose_hoover.f
!------------------------------------------------------------------------------
!
! Integrator of Nose'-Hoover chain of thermostat
!
! Ref: Jang & Voth, J. Chem. Phys. 107, 9514 (1997)
!
! Last update: 03/09/2006
!
!------------------------------------------------------------------------------
#include "../include/dprec.fh"


!==============================================================================
module nose_hoover_module
!==============================================================================

  use abfqmmm_module, only: abfqmmm_param

  implicit none

  character(*), parameter :: module_name = "nose_hoover_module:"
  logical, save :: module_init = .false.
  integer, parameter :: dp = 8   !! double precision

  integer, parameter :: M = 9   !! max. chain length

  type SystemCoordinate_type
     real(dp) :: mass
     real(dp), pointer :: vel
     real(dp), pointer :: frc
  endtype

  type Thermostat_type
     real(dp) :: eta    ( M )   !! position (heat bath variable)
     real(dp) :: eta_old( M )   !! old position for Leap-frog
     real(dp) :: v      ( M )   !! velocity (thermostat variable)
     real(dp) :: v_old  ( M )   !! old velocity for Leap-frog
     real(dp) :: a      ( M )   !! acceleration (without friction term)
     real(dp) :: Q      ( M )   !! mass
     real(dp) :: Q_qm   ( M )   !! lam81 (what is this?)
     real(dp) :: Q_inv  ( M )   !! inverse mass
     real(dp) :: Q_inv_qm  ( M ) !! lam81 (what is this?)
     real(dp) :: gamma_ln     !! APJ: Langevin parameter for chains 
     real(dp) :: gamma_ln_qm  ! lam81 (what is this?)
     real(dp) :: E_ln        !! APJ: tracks changes to energy / extended energy 
                             !! APJ: from the langevin process jumping from one 
                             !! APJ: deterministic trajectory to another 
     real(dp) :: kT, Ndof_kT   !! kB * temperature
     real(dp) :: Ekin2       !! APJAPJ: kinetic energy (x2) ON THE FULL STEP
     type( SystemCoordinate_type ), pointer :: system_coords( : )
     integer :: num_system_coords
     integer :: system_coord_id = 0
     logical :: activated   !! flag to activate the thermostat
     logical :: initialized = .false.
     logical :: langevin = .false.
     logical :: adaptive = .false.
     integer :: iatom
  endtype

  !..................................................

  integer, save :: print_level = 0   !! output is suppressed
!!  integer :: print_level = 1

  !..................................................

!! Random number generator from Numerical Recipes.

  integer, save :: random_number_seed = 100
  integer, save :: idum_ran2

  !..................................................

  type(Thermostat_type), save :: thermo_lnv

  _REAL_, save :: c2_lnv, mass_lnv, v_lnv, f_lnv_v, f_lnv_p, x_lnv, x_lnv_old


!! Nose Hoover vars (originally in nose_hoover_vars.F90 module)
  logical :: use_nose_hoover = .false.

  integer :: file_nhc = 1001

  integer :: nchain, nthermo

  type( Thermostat_type ), allocatable :: thermo( :, : )

  _REAL_ :: Econserved = 0.d0

  _REAL_ :: tau  !! characteristic time scale of the system
  _REAL_ :: tau_qm
 
contains

!------------------------------------------------------------------------------
subroutine nose_hoover_module_init
!------------------------------------------------------------------------------

  implicit none


  if ( module_init ) return

  if ( print_level > 0 ) then
     write(6,*)
     write(6,*) "Initializing ", module_name
  endif

  module_init = .true.

!!
!! Output.
!!

  if ( print_level > 0 ) then
     write(6,*)
     write(6,*) module_name
     write(6,*)
     write(6,*) "   chain length = ", M
  endif

end subroutine

!------------------------------------------------------------------------------
subroutine Thermostat_init  &
     ( nchain, thermo, num_system_coords, num_constraints, kT, tau, &
       langevin, adaptive, pos_init, vel_init, activate, tau_qm, iatom )  ! APJ lam81

! initializes a Thermostat object.
!------------------------------------------------------------------------------

  use random, only : gauss

  implicit none

  type( Thermostat_type ) :: thermo
  integer, intent(in) :: nchain  !! number of oscillators in each chain
  integer, intent(in) :: num_system_coords   !! number of system coordinates
                                             !! to be thermostatted
  integer, intent(in) :: num_constraints   !! number of geometrical constraints
                                           !! in the system
  real(dp), intent(in) :: kT   !! external temperature
  real(dp), intent(in) :: tau   !! characteristic timescale of the system
  logical, intent(in) :: langevin !! Nose-Hoover-Chain -> NHCLangevin?
                                       !! ( NHL = NHCL with nchain=1 )
  logical, intent(in) :: adaptive !! extra 'adaptive' thermostat ?
                                  !! ( n+1-th thermostat connects parallel to 1st )
  real(dp), intent(in), optional :: pos_init, vel_init   !! initial condition
                                                         !! of the thermostat

  logical, intent(in), optional :: activate   !! flag to activate the thermostat
  real(dp), intent(in), optional :: tau_qm
  integer, intent(in), optional :: iatom

  integer :: Ndof, j

  if ( .not. module_init ) call nose_hoover_module_init

  if (adaptive .and. (nchain+1) > (M-1) ) then
     write(6,*) "Stop. Number of oscillators too large."
     write(6,*) "In adaptive case, nchain must be <= ", (M-2)
     stop
  else if ( nchain > (M-1) ) then
     write(6,*) "Stop. Number of oscillators too large."
     write(6,*) "nchain must be <= ", (M-1)
     stop
  end if

  thermo%initialized = .true.

  allocate( thermo%system_coords ( num_system_coords ) )
  thermo%num_system_coords = num_system_coords

  Ndof = num_system_coords - num_constraints

!! Set up thermostat masses / coupling constants, target temperature.

  thermo%Q(    1     ) = kT * tau**2 * dble( Ndof )
  thermo%Q( 2:nchain ) = kT * tau**2
  thermo%Q_inv( 1:nchain ) = 1.0d0 / thermo%Q( 1:nchain )

  thermo%Q( nchain+1:M ) = 0.d0
  thermo%Q_inv( nchain+1:M ) = 0.d0

  if(abfqmmm_param%abfqmmm == 1) then
    thermo%iatom = iatom

    thermo%Q_qm(    1     ) = kT * tau_qm**2 * dble( Ndof )
    thermo%Q_qm( 2:nchain ) = kT * tau_qm**2
    thermo%Q_inv_qm( 1:nchain ) = 1.0d0 / thermo%Q_qm( 1:nchain )

    thermo%Q_qm( nchain+1:M ) = 0.d0
    thermo%Q_inv_qm( nchain+1:M ) = 0.d0
  end if

  if (adaptive) then
    thermo%Q( nchain+1 ) = kT * tau**2
    thermo%Q_inv( nchain+1 ) = 1.0d0 / thermo%Q( nchain+1 )
    if(abfqmmm_param%abfqmmm == 1) then
      thermo%Q_qm( nchain+1 ) = kT * tau_qm**2
      thermo%Q_inv_qm( nchain+1 ) = 1.0d0 / thermo%Q_qm( nchain+1 )
    end if
  endif

  if (langevin) then
    thermo%gamma_ln = 1.0d0 / tau
    if(abfqmmm_param%abfqmmm == 1) thermo%gamma_ln_qm = 1.0d0 / tau_qm
    thermo%E_ln = 0.d0
  endif

  thermo%kT = kT
  thermo%Ndof_kT = dble( Ndof ) * kT
  thermo%Ekin2 = 0.d0 ! APJAPJ

!!
!! Initial position, velocity, and acceleration.
!!

  thermo%eta( : ) = 0.0d0
  thermo%eta_old( : ) = 0.0d0

  thermo%v( : ) = 0.0d0
  thermo%v_old( : ) = 0.0d0

  thermo%a( : ) = 0.0d0
      
  do j = 1, nchain
    if(abfqmmm_param%abfqmmm /= 1) then
      call gauss( 0.d0, sqrt(kT/thermo%Q(j)), thermo%v(j) )
    else
      if(abfqmmm_param%id(thermo%iatom) <= 4) then
        call gauss( 0.d0, sqrt(kT/thermo%Q_qm(j)), thermo%v(j) )
      else
        call gauss( 0.d0, sqrt(kT/thermo%Q(j)), thermo%v(j) )
      end if
    end if
    thermo%v_old( j ) = thermo%v( j )
  enddo
  if (adaptive) then
    j = nchain + 1
    if(abfqmmm_param%abfqmmm /= 1) then
      call gauss( 0.d0, sqrt(kT/thermo%Q(j)), thermo%v(j) )
    else
      if(abfqmmm_param%id(thermo%iatom) <= 4) then
        call gauss( 0.d0, sqrt(kT/thermo%Q_qm(j)), thermo%v(j) )
      else
        call gauss( 0.d0, sqrt(kT/thermo%Q(j)), thermo%v(j) )
      end if
    end if
    thermo%v_old( j ) = thermo%v( j )
  endif

  if ( present( pos_init ) ) thermo%eta( : ) = pos_init
  if ( present( vel_init ) ) thermo%v  ( : ) = vel_init

  do j = 2, nchain
    if(abfqmmm_param%abfqmmm /= 1) then
      thermo%a( j ) = thermo%Q_inv( j )  &
         * ( thermo%Q( j-1 ) * thermo%v( j-1 )**2 - thermo%kT )
    !! Note: thermo%a( 1 ) is computed in Thermostat_integrate().
    else
      if(abfqmmm_param%id(thermo%iatom) <= 4) then
        thermo%a( j ) = thermo%Q_inv_qm( j )  &
           * ( thermo%Q_qm( j-1 ) * thermo%v( j-1 )**2 - thermo%kT )
      else
        thermo%a( j ) = thermo%Q_inv( j )  &
           * ( thermo%Q( j-1 ) * thermo%v( j-1 )**2 - thermo%kT )
      end if
    end if
  enddo

!! APJ: Special flags
  thermo%langevin = langevin
  thermo%adaptive = adaptive

!! Flag to activate the thermostat.

  if ( present( activate ) ) then
     thermo%activated = activate
  else
     thermo%activated = .true.
  endif

end subroutine

!------------------------------------------------------------------------------
subroutine Thermostat_link ( thermo, system_mass, system_velocity, system_force)

! establishes the link to system velocicities via pointer assignment.
!------------------------------------------------------------------------------

  implicit none

  type( Thermostat_type ) :: thermo
  real(dp), intent(in) :: system_mass
  real(dp), intent(in), target :: system_velocity
  real(dp), intent(in), target :: system_force

  integer :: id

  if ( .not. thermo%initialized ) then
     write(6,*) module_name, " thermostat is not initialized. stop."
     stop
  endif

  thermo%system_coord_id = thermo%system_coord_id + 1
  id = thermo%system_coord_id
  if ( id <= thermo%num_system_coords ) then
     thermo%system_coords( id )%mass = system_mass
     thermo%system_coords( id )%vel => system_velocity   !! establish link
     thermo%system_coords( id )%frc => system_force  ! APJ 
     thermo%Ekin2 = thermo%Ekin2 + system_mass*system_velocity**2 ! APJAPJ
  else
     write(6,*) module_name, " thermo%num_system_coords is too small. stop."
     stop
  endif

end subroutine

!------------------------------------------------------------------------------
subroutine Thermostat_switch ( thermo, activate )
!
! Switch thermostat on or off.
!------------------------------------------------------------------------------

  implicit none

  type( Thermostat_type ) :: thermo
  logical, intent(in) :: activate

  thermo%activated = activate
  thermo%v(:) = 0.d0
  thermo%v_old(:) = 0.d0

end subroutine

!------------------------------------------------------------------------------
subroutine Thermostat_integrate_1  &
           ( nchain, thermostats, num_thermostats, timestep, ntp )
!------------------------------------------------------------------------------

  implicit none

  integer, intent(in) :: nchain
  integer, intent(in) :: num_thermostats
  type( Thermostat_type ), target :: thermostats( num_thermostats )
  real(dp), intent(in) :: timestep

  integer :: ithermo, j, k, ntp, istart3,iend3
  real(dp) :: dt, hdt, kT, Ndof_kT, Ekin2, exp1, exp2
  type( Thermostat_type ), pointer :: thermo
  type( SystemCoordinate_type ), pointer :: system( : )
#ifdef MPI
# include "parallel.h"
#endif


  dt = timestep
  hdt = 0.5d0 * dt

  !------------------------------
#ifdef MPI
  istart3 = iparpt3(mytaskid)+1
  iend3   = iparpt3(mytaskid+1)
#else
  istart3 = 1
  iend3   = num_thermostats
#endif

  do ithermo = 1, num_thermostats
     thermo => thermostats( ithermo )
     if ( .not.thermo%activated ) cycle
     !! Error check.
     if ( .not. thermo%initialized ) goto 9000
     if ( thermo%num_system_coords /= thermo%system_coord_id ) goto 9100
     system => thermo%system_coords( : )
     kT      = thermo%kT
     Ndof_kT = thermo%Ndof_kT
     !------------------------------
     !! Update coordinates for oscillators of Nose'-Hoover chains.
     !! - odd positions and forces.
     !! - even velocities.
     do j = 1, nchain, 2
        thermo%eta_old( j ) = thermo%eta( j )
        thermo%eta( j ) = thermo%eta( j ) + dt * thermo%v( j )
     end do
     do j = 2, nchain, 2
        thermo%v_old( j ) = thermo%v( j )
        exp1 = 1.d0
        if ( j < nchain ) exp1 = exp( - hdt * thermo%v( j+1 ) )
        exp2 = exp1*exp1
        thermo%v( j ) = thermo%v( j ) * exp2  &
           + dt * thermo%a( j ) * exp1
     end do
     do j = 3, nchain, 2
        thermo%a( j ) = thermo%Q_inv( j )  &
           * ( thermo%Q( j-1 ) * thermo%v( j-1 )**2 - kT )
     end do
     !! Update force acting on the first thermostat oscillator.
     Ekin2 = 0.0d0

     do k = 1, thermo%num_system_coords
        Ekin2 = Ekin2 + system( k )%mass * system( k )%vel**2
     enddo
     thermo%a( 1 ) = thermo%Q_inv( 1 ) * ( Ekin2 - Ndof_kT )
  enddo  !! loop over thermostats

  if( ntp>0 ) then
  !------------------------------
     kT = thermo_lnv%kT
     do j = 1, nchain, 2
        thermo_lnv%eta_old( j ) = thermo_lnv%eta( j )
        thermo_lnv%eta( j ) = thermo_lnv%eta( j ) + dt * thermo_lnv%v( j )
     end do

     do j = 2, nchain, 2
        thermo_lnv%v_old( j ) = thermo_lnv%v( j )
        exp1 = 1.d0
        if ( j < nchain ) exp1 = exp( -hdt * thermo_lnv%v( j+1 ) )
        exp2 = exp1*exp1
        thermo_lnv%v( j ) = thermo_lnv%v( j ) * exp2  &
           + dt * thermo_lnv%a( j ) * exp1
     end do

     do j = 3, nchain, 2
        thermo_lnv%a( j ) = thermo_lnv%Q_inv( j )  &
           * ( thermo_lnv%Q( j-1 ) * thermo_lnv%v( j-1 )**2 - kT )
     end do
     !! Update force acting on the first thermostat oscillator.
     Ekin2 = 0.0d0
     thermo_lnv%a(1) = thermo_lnv%Q_inv(1) * (mass_lnv*v_lnv*v_lnv-kT )
  end if

  return

!!
!! Error handling.
!!

9000 continue
  write(6,*) module_name, " thermostat is not initialized. stop."
  stop

9100 continue
  write(6,*) module_name, " link to system velocities is incomplete. stop."
  stop

end subroutine

!------------------------------------------------------------------------------
subroutine Thermostat_integrate_2  &
           ( nchain, thermostats, num_thermostats, timestep, ntp )
!------------------------------------------------------------------------------

  implicit none

  integer, intent(in) :: nchain,ntp
  integer, intent(in) :: num_thermostats
  type( Thermostat_type ), target :: thermostats( num_thermostats )
  real(dp), intent(in) :: timestep

  integer :: ithermo, j, istart3,iend3
  real(dp) :: dt, hdt, kT, Ndof_kT, exp1, exp2
  type( Thermostat_type ), pointer :: thermo
  type( SystemCoordinate_type ), pointer :: system( : )

#ifdef MPI
# include "parallel.h"
#endif

  dt = timestep
  hdt = 0.5d0 * dt
  thermo => thermostats(1)

  !------------------------------
#ifdef MPI
  istart3 = iparpt3(mytaskid)+1
  iend3   = iparpt3(mytaskid+1)
#else
  istart3 = 1
  iend3   = num_thermostats
#endif

  do ithermo = 1, num_thermostats
     thermo => thermostats( ithermo )
     if ( .not.thermo%activated ) cycle
     system => thermo%system_coords( : )
     kT      = thermo%kT
     Ndof_kT = thermo%Ndof_kT

     !------------------------------

     !! Update coordinates for oscillators of Nose'-Hoover chains.
     !! - even positions and forces.
     !! - odd velocities.

     do j = 2, nchain, 2
        thermo%eta( j ) = thermo%eta( j ) + dt * thermo%v( j )
        thermo%eta_old( j ) = thermo%eta( j )
     end do

     do j = 1, nchain, 2
        exp1 = 1.d0
        if ( j < nchain ) exp1 = exp( - hdt * thermo%v( j+1 ) )
        exp2 = exp1**2
        thermo%v( j ) = thermo%v( j ) * exp2  &
           + dt * thermo%a( j ) * exp1
        thermo%v_old( j ) = thermo%v( j )
     end do

     do j = 2, nchain, 2
        thermo%a( j ) = thermo%Q_inv( j )  &
           * ( thermo%Q( j-1 ) * thermo%v( j-1 )**2 - kT )
     end do

  enddo  !! loop over thermostats

  if( ntp > 0 ) then
     kT      = thermo%kT
     do j = 2, nchain, 2
        thermo_lnv%eta( j ) = thermo_lnv%eta( j ) + dt * thermo_lnv%v( j )
        thermo_lnv%eta_old( j ) = thermo_lnv%eta( j )
     end do

     do j = 1, nchain, 2
        exp1 = 1.d0
        if ( j < nchain ) exp1 = exp( - hdt * thermo_lnv%v( j+1 ) )
        exp2 = exp1**2

        thermo_lnv%v( j ) = thermo_lnv%v( j ) * exp2  &
           + dt * thermo_lnv%a( j ) * exp1

        thermo_lnv%v_old( j ) = thermo_lnv%v( j )
     end do

     do j = 2, nchain, 2
        thermo_lnv%a( j ) = thermo_lnv%Q_inv( j )  &
           * ( thermo_lnv%Q( j-1 )*thermo_lnv%v( j-1 )**2 - kT )
     end do

  end if
  !------------------------------
end subroutine

!------------------------------------------------------------------------------
!  Update the thermostats (numerically integrate equations of motion)
!------------------------------------------------------------------------------
subroutine Adaptive_Thermostat_integrate &
           ( nchain, thermostats, num_thermostats, timestep, ntp, step)
!------------------------------------------------------------------------------

  use random, only : gauss

  implicit none

  integer, intent(in) :: nchain
  type( Thermostat_type ), target :: thermostats( num_thermostats )
  integer, intent(in) :: num_thermostats
  real(dp), intent(in) :: timestep
  integer, intent(in) :: ntp
  integer, intent(in) :: step

  integer :: ithermo, j, k,  istart3,iend3, ibth,ithm,iacc
  real(dp) :: dt, hdt, qdt, kT, Ndof_kT, Ekin2, exp1, acc,thm
  real(dp) :: exp_ln,sig_ln,gss_ln, Ekin2_old,Ekin2_new
  logical :: langevin, adaptive
  type( Thermostat_type ), pointer :: thermo
  type( SystemCoordinate_type ), pointer :: system( : )
#ifdef MPI
# include "parallel.h"
#endif


  dt = timestep
  hdt = 0.5d0 * dt
  qdt = 0.25d0 * dt

  !------------------------------
#ifdef MPI
  istart3 = iparpt3(mytaskid)+1
  iend3   = iparpt3(mytaskid+1)
#else
  istart3 = 1
  iend3   = num_thermostats
#endif

  istart3 = 1
  iend3   = num_thermostats

  !------------------------------
  ! These move with the Leapfrog velocities:
  !  even-numbered thermostats
  !  odd-numbered baths and accelerations
  if(step.eq.1) then
    ithm = 2
    iacc = 3 ! 1st thermo is 'accelerated' by system,(effectively the 'zero-th thermo')
    ibth = 1
  !------------------------------
  ! These move with the Leapfrog positions:
  !  odd-numbered thermostats
  !  even-numbered baths and accelerations
  else if (step.eq.2) then
    ithm = 1
    iacc = 2
    ibth = 2
  else
    write(6,*) "bad call to Thermostat_integrate. stop."
    stop
  end if

  do ithermo = istart3, iend3
     thermo => thermostats( ithermo )
     system => thermo%system_coords( : )
     !! Error check.
     if ( .not. thermo%initialized ) goto 9000
     if ( thermo%num_system_coords /= thermo%system_coord_id ) goto 9100
     kT      = thermo%kT
     Ndof_kT = thermo%Ndof_kT
     adaptive = thermo%adaptive
     langevin = thermo%langevin

     if(step.eq.1) then
        !------------------------------
        !! If inactivated, update the velocities only with Newtonian forces.
        if ( .not. thermo%activated ) then
           do k = 1, thermo%num_system_coords
              acc = system( k )%frc / system( k )%mass
              system( k )%vel = system(k)%vel + dt*acc
           enddo
           cycle   ! skip to next thermostat
        !------------------------------
        !! Simple Langevin, Adaptive Langevin.
        else if(nchain.eq.0 .and. langevin) then
           Ekin2_old = 0.d0
           Ekin2_new = 0.d0
!APJAPJ    Ekin2 = 0.d0
           do k = 1, thermo%num_system_coords
              acc = system( k )%frc / system( k )%mass
              if(abfqmmm_param%abfqmmm /= 1) then
                exp_ln = exp(-dt*thermo%gamma_ln)
              else
                if(abfqmmm_param%id(thermo%iatom) <= 4) then
                  exp_ln = exp(-dt*thermo%gamma_ln_qm)
                else
                  exp_ln = exp(-dt*thermo%gamma_ln)
                end if
              end if
              sig_ln = sqrt( (1.d0-exp_ln*exp_ln)*kT/system(k)%mass )
              call gauss(0.d0,sig_ln,gss_ln)
              exp1 = 1.d0
              if(adaptive) exp1 = exp(-hdt*thermo%v(nchain+1))
!APJAPJ       system( k )%vel = ( ( system(k)%vel*exp1 ) + hdt*acc )*exp1  ! half a step on old trajectory
              system( k )%vel = system(k)%vel + hdt*acc                    ! half timestep of Newtonian dynamics ! APJAPJ
              system( k )%vel = system(k)%vel*exp1                         ! half timestep of thermostatting on old trajectory ! APJAPJ 
              Ekin2_old = Ekin2_old + system(k)%mass*system(k)%vel**2      ! save old energy
              system( k )%vel = system(k)%vel*exp_ln + gss_ln              ! langevin leap to new trajectory
              Ekin2_new = Ekin2_new + system(k)%mass*system(k)%vel**2      ! save new energy
!APJAPJ       system( k )%vel = ( ( system(k)%vel*exp1 ) + hdt*acc )*exp1  ! half a step on new trajectory
!APJAPJ       Ekin2 = Ekin2 + system(k)%mass*system(k)%vel**2  ! save final energy
              system( k )%vel = system(k)%vel*exp1                         ! half timestep of thermostatting on new trajectory ! APJAPJ
              system( k )%vel = system(k)%vel + hdt*acc                    ! half timestep of Newtonian dynamics ! APJAPJ
           enddo
           thermo%Ekin2 = 0.5d0*(Ekin2_old+Ekin2_new)  ! APJAPJ
           if(abfqmmm_param%abfqmmm /= 1) then
             if(adaptive) thermo%a( nchain+1 ) = thermo%Q_inv( nchain+1 ) * ( thermo%Ekin2 - Ndof_kT )        ! APJ ! APJAPJ
           else
             if(abfqmmm_param%id(thermo%iatom) <= 4) then
               if(adaptive) thermo%a( nchain+1 ) = thermo%Q_inv_qm( nchain+1 ) * ( thermo%Ekin2 - Ndof_kT ) ! APJAPJ
             else
               if(adaptive) thermo%a( nchain+1 ) = thermo%Q_inv( nchain+1 ) * ( thermo%Ekin2 - Ndof_kT ) ! APJAPJ
             end if 
           end if
           thermo%E_ln = thermo%E_ln + 0.5d0*( Ekin2_new - Ekin2_old )
        else
           !------------------------------
           !! If activated, and not direct Langevin,
           !! update the velocities with Newtonian forces
           !! and chain thermostat. Update the feedback.
           thm = 0.d0
           if(nchain.ge.1) thm = thm + thermo%v(1)
           if(adaptive) thm = thm + thermo%v(nchain+1)
           if(nchain<1 .and. .not.adaptive) then
             write(6,*) 'One of the thermostats isnt doing anything. Stopping ...'
             stop
           endif
           exp1 = exp( - hdt*thm )
           Ekin2 = 0.0d0
           do k = 1, thermo%num_system_coords
              acc = system( k )%frc / system( k )%mass
! APJAPJ      system( k )%vel = ( ( system(k)%vel*exp1 ) + dt*acc )*exp1
              system( k )%vel = system(k)%vel + hdt*acc              ! half timestep of Newtonian dynamics ! APJAPJ
              system( k )%vel = system(k)%vel * exp1                 ! half timestep of thermostatting     ! APJAPJ
              Ekin2 = Ekin2 + system( k )%mass * system( k )%vel**2
              system( k )%vel = system(k)%vel * exp1                 ! half timestep of thermostatting     ! APJAPJ
              system( k )%vel = system(k)%vel + hdt*acc              ! half timestep of Newtonian dynamics ! APJAPJ
           enddo
           thermo%Ekin2 = Ekin2 ! APJAPJ
           if(abfqmmm_param%abfqmmm /= 1) then
             thermo%a( 1 ) = thermo%Q_inv( 1 ) * ( Ekin2 - Ndof_kT )
           else
             if(abfqmmm_param%id(thermo%iatom) <= 4) then
               thermo%a( 1 ) = thermo%Q_inv_qm( 1 ) * ( Ekin2 - Ndof_kT )
             else
               thermo%a( 1 ) = thermo%Q_inv( 1 ) * ( Ekin2 - Ndof_kT )
             end if
           end if
           if(abfqmmm_param%abfqmmm /= 1) then
             if(adaptive) thermo%a( nchain+1 ) = thermo%Q_inv( nchain+1 ) * ( Ekin2 - Ndof_kT )
           else
             if(abfqmmm_param%id(thermo%iatom) <= 4) then
               if(adaptive) thermo%a( nchain+1 ) = thermo%Q_inv_qm( nchain+1 ) * ( Ekin2 - Ndof_kT )
             else
               if(adaptive) thermo%a( nchain+1 ) = thermo%Q_inv( nchain+1 ) * ( Ekin2 - Ndof_kT )
             end if
           end if
        endif
     endif

     !------------------------------
     !! Update thermos and heat-baths of Nose'-Hoover chains.
     do j = ithm, nchain, 2
        thermo%v_old( j ) = thermo%v( j )
        if ( j < nchain ) then
           exp1 = exp( - hdt * thermo%v(j+1) )
           thermo%v( j ) = ( ( thermo%v(j)*exp1 ) + dt*thermo%a(j) )* exp1
        else if ( langevin .and. j==nchain ) then
          if(abfqmmm_param%abfqmmm /= 1) then
            exp_ln = exp( -dt*thermo%gamma_ln )
            call gauss( 0.d0, sqrt((1.d0-exp_ln*exp_ln)*kT/thermo%Q(j)), gss_ln )
            thermo%v( j ) = thermo%v(j) + hdt*thermo%a(j)  ! half a step
            Ekin2_old = thermo%Q(j)*thermo%v(j)**2         ! save energy
            thermo%v( j ) = thermo%v( j )*exp_ln + gss_ln  ! langevin jump
            Ekin2_new = thermo%Q(j)*thermo%v(j)**2         ! save energy
            thermo%v( j ) = thermo%v(j) + hdt*thermo%a(j)  ! half a step
            thermo%E_ln = thermo%E_ln + 0.5d0*( Ekin2_new - Ekin2_old )
          else
            if(abfqmmm_param%id(thermo%iatom) <= 4) then
              exp_ln = exp( -dt*thermo%gamma_ln_qm )
              call gauss( 0.d0, sqrt((1.d0-exp_ln*exp_ln)*kT/thermo%Q_qm(j)), gss_ln )
              thermo%v( j ) = thermo%v(j) + hdt*thermo%a(j)  ! half a step
              Ekin2_old = thermo%Q_qm(j)*thermo%v(j)**2      ! save energy
              thermo%v( j ) = thermo%v( j )*exp_ln + gss_ln  ! langevin jump
              Ekin2_new = thermo%Q_qm(j)*thermo%v(j)**2      ! save energy
              thermo%v( j ) = thermo%v(j) + hdt*thermo%a(j)  ! half a step
              thermo%E_ln = thermo%E_ln + 0.5d0*( Ekin2_new - Ekin2_old )
            else
              exp_ln = exp( -dt*thermo%gamma_ln )
              call gauss( 0.d0, sqrt((1.d0-exp_ln*exp_ln)*kT/thermo%Q(j)), gss_ln )
              thermo%v( j ) = thermo%v(j) + hdt*thermo%a(j)  ! half a step
              Ekin2_old = thermo%Q(j)*thermo%v(j)**2         ! save energy
              thermo%v( j ) = thermo%v( j )*exp_ln + gss_ln  ! langevin jump
              Ekin2_new = thermo%Q(j)*thermo%v(j)**2         ! save energy
              thermo%v( j ) = thermo%v(j) + hdt*thermo%a(j)  ! half a step
              thermo%E_ln = thermo%E_ln + 0.5d0*( Ekin2_new - Ekin2_old )
            end if
          end if
        else
          thermo%v( j ) = thermo%v(j) + dt*thermo%a(j)
        endif
     end do
     do j = iacc, nchain, 2
      if(abfqmmm_param%abfqmmm /= 1) then
        thermo%a( j ) = thermo%Q_inv( j ) * ( thermo%Q( j-1 ) * thermo%v( j-1 )**2 - kT )
      else
       if(abfqmmm_param%id(thermo%iatom) <= 4) then
        thermo%a( j ) = thermo%Q_inv_qm( j ) * ( thermo%Q_qm( j-1 ) * thermo%v( j-1 )**2 - kT )
       else
        thermo%a( j ) = thermo%Q_inv( j ) * ( thermo%Q( j-1 ) * thermo%v( j-1 )**2 - kT )
       end if
      end if
     end do
     do j = ibth, nchain, 2
        thermo%eta_old( j ) = thermo%eta( j )
        thermo%eta( j ) = thermo%eta( j ) + dt * thermo%v( j )
     end do
     if(adaptive) then
        j = nchain+1
        if(step==1) then
           thermo%eta_old( j ) = thermo%eta( j )
           thermo%eta( j ) = thermo%eta( j ) + dt * thermo%v( j )
        else if(step==2) then
           thermo%v_old( j ) = thermo%v( j )
           thermo%v( j ) = thermo%v(j) + dt*thermo%a(j)
        endif
     endif

  enddo  !! loop over thermostats

  if( ntp>0 ) then 
  !------------------------------
     kT = thermo_lnv%kT

     if(step.eq.1) then
        !------------------------------
        !! Update the velocities with Newtonian forces and thermostat forces
        !! Update the force acting on the first thermostat oscillator.
        exp1 = exp( -hdt * thermo_lnv%v(1) )
        acc = ( f_lnv_v + f_lnv_p ) / mass_lnv
        v_lnv = ( ( v_lnv*exp1 ) + dt*acc )*exp1
        thermo_lnv%a( 1 ) = thermo_lnv%Q_inv( 1 ) * (mass_lnv*v_lnv*v_lnv-kT ) ! APJ    
     endif

     do j = ithm, nchain, 2
        thermo_lnv%v_old( j ) = thermo_lnv%v( j )
        exp1 = 1.d0
        if ( j < nchain ) exp1 = exp( -hdt * thermo_lnv%v( j+1 ) )
        thermo_lnv%v( j ) = ( ( thermo_lnv%v(j)*exp1 ) + dt*thermo_lnv%a(j) )*exp1
     end do

     do j = iacc, nchain, 2
        thermo_lnv%a( j ) = thermo_lnv%Q_inv( j ) * ( thermo_lnv%Q( j-1 ) * thermo_lnv%v( j-1 )**2 - kT )
     end do

     do j = ibth, nchain, 2
        thermo_lnv%eta_old( j ) = thermo_lnv%eta( j )
        thermo_lnv%eta( j ) = thermo_lnv%eta( j ) + dt * thermo_lnv%v( j )
     end do

  end if

  return

!!
!! Error handling.
!!

9000 continue
  write(6,*) module_name, " thermostat is not initialized. stop."
  stop

9100 continue
  write(6,*) module_name, " link to system velocities is incomplete. stop."
  stop

end subroutine

!------------------------------------------------------------------------------
function Thermostat_hamiltonian ( nchain, thermostats, num_thermostats ) result ( E )

! computes the thermostat terms in the extended Hamiltonian.
!------------------------------------------------------------------------------

  implicit none

  integer, intent(in) :: nchain
  integer, intent(in) :: num_thermostats
  type( Thermostat_type ), intent(in), target :: thermostats( num_thermostats )
  real(dp) :: E, av_v, av_eta
  integer :: j, ithermo
  type( Thermostat_type ), pointer :: thermo

#ifdef MPI
   include 'mpif.h'
#  include "parallel.h"
  integer :: ierr
  _REAL_ :: Etmp
#endif
  integer :: istart3, iend3
#include "../include/md.h"
#include "../include/memory.h"

  E = 0.0d0

#ifdef MPI
  if( mpi_orig) then
     istart3 = 1
     iend3 = natom*3
  else
     istart3 = 3*iparpt(mytaskid)+1
     iend3 = 3*iparpt(mytaskid+1)
  end if
#else
  istart3 = 1
  iend3 = 3*natom
#endif

  do ithermo = istart3, iend3
     thermo => thermostats( ithermo )
     av_v   = 0.5d0 * ( thermo%v_old( 1 ) + thermo%v( 1 ) )
     av_eta = 0.5d0 * ( thermo%eta_old( 1 ) + thermo%eta( 1 ) )
     E = E + 0.5d0 * thermo%Q( 1 ) * av_v**2 + thermo%Ndof_kT * av_eta
     do j = 2, nchain
        av_v   = 0.5d0 * ( thermo%v_old( j ) + thermo%v( j ) )
        av_eta = 0.5d0 * ( thermo%eta_old( j ) + thermo%eta( j ) )
        E = E + 0.5d0 * thermo%Q( j ) * av_v**2 + thermo%kT * av_eta
     enddo
  enddo

#ifdef MPI
  call mpi_allreduce(E,Etmp,1,MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
  E=Etmp
#endif

  if(ntp>0) then
     do j=1,nchain
        av_v   = 0.5d0 * ( thermo_lnv%v_old(j) + thermo_lnv%v(j) )
        av_eta = 0.5d0 * ( thermo_lnv%eta_old(j) + thermo_lnv%eta_old(j) )
        E=E+0.5*thermo_lnv%Q(j)*av_v*av_v + thermo_lnv%kT*av_eta
     end do
  end if


end function

!------------------------------------------------------------------------------
function Adaptive_Thermostat_hamiltonian ( nchain, thermostats, num_thermostats ) result ( E )

! computes the thermostat terms in the extended Hamiltonian.
!------------------------------------------------------------------------------

  implicit none

  integer, intent(in) :: nchain
  integer, intent(in) :: num_thermostats
  type( Thermostat_type ), intent(in), target :: thermostats( num_thermostats )
  real(dp) :: E, v, eta
  integer :: j, ithermo
  type( Thermostat_type ), pointer :: thermo
  logical :: adaptive,langevin

#ifdef MPI
   include 'mpif.h'
#  include "parallel.h"
  integer :: ierr
  _REAL_ :: Etmp
#endif
  integer :: istart3, iend3
#include "../include/md.h"
#include "../include/memory.h"

  E = 0.0d0

#ifdef MPI
  if( mpi_orig) then
     istart3 = 1
     iend3 = natom*3
  else
     istart3 = 3*iparpt(mytaskid)+1
     iend3 = 3*iparpt(mytaskid+1) 
  end if
#else
  istart3 = 1
  iend3 = 3*natom
#endif

! istart3 = 1
! iend3   = num_thermostats

  !! APJ: Leapfrog is update x, update f, update v
  !! APJ: AMBER uses            update f, update v, update x
  !! APJ: Thus, x and all variables that move with x, e.g. therm%v(1),
  !! APJ: are 1 step ahead of the forces that were calculated
  !! APJ: Therefore, part of the thermostat energy is one step ahead,
  !! APJ: unless we use the old positions. 
  !! APJ: For variables that move with v, we need an average of new and old
  do ithermo = istart3, iend3
     thermo => thermostats( ithermo )
     langevin = thermo%langevin
     adaptive = thermo%adaptive
     v   = thermo%v_old( 1 )  ! APJ: correct is thermo%v_old( 1 )
     eta = 0.5d0 * ( thermo%eta_old( 1 ) + thermo%eta( 1 ) )
     if(abfqmmm_param%abfqmmm /= 1) then
       E = E + 0.5d0 * thermo%Q(1)*v**2 + thermo%Ndof_kT * eta
     else
       if(abfqmmm_param%id(thermo%iatom) <= 4) then
         E = E + 0.5d0 * thermo%Q_qm(1)*v**2 + thermo%Ndof_kT * eta
       else
         E = E + 0.5d0 * thermo%Q(1)*v**2 + thermo%Ndof_kT * eta
       end if
     end if
     do j = 3, nchain,2
        v   = thermo%v_old( j ) ! APJ: correct is thermo%v_old( j )
        eta = 0.5d0 * ( thermo%eta_old( j ) + thermo%eta( j ) )
      if(abfqmmm_param%abfqmmm /= 1) then
        E = E + 0.5d0 * thermo%Q(j)*v**2 + thermo%kT * eta
      else
        if(abfqmmm_param%id(thermo%iatom) <= 4) then
          E = E + 0.5d0 * thermo%Q_qm(j)*v**2 + thermo%kT * eta
        else
          E = E + 0.5d0 * thermo%Q(j)*v**2 + thermo%kT * eta
        end if
      end if
     enddo
     do j = 2, nchain,2
        v   = 0.5d0 * ( thermo%v_old( j ) + thermo%v( j ) )
        eta = thermo%eta_old( j ) ! APJ: correct is thermo%eta_old( j )
      if(abfqmmm_param%abfqmmm /= 1) then
        E = E + 0.5d0 * thermo%Q(j)*v**2 + thermo%kT * eta            ! APJ 
      else
        if(abfqmmm_param%id(thermo%iatom) <= 4) then
          E = E + 0.5d0 * thermo%Q_qm(j)*v**2 + thermo%kT * eta
        else
          E = E + 0.5d0 * thermo%Q(j)*v**2 + thermo%kT * eta
        end if
      end if
     enddo
     if (adaptive) then
       j=nchain+1
       v   = thermo%v_old(j) ! APJ: correct is thermo%v_old( 1 )
       eta = thermo%eta_old(j)
       if(abfqmmm_param%abfqmmm /= 1) then
         E = E + 0.5d0 * thermo%Q(j)*v**2 + thermo%Ndof_kT * eta
       else
         if(abfqmmm_param%id(thermo%iatom) <= 4) then
           E = E + 0.5d0 * thermo%Q_qm(j)*v**2 + thermo%Ndof_kT * eta
         else
           E = E + 0.5d0 * thermo%Q(j)*v**2 + thermo%Ndof_kT * eta
         end if
       end if
     endif
     if (langevin) E = E - thermo%E_ln
  enddo

  

#ifdef MPI
  call mpi_allreduce(E,Etmp,1,MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
  E=Etmp
#endif

  if(ntp>0) then
     do j=1,nchain,2
        v   = thermo_lnv%v_old(j)    ! APJ: should be thermo_lnv%v_old(j)
        eta = 0.5d0 * ( thermo_lnv%eta_old(j) + thermo_lnv%eta_old(j) )
        E=E+0.5*thermo_lnv%Q(j)*v**2 + thermo_lnv%kT*eta
     end do
     do j=2,nchain,2
        v   = 0.5d0 * ( thermo_lnv%v_old(j) + thermo_lnv%v(j) )
        eta = thermo_lnv%eta_old(j)   ! APJ: should be thermo_lnv%eta_old(j)
        E=E+0.5*thermo_lnv%Q(j)*v**2 + thermo_lnv%kT*eta
     end do
  end if                      

                                                          
end function

end module
