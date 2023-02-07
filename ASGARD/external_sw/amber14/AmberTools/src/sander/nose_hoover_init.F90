subroutine nose_hoover_init ( mass, v, f)
! initialize Nose-Hoover chain thermostats
  use pimd_vars, only: nbead,NMPIMD,CMD,ipimd
  use constants, only: pi, hbar, kB
  use full_pimd_vars, only : mybeadid
  use nose_hoover_module, only: nose_hoover_module_init,  &
                                Thermostat_init,  &
                                Thermostat_link,  &
                                Thermostat_type
  use nose_hoover_module, only: file_nhc, nchain, thermo, nthermo, tau, tau_qm

  use cmd_vars, only: omega_nmode
  
  use abfqmmm_module, only: abfqmmm_param

  implicit none
# include "../include/memory.h"
# include "../include/md.h"
# include "parallel.h"
  _REAL_, intent(in) :: mass( natom )
  _REAL_ :: v( 3, natom )
  _REAL_ :: f( 3, natom )
  _REAL_ :: kT, beta
  integer :: idim, iatom
  integer :: j
  logical :: langevin,adaptive
  type(Thermostat_type) :: thermo_dummy

  if ( worldrank.eq.0 ) then
     open( file_nhc, file = "NHC.dat" )
  endif

  kT = kB * temp0
  beta = 1.d0 / kT

  call nose_hoover_module_init
  nthermo = 3 * natom

  allocate(thermo(3,natom))

  if (ipimd==NMPIMD.or.ipimd==CMD) then
     if (irest==1) then
        call trans_vel_from_cart_to_nmode(v)
     endif
  endif

  if (ipimd > 0) then
     if (ipimd.eq.NMPIMD.or.ipimd.eq.CMD) then
        tau = 1.d0/omega_nmode(mybeadid)
     else
        tau = hbar*beta/sqrt(dble(nbead))
     endif
  else
     tau = 1.0/gamma_ln*20.455d0
     if(abfqmmm_param%abfqmmm == 1) then
        if(abfqmmm_param%gamma_ln_qm == 0.0d0) abfqmmm_param%gamma_ln_qm = gamma_ln
        tau_qm = 1.0/abfqmmm_param%gamma_ln_qm*20.455d0
     end if
  end if

  if (ipimd>0) then
     do j=1,(mybeadid-1)*natom*3
        call Thermostat_init(nchain,thermo_dummy,1,0,kT,tau,.false.,&
                             .false.,tau_qm=tau_qm,iatom=iatom)
     enddo
  endif

  langevin = .false.
  adaptive = .false.
  if( ntt==3 .or. ntt==5 .or. ntt==6 .or. ntt ==8 ) langevin = .true.
  if( ntt==6 .or. ntt==7 .or. ntt ==8 )             adaptive = .true.
  do iatom = 1, natom
     do idim  = 1, 3
        call Thermostat_init(nchain,thermo(idim,iatom),1,0,kT, &
            tau,langevin,adaptive,tau_qm=tau_qm,iatom=iatom)
        call Thermostat_link(thermo(idim,iatom),mass(iatom), &
                             v(idim,iatom),f(idim,iatom))
     enddo
  enddo

  if (ipimd>0) then
     do j=mybeadid*natom*3+1,nbead*natom*3
        call Thermostat_init(nchain,thermo_dummy,1,0,kT,tau,.false.,.false.,tau_qm=tau_qm,iatom=iatom)
     enddo
  endif

end subroutine nose_hoover_init


