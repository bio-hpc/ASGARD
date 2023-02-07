!<compile=optimized>

!************************************************************************
!                              AMBER                                   **
!                                                                      **
!                        Copyright (c) 2008                            **
!                Regents of the University of California               **
!                       All Rights Reserved.                           **
!                                                                      **
!  This software provided pursuant to a license agreement containing   **
!  restrictions on its disclosure, duplication, and use. This software **
!  contains confidential and proprietary information, and may not be   **
!  extracted or distributed, in whole or in part, for any purpose      **
!  whatsoever, without the express written permission of the authors.  **
!  This notice, and the associated author list, must be attached to    **
!  all copies, or extracts, of this software. Any additional           **
!  restrictions set forth in the license agreement also apply to this  **
!  software.                                                           **
!************************************************************************



!+++++++++++++++++++++++++++++++++++++++++++++
!This module contains an interface and some state for an EMIL free energy calc.
!The meaningful work is done outside pmemd, in the EMIL library code  
!(kept in ../../src/emil at the time of writing).
!+++++++++++++++++++++++++++++++++++++++++++++
module emil_mod
   implicit none

   integer, public, save             :: emil_do_calc

#ifdef EMIL
   !Bring in some state variables
#include "../include/md.h"
#include "ew_cntrl.h"
#include "../include/memory.h"
#ifdef MPI
#include "parallel.h"
#endif

   private

   !state variables for this module
   double precision, public, save    :: emil_lambdamol, emil_lambdarest
   double precision, public, save    :: emil_softcore_scale, emil_dvdl

   integer, public, save             :: emil_mc_done, emil_n_atoms_reverted 
   
   !!flag to indicate if a new forces calc is needed.
   integer, public, save             :: emil_calc_AMBER
   integer, public, save             :: emil_save_pme, emil_save_gb

   !temporary storage for energies 
   double precision, public          :: emil_base_epot, emil_base_evdw, emil_base_eelec 

   !!Stuff for Metropolis-Hastings... currently does nothing
   integer, dimension(:), allocatable, save :: emil_atoms_reverted

   !subroutines
   public  emil_init
   public  emil_step


 contains

subroutine emil_init (atm_cnt, beta, mass, crd, frc, vel, pbc_box)
  
  !import some state variables
  !use mdread_mod,          only : nstlim, using_gb_potential, vdw_cutoff
  !use prmtop_dat_mod,      only : atm_nsp, atm_qterm, gbl_cn1, gbl_cn2, &
  !                                gbl_labres, gbl_res_atms, nres, ntypes, typ_ico, &
  !                                atm_numex, gbl_natex
  !use mol_list_mod,        only : gbl_mol_cnt 
  !use pbc_mod,             only : pbc_box
  
  use mdin_emil_dat_mod,   only :  emil_paramfile, emil_logfile, emil_model_infile, &
                                   emil_model_outfile

  use memory_module, only        : atom_type_index, atom_name, atoms_per_molecule, &
                                   number_excluded_atoms, excluded_atoms_list, charge, &
                                   residue_label, residue_pointer, nonbonded_parm_index

  use parms, only                : cn1, cn2
  use nblist,               only : cutoffnb


  implicit none

  !Arguments
  integer, intent(in)               :: atm_cnt
  double precision, intent(in)      :: beta, pbc_box(3) 
  double precision, intent(in)      :: mass(atm_cnt)
  double precision, intent(inout)   :: crd(3, atm_cnt), frc(3, atm_cnt), vel(3, atm_cnt)

  emil_lambdamol      = 0.0
  emil_lambdarest     = 0.0
  emil_softcore_scale = 0.0

  emil_save_pme       = 0
  emil_save_gb        = 0

 !  if ( mytaskid .eq. 1 ) then
 !   open(unit=99,file='a1',status='old', access='append')
 !   write(99,*) "again ",mytaskid, numtasks, commsander;
 !   close(unit=99); 
 !  else
 !   open(unit=995,file='a0',status='old', access='append')
 !   write(995,*) "again ",mytaskid, numtasks, commsander;
 !   close(unit=995); 
 !  end if



  ! Setup emil - position wells, open outfiles, all that
#ifdef MPI
  call c_emil_init( atm_cnt, nres, nspm, residue_pointer, atoms_per_molecule, mass, atom_name, &
             residue_label, ntypes, cutoffnb, atom_type_index, nonbonded_parm_index, charge, &
             crd, frc, vel, pbc_box, cn1, cn2, &
             emil_lambdamol, emil_lambdarest,  emil_softcore_scale, &
             beta,  mytaskid, numtasks, commsander, nstlim, &
             number_excluded_atoms, excluded_atoms_list, & 
             trim(emil_paramfile)//char(0), &
             trim(emil_logfile)  //char(0), &
             trim(emil_model_infile)//char(0), &
             trim(emil_model_outfile)//char(0), 0  ) !C strings need a null-terminator
#else
  call c_emil_init( atm_cnt, nres, nspm, residue_pointer, atoms_per_molecule, mass, atom_name, &
             residue_label, ntypes, cutoffnb, atom_type_index, nonbonded_parm_index, charge, &
             crd, frc, vel, pbc_box, cn1, cn2, &
             emil_lambdamol, emil_lambdarest,  emil_softcore_scale, &
             beta,  0, 1, nstlim, &
             number_excluded_atoms, excluded_atoms_list, &
             trim(emil_paramfile)//char(0), &
             trim(emil_logfile)  //char(0), &
             trim(emil_model_infile)//char(0) , &
             trim(emil_model_outfile)//char(0), 0 ) !C strings need a null-terminator
#endif

  return !End subroutine emil_init
  end subroutine emil_init

subroutine emil_step (atm_cnt, nstep, beta, crd, frc, vel, gb_pot_ene, pme_pot_ene, pbc_box)
 
  !import derived datatypes
  !use gb_force_mod,        only : gb_pot_ene_rec, null_gb_pot_ene_rec
  !use pme_force_mod,       only : pme_pot_ene_rec, null_pme_pot_ene_rec

  !Import the energy record derived datatype
  use  state,               only : potential_energy_rec, zero_pot_energy

  !use state_info_mod

  !import state
  !use mdin_ctrl_dat_mod,   only : using_gb_potential, using_pme_potential
  
  implicit none

  !Arguments
  integer, intent(in)                  :: atm_cnt
  integer, intent(in)                  :: nstep
  double precision, intent(in)         :: beta 
  double precision, intent(in)         :: pbc_box(3)
  double precision, intent(inout)      :: crd(3, atm_cnt), frc(3, atm_cnt), vel(3, atm_cnt)
  type(potential_energy_rec), intent(inout)   :: gb_pot_ene
  type(potential_energy_rec), intent(inout)   :: pme_pot_ene
  
  !Local variables
  double precision      :: potential_energy

#ifdef MPI
  integer               :: my_first_atom, my_last_atom
  
  my_first_atom = iparpt(mytaskid)+1
  my_last_atom  = iparpt(mytaskid+1)

#endif

  !Where do we find the total energy of the system?
  if ( igb /= 0 )  then
         potential_energy = gb_pot_ene%tot
  else if ( use_pme /= 0 ) then
         potential_energy = pme_pot_ene%tot
  else 
         potential_energy = 0.0
  endif

  !Try some MC moves
  call c_emil_mcattempt( crd, frc, vel, pbc_box, emil_base_epot, emil_mc_done, &
                           beta, emil_base_evdw, emil_base_eelec )

  !determine forces and energies due to the fancy emil restraining potentials
#ifdef MPI      
  call c_emil_forces( crd, frc, pbc_box, emil_dvdl, potential_energy, &
                        0.0, beta, nstep,                    &
                        my_first_atom, my_last_atom,         &
                        emil_calc_AMBER )
#else
  call c_emil_forces( crd, frc, pbc_box, emil_dvdl, potential_energy, &
                        0.0, beta, nstep, &
                        0, atm_cnt,       &
                        emil_calc_AMBER )
#endif

  !Is the AMBER Hamiltonian in play at all?  If not, can save time:
  ! This block should switch OFF whichever of pme or gb is ON
  !if ( emil_calc_AMBER .eq. 0 ) then
  if( 1 .eq. 0 ) then
     if ( use_pme /= 0 ) then
        emil_save_pme       = use_pme
        use_pme             = 0
        
        !reset for tidier logging
        call zero_pot_energy(pme_pot_ene)

     else if ( igb /= 0 ) then
        emil_save_gb        = igb
        igb                 = 0
   
       !reset for tidier logging
        call zero_pot_energy(gb_pot_ene)

     end if
  else 
     !This block should switch ON whichever of them was previously OFF
     if ( emil_save_pme /= 0) then
        use_pme             = emil_save_pme
        emil_save_pme       = 0
     else if ( emil_save_gb /= 0) then
        igb                 = emil_save_gb
        emil_save_gb        = 0
     end if
  end if

  return !End subroutine emil_step
  end subroutine emil_step

#endif 

end module emil_mod
!===================================================================================================
