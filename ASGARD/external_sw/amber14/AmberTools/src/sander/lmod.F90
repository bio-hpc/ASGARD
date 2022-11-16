#include "../include/assert.fh"
#include "../include/dprec.fh"


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Interface to Istvan Kolossvary's LMOD libraries.

module lmod_driver

   ! Description:
   ! Module containing the ARPACK reverse communication driver routines
   ! for the three LMOD packages of Istvan Kolossvary: XMIN, LMOD,
   ! and NVIB.
   ! sander.f contains examples of the usage of this module.

   ! The current set of public routines includes:
   !   subroutine read_lmod_namelist( )
   !          Read the lmod namelist.
   !   subroutine run_lmod( xx, ix, ih, ipairs,  &
   !      coordinates, forces, energies )
   !          Reverse communication driver for LMOD.
   !   subroutine run_xmin( xx, ix, ih, ipairs,  &
   !      coordinates, forces, energies )
   !          Reverse communication driver for XMIN.
   !   subroutine write_lmod_namelist( )
   !          Write the lmod namelist.

   ! The current set of public parameters includes:
   !     LMOD_NTMIN_LMOD
   !     LMOD_NTMIN_XMIN
   ! which are labels for some integral values of input option ntmin.

   ! History:
   ! $Id: lmod.f,v 10.8 2010/03/19 18:04:21 mjw Exp $

   ! Code Description:
   !   Language:           Fortran 90.
   !   Software Standards: Internal Amber Standards.

   ! **********************************************************************
   !  Copyright 2003                                                      *
   !                                                                      *
   !   Modified BSD license                                               *
   !                                                                      *
   !   Redistribution and use in source and binary forms, with or without *
   !   modification, are permitted provided that the following conditions *
   !   are met:                                                           *
   !                                                                      *
   !    1.Redistributions of source code must retain the above copyright  *
   !      notice, this list of conditions and the following disclaimer.   *
   !    2.Redistributions in binary form must reproduce the above         *
   !      copyright notice, this list of conditions and the following     *
   !      disclaimer in the documentation and/or other materials provided *
   !      with the distribution.                                          *
   !    3.The name of the author may not be used to endorse or promote    *
   !      products derived from this software without specific prior      *
   !      written permission.                                             *
   !                                                                      *
   !   THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS''                  *
   !   AND ANY EXPRESS OR IMPLIED WARRANTIES,                             *
   !   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED                         *
   !   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR                      *
   !   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT                   *
   !   SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,                         *
   !   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR                       *
   !   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT                          *
   !   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR                     *
   !   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR                        *
   !   BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON                       *
   !   ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,                      *
   !   STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE                    *
   !   OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE                    *
   !   OF THIS SOFTWARE, EVEN IF ADVISED OF THE                           *
   !   POSSIBILITY OF SUCH DAMAGE.                                        *
   !                                                                      *
   !  To report bugs, suggest enhancements, etc., contact:                *
   !    Scott Brozell                                                     *
   !                  send email to sbrozell@chemistry.ohio-state.edu     *
   !                                                                      *
   ! **********************************************************************

   implicit none

   private
   public :: LMOD_NTMIN_LMOD
   public :: LMOD_NTMIN_XMIN
   public :: read_lmod_namelist
   public :: run_lmod
   public :: run_xmin
   public :: write_lmod_namelist


   ! I/O unit for conflib
   integer,      parameter :: CONFLIB_UNIT   = 69
   ! I/O unit for lmod_trajectory
   integer,      parameter :: LMOD_TRAJECTORY_UNIT   = 96
   ! ntmin value to specify LMOD lmod minimization.
   integer,      parameter :: LMOD_NTMIN_LMOD   = 4
   ! ntmin value to specify LMOD xmin minimization.
   integer,      parameter :: LMOD_NTMIN_XMIN   = 3


   ! lmod namelist definitions in alphabetical order

   integer                 :: arnoldi_dimension = 0
   ! The dimension of the ARPACK Arnoldi factorization.
   ! Zero specifies the whole space, ie, three times the number of atoms.
   ! Renaming of Kolossvary's ndim_arnoldi.
   namelist /lmod/ arnoldi_dimension

   integer                 :: conflib_size = 3
   ! The number of conformations to store in conflib.
   ! Renaming of Kolossvary's nconf.
   namelist /lmod/ conflib_size

   _REAL_                  :: energy_window = 0
   ! The energy window for conformation storage; the energy of a stored
   ! structure will be in the interval [global_min, global_min + energy_window].
   namelist /lmod/ energy_window

   integer                 :: explored_low_modes = 3
   ! The number of low frequency vibrational modes used per LMOD iteration.
   ! Renaming of Kolossvary's kmod.
   namelist /lmod/ explored_low_modes

   integer                 :: frequency_eigenvector_recalc = 3
   ! The frequency, in LMOD iterations, of the recalculation of eigenvectors.
   ! Renaming of Kolossvary's eig_recalc.
   namelist /lmod/ frequency_eigenvector_recalc

   integer                 :: frequency_ligand_rotrans = 1
   ! The frequency, in LMOD iterations, of the rotation/translation of
   ! ligand(s).
   ! Renaming of Kolossvary's apply_rigdock.
   namelist /lmod/ frequency_ligand_rotrans

   integer                 :: lbfgs_memory_depth = 3
   ! The depth of the LBFGS memory for XMIN's LBFGS minimization or TNCG
   ! preconditioning.
   ! The value 0 turns off preconditioning in TNCG minimization.
   ! Renaming of Kolossvary's m_lbfgs.
   namelist /lmod/ lbfgs_memory_depth

   _REAL_                  :: lmod_minimize_grms = 0.1
   ! The gradient RMS convergence criterion of structure minimization.
   ! Renaming of Kolossvary's minim_grms.
   namelist /lmod/ lmod_minimize_grms

   _REAL_                  :: lmod_relax_grms    = 1.0
   ! The gradient RMS convergence criterion of structure relaxation.
   namelist /lmod/ lmod_relax_grms

   integer                 :: lmod_restart_frequency = 5
   ! The frequency, in LMOD iterations, of conflib updating and LMOD
   ! restarting with a randomly chosen structure from the pool.
   ! Renaming of Kolossvary's lmod_restart.
   namelist /lmod/ lmod_restart_frequency

   _REAL_                  :: lmod_step_size_max = 5.0
   ! The maximum length of a single LMOD ZIG move.
   namelist /lmod/ lmod_step_size_max

   _REAL_                  :: lmod_step_size_min = 2.0
   ! The minimum length of a single LMOD ZIG move.
   namelist /lmod/ lmod_step_size_min

   ! The verbosity of the internal status output from the LMOD package:
   ! 0 = none, 1 = some details, 2 = more details, 3 = everything
   ! including ARPACK.
   ! Renaming of Kolossvary's print_level.
   integer,      parameter :: MAXIMUM_LMOD_VERBOSITY = 3
   integer,      parameter :: MINIMUM_LMOD_VERBOSITY = 0
   integer                 :: lmod_verbosity = MINIMUM_LMOD_VERBOSITY
   namelist /lmod/ lmod_verbosity

   ! XMIN's finite difference Hv matrix-vector product method.
   ! Renaming of Kolossvary's numdiff.
   character(len=*), parameter :: MVPM_CENTRAL      = 'central'
   integer,          parameter :: MVPM_CENTRAL_CODE = 2
   character(len=*), parameter :: MVPM_FORWARD      = 'forward'
   integer,          parameter :: MVPM_FORWARD_CODE = 1
   !character(len=len(MVPM_FORWARD))::matrix_vector_product_method=MVPM_FORWARD
   character(len=7) :: matrix_vector_product_method=MVPM_FORWARD ! gfortran hack
   integer                     :: mvpm_code = MVPM_FORWARD_CODE
   namelist /lmod/ matrix_vector_product_method

   ! The Monte Carlo method.
   ! Renaming of Kolossvary's mc_option.
   character(len=*), parameter :: MC_METHOD_METROPOLIS        = 'Metropolis'
   integer,          parameter :: MC_METHOD_METROPOLIS_CODE   = 1
   character(len=*), parameter :: MC_METHOD_TOTAL_QUENCH      = 'Total_Quench'
   integer,          parameter :: MC_METHOD_TOTAL_QUENCH_CODE = 2
   character(len=*), parameter :: MC_METHOD_QUICK_QUENCH      = 'Quick_Quench'
   integer,          parameter :: MC_METHOD_QUICK_QUENCH_CODE = 3
   !character(len=len(MC_METHOD_QUICK_QUENCH)) :: Monte_Carlo_method &
   !                                                  = MC_METHOD_METROPOLIS
   character(len=12) :: Monte_Carlo_method = MC_METHOD_METROPOLIS
   integer                :: Monte_Carlo_method_code = MC_METHOD_METROPOLIS_CODE
   namelist /lmod/ Monte_Carlo_method

   integer                 :: number_free_rotrans_modes =  6  ! no frozen atoms
   ! The number of rotational and translational degrees of freedom.
   ! This is related to the number of frozen/tethered atoms in the system:
   ! 0 atoms dof=6, 1 atom dof=3, 2 atoms dof=1, >=3 atoms dof=0.
   ! Renaming of Kolossvary's nrotran_dof.
   namelist /lmod/ number_free_rotrans_modes

   integer                 :: number_ligand_rotrans = 0
   ! The number of rotational and translational motions applied to
   ! the ligand(s).
   ! Renaming of Kolossvary's nof_poses_to_try.
   namelist /lmod/ number_ligand_rotrans

   integer                 :: number_ligands = 0
   ! The number of ligands for flexible docking.
   ! Renaming of Kolossvary's nlig.
   namelist /lmod/ number_ligands

   integer                 :: number_lmod_iterations = 10
   ! The number of LMOD iterations.
   ! Note that LMOD iterations do not correspond to gradient evaluations.
   ! Renaming of Kolossvary's niter.
   namelist /lmod/ number_lmod_iterations

   integer                 :: number_lmod_moves = 0
   ! The number of LMOD ZIG-ZAG moves.
   ! Zero specifies an unlimited number of ZIG-ZAG moves.
   ! Renaming of Kolossvary's nof_lmod_steps.
   namelist /lmod/ number_lmod_moves

   integer                 :: random_seed = 314159
   ! The seed of the random number generator.
   namelist /lmod/ random_seed

   integer                 :: restart_pool_size = 3
   ! The size of the pool of lowest-energy structures to be used for restarting.
   ! Renaming of Kolossvary's n_best_struct.
   namelist /lmod/ restart_pool_size

   _REAL_                  :: rtemperature = 1.5
   ! The value of RT in Amber energy units.
   ! rtemperature is utilized in the Metropolis criterion.
   ! Renaming of Kolossvary's rtemp.
   namelist /lmod/ rtemperature

   integer,      parameter :: DEFAULT_TOTAL_LOW_MODES = 10
   integer                 :: total_low_modes = DEFAULT_TOTAL_LOW_MODES
   ! The total number of low frequency vibrational modes used.
   ! Renaming of Kolossvary's nmod.
   namelist /lmod/ total_low_modes

   ! XMIN minimization method.
   character(len=*), parameter :: XMIN_METHOD_LBFGS      = 'LBFGS'
   integer,          parameter :: XMIN_METHOD_LBFGS_CODE = 2
   character(len=*), parameter :: XMIN_METHOD_PRCG       = 'PRCG'
   integer,          parameter :: XMIN_METHOD_PRCG_CODE  = 1
   character(len=*), parameter :: XMIN_METHOD_TNCG       = 'TNCG'
   integer,          parameter :: XMIN_METHOD_TNCG_CODE  = 3
   !character(len=len(XMIN_METHOD_LBFGS)) :: xmin_method = XMIN_METHOD_LBFGS
   character(len=5) :: xmin_method = XMIN_METHOD_LBFGS
   integer                     :: xmin_method_code = XMIN_METHOD_LBFGS_CODE
   namelist /lmod/ xmin_method

   ! Verbosity of the internal status output from the XMIN package:
   ! 0 = none, 1 = minimization details, 2 = minimization and
   ! line search details plus CG details in TNCG.
   integer,      parameter :: MAXIMUM_XMIN_VERBOSITY = 2
   integer,      parameter :: MINIMUM_XMIN_VERBOSITY = 0
   integer                 :: xmin_verbosity = MINIMUM_XMIN_VERBOSITY
   namelist /lmod/ xmin_verbosity

   character(len=80) :: conflib_filename = 'conflib'
   namelist /lmod/ conflib_filename

   character(len=80) :: lmod_trajectory_filename = 'lmod_trajectory'
   namelist /lmod/ lmod_trajectory_filename

   character(len=80) :: lmod_job_title = 'job_title_goes_here'
   namelist /lmod/ lmod_job_title


   !Array inputs for (multiple) ligand docking run
   !arrays must be in format [value1,value2,...,valuen] where n=number_ligands
   !ligstart, ligend, ligcent are integer arrays so values must be integer
   !rotmax, rotmin, trmax, trmin are arrays of real, values can be any type of floating point format
   !eg. 1.1 or 1e+3 or 1e-3
   character(len=80) :: ligcent_list = ''
   namelist /lmod/ ligcent_list
   character(len=80) :: ligend_list = ''
   namelist /lmod/ ligend_list
   character(len=80) :: ligstart_list = ''
   namelist /lmod/ ligstart_list
   character(len=80) :: rotmax_list = ''
   namelist /lmod/ rotmax_list
   character(len=80) :: rotmin_list = ''
   namelist /lmod/ rotmin_list
   character(len=80) :: trmax_list = ''
   namelist /lmod/ trmax_list
   character(len=80) :: trmin_list = ''
   namelist /lmod/ trmin_list

contains
! public routines in alphabetical order


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read lmod namelist, validate input, translate strings to codes.

subroutine read_lmod_namelist( )

   implicit none
#include "../include/memory.h"

   integer :: ios
   read( 5, nml=lmod, iostat=ios )
   if ( ios > 0 ) then
      write(6, '(1x,a,i5/)') 'Error reading lmod namelist, iostat= ', ios
      call mexit(6,1)
   end if

   if ( lmod_verbosity > MAXIMUM_LMOD_VERBOSITY ) then
      lmod_verbosity = MAXIMUM_LMOD_VERBOSITY
   else if ( lmod_verbosity < MINIMUM_LMOD_VERBOSITY ) then
      lmod_verbosity = MINIMUM_LMOD_VERBOSITY
   end if

   ! lowercase( matrix_vector_product_method )
   select case ( matrix_vector_product_method )
   ! XMIN's finite difference Hv matrix-vector product method:
   case ( MVPM_CENTRAL )
      mvpm_code = MVPM_CENTRAL_CODE
   case ( MVPM_FORWARD )
      mvpm_code = MVPM_FORWARD_CODE
   case default
      ! invalid matrix_vector_product_method
      write(6,'(/2x,a,a,a)') 'Error: Invalid MATRIX_VECTOR_PRODUCT_METHOD (' &
            , matrix_vector_product_method, ').'
      call mexit(6,1)
   end select

   ! lowercase( Monte_Carlo_method )
   select case ( Monte_Carlo_method )
   ! LMOD's Monte Carlo method:
   case ( MC_METHOD_METROPOLIS )
      Monte_Carlo_method_code = MC_METHOD_METROPOLIS_CODE
   case ( MC_METHOD_TOTAL_QUENCH )
      Monte_Carlo_method_code = MC_METHOD_TOTAL_QUENCH_CODE
   case ( MC_METHOD_QUICK_QUENCH )
      Monte_Carlo_method_code = MC_METHOD_QUICK_QUENCH_CODE
   case default
      ! invalid Monte_Carlo_method
      write(6,'(/2x,a,a,a)') 'Error: Invalid MONTE_CARLO_METHOD (', &
            Monte_Carlo_method, ').'
      call mexit(6,1)
   end select

   ! lowercase( xmin_method )
   select case ( xmin_method )
   ! XMIN minimization method:
   case ( XMIN_METHOD_LBFGS )
      xmin_method_code = XMIN_METHOD_LBFGS_CODE
   case ( XMIN_METHOD_PRCG )
      xmin_method_code = XMIN_METHOD_PRCG_CODE
   case ( XMIN_METHOD_TNCG )
      xmin_method_code = XMIN_METHOD_TNCG_CODE
   case default
      ! invalid xmin_method
      write(6,'(/2x,a,a,a)') 'Error: Invalid XMIN_METHOD (', xmin_method, ').'
      call mexit(6,1)
   end select

   if ( xmin_verbosity > MAXIMUM_XMIN_VERBOSITY ) then
      xmin_verbosity = MAXIMUM_XMIN_VERBOSITY
   else if ( xmin_verbosity < MINIMUM_XMIN_VERBOSITY ) then
      xmin_verbosity = MINIMUM_XMIN_VERBOSITY
   end if

end subroutine read_lmod_namelist


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver routine for LMOD LMOD minimization.
!-----------------------------------------------------------------------
! First crude implementation by Scott Brozell.
! Based on run_xmin and Istvan Kolossvary's NAB implementation.
!
! LMOD LMOD reference.
! This package can be used with any program that can calculate the energy
! and the gradient for a particular [x, y, z] atomic configuration.
! There is no size limit, but the xyz[], grad[], conflib[],
! lmod_trajectory[], lig_start[], lig_end[], lig_cent[], tr_min[], tr_max[],
! rot_min[], and rot_max[] arrays must be allocated by the calling program.
!
! Input arguments:  Number of atoms, xyz[] and grad[] arrays,
! conflib[] and lmod_trajectory[] arrays, number of LMOD simulation steps,
! number of low-frequency vibrational modes used by LMOD,
! number of conformations to be stored within a user-specified energy window,
! the number of ligands to be docked with their identification lig_start[]
! and lig_end[] as well as rotran parameters tr_min[], tr_max[], rot_min[]
! and rot_max[] to specify how the ligand(s) move around during the
! simulation (see details below).
!
! Output arguments: CPU time spent in the lmod() routine, and possibly
! an error message (see below).
!
! On exit, the low-energy conflib library will be loaded into the conflib[]
! array and the LMOD simulation path (trajectory) is loaded into
! lmod_trajectory[]. The conflib libray is also stored and periodically
! updated in a binary file called "conflib.dat". Each block in conflib.dat
! represents a single conformation sorted by increasing energy
! (glob. min. first). A single block consists of a header followed by
! the x, y, z coordinates of the atoms. The header holds three numbers:
! float (8 bytes) energy of the conformation, float (8 bytes) radius of
! gyration and int (4 bytes) multiplicity, i.e., the number of times that
! particular conformation was found during the LMOD simulation.
! The return value of call_lmod() is the global minimum energy and the glob.
! min. structure is loaded into xyz[].
!
! C declaration.
! float lmod( int niter, int nmod, int kmod, int nrotran_dof, int natm,
!    float xyz, float energy, float grad, int nconf, float energy_window,
!    float conflib, float lmod_trajectory, int eig_recalc, int ndim_arnoldi,
!    int lmod_restart, int n_best_struct, int mc_option, float rtemp,
!    float lmod_step_size_min, float lmod_step_size_max, int nof_lmod_steps,
!    int nlig, int lig_start, int lig_end, int apply_rigdock,
!    int nof_poses_to_try, float tr_min, float tr_max, int lig_cent,
!    float rot_min, float rot_max, int random_seed, int print_level,
!    float lmod_time, float aux_time, int return_flag, int status_flag );
!
! Arguments:
!
! apply_rigdock   Frequency by which lmod() applies "nof_poses_to_try"
!                 rigid-body explicit rot/trans to the ligand(s).
!                 apply_rigdock=1 means that such rot/trans takes place at
!                 every LMOD iteration, apply_rigdock=5 means rot/trans is
!                 applied at every five iterations, etc.
! aux_time        Total time spent in the calling program for calculating
!                 energies/gradients and do energy minimizations when
!                 requested by lmod().
! conflib         Allocated array of (x, y, z) atomic coordinates of conflib
!                 structures (dimension = nconf).
! eig_recalc      Frequency by which lmod() recalculates the low-mode
!                 eigenvectors. eig_recalc=1 means that eigenvectors are
!                 recalculated at every LMOD iteration, eig_recalc=5 means
!                 they are only recalculated at every five iterations, etc.
! energy          Energy value.
! energy_window   Max energy gap above global minimum for stored structures.
!                 Energy unit defined in calling program.
! grad            Allocated array of the gradient (dimension=3*natm).
! kmod            Number of modes out of nmod to be explored in each LMOD
!                 iteration.
! lig_start       It is assumed that ligand atoms form a consecutive list
! & lig_end       (no gaps). lig_start[i] is the lowest atom number in
!                 ligand 'i' and lig_end[i] is the highest atom number in
!                 ligand 'i'.  lig_start[] and lig_end[] must be allocated
!                 in the calling program (dimension = nlig).
! lig_cent        lig_cent[i] defines the center of rotation of ligand 'i'.
!                 lig_start[i] <= lig_cent[i] <= lig_end[i] specifies a
!                 particular atom as the center of rotation whereas
!                 lig_cent[i] = 0 means that the center of rotation will
!                 be the center of gravity (geometric centroid) of the
!                 ligand. lig_cent[] should also be allocated by the
!                 calling program to nlig dimensions.
! lmod_relax_grms This is the endpoint criterion of the ZIG
!                 relaxation/minimization in terms of gradient RMS.
!                 lmod_relax_grms = 1 is generally a good value for
!                 crossing the closest barrier. For larger moves a less
!                 stringent criterion can be adequate.
! lmod_restart    Frequency by which conflib[] is updated and LMOD
!                 simulation is restarted with a randomly chosen
!                 structure among the n_best_struct lowest energy
!                 structures found so far. A good value for lmod_restart
!                 is 10, which means that conflib[] is updated
!                 (and written to the binary file conflib.dat) after every
!                 10th LMOD iteration and the simulation is restarted
!                 with one of the lowest-energy structures.
!                 If lmod_restart >= niter (number of LMOD iterations),
!                 the simulation will never restart.
! lmod_step_size_min  Minimum lengths of a single LMOD ZIG move specified
!                     in the distance unit used in the calling program.
!                     A generally good value is 2 Angs.
! lmod_step_size_max  Maximum lengths of a single LMOD ZIG move specified
!                     in the distance unit used in the calling program.
!                     A generally good value is 5 Angs.
!                     The actual length of the ZIG move will be chosen
!                     randomly between the min and max values.
! lmod_time       Total time spent in the lmod() routine in CPU sec.
! lmod_trajectory Allocated array of (x, y, z) atomic coordinates of
!                 consecutive structures visited along the LMOD
!                 path/trajectory (dimension = niter+1).
! mc_option       Monte Carlo option. Allowed values:
!                 '1' Metropolis Monte Carlo, '2' "total quenching"
!                 (the LMOD trajectory always proceeds towards the
!                 lowest lying neighbor of a particular energy well
!                 found after exhaustive search along all of the low modes),
!                 and '3' "quick quenching" (the LMOD trajectory proceeds
!                 towards the first neighbor found, which is lower in
!                 energy than the current point on the path, without
!                 exploring the remaining modes).
! minim_grms      RMS gradient convergence criterion for minimization.
!                 Should be <= 0.1.
! natm            Number of atoms.
! nconf           Number of conformations or docking modes/poses to be
!                 stored in conflib[].
! ndim_arnoldi    This is the dimension of the Arnoldi factorization.
!                 Basically, the ARPACK package used for the eigenvector
!                 calculations solves multiple "small" eigenproblems instead
!                 of a single "large" problem, which is the diagonalization
!                 of the 3xnatm by 3xnatm Hessian matrix. This parameter
!                 is the user-specified dimension of the "small" problem.
!                 Allowed range is nmod+1 <= ndim_arnoldi <= 3xnatm.
!                 ndim_arnoldi=0 translates to ndim_arnoldi=3xnatm which
!                 means that the "small" problem and the "large" problem
!                 are identical. This is the preferred/fastest calculation
!                 for small to medium size systems, because ARPACK is
!                 guaranteed to converge in a single iteration.  The ARPACK
!                 calculation scales with 3xnatm x ndim_arnoldi^2 and,
!                 therefore, for larger molecules there is an optimal
!                 ndim_arnoldi << 3xnatm that converges much faster in
!                 multiple iterations (possibly thousands or
!                 tens of thousands of iterations). For proteins,
!                 ndim_arnoldi=1000 is generally a good value.
! niter           Number of LMOD iterations.
! nlig            Number of ligands considered for flexible docking.
! nmod            Total number of low-frequency vibrational modes utilized
!                 by LMOD.
! nof_lmod_steps  Total number of ZIG-ZAG moves. nof_lmod_steps=0 means
!                 that the number of ZIG-ZAG moves is not pre-defined,
!                 instead LMOD will attempt to cross the barrier in as
!                 many ZIG-ZAG moves as it is necessary. The criterion
!                 of crossing a barrier is as follows:
!
!                     IF the following Boolean holds:
!
!                     1.  The current endpoint of the zigzag trajectory
!                         is lower than the bottom of the current energy
!                         well from where the zigzag trajectory starts.
!                             - OR -
!                     2.  The endpoint is at least lower than it was in
!                         the previous zigzag iteration step.
!                              - AND -
!                         The molecule has also moved farther away from
!                         the starting point in terms of all-atom
!                         superposition RMS.
!
!                     THEN
!
!                         The LMOD ZIG-ZAG trajectory has crossed an energy
!                         barrier :-)
!
!                 nof_lmod_steps > 0 means that multiple barriers may be
!                 crossed and LMOD can carry the molecule to a large
!                 distance on the PES without severely distorting the geometry.
! nof_poses_to_try Number of explicit rotations/translations applied to the
!                 ligand(s) after each apply_rigdock-th LMOD iteration.
! nrotran_dof     Number of external/trivial/rotranslational degrees of
!                 freedom. Allowed values depend on the number of
!                 frozen/tethered atoms in the system: 0 atoms dof=6,
!                 1 atom dof=3, 2 atoms dof=1, >=3 atoms dof=0.
! n_best_struct   Number of the lowest-energy structures found so far at
!                 a particular LMOD restart point.  The structure to be
!                 used for the restart will be chosen randomly from this pool.
!                 n_best_struct = 10 is generally a good choice.
!                 n_best_struct = 1 allows the user to explore the
!                 neighborhood of the then current global minimum.
! print_level     Verbosity: 0= none, 1= some details, 2= more details,
!                 3= everything incl. ARPACK.
! random_seed     Random seed used to reproduce previous runs.
!                 random_seed = 0 uses a hardware seed.
! rot_min         Range of random rotation of ligand 'i' between a minimum
! & rot_max       angle rot_min[i] and a maximum angle rot_max[i] about
!                 the origin specified by lig_cent[i]. The angle is given
!                 in degrees.  rot_min[] and rot_max[] also must be allocated
!                 in the calling program (dimension = nlig).
! rtemp           The value of RT at a particular, user-defined temperature,
!                 given in the energy unit used in the calling program.
!                 rtemp is utilized in the Metropolis criterion.
!
!                 The basic tenet of LMOD is climbing energy barriers with
!                 ease. In libLMOD this is done by utilizing the LMOD
!                 ZIG-ZAG algorithm:
!
!                 A single LMOD move inherently involves excessive bond
!                 stretching and bond angle bending in Cartesian space.
!                 Therefore the primarily torsional trajectory drawn by
!                 the low modes of vibration on the PES is severely
!                 contaminated by this naive, linear approximation and,
!                 therefore, the actual Cartesian LMOD trajectory often
!                 misses its target by climbing walls rather than crossing
!                 over into neighboring valleys at not too high altitudes.
!                 The ZIG-ZAG algorithm consists of a series of alternating
!                 short LMOD moves along the low-mode eigenvector (ZIG)
!                 followed by a few steps of minimization (ZAG), which is
!                 expected to relax excessive stretches and bends more
!                 than reversing the torsional move.  Therefore, it is
!                 expected that such a ZIG-ZAG trajectory will eventually
!                 be dominated by concerted torsional movements and will
!                 carry the molecule over the energy barrier in a way
!                 that is not too different from finding a saddle point
!                 and crossing over into the next valley like passing
!                 through a mountain pass.
!
! tr_min          Range of random translation of ligand 'i' between a
! & tr_max        minimum distance tr_min[i] and a maximum distance tr_max[i].
!                 The distance must be specified in the distance units
!                 used in the calling program. tr_min[] and tr_max[] must
!                 be allocated in the calling program (dimension = nlig).
! xyz             Allocated array of (x,y,z) atomic coordinates (dimension=3*natm).

subroutine run_lmod( xx, ix, ih, ipairs, &
      coordinates, forces, energies, qsetup )

   use state
   use file_io_dat
   implicit none

   _REAL_,  intent(inout) :: xx(*)          ! real dynamic memory
   integer, intent(inout) :: ix(*)          ! integer dynamic memory
   character(len=4), intent(inout) :: ih(*) ! hollerith dynamic memory
   integer, intent(inout) :: ipairs(*)      ! nonbond pair list dynamic memory
   _REAL_,  intent(inout) :: coordinates(*)
   _REAL_,  intent(inout) :: forces(*)
   type(state_rec),  intent(inout) :: energies

   ! ------ External functions -----------------
   _REAL_   lmodC
   external lmodC

#include "../include/md.h"
#include "../include/memory.h"
#include "xmin.h"
#include "parallel.h"

   ! ------ local variables --------------------
! Moved to namelist
!   character(len=80) :: conflib_filename = 'conflib'
!   character(len=80) :: lmod_trajectory_filename = 'lmod_trajectory'
   integer :: i
   integer :: ier
   logical :: is_dump_formatted
   logical :: is_error
   logical :: is_lmod_done
   logical :: qsetup
   integer :: n_force_calls
   integer :: xmin_iter
   integer :: amber_xmin_print_level
   _REAL_  :: aux_time
   _REAL_  :: lmod_time
   integer :: maxiter
   _REAL_  :: grms_tol
   _REAL_  :: minimum_energy
   integer :: natm
   _REAL_  :: grms
   !integer :: iter
   integer :: return_flag
   integer :: status_flag

_REAL_ , dimension(:, :, :), allocatable :: conflib
_REAL_ , dimension(:, :, :), allocatable :: lmod_trajectory
integer, dimension(:), allocatable :: lig_cent
integer, dimension(:), allocatable :: lig_end
integer, dimension(:), allocatable :: lig_start
_REAL_ , dimension(:), allocatable :: rot_max
_REAL_ , dimension(:), allocatable :: rot_min
_REAL_ , dimension(:), allocatable :: tr_max
_REAL_ , dimension(:), allocatable :: tr_min

   ! As of 3/1/14, the lmod capability in sander is still broken:
   ! At least let users know of the sad news:

   write(0,*) 'The lmod capability in sander is broken!'
   write(0,*) '(You might be able to use the NAB interface)'
   write(0,*) 'Exiting....'
   write(6,*) 'The lmod capability in sander is broken!'
   write(6,*) '(You might be able to use the NAB interface)'
   write(6,*) 'Exiting....'
   call mexit(6,1)

   ! Zero the state type as done in runmd()
   energies = null_state_rec

   ! input validation
   if ( total_low_modes > 3*natom ) then
      write(6,'(/2x,a,i5,a)') 'Warning: Invalid total_low_modes (' &
            , total_low_modes,').'
      total_low_modes = min( DEFAULT_TOTAL_LOW_MODES, &
                             3*natom - number_free_rotrans_modes )
      write(6,'(/2x,a,i5,a)') '    Redefined as ' &
            , total_low_modes, '.'
   end if

   ! maxcyc is the maximum number of gradient evaluations.
   ! Since lmod iterations do not correspond to gradient evaluations,
   ! we implement our own cycle control and guarantee that lmod never
   ! exceeds its maximum number of iterations.
   maxiter            = maxcyc + 1
   grms_tol           = drms
   natm               = natom

   allocate( conflib(1:3, natm, conflib_size), stat = ier )
   if( ier /= 0 ) then
      write(6,*) 'Failed to allocate memory for conflib:', 3*natm*conflib_size
      call mexit(6,1)
   end if
   allocate( lmod_trajectory(1:3, natm, number_lmod_iterations + 1), stat = ier)
   if( ier /= 0 ) then
      write(6,*) 'Failed to allocate memory for lmod_trajectory:', 3*natm*(number_lmod_iterations + 1)
      call mexit(6,1)
   end if

!  Case of docking
   if( number_ligands > 0 ) then
      allocate( lig_start(number_ligands), stat = ier )
      if( ier /= 0 ) then
        write(6,*) 'Failed to allocate memory for lig_start:', 1*number_ligands
        call mexit(6,1)
      end if
      allocate( lig_end(number_ligands), stat = ier )
      if( ier /= 0 ) then
       write(6,*) 'Failed to allocate memory for lig_end:', 1*number_ligands
       call mexit(6,1)
      end if
      allocate( rot_min(number_ligands), stat = ier )
      if( ier /= 0 ) then
       write(6,*) 'Failed to allocate memory for rot_min:', 1*number_ligands
       call mexit(6,1)
      end if
      allocate( rot_max(number_ligands), stat = ier )
      if( ier /= 0 ) then
       write(6,*) 'Failed to allocate memory for rot_max:', 1*number_ligands
       call mexit(6,1)
      end if
      allocate( tr_min(number_ligands), stat = ier )
      if( ier /= 0 ) then
       write(6,*) 'Failed to allocate memory for tr_min:', 1*number_ligands
       call mexit(6,1)
      end if
      allocate( tr_max(number_ligands), stat = ier )
      if( ier /= 0 ) then
       write(6,*) 'Failed to allocate memory for tr_max:', 1*number_ligands
       call mexit(6,1)
      end if
      allocate( lig_cent(number_ligands), stat = ier )
      if( ier /= 0 ) then
       write(6,*) 'Failed to allocate memory for lig_cent:', 1*number_ligands
       call mexit(6,1)
      end if
      call split_list_int(ligstart_list, number_ligands, lig_start)
      call split_list_int(ligend_list, number_ligands, lig_end)
      call split_list_real(rotmin_list, number_ligands, rot_min)
      call split_list_real(rotmax_list, number_ligands, rot_max)
      call split_list_real(trmin_list, number_ligands, tr_min)
      call split_list_real(trmax_list, number_ligands, tr_max)
      call split_list_int(ligcent_list, number_ligands, lig_cent)

!     Just to see the arrays are OK
!      do i=1, number_ligands
!        print *, ">>ligstart(",i,")",lig_start(i),"<<"
!        print *, ">>ligend  (",i,")",lig_end(i),"<<"
!        print *, ">>rotmin  (",i,")",rot_min(i),"<<"
!        print *, ">>rotmax  (",i,")",rot_max(i),"<<"
!        print *, ">>trmin   (",i,")",tr_min(i),"<<"
!        print *, ">>trmax   (",i,")",tr_max(i),"<<"
!        print *, ">>ligcent (",i,")",lig_cent(i),"<<"
!      enddo
   end if


   status_flag        = 0
   n_force_calls = 0
   is_error = .false.
   is_lmod_done = .false.
   forces(1:3*natom) = 1.d0   !  keep from thinking atoms are frozen
   do while ( .not. is_lmod_done .and. .not. is_error )
      forces(1:3*natom) = -forces(1:3*natom)

      minimum_energy = lmodC( number_lmod_iterations, total_low_modes, &
            explored_low_modes, number_free_rotrans_modes, natm, &
            coordinates, energies%tot, forces, conflib_size, energy_window, &
            conflib, lmod_trajectory, frequency_eigenvector_recalc, &
            arnoldi_dimension, lmod_restart_frequency, restart_pool_size, &
            Monte_Carlo_method_code, rtemperature, &
            lmod_step_size_min, lmod_step_size_max, number_lmod_moves, &
            number_ligands, lig_start, lig_end, frequency_ligand_rotrans, &
            number_ligand_rotrans, tr_min, tr_max, lig_cent, &
            rot_min, rot_max, random_seed, lmod_verbosity, &
            lmod_time, aux_time, return_flag, status_flag )

      forces(1:3*natom) = -forces(1:3*natom)

      select case ( return_flag )
      case ( DONE )
         ! Finished minimization.
         is_lmod_done = .true.
         is_error = status_flag < 0
      case ( UNIMPLEMENTED )
         ! Normal Amber control of NB list updates.
         is_error = status_flag < 0
         if ( .not. is_error ) then
            amber_xmin_print_level = ntpr
            ntpr = 999999
!            write(6,'(/a)') 'ARPACK requested grad calc; ignore'
            call gradient_calc( xx, ix, ih, ipairs, &
                  coordinates, forces, energies, NBL_UPDATE_NORMAL, &
                  n_force_calls,qsetup )
            ntpr = amber_xmin_print_level
         end if
      case ( 1, 2 )
         ! Coerce a NB list update.
         is_error = status_flag < 0
         if ( .not. is_error ) then
            amber_xmin_print_level = ntpr
            ntpr = 999999
!            write(6,'(/a)') 'ARPACK requested grad calc; ignore'
            call gradient_calc( xx, ix, ih, ipairs, &
                  coordinates, forces, energies, NBL_UPDATE_ALWAYS, &
                  n_force_calls,qsetup )
            ntpr = amber_xmin_print_level
         end if
      case ( 3, 4 )
         ! Prevent a NB list update.
         is_error = status_flag < 0
         if ( .not. is_error ) then
            amber_xmin_print_level = ntpr
            ntpr = 999999
!            write(6,'(/a)') 'ARPACK requested grad calc; ignore'
            call gradient_calc( xx, ix, ih, ipairs, &
                  coordinates, forces, energies, NBL_UPDATE_NEVER, &
                  n_force_calls,qsetup )
            ntpr = amber_xmin_print_level
         end if
      case ( MINIMIZE )
         ! Minimize the current structure with xmin.
         is_error = status_flag < 0
         if ( .not. is_error ) then
            write(6,'(a)') 'Start of LMOD requested XMIN minimization'
            grms = drms  ! drms is an argument to run_xmin via common
            drms = lmod_minimize_grms
            amber_xmin_print_level = ntpr  ! ntpr is a common argument to run_xmin
            ntpr = maxcyc  ! turn off run_xmin printing
            xmin_iter = 0
            call run_xmin( xx, ix, ih, ipairs, &
                  coordinates, forces, energies,qsetup, xmin_iter, ntpr )
            drms = grms  ! restore user specified drms
            ntpr = amber_xmin_print_level  ! restore user specified ntpr
            write(6,'(/a)') 'End of LMOD requested XMIN minimization'
         end if
      case ( RELAX )
         ! Relax the current structure with xmin.
         is_error = status_flag < 0
         if ( .not. is_error ) then
            write(6,'(a)') 'Start of LMOD requested XMIN relaxation'
            grms = drms  ! drms is an argument to run_xmin via common
            drms = lmod_relax_grms
            amber_xmin_print_level = ntpr  ! ntpr is a common argument to run_xmin
            ntpr = maxcyc  ! turn off run_xmin printing
            xmin_iter = 0
            call run_xmin( xx, ix, ih, ipairs, &
                  coordinates, forces, energies, qsetup, xmin_iter, ntpr )
            drms = grms  ! restore user specified drms
            ntpr = amber_xmin_print_level  ! restore user specified ntpr
            write(6,'(/a)') 'End of LMOD requested XMIN relaxation'
         end if
      case default
         ! error from LMOD or the return_flag is corrupted.
         is_error = status_flag < 0
         ASSERT( is_error )
      end select

!      ERRONEOUS ASSUMPTION. LMOD CALLS GRADIENT NOT ONLY FOR MINIMIZATION
!                            BUT FOR ARPACK CALCULATIONS AS WELL.
!      if ( n_force_calls >= maxcyc ) then
!         ! The number of lmod iterations is <= to the number of gradient calcs.
!         ! So ensure termination of the minimization
!         ! maxiter = iter - 1  ! this is sometimes insufficient
!         is_lmod_done = .true.
!         write(6,'(//,a)') '  Maximum number of minimization cycles reached.'
!      end if
   end do

   if ( is_error ) then
      write(6,'(a,i4)') '  LMOD ERROR: Status is ', status_flag
      call mexit(6,1)
   end if

   call report_min_results( xmin_iter, grms, coordinates, &
         forces, energies, ih(m04), xx, ix, ih )  ! ih(m04) = atom names

   ! dump conflib and lmod_trajectory
   is_dump_formatted = ioutfm <= 0
   call amopen( CONFLIB_UNIT, conflib_filename, owrite, 'F', 'W' )
   write( CONFLIB_UNIT, '(a80)' ) lmod_job_title
   do i = 1, conflib_size
      call corpac( conflib, (i-1)*3*natm + 1, i*3*natm, CONFLIB_UNIT, &
            is_dump_formatted )
   end do
   call amopen( LMOD_TRAJECTORY_UNIT, lmod_trajectory_filename, owrite, 'F', 'W' )
   write( LMOD_TRAJECTORY_UNIT, '(a80)' ) lmod_job_title
   do i = 1, number_lmod_iterations
      call corpac( lmod_trajectory, (i-1)*3*natm + 1, i*3*natm, &
            LMOD_TRAJECTORY_UNIT, is_dump_formatted )
   end do
   return

end subroutine run_lmod


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write the lmod namelist if the methods are requested.

subroutine write_lmod_namelist( )

   implicit none
#include "../include/md.h"

   write(6,'(/a)') 'LMOD:'
   if ( ntmin /= LMOD_NTMIN_XMIN .and. ntmin /= LMOD_NTMIN_LMOD ) then
      write(6,'(5x,2(a))') &
            'Warning: namelist present but methods unrequested.'
   end if

   ! always emit XMIN options since LMOD uses XMIN
   write(6,'(5x,2(a))') &
         'xmin_method                  = ', xmin_method
   select case ( xmin_method )
   case ( XMIN_METHOD_LBFGS, XMIN_METHOD_TNCG )
      write(6,'(5x,a,i7)') &
         'lbfgs_memory_depth           = ', lbfgs_memory_depth
   case ( XMIN_METHOD_PRCG )
      ! lbfgs_memory_depth is irrelevant
   case default
      ! xmin_method input validation occurred in read_lmod_namelist
      ASSERT( .false. )
   end select
   write(6,'(5x,3(a))') &
         'matrix_vector_product_method = ', matrix_vector_product_method, &
         ' finite difference'
   write(6,'(5x,a,i7)') &
         'xmin_verbosity               = ', xmin_verbosity

   ! conditionally emit LMOD options
   if ( ntmin == LMOD_NTMIN_LMOD ) then
      write(6,'(5x,a,i7)') &
         'arnoldi_dimension            = ', arnoldi_dimension
      write(6,'(5x,a,i7)') &
         'conflib_size                 = ', conflib_size
      write(6,'(5x,a,g11.5)') &
         'energy_window                = ', energy_window
      write(6,'(5x,a,i7)') &
         'explored_low_modes           = ', explored_low_modes
      write(6,'(5x,a,i7)') &
         'frequency_eigenvector_recalc = ', frequency_eigenvector_recalc
      write(6,'(5x,a,i7)') &
         'frequency_ligand_rotrans     = ', frequency_ligand_rotrans
      write(6,'(5x,a,g11.5)') &
         'lmod_minimize_grms           = ', lmod_minimize_grms
      write(6,'(5x,a,g11.5)') &
         'lmod_relax_grms              = ', lmod_relax_grms
      write(6,'(5x,a,i7)') &
         'lmod_restart_frequency       = ', lmod_restart_frequency
      write(6,'(5x,a,g11.5)') &
         'lmod_step_size_max           = ', lmod_step_size_max
      write(6,'(5x,a,g11.5)') &
         'lmod_step_size_min           = ', lmod_step_size_min
      write(6,'(5x,a,i7)') &
         'lmod_verbosity               = ', lmod_verbosity
      write(6,'(5x,2(a))') &
         'Monte_Carlo_method           = ', Monte_Carlo_method
      write(6,'(5x,a,i7)') &
         'number_free_rotrans_modes    = ', number_free_rotrans_modes
      write(6,'(5x,a,i7)') &
         'number_ligand_rotrans        = ', number_ligand_rotrans
      write(6,'(5x,a,i7)') &
         'number_ligands               = ', number_ligands
      write(6,'(5x,a,i7)') &
         'number_lmod_iterations       = ', number_lmod_iterations
      write(6,'(5x,a,i7)') &
         'number_lmod_moves            = ', number_lmod_moves
      write(6,'(5x,a,i7)') &
         'random_seed                  = ', random_seed
      write(6,'(5x,a,i7)') &
         'restart_pool_size            = ', restart_pool_size
      write(6,'(5x,a,g11.5)') &
         'rtemperature                 = ', rtemperature
      write(6,'(5x,a,i7)') &
         'total_low_modes              = ', total_low_modes
      write(6,'(5x,2(a))') &
         'conflib_filename             = ', conflib_filename
      write(6,'(5x,2(a))') &
         'lmod_trajectory_filename     = ', lmod_trajectory_filename
      write(6,'(5x,2(a))') &
         'lmod_job_title               = ', lmod_job_title
      write(6,'(5x,2(a))') &
         'ligstart_list                = ', ligstart_list
      write(6,'(5x,2(a))') &
         'ligend_list                  = ', ligend_list
      write(6,'(5x,2(a))') &
         'ligcent_list                 = ', ligcent_list
      write(6,'(5x,2(a))') &
         'rotmin_list                  = ', rotmin_list
      write(6,'(5x,2(a))') &
         'rotmax_list                  = ', rotmax_list
      write(6,'(5x,2(a))') &
         'trmin_list                   = ', trmin_list
      write(6,'(5x,2(a))') &
         'trmax_list                   = ', trmax_list
   end if

end subroutine write_lmod_namelist


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver routine for LMOD XMIN minimization.
!-----------------------------------------------------------------------
! First crude implementation by Istvan Kolossvary and Scott Brozell.
! Based on Istvan's NAB implementation.
!
!
! LMOD XMIN reference:
! This package can be used with any program that can calculate the energy
! and the gradient for a particular [x, y, z] atomic configuration.
! There is no size limit, but the xyz[] and grad[] arrays must be allocated
! by the calling program.
!
! Input arguments:  Number of atoms, xyz[] and grad[] arrays, minimization
! method, max number of minimization steps, convergence criterion (see below).
!
! Output arguments: Number of minimization steps completed, final energy and
! the RMS of the gradient, the CPU time spent in the xmin() routine, and
! possibly an error message (see below).
! On exit, xmin() loads the minimized coordinates into xyz[] and grad[] will
! be up-to-date, too.
!
! C declaration.
! float xmin( int method, int maxiter, float grms_tol, int natm, int m_lbfgs,
!    int numdiff, float xyz, float energy, float grad, float grms, int iter,
!    float xmin_time, int print_level, int return_flag, int status_flag );
!
! Arguments:
!
! energy      Energy value.
! grad        Allocated array of the gradient.
! grms        RMS of the gradient.
! grms_tol    Convergence criterion in terms of the RMS of the gradient.
! iter        Number of minimization steps completed.
! maxiter     Max number of minimization steps.
! method      Minimization method: 1= PRCG, 2= LBFGS,
!             3= LBFGS-preconditioned TNCG.
! m_lbfgs     Depth of LBFGS memory for LBFGS minimization or TNCG
!             preconditioning.  Suggested value 5. m_lbfgs=0 with
!             TNCG minimization turns off preconditioning.
! natm        Number of atoms.
! numdiff     Method used in finite difference Hv matrix-vector products:
!             1= forward difference, 2= central difference.
! print_level Verbosity: 0= none, 1= minim details, 2= minim and
!             line search details plus CG details in TNCG
! xmin_time   Time spent in the xmin() routine in CPU sec.
! xyz         Allocated array of (x, y, z) atomic coordinates.

subroutine run_xmin( xx, ix, ih, ipairs, &
      coordinates, forces, energies, qsetup, xmin_iter, ntpr2)

   use constants, only : zero
   use state
   use bintraj, only: end_binary_frame
   use file_io_dat
   implicit none

   _REAL_,  intent(inout) :: xx(*)          ! real dynamic memory
   integer, intent(inout) :: ix(*)          ! integer dynamic memory
   character(len=4), intent(inout) :: ih(*) ! hollerith dynamic memory
   integer, intent(inout) :: ipairs(*)      ! nonbond pair list dynamic memory
   _REAL_,  intent(inout) :: coordinates(*)
   _REAL_,  intent(inout) :: forces(*)
   !_REAL_,  intent(inout) :: energies(*)
   type(state_rec),  intent(inout) :: energies
   integer, intent(inout) :: xmin_iter
   ! BPR :: ntpr now comes in from files.h,
   ! BPR :: so "ntpr2" is used to avoid confusion
   integer, intent(in)    :: ntpr2          ! print frequency
   ! ------ External functions -----------------
   _REAL_   xminc,ddot
   external xminc,ddot

#include "extra.h"
#include "../include/md.h"
#include "../include/memory.h"
#include "xmin.h"
#include "parallel.h"

   ! ------ local variables --------------------
   _REAL_  :: grms            = ZERO
   _REAL_  :: grms_tol
   logical :: is_error
   logical :: ixdump, itdump, iprint
   logical :: loutfm
   logical :: qsetup
   integer :: iter
   integer :: maxiter
   _REAL_  :: minimum_energy
   integer :: n_force_calls
   integer :: nrx
   integer :: return_flag
   integer :: status_flag
   _REAL_  :: xmin_time
   integer ::  ls_method, ls_maxiter, ls_iter, xyz_min
   _REAL_  :: ls_maxatmov, beta_armijo, c_armijo, mu_armijo, ftol_wolfe, &
              gtol_wolfe
   integer :: NBL_CASE

   ! Zero the state type as done in runmd()
   energies = null_state_rec


   ! Open the mdinfo file
   if (master) call amopen(7,mdinfo,'U','F','W')

   ! maxcyc is Amber's maximum number of gradient evaluations.
   ! Since xmin iterations do not correspond to gradient evaluations,
   ! we implement our own cycle control and guarantee that xmin never
   ! exceeds its maximum number of iterations.
   maxiter  = maxcyc
   ! drms is Amber's RMS of the energy gradient convergence criterion.
   ! We use XMIN's name since it is more readable.
   grms_tol = drms
   status_flag = 0
   xyz_min = 1
   ls_method = 2
   ls_maxiter = 20
   ls_maxatmov = 0.2
   beta_armijo = 0.5
   c_armijo = 0.4
   mu_armijo = 1.0
   ftol_wolfe = 0.0001
   gtol_wolfe = 0.9

   n_force_calls = 0
   is_error = .false.

   ! Set an iteration counter for printing and saving purposes
   iter = 0

   ! Ben Roberts: Enable writing to a trajectory file if
   ! requested.
   ! Number of atoms to write to the trajectory
   ! If NTWPRT.NE.0, only print the atoms up to this value
   nrx = 3*nrp
   if (ntwprt > 0) nrx = ntwprt*3
   ! Trajectory format
   loutfm = (ioutfm <= 0)

   ! set initial gradient to non-zero value so that xminC() won't think
   !   that everything is frozen:
   forces(1:3*nrp) = 1.d0

   do while ( .not. is_error ) !Will exit loop if done or encounter an error
      ! mjw Feb 2010
      ! Weird bug here, according to the old code, current_energy get its
      ! value from energies(1), which is energies%tot. This does not work,
      ! its need to get its value from energies%pot%tot  ( energies(23) )

      !MC June 2011 current_energy was only used to pass energies%pot%total to
      !xminc() and has been removed

      forces(1:3*nrp) = -forces(1:3*nrp)

      minimum_energy = xminc( xyz_min, xmin_method_code, maxiter, grms_tol, &
            nrp, lbfgs_memory_depth, mvpm_code, &
            coordinates, energies%pot%tot, forces, grms, xmin_iter, xmin_time, &
            xmin_verbosity, ls_method, ls_maxiter, ls_iter, ls_maxatmov, &
            beta_armijo, c_armijo, mu_armijo, ftol_wolfe, gtol_wolfe,  &
            return_flag, status_flag )

      forces(1:3*nrp) = -forces(1:3*nrp)
      is_error = status_flag < 0

      ! Here if we are starting a new step or if we are done;
      ! print progress if requested
      ! AWG: print if
      !      1) starting geometry (step 0, n_force_calls = 1)
      !      2) if xmin made a step (xmin_iter has increased)
      !      3) when we are done (in which case xmin did not increase xmin_iter)
      iprint = .false.
      ixdump = .false.
      itdump = .false.
      if( .not. is_error .and. (n_force_calls==1) ) then 
         iprint = .true.
         if (ntwx > 0) itdump = .true.
      end if
      if( .not. is_error .and. (xmin_iter/=iter+1)) then
         iter = iter + 1
         if ( mod(iter,ntpr2) == 0 ) iprint = .true.
         if ( (ntwr /= 0) .and. (mod(iter,ntwr) == 0) ) ixdump = .true.
         if ( (ntwx /= 0) .and. (mod(iter,ntwx) == 0) .and. (imin /= 5) ) itdump = .true.
      end if
      if( .not. is_error .and. (return_flag==DONE)) then 
         iter = iter + 1
         iprint = .true.
         ixdump = .true.
         if (ntwx > 0) itdump = .true.
      end if
      
      ! Print energies / gradients
      if ( iprint ) then
         call rmsgrd(forces,grms)
         call report_min_progress( iter, grms, forces, energies, &
               ih(m04), xx(l15) )  ! ih(m04) = atom names, xx(l15) = charges
      end if

      if (master) then
         ! Write a restart file if appropriate
         if (ixdump) call minrit(iter,nrp,ntxo,coordinates)
         ! Write to the trajectory if appropriate
         if (itdump) then
            call corpac(coordinates,1,nrx,MDCRD_UNIT,loutfm)
            if (.not. loutfm) call end_binary_frame(MDCRD_UNIT)
         end if
      end if

      select case ( return_flag )
      case ( DONE ) 
         ! Finished minimization.
         exit
      case ( CALCENRG, CALCGRAD, CALCBOTH )
         ! Normal Amber control of NB list updates.
         NBL_CASE=NBL_UPDATE_NORMAL
      case ( CALCENRG_NEWNBL, CALCGRAD_NEWNBL, CALCBOTH_NEWNBL )
         ! Coerce a NB list update.
         NBL_CASE=NBL_UPDATE_ALWAYS
      case ( CALCENRG_OLDNBL, CALCGRAD_OLDNBL, CALCBOTH_OLDNBL )
         ! Prevent a NB list update.
         NBL_CASE=NBL_UPDATE_NEVER 
      case default
         ! error from XMIN or the return_flag is corrupted.
         ASSERT( is_error )
      end select

      if ( .not. is_error ) then
         call gradient_calc( xx, ix, ih, ipairs, coordinates, forces, & 
              energies, NBL_CASE, xmin_iter, qsetup )
         n_force_calls = n_force_calls + 1
     end if
     
  end do
   
   if ( is_error ) then
      write(6,'(a,i4)') '  XMIN ERROR: Status is ', status_flag
      call mexit(6,1)
   end if

   return

end subroutine run_xmin


! private routines in alphabetical order


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Do processing associated with gradient computation for minimization.
!-----------------------------------------------------------------------
! Maintain bookkeeping for calls to force; call force, ie, calculate
! the gradient; this is based on runmin.  Shake is not performed.
! Various types of nonbond pair list update control are possible.
! Emit a minimization progress report.

subroutine gradient_calc( xx, ix, ih, ipairs, coordinates, &
      forces, energies, list_control, xmin_iter, qsetup)

   use state
   use file_io_dat
   implicit none

   _REAL_,  intent(inout) :: xx(*)          ! real dynamic memory
   integer, intent(inout) :: ix(*)          ! integer dynamic memory
   character(len=4), intent(inout) :: ih(*) ! hollerith dynamic memory
   integer, intent(inout) :: ipairs(*)      ! nonbond pair list dynamic memory
   _REAL_,  intent(inout) :: coordinates(*)
   _REAL_,  intent(inout) :: forces(*)
   !_REAL_,  intent(inout) :: energies(*)
   type(state_rec), intent(inout) :: energies
   integer, intent(in)    :: list_control   ! nonbond pair list update control
   integer, intent(in) :: xmin_iter
   logical :: do_list_update=.false.

#   include "../include/md.h"
#   include "../include/memory.h"
#   include "nmr.h"
#   include "xmin.h"

   logical      :: qsetup
   ! Only used by subroutine force when PSANDER is defined.
   ! By default for minimization list updating is controlled the old
   ! way, ie, nbflag is 0 and updating occurs every nsnb steps;
   ! ntnb is the actual variable that is tested by subroutine nonbond_list.

   _REAL_        :: virials(4)

   iprint = 0
   if ( xmin_iter == maxcyc .or. xmin_iter == 1 ) then
      iprint = 1
   end if

   select case ( list_control )
   case ( NBL_UPDATE_ALWAYS )
      ntnb = 1
   case ( NBL_UPDATE_NEVER )
      ntnb = 0
   case ( NBL_UPDATE_NORMAL )
      if ( mod(xmin_iter,nsnb) == 0 ) then
         ntnb = 1
      end if
   case default
      ! invalid list_control
      ASSERT( .false. )
   end select

   !  no shake for xmin

   call force( xx,ix,ih,ipairs,coordinates,forces,energies,virials, &
         xx(l96),xx(l97),xx(l98), xx(l99), qsetup, &
         do_list_update, xmin_iter)

end subroutine gradient_calc

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine split_list_int(numlist,size,target)
! examples of numlists: [1,5,6,8,14]
   implicit none

   character(*) numlist
   integer size, target(*)
   character(80) buffer
   character(1)  symbol
   integer p, i, inplen, res1, ios, j

   i = 1
   j = 1
   res1 = 1
   inplen = index(numlist,']')-1
   do p=1,inplen
   if ( j.le.size) then
      symbol = numlist(p:p)
      if (symbol.ge.'0'.and.symbol.le.'9') then
         buffer(i:i) = symbol
         i = i + 1
      end if
      if (symbol.eq.','.or. p == inplen) then
!            print *, "res1: >>",buffer(1:i-1),"<<"
            read(buffer(1:i-1),*,iostat=ios) res1
!            write(6,'(5x,a,i7)') &
!            'res1                         = ', res1
            i = 1
            target(j) = res1
            j = j + 1
      end if
    endif
    end do

    if ( j.le.size ) then
       write(6,*) 'Too few parameter at one of LMOD parameter array:', size-j
       call mexit(6,1)
    endif

end subroutine split_list_int

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine split_list_real(numlist,size,target)
! examples of numlists: [1.1,5.6,6.7,8.7,14.11,1e-3,1e+2]
   implicit none

   character(*) numlist
   integer size
   _REAL_ target(*)
   character(80) buffer
   character(1)  symbol
   integer p, i, inplen, ios, j
   _REAL_ res1
   logical switch
   i = 1
   j = 1
   res1 = 1
   inplen = index(numlist,']')-1
   do p=1,inplen
   if ( j.le.size) then
      symbol = numlist(p:p)
      switch = ((symbol.ge.'0'.and.symbol.le.'9').or.(symbol.eq.'.'))
      switch = switch.or.(symbol.eq.'E').or.(symbol.eq.'e')
      switch = switch.or.(symbol.eq.'+').or.(symbol.eq.'-')
      if (switch) then
         buffer(i:i) = symbol
         i = i + 1
      end if
      if (symbol.eq.','.or. p == inplen) then
!            print *, "res1: >>",buffer(1:i-1),"<<"
            read(buffer(1:i-1),*,iostat=ios) res1
!            write(6,'(5x,a,i7)') &
!            'res1                         = ', res1
            i = 1
            target(j) = res1
            j = j + 1
      end if
    endif
    end do

    if ( j.le.size ) then
       write(6,*) 'Too few parameter at one of LMOD parameter array:', size-j
       call mexit(6,1)
    endif

end subroutine split_list_real

end module lmod_driver
