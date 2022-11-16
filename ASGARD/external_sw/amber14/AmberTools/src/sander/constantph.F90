! <compile=optimized>
#include "copyright.h"
#include "../include/dprec.fh"
#include "../include/assert.fh"

! This module handles all of the constant pH MD capabilities in Amber.
module constantph

use random

! Private by default
private

#  include "dynph.h"

! Public variable declaration
_REAL_, public  :: chrgdat(0:ATOM_CHRG_C-1)
_REAL_, public  :: cph_intdiel
integer, public :: cph_igb
integer, public :: cphfirst_sol
integer, public :: relaxations
integer, public :: tot_relax
integer, public :: mccycles
logical, public :: write_dump

! Public for pH REM
_REAL_, public  :: target_ph

! Private variable declaration
type(const_ph_info)  :: stateinf(0:TITR_RES_C-1)
type(rand_gen_state), save :: cph_gen

character(len=40)   :: resname(0:TITR_RES_C)

_REAL_  :: statene(0:TITR_STATES_C-1)
integer :: hindex(0:TITR_RES_C*(MAX_H_COUNT+1)-1)
integer :: tpair(0:TITR_RES_C*4)
integer :: resstate(0:TITR_RES_C-1)
integer :: iselres(0:TITR_RES_C)
integer :: iselstat(0:2)
integer :: protcnt(0:TITR_STATES_C-1)
integer :: trescnt
logical :: cph_success
logical :: first

integer, allocatable :: mobile_atoms(:)
_REAL_, allocatable  :: frac_pop(:,:)

!  Variable descriptions
!
!  cph_gen        : random number generator for constant pH MC
!  stateinf       : information about titratable residues. See dynph.h
!  chrgdat        : partial atomic charges for every state of every titr. res.
!  statene        : relative energy of every state for each titr. res.
!  hindex         : location of each titratable hydrogen; see below
!  tpair          : neighborlist for interacting titr. residues; see below
!  resstate       : list of current protonation states for each titr. res.
!  iselres        : array of selected, titrating residues
!  iselstat       : array of selected, proposed, states for selected titr. res.
!  protcnt        : number of protons for each state of each titr. res.
!  trescnt        : number of titrating residues
!  cphfirst_sol   : atom number of the first solvent residue
!  cph_igb        : igb model to use for MC evals. in explicit solvent pH MD
!  cph_success    : did the proton exchange attempt succeed or not
!  cph_intdiel    : Internal dielectric to use for MC moves in explicit solvent
!  first          : is this the first pass through? (for CPOUT printing)
!  relaxations    : how many times we have run water relaxations (for timing)
!  tot_relax      : how many relaxations we've done in total (for timing)
!  mccycles       : The number of MC steps that will be attempted at once in
!                   explicit solvent CpHMD
!  write_dump     : Determines if we write a full account of the titration
!                   behavior for explicit solvent CpHMD
!  frac_pop       : Track the populations of each state for each residue when
!                   statistics for a cycle of MC attempts in explicit solvent

! tpair and hindex arrays: 
!  These arrays contain index information for each titr. res., and have similar
!  structures. The head of each array is a pointer section, from indices 0 to 
!  trescnt.  The residue number (0 to trescnt-1) is an index that shows the 
!  location of the FIRST element belonging to that residue (be it neighbors or
!  titrating protons). You can get the *number* of elements belonging to residue
!  "x" by subtracting array(x+1) - array(x). The trescnt-th element of each
!  array is a tail index so that trick also works for the last titr. res.
!  Let me show you an example using hindex: Consider 2 titr. res.  The 0th res
!  has titrating proton numbers 12, 15, 16, and 17. The 1st res has titrating
!  proton numbers 22 and 23. The hindex array would then be arranged as follows:
!  hindex = 3, 7, 9, 12, 15, 16, 17, 22, 23

! Make necessary subroutines public
public cnstphinit, cnstphupdatepairs, cnstphbeginstep, cnstphendstep,  &
       cnstphwriterestart, cnstphwrite, cnstph_explicitmd, cnstphread, &
       total_protonation, cnstph_finalize, &
#ifdef MPI
       cnstph_bcast, &
#endif
       cnstph_zero

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Zero all of the constant pH arrays
subroutine cnstph_zero()
   use constants, only : ZERO, ONE
   implicit none

   ! I zero out all of the arrays at first because when the cprestrt is written,
   ! it simply writes the whole namelist, but unused values would not otherwise
   ! be initialized.

   ! NULL_CPH_INFO is a zero-filled const_ph_info type. See dynph.h

   stateinf(0:TITR_RES_C-1)                  = NULL_CPH_INFO
   chrgdat(0:ATOM_CHRG_C-1)                  = 0
   hindex(0:TITR_RES_C*(MAX_H_COUNT+1)-1)    = 0
   tpair(0:TITR_RES_C*4-1)                   = 0
   resstate(0:TITR_RES_C-1)                  = 0
   iselres(0:2)                              = 0
   iselstat(0:2)                             = 0
   statene(0:TITR_STATES_C-1)                = ZERO
   protcnt(0:TITR_STATES_C-1)                = 0
   resname(0:TITR_RES_C)                     = ' '
   relaxations                               = 0
   tot_relax                                 = 0
   cph_igb                                   = 0
   cph_intdiel                               = ONE
   cphfirst_sol                              = 0

end subroutine cnstph_zero

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read the CPIN file and initialize arrays
subroutine cnstphread(charge)
   use constants, only : AMBER_ELECTROSTATIC, ZERO, ONE
   use file_io_dat, only : CNSTPH_UNIT, cpin, CPH_DUMP_UNIT, cph_dump, owrite
   use memory_module, only : natom
   implicit none
#  include "extra.h"
#  include "../include/md.h"

   ! Variable Descriptions:
   !
   ! Passed:
   !  charge         : partial charge array
   !
   ! Internal:
   !  itres          : residue loop pointer
   !  res            : residue pointer
   !  iatom          : atom pointer
   !  icumstat       : sum of states
   !  icumchrg       : sum of charge
   !
   ! Common memory:
   !  gbgamma        : (md.h) OBC GB parameter
   !  gbbeta         : (md.h) OBC GB parameter
   !  gbalpha        : (md.h) OBC GB parameter

   _REAL_, intent(inout) :: charge(*)

   integer  :: itres, res, iatom, icumstat, icumchrg, ier, i

   ! This subroutine reads in and initializes the constant pH data. It 
   ! multiplies in the AMBER_ELECTROSTATIC factor and adjusts the initial charge
   ! array for the initial protonation states.

   namelist /cnstph/ stateinf, resstate, protcnt, chrgdat, statene, &
                     trescnt, resname, cphfirst_sol, cph_igb, cph_intdiel

   icumstat = -1
   icumchrg = 0
   res = 0

   write(6,'(a,a)') '|reading charge increments from file: ', cpin
   call amopen(CNSTPH_UNIT, cpin, 'O', 'F', 'R')
   read(CNSTPH_UNIT, nml=cnstph)
   do iatom = 0, ATOM_CHRG_C-1
      chrgdat(iatom) = chrgdat(iatom) * AMBER_ELECTROSTATIC
   end do
   close(CNSTPH_UNIT)

   ! Set initial charges
   do itres = 0, trescnt - 1
      do iatom = 0, stateinf(itres)%num_atoms - 1
         charge(iatom + stateinf(itres)%first_atom) &
             = chrgdat(stateinf(itres)%first_charge + iatom + resstate(itres)* &
               stateinf(itres)%num_atoms)
      end do
   end do

   ! Set proper GB stuff
   if ( cph_igb .eq. 2 ) then
      gbgamma = 2.90912499999d0
      gbbeta  = ZERO
      gbalpha = 0.8d0
   else if ( cph_igb .eq. 5 ) then
      gbgamma = 4.85d0
      gbbeta  = 0.8d0
      gbalpha = ONE
   end if

   ! Overflow checking
   if (trescnt > TITR_RES_C) then
      write(6,*) 'Too many titrating residues; alter dynph.h and recompile'
      call mexit(6,1)
   end if
   do itres = 0, trescnt - 1
      if (stateinf(itres)%first_state > icumstat) then
         icumstat = stateinf(itres)%first_state
         res = itres
      end if
   end do

   icumstat = stateinf(res)%first_state + stateinf(res)%num_states
   icumchrg = stateinf(res)%first_charge + stateinf(res)%num_atoms * &
              stateinf(res)%num_states

   if (icumstat > TITR_STATES_C) then
      write(6,*) 'Too many titrating states; alter dynph, recompile'
      call mexit(6,1)
   end if
   if (icumchrg > ATOM_CHRG_C) then
      write(6,*) 'Too much charge data; alter dynph.h, recompile'
      call mexit(6,1)
   end if

   if (icnstph > 1) then
      write(6,'(a,i4,a,i4,a)') '| Attempting ', trescnt, &
               ' MC protonation changes every ', ntcnstph, ' steps.'
      allocate(mobile_atoms(natom), stat=ier)
      REQUIRE (ier == 0)
      mobile_atoms(1:cphfirst_sol-1) = 0
      mobile_atoms(cphfirst_sol:natom) = 1
      ! Only dump to separate file if we have more than 1 MC cycle
      if (mccycles == 1 .or. .not. master) write_dump = .false.
      if (write_dump) then
         call setup_pop_arry ! only master needs this
         call amopen(CPH_DUMP_UNIT, cph_dump, owrite, 'F', 'W')
         write(CPH_DUMP_UNIT, '("Residue Number",6x)', advance='NO')
         do i = 0, max_states() - 1
            write(CPH_DUMP_UNIT, '("State",x,i2,2x)', advance='NO') i
         end do
         write(CPH_DUMP_UNIT, '()')
         write(CPH_DUMP_UNIT, '()')
      end if
   else
      ! Never write all data in implicit solvent (all data is already written to
      ! cpout file)
      write_dump = .false.
   end if

   target_ph = solvph

   ! Set intdiel here. It's not used anywhere except GB, so it doesn't matter if
   ! we set it here

   intdiel = cph_intdiel

end subroutine cnstphread
              
#ifdef MPI
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Broadcasts all of the constant pH information to each thread
subroutine cnstph_bcast(ierr)
   use memory_module, only : natom
   implicit none
#include "parallel.h"
#include "../include/md.h"
 include 'mpif.h'
   integer, intent(out) :: ierr

   ! Broadcast everything
   call mpi_bcast(stateinf, TITR_RES_C*STATEINF_FLD_C, MPI_INTEGER, 0, &
                  commsander, ierr)
   call mpi_bcast(trescnt, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(chrgdat, ATOM_CHRG_C, MPI_DOUBLE_PRECISION, 0, &
                  commsander, ierr)
   call mpi_bcast(resstate, TITR_RES_C, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(statene, TITR_STATES_C, MPI_DOUBLE_PRECISION, 0, &
                  commsander, ierr)
   call mpi_bcast(protcnt, TITR_STATES_C, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(resname, 40*(TITR_RES_C+1), MPI_CHARACTER, 0, &
                  commsander, ierr)
   call mpi_bcast(cphfirst_sol, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(cph_igb, 1, MPI_INTEGER, 0, commsander, ierr)
   call mpi_bcast(mccycles, 1, MPI_INTEGER, 0, commsander, ierr)
   ! Now we have to allocate and broadcast our mobile atoms array
   if (icnstph > 1) then
      if (.not. allocated(mobile_atoms)) &
         allocate(mobile_atoms(natom), stat=ierr)
      REQUIRE (ierr == 0)
      call mpi_bcast(mobile_atoms, natom, MPI_INTEGER, 0, commsander, ierr)
   end if

end subroutine cnstph_bcast
#endif /* MPI */

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Set up an array to track populations of ea state of ea residue for mc cycles
subroutine setup_pop_arry

   implicit none

   integer :: ierror
   integer :: maxstates
   
   ! Make this subroutine re-usable---deallocate the array if it's already
   ! allocated
   if (allocated(frac_pop)) deallocate(frac_pop)

   maxstates = max_states()

   ! Index from zero, like the rest of our arrays here
   allocate(frac_pop(0:trescnt-1, 0:maxstates-1), stat=ierror)

   REQUIRE( ierror == 0 )

   return

end subroutine setup_pop_arry

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Find the most states in any titratable residue
integer function max_states()
   
   implicit none

   integer i

   max_states = stateinf(0)%num_states

   do i = 1, trescnt - 1
      if (stateinf(i)%num_states > max_states) &
         max_states = stateinf(i)%num_states
   end do

   return

end function max_states

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Find all titr. protons in each titr. res. for cnstphupdatepairs
subroutine cnstphinit(x, ig)
   implicit none

#ifdef MPI
   include 'mpif.h'
#  include "parallel.h"
#endif

   ! Variable Descriptions:
   !
   ! Passed:
   !  x              : global position array
   !  ig             : random seed -- used for initializing cph_gen
   !
   ! Local:
   !  res, state     : descriptive pointers in hindex initialization loops
   !  atom, i        : (descriptive) pointers in hindex initialization loops
   !  index_pointer  : index for initializing portions of hindex array
   !  used_atoms     : memory for which titr. protons have been found
   !  is_used        : whether the current proton has been found yet or not
   

   _REAL_, intent(in)  :: x(*)
   integer, intent(in) :: ig
   
   integer :: res, state, atom, i, index_pointer
   logical :: is_used
   integer, dimension(MAX_H_COUNT) :: used_atoms
   integer :: seed
#ifdef MPI
   integer :: ierr
#endif

   ! This subroutine will initialize the hindex array according to the format
   ! described above. We must find each titr. proton in each res. exactly once.
   ! We will use the fact that a titr. proton must have a charge of 0 in at
   ! least 1 state to find them. However, it can be 0 in multiple states, so
   ! take care not to double-count. At the end, initialize the tpair array for
   ! the first time. It also initializes the cph_gen random generator

   ! Every thread in parallel MUST have an identical random number stream.
   ! However, setting ig=-1 in the input file forces desynchronization of the
   ! random numbers. For simplicity's sake, we'll just broadcast the master's
   ! "ig" value (but we don't want to clobber it, so copy it here)

   seed = ig
#ifdef MPI
   call mpi_bcast(seed, 1, mpi_integer, 0, commsander, ierr)
#endif
   call amrset_gen(cph_gen, seed)

   ! The first element points past the header section
   index_pointer = trescnt + 1
   
   do res = 0, trescnt - 1
      hindex(res) = index_pointer
      used_atoms(1:MAX_H_COUNT) = 0

      do atom = 0, stateinf(res)%num_atoms - 1
         do state = 0, stateinf(res)%num_states - 1
            if (chrgdat(stateinf(res)%first_charge + &
                state * stateinf(res)%num_atoms + atom) .eq. 0.d0) then

               is_used = .false.
               do i = 1, MAX_H_COUNT
                  if (used_atoms(i) .eq. (stateinf(res)%first_atom + atom)) &
                     is_used = .true.
               end do ! i = 1, MAX_H_COUNT

               if (.not. is_used) then
                  hindex(index_pointer) = stateinf(res)%first_atom + atom
                  index_pointer = index_pointer + 1

                  do i = 1, MAX_H_COUNT
                     if (used_atoms(i) .eq. 0) then
                        used_atoms(i) = stateinf(res)%first_atom + atom
                        exit
                     end if ! used_atoms(i) .eq. 0
                  end do ! i = 1, MAX_H_COUNT
               
               end if ! .not. is_used
            end if ! chrgdat(stateinf(res)%first_charge ...) .eq. 0.d0
         end do ! state = 0, ...
      end do ! atom = 0, ...
   end do ! res = 0, ...

   ! set tail index
   hindex(trescnt) = index_pointer 

   ! this is our first pass through
   first = .true.

   call cnstphupdatepairs(x)

end subroutine cnstphinit

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Sets the neighborlist (tpair) for "coupled" titrating residues
subroutine cnstphupdatepairs(x)
   implicit none

   ! Variable descriptions:
   !
   ! Passed:
   !  x              : global coordinate array
   !
   ! Local:
   !  resi, resj     : residue loop counters
   !  hi, hj         : proton loop counters for residues i, j
   !  tpair_index    : index for tpair array initialization
   !  ti, tj         : pointer for atom i, j in global position arrays
   !  xij, yij, zij  : xj - xi, yj - yi, zj - zi
   !  CUTOFF         : parameter that determines "neighbor" cutoff criteria

   _REAL_, intent(in)   :: x(*)

   integer  :: resi, resj, hi, hj, ti, tj, tpair_index
   _REAL_   :: xij, yij, zij

   _REAL_, parameter :: CUTOFF = 4.0d0

   ! This subroutine sets the tpair array according to the description above.
   ! Two residues are considered "neighbors" if the square of their distance is
   ! less than or equal to CUTOFF defined above

   ! First element points beyond the header
   tpair_index = trescnt + 1

   do resi = 0, trescnt - 1
      tpair(resi) = tpair_index
      do resj = 0, trescnt - 1
         if (resi .eq. resj) cycle
         
         hloop: do hi = 0, hindex(resi + 1) - hindex(resi) - 1
            ti = 3 * hindex(hi + hindex(resi))

            do hj = 0, hindex(resj + 1) - hindex(resj) - 1
               tj = 3 * hindex(hj + hindex(resj))
               xij = x(ti-2) - x(tj-2)
               yij = x(ti-1) - x(tj-1)
               zij = x(ti)   - x(tj)
               if (4.0d0 .gt. xij * xij + yij * yij + zij * zij) then
                  if (tpair_index < TITR_STATES_C * 4) then
                     tpair(tpair_index) = resj
                     tpair_index = tpair_index + 1
                     exit hloop
                  else
                     write(6,*) "Constant pH pair list overflow. &
                                &Increase TITR_RES_C in dynph.h; recompile"
                     call mexit(6,1)
                  end if 
               end if
            end do
         end do hloop
      end do ! resj
   end do ! resi

   tpair(trescnt) = tpair_index

end subroutine cnstphupdatepairs

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Select random new state for random titrating residue, set charges
subroutine cnstphbeginstep(dcharge)
   implicit none

   ! Variable descriptions:
   !
   ! Passed:
   !  dcharge        : charge array for the proposed state
   !
   ! Internal:
   !  neighborcount  : how many found neighbors we have
   !  randval        : random value from random number generator

   _REAL_, intent(inout) :: dcharge(1:*)

   integer  :: neighborcount
   _REAL_   :: randval

   ! This subroutine will randomly select what to titrate. There is a 25% chance
   ! that we will attempt a multi-site move among neighboring residues. This is
   ! done so we can at times observe a proton transfer when that would otherwise
   ! be a very very rare event. The code only currently supports a 2-site move,
   ! but iselres/iselstat and associated loops have been generalized to an
   ! arbitrary number of sites.
   ! First determine if we will do a multisite move. Then randomly select a res
   ! to titrate. If we tried a multisite move, see if our selected site has any
   ! neighbors, then select from among those to titrate also. Update the dcharge
   ! array for the proposed states. If we have no neighbors, don't do multisite.

   call amrand_gen(cph_gen, randval)
   if (randval .lt. 0.25d0) then
      iselres(0) = 2
   else
      iselres(0) = 1
   end if

   ! select a random residue
   call amrand_gen(cph_gen, randval)
   iselres(1) = int((randval * 0.9999999d0)*trescnt)

   ! select a different, random state
   call amrand_gen(cph_gen, randval)
   iselstat(1) = int((randval * 0.9999999d0) * &
                     (stateinf(iselres(1))%num_states-1))
   if (iselstat(1) .ge. resstate(iselres(1))) &
      iselstat(1) = iselstat(1)+1
   call cnstphupdatechrg(dcharge, iselres(1), iselstat(1))

   if (iselres(0) .eq. 2) then
      neighborcount = tpair(iselres(1)+1) - tpair(iselres(1))
      if (neighborcount .gt. 0) then
         call amrand_gen(cph_gen, randval)
         iselres(2) = tpair(tpair(iselres(1)) + int(randval * 0.9999999d0 * &
                                                    neighborcount)            )
         call amrand_gen(cph_gen, randval)
         iselstat(2) = int((randval*0.9999999d0) * &
                           (stateinf(iselres(2))%num_states - 1))
         if (iselstat(2) .ge. resstate(iselres(2))) iselstat(2) = iselstat(2)+1
         call cnstphupdatechrg(dcharge, iselres(2), iselstat(2))
      else
         iselres(0) = 1
      end if
   end if

end subroutine cnstphbeginstep

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Same as cnstphbeginstep, but for explicit solvent, so we can do multi-MC
subroutine excnstphbeginstep(dcharge, done_res, num_done)
   implicit none

   ! Variable descriptions:
   !
   ! Passed:
   !  dcharge        : charge array for the proposed state
   !  done_res       : array for which residues have been done
   !  num_done       : how many residues have been finished
   !
   ! Internal:
   !  i              : loop counters
   !  neighborcount  : how many found neighbors we have
   !  randval        : random value from random number generator
   !  holder_res     : holder variable so we insert selected residue in order

   _REAL_, intent(inout)  :: dcharge(1:*)
   integer, intent(inout) :: done_res(1:TITR_RES_C)
   integer, intent(in)    :: num_done

   integer  :: i, neighborcount, holder_res(1:TITR_RES_C)
   _REAL_   :: randval

   ! This subroutine will randomly select what to titrate, as long as it has. 
   ! not yet been titrated. That is the only difference from cnstphbeginstep

   call amrand_gen(cph_gen, randval)
   if (randval .lt. 0.25d0) then
      iselres(0) = 2
   else
      iselres(0) = 1
   end if

   ! select a random residue
   call amrand_gen(cph_gen, randval)
   iselres(1) = int((randval * 0.9999999d0) * (trescnt-num_done))
   
   if (num_done .eq. 0) then
      done_res(1) = iselres(1)
   else
      ! now move the selected residue off of those that have been selected already
      ! THIS ARRAY SHOULD ALREADY BE SORTED (see below)
      do i = 1, num_done
         if (iselres(1) .ge. done_res(i)) iselres(1) = iselres(1) + 1
      end do

      ! If the residue we chose was the larger than our last number, add it to
      ! the end
      if (iselres(1) .gt. done_res(num_done)) &
         done_res(num_done + 1) = iselres(1)
   end if ! (num_done .eq. 0)

   ! now add this to the list of residues that we've titrated already, but
   ! keep this list in order!
   do i = 1, num_done
      if (iselres(1) .lt. done_res(i)) then
         holder_res(i:num_done) = done_res(i:num_done)
         done_res(i) = iselres(1)
         done_res(i+1:num_done+1) = holder_res(i:num_done)
         exit
      end if
   end do

   ! select a different, random state
   call amrand_gen(cph_gen, randval)
   iselstat(1) = int((randval * 0.9999999d0) * &
                     (stateinf(iselres(1))%num_states-1))
   if (iselstat(1) .ge. resstate(iselres(1))) &
      iselstat(1) = iselstat(1)+1
   call cnstphupdatechrg(dcharge, iselres(1), iselstat(1))

   if (iselres(0) .eq. 2) then
      neighborcount = tpair(iselres(1)+1) - tpair(iselres(1))
      if (neighborcount .gt. 0) then
         call amrand_gen(cph_gen, randval)
         iselres(2) = tpair(tpair(iselres(1)) + int(randval * 0.9999999d0 * &
                                                    neighborcount)            )
         call amrand_gen(cph_gen, randval)
         iselstat(2) = int((randval*0.9999999d0) * &
                           (stateinf(iselres(2))%num_states - 1))
         if (iselstat(2) .ge. resstate(iselres(2))) iselstat(2) = iselstat(2)+1
         call cnstphupdatechrg(dcharge, iselres(2), iselstat(2))
      else
         iselres(0) = 1
      end if
   end if

end subroutine excnstphbeginstep

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Select random new state for random titrating residue, set charges
subroutine cnstphendstep(dcharge, charge, dvdl, temp, solvph)
   use file_io_dat, only : CPH_DUMP_UNIT
   use constants, only : KB, LN_TO_LOG
   implicit none
#include "extra.h"

   ! Variable descriptions
   !
   ! Passed:
   !  dcharge        : charge array for the proposed state
   !  charge         : charge array for the current state
   !  dvdl           : change in energy between the two states
   !  temp           : system temperature
   !  solvph         : solvent pH
   !
   ! Internal:
   !  randval        : random value from random number generator
   !  deltae         : Non-electrostatic factors adjusting dvdl
   !  i              : loop counter
   !  statebase      : the first state of the selected res. in state arrays

   _REAL_, intent(in)    :: temp, solvph
   _REAL_, intent(inout) :: dcharge(1:*), charge(1:*)
   _REAL_, intent(inout) :: dvdl

   _REAL_   :: randval, deltae
   integer  :: i, statebase
   integer  :: orig_states(2) ! 2 is max # of residues changed in multi-site MC

   ! This subroutine adjusts the delta E passed from electrostatic calculations
   ! for non-EEL factors, and uses the final result to evaluate a MC transition.
   ! If successful, it will call cnstphupdatechrg charge to modify the charge
   ! array to be the same as the dcharge array and update the resstate to be
   ! the newly updated residue state. If it fails, it will reset the
   ! dcharge array to be the same as the charge array. Then it will modify the
   ! cph_success variable based on whether or not we succeeded.

   ! Back up original states
   orig_states(1) = iselstat(1)
   if (iselres(0) == 2) orig_states(2) = iselstat(2)

   deltae = 0.d0
   do i = 1, iselres(0)
      statebase = stateinf(iselres(i))%first_state
      !  deltae = E(proposed state) - E(current state)
      deltae = deltae + statene(iselstat(i)+statebase) - &
                        statene(resstate(iselres(i))+statebase)
      !  correct for pH (delta protons * pH * LN_TO_LOG * KB * T)
      deltae = deltae - (protcnt(iselstat(i)+statebase) - &
                         protcnt(resstate(iselres(i))+statebase)) * solvph * &
                         LN_TO_LOG * KB * temp
   end do
   call amrand_gen(cph_gen, randval)

   dvdl = dvdl - deltae

   if ((dvdl .lt. 0.d0) .or. (randval .le. exp(-dvdl/(KB * temp)))) then
      cph_success = .true.
      do i = 1, iselres(0)
         call cnstphupdatechrg(charge, iselres(i), iselstat(i))
         resstate(iselres(i)) = iselstat(i)
      end do
   else
      do i = 1, iselres(0)
         call cnstphupdatechrg(dcharge, iselres(i), resstate(iselres(i)))
      end do
   end if

   return

end subroutine cnstphendstep

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Update charges to reflect new state
subroutine cnstphupdatechrg(charge, res, inewstate)
   implicit none

   ! Variables
   !
   ! Passed:
   !  charge         : charge array to be updated
   !  res            : the titr. res. to change the charges of
   !  inewstate      : the charge state to change to
   !
   ! Internal:
   !  i              : loop counter

   _REAL_, intent(inout) :: charge(*)
   integer, intent(in)   :: res, inewstate

   integer :: i

   do i = 0, stateinf(res)%num_atoms - 1
      charge(i+stateinf(res)%first_atom) &
            = chrgdat(stateinf(res)%first_charge &
              + inewstate * stateinf(res)%num_atoms + i)
   end do

end subroutine cnstphupdatechrg

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write the cprestrt file
subroutine cnstphwriterestart(inchrgdat)
   use constants, only : AMBER_ELECTROSTATIC
   use file_io_dat, only : CNSTPH_UNIT, cprestrt, owrite, CPOUT_UNIT
   implicit none

   ! Variable descriptions
   !
   ! Passed:
   !  inchrgdat      : chrgdat array - sent as a copy to multiply EEL factor
   !
   ! Internal:
   !  chrgdat        : chgrdat array adjusted by AMBER_ELECTROSTATIC
   !  iatom          : atom counter
   !  stat           : file status
   !  first          : is this the first pass through?

   _REAL_, intent(in)   :: inchrgdat(0:ATOM_CHRG_C-1)

   _REAL_            :: chrgdat(0:ATOM_CHRG_C-1)
   integer           :: iatom
   character(len=7)  :: stat
   logical, save     :: first_cprestrt = .true.

   namelist /cnstph/ stateinf, resstate, protcnt, chrgdat, statene, &
                     trescnt, resname, cphfirst_sol, cph_igb, cph_intdiel

   do iatom = 0, ATOM_CHRG_C - 1
      chrgdat(iatom) = inchrgdat(iatom) / AMBER_ELECTROSTATIC
   end do

   if (first_cprestrt) then
      if (owrite == 'N') then
         stat = 'NEW'
      else if (owrite == 'O') then
         stat = 'OLD'
      else if (owrite == 'R') then
         stat = 'REPLACE'
      else if (owrite == 'U') then
         stat = 'UNKNOWN'
      end if
      open(unit=CNSTPH_UNIT, file=cprestrt, status=stat, form='FORMATTED', &
           delim='APOSTROPHE')
      first_cprestrt = .false.
   else
      open(unit=CNSTPH_UNIT, file=cprestrt, status='OLD', form='FORMATTED', &
           delim='APOSTROPHE')
   end if

   write(CNSTPH_UNIT, nml=cnstph)
   close(CNSTPH_UNIT)

   ! flush all cpout data
   call amflsh(CPOUT_UNIT)

end subroutine cnstphwriterestart

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Write data to the cpout file
subroutine cnstphwrite(rem)
   use file_io_dat, only : CPOUT_UNIT, ntwx
   implicit none
#  include "../include/md.h"

   ! Variables
   !
   ! Passed:
   !  rem            : remd method used
   !
   ! Internal:
   !  full           : do we write a full record or not
   !  i              : loop counter
   !  
   ! Common memory:
   !  irespa         : (md.h) step counter for updating slow-varying terms
   !  nstlim         : (md.h) limit of how many steps we run
   !  solvph         : (md.h) solvent pH
   !  ntcnstph       : (md.h) how often MC protonation jumps are attempted
   !  t              : (md.h) time

   integer, intent(in) :: rem

   logical  :: full
   integer  :: i
   integer  :: nstlim_

   nstlim_ = nstlim
   if (numexchg > 0) &
      nstlim_ = nstlim * numexchg

   if (ntwx .gt. 0) then
      full = (first .or. irespa .eq. nstlim_ .or. mod(irespa,ntwx) .eq. 0)
   else
      full = (first .or. irespa .eq. nstlim_)
   end if

   if (full) then
      write(CPOUT_UNIT, '(a,f8.5)') 'Solvent pH: ', solvph
      write(CPOUT_UNIT, '(a,i8)')   'Monte Carlo step size: ', ntcnstph
      write(CPOUT_UNIT, '(a,i8)')   'Time step: ', irespa
      write(CPOUT_UNIT, '(a,f10.3)') 'Time: ', t
      if (rem == 4 .or. rem == -1) then
         do i = 0, trescnt - 1
            write(CPOUT_UNIT, '(a,i4,a,i2,a,f7.3)') &
               'Residue ', i, ' State: ', resstate(i), ' pH: ', solvph
         end do
      else
         do i = 0, trescnt - 1
            write(CPOUT_UNIT, '(a,i4,a,i2)') 'Residue ', i, ' State: ', resstate(i)
         end do
      end if
   else
      if (rem == 4 .or. rem == -1) then
         do i = 1, iselres(0)
            write(CPOUT_UNIT, '(a,i4,a,i2,a,f7.3)') 'Residue ', iselres(i), &
               ' State: ', resstate(iselres(i)), ' pH: ', solvph
         end do
      else
         do i = 1, iselres(0)
            write(CPOUT_UNIT, '(a,i4,a,i2)') 'Residue ', iselres(i), ' State: ', &
                  resstate(iselres(i))
         end do
      end if
   end if
   write(CPOUT_UNIT, '()')
   first = .false.

end subroutine cnstphwrite

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Calculates the total number of "active" protons
subroutine total_protonation(tot_prot)
   
   implicit none

! Passed Variables
   integer, intent(out) :: tot_prot

! Local Variables
   integer              :: i ! counter

   tot_prot = 0

   do i = 0, trescnt - 1
      tot_prot = tot_prot + protcnt(stateinf(i)%first_state + resstate(i))
   end do

end subroutine total_protonation

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Sets up explicit constant pH transition and water relaxation
subroutine cnstph_explicitmd(xx,ix,ih,ipairs,x,winv,amass,f,v,vold, &
                             xr,xc,conp,skip,nsp,tma,erstop,qsetup, &
                             do_list_update,rem)
   use amd_mod, only : iamd
   use file_io_dat, only : CPH_DUMP_UNIT
   use memory_module, only : natom, l15, l190, l96, l97, l98, l99, ibellygp
   use nbips, only: ips
   use state

   implicit none

#include "box.h"
#include "extra.h"
#include "../include/md.h"
#include "nmr.h"

!  Passed Variables
!
!  XX    : Global real array
!  X     : Global coordinate array
!  WINV  : Inverted masses
!  AMASS : Atomic masses
!  F     : Force array
!  V     : Velocity array
!  VOLD  : Old velocity array (for velocity verlet integrator)
!  XR    : Coords relative to COM of molecule
!  XC    : Restraint reference coords
!  CONP  : Bond parameter for SHAKE
!  SKIP  : Skip array for shake with logical value for each atom
!  NSP   : Submolecule index array (?)
!  TMA   : Submolecular weight array (?)
!  ERSTOP: Should we stop in error (?)
!  QSETUP: Not quite sure...
!  DO_LIST_SETUP : logical that indicates if we should rebuild pairlist (?)

!  Local Variables
!
!  holder_cut     : temporarily store the &cntrl CUT value while we do GB
!  holder_ntb     : temporarily store the &cntrl NTB value while we do GB
!  holder_natom   : temporarily store the true number of atoms
!  holder_ips     : temporarily store the ips value
!  holder_iamd    : temporarily store the iamd value
!  last_solute    : last atom of the solute (first atom of solvent - 1)
!  cph_ener       : record for storing energies
!  iselres        : constant pH residue selection (change prot state)
!  iselstat       : constant pH state selection for chosen residues
!  vtemp          : temporary array holder for the velocities
!  vtemp          : holder array for velocities
!  natom3         : natom * 3
!  selres_holder  : residues that have been chosen to titrate

!  Other used variables from common memory
!  
!  igb         : (md.h) -- GB solvent model to use
!  ntb         : (box.h) -- PBC (1 -> V, 2 -> P, 0 -> none)
!  ntp         : (md.h) -- constant pressure control
!  cut         : (box.h) -- value of cutoff
!  master      : (extra.h) -- True if this is PID 0, False otherwise
!  ntt         : (md.h) -- thermostat
!
!  NOTE: All pointer variables into IX and XX are defined in memory.h and
!        described in locmem.f, where they're used during the allocation
!        of the global integer, real, and hollerith arrays.

   ! Passed arguments
   _REAL_  :: xx(*), x(*), winv(*), amass(*), f(*), v(*), vold(*), &
              xr(*), xc(*), conp(*), tma(*)
   logical :: skip(*), erstop, qsetup, do_list_update
   integer :: ix(*), ipairs(*), nsp(*), rem
   character (len=4) :: ih(*)

   ! local variable declarations
   integer          :: holder_ntb, holder_natom, j,  &
                       i, k, natom3, selres_holder(1:TITR_RES_C),&
                       holder_ntp, holder_ips, holder_iamd
   type (state_rec) :: cph_ener
   _REAL_, dimension(natom*3) :: vtemp
   _REAL_           :: holder_cut

   ! This subroutine is the main driver for running constant pH MD in explicit
   ! solvent. The first thing it does is call cnstphbeginstep to set up the MC.
   ! Then it removes the periodic boundary conditions (PBC), selects a GB model,
   ! then calls force with nstep = 0 which will force oncpstep to be .true., so
   ! we'll get dvdl back. This is then passed to cnstphendstep to evaluate the
   ! transition. If it fails, then we restore the PBC, turn off the GB model,
   ! restore the CUToff, and return to the calling routine. If it's successful,
   ! we still return the above variables, but then we also turn on belly and fix
   ! the protein while we call runmd to relax the solvent for ntrelax steps.
   ! Some important things we do are:
   !  o  turn off T-coupling, since ntt=3 doesn't work with belly
   !  o  turn off trajectory/energy printing
   !  o  zero-out the v and vold arrays for the solute
   !  o  turn on belly and assign belly arrays

   ! Local variable assignments
   natom3          = natom * 3
   holder_cut      = cut
   holder_natom    = natom
   holder_ntb      = ntb
   holder_ntp      = ntp
   holder_ips      = ips

   cut   = 9801.0d0    ! 99 * 99
   natom = cphfirst_sol - 1
   ntb   = 0
   igb   = cph_igb

   cph_success = .false.

   if (write_dump) &
      frac_pop(:,:) = 0

   do j = 1, mccycles
      ! Zero-out the holder for the selected residues
      selres_holder(1:trescnt) = 0

      ! Initiate trial moves for protonation states
      do i = 1, trescnt
         call excnstphbeginstep(xx(l190), selres_holder, i-1)

         ! call force to get energies
         call force(xx,ix,ih,ipairs,x,f,cph_ener,cph_ener%vir, &
                 xx(l96),xx(l97),xx(l98),xx(l99),qsetup, &
                 do_list_update,0)

         ! end constant pH stuff
         call cnstphendstep(xx(l190), xx(l15), cph_ener%pot%dvdl, temp0, solvph)

         ! We do _not_ want this force call to count toward our nmropt counter,
         ! so we decrement it here if nmropt is on
         if (nmropt /= 0) call nmrdcp()
      end do ! i = 1, trescnt

      if (write_dump) then
         do i = 0, trescnt - 1
            frac_pop(i, resstate(i)) = frac_pop(i, resstate(i)) + 1
         end do
      end if

   end do ! j = 1, mccycles

   if (write_dump) then
      write(CPH_DUMP_UNIT, '(29("="),x,"pH",f8.3,x,29("="))') solvph
      do i = 0, trescnt - 1

         write(CPH_DUMP_UNIT, '(("Residue "),i3,9x)', advance='NO') i

         do k = 0, stateinf(i)%num_states - 1
            write(CPH_DUMP_UNIT, '(1f8.6,2x)', advance='NO') &
               dble(frac_pop(i, k)) / dble(mccycles)
         end do

         ! now add the newline
         write(CPH_DUMP_UNIT, '()')

      end do
   end if

   ! restore PBC variables
   cut   = holder_cut
   natom = holder_natom
   ntb   = holder_ntb
   igb   = 0
   ntp   = holder_ntp
   ips   = holder_ips

   ! set iselres so that every residue is printed out
   iselres(0) = trescnt
   do i = 1, trescnt
      iselres(i) = i - 1
   end do

   ! write out the results to the cpout file
   if (master) call cnstphwrite(rem)

   ! if no exchange attempts succeeded, just return and continue dynamics
   if ( .not. cph_success ) return

   ! increase the relaxation count
   relaxations = relaxations + 1
   tot_relax = tot_relax + 1

   ! Disable AMD for relaxation dynamics
   holder_iamd = iamd
   iamd = 0

   do i = 1, cphfirst_sol * 3 - 3
      vtemp(i) = v(i)
      v(i)     = 0
   end do

   ! call runmd to relax the waters
   call relaxmd(xx,ix,ih,ipairs,x,winv,amass,f, &
      v,vold,xr,xc,conp,skip,nsp,tma,erstop, qsetup, &
      ntrelax, mobile_atoms, .false.)

   ! restore the original values
   do i = 1, cphfirst_sol * 3 - 3
      v(i) = vtemp(i)
   end do
   
   ! restore iamd
   iamd = holder_iamd

end subroutine cnstph_explicitmd

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Closes all relevant file units and does necessary deallocations
subroutine cnstph_finalize

   use file_io_dat, only : CPH_DUMP_UNIT
   implicit none
#  include "extra.h"

   if (allocated(mobile_atoms)) deallocate(mobile_atoms)
   if (allocated(frac_pop)) deallocate(frac_pop)

   if (write_dump) close(CPH_DUMP_UNIT)

end subroutine cnstph_finalize

end module constantph
