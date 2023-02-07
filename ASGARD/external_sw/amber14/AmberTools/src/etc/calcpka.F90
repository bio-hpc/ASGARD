!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!                                                                              :
!      This program will calculate the pKa of each residue based on the        :
!      given cpin and cpout files.                                             :
!                                                                              :
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

program calcpka

   implicit none

!  Variables:

!     STATEINF_FLD_C     : How many fields are in the type const_ph_info
!     TITR_RES_C         : How many residues may be titrated
!     TITR_STATES_C      : How many states may exist
!     ATOM_CHRG_C        : How many charges may exist
!     const_ph_info      : Type for constant pH info
!     max_nstate         : Maximum number of states in titratable residues

!     cpin               : Name of cpin specified on command line (CL)
!     cpout(:)           : Array holding all cpout file names specified on CL
!     output_file        : File where final statistics are printed
!     alternative_output : Output file for dumping of chunks or running sum
!     population_output  : Output file with detailed populations for each state
!     arg_holder         : holds a command-line argument
!     line_holder        : holds a parsed line from the cpout

!     cpin_unit          : Fortran file unit for the cpin file
!     cpout_unit         : Fortran file unit for the cpout file
!     output_unit        : Fortran file unit for the output file
!     alternative_unit   : Fortran file unit for the time-dumped output file
!     population_unit    : Fortran file unit for the population data file
!     ios                : stat for opening file -- catch errors

!     narg               : Number of CL arguments
!     num_cpout          : Number of cpout files to parse
!     dump_interval      : How frequently to dump pKa chunks or cumulative sums.
!                          Literally, how many ntcnstph steps between dumping pKas
!     i, j               : Counters in do loops
!     cpout_done         : Logical to see if we're done finding cpout files on the CL

!     solvph             : Solvent pH in cpout file
!     step_size          : Monte carlo step size

!     stateinf           : Holds info for each titrating residue
!     resstate           : Array for the state of each residue
!     protcnt            : Array for number of protons in each state
!     statene            : State energy array for each state
!     chrgdat            : Array with all partial atomic charges
!     trescnt            : Number of titratable residues
!     resname            : Array holding all residue names
!     cph_igb            : GB model used for hybrid explicit CpHMD
!     cphfirst_sol       : Atom number of first bulk solvent atom
!     avg_prot           : Average protonation count -- running avg
!     avg_prot_chunk     : Average protonation count reset for chunks
!     frames             : counter for how many frames we have
!     frames_chunk       : counter for how many frames since the last output dump

!     protonations       : Array holding populations of every protonation state
!     protonations_chunk : same as above, but reset every interval steps
!     updating           : Logical to see if we are still updating protonations
!     res_holder         : Holder variable for which residue was selected
!     state_holder       : Holder for which state was selected for the above residue
!     onstep             : Which step we are currently on
!     protonated         : array that holds the protcnt for protonated for each residue
!     transitions        : array that holds the number of transitions that each residue has made

   ! From dynph.h
   integer, parameter :: STATEINF_FLD_C = 5
   integer, parameter :: TITR_RES_C     = 50
   integer, parameter :: TITR_STATES_C  = 200
   integer, parameter :: ATOM_CHRG_C    = 1000

   integer :: max_nstate = 0

   type :: const_ph_info
      sequence
      integer :: num_states, first_atom, num_atoms, first_state, first_charge
   end type const_ph_info
   ! end dynph.h

   character (len=256)              :: cpin
   character (len=256), allocatable :: cpout(:)
   character (len=256)              :: output_file = 'none'
   character (len=256)              :: alternative_output = "pKa_evolution.dat"
   character (len=256)              :: population_output = "populations.dat"
   character (len=256)              :: arg_holder
   character (len=80)               :: line_holder

   integer :: narg, iargc
   integer :: num_cpout = 1
   integer :: dump_interval = 0
   integer :: i, j
   integer :: ios
   integer :: frames = 0
   integer :: frames_chunk = 0
   real    :: avg_prot = 0.0d0
   real    :: avg_prot_chunk = 0.0d0
   logical :: cpout_done = .false.
   logical :: rem = .false.

   real    :: solvph
   integer :: step_size

   integer, parameter :: cpin_unit = 10
   integer, parameter :: cpout_unit = 20
   integer, parameter :: output_unit = 6
   integer, parameter :: alternative_unit = 30
   integer, parameter :: population_unit = 40

   type (const_ph_info) :: stateinf(0:TITR_RES_C-1)
   integer              :: resstate(0:TITR_RES_C-1)
   integer              :: protcnt(0:TITR_STATES_C-1)
   real                 :: statene(0:TITR_STATES_C-1)
   real                 :: chrgdat(0:ATOM_CHRG_C-1)
   real                 :: cph_intdiel
   integer              :: trescnt
   character (len=40)   :: resname(0:TITR_RES_C)
   integer              :: cph_igb, cphfirst_sol

   integer, allocatable :: protonations(:,:), protonations_chunk(:,:)
   logical              :: updating
   integer              :: res_holder, state_holder
   integer              :: onstep = 0
   integer, allocatable :: protonated(:), transitions(:)

   namelist /cnstph/     stateinf, resstate, protcnt,  statene,  &
                         chrgdat,  trescnt,  resname,  cph_igb,  &
                         cph_intdiel, cphfirst_sol

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! END VARIABLE DECLARATIONS !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   narg = iargc() ! load number of arguments into narg

   ! Make sure that enough arguments are given
   if (narg .lt. 2) then
      call usage()
      call exit(1)
   end if

   ! get the cpin
   call getarg(1, cpin) 

   ! answer cries for help
   if (cpin .eq. "-h" .or. cpin .eq. "--help") then
      call usage()
      call exit(0)
   end if

   ! get the rest of the CL arguments
   i = 3
   do while (i .le. narg)
      call getarg(i, arg_holder)

      if (arg_holder(1:1) .eq. '-') then
         cpout_done = .true.
      end if

      if (arg_holder .eq. '-o') then
         i = i + 1
         call getarg(i, output_file)
      else if (arg_holder .eq. '-t') then
         i = i + 1
         call getarg(i, arg_holder)
         read (arg_holder,FMT='(I10)') dump_interval
      else if (arg_holder .eq. '-ao') then
         i = i + 1
         call getarg(i, alternative_output)
      else if (arg_holder .eq. '-po') then
         i = i + 1
         call getarg(i, population_output)
      else if (cpout_done) then
         call usage()
         call exit(1)
      else
         num_cpout = num_cpout + 1
      end if

      i = i + 1
   end do

   ! allocate cpout name array and fill it
   allocate(cpout(num_cpout))
   do i = 1, num_cpout
      call getarg(i+1,cpout(i))
   end do

   ! open the cpin file and parse it
   open(unit=cpin_unit, file=cpin, status='OLD', iostat=ios)
   if (ios .ne. 0) then
      write(0,*) 'Error: CPIN file ', trim(cpin), ' cannot be opened!'
      call usage()
      call exit(1)
   end if
   read(cpin_unit, nml=cnstph)
   close(cpin_unit)

   ! open output_file to stdout if it's specified
   if (output_file .ne. 'none') then
      open(unit=output_unit, file=output_file, status='REPLACE')
   end if

   ! open the alternative output file if it's specified
   if (dump_interval .ne. 0) then
      open(unit=alternative_unit, file=alternative_output, status='REPLACE')
   end if

   ! open the population output file
   open(unit=population_unit, file=population_output, status='REPLACE')

   ! open up the first cpout file and get the necessary information
   open(unit=cpout_unit, file=cpout(1), status='OLD',iostat=ios)
   if (ios .ne. 0) then
      write(0,fmt="(a,a,a)") "Error: CPOUT ", trim(cpout(1)), " does not exist!"
      call usage()
      stop 1
   end if

   read(cpout_unit,fmt='(1a80)') line_holder
   if (line_holder(1:11) .ne. "Solvent pH:") then
      write(0,'(a,a,a)') 'Error: ', trim(cpout(1)), ' is an invalid CPOUT file!'
      call usage()
      stop 1
   end if
   read(line_holder(13:),'(f8.5)') solvph

   read(cpout_unit,fmt='(1a80)') line_holder
   read(line_holder(24:),fmt='(i8)') step_size
   close(cpout_unit)

   ! allocate storage for all of the protonations and zero it out
   allocate(transitions(trescnt))
   allocate(protonated(trescnt))
   do i = 1, trescnt
      transitions(i) = 0
      protonated(i) = 0
      do j = 0, stateinf(i-1)%num_states-1
         if (protcnt(stateinf(i-1)%first_state + j) .gt. protonated(i)) then
            protonated(i) = protcnt(stateinf(i-1)%first_state + j)
         end if
      end do
   end do

   ! Find the maximum number of states
   do i = 1, trescnt
      if (stateinf(i-1)%num_states > max_nstate) then
         max_nstate = stateinf(i-1)%num_states
      end if
   end do

   allocate(protonations(trescnt, max_nstate)) 
   allocate(protonations_chunk(trescnt, max_nstate))
   call empty_protonation(protonations, trescnt, max_nstate)
   call empty_protonation(protonations_chunk, trescnt, max_nstate)

   ! Now loop through every cpout file to calculate the pKas.
   i = 1
   do while (i <= num_cpout)

      ! open up the i-th cpout file, but go to the next if it's not valid
      open(unit=cpout_unit, file=cpout(i), status='OLD', iostat=ios)
      if (ios .ne. 0) then
         write(0,'(a,a,a)') 'Error: CPOUT file ', trim(cpout(i)), ' does not exist!'
         call exit(1)
      end if

      ! read the cpout file
      do while(.true.)
         read(unit=cpout_unit,fmt='(1a80)',end=9) line_holder

         ! Now update the protonations
         if (line_holder(1:8) .eq. "Residue ") then
            frames = frames + 1
            frames_chunk = frames_chunk + 1
            updating = .true.
            do while(updating)
               read(line_holder(9:12),'(i4)') res_holder
               read(line_holder(21:22),'(i2)') state_holder
               if (protcnt(stateinf(res_holder)%first_state + resstate(res_holder)) .ne. &
                   protcnt(stateinf(res_holder)%first_state + state_holder)) then
                  transitions(res_holder + 1) = transitions(res_holder + 1) + 1
               end if
               resstate(res_holder) = state_holder
               read(cpout_unit,fmt='(1a80)',end=9) line_holder
               if (line_holder(1:8) .ne. "Residue ") then 
                  updating = .false.
               end if
            end do !while(updating)

            onstep = onstep + step_size

            call average_protonations(avg_prot, resstate, protcnt, stateinf, trescnt)
            call average_protonations(avg_prot_chunk, resstate, protcnt, stateinf, trescnt)
            
            if (dump_interval .eq. 0) then
               do j = 1, trescnt
                  protonations(j,resstate(j-1)+1) = protonations(j, resstate(j-1)+1 ) + 1
               end do ! j = 1, trescnt
            else
               do j = 1, trescnt
                  protonations(j,resstate(j-1)+1) = protonations(j,resstate(j-1)+1 ) + 1
                  protonations_chunk(j,resstate(j-1)+1) = protonations_chunk(j,resstate(j-1)+1 ) + 1
               end do
               if (mod(onstep, dump_interval) .eq. 0) then
                  write(unit=alternative_unit,fmt='(a)') "========================== CUMULATIVE ========================"
                  call calculate_pKas(protonations, solvph, stateinf, protcnt, trescnt, max_nstate, alternative_unit, &
                                      protonated, resname, transitions, avg_prot, frames)
                  write(unit=alternative_unit,fmt='(a)') "============================ CHUNK ==========================="
                  call calculate_pKas(protonations_chunk, solvph, stateinf, protcnt, trescnt, max_nstate, alternative_unit, &
                                      protonated, resname, transitions, avg_prot_chunk, frames_chunk)
                  call empty_protonation(protonations_chunk, trescnt, max_nstate)
                  avg_prot_chunk = 0
                  frames_chunk = 0
               end if ! mod(onstep
            end if ! dump_interval
         end if ! line_holder(1:8)

      end do ! while(.true.)
      

9  close(cpout_unit)
   i = i + 1
   end do

   call calculate_pKas(protonations, solvph, stateinf, protcnt, trescnt, max_nstate, &
                       output_unit, protonated, resname, transitions, avg_prot, frames)

   write(population_unit,*) "Populations: "
   write(population_unit,*)

   call dump_protonations(protonations, stateinf, protcnt, trescnt, max_nstate, &
                          population_unit, resname)

   close(output_unit)
   close(alternative_unit) 
   close(population_unit)

end program calcpka

subroutine usage()

   implicit none

   write(0,*) 'Usage: calcpka <cpin> <cpout1> {<cpout2> ... <cpoutN>} {-o output} \\'
   write(0,*) '               {-t dump_interval} {-ao dump_output} {-po population_output}'

end subroutine usage

subroutine empty_protonation(array, index1, index2)

   integer, intent (in)  :: index1, index2
   integer, intent (out) :: array(index1, index2)
   integer :: i, j

   do i = 1, index1
      do j = 1, index2
         array(i,j) = 0
      end do
   end do

end subroutine empty_protonation

subroutine calculate_pKas(protonations, solvph, stateinf, protcnt, trescnt, &
                          max_nstate, fileno, protonated, resname, transitions, &
                          avg_prot, nsteps)

   type :: const_ph_info ! from dynph.h
      sequence
      integer :: num_states, first_atom, num_atoms, first_state, first_charge
   end type const_ph_info

   integer, intent (in)             :: max_nstate, trescnt, fileno
   integer, intent (in)             :: protonations(1:trescnt, 1:max_nstate)
   type (const_ph_info), intent(in) :: stateinf(0:*)
   integer, intent (in)             :: protcnt(0:*)
   integer, intent(in)              :: protonated(0:*)
   real, intent(in)                 :: solvph
   character (len=40), intent(in)   :: resname(0:*)
   integer, intent(in)              :: transitions(*)
   real, intent(in)                 :: avg_prot
   integer, intent(in)              :: nsteps

   integer                          :: i, j ! counters
   real, dimension(trescnt)         :: pkas, fracprot
   real                             :: numprot, numdep
   character (len=40)               :: rnm

   do i = 0, trescnt-1
      numprot = 0
      numdep = 0
      do j = 0, stateinf(i)%num_states-1
         if (protcnt(stateinf(i)%first_state + j) .eq. protonated(i)) then
            numprot = numprot + protonations(i+1,j+1)
         else
            numdep = numdep + protonations(i+1,j+1)
         end if
      end do
      pkas(i+1) = solvph - log10(numdep / numprot)
      fracprot(i+1) = numprot / (numdep + numprot)
   end do

   write(fileno, '(a,f8.3)') "Solvent pH is ", solvph
   do i = 1, trescnt
      rnm = resname(i)
      write(fileno, '(a,a,f6.3,a,f6.3,a,f5.3,a,i9)') rnm(10:17), ': Offset ', &
               pkas(i) - solvph, "  Pred ", pkas(i), "  Frac Prot ", fracprot(i), "  Transitions ", &
               transitions(i)
   end do

   write(fileno, '()')
   write(fileno, '(a,f7.3)') "Average total molecular protonation: ", avg_prot / nsteps

end subroutine calculate_pKas

subroutine dump_protonations(protonations, stateinf, protcnt, trescnt, max_nstate, &
                             fileno, resname)

   type :: const_ph_info ! from dynph.h
      sequence
      integer :: num_states, first_atom, num_atoms, first_state, first_charge
   end type const_ph_info

   integer, intent (in)             :: max_nstate, trescnt, fileno
   integer, intent (in)             :: protonations(1:trescnt, 1:max_nstate)
   type (const_ph_info), intent(in) :: stateinf(0:*)
   integer, intent (in)             :: protcnt(0:*)
   character (len=40), intent(in)   :: resname(0:*)

   integer                          :: i, j, linelen ! counters
   real                             :: totpts = 0
   character (len=256)              :: line
   character (len=13)               :: holder

   write(line,'(1a17,x)') "Residue Number"
   linelen = 18

   do i = 0, max_nstate - 1 ! write the header
      write(holder, '(1a9,i3,x)') "    State", i

      line = line(1:linelen) // holder
      linelen = linelen + 13
   end do

   write(fileno,'(a)') trim(line) ! print the header to the screen
   line = '------------------'
   do i = 1, max_nstate
      line = trim(line) // '-------------'
   end do
   
   write(fileno,'(a)') trim(line)

   do i = 0, trescnt - 1 ! write out all of the populations
      line = ''
      totpts = 0
      write(line, '(1a17)') resname(i+1)
      do j = 0, stateinf(i)%num_states - 1
         totpts = totpts + protonations(i+1,j+1)
      end do
      linelen = 18
      do j = 0, stateinf(i)%num_states - 1
         write(holder,'(f8.6,1a2,i1,1a2)') protonations(i+1,j+1)/totpts, ' (', &
            protcnt(stateinf(i)%first_state+j), ') '
         line = line(1:linelen) // holder
         linelen = linelen + 13
      end do
      write(fileno,'(a)') trim(line)
   end do

end subroutine dump_protonations

subroutine average_protonations(avg_prot, resstate, protcnt, stateinf, trescnt)
   implicit none

   ! stateinf type
   type :: const_ph_info
      sequence
      integer :: num_states, first_atom, num_atoms, first_state, first_charge
   end type const_ph_info

   ! Variables

   ! Passed:
   !  avg_prot    : average protonations that we're adding to
   !  resstate    : current state of each titratable residue
   !  protcnt     : number of protons in each state of each residue
   !  stateinf    : titratable residue state information
   !  trescnt     : titratable residue count

   ! Local:
   !  i           : loop counter

   real, intent(inout)              :: avg_prot
   integer, intent(in)              :: resstate(0:*)
   integer, intent(in)              :: protcnt(0:*)
   integer, intent(in)              :: trescnt
   type (const_ph_info), intent(in) :: stateinf(0:*)

   integer :: i

   do i = 0, trescnt - 1
      avg_prot = avg_prot + protcnt(stateinf(i)%first_state + resstate(i))
   end do

end subroutine average_protonations
