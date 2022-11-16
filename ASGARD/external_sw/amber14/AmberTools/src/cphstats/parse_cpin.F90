
#define STATEINF_FLD_C 5
#define TITR_RES_C 50
#define TITR_STATES_C 200
#define ATOM_CHRG_C 1000
#define MAX_H_COUNT 4
#define FN_LEN 256

subroutine parse_cpin(trescnt, protcnt, stateinf, resname, cpin_name, ierr)

   implicit none

   ! The stateinf struct
   type :: const_ph_info
      sequence
      integer :: num_states
      integer :: first_atom
      integer :: num_atoms
      integer :: first_state
      integer :: first_charge
   end type const_ph_info

   ! The namelist variables

   integer             :: trescnt
   integer             :: cphfirst_sol
   integer             :: cph_igb
   integer             :: protcnt(0:TITR_STATES_C-1)
   integer             :: resstate(0:TITR_RES_C-1)
   integer             :: ierr

   double precision    :: cph_intdiel
   double precision    :: statene(0:TITR_STATES_C-1)
   double precision    :: chrgdat(0:ATOM_CHRG_C-1)

   character(len=40)   :: resname(0:TITR_RES_C)

   type(const_ph_info) :: stateinf(0:TITR_RES_C-1)
   type(const_ph_info) :: null_cnstph_info = const_ph_info(0,0,0,0,0)

   ! Is our cpin file read yet?

   logical             :: is_read = .false.

   ! File unit

   integer, parameter  :: CPIN_UNIT = 10

   ! The cpin name

   character(len=FN_LEN), intent(in) :: cpin_name

   ! The public functions


   ! We read it as a namelist
   namelist /cnstph/ stateinf, resstate, protcnt, chrgdat, statene, &
                     trescnt, resname, cphfirst_sol, cph_igb, cph_intdiel

   ! Initialize the namelist variables
   trescnt = 0
   resstate(:) = 0
   protcnt(:) = 0
   chrgdat(:) = 0.d0
   statene(:) = 0.d0
   resname(:) = ' '
   stateinf(:) = null_cnstph_info
   cphfirst_sol = 0
   cph_igb = 0
   cph_intdiel = 0.d0

   ierr = 0

   ! Open the unit, bailing on error
   open(unit=CPIN_UNIT, file=cpin_name, status='OLD', iostat=ierr)
   if (ierr .ne. 0) return

   ! Read the namelist, bailing on error
   read(CPIN_UNIT, nml=cnstph, iostat=ierr)
   if (ierr .ne. 0) return

   ! If we got this far, then our file is read
   is_read = .true.

   return

end subroutine parse_cpin
