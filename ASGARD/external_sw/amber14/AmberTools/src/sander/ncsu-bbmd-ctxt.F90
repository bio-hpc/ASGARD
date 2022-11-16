#include "ncsu-utils.h"
#include "ncsu-config.h"


module ncsu_bbmd_ctxt

#ifdef NCSU_ENABLE_BBMD

use ncsu_constants, only : SL => STRING_LENGTH, BBMD_MONITOR_UNIT

use ncsu_umbrella, only : umbrella_t, &
   MAX_NUMBER_OF_COLVARS => UMBRELLA_MAX_NEXTENTS

use ncsu_colvar_type, only : colvar_t

implicit none

private

integer, private, parameter :: MONITOR_UNIT = BBMD_MONITOR_UNIT

character(*), private, parameter :: SECTION = 'ncsu_bbmd'

character(*), private, parameter :: &
   DEFAULT_MONITOR_FILE = 'ncsu-bbmd-monitor', &
   DEFAULT_UMBRELLA_FILE = 'ncsu-bbmd-umbrella', &
   DEFAULT_SNAPSHOTS_BASENAME = 'ncsu-bbmd-umbrella-snapshot'

integer, private, parameter :: MODE_NONE = 5432

integer, private, parameter :: MODE_ANALYSIS = 1234
integer, private, parameter :: MODE_UMBRELLA = 2345
integer, private, parameter :: MODE_FLOODING = 3456

!-------------------------------------------------------------------------------

type, public :: bbmd_ctxt_t

   character(SL) :: mdout ! master only
   character(SL) :: restrt ! master only
   character(SL) :: mdvel ! master only
   character(SL) :: mden ! master only
   character(SL) :: mdcrd ! master only
   character(SL) :: mdinfo ! master only

   integer :: ioutfm
   integer :: ntpr
   integer :: ntwr
   integer :: ntwx

   character(SL) :: monitor_file ! master only
   character(SL) :: umbrella_file ! master only
   character(SL) :: snapshots_basename ! master only

   character(SL) :: monitor_fmt ! master only

   integer :: monitor_freq
   integer :: snapshots_freq ! master only

   NCSU_REAL :: timescale ! master only

   integer :: mode
   integer :: ncolvars

   type(colvar_t)   :: colvars(MAX_NUMBER_OF_COLVARS)
   type(umbrella_t) :: umbrella ! master only (sanderrank.eq.0)

   logical :: umbrella_file_existed ! home-master only
   logical :: umbrella_discretization_changed  ! home-master only

#ifndef NCSU_DISABLE_ASERT
   logical :: initialized = .false.
#endif /* NCSU_DISABLE_ASSERT */

end type bbmd_ctxt_t

public :: ctxt_init
public :: ctxt_fini

public :: ctxt_print
public :: ctxt_bcast

public :: ctxt_on_force
public :: ctxt_on_mdwrit

public :: ctxt_Um
public :: ctxt_Uo

public :: ctxt_send
public :: ctxt_recv

public :: ctxt_close_units
public :: ctxt_open_units

!-------------------------------------------------------------------------------

NCSU_REAL, parameter, private :: TINY = 0.00001d0

!-------------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------------

subroutine ctxt_init(self, root, amass)

   use ncsu_utils
   use ncsu_cftree
   use ncsu_colvar
   use ncsu_constants
   use ncsu_umbrella
   use ncsu_sander_proxy
   use file_io_dat

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self
   type(node_t), pointer :: root
   NCSU_REAL, intent(in) :: amass(*)

#  include "ncsu-mpi.h"

   type(child_t), pointer :: child

   logical :: umbrella_file_exists
   type(umbrella_t) :: umbrella_from_file

   logical :: found, do_transfer
   integer :: n, error

   integer :: cv_extents(UMBRELLA_MAX_NEXTENTS)
   logical :: cv_periodicity(UMBRELLA_MAX_NEXTENTS)

   NCSU_REAL :: cv_origin(UMBRELLA_MAX_NEXTENTS)
   NCSU_REAL :: cv_spacing(UMBRELLA_MAX_NEXTENTS)

   NCSU_REAL :: tmp

   character(len = SL) :: astring

   ncsu_assert(.not.self%initialized)

   if (sanderrank.gt.0) &
      goto 1

   ncsu_assert(associated(root))

   ! store SANDER's filenames

   self%mdout  = mdout
   self%restrt = restrt
   self%mdvel  = mdvel
   self%mden   = mden
   self%mdcrd  = mdcrd
   self%mdinfo = mdinfo

   self%ioutfm = ioutfm
   self%ntpr = ntpr
   self%ntwr = ntwr
   self%ntwx = ntwx

   ! setup defaults

   write (unit = self%monitor_file, fmt = '(a,a,i3.3)') &
      DEFAULT_MONITOR_FILE, '-', (masterrank + 1)
   write (unit = self%umbrella_file, fmt = '(a,a,i3.3,a)') &
      DEFAULT_UMBRELLA_FILE, '-', (masterrank + 1), '.nc'
   write (unit = self%snapshots_basename, fmt = '(a,a,i3.3)') &
      DEFAULT_SNAPSHOTS_BASENAME, '-', (masterrank + 1)

   self%ncolvars = 0

   ! discover the run-mode

   found = node_lookup_string(root, 'mode', astring)
   if (.not.found) &
      call fatal('could not find ''mode'' in the '''//SECTION//''' section')

   if (astring == 'NONE') then
      self%mode = MODE_NONE
      goto 1
   else if (astring == 'ANALYSIS') then
      self%mode = MODE_ANALYSIS
   else if (astring == 'UMBRELLA') then
      self%mode = MODE_UMBRELLA
   else if (astring == 'FLOODING') then
      self%mode = MODE_FLOODING
   else
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
         NCSU_ERROR, 'unknown mode ''', trim(astring), ''''
      call terminate()
   end if

   found = node_lookup_string(root, 'umbrella_file', astring)
   if (found) &
      self%umbrella_file = astring

#ifdef NCSU_NO_NETCDF
   umbrella_file_exists = .false.
   write (unit = ERR_UNIT, fmt = '(a,a)') NCSU_WARNING, &
      'netCDF is not available (try ''-bintraj'' configure option)'
#else
   inquire (file = self%umbrella_file, exist = umbrella_file_exists)
#endif /* NCSU_NO_NETCDF */

   if (.not.umbrella_file_exists.and.self%mode.eq.MODE_UMBRELLA) then
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') NCSU_ERROR, '''', &
      trim(self%umbrella_file), ''' does not exist (required for UMBRELLA mode)'
      call terminate()
   end if

   if (self%mode.eq.MODE_ANALYSIS) &
      umbrella_file_exists = .false.

#ifndef NCSU_NO_NETCDF
   if (umbrella_file_exists) &
      call umbrella_load(umbrella_from_file, self%umbrella_file)
#endif /* NCSU_NO_NETCDF */

   ! collective variables
   ncsu_assert(self%ncolvars.eq.0)

   child => node_children(root)
   do while (associated(child))
      if (node_title(child%node) == 'variable') &
         self%ncolvars = self%ncolvars + 1
      child => child%next
   end do

   if (self%ncolvars.eq.0) &
      call fatal('no variable(s) in the '''//SECTION//''' section')

   if (self%ncolvars.gt.MAX_NUMBER_OF_COLVARS) &
      call fatal('too many variables in the '''//SECTION//''' section')

   if (umbrella_file_exists) then
      if(umbrella_nextents(umbrella_from_file).ne.self%ncolvars) &
         call fatal('number of variables in the '''//SECTION//''' does not &
                 &match with the number of extents found in the umbrella_file')
   end if ! umbrella_file_exists

   n = 1
   child => node_children(root)

   do while (associated(child))
      if (node_title(child%node) == 'variable') then

         ncsu_assert(n <= self%ncolvars)
         call colvar_mdread(self%colvars(n), child%node, n)

         if (self%mode.eq.MODE_FLOODING) then
            found = node_lookup_positive_real(child%node, &
                        'resolution', cv_spacing(n))
            cv_spacing(n) = cv_spacing(n)/4
            if (.not.found.and..not.umbrella_file_exists) then
               write (unit = ERR_UNIT, fmt = '(/a,a,i1/)') NCSU_ERROR, &
                  'could not determine ''resolution'' for CV #', n
               call terminate()
            end if

            if (.not.found) &
               cv_spacing(n) = umbrella_spacing(umbrella_from_file, n)

            cv_periodicity(n) = colvar_is_periodic(self%colvars(n))

            if (cv_periodicity(n)) then

               ncsu_assert(colvar_has_min(self%colvars(n)))
               ncsu_assert(colvar_has_max(self%colvars(n)))

               cv_origin(n) = colvar_min(self%colvars(n))

               ncsu_assert(cv_spacing(n).gt.ZERO)
               ncsu_assert(colvar_max(self%colvars(n)).gt.cv_origin(n))

               cv_extents(n) = &
               int((colvar_max(self%colvars(n)) - cv_origin(n))/cv_spacing(n))

               if (cv_extents(n).lt.UMBRELLA_MIN_EXTENT) then
                  write (unit = ERR_UNIT, fmt = '(/a,a,i1,a/)') NCSU_ERROR, &
                     'CV #', n, ' : ''resolution'' is too big'
                  call terminate()
               end if

               cv_spacing(n) = &
                  (colvar_max(self%colvars(n)) - cv_origin(n))/cv_extents(n)

            else ! .not.periodic

               found = node_lookup_real(child%node, 'min', cv_origin(n))
               if (.not.found) then
                  if (umbrella_file_exists) then
                     cv_origin(n) = umbrella_origin(umbrella_from_file, n)
                  else if (colvar_has_min(self%colvars(n))) then
                     cv_origin(n) = colvar_min(self%colvars(n))
                  else
                     write (unit = ERR_UNIT, fmt = '(/a,a,i1/)') NCSU_ERROR, &
                        'could not determine ''min'' for CV #', n
                     call terminate()
                  end if
               end if

               found = node_lookup_real(child%node, 'max', tmp)
               if (.not.found) then
                  if (umbrella_file_exists) then
                     tmp = umbrella_origin(umbrella_from_file, n) &
                        + umbrella_spacing(umbrella_from_file, n) &
                        * (umbrella_extent(umbrella_from_file, n) - 1)
                  else if (colvar_has_max(self%colvars(n))) then
                     tmp = colvar_max(self%colvars(n))
                  else
                     write (unit = ERR_UNIT, fmt = '(/a,a,i1/)') NCSU_ERROR, &
                        'could not determine ''max'' for CV #', n
                     call terminate()
                  end if
               end if

               if (cv_origin(n).ge.tmp) then
                  write (unit = ERR_UNIT, fmt = '(/a,a,i1/)') NCSU_ERROR, &
                     'min.ge.max for CV #', n
                  call terminate()
               end if

               ncsu_assert(cv_spacing(n).gt.ZERO)
               cv_extents(n) = 1 &
                  + int((tmp - cv_origin(n))/cv_spacing(n))

               if (cv_extents(n).lt.UMBRELLA_MIN_EXTENT) then
                  write (unit = ERR_UNIT, fmt = '(/a,a,i1,a/)') NCSU_ERROR, &
                     'CV #', n, ' : the ''resolution'' is too big'
                  call terminate()
               end if

               cv_spacing(n) = (tmp - cv_origin(n))/(cv_extents(n) - 1)

            end if ! cv_periodicity(n)
         end if ! mode.eq.MODE_FLOODING
         n = n + 1
      end if ! title == 'variable'
      child => child%next
   end do

   ! monitor
   found = node_lookup_string(root, 'monitor_file', astring)
   if (found) &
      self%monitor_file = astring

   found = node_lookup_positive_integer(root, 'monitor_freq', self%monitor_freq)
   if (.not.found) &
      self%monitor_freq = 50

   self%monitor_freq = min(self%monitor_freq, sander_nstlim())
   self%monitor_freq = max(1, self%monitor_freq)

   ! umbrella snapshots
   found = node_lookup_string(root, 'snapshots_basename', astring)
   if (found) &
      self%snapshots_basename = astring

   found = node_lookup_integer(root, 'snapshots_freq', self%snapshots_freq)
   if (.not.found) &
      self%snapshots_freq = -1 ! no snapshots

   if (self%mode.eq.MODE_FLOODING) then
      found = node_lookup_positive_real(root, 'timescale', self%timescale)
      if (.not.found) &
         call fatal('could not find ''timescale'' &
            &in the '''//SECTION//''' section')

   end if ! mode.eq.MODE_FLOODING

1  continue ! sanderrank.gt.0 jumps here

   call mpi_bcast(self%mode, 1, MPI_INTEGER, 0, commsander, error)
   ncsu_assert(error.eq.0)

   if (self%mode.eq.MODE_NONE) &
      goto 2

   ncsu_assert(.not.is_master().or.self%ncolvars.gt.0)

   call mpi_bcast(self%ncolvars, 1, MPI_INTEGER, 0, commsander, error)
   ncsu_assert(error.eq.0)

   call mpi_bcast(self%monitor_freq, 1, MPI_INTEGER, 0, commsander, error)
   ncsu_assert(error.eq.0)

   ncsu_assert(self%ncolvars.gt.0)
   ncsu_assert(self%ncolvars.le.MAX_NUMBER_OF_COLVARS)

   if (multisander_numwatkeep().gt.0) &
      call fatal('numwatkeep.gt.0 is not supported')

   if (sander_imin().ne.0) &
      call fatal('imin.ne.0 is not supported')

   do n = 1, self%ncolvars
      call colvar_bootstrap(self%colvars(n), n, amass)
   end do

   if (sanderrank.gt.0) &
      goto 2

   do_transfer = .false.

   if (self%mode.eq.MODE_UMBRELLA) then
      ncsu_assert(umbrella_file_exists)
      do_transfer = .false.
      call umbrella_swap(self%umbrella, umbrella_from_file)
   else if (self%mode.eq.MODE_FLOODING) then
      if (umbrella_file_exists) then
         do_transfer = .false.
         do n = 1, self%ncolvars
            do_transfer = do_transfer &
               .or.(cv_extents(n).ne.umbrella_extent(umbrella_from_file, n))
            do_transfer = do_transfer &
               .or.(cv_periodicity(n).neqv.&
                  umbrella_periodicity(umbrella_from_file, n))
            do_transfer = do_transfer &
               .or.(abs(cv_origin(n) - umbrella_origin(umbrella_from_file, n)) &
                  .gt.TINY)
            do_transfer = do_transfer &
               .or.(abs(cv_spacing(n) - umbrella_spacing(umbrella_from_file, &
                  n)).gt.TINY)
            if (do_transfer) &
               exit
         end do
         if (do_transfer) then
            call umbrella_init(self%umbrella, self%ncolvars, cv_extents, &
                               cv_origin, cv_spacing, cv_periodicity)
            call umbrella_transfer(self%umbrella, umbrella_from_file)
            call umbrella_fini(umbrella_from_file)
         else
            call umbrella_swap(self%umbrella, umbrella_from_file)
         end if ! do_transfer
      else
         call umbrella_init(self%umbrella, self%ncolvars, cv_extents, &
                            cv_origin, cv_spacing, cv_periodicity)
      end if ! umbrella_file_exits
   end if ! self%mode.eq.MODE_FLOODING

   self%umbrella_file_existed = umbrella_file_exists
   self%umbrella_discretization_changed = do_transfer

   ! prepare monitor_fmt & open MONITOR_UNIT

   open (unit = MONITOR_UNIT, file = self%monitor_file, iostat = error, &
         form = 'FORMATTED', action = 'WRITE', status = 'REPLACE')

   if (error.ne.0) then
      write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
         NCSU_ERROR, 'failed to open ''', trim(self%monitor_file), &
         ''' for writing'
      call terminate()
   end if

   write (unit = MONITOR_UNIT, fmt = '(a,/a)', advance = 'NO') &
      '#', '# MD time (ps), '
   do n = 1, self%ncolvars - 1
      write (unit = MONITOR_UNIT, fmt = '(a,i1,a)', advance = 'NO') &
         'CV #', n, ', '
   end do

   if (self%mode == MODE_FLOODING) then
      write (unit = MONITOR_UNIT, fmt = '(a,i1,a,/a)') &
         'CV #', self%ncolvars, ', E_{bias} (kcal/mol)', '#'
      write (unit = self%monitor_fmt, fmt = '(a,i1,a)') &
         '(f12.4,', self%ncolvars, '(1x,f16.10),1x,f16.10)'
   else
      write (unit = MONITOR_UNIT, fmt = '(a,i1,/a)') &
         'CV #', self%ncolvars, '#'
      write (unit = self%monitor_fmt, fmt = '(a,i1,a)') &
         '(f12.4,', self%ncolvars, '(1x,f16.10))'
   end if

   call flush_UNIT(MONITOR_UNIT)

2  continue ! sanderrank.gt.0 jump here (or mode.eq.MODE_NONE)

#  ifndef NCSU_DISABLE_ASSERT
   self%initialized = .true.
#  endif /* NCSU_DISABLE_ASSERT */

end subroutine ctxt_init

!-------------------------------------------------------------------------------

subroutine ctxt_fini(self)

   NCSU_USE_AFAILED

   use ncsu_colvar
   use ncsu_umbrella
   use ncsu_sander_proxy

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self

#  include "ncsu-mpi.h"

   integer :: n

   ncsu_assert(self%initialized)

   if (self%mode.ne.MODE_NONE) then
      ncsu_assert(self%ncolvars.gt.0)
      do n = 1, self%ncolvars
         call colvar_cleanup(self%colvars(n))
      end do
   end if

   if (sanderrank.eq.0) then
      if (self%mode.eq.MODE_FLOODING.or.self%mode.eq.MODE_UMBRELLA) &
         call umbrella_fini(self%umbrella)
   end if ! sanderrank.eq.0

   self%mode = MODE_NONE

#  ifndef NCSU_DISABLE_ASSERT
   self%initialized = .false.
#  endif /* NCSU_DISABLE_ASSERT */

end subroutine ctxt_fini

!-------------------------------------------------------------------------------

subroutine ctxt_print(self, lun)

   use ncsu_utils
   use ncsu_colvar
   use ncsu_umbrella
   use ncsu_constants
   use ncsu_sander_proxy

   implicit none

   type(bbmd_ctxt_t), intent(in) :: self
   integer, intent(in) :: lun

   integer :: n
   NCSU_REAL :: tmp

   ncsu_assert(self%initialized)

   write (unit = lun, fmt = '(a,a)', advance = 'NO') NCSU_INFO, 'mode = '

   select case(self%mode)
      case(MODE_NONE)
         write (unit = lun, fmt = '(a)') 'NONE'
         goto 1
      case(MODE_ANALYSIS)
         write (unit = lun, fmt = '(a)') 'ANALYSIS'
      case(MODE_UMBRELLA)
         write (unit = lun, fmt = '(a)') 'UMBRELLA'
      case(MODE_FLOODING)
         write (unit = lun, fmt = '(a)') 'FLOODING'
      case default
         ncsu_assert_not_reached()
         continue
   end select

   write (unit = lun, fmt = '(a)') NCSU_INFO

   do n = 1, self%ncolvars
      write (unit = lun, fmt = '(a,a,i1)') NCSU_INFO, 'CV #', n
      call colvar_print(self%colvars(n), lun)
      write (unit = lun, fmt = '(a)') NCSU_INFO
   end do

   write (unit = lun, fmt = '(a,a,a)') NCSU_INFO, &
      'monitor_file = ', trim(self%monitor_file)
   write (unit = lun, fmt = '(a,a,'//pfmt &
      (self%monitor_freq)//',a,'//pfmt &
      (self%monitor_freq*sander_timestep(), 4)//',a)') NCSU_INFO, &
      'monitor_freq = ', self%monitor_freq, ' (', &
      self%monitor_freq*sander_timestep(), ' ps)'

   if (self%mode.eq.MODE_ANALYSIS) &
      goto 1

   write (unit = lun, fmt = '(a,a,a,a)', advance = 'NO') NCSU_INFO, &
      'umbrella_file = ', trim(self%umbrella_file), ' ('

   if (self%umbrella_file_existed) then
      write (unit = lun, fmt = '(a)') 'loaded)'
   else
      write (unit = lun, fmt = '(a)') 'not found)'
   end if

   write (unit = lun, fmt = '(a)') NCSU_INFO
   write (unit = lun, fmt = '(a,a)', advance = 'NO') NCSU_INFO, &
      'umbrella discretization '

   if (self%umbrella_file_existed) then
      if (self%umbrella_discretization_changed) then
         write (unit = lun, fmt = '(a)') '(modified) :'
      else
         write (unit = lun, fmt = '(a)') '(unchanged) :'
      end if
   else
      write (unit = lun, fmt = '(a)') '(new) :'
   end if

   do n = 1, self%ncolvars
      write (unit = lun, fmt = '(a,a,i1)', advance = 'NO') &
         NCSU_INFO, 'CV #', n
      if (umbrella_periodicity(self%umbrella, n)) then
         write (unit = lun, fmt = '(a)', advance = 'NO') ' periodic, '
         tmp = umbrella_origin(self%umbrella, n) &
         + umbrella_spacing(self%umbrella, n)*umbrella_extent(self%umbrella, n)
      else
         write (unit = lun, fmt = '(a)', advance = 'NO') ' not periodic, '
         tmp = umbrella_origin(self%umbrella, n) &
         + umbrella_spacing(self%umbrella, n)*(umbrella_extent(self%umbrella, n) - 1)
      end if

      write (unit = lun, &
         fmt = '('//pfmt(umbrella_extent(self%umbrella, n))//',a,'//pfmt &
         (umbrella_origin(self%umbrella, n), 6)//',a,'//pfmt(tmp, 6)//')') &
         umbrella_extent(self%umbrella, n), ' points, min/max = ', &
         umbrella_origin(self%umbrella, n), '/', tmp
   end do

   if (self%mode.eq.MODE_UMBRELLA) &
      goto 1

   write (unit = lun, fmt = '(a/,a,a,'//pfmt(self%timescale, 3)//',a)') &
      NCSU_INFO, NCSU_INFO, 'flooding timescale = ', self%timescale, ' ps'

   if (self%snapshots_freq.gt.0) then
      write (unit = lun, fmt = '(a,a,a)') NCSU_INFO, &
         'snapshots_basename = ', trim(self%snapshots_basename)
      write (unit = lun, &
         fmt = '(a,a,'//pfmt(self%snapshots_freq)//',a,'//pfmt &
         (self%snapshots_freq*sander_timestep(), 4)//',a)') NCSU_INFO, &
         'snapshots_freq = ', self%snapshots_freq, ' (', &
         self%snapshots_freq*sander_timestep(), ' ps)'
   end if

1  call flush_UNIT(lun)

end subroutine ctxt_print

!-------------------------------------------------------------------------------

subroutine ctxt_bcast(self, masterroot, amass)

   use ncsu_utils
   use ncsu_umbrella, only : umbrella_bcast

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self
   integer, intent(in) :: masterroot
   NCSU_REAL, intent(in) :: amass(*)

#  include "ncsu-mpi.h"

   integer :: n, error

   ncsu_assert(masterroot.ge.0.and.masterroot.lt.mastersize)

#  ifndef NCSU_DISABLE_ASSERT
   if (masterroot.eq.masterrank) then
      ncsu_assert(self%initialized)
   else
      self%initialized = .true.
   end if ! masterroot.eq.masterrank
#  endif /* NCSU_DISABLE_ASSERT */

   ! SANDER's files (no matter what the mode is)

   if (sanderrank.eq.0) then
      ncsu_assert(commmaster.ne.MPI_COMM_NULL)

      call mpi_bcast(self%mdout, len(self%mdout), MPI_CHARACTER, &
                     masterroot, commmaster, error)
      ncsu_assert(error.eq.0)

      call mpi_bcast(self%restrt, len(self%restrt), MPI_CHARACTER, &
                     masterroot, commmaster, error)
      ncsu_assert(error.eq.0)

      call mpi_bcast(self%mdvel, len(self%mdvel), MPI_CHARACTER, &
                     masterroot, commmaster, error)
      ncsu_assert(error.eq.0)

      call mpi_bcast(self%mden, len(self%mden), MPI_CHARACTER, &
                     masterroot, commmaster, error)
      ncsu_assert(error.eq.0)

      call mpi_bcast(self%mdcrd, len(self%mdcrd), MPI_CHARACTER, &
                     masterroot, commmaster, error)
      ncsu_assert(error.eq.0)

      call mpi_bcast(self%mdinfo, len(self%mdinfo), MPI_CHARACTER, &
                     masterroot, commmaster, error)
      ncsu_assert(error.eq.0)

      call mpi_bcast(self%ioutfm, 1, MPI_INTEGER, masterroot, commmaster, error)
      ncsu_assert(error.eq.0)

      call mpi_bcast(self%ntpr, 1, MPI_INTEGER, masterroot, commmaster, error)
      ncsu_assert(error.eq.0)

      call mpi_bcast(self%ntwr, 1, MPI_INTEGER, masterroot, commmaster, error)
      ncsu_assert(error.eq.0)

      call mpi_bcast(self%ntwx, 1, MPI_INTEGER, masterroot, commmaster, error)
      ncsu_assert(error.eq.0)

      call mpi_bcast(self%mode, 1, MPI_INTEGER, masterroot, commmaster, error)
      ncsu_assert(error.eq.0)

   end if ! sanderrank.eq.0

   if (masterrank.ne.masterroot) then
      call mpi_bcast(self%mode, 1, MPI_INTEGER, 0, commsander, error)
      ncsu_assert(error.eq.0)
   end if ! masterrank.ne.masterroot

   if (self%mode.eq.MODE_NONE) &
      return

   ! CVs

   ncsu_assert(self%ncolvars.gt.0.or.masterrank.ne.masterroot)

   if (sanderrank.eq.0) then
      call mpi_bcast(self%ncolvars, 1, &
                     MPI_INTEGER, masterroot, commmaster, error)
      ncsu_assert(error.eq.0)
   end if ! sanderrank.eq.0

   if (masterrank.ne.masterroot) then
      call mpi_bcast(self%ncolvars, 1, MPI_INTEGER, 0, commsander, error)
      ncsu_assert(error.eq.0)
   end if ! masterrank.ne.masterroot

   ncsu_assert(self%ncolvars.gt.0)

   do n = 1, self%ncolvars
      call bcast_colvar(self%colvars(n), n + 10*masterrank)
   end do

   ! mode = ANALYSIS

   if (sanderrank.eq.0) then
      call mpi_bcast(self%monitor_file, len(self%monitor_file), &
                     MPI_CHARACTER, masterroot, commmaster, error)
      ncsu_assert(error.eq.0)

      call mpi_bcast(self%monitor_fmt, len(self%monitor_fmt), &
                     MPI_CHARACTER, masterroot, commmaster, error)
      ncsu_assert(error.eq.0)

      call mpi_bcast(self%monitor_freq, 1, &
                     MPI_INTEGER, masterroot, commmaster, error)
      ncsu_assert(error.eq.0)
   end if ! sanderrank.eq.0

   if (masterrank.ne.masterroot) then
      call mpi_bcast(self%monitor_freq, 1, MPI_INTEGER, 0, commsander, error)
      ncsu_assert(error.eq.0)
   end if ! masterrank.ne.masterroot

   if (self%mode.eq.MODE_ANALYSIS) &
      return

   ! mode = UMBRELLA | FLOODING (these are on masters only)

   if (sanderrank.eq.0) then
      call mpi_bcast(self%umbrella_file, len(self%umbrella_file), &
                     MPI_CHARACTER, masterroot, commmaster, error)
      ncsu_assert(error.eq.0)

      call mpi_bcast(self%snapshots_basename, len(self%snapshots_basename), &
                     MPI_CHARACTER, masterroot, commmaster, error)
      ncsu_assert(error.eq.0)

      call mpi_bcast(self%snapshots_freq, 1, &
                     MPI_INTEGER, masterroot, commmaster, error)
      ncsu_assert(error.eq.0)

      call mpi_bcast(self%timescale, 1, &
                     MPI_DOUBLE_PRECISION, masterroot, commmaster, error)
      ncsu_assert(error.eq.0)

      call umbrella_bcast(self%umbrella, commmaster, masterroot)

   end if ! sanderrank.eq.0

contains

subroutine bcast_colvar(cv, cvno)

   use ncsu_colvar, only : colvar_bootstrap

   implicit none

   type(colvar_t), intent(inout) :: cv
   integer, intent(in) :: cvno

   integer :: bcastdata(3)

   if (sanderrank.eq.0) then

      if (masterrank.eq.masterroot) then

         ncsu_assert(cv%type.gt.0)

         bcastdata(1) = cv%type

         bcastdata(2) = 0
         if (associated(cv%i)) &
            bcastdata(2) = size(cv%i)

         bcastdata(3) = 0
         if (associated(cv%r)) &
            bcastdata(3) = size(cv%r)

      end if ! masterrank.eq.masterroot

      call mpi_bcast(bcastdata, 3, MPI_INTEGER, masterroot, commmaster, error)
      ncsu_assert(error.eq.0)

      if (masterrank.ne.masterroot) then
         cv%type = bcastdata(1)

         if (bcastdata(2).gt.0) then
            allocate(cv%i(bcastdata(2)), stat = error)
            if (error.ne.0) &
               NCSU_OUT_OF_MEMORY
         end if ! bcastdata(2).gt.0

         if (bcastdata(3).gt.0) then
            allocate(cv%r(bcastdata(3)), stat = error)
            if (error.ne.0) &
               NCSU_OUT_OF_MEMORY
         end if ! bcastdata(3).gt.0

      end if ! masterrank.ne.masterroot

      if (bcastdata(2).gt.0) then
         call mpi_bcast(cv%i, bcastdata(2), MPI_INTEGER, &
                        masterroot, commmaster, error)
         ncsu_assert(error.eq.0)
      end if ! bcastdata(2).gt.0

      if (bcastdata(3).gt.0) then
         call mpi_bcast(cv%r, bcastdata(3), MPI_DOUBLE_PRECISION, &
                        masterroot, commmaster, error)
         ncsu_assert(error.eq.0)
      end if ! bcastdata(3).gt.0

   end if ! sanderrank.eq.0

   if (masterrank.ne.masterroot) &
      call colvar_bootstrap(cv, cvno, amass)

end subroutine bcast_colvar

end subroutine ctxt_bcast

!-------------------------------------------------------------------------------

subroutine ctxt_on_force(self, x, f, mdstep)

   use ncsu_utils
   use ncsu_colvar
   use ncsu_umbrella
   use ncsu_constants
   use ncsu_sander_proxy

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self

   NCSU_REAL, intent(in) :: x(*)

   NCSU_REAL, intent(inout) :: f(*)

   integer, intent(in) :: mdstep

#  include "ncsu-mpi.h"

   NCSU_REAL :: u_derivative(UMBRELLA_MAX_NEXTENTS)
   NCSU_REAL :: instantaneous(UMBRELLA_MAX_NEXTENTS)
   NCSU_REAL :: alt, u_value

   character(len = SL + 16) :: snapshot

   integer :: n, error
   logical :: real_mdstep

   ncsu_assert(self%initialized)

   if (self%mode.eq.MODE_NONE) &
      return

   real_mdstep = (sander_init().eq.4)

   ncsu_assert(self%ncolvars.gt.0)

   if (self%mode.eq.MODE_ANALYSIS) then
      if (real_mdstep.and.mod(mdstep, self%monitor_freq).eq.0) then
         do n = 1, self%ncolvars
            instantaneous(n) = colvar_value(self%colvars(n), x)
         end do

         if (sanderrank.eq.0) then
            write (unit = MONITOR_UNIT, fmt = self%monitor_fmt) &
               sander_mdtime(), instantaneous(1:self%ncolvars)
            call flush_UNIT(MONITOR_UNIT)
         end if ! sanderrank.eq.0
      end if

      return
   end if ! self%mode.eq.MODE_ANALYSIS

   !
   ! either UMBRELLA or FLOODING
   !

   do n = 1, self%ncolvars
      instantaneous(n) = colvar_value(self%colvars(n), x)
   end do

   if (sanderrank.eq.0) &
      call umbrella_eval_vdv(self%umbrella, instantaneous, &
                             u_value, u_derivative)

   call mpi_bcast(u_derivative, self%ncolvars, &
      MPI_DOUBLE_PRECISION, 0, commsander, error)
   ncsu_assert(error.eq.0)

   ! FIXME: virial
   do n = 1, self%ncolvars
      call colvar_force(self%colvars(n), x, -u_derivative(n), f)
   end do

   if (.not.real_mdstep.or.sanderrank.ne.0) &
      return

   ncsu_assert(self%mode.eq.MODE_UMBRELLA.or.self%mode.eq.MODE_FLOODING)

   if (mod(mdstep, self%monitor_freq).eq.0) then
      if (self%mode.eq.MODE_FLOODING) then
         write (unit = MONITOR_UNIT, fmt = self%monitor_fmt) &
            sander_mdtime(), instantaneous(1:self%ncolvars), u_value
      else
         write (unit = MONITOR_UNIT, fmt = self%monitor_fmt) &
            sander_mdtime(), instantaneous(1:self%ncolvars)
      end if
      call flush_UNIT(MONITOR_UNIT)
   end if

   if (self%mode.eq.MODE_FLOODING) then
#  ifndef NCSU_NO_NETCDF
      if (self%snapshots_freq.gt.0 &
          .and.mod(mdstep, self%snapshots_freq).eq.0) then
         write (unit = snapshot, fmt = '(a,a,i10.10,a)') &
            trim(self%snapshots_basename), '.', mdstep, '.nc'
         call umbrella_save(self%umbrella, snapshot)
         write (unit = OUT_UNIT, fmt = '(/a,a,'//pfmt &
            (sander_mdtime(), 4)//',a,/a,a,a,a)') &
            NCSU_INFO, 'biasing potential snapshot at t = ', &
            sander_mdtime(), ' ps', NCSU_INFO, 'saved as ''', &
            trim(snapshot), ''''
      end if
#  endif /* NCSU_NO_NETCDF */
      alt = sander_timestep()/self%timescale
      call umbrella_hill(self%umbrella, instantaneous, alt)
   end if ! self%mode.eq.MODE_FLOODING

end subroutine ctxt_on_force

!-------------------------------------------------------------------------------

subroutine ctxt_on_mdwrit(self)

   use ncsu_utils
   use ncsu_umbrella
   use ncsu_sander_proxy

   implicit none

   type(bbmd_ctxt_t), intent(in) :: self

   ncsu_assert(self%initialized)

#ifndef NCSU_NO_NETCDF
   if (self%mode.eq.MODE_FLOODING) then
      ncsu_assert(is_master())
      call umbrella_save(self%umbrella, self%umbrella_file)
   end if
#endif /* NCSU_NO_NETCDF */

end subroutine ctxt_on_mdwrit

!-------------------------------------------------------------------------------

!
! U_m? are valid for masterrank.lt.r_masterrank
!

! assumes that local self is up to date
subroutine ctxt_Um(self, r_masterrank, x, U_mm, U_mo)

   NCSU_USE_AFAILED

   use ncsu_colvar
   use ncsu_umbrella
   use ncsu_constants

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self

   integer, intent(in) :: r_masterrank
   NCSU_REAL, intent(in) :: x(*)

   NCSU_REAL, intent(out) :: U_mm, U_mo

#  include "ncsu-mpi.h"

   NCSU_REAL :: local_inst(UMBRELLA_MAX_NEXTENTS)
   NCSU_REAL :: remote_inst(UMBRELLA_MAX_NEXTENTS)

   integer :: n, error
   ncsu_assert(self%initialized)

   U_mm = ZERO
   U_mo = ZERO

   if (self%mode.eq.MODE_NONE.or.self%mode.eq.MODE_ANALYSIS) &
      return

   ncsu_assert(self%ncolvars.gt.0)
   do n = 1, self%ncolvars
      local_inst(n) = colvar_value(self%colvars(n), x)
   end do

   if (sanderrank.gt.0) &
      return

   ncsu_assert(self%mode.eq.MODE_UMBRELLA.or.self%mode.eq.MODE_FLOODING)

   if (masterrank.lt.r_masterrank) then
      call mpi_recv(remote_inst, self%ncolvars, MPI_DOUBLE_PRECISION, &
                    r_masterrank, 0, commmaster, MPI_STATUS_IGNORE, error)
      ncsu_assert(error.eq.0)
      U_mm = umbrella_eval_v(self%umbrella, local_inst)
      U_mo = umbrella_eval_v(self%umbrella, remote_inst)
   else
      call mpi_send(local_inst, self%ncolvars, MPI_DOUBLE_PRECISION, &
                    r_masterrank, 0, commmaster, error)
      ncsu_assert(error.eq.0)
   end if ! masterrank.lt.r_masterrank

end subroutine ctxt_Um

! assumes that remote self is up to date (if mode.eq.MODE_FLOODING)
subroutine ctxt_Uo(self, r_masterrank, x, U_om, U_oo)

   NCSU_USE_AFAILED

   use ncsu_colvar
   use ncsu_umbrella
   use ncsu_constants

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self

   integer, intent(in) :: r_masterrank
   NCSU_REAL, intent(in) :: x(*)

   NCSU_REAL, intent(out) :: U_om, U_oo

#  include "ncsu-mpi.h"

   NCSU_REAL :: local_inst(UMBRELLA_MAX_NEXTENTS)
   NCSU_REAL :: remote_inst(UMBRELLA_MAX_NEXTENTS)

   NCSU_REAL :: tmp(2)

   integer :: n, error

   ncsu_assert(self%initialized)

   U_om = ZERO
   U_oo = ZERO

   if (self%mode.eq.MODE_NONE.or.self%mode.eq.MODE_ANALYSIS) &
      return

   ncsu_assert(self%ncolvars.gt.0)
   do n = 1, self%ncolvars
      local_inst(n) = colvar_value(self%colvars(n), x)
   end do

   if (sanderrank.gt.0) &
      return

   if (self%mode.eq.MODE_UMBRELLA) then
      ! remote umbrella is same as local
      if (masterrank.lt.r_masterrank) then
         call mpi_recv(remote_inst, self%ncolvars, MPI_DOUBLE_PRECISION, &
                       r_masterrank, 0, commmaster, MPI_STATUS_IGNORE, error)
         ncsu_assert(error.eq.0)
         U_om = umbrella_eval_v(self%umbrella, local_inst)
         U_oo = umbrella_eval_v(self%umbrella, remote_inst)
      else
         call mpi_send(local_inst, self%ncolvars, MPI_DOUBLE_PRECISION, &
                       r_masterrank, 0, commmaster, error)
         ncsu_assert(error.eq.0)
      end if ! masterrank.lt.r_masterrank
   else
      ncsu_assert(self%mode.eq.MODE_FLOODING)
      if (masterrank.gt.r_masterrank) then
         call mpi_recv(remote_inst, self%ncolvars, MPI_DOUBLE_PRECISION, &
                       r_masterrank, 0, commmaster, MPI_STATUS_IGNORE, error)
         ncsu_assert(error.eq.0)
         tmp(1) = umbrella_eval_v(self%umbrella, local_inst)
         tmp(2) = umbrella_eval_v(self%umbrella, remote_inst)
         call mpi_send(tmp, 2, MPI_DOUBLE_PRECISION, &
                       r_masterrank, 1, commmaster, error)
         ncsu_assert(error.eq.0)
      else
         call mpi_send(local_inst, self%ncolvars, MPI_DOUBLE_PRECISION, &
                       r_masterrank, 0, commmaster, error)
         ncsu_assert(error.eq.0)
         call mpi_recv(tmp, 2, MPI_DOUBLE_PRECISION, &
                       r_masterrank, 1, commmaster, MPI_STATUS_IGNORE, error)
         ncsu_assert(error.eq.0)
         U_om = tmp(2)
         U_oo = tmp(1)
      end if ! masterrank.gt.r_masterrank
   end if ! self%mode.eq.MODE_UMBRELLA

end subroutine ctxt_Uo

!-------------------------------------------------------------------------------

subroutine ctxt_send(self, dst_masterrank)

   NCSU_USE_AFAILED

   use ncsu_umbrella, only : umbrella_send_coeffs

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self
   integer, intent(in) :: dst_masterrank

   integer :: error

#  include "ncsu-mpi.h"

   ncsu_assert(sanderrank.eq.0)
   ncsu_assert(self%initialized)

   if (self%mode.eq.MODE_FLOODING) then
      call umbrella_send_coeffs(self%umbrella, dst_masterrank, commmaster)
   else
      ! for synchronization
      call mpi_send(masterrank, 1, MPI_INTEGER, &
                    dst_masterrank, 8, commmaster, error)
      ncsu_assert(error.eq.0)
   end if

end subroutine ctxt_send

!-------------------------------------------------------------------------------

subroutine ctxt_recv(self, src_masterrank)

   NCSU_USE_AFAILED

   use ncsu_umbrella, only : umbrella_recv_coeffs

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self
   integer, intent(in) :: src_masterrank

#  include "ncsu-mpi.h"

   integer :: error, itmp

   ncsu_assert(sanderrank.eq.0)
   ncsu_assert(self%initialized)

   if (self%mode.eq.MODE_FLOODING) then
      call umbrella_recv_coeffs(self%umbrella, src_masterrank, commmaster)
   else
      ! for synchronization
      call mpi_recv(itmp, 1, MPI_INTEGER, src_masterrank, 8, &
                    commmaster, MPI_STATUS_IGNORE, error)
      ncsu_assert(error.eq.0)
      ncsu_assert(itmp.eq.src_masterrank)
   end if

end subroutine ctxt_recv

!-------------------------------------------------------------------------------

subroutine ctxt_close_units(self)

   NCSU_USE_AFAILED

   use ncsu_constants
   use ncsu_sander_proxy

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self

#  include "ncsu-mpi.h"

   ncsu_assert(self%initialized)
   ncsu_assert(sanderrank.eq.0)

   if (self%mode.ne.MODE_NONE) &
      close (MONITOR_UNIT)

   if (self%mdout.ne.'stdout') &
      close (OUT_UNIT)

   call close_dump_files()

end subroutine ctxt_close_units

!-------------------------------------------------------------------------------

subroutine ctxt_open_units(self)

   NCSU_USE_AFAILED

   use ncsu_constants
   use ncsu_sander_proxy
   use file_io_dat

   implicit none

   type(bbmd_ctxt_t), intent(inout) :: self

#  include "ncsu-mpi.h"

   integer :: error

   ncsu_assert(self%initialized)
   ncsu_assert(sanderrank.eq.0)

   if (self%mode.ne.MODE_NONE) then
      open (unit = MONITOR_UNIT, file = self%monitor_file, iostat = error, &
        form = 'FORMATTED', action = 'WRITE', status = 'OLD', position = 'APPEND')
      if (error.ne.0) then
         write (unit = ERR_UNIT, fmt = '(/a,a,a,a/)') &
            NCSU_ERROR, 'failed to open ''', trim(self%monitor_file), &
            ''' for appending'
         call terminate()
      end if
   end if ! self%mode.ne.MODE_NONE

   mdout = self%mdout
   mdinfo = self%mdinfo
   restrt = self%restrt
   mdvel = self%mdvel
   mden = self%mden
   mdcrd = self%mdcrd

   ioutfm = self%ioutfm

   ntpr = self%ntpr
   ntwr = self%ntwr
   ntwx = self%ntwx

   if (self%mdout.ne.'stdout') &
      call amopen(OUT_UNIT, mdout, 'O', 'F', 'A')

   call amopen(MDINFO_UNIT, mdinfo, 'U', 'F', 'W')

   facc = 'A'
   owrite = 'U'
   call open_dump_files()
   facc = 'W'

end subroutine ctxt_open_units

!-------------------------------------------------------------------------------

#endif /* NCSU_ENABLE_BBMD */

end module ncsu_bbmd_ctxt
