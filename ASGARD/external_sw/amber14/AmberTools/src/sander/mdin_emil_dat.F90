!*******************************************************************************
!
! Module: mdin_emil_dat_mod
!
! Description: read in (mostly optional) extra config information for EMIL calcs
!              
!*******************************************************************************


module mdin_emil_dat_mod

  use file_io_dat,      only : mdin, mdout
  implicit none

  !these strings pasted in from pmemd, only EMIL seems to use them currently
  !however it might cut the memory footprint minutely if other parts of sander also did?
  character(11), parameter, public :: extra_line_hdr = '|          '
  character(11), parameter, public :: warn_hdr =       '| WARNING: '
  character(11), parameter, public :: info_hdr =       '| INFO:    '
  character(11), parameter, public :: error_hdr =      '| ERROR:   '

#ifdef EMIL

  private


  ! Public variables
  character(256), public   :: emil_paramfile, emil_logfile, emil_model_infile, &
                              emil_model_outfile

  ! mdin_unit seems to be hardcoded as 5 elsewhere in sander, but without a variable name.
  integer                          :: mdin_unit = 5

  ! Subroutines
  public   init_emil_dat
  private  validate_and_log_filename 

  contains

  
!*******************************************************************************
!
! Subroutine:  init_emil_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine init_emil_dat(mdin_unit, mdout_unit)
   implicit none

   ! Arguments
   integer, intent(in) :: mdin_unit, mdout_unit

   ! Local variables:
   integer             :: ifind
   namelist /emil_cntrl/       emil_paramfile, emil_logfile, emil_model_infile, &
                               emil_model_outfile

   !let the emil module set the defaults if we have no preferences ourselves
   emil_paramfile     = ''
   emil_logfile       = ''
   emil_model_infile  = ''
   emil_model_outfile = ''


   ! Read input in namelist format:
   rewind(mdin_unit)                  ! Search the mdin file from the beginning
   call nmlsrc('emil_cntrl', mdin_unit, ifind)

   if (ifind .ne. 0) then        

      ! Namelist found. Read it:
      read(mdin_unit, nml = emil_cntrl)
      write(mdout_unit, '(a, a)') info_hdr,       'Found an "emil_cntrl" namelist'
      call validate_and_log_filename(mdout_unit, emil_paramfile, "emil_paramfile")
      call validate_and_log_filename(mdout_unit, emil_logfile,   "emil_logfile")
      call validate_and_log_filename(mdout_unit, emil_model_infile, "emil_model_infile")
      call validate_and_log_filename(mdout_unit, emil_model_outfile, "emil_model_outfile")
      write(mdout_unit, '(a)') ''
   else
      write(mdout_unit, '(a, a)') info_hdr,  &
        'Did not find an "emil_cntrl" namelist: assuming default emil i/o files'
      
   end if

end subroutine init_emil_dat

subroutine validate_and_log_filename(mdout_unit, filename, id_string)

   implicit none

   !arguments
   integer, intent(in)      :: mdout_unit
   character(*), intent(in) :: id_string
   character(*), intent(in) :: filename
 
   !local variables
   integer                  :: string_length

   string_length = len_trim(filename)
   if( string_length .ge. 255 ) then
           write(mdout_unit,'(a,a)') warn_hdr, &
                        id_string // ' is >= max length! Probably cropped.'
   end if
   if( string_length .gt. 0 ) then
           write(mdout_unit,'(a,a)') extra_line_hdr, &
                        id_string // ' set to: ' // trim(filename)
   else
           write(mdout_unit,'(a,a)') extra_line_hdr, &
                        id_string // ' not set in namelist, using default.'
   end if
       
end subroutine validate_and_log_filename

#endif
end module mdin_emil_dat_mod


