#include "copyright.h"

module commandline_module

    private

    public :: mdfil

contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Reads the command-line and fills the corresponding variables
subroutine mdfil(VERSION, version_requested)
   
   !     Author: George Seibel; many subsequent modifications by many others.
   use constantph, only : write_dump
   use file_io_dat
#ifdef MPI
   use remd, only : rem, rremd 
#endif
#ifdef RISMSANDER
   use sander_rism_interface, only : xvvfile, guvfile, huvfile, cuvfile, uuvfile,&
        asympfile, quvfile, chgdistfile, exchemfile, solvenefile, entropyfile,&
        exchemGFfile, solveneGFfile, entropyGFfile,&
        exchemUCfile, solveneUCfile, entropyUCfile,&
        potUVfile
#endif
   use cns_xref, only : is_xref_on

   implicit none
   
   ! If we are given the VERSION flag, then it's the first pass through mdfil

   character(4), optional, intent(in) :: VERSION
   logical, optional, intent(inout)   :: version_requested

   !     Modified for multisander to allow reading command line from a string
   !     rather than the command line.  Use iargc_wrap and getarg_wrap instead
   !     of the intrinsics.
   
   !     OUTPUT: (to common)
   
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#  include "../include/dprec.fh"
#endif
   
   !     INTERNAL:
   
   character(len=MAX_FN_LEN) arg
   character(len=MAX_FN_LEN) outfile_suffix
   character(len=MAX_FN_LEN) groupfile_holder
   !         temp for each of the whitespace delimited command line arguments
   integer iarg
   !         index of the current argument
   integer iargc_wrap
   !         wrapper to intrinsic that returns the index of the last argument
   !         from either the command line or a string
   integer last_arg_index
   !         index of the last argument
   logical :: suffix_specified   = .false.
   logical :: always_add_suffix  = .false.
   logical :: mdout_specified    = .false.
   logical :: restrt_specified   = .false.
   logical :: mdinfo_specified   = .false.
   logical :: mdcrd_specified    = .false.
   logical :: mdvel_specified    = .false.
   logical :: mdfrc_specified    = .false.
   logical :: mden_specified     = .false.
   logical :: cpout_specified    = .false.
   logical :: cprestrt_specified = .false.
   logical :: evbout_specified   = .false.
   logical :: amdlog_specified   = .false.
   logical :: scaledMDlog_specified   = .false.
   logical :: cph_dump_specified = .false.
   ! We need to know if these files have been specified on the command-line so 
   ! we know whether or not to add suffixes to them.

   !     --- default file names ---

   
   mdin   = 'mdin'
   mdout  = 'mdout'
   inpcrd = 'inpcrd'
   parm   = 'prmtop'
   restrt = 'restrt'
   refc   = 'refc'
   mdvel  = 'mdvel'
   mdfrc  = 'mdfrc'
   mden   = 'mden'
   mdcrd  = 'mdcrd'
   inptraj = 'inptraj'
   mdinfo = 'mdinfo'
   mtmd   = 'mtmd'
   vecs   = 'vecs'
   freqe  = 'dummy'
   rstdip = 'rstdip'
   inpdip = 'inpdip'
   mddip  = 'mddip'
   radii  = 'radii'
   cpin   = 'cpin'
   cpout  = 'cpout'
   cprestrt = 'cprestrt'
   evbin  = 'evbin'                                        ! EVB input file
   evbout = 'evbout'                                       ! EVB output file
   pimdout = 'pimdout'
   cph_dump = 'explicit_titration.dat'
#ifdef RISMSANDER
   xvvfile       = ''
   guvfile       = ''
   huvfile       = ''
   cuvfile       = ''
   uuvfile       = ''
   asympfile     = ''
   quvfile       = ''
   chgdistfile   = ''
   exchemfile    = ''
   solvenefile   = ''
   entropyfile   = ''
   exchemGFfile    = ''
   solveneGFfile   = ''
   entropyGFfile   = ''
   exchemUCfile    = ''
   solveneUCfile   = ''
   entropyUCfile   = ''
   potUVfile     = ''
   rismcrdfil    = ''
   rismfrcfil    = ''
   rismcrdrstfil = ''
   rismfrcrstfil = ''
#endif
   if (numgroup == 1) then
      groups(:) = ' '
      groupfile_holder(:) = ' '
#ifdef MPI
   else
      suffix_specified = .true.
      if (masterrank < 10) then
         write(outfile_suffix, '(a,i1)') '00', masterrank

      else if (masterrank < 100) then
         write(outfile_suffix, '(a,i2)') '0', masterrank

      else if (masterrank < 1000) then
         write(outfile_suffix, '(i3)') masterrank

      else if (masterrank < 10000) then
         write(outfile_suffix, '(i4)') masterrank

      else if (masterrank < 100000) then
         write(outfile_suffix, '(i5)') masterrank

      end if
#endif
   end if ! numgroup == 1

   !     --- default status of output: New
   
   owrite = 'N'

   facc = 'W'           ! default: overwriting
   is_xref_on = .false.   !  default for cns_xref

!AMD log file
   amdlog = 'amd.log'   ! default log file name for AMD

!scaledMD log file
   scaledMDlog = 'scaledMD.log'   ! default log file name for scaledMD

   !     --- get command line arguments ---
   
   iarg = 0
   last_arg_index = iargc_wrap()
   do while (iarg < last_arg_index)
      iarg = iarg + 1

      call getarg_wrap(iarg,arg)

      ! Look for a --version or -V flag and dump out the version of sander
      if (arg == '-V' .or. arg == '--version') then
         ! If VERSION is present, print out our program name and the version
         if (present(VERSION)) then
            call getarg_wrap(0, arg)
            call basename(arg) ! gets rid of the preceding path
            write(6, '(3a)') trim(arg), ': Version ', trim(VERSION)
            version_requested = .true.
            return
         else
            write(0, '(2a)') trim(arg), ' is a command-line flag! Not for use &
                                 &in the groupfile'
         end if
      else if (arg == '-O') then
         owrite = 'R' !      status of output: Replace
      else if (arg == '-A') then
         owrite = 'U' !      status of output: unknown
         facc = 'A'
      else if (arg == '-i') then
         iarg = iarg + 1
         call getarg_wrap(iarg,mdin)
      else if (arg == '-o') then
         mdout_specified = .true.
         iarg = iarg + 1
         call getarg_wrap(iarg,mdout)
      else if (arg == '-p') then
         iarg = iarg + 1
         call getarg_wrap(iarg,parm)
      else if (arg == '-c') then
         iarg = iarg + 1
         call getarg_wrap(iarg,inpcrd)
      else if (arg == '-vecs') then
         iarg = iarg + 1
         call getarg_wrap(iarg,vecs)
      else if (arg == '-radii') then
         iarg = iarg + 1
         call getarg_wrap(iarg,radii)
      else if (arg == '-f') then
         iarg = iarg + 1
         call getarg_wrap(iarg,freqe)
      else if (arg == '-r') then
         restrt_specified = .true.
         iarg = iarg + 1
         call getarg_wrap(iarg,restrt)
      else if (arg == '-ref' .or. arg == '-z') then
         iarg = iarg + 1
         call getarg_wrap(iarg,refc)
      else if (arg == '-e') then
         mden_specified = .true.
         iarg = iarg + 1
         call getarg_wrap(iarg,mden)
      else if (arg == '-v') then
         mdvel_specified = .true.
         iarg = iarg + 1
         call getarg_wrap(iarg,mdvel)
      else if (arg == '-x'.or.arg == '-t') then
         mdcrd_specified = .true.
         iarg = iarg + 1
         call getarg_wrap(iarg,mdcrd)
      else if (arg == '-frc') then
         mdfrc_specified = .true.
         iarg = iarg + 1
         call getarg_wrap(iarg,mdfrc)
      else if (arg == '-y') then
         iarg = iarg + 1
         call getarg_wrap(iarg,inptraj)
      else if (arg == '-inf') then
         mdinfo_specified = .true.
         iarg = iarg + 1
         call getarg_wrap(iarg,mdinfo)
      else if (arg == '-mtmd') then
         iarg = iarg + 1
         call getarg_wrap(iarg,mtmd)
      else if (arg == '-idip') then
         iarg = iarg + 1
         call getarg_wrap(iarg,inpdip)
      else if (arg == '-rdip') then
         iarg = iarg + 1
         call getarg_wrap(iarg,rstdip)
      else if (arg == '-mdip') then
         iarg = iarg + 1
         call getarg_wrap(iarg,mddip)
      else if (arg == '-cpin') then
         iarg = iarg + 1
         call getarg_wrap(iarg,cpin)
      else if (arg == '-cpout') then
         cpout_specified = .true.
         iarg = iarg + 1
         call getarg_wrap(iarg,cpout)
      else if (arg == '-cprestrt') then
         cprestrt_specified = .true.
         iarg = iarg + 1
         call getarg_wrap(iarg,cprestrt)
      else if (arg == '-cph-data') then
         write_dump = .true.
         cph_dump_specified = .true.
         iarg = iarg + 1
         call getarg_wrap(iarg, cph_dump)
      else if (arg == '-evbin') then                       ! EVB input file
         iarg = iarg + 1
         call getarg_wrap(iarg,evbin)
      else if (arg == '-evbout') then                      ! EVB output file
         evbout_specified = .true.
         iarg = iarg + 1
         call getarg_wrap(iarg,evbout)
!AMD log file
      else if (arg == '-amd') then
         amdlog_specified = .true.
         iarg = iarg + 1
         call getarg_wrap(iarg, amdlog)
!scaledMD log file
      else if (arg == '-scaledMD') then
         scaledMDlog_specified = .true.
         iarg = iarg + 1
         call getarg_wrap(iarg, scaledMDlog)

#ifdef RISMSANDER
      else if (arg == '-xvv') then
         iarg = iarg + 1
         call getarg_wrap(iarg,xvvfile)
      else if (arg == '-guv') then
         iarg = iarg + 1
         call getarg_wrap(iarg,guvfile)
      else if (arg == '-huv') then
         iarg = iarg + 1
         call getarg_wrap(iarg,huvfile)
      else if (arg == '-cuv') then
         iarg = iarg + 1
         call getarg_wrap(iarg,cuvfile)
      else if (arg == '-uuv') then
         iarg = iarg + 1
         call getarg_wrap(iarg,uuvfile)
      else if (arg == '-asymp') then
         iarg = iarg + 1
         call getarg_wrap(iarg,asympfile)
      else if (arg == '-quv') then
         iarg = iarg + 1
         call getarg_wrap(iarg,quvfile)
      else if (arg == '-chgdist') then
         iarg = iarg + 1
         call getarg_wrap(iarg,chgdistfile)
      else if (arg == '-exchem') then
         iarg = iarg + 1
         call getarg_wrap(iarg,exchemfile)
      else if (arg == '-solvene') then
         iarg = iarg + 1
         call getarg_wrap(iarg,solvenefile)
      else if (arg == '-entropy') then
         iarg = iarg + 1
         call getarg_wrap(iarg,entropyfile)
      else if (arg == '-exchemGF') then
         iarg = iarg + 1
         call getarg_wrap(iarg,exchemGFfile)
      else if (arg == '-solveneGF') then
         iarg = iarg + 1
         call getarg_wrap(iarg,solveneGFfile)
      else if (arg == '-entropyGF') then
         iarg = iarg + 1
         call getarg_wrap(iarg,entropyGFfile)
      else if (arg == '-exchemUC') then
         iarg = iarg + 1
         call getarg_wrap(iarg,exchemUCfile)
      else if (arg == '-solveneUC') then
         iarg = iarg + 1
         call getarg_wrap(iarg,solveneUCfile)
      else if (arg == '-entropyUC') then
         iarg = iarg + 1
         call getarg_wrap(iarg,entropyUCfile)
      else if (arg == '-potUV') then
         iarg = iarg + 1
         call getarg_wrap(iarg,potUVfile)
      else if (arg == '-rismcrd') then
         iarg = iarg + 1
         call getarg_wrap(iarg,rismcrdfil)
      else if (arg == '-rismfrc') then
         iarg = iarg + 1
         call getarg_wrap(iarg,rismfrcfil)
      else if (arg == '-rismcrdrst') then
         iarg = iarg + 1
         call getarg_wrap(iarg,rismcrdrstfil)
      else if (arg == '-rismfrcrst') then
         iarg = iarg + 1
         call getarg_wrap(iarg,rismfrcrstfil)
#endif

#ifdef MPI
      else if (arg(1:3) == '-p4') then
         iarg = iarg+1
      else if (arg == '-np') then
         iarg = iarg+1
      else if (arg == '-mpedbg') then
         continue
      else if (arg == '-dbx') then
         continue
      else if (arg == '-gdb') then
         continue

      else if (arg == '-nrecip') then
         iarg = iarg + 1
         call getarg_wrap(iarg,arg)
         read(arg,'(i5)',err=91) num_recip
         if (num_recip == numtasks) then
            num_direct=numtasks
         else
            num_direct=numtasks-num_recip
         end if
         
         !     Parse input options for multisander
         
      else if (arg == '-ng') then
         iarg = iarg + 1
         call getarg_wrap(iarg,arg)
         read(arg,'(i5)',err=91) numgroup

      else if (arg == '-ng-nonsequential') then
         ng_sequential = .false.

      else if (arg == '-groupfile') then
         iarg = iarg + 1
         call getarg_wrap(iarg,groupfile_holder)

      else if (arg == '-gpes') then
         iarg = iarg + 1
         call getarg_wrap(iarg,gpes)

      else if (arg == '-rem') then
         iarg = iarg + 1
         call getarg_wrap(iarg, arg)
         read(arg, '(i5)', err=91) rem

         if (rem /= -1 .and. len_trim(remd_dimension_file) /= 0) then
            write(6, *) 'Error: -rem cannot be set when -remd-file is used!'
            call mexit(6, 1)
         end if

      else if (arg == '-nslice') then                      ! # of PIMD slices
         iarg = iarg + 1
         call getarg_wrap(iarg,arg)
         read(arg, '(i5)', err=91) nslice

      else if (arg == '-remlog') then
         iarg = iarg + 1
         call getarg_wrap(iarg, remlog)

      !     RREMD Options
      else if (arg == '-rremd') then
         iarg = iarg + 1
         call getarg_wrap(iarg,arg)
         read(arg, '(i5)',err=91) rremd

      else if (arg == '-saveene') then
         iarg = iarg + 1
         call getarg_wrap(iarg, saveenefile)

      else if (arg == '-clusterinfo') then
         iarg = iarg + 1
         call getarg_wrap(iarg, clusterinfofile)

      else if (arg == '-reservoir') then
         iarg = iarg + 1
         call getarg_wrap(iarg, reservoirname)
      !     End RREMD options
      else if (arg == '-remtype') then
         iarg = iarg + 1
         call getarg_wrap(iarg, remtype)

      else if (arg == '-hybridtraj') then
         iarg = iarg + 1
         call getarg_wrap(iarg, remstripcoord)

      else if (arg == '-remd-file') then
         iarg = iarg + 1
         call getarg_wrap(iarg, remd_dimension_file)
         if (rem > 0) then
            write(6, *) 'Error: -rem cannot be set with -remd-file!'
            call mexit(6, 1)
         end if
         rem = -1
      ! End REMD Options

#endif
#ifdef PUPIL_SUPPORT
      else if ((arg == '-ORBInitialPort') .or. &
               (arg == '-ORBInitialHost') .or. &
               (arg == '-jxms'          ) .or. &
               (arg == '-jxmx'          ) .or. &
               (arg == '-jxss'          ) .or. &
               (arg == '-OptPrint')) then
        iarg = iarg + 1
#endif /*PUPIL_SUPPORT*/
      else if (arg == '-pimdout') then
         iarg = iarg + 1
         call getarg_wrap(iarg,pimdout)
      
      else if (arg == '-suffix') then
         ! suffix_specified is already true for REMD. If we specify it INSIDE
         ! REMD, that means we want to add that file suffix to even files we
         ! specify
         if (suffix_specified) then
            always_add_suffix = .true.
         else
            suffix_specified = .true.
         end if
         iarg = iarg + 1
         call getarg_wrap(iarg, outfile_suffix)
      else if (arg == '-cns') then
         is_xref_on = .true.
      else if (arg == ' ') then
         continue
      else
         write(6,'(/,5x,a,a)') 'mdfil: Error unknown flag: ',arg
         write(6,9000)
         call mexit(6, 1)
      end if 
   end do  !  while (iarg < last_arg_index)
   
   ! Now it may be time to add the suffixes to all of the output files. File
   ! suffixes are added if:
   !
   ! 1. All unspecified files IF a groupfile is used OR -suffix is supplied
   !
   ! 2. ALL output files (regardless of if it was specified) if a groupfile is
   !    used AND -suffix is used on that line
   
   if (suffix_specified) then
      
      if (.not. mdout_specified .or. always_add_suffix) &
         call add_suffix(mdout, outfile_suffix)

      if (.not. restrt_specified .or. always_add_suffix) &
         call add_suffix(restrt, outfile_suffix)

      if (.not. mdinfo_specified .or. always_add_suffix) &
         call add_suffix(mdinfo, outfile_suffix)

      if (.not. mdcrd_specified .or. always_add_suffix) &
         call add_suffix(mdcrd, outfile_suffix)

      if (.not. mdvel_specified .or. always_add_suffix) &
         call add_suffix(mdvel, outfile_suffix)
      
      if (.not. mdfrc_specified .or. always_add_suffix) &
         call add_suffix(mdfrc, outfile_suffix)

      if (.not. mden_specified .or. always_add_suffix) &
         call add_suffix(mden, outfile_suffix)

      if (.not. cpout_specified .or. always_add_suffix) &
         call add_suffix(cpout, outfile_suffix)

      if (.not. cprestrt_specified .or. always_add_suffix) &
         call add_suffix(cprestrt, outfile_suffix)

      if (.not. amdlog_specified .or. always_add_suffix) &
         call add_suffix(amdlog, outfile_suffix)

      if (.not. scaledMDlog_specified .or. always_add_suffix) &
         call add_suffix(scaledMDlog, outfile_suffix)

      if (.not. cph_dump_specified .or. always_add_suffix) &
         call add_suffix(cph_dump, outfile_suffix)

   end if ! suffix_specified

   ! Load the groupfile_holder into the main groupfile
   groups = groupfile_holder

   return
   
#ifdef MPI
   91 write(6,*) 'mdfil: Error "-nrecip", "-rem" and "-ng"', &
                 'require integer arguments'
   write(6,*)'                   '
   call mexit(6, 1)
#endif
   9000 format(/,5x, &
         'usage: sander  [-O|A] -i mdin -o mdout -p prmtop -c inpcrd ', &
         '-r restrt',/19x,'[-ref refc -x mdcrd -v mdvel -e mden -frc mdfrc ', &
         '-idip inpdip -rdip rstdip -mdip mddip ', /19x, &
#ifdef MPI
         '-ng numgroup -remlog rem.log -remtype rem.type -rem [0|1|2] ', &
         '-rremd [0|1|2|3] -saveene saveene ',/19x, &
         '-clusterinfo cluster.info ', &
         '-reservoir reserv/frame -hybridtraj hybrid.strip.crd', /19x, &
#endif
         '-inf mdinfo -radii radii -y inptraj -amd amd.log -scaledMD scaledMD.log] -cph-data <file>' &
         , /, 'Consult the manual for additional options.')
end subroutine mdfil 

end module commandline_module

!        '-O                Overwrite existing files.',
!        '-A                Append existing files.',
!        '-i MDIN           Namelist control input file',
!        '-o MDOUT          Output file.',
!        '-p PARM           ParmTop file.',
!        '-c INPCRD         Coordinate file.',
!        '-vecs VECS        ???',
!        '-radii RADII      ???',
!        '-f FREQE          ???',
!        '-r RESTRT         ???',
!        '-ref REFC         ???',
!        '-z REFC           alias for -ref.',
!        '-e MDEN           ???',
!        '-v MDVEL          ???',
!        '-x MDCRD          ???',
!        '-t MDCRD          alias for -x MDCRD',
!        '-inf MDINFO       ???',
!        '-y INPCRD         input trajectory for imin==5',
!        '-idip INPDIP      ???',
!        '-rdip RSTDIP      ???',
!        '-mdip MDDIP       ???',
!        '-ng NUMGROUP      number of separate sander groups',
!        '-rem [0|1|2]      type of REM simulation',
!        '-remlog REMLOG    the filename for REM log file',
!        '-rremd [0|1|2]    type of REMD reservoir',
!        '-remtype REMTYPE  the filename for REM simulation type output file',
!        '-hybridtraj TRAJOUT    the filename for hybrid stripped output traj',
!        '-saveene SAVEENE  the filename for energies of reservoir structures',
!        '-clusterinfo CLUSTERINFO    the filename with dihedral cluster info',
!        '-reservoir RESNAME the filename of reservoir structures',
!        '-cpin CPDAT       Constant pH state information ',
!        '-cprstrt CPDAT    Constant pH state restart information',
!        '-cpout CPOUT      Constant pH protonation output
!        '-evbin  EVBIN     EVB input file ',
!        '-evbout EVBOUT    EVB output file ',
!        '-amd AMDLOG    the filename for AMD log file',
!        '-scaledMD scaledMDLOG    the filename for scaledMD log file',
!#ifdef MPI
!        '-nrecip N     Set number of reciprocal tasks to N;'
!        '              if < numtasks, remaining tasks go to direct.'
!        '-p4 -np -mpedbg -dbx -gdb -- flags read by MPI'
!#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Wrapper for IARGC to support reading command line input from a string
integer function iargc_wrap()

#ifdef MPI
   use file_io_dat
#endif

   implicit none
   integer iargc
#ifdef MPI
#  include "parallel.h"

   integer istart, iend
   integer ia, ie

   if (numgroup > 1 .and. groups(1:1) /= ' ') then
      ia = 0
      istart = 1
      iend = len(groupbuffer)

      do while (istart <= iend)
         if ( groupbuffer(istart:istart) == ' ' ) then
            istart = istart + 1
         else
            do ie = istart, iend
              if (groupbuffer(ie:ie) == ' ') exit
            end do
            ia = ia + 1
            istart = ie
         end if
      end do
      iargc_wrap = ia
   else
#endif

      iargc_wrap = iargc()

#ifdef MPI
   end if
#endif
end function iargc_wrap 



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Wrapper for GETARG to support grabbing command line arguments from a string
subroutine getarg_wrap(iarg, arg)

#ifdef MPI
   use file_io_dat
#endif
   implicit none
   integer iarg
   character(len=*) arg

#ifdef MPI
   integer istart, iend
   integer ia, ie

   if (groups(1:1) /= ' ') then

      ia = 0
      istart = 1
      iend = len(groupbuffer)

      do while (istart <= iend)
         if ( groupbuffer(istart:istart) == ' ' ) then
            istart = istart + 1
         else
            do ie = istart, iend
              if (groupbuffer(ie:ie) == ' ') exit
            end do
            ia = ia + 1

            if (iarg == ia) then
               arg = groupbuffer(istart:ie)
               return
            end if
            istart = ie
         end if
      end do

   else
#endif

      ! Intrinsic getarg requires a 4 byte integer argument; 
      ! this guards the argument for builds with default 8 byte integers.
      call getarg(int(iarg,4), arg)

#ifdef MPI
   end if
#endif
end subroutine getarg_wrap 

subroutine add_suffix(filename, suffix)

   use file_io_dat

   implicit none

   character (len=MAX_FN_LEN), intent(in out) :: filename
   character (len=MAX_FN_LEN), intent(in)     :: suffix

   if (len_trim(filename) + len_trim(suffix) .gt. MAX_FN_LEN) then
      write(0,*) 'Filename + suffix too long!'
      call mexit(6, 1)
   else
      filename = trim(filename) // '.' // trim(suffix)
   end if

end subroutine add_suffix

subroutine basename(filename)

   implicit none

   ! Passed variables
   character(*), intent(in out)  :: filename

   ! Local variables
   integer :: i, pos

   pos = 1
   do i = 1, len_trim(filename)
      if (filename(i:i) .eq. '/') &
         pos = i + 1
   end do

   filename = filename(pos:len_trim(filename))

   return

end subroutine basename
