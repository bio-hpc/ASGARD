!<compile=optimized>

!The 3D-RISM-KH software found here is copyright (c) 2011-2012 by Tyler Luchko 
!and David A. Case.
!
!This program is free software: you can redistribute it and/or modify it
!under the terms of the GNU General Public License as published by the Free
!Software Foundation, either version 3 of the License, or (at your option)
!any later version.
!
!This program is distributed in the hope that it will be useful, but
!WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
!or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!for more details.
!
!You should have received a copy of the GNU General Public License in the
!../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
!
!Users of the 3D-RISM capability found here are requested to acknowledge
!use of the software in reports and publications.  Such acknowledgement
!should include the following citations:
!
!1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999); 
!ibid. 112:10391-10417 (2000).   
!
!2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by  
!F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.  
!
!3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
!and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010). 

#include "../include/dprec.fh"

!This module is used to read command line arguments with a simple interface.
!The module is used in three steps:
!
!1)Initialize (add variables)
!For each key-value the user can set on the commandline, used the getopt_add
!subroutine to add the keyname and default value(s).  E.g.,
!
!  call getopt_add("longkey", default, "shortkey", maxvalues, minvalues)
!
!shortkey, minvalues and maxvlues are all optional.
!longkey  - multicharacter key
!default  - default value
!shortkey - single character key
!minvalues - minimum number of values the user must provide (default zero)
!maxvalues - maximum number of values the user can provide (default huge)
!
!Will allow the user to use on the command line 
!
!  executable --longkey value1 -shortkey value2 value3
!
!Of course, 'shortkey' is really a single character.  'value' can be
!of type character, integer, _REAL_ or logical.  The user can specify
!multiple values for each key (all will be stored).  Logicals can only
!have one value each.  The user can use -nokey to set it to .false. or
!-key to set it to .true..
!
!2)Process commandline
!
!  call getopt_process()
!
!3)Get values
!Values for specific keys can be obtained individually or as an array. If 
!nothing has been set then the default value(s) are used. E.g., get the first 
!value:
!
!  call getopt_get("key",1, value)
!
!where value is a scalar.  For the full array this would be
!
!  value => getopt_get("key",value)
!
!where value is a pointer of the appropriate type.
!
!3b)Get leftover values
!If not all command line arguments map to specific keys (such as input/output
!file names) they will be leftover at the end.  These can be obtained as a character
!pointer.  E.g.
!
!  extraArgs => getopt_unread()
!
!KNOWN ISSUES:
!-Not fully POSIX compliant: http://www.opengroup.org/onlinepubs/009695399/basedefs/xbd_chap12.html#tag_12_02
!  +Guideline 3:  multicharacter options are allowed.  '-W' is not
!                 reserved.
!  +Guideline 4:  multicharacter options are proceeded by '--'
!  +Guideline 5:  option grouping is not supported
!  +Guideline 8:  only space separated multi-argument options allowed.  
!  +Guildline 9:  operands may appear anywhere.  However, as operands are generally
!                 strings, they should not follow a string option
module getopts_c
  use safemem
  use array_util, only : array_index
  use rism_report_c
  implicit none
  
  !Error enumeration
  integer, public, parameter :: getopts_OK=0, getopts_TOOFEW=-1, getopts_TOOMANY=-2

  integer, private,  parameter :: CLEN =1024
  !store character keys
  character(len=CLEN), private, pointer :: clongkey(:)=>NULL()
  !store single character keys
  character, private, pointer :: cshortkey(:)=>NULL()
  !store character values
  character(len=CLEN), private, pointer :: cval(:,:)=>NULL()
  !store character default values
  character(len=CLEN), private, pointer :: cdef(:,:)=>NULL()
  !minmum and maximum values for each character key
  integer, private, pointer :: cminmax(:,:)=>NULL()
  !store integer keys
  character(len=CLEN), private, pointer :: ilongkey(:)=>NULL()
  !store single character integer keys
  character(len=CLEN), private, pointer :: ishortkey(:)=>NULL()
  !store integer values
  integer, private, pointer :: ival(:,:)=>NULL()
  !store integer default values
  integer, private, pointer :: idef(:,:)=>NULL()
  !minmum and maximum values for each integer key
  integer, private, pointer :: iminmax(:,:)=>NULL()
  !store real keys
  character(len=CLEN), private, pointer :: rlongkey(:)=>NULL()
  !store single character real keys
  character(len=CLEN), private, pointer :: rshortkey(:)=>NULL()
  !store real values
  _REAL_, private, pointer :: rval(:,:)=>NULL()
  !store real default values
  _REAL_, private, pointer :: rdef(:,:)=>NULL()
  !minmum and maximum values for each real key
  integer, private, pointer :: rminmax(:,:)=>NULL()
  !store logical keys
  character(len=CLEN), private, pointer :: llongkey(:)=>NULL()
  !store single character logical keys
  character(len=CLEN), private, pointer :: lshortkey(:)=>NULL()
  !store logical values
  logical, private, pointer :: lval(:)=>NULL()

  !list of command line argument indices that have been read
  integer, private, pointer :: read_arg(:)=>NULL()

  !list of all command line arguments
  character(len=CLEN), private, pointer :: args(:)=>NULL()
  !number of command line argumens
  integer, private :: nargs=0

  public getopts_add, getopts_process, getopts_get, getopts_narg, getopts_unread, &
       getopts_cleanup

  private add_c, add_i, add_r, add_l,&
       findkey, read_c, read_i, read_r, read_l, sanitycheck, &
       get_c, get_i, get_r, get_l, get_c_pt, get_i_pt, get_r_pt

  interface getopts_add
     module procedure add_c, add_i, add_r, add_l
  end interface

  interface getopts_get
     module procedure get_c, get_i, get_r, get_l
  end interface

  interface getopts_getAll
     module procedure get_c_pt, get_i_pt, get_r_pt
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Adds a character key and default value. Optionally, adds a single character
!!!key and min and max number of values.
!!!IN:
!!!   longkey : name of the key and argument as it will appear on the command line  
!!!   default : default value associated with this key  
!!!   shortkey : (optional) single character name of the key and argument as it
!!!              will appear on the command line  
!!!   min : (optional) minmum number of values
!!!   max : (optional) maximum number of values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine add_c(longkey,default,shortkey,max,min)
    use safemem
    implicit none
    character(len=*), intent(in) :: longkey,default
    character, optional, intent(in) :: shortkey
    integer, optional, intent(in) :: min,max
    clongkey   => safemem_realloc(clongkey,len(clongkey),ubound(clongkey,1)+1)
    cshortkey   => safemem_realloc(cshortkey,len(cshortkey),ubound(cshortkey,1)+1)
    cdef    => safemem_realloc(cdef,len(cdef),ubound(cdef,1)+1,1)
    cminmax => safemem_realloc(cminmax,ubound(cminmax,1)+1,2)
    clongkey(ubound(clongkey,1)) = longkey
    cdef(ubound(cdef,1),1) = default
    if(present(shortkey))then
       cshortkey(ubound(cshortkey,1)) = shortkey
    else
       cshortkey(ubound(cshortkey,1)) = ''
    end if
    if(present(min))then
       cminmax(ubound(cminmax,1),1) = min
    else
       cminmax(ubound(cminmax,1),1) = 0
    end if
    if(present(max))then
       cminmax(ubound(cminmax,1),2) = max
    else
       cminmax(ubound(cminmax,1),2) = huge(1)
    end if
  end subroutine add_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Adds an integer key and default value. Optionally, adds a single character
!!!key and min and max number of values.
!!!IN:
!!!   longkey : name of the key and argument as it will appear on the command line  
!!!   default : default value associated with this key  
!!!   shortkey : (optional) single character name of the key and argument as it
!!!              will appear on the command line  
!!!   min : (optional) minmum number of values
!!!   max : (optional) maximum number of values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine add_i(longkey,default,shortkey,max,min)
    use safemem
    implicit none
    character(len=*), intent(in) :: longkey
    integer, intent(in) :: default
    character, optional, intent(in) :: shortkey
    integer, optional, intent(in) :: min,max
    ilongkey => safemem_realloc(ilongkey,len(ilongkey),ubound(ilongkey,1)+1)
    ishortkey   => safemem_realloc(ishortkey,len(ishortkey),ubound(ishortkey,1)+1)
    idef  => safemem_realloc(idef,ubound(idef,1)+1,1)
    iminmax => safemem_realloc(iminmax,ubound(iminmax,1)+1,2)
    ilongkey(ubound(ilongkey,1)) = longkey
    idef(ubound(idef,1),1) = default
    if(present(shortkey))then
       ishortkey(ubound(ishortkey,1)) = shortkey
    else
       ishortkey(ubound(ishortkey,1)) = ''
    end if
    if(present(min))then
       iminmax(ubound(iminmax,1),1) = min
    else
       iminmax(ubound(iminmax,1),1) = 0
    end if
    if(present(max))then
       iminmax(ubound(iminmax,1),2) = max
    else
       iminmax(ubound(iminmax,1),2) = huge(1)
    end if
  end subroutine add_i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Adds an _REAL_ key and default value. Optionally, adds a single character
!!!key and min and max number of values.
!!!IN:
!!!   longkey : name of the key and argument as it will appear on the command line  
!!!   default : default value associated with this key  
!!!   shortkey : (optional) single character name of the key and argument as it
!!!              will appear on the command line  
!!!   min : (optional) minmum number of values
!!!   max : (optional) maximum number of values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine add_r(longkey,default,shortkey,max,min)
    use safemem
    implicit none
    character(len=*), intent(in) :: longkey
    _REAL_, intent(in) :: default
    character, optional, intent(in) :: shortkey
    integer, optional, intent(in) :: min,max
    rlongkey => safemem_realloc(rlongkey,len(rlongkey),ubound(rlongkey,1)+1)
    rshortkey => safemem_realloc(rshortkey,len(rshortkey),ubound(rshortkey,1)+1)
    rdef  => safemem_realloc(rdef,ubound(rdef,1)+1,1)
    rminmax => safemem_realloc(rminmax,ubound(rminmax,1)+1,2)
    rlongkey(ubound(rlongkey,1)) = longkey
    rdef(ubound(rdef,1),1) = default
    if(present(shortkey))then
       rshortkey(ubound(rshortkey,1)) = shortkey
    else
       rshortkey(ubound(rshortkey,1)) = ''
    end if
    if(present(min))then
       rminmax(ubound(rminmax,1),1) = min
    else
       rminmax(ubound(rminmax,1),1) = 0
    end if
    if(present(max))then
       rminmax(ubound(rminmax,1),2) = max
    else
       rminmax(ubound(rminmax,1),2) = huge(1)
    end if
  end subroutine add_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Adds an logical key and default value. Optionally, adds a single character
!!!key.
!!!IN:
!!!   longkey : name of the key and argument as it will appear on the command line  
!!!   default : default value associated with this key  
!!!   shortkey : (optional) single character name of the key and argument as it
!!!              will appear on the command line  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine add_l(longkey,default,shortkey)
    use safemem
    implicit none
    character(len=*), intent(in) :: longkey
    logical, intent(in) :: default
    character, optional, intent(in) :: shortkey
    llongkey => safemem_realloc(llongkey,len(llongkey),ubound(llongkey,1)+1)
    lshortkey   => safemem_realloc(lshortkey,len(lshortkey),ubound(lshortkey,1)+1)
    lval  => safemem_realloc(lval,ubound(lval,1)+1)
    llongkey(ubound(llongkey,1)) = longkey
    lval(ubound(lval,1)) = default
    if(present(shortkey))then
       lshortkey(ubound(lshortkey,1)) = shortkey
    else
       lshortkey(ubound(lshortkey,1)) = ''
    end if
  end subroutine add_l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Processes command line arguments using specified keys
!!!OUT:
!!!    number of tokens on the command line if everything is ok.  Otherwise an 
!!!    error number is returned
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function getopts_process() result(error)
    implicit none
    integer :: i, error,j, iarg, iargc
    error = getopts_OK
    !fill arguments array
    !NOTE: this is the only place that iargc() and getarg() are used.  Upgrade to
    !      F2003 here.
    args => safemem_realloc(args,len(args),0)
    read_arg => safemem_realloc(read_arg,0)
    do iarg=1, iargc()
       args => safemem_realloc(args,len(args),iarg)
       call getarg(iarg,args(iarg))
    end do
    nargs = ubound(args,1)

    !assign character strings
    do i=1,ubound(clongkey,1)
       iarg=1
       do while(iarg < nargs)
          if(findkey(clongkey(i),cshortkey(i),iarg))then
             call read_c(i,iarg)
          else
             exit
          end if
       end do
    end do

    !assign integers
    do i=1,ubound(ilongkey,1)
       iarg=1
       do while(iarg < nargs)
          if(findkey(ilongkey(i),ishortkey(i),iarg))then
             call read_i(i,iarg)
          else
             exit
          end if
       end do
    end do

    !assign _REAL_s
    do i=1,ubound(rlongkey,1)
       iarg=1
       do while(iarg < nargs)
          if(findkey(rlongkey(i),rshortkey(i),iarg))then
             call read_r(i,iarg)
          else
             exit
          end if
       end do
    end do
    
    !assign logicals
    do i=1,ubound(llongkey,1)
       iarg=1
       do while(iarg <= nargs)
          if(findkey(llongkey(i),lshortkey(i),iarg,.true.))then
             call read_l(i,iarg)
          else
             exit
          end if
       end do
    end do
    error = sanitycheck()
    if(error == getopts_OK) error = nargs
  end function getopts_process

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns an allocated array of unread arguments.  These are in the order they 
!!!appear in the command line.
!!!OUT:
!!!    Pointer array of unprocessed arguments.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function getopts_unread() result(argv)
    implicit none
    character(len=CLEN), pointer :: argv(:)
    integer :: i
    nullify(argv)
    argv => safemem_realloc(argv,CLEN,0)
    do i = 1, nargs
       if(array_index(read_arg,i) >0)cycle
       if(trim(args(i)) .eq. "--") cycle
       argv => safemem_realloc(argv,CLEN,ubound(argv,1)+1)
       argv(ubound(argv,1)) = trim(args(i))
    end do
  end function getopts_unread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Writes a summary of read in options to the specified unit.
!!!IN:
!!!   unit : location to print to
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getopts_summary(unit)
    implicit none
    integer, intent(in) :: unit
    integer ::  i,j
    write(unit,'(a)') "Character Options:"
    write(unit,'(a)') "------------------"
    do i=1,ubound(clongkey,1)
       write(unit,'(a)',advance="no") trim(clongkey(i))//","//cshortkey(i)//":"
       do j=1,getopts_narg(clongkey(i))
          write(unit,'(a)',advance="no") " "//trim(cval(i,j))
       end do
       write(unit,'(a)')
    end do
    write(unit,'(a)') "Integer Options:"
    write(unit,'(a)') "------------------"
    do i=1,ubound(ilongkey,1)
       write(unit,'(a)',advance="no") trim(ilongkey(i))//","//ishortkey(i)//":"
       do j=1,getopts_narg(ilongkey(i))
          write(unit,'(a,i8)',advance="no") " ",ival(i,j)
       end do
       write(unit,'(a)')
    end do
    write(unit,'(a)') "_REAL_ Options:"
    write(unit,'(a)') "------------------"
    do i=1,ubound(rlongkey,1)
       write(unit,'(a)',advance="no") trim(rlongkey(i))//","//rshortkey(i)//":"
       do j=1,getopts_narg(rlongkey(i))
          write(unit,'(a,g12.4)',advance="no") " ",rval(i,j)
       end do
       write(unit,'(a)')
    end do
    write(unit,'(a)') "Logical Options:"
    write(unit,'(a)') "------------------"
    do i=1,ubound(llongkey,1)
       write(unit,'(a)',advance="no") trim(llongkey(i))//","//lshortkey(i)//":"
       write(unit,'(a,l)',advance="no") " ",lval(i)
       write(unit,'(a)')
    end do
  end subroutine getopts_summary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the i'th argument of the given string option
!!!IN:
!!!   key : name of the option (must exist)
!!!   i   : index number.
!!!   value : the value at this index.  If the key or index does not exist, an empty 
!!!           string is returned.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_c(key,i, value)
    implicit none
    character(len=*), intent(in) :: key
    integer, intent(in) :: i
    character(len=*) :: value
    integer :: ikey, nval
    value = ''
    ikey = array_index(clongkey,key)
    if(ikey < 0)then
       ikey = array_index(cshortkey,key)
       if(ikey < 0) return
    end if
    if(associated(cval))then
       nval =  array_index(cval(ikey,:),'')-1
       if(nval==-2) nval = ubound(cval,2)
    else
       nval = 0
    end if
    if(i > nval)then
       !if i==1 and no option has been set, use default
       if(i==1) value = trim(cdef(ikey,i))
       return
    end if
    value = trim(cval(ikey,i))
  end subroutine get_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the i'th argument of the given integer option
!!!IN:
!!!   key : name of the option (must exist)
!!!   i   : index number.
!!!   value : the value at this index.  If the key or index does not exist, huge(1) is
!!!           returned
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_i(key,i, value)
    implicit none
    character(len=*), intent(in) :: key
    integer, intent(in) :: i
    integer :: value
    integer :: ikey,nval
    value=huge(value)
    ikey = array_index(ilongkey,key)
    if(ikey < 0)then
       ikey = array_index(ishortkey,key)
       if(ikey < 0) return
    end if
    if(associated(ival))then
       nval =  array_index(ival(ikey,:),huge(1))-1
       if(nval==-2) nval = ubound(ival,2)
    else
       nval = 0
    end if
    if(i > nval)then
       !if i==1 and no option has been set, use default
       if(i==1) value = idef(ikey,i)
       return
    end if
    value = ival(ikey,i)
  end subroutine get_i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the i'th argument of the given real option
!!!IN:
!!!   key : name of the option (must exist)
!!!   i   : index number.
!!!   value : the value at this index.  If the key or index does not exist, 
!!!           huge(1d0) or huge(1e0) is returned
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_r(key,i, value)
    implicit none
    character(len=*), intent(in) :: key
    integer, intent(in) :: i
    _REAL_ :: value
    integer :: ikey, nval
    value=huge(value)
    ikey = array_index(rlongkey,key)
    if(ikey < 0)then
       ikey = array_index(rshortkey,key)
       if(ikey < 0) return
    end if
    if(associated(rval))then
       nval =  array_index(rval(ikey,:),huge(1d0))-1
       if(nval==-2) nval = ubound(rval,2)
    else
       nval = 0
    end if
    if(i > nval)then
       !if i==1 and no option has been set, use default
       if(i==1) value = rdef(ikey,i)
       return
    end if
    value = rval(ikey,i)
  end subroutine get_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns value of the given logical option
!!!IN:
!!!   key : name of the option (must exist)
!!!   value : the value at this index.  If the key does not exist, 
!!!           .false. is returned
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_l(key, value)
    implicit none
    character(len=*), intent(in) :: key
    logical :: value
    integer :: ikey
    value=.false.
    ikey = array_index(llongkey,key)
    if(ikey < 0)then
       ikey = array_index(lshortkey,key)
       if(ikey < 0) return
    end if
    value = lval(ikey)
  end subroutine get_l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the full array of arguments for the given string option.  The array is
!!!allocated and the user is responsible for returning the memory
!!!IN:
!!!   key : name of the option (must exist)
!!!   type : the pointer destination.  Should be the same as value
!!!OUT:
!!!   value : the allocated pointer of values.  If the key does not exist, a 
!!!           NULL() pointer is return.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function get_c_pt(key, type,length)! result(value)
    implicit none
    character(len=*), intent(in) :: key
    character(len=length),pointer,dimension(:) :: type
    character(len=length),pointer,dimension(:) :: get_c_pt
    integer, intent(in) :: length
    integer :: ikey, nval
    get_c_pt=>NULL()
    ikey = array_index(clongkey,key)
    if(ikey < 0)then
       ikey = array_index(cshortkey,key)
       if(ikey < 0)then
          if(associated(get_c_pt)) deallocate(get_c_pt)
          get_c_pt=>NULL()
          return
       end if
    end if
    if(.not.associated(cval))then
       get_c_pt=> safemem_realloc(get_c_pt,length,ubound(cdef,2))
       get_c_pt= cdef(ikey,:)
       return
    end if
    nval = array_index(cval(ikey,:),'')-1
    if(nval<0) nval = ubound(cval,2)
    get_c_pt=> safemem_realloc(get_c_pt,length,nval)
    get_c_pt= cval(ikey,:nval)
  end function get_c_pt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the full array of arguments for the given integer option.  The array is
!!!allocated and the user is responsible for returning the memory
!!!IN:
!!!   key : name of the option (must exist)
!!!   type : the pointer destination.  Should be the same as value
!!!OUT:
!!!   value : the allocated pointer of values.  If the key does not exist, a 
!!!           NULL() pointer is return.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function get_i_pt(key, type)! result(value)
    implicit none
    character(len=*), intent(in) :: key
    integer,pointer,dimension(:) :: type
    integer,pointer,dimension(:) :: get_i_pt
    integer :: ikey,nval
    get_i_pt=>NULL()
    ikey = array_index(ilongkey,key)
    if(ikey < 0)then
       ikey = array_index(ishortkey,key)
       if(ikey < 0)then
          if(associated(get_i_pt)) deallocate(get_i_pt)
          get_i_pt=>NULL()
          return
       end if
    end if
    if(.not.associated(ival))then
       get_i_pt=> safemem_realloc(get_i_pt,ubound(idef,2))
       get_i_pt= idef(ikey,:)
       return
    end if
    nval = array_index(ival(ikey,:),huge(1))-1
    if(nval<0) nval = ubound(ival,2)
    get_i_pt=> safemem_realloc(get_i_pt,nval)
    get_i_pt= ival(ikey,:nval)
  end function get_i_pt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the full array of arguments for the given _REAL_ option.  The array is
!!!allocated and the user is responsible for returning the memory
!!!IN:
!!!   key : name of the option (must exist)
!!!   type : the pointer destination.  Should be the same as value
!!!OUT:
!!!   value : the allocated pointer of values.  If the key does not exist, a 
!!!           NULL() pointer is return.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function get_r_pt(key, type)! result(value)
    implicit none
    character(len=*), intent(in) :: key
    _REAL_,pointer,dimension(:) :: type
    _REAL_,pointer,dimension(:) :: get_r_pt
    integer :: ikey,nval
    get_r_pt=>NULL()
    ikey = array_index(rlongkey,key)
    if(ikey < 0)then
       ikey = array_index(rshortkey,key)
       if(ikey < 0)then
          if(associated(get_r_pt)) deallocate(get_r_pt)
          get_r_pt=>NULL()
          return
       end if
    end if
    if(.not.associated(rval))then
       get_r_pt=> safemem_realloc(get_r_pt,ubound(rdef,2))
       get_r_pt= rdef(ikey,:)
       return
    end if
    nval = array_index(rval(ikey,:),huge(1d0))-1
    if(nval<0) nval = ubound(rval,2)
    get_r_pt=> safemem_realloc(get_r_pt,nval)
    get_r_pt= rval(ikey,:nval)
  end function get_r_pt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the numbers of arguments supplied for this key.  Does not include 
!!!defaults.
!!!IN:
!!!   key : name of the key
!!!OUT:
!!!    number of values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function getopts_narg(key) result(narg)
    implicit none
    character(len=*), intent(in) :: key
    integer :: narg, ikey
    narg=0
    !character key
    ikey = array_index(clongkey,key)
    if(ikey>0)then
       if(associated(cval))then
          narg = array_index(cval(ikey,:),'')-1
          if(narg <-1)then
             narg = ubound(cval,2)
          end if
       end if
       return
    end if

    !integer key
    ikey = array_index(ilongkey,key)
    if(ikey>0)then
       if(associated(ival))then
          narg = array_index(ival(ikey,:),huge(1))-1
          if(narg <-1)then
             narg = ubound(ival,2)
          end if
       end if
       return
    end if

    !_REAL_ key
    ikey = array_index(rlongkey,key)
    if(ikey>0)then
       if(associated(rval))then
          narg = array_index(rval(ikey,:),huge(rval(1,1)))-1
          if(narg <-1)then
             narg = ubound(rval,2)
          end if
       end if
       return
    end if
  end function getopts_narg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!frees all associated memory
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine getopts_cleanup
    implicit none
    if(safemem_dealloc(clongkey)/=0) call rism_report_error("GETOPTS: Failed to deallocate clongkey")
    if(safemem_dealloc(cshortkey)/=0) call rism_report_error("GETOPTS: Failed to deallocate cshortkey")
    if(safemem_dealloc(cval)/=0) call rism_report_error("GETOPTS: Failed to deallocate cval")
    if(safemem_dealloc(cdef)/=0) call rism_report_error("GETOPTS: Failed to deallocate cdef")
    if(safemem_dealloc(cminmax)/=0) call rism_report_error("GETOPTS: Failed to deallocate cminmax")

    if(safemem_dealloc(ilongkey)/=0) call rism_report_error("GETOPTS: Failed to deallocate ilongkey")
    if(safemem_dealloc(ishortkey)/=0) call rism_report_error("GETOPTS: Failed to deallocate ishortkey")
    if(safemem_dealloc(ival)/=0) call rism_report_error("GETOPTS: Failed to deallocate ival")
    if(safemem_dealloc(idef)/=0) call rism_report_error("GETOPTS: Failed to deallocate idef")
    if(safemem_dealloc(iminmax)/=0) call rism_report_error("GETOPTS: Failed to deallocate iminmax")

    if(safemem_dealloc(rlongkey)/=0) call rism_report_error("GETOPTS: Failed to deallocate rlongkey")
    if(safemem_dealloc(rshortkey)/=0) call rism_report_error("GETOPTS: Failed to deallocate rshortkey")
    if(safemem_dealloc(rval)/=0) call rism_report_error("GETOPTS: Failed to deallocate rval")
    if(safemem_dealloc(rdef)/=0) call rism_report_error("GETOPTS: Failed to deallocate rdef")
    if(safemem_dealloc(rminmax)/=0) call rism_report_error("GETOPTS: Failed to deallocate rminmax")

    if(safemem_dealloc(llongkey)/=0) call rism_report_error("GETOPTS: Failed to deallocate llongkey")
    if(safemem_dealloc(lshortkey)/=0) call rism_report_error("GETOPTS: Failed to deallocate lshortkey")
    if(safemem_dealloc(lval)/=0) call rism_report_error("GETOPTS: Failed to deallocate lval")

    if(safemem_dealloc(read_arg)/=0) call rism_report_error("GETOPTS: Failed to deallocate read_arg")
    if(safemem_dealloc(args)/=0) call rism_report_error("GETOPTS: Failed to deallocate args")
    nargs=0
  end subroutine getopts_cleanup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Searches the command line arguments for the give key.  If found, iarg is set 
!!!to that index and .true. is returned.  Otherwise .false. is returned and iarg
!!!is left untouched.
!!!IN:
!!!   key  : option we are searching for
!!!   iarg : position in list of all arguments
!!!   bool : (OPTIONAL) the key is a logical (boolean) value.  Accept 'no' as a 
!!!          option prefix to indicate .false.
!!!OUT:
!!!    .true. if the key is found. .false. otherwise.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function findkey(key,shortkey,iarg,bool) result(found)
    implicit none
    character(len=*), intent(in) :: key,shortkey
    integer, intent(inout) :: iarg
    logical, optional, intent(in) :: bool
    logical :: found, allow_no 
    integer :: iargOld
    found = .false.
    allow_no = .false.
    iargOld = iarg
    if(present(bool)) allow_no = bool
    do iarg = iarg, nargs
       if("--"//key .eq. args(iarg) .or.(allow_no .and. "--no"//key .eq. args(iarg)) &
          .or. "-"//shortkey .eq. args(iarg) .or.(allow_no .and. "-no"//shortkey .eq. args(iarg))  )then
          found = .true.
          read_arg => safemem_realloc(read_arg,ubound(read_arg,1)+1)
          read_arg(ubound(read_arg,1)) = iarg
          return
       end if
    end do
    iarg = iargOld
  end function findkey

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!checks if this is a registered key.
!!!IN:
!!!   key : key name. long or short
!!!OUT:
!!!   .true. if it is a key
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function iskey(key)
    implicit none
    character(len=*), intent(in) :: key
    logical :: iskey
    iskey = .false.

    !character key
    if(array_index(cshortkey,key) >= 1 .or. array_index(clongkey,key) >= 1)then
       iskey = .true.
       return
    end if

    !integer key
    if(array_index(ishortkey,key) >= 1 .or. array_index(ilongkey,key) >= 1)then
       iskey = .true.
       return
    end if

    !_REAL_ key
    if(array_index(rshortkey,key) >= 1 .or. array_index(rlongkey,key) >= 1)then
       iskey = .true.
       return
    end if

    !logical key
    if(array_index(lshortkey,key) >= 1 .or. array_index(llongkey,key) >= 1)then
       iskey = .true.
       return
    end if
  end function iskey

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Walks through the command line arguments, updating iarg along the way, and
!!!adds values to cval, reallocating memory as necessary.  The subroutine returns
!!!when it finds a value that can be considered a key.  It 
!!!is assumed that no values have been added for this key.  iarg points to the 
!!!next unread argument.
!!!IN:
!!!   i : index of cval to insert values.  I.e., corresponds to the correct key.
!!!   iarg : position in list of all arguments to start from
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_c(i,iarg)
    implicit none
    integer, intent(in) :: i
    integer, intent(inout) :: iarg
    integer :: iostat, itest, index
    _REAL_ :: rtest
    do index = 1,ubound(cval,2)
       if(len_trim(cval(i,index)) == 0) exit
    end do
    index = index-1
    do while(iarg < nargs)
       iarg = iarg+1
       index = index+1
       !don't bother checking if this is a 'string'.  Even numbers can be strings.
!!$       !check if this is an integer
!!$       read(args(iarg),*,iostat=iostat) itest
!!$       if(iostat == 0) return
!!$
!!$       !check if this ia a _REAL_
!!$       read(args(iarg),*,iostat=iostat) rtest
!!$       write(0,*) "REAL ", trim(args(iarg)), iostat
!!$       if(iostat == 0) then
!!$          !file names with some characters can be confused for reals
!!$          !Try looking non numeric characters instead
!!$          iostat = scan(args(iarg),"AaBbCcFfHhIiJjKkLlMmNnOoPpQqRrSsTtVvUuWwXxYyZz!@#$%^&*()_={}|\:;'<,>?/`~")
!!$          iostat = iostat + scan(args(iarg),'"')
!!$          if(iostat==0)return
!!$       end if

       !check if this is a key
       if(args(iarg)(1:2) .eq. '--')then
          if(iskey(trim(args(iarg)(3:)))) &
               return
       elseif(args(iarg)(1:1) .eq. '-')then
          if(iskey(trim(args(iarg)(2:)))) &
               return
       end if
       

       !this is probably a string.  
       !Allocate memory if necessary
       if(ubound(cval,2) < index)then
          cval => safemem_realloc(cval,CLEN,ubound(clongkey,1),ubound(cval,2)+1)
          cval(:,ubound(cval,2)) = ''
       end if
       !add value
       cval(i,index) = trim(args(iarg))
       !mark argument as read
       read_arg => safemem_realloc(read_arg,ubound(read_arg,1)+1)
       read_arg(ubound(read_arg,1)) = iarg
    end do
  end subroutine read_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Walks through the command line arguments, updating iarg along the way, and
!!!adds values to ival, reallocating memory as necessary.  The subroutine returns
!!!when it finds a value that in not an integer.  It is assumed that no values 
!!!have been added for this key. iarg points to the next unread argument.
!!!IN:
!!!   i : index of ival to insert values.  I.e., corresponds to the correct key.
!!!   iarg : position in list of all arguments to start from
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_i(i,iarg)
    implicit none
    integer, intent(in) :: i
    integer, intent(inout) :: iarg
    integer :: iostat, itest, index
    do index = 1,ubound(ival,2)
       if(ival(i,index) == huge(1)) exit
    end do
    index = index-1
    do while(iarg < nargs)
       iarg = iarg+1
       index = index+1
       !check if this is an integer
       read(args(iarg),*,iostat=iostat) itest
       if(iostat /= 0) return

       !Allocate memory if necessary
       if(ubound(ival,2) < index)then
          ival => safemem_realloc(ival,ubound(ilongkey,1),ubound(ival,2)+1)
          ival(:,ubound(ival,2)) = huge(1)
       end if
       !add value
       ival(i,index) = itest
       !mark argument as read
       read_arg => safemem_realloc(read_arg,ubound(read_arg,1)+1)
       read_arg(ubound(read_arg,1)) = iarg
    end do
  end subroutine read_i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Walks through the command line arguments, updating iarg along the way, and
!!!adds values to rval, reallocating memory as necessary.  The subroutine returns
!!!when it finds a value that in not an integer.  It is assumed that no values 
!!!have been added for this key. iarg points to the next unread argument.
!!!IN:
!!!   i : index of rval to insert values.  I.e., corresponds to the correct key.
!!!   iarg : position in list of all arguments to start from
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_r(i,iarg)
    implicit none
    integer, intent(in) :: i
    integer, intent(inout) :: iarg
    integer :: iostat, index
    _REAL_ :: rtest
    do index = 1,ubound(rval,2)
       if(rval(i,index) == huge(1d0)) exit
    end do
    index = index-1
    do while(iarg < nargs)
       iarg = iarg+1
       index = index+1
       !check if this is a _REAL_
       read(args(iarg),*,iostat=iostat) rtest
       if(iostat /= 0) return

       !Allocate memory if necessary
       if(ubound(rval,2) < index)then
          rval => safemem_realloc(rval,ubound(rlongkey,1),ubound(rval,2)+1)
          rval(:,ubound(rval,2)) = huge(1d0)
       end if
       !add value
       rval(i,index) = rtest
       !mark argument as read
       read_arg => safemem_realloc(read_arg,ubound(read_arg,1)+1)
       read_arg(ubound(read_arg,1)) = iarg
    end do
  end subroutine read_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Reads the logical key and sets the appropriate index in lval.  Assumes that
!!!iarg is pointing at the key and iarg has been labeled as read.  iarg points 
!!!to the next unread argument.
!!!IN:
!!!   i : index of lval to modify.  I.e., corresponds to the correct key.
!!!   iarg : position in list of all arguments to start from
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_l(i,iarg)
    implicit none
    integer, intent(in) :: i
    integer, intent(inout) :: iarg
    integer :: iostat
    lval(i) = args(iarg)(2:3) .ne. 'no'
    iarg = iarg+1
  end subroutine read_l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Check that everything is ok with the input
!!!OUT:
!!!    see error value enumeration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function sanitycheck() result(error)
    implicit none
    integer :: error
    integer :: i, nval
    error = getopts_OK
    !
    !check min/max for each key
    !

    !character
    if(.not.associated(clongkey)) clongkey=>safemem_realloc(clongkey,CLEN,0)
    if(.not.associated(cshortkey)) cshortkey=>safemem_realloc(cshortkey,len(cshortkey),0)
    do i = 1, ubound(clongkey,1)
       nval = getopts_narg(clongkey(i))
       call flush(0)
       if(nval < cminmax(i,1))then
          error=getopts_TOOFEW
          call rism_report_warn("Too few values for '--"//trim(clongkey(i))//"'")
          return
       elseif(nval > cminmax(i,2))then
          error=getopts_TOOMANY
          call rism_report_warn("Too many values for '--"//trim(clongkey(i))//"'")
          return
       end if
    end do

    !integer
    if(.not.associated(ilongkey)) ilongkey=>safemem_realloc(ilongkey,CLEN,0)
    if(.not.associated(ishortkey)) ishortkey=>safemem_realloc(ishortkey,len(ishortkey),0)
    do i = 1, ubound(ilongkey,1)
       nval = getopts_narg(ilongkey(i))
       if(nval < iminmax(i,1))then
          error=getopts_TOOFEW
          call rism_report_warn("Too few values for '--"//trim(ilongkey(i))//"'")
          return
       elseif(nval > iminmax(i,2))then
          error=getopts_TOOMANY
          call rism_report_warn("Too many values for '--"//trim(ilongkey(i))//"'")
          return
       end if
    end do

    !_REAL_
    if(.not.associated(rlongkey)) rlongkey=>safemem_realloc(rlongkey,CLEN,0)
    if(.not.associated(rshortkey)) rshortkey=>safemem_realloc(rshortkey,len(rshortkey),0)
    do i = 1, ubound(rlongkey,1)
       nval = getopts_narg(rlongkey(i))
       if(nval < rminmax(i,1))then
          error=getopts_TOOFEW
          call rism_report_warn("Too few values for '--"//trim(rlongkey(i))//"'")
          return
       elseif(nval > rminmax(i,2))then
          error=getopts_TOOMANY
          call rism_report_warn("Too many values for '--"//trim(rlongkey(i))//"'")
          return
       end if
    end do
    
  end function sanitycheck
end module getopts_c
