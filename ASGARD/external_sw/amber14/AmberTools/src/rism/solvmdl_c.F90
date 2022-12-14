#include "../include/dprec.fh"

!The 3D-RISM-KH software found here is copyright (c) 2011-2012 by 
!Tyler Luchko and David A. Case.
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

!#include "precision.fh"
module solvMDL_c
  use rism_report_c
  implicit none
  type solvMDL
     !title : description of SolvMDL
     character(len=80) :: title
     !nsite : number of sites in the solvMDL 
     integer :: nsite
     !ntype : number of atom types in the solvMDL (maybe more than one of each atom type)
     integer :: ntype
     !atmtyp : atom type number
     integer, pointer ::  atmtyp(:)=>NULL()
     !atmname(ntype) : character name of the atom types
     character(len=4), pointer :: atmname(:)=>NULL()
     !mass(ntype) : mass of each type (g/mol)
     _REAL_, pointer :: mass(:)=>NULL()
     !charge(ntype) : charge on each type (e*18.2223) (Amber units)
     _REAL_, pointer :: charge(:)=>NULL()
     !rmin(ntype) : lennard-jones rmin parameter
     _REAL_, pointer :: rmin(:)=>NULL()
     !epsilon(ntype) : lennard-jones epsilon parameter (kcal/mol)
     _REAL_, pointer :: epsilon(:)=>NULL()
     !multi(ntype) : multiplicity of each type
     integer, pointer :: multi(:)=>NULL()
     !coord(3,max(multi),ntype) : coodinamtes of each site
!     _REAL_, pointer :: coord(:,:,:)=>NULL()
     _REAL_, pointer :: coord(:,:)=>NULL()
  end type solvMDL
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!creates an new solvMDL object from an MDL file
!!!IN:
!!!   this : new SolvMDL object
!!!   file : file name to create it from
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine solvMDL_new(this,file)
    use rism_util
    use rism_parm
    use constants, only : COULOMB_CONST_E, KB, AMBER_ELECTROSTATIC, AVOGADRO
    implicit none
    type(solvMDL),intent(inout) :: this
    character(len=*), intent(in) :: file
    integer :: unit, iostat,iok
    character(len=256) :: fmt
    _REAL_ :: chgErr
    _REAL_ :: unitChg
    unitChg=sqrt(COULOMB_CONST_E/KB)
    call rism_parm_nxtsec_init()
    unit = freeUnit()
    open(unit=unit, file=file,status='old', iostat=iostat)
    if(iostat /= 0) &
       call rism_report_error('(a,i5)',"could not open "//trim(file)//":",iostat)

    call rism_parm_nxtsec(unit,6,0,'','TITLE',fmt,iok)
    read(unit, fmt) this%title

    call rism_parm_nxtsec(unit,6,0,'','POINTERS',fmt,iok)
    read(unit, fmt) this%nsite, this%ntype

    call solvMDL_allocate(this,this%nsite,this%ntype)

    call rism_parm_nxtsec(unit,6,0,'','MULTI',fmt,iok)
    read(unit, fmt) this%multi

    call rism_parm_nxtsec(unit,6,0,'','ATMTYP',fmt,iok)
    read(unit, fmt) this%atmtyp

    call rism_parm_nxtsec(unit,6,0,'','ATMNAME',fmt,iok)
    read(unit, fmt) this%atmname

    call rism_parm_nxtsec(unit,6,0,'','MASS',fmt,iok)
    read(unit, fmt) this%mass

    !read in charge and ensure that it is a integer number of unit charges
    call rism_parm_nxtsec(unit,6,0,'','CHG',fmt,iok)
    read(unit, fmt) this%charge
    !convert to k A / e (we don't know T)
    this%charge = this%charge/AMBER_ELECTROSTATIC*sqrt(COULOMB_CONST_E/KB)

    chgErr = sum(this%charge*this%multi)/unitChg
    chgErr =  (chgErr - nint(chgErr))/this%nsite
    this%charge = this%charge - chgErr*unitChg

    call rism_parm_nxtsec(unit,6,0,'','LJEPSILON',fmt,iok)
    read(unit, fmt) this%epsilon
    !convert to kT (we don't know T)
    this%epsilon = this%epsilon/KB

    call rism_parm_nxtsec(unit,6,0,'','LJSIGMA',fmt,iok)
    read(unit, fmt) this%rmin
    this%rmin = this%rmin*2d0

    call rism_parm_nxtsec(unit,6,0,'','COORD',fmt,iok)
    read(unit, fmt) this%coord

    close(unit)
  end subroutine solvMDL_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!Allocate memory for a solvMDL object
!!!IN:
!!!   this : solvMDL object
!!!   nsite : number of sites in the solvMDL 
!!!   ntype : number of atom types in the solvMDL (maybe more than one of each 
!!!           atom type) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine solvMDL_allocate(this,nsite,ntype)
    use safemem
    implicit none
    type(solvMDL), intent(inout) :: this
    integer, intent(in) :: nsite, ntype
    this%atmtyp => safemem_realloc(this%atmtyp,ntype)
    this%atmname => safemem_realloc(this%atmname,len(this%atmname),ntype)
    this%mass => safemem_realloc(this%mass,ntype)
    this%charge => safemem_realloc(this%charge,ntype)
    this%rmin => safemem_realloc(this%rmin,ntype)
    this%epsilon => safemem_realloc(this%epsilon,ntype)
    this%multi => safemem_realloc(this%multi,ntype)
    this%coord => safemem_realloc(this%coord,3,nsite)
  end subroutine solvMDL_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!Deallocates memory for a solvMDL object
!!!IN:
!!!   this : solvMDL object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  subroutine solvMDL_destroy(this)
    use safemem
    implicit none
    type(solvMDL), intent(inout) :: this
    integer :: err=0
    err = err+safemem_dealloc(this%atmtyp)
    err = err+safemem_dealloc(this%atmname)
    err = err+safemem_dealloc(this%mass)
    err = err+safemem_dealloc(this%charge)
    err = err+safemem_dealloc(this%rmin)
    err = err+safemem_dealloc(this%epsilon)
    err = err+safemem_dealloc(this%multi)
    err = err+safemem_dealloc(this%coord)
    if(err /= 0)&
         call rism_report_error("failure deallocating SolvMDL")
  end subroutine solvMDL_destroy
end module solvMDL_c
