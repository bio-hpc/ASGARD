!<compile=optimized>

!The 3D-RISM-KH software found here is copyright (c) 2010-2012 by 
!Andriy Kovalenko, Tyler Luchko and David A. Case.
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Object for solute information coming from 1D-RISM to be used in 3D-RISM.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module rism3d_solu_c
  use safemem
  use rism_report_c
  implicit none
  type rism3d_solu
!     sequence
     !natom :: number of solute atoms
     integer :: natom=0
     !mass :: mass of each solute atom [au]
     !charge:: partial charge of each solute atom [sqrt(kT A)]
     !origCharge:: Holds the orignal values in charge() when _unsetcharges() 
     !             is used.  Otherwise it is not used. [sqrt(kT A)]
     !ratu :: cartesian cooridinates of each solute atom (dim,atom) [A] 
     !sig_s:: Lennard-Jones parameter sigma* or r_min [A]
     !eps  :: Lennard-Jones parameter epsilon [kT]
     !    U_LJ = eps( (sig_s/r)**12 - 2*(sig_s/r)**6)
     _REAL_,pointer ::mass(:)=>NULL(),charge(:)=>NULL(),origCharge(:)=>NULL(),ratu(:,:)=>NULL(),&
          sig_s(:)=>NULL(),eps(:)=>NULL()

     !totCharge : solute total charge
     _REAL_ :: totCharge

     !charged  :: .false. if all partial charges are zero
     logical :: charged
  end type rism3d_solu

interface rism3d_solu_new
   module procedure rism3d_solu_new_internal, rism3d_solu_new_sander, rism3d_solu_new_parm7
end interface

interface rism3d_solu_setCoord
   module procedure rism3d_solu_setCoord_ratu, rism3d_solu_setCoord_crd
end interface rism3d_solu_setCoord

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!private subroutines and functions  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
private :: allocate_solu, setCharge

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Constructor using internal data representation
!!!IN:
!!!   this   :: solute object
!!!   natom  :: number of atoms
!!!   mass   :: mass of each atom [au]
!!!   charge :: charge of each atom [sqrt(kT A)]
!!!   ratu   :: atom coordinates [A]
!!!   sig_s  :: r_min/2 [A]
!!!   eps    :: LJ epsilon [kT]
!!!   o_mpicomm :: (optional) MPI communicator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_solu_new_internal (this,natom,mass,charge,ratu,sig_s,eps, o_mpicomm)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif /*MPI*/    
    type(rism3d_solu),intent(inout) :: this
    integer,intent(in) :: natom
    _REAL_,intent(in) :: mass(natom),charge(natom),ratu(3,natom),sig_s(natom),eps(natom)
    integer, optional, intent(in) :: o_mpicomm
    integer :: mpicomm, mpirank, mpisize, err
    mpicomm = 0
    mpirank = 0
    mpisize = 1
#ifdef MPI
    if(present(o_mpicomm)) then
       mpicomm = o_mpicomm
       if(mpicomm == MPI_COMM_NULL)&
            call rism_report_error("RISM3D_SOLU: received NULL MPI communicator")
       call mpi_comm_rank(mpicomm,mpirank,err)
       if(err /=0) call rism_report_error&
            ("(a,i8)","RISM3D SOLU: could not get MPI rank for communicator ",mpicomm)
       call mpi_comm_size(mpicomm,mpisize,err)
       if(err /=0) call rism_report_error&
            ("(a,i8)","RISM3D SOLU: could not get MPI size for communicator ",mpisize)
    end if
#endif
    if(mpirank == 0)then
       this%natom = natom
       call allocate_solu(this)
       call setCharge(this,charge)
       this%mass = mass
       this%sig_s = sig_s
       this%eps = eps
       call rism3d_solu_setCoord(this,ratu)
    end if
#ifdef MPI
    if(mpisize /=1)&
         call rism3d_solu_mpi_clone(this,mpirank,mpicomm)
#endif /*MPI*/
  end subroutine rism3d_solu_new_internal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Constructor for dealing with Amber's SANDER internal molecular representation
!!!IN:
!!!   this :: solute object
!!!   natom  :: number of atoms
!!!   ntype  :: number of atom types
!!!   iac    ::
!!!   ioc    ::
!!!   charge :: charge of each atom [e*18.2223]
!!!   cn1    :: LJ A parameter [kcal/mol*A^12]
!!!   cn2    :: LJ B parameter [kcal/mol*A^6]
!!!   mass   :: mass of each atom [au]
!!!   temperature :: temperature of the solvent
!!!   o_mpicomm :: (optional) MPI communicator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_solu_new_sander (this,natom,ntype,iac,ico,charge,cn1,cn2,&
       mass, temperature, o_mpicomm)
    use constants, only : KB, COULOMB_CONST_E, BOLTZMANN, AVOGADRO, AMBER_ELECTROSTATIC
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif /*MPI*/    
    type(rism3d_solu),intent(inout) :: this
    integer,intent(in) :: natom, iac(natom), ico(ntype**2),ntype
    _REAL_,intent(in) :: charge(natom), cn1(NTYPE*(NTYPE+1)/2), &
         cn2(NTYPE*(NTYPE+1)/2), mass(natom), temperature
    integer, optional, intent(in) :: o_mpicomm
    
    integer :: id,iu,iv,ir
    !sigu :: LJ sigma for solute
    !epsu :: LJ epsilon for solute
    _REAL_ :: sigu(natom), epsu(natom)

    integer :: mpicomm, mpirank, mpisize, err
    mpicomm = 0
    mpirank = 0
    mpisize = 1
#ifdef MPI
    if(present(o_mpicomm)) then 
       mpicomm = o_mpicomm
       if(mpicomm == MPI_COMM_NULL)&
            call rism_report_error("RISM3D_SOLU: received NULL MPI communicator")
       call mpi_comm_rank(mpicomm,mpirank,err)
       if(err /=0) call rism_report_error&
            ("(a,i8)","RISM3D SOLU: could not get MPI rank for communicator ",mpicomm)
       call mpi_comm_size(mpicomm,mpisize,err)
       if(err /=0) call rism_report_error&
            ("(a,i8)","RISM3D SOLU: could not get MPI size for communicator ",mpisize)
    end if
#endif
    if(mpirank == 0)then
    
       this%natom = natom
       call allocate_solu(this)
       !UNIT CONVERSION 
    !from [e*18.2223] to [sqrt(kT A)]
       call setCharge(this,charge(1:natom)/AMBER_ELECTROSTATIC&
            *sqrt(COULOMB_CONST_E/(KB *temperature)))
       this%mass = mass
!!$    this%sig_s(1)=1.77670000d+00
!!$    this%sig_s(2)=2.24492410d-01
!!$    this%sig_s(3)=2.24492410d-01
       do iu=1,natom
          id = ico(ntype*(iac(iu)-1)+iac(iu))
          !!FIX - is this the proper way of handling no LJ?  No.  The best is to find the volume of the 
          !!      enclosing atom and calculating a radius for the enclosed atom that coincides with this
          !!      and a well depth of 1/10 the parent atom
          if(cn2(id)==0d0)then
             this%sig_s(iu)=0.7d0
!!$          this%sig_s(iu)=2.24492410e-01
          else
             this%sig_s(iu) = (2*cn1(id)/cn2(id))**(1.d0/6.d0)/2d0
          endif
          !!endfix
!!$       write(0,*) "SIGMA", iu, this%sig_s(iu)
       enddo
!!$    this%eps(1)=1.55300000d-01
!!$    this%eps(2)=4.60325047d-02
!!$    this%eps(3)=4.60325047d-02
       do iu=1,natom
          id = ico(ntype*(iac(iu)-1)+iac(iu))
          if(cn1(id)==0d0)then
             this%eps(iu)=1d-2
!!$          this%eps(iu)=4.48984819d-01
          else
             this%eps(iu) = cn2(id)**2/(4*cn1(id))
             this%eps(iu) = this%eps(iu)
          endif
          this%eps(iu) = this%eps(iu)/(KB*temperature)
!!$       write(0,*) "EPS", iu, this%eps(iu)
       enddo
    end if
#ifdef MPI
    if(mpisize /=1)&
       call rism3d_solu_mpi_clone(this,mpirank,mpicomm)
#endif /*MPI*/
  end subroutine rism3d_solu_new_sander

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Constructor for dealing with Amber parm7 and crd7 files
!!!IN:
!!!   this :: solute object
!!!   parm :: parameter file name
!!!   temperature :: temperature of the solvent [K]
!!!   o_mpicomm :: (optional) MPI communicator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_solu_new_parm7(this, parm, temperature, o_mpicomm)
    use rism_util, only : freeUnit
    use rism_parm
    implicit none
    type(rism3d_solu),intent(inout) :: this
    character(len=*), intent(in) :: parm
    _REAL_, intent(in) :: temperature
    integer, optional, intent(in) :: o_mpicomm

    character(len=80) fmt
    character(len=80) fmtin,ifmt,afmt,rfmt,type, line
    integer :: i, iok, unit
    integer :: natom,ntype,nbonh,mbona,ntheth,mtheta,nphih,mphia, &
         nhparm,nparm,nnb,nres,nbona,ntheta,nphia, &
         numbnd,numang,nptra,natyp,nphb,ifpert,nbper,ngper,ndper, &
         mbper,mgper,mdper,ifbox,nmxrs,ifcap,numextra,ncopy
    integer,pointer :: iac(:)=>NULL(), ico(:)=>NULL()
    _REAL_,pointer :: charge(:)=>NULL(), cn1(:)=>NULL(), &
         cn2(:)=>NULL(), mass(:)=>NULL()
    integer :: mpicomm
    mpicomm = 0
    if(present(o_mpicomm)) mpicomm = o_mpicomm
    
    !this liberally uses code from sander's rdparm.f
    unit = freeUnit()
    open (unit,file=parm,status='old')

    ifmt = '(12I6)'
    afmt = '(20A4)'
    rfmt = '(5E16.8)'
    call rism_parm_nxtsec_init()

    fmtin = ifmt
    type = 'POINTERS'
    call rism_parm_nxtsec(unit,  6,  0,fmtin,  type,  fmt,  iok)
    read(unit,fmt) natom,ntype,nbonh,mbona,ntheth,mtheta,nphih,mphia, &
         nhparm,nparm,nnb,nres,nbona,ntheta,nphia, &
         numbnd,numang,nptra,natyp,nphb,ifpert,nbper,ngper,ndper, &
         mbper,mgper,mdper,ifbox,nmxrs,ifcap,numextra,ncopy

    !allocate some temporary storage.  This uses more memory than we
    !would like but we can use the Sander constructor directly this
    !way
    iac => safemem_realloc(iac,natom,.false.)
    ico => safemem_realloc(ico,ntype**2,.false.)
    charge => safemem_realloc(charge,natom,.false.)
    cn1 => safemem_realloc(cn1,Ntype*(Ntype+1)/2,.false.)
    cn2 => safemem_realloc(cn2,Ntype*(Ntype+1)/2,.false.)
    mass => safemem_realloc(mass,natom,.false.)

    fmtin = ifmt
    type = 'ATOM_TYPE_INDEX'
    call rism_parm_nxtsec(unit,  6,  0,fmtin,  type,  fmt,  iok)
    read(unit,fmt) iac

    fmtin = ifmt
    type = 'NONBONDED_PARM_INDEX'
    call rism_parm_nxtsec(unit,  6,  0,fmtin,  type,  fmt,  iok)
    read(unit,fmt) ico

    fmtin = rfmt
    type = 'CHARGE'
    call rism_parm_nxtsec(unit,  6,  0,fmtin,  type,  fmt,  iok)
    read(unit,fmt) charge

    fmtin = rfmt
    type = 'LENNARD_JONES_ACOEF'
    call rism_parm_nxtsec(unit,  6,  0,fmtin,  type,  fmt,  iok)
    read(unit,fmt) cn1

    fmtin = rfmt
    type = 'LENNARD_JONES_BCOEF'
    call rism_parm_nxtsec(unit,  6,  0,fmtin,  type,  fmt,  iok)
    read(unit,fmt) cn2

    fmtin = rfmt
    type = 'MASS'
    call rism_parm_nxtsec(unit,  6,  0,fmtin,  type,  fmt,  iok)
    read(unit,fmt) mass
    
    close(unit)    
    call rism3d_solu_new_sander (this,natom,ntype,iac,ico,charge,cn1,cn2,&
       mass, temperature,mpicomm)

    if(safemem_dealloc(iac) /= 0) call rism_report_error("RISM3D_SOLU: Failed to deallocate IAC")
    if(safemem_dealloc(ico) /= 0) call rism_report_error("RISM3D_SOLU: Failed to deallocate ICO")
    if(safemem_dealloc(charge) /= 0) call rism_report_error("RISM3D_SOLU: Failed to deallocate CHARGE")
    if(safemem_dealloc(cn1) /= 0) call rism_report_error("RISM3D_SOLU: Failed to deallocate CN1")
    if(safemem_dealloc(cn2) /= 0) call rism_report_error("RISM3D_SOLU: Failed to deallocate CN2")
    if(safemem_dealloc(mass) /= 0) call rism_report_error("RISM3D_SOLU: Failed to deallocate MASS")

  end subroutine rism3d_solu_new_parm7

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Clone constructor using internal data representation
!!!IN:
!!!   this :: object to be copied
!!!   clone :: clone of the object.  No memory space is shared.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_solu_clone (this,clone)
    implicit none
    type(rism3d_solu),intent(inout) :: clone
    type(rism3d_solu),intent(in) :: this
    
    call rism3d_solu_new(clone,this%natom,this%mass,this%charge,this%ratu,&
         this%sig_s,this%eps)
  end subroutine rism3d_solu_clone

#ifdef MPI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Allocates memory on non-master nodes and then distributes information out 
!!!from the master.  It is assumed that the object on the master already exists.
!!!IN:
!!!   this :: solute object
!!!   rank :: MPI rank
!!!   comm :: MPI communicator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_solu_mpi_clone(this,rank,comm)
    implicit none
    type(rism3d_solu),intent(inout) :: this
    integer, intent(in) :: rank,comm
    integer :: err
    include 'mpif.h'
    !first distribute the pieces of information needed to allocate memory
    call mpi_bcast(this%natom,1,mpi_integer,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLU: could not broadcast NATOM")

    !non-master processes should now allocate memory
    if(rank /= 0) then
       call allocate_solu(this)
    end if

    !now distribute the arrays to the non-master processes
    call mpi_bcast(this%charge,size(this%charge),mpi_double_precision,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLU: could not broadcast CHARGE")
    call mpi_bcast(this%mass,size(this%mass),mpi_double_precision,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLU: could not broadcast MASS")
    call mpi_bcast(this%ratu,size(this%ratu,1)*size(this%ratu,2),mpi_double_precision,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLU: could not broadcast RATU")
    call mpi_bcast(this%sig_s,size(this%sig_s),mpi_double_precision,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLU: could not broadcast SIG_S")
    call mpi_bcast(this%eps,size(this%eps),mpi_double_precision,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLU: could not broadcast EPS")
    
    if(rank /=0) call setCharge(this,this%charge)
  end subroutine rism3d_solu_mpi_clone
#endif /*MPI*/

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Set the solute coordinates
!!!IN:
!!!   this :: solute object
!!!   ratu :: coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_solu_setCoord_ratu(this,ratu)
    implicit none
    type(rism3d_solu), intent(inout) :: this
    _REAL_, intent(in) :: ratu(:,:)
    if(ubound(ratu,1) /= 3 .and. ubound(ratu,2) /= this%natom)then
       call rism_report_error("solute size and coordinate array do not match")
       stop
    end if
    this%ratu = ratu
  end subroutine rism3d_solu_setCoord_ratu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Set the solute coordinates
!!!IN:
!!!   this :: solute object
!!!   crd :: coordinate file. Amber crd/rst or PDB 3.2 file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_solu_setCoord_crd(this, crd)
    use rism_util, only : freeUnit
    implicit none
    type(rism3d_solu), intent(inout) :: this
    character(len=*), intent(in) :: crd
    character(len=128) :: test
    integer :: unit, natom, iostat
    unit = freeUnit()
    open (unit,file=crd,status='old')

    !
    !first try Amber format
    !
    !read off title and discard
    read(unit,*)
    !read in the number of atoms, we can discard the time if it is there
    read(unit,*) test
    read(test,*,iostat=iostat) natom
    if(natom==this%natom .and. iostat==0) then
       read(unit,'(6F12.7)') this%ratu
    else
       !
       !try PDB format
       !
       natom = 0
       rewind(unit)
       iostat=0
       do 
          read(unit,'(a)',iostat=iostat) test
          if(iostat <0) exit
          if(iostat >0) call rism_report_error("RISM3D_SOLU: failed on reading PDB")
          if(test(1:6) .ne. "ATOM  " .and. test(1:6) .ne. "HETATM") cycle
          natom = natom+1
          if(natom > this%natom) &
               call rism_report_error('(a,i6,a)',&
               "RISM3D_SOLU: more PDB atoms than present in parameter file (",this%natom,")")
          read(test(31:54),'(3(f8.3))') this%ratu(:,natom)
       end do
       if(natom /= this%natom)&
            call rism_report_error("RISM3D_SOLU: failed reading coordinate file")
    end if
    close(unit)
  end subroutine rism3d_solu_setCoord_crd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Indicates if the partial charges are on or not
!!!OUT:
!!!    Returns true if partial charges are set, false if they are no
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism3d_solu_charged(this) result(charged)
    implicit none
    type(rism3d_solu), intent(inout) :: this
    logical :: charged
    charged=.not.associated(this%origCharge)
  end function rism3d_solu_charged

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets all solute partial charges to zero
!!!IN:
!!!   this :: solute object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_solu_unsetCharges(this)
      implicit none
      type(rism3d_solu), intent(inout) :: this
      this%origCharge => safemem_realloc(this%origCharge,ubound(this%charge,1))
      this%origCharge = this%charge
      call setCharge(this,(/0d0/))
    end subroutine rism3d_solu_unsetCharges

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets all solute partial charges to to their original
!!!values. (Undoes rism3d_solu_unsetCharge().)
!!!IN:
!!!   this :: solute object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_solu_resetCharges(this)
      implicit none
      type(rism3d_solu), intent(inout) :: this
      call setCharge(this,this%origCharge)
      if(safemem_dealloc(this%origCharge)/=0)&
           call rism_report_error("unable to deallocate 'origCharge'")
    end subroutine rism3d_solu_resetCharges

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!deconstructor
!!!IN:
!!!   this :: solute object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_solu_destroy(this)
    implicit none
    type(rism3d_solu) :: this
    if(safemem_dealloc(this%mass) /=0 )&
         call rism_report_error("RISM3D_SOLU: could not deallocate MASS")
    if(safemem_dealloc(this%charge) /=0 )&
         call rism_report_error("RISM3D_SOLU: could not deallocate CHARGE")
    if(safemem_dealloc(this%origCharge) /=0 )&
         call rism_report_error("RISM3D_SOLU: could not deallocate ORIGCHARGE")
    if(safemem_dealloc(this%ratu) /=0 )&
         call rism_report_error("RISM3D_SOLU: could not deallocate RATU")
    if(safemem_dealloc(this%sig_s) /=0 )&
         call rism_report_error("RISM3D_SOLU: could not deallocate SIG_s")
    if(safemem_dealloc(this%eps) /=0 )&
         call rism_report_error("RISM3D_SOLU: could not deallocate EPS")
  end subroutine rism3d_solu_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!private subroutines and functions  
!!!IN:
!!!   this :: solute object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine allocate_solu(this)
    implicit none
    type(rism3d_solu),intent(inout) :: this
    this%mass => safemem_realloc(this%mass,this%natom,.false.)
    this%charge => safemem_realloc(this%charge,this%natom,.false.)
    this%sig_s => safemem_realloc(this%sig_s,this%natom,.false.)
    this%eps => safemem_realloc(this%eps,this%natom,.false.)
    this%ratu => safemem_realloc(this%ratu,3,this%natom,.false.)
    !ensure that ratu is initialized in case we have to do an MPI broadcast
    !before getting coordinates
    this%ratu = 0d0
  end subroutine allocate_solu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!private subroutine to take care of book keeping when setting charges
!!!IN:
!!!   this :: solute object
!!!   charge :: array of charges, one per charge site or a single
!!!             value may be passed which will be applied to all sites
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine setCharge(this,charge)
    implicit none
    type(rism3d_solu),intent(inout) :: this
    _REAL_, intent(in) :: charge(:)
    if(ubound(charge,1) /= ubound(this%charge,1) .and.&
         ubound(charge,1)/=1)&
         call rism_report_error("Bad solute charge set")
    if(ubound(charge,1)/=1)then
       this%charge = charge
    else
       this%charge = charge(1)
    end if
    this%totCharge = sum(this%charge)
    if(sum(abs(this%charge))>0d0)then
       this%charged=.true.
    else
       this%charged=.false.
    end if
  end subroutine setCharge
end module rism3d_solu_c
