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
!!!Object for solvent information coming from 1D-RISM to be used in 3D-RISM.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module rism3d_solv_c
  use rism_report_c
  use safemem
  implicit none
  type rism3d_solv
!     sequence
     !dr          :: grid spacing from RISM1D calculation [A]
     !dt          :: Fourier grid spacing from RISM1D calculation dt = pi/nr/dr. [1/A]
     !               Currently this is used only to set up the solvent calculations.
     !temperature :: solvent temperature [K]
     !dielconst   :: solvent dielectric constant
     !xappa       :: inverse Debye length [1/A]
     !xikt        :: compressibility [A^3]
     !smear       :: charge smear parameter required for long range electrostatics [A]
     !xikt_dT     :: compressibility temperature derivative [A^3/K]
     _REAL_ :: dr=0, dt=0, temperature=0, dielconst=0, xappa=0, xikt=0, smear=0,&
          xikt_dT
     !natom       :: number of solvent atom types
     !nspecies    :: number of molecular species
     integer ::  natom=0, nspecies
     !nr :: the number of points in the solvent-solvent RDF
     integer :: nr=0


     !mult      :: the multiplicity of each atom in the solvent model
     !natomspecies :: number of solvent sites for each species
     integer,pointer :: mult(:)=>NULL(), natomspecies(:)=>NULL()

     !atomname :: name of the solvent atoms.  Used for the GUV,CUV and HUV output files.
     character(len=4),pointer :: atomname(:)=>NULL()

     !xvv        :: original solvent chi
     !fourier_tbl   :: Precalculated table of the radial spacing for 1D-RISM in Fourier space.
     !              dt=pi/(nr*dr) This is then used for interpolating Xvv for the needs of 3D-RISM
     !charge         :: solvent atom charge by type [sqrt(kT A)]
     !charge_sp       :: total solvent species charge for each atom by type [sqrt(kT A)]
     !rho       :: solvent density by atom type [#/A^3]
     !rho_sp       :: solvent density by molecular species [#/A^3]
     !sig_s       :: solvent LJ-sigma by atom type [A]
     !eps       :: solvent LJ-epsilon by atom type [kT]
     !delhv0     :: -Lim_k->0 ( Sum_v1 Qv1*Xv1v2(k)4pi/k^2 - hlkv0 ) from 1D-RISM
     !delhv0_dT     :: temperature derivative analogue
     _REAL_,pointer :: fourier_tbl(:)=>NULL(), xvv(:,:,:)=>NULL(),&
          charge(:)=>NULL(),charge_sp(:)=>NULL(),rho(:)=>NULL(),rho_sp(:)=>NULL(),&
          sig_s(:)=>NULL(),eps(:)=>NULL(),delhv0(:)=>NULL(),&
          delhv0_dT(:)=>NULL(), xvv_dT(:,:,:)=>NULL()

     !ionic :: .true. if the any of the solvent species has non-zero charge
     logical :: ionic

     real :: xvv_version
  end type rism3d_solv

  interface rism3d_solv_new
     module procedure rism3d_solv_new_all, rism3d_solv_new_readxvv
  end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!public subroutines and functions  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  public :: rism3d_solv_new, rism3d_solv_clone, rism3d_solv_destroy
#ifdef MPI
  public :: rism3d_solv_mpi_clone
#endif /*MPI*/

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!private subroutines and functions  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  private :: allocate_solv, calc_fourier_tbl, readxvv1, readxvv2

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!constructor.  Creates new solvent object be specifying all data. if
!!!this is an MPI run, only the rank 0 process needs valid arguments (other than
!!!the communicator).
!!!IN:
!!!   this       :: the new object
!!!   dr         :: grid spacing from RISM1D calculation [A]
!!!   dt         :: Fourier grid spacing from RISM1D calculation dt = pi/nr/dr. [1/A]
!!!                 Currently this is used only to set up the solvent calculations.
!!!   temperature:: solvent temperature [K]
!!!   dielconst  :: solvent dielectric constant
!!!   xappa      :: inverse Debye length [1/A]
!!!   xikt       :: compressibility [A^3]
!!!   xikt_dT    :: compressibility temperature derivative [A^3/K]
!!!   smear      :: ???
!!!   natom      :: number of solvent atom types
!!!   nspecies   :: number of molecular species
!!!   natomspecies :: number of solvent sites for each molecular species
!!!   nr         :: the number of points in the solvent-solvent RDF
!!!   mult      :: the multiplicity of each atom in the solvent model
!!!   atomname    :: name of the solvent atoms.  Used for the GUV,CUV and HUV output files.
!!!   xvv        :: original solvent chi
!!!   fourier_tbl:: Precalculated table of the radial spacing for 1D-RISM in Fourier space.
!!!                 dt=pi/(nr*dr) This is then used for interpolating Xvv for the needs of 3D-RISM
!!!   charge         :: solvent atom charge by type [sqrt(kT A)]
!!!   charge_sp       :: total solvent species charge for each atom by type [sqrt(kT A)]
!!!   rho       :: solvent density by atom type [#/A^3]
!!!   rho_sp       :: solvent density by molecular species [#/A^3]
!!!   sig_s       :: solvent LJ-sigma by atom type [A]
!!!   eps       :: solvent LJ-epsilon by atom type [kT]
!!!   delhv0     :: ???
!!!   o_mpicomm :: (optional) MPI communicator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_solv_new_all(this,dr, dt, temperature, dielconst, xappa, &
       xikt, xikt_dT, smear,natom,nspecies, natomspecies, nr,mult, atomname,&
       fourier_tbl,xvv, charge,charge_sp,rho,rho_sp,sig_s,eps,delhv0,&
       delhv0_dT,xvv_dT,&
       o_mpicomm)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    type(rism3d_solv),intent(inout) :: this
    _REAL_, intent(in) :: dr, dt, temperature, dielconst, xappa, xikt, xikt_dT, smear
    integer, intent(in) ::  natom,nspecies,nr
    integer, intent(in) :: mult(natom),natomspecies(nspecies)
    character(len=4), intent(in) :: atomname(natom)
    _REAL_,intent(in) :: fourier_tbl(nr), xvv(nr,natom,natom),&
         charge(natom),charge_sp(natom),rho(natom),rho_sp(nspecies),&
         sig_s(natom),eps(natom),delhv0(natom),&
         delhv0_dT(natom),xvv_dT(nr,natom,natom)
    integer, optional, intent(in) :: o_mpicomm
    integer :: mpicomm, mpirank, mpisize, err


    mpicomm =0
    mpirank = 0
    mpisize = 1
#ifdef MPI
    if(present(o_mpicomm)) then
       mpicomm = o_mpicomm
       if(mpicomm == MPI_COMM_NULL)&
            call rism_report_error("RISM3D_SOLV: received NULL MPI communicator")
       call mpi_comm_rank(mpicomm,mpirank,err)
       if(err /=0) call rism_report_error&
            ("(a,i8)","RISM3D INIT: could not get MPI rank for communicator ",mpicomm)
       call mpi_comm_size(mpicomm,mpisize,err)
       if(err /=0) call rism_report_error&
            ("(a,i8)","RISM3D SOLV: could not get MPI size for communicator ",mpisize)
    end if
#endif /*MPI*/
    if(mpirank == 0) then
       this%dr = dr
       this%dt = dt
       this%temperature = temperature
       this%dielconst = dielconst
       this%xappa = xappa
       this%xikt = xikt
       this%smear = smear
       this%natom = natom
       this%nspecies = nspecies
       this%nr = nr
       call allocate_solv(this)
       this%xikt_dT = xikt_dT !after allocate so it is not overwritten with HUGE
       this%mult = mult
       this%natomspecies = natomspecies
       this%atomname = atomname
       this%fourier_tbl = fourier_tbl
       this%xvv = xvv
       this%charge = charge
       this%charge_sp = charge_sp
       this%rho = rho
       this%rho_sp = rho_sp
       this%sig_s = sig_s
       this%eps = eps
       this%delhv0 = delhv0

       this%delhv0_dT = delhv0_dT

       call calc_fourier_tbl(this)    
    end if
#ifdef MPI
    if(mpisize /=1)&
         call rism3d_solv_mpi_clone(this, mpirank,mpicomm)
#endif /*MPI*/    
       !optional - may not be defined
       if(all(this%delhv0_dT/=huge(1d0)))then
          this%xvv_dT = xvv_dT
       else
          this%delhv0_dT = huge(1d0)
          if(safemem_dealloc(this%xvv_dT)/=0)&
               call rism_report_error("Failed to deallocate xvv_dT")
       end if

       if(sum(abs(this%charge_sp)) > 0d0)then
          this%ionic = .true.
       else
          this%ionic=.false.
       end if
  end subroutine rism3d_solv_new_all

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!constructor.  Creates new solvent object directly from an Amber7 format XVV 
!!!file. if this is an MPI run, only the rank 0 process needs a valid Xvv file.
!!!IN:
!!!   this :: the new object
!!!   xvvfile :: name of the XVV file to read
!!!   o_mpicomm :: (optional) MPI communicator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_solv_new_readxvv(this,xvvfile, o_mpicomm)
    use rism_util, only : freeUnit
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    type(rism3d_solv),intent(inout) :: this
    character(*), intent(in) :: xvvfile
    integer, optional, intent(in) :: o_mpicomm
    integer :: mpicomm, mpirank, mpisize, err, unit
    integer :: iatom
    mpicomm =0
    mpirank = 0
    mpisize = 1
#ifdef MPI
    if(present(o_mpicomm)) then
       mpicomm = o_mpicomm
       if(mpicomm == MPI_COMM_NULL)&
            call rism_report_error("RISM3D_SOLV: received NULL MPI communicator")
       call mpi_comm_rank(mpicomm,mpirank,err)
       if(err /=0) call rism_report_error&
            ("(a,i8)","RISM3D INIT: could not get MPI rank for communicator ",mpicomm)
       call mpi_comm_size(mpicomm,mpisize,err)
       if(err /=0) call rism_report_error&
            ("(a,i8)","RISM3D SOLV: could not get MPI size for communicator ",mpisize)
    end if
#endif /*MPI*/
    if(mpirank == 0) then

       unit = freeUnit()
       !unit conversion done in readxvv()
       open (unit,file=xvvfile,status='old')
       call readxvv1(this,unit,this%natom,this%nspecies,this%nr)
       call allocate_solv(this)
       call readxvv2(this,unit)
       close(unit)
       call calc_fourier_tbl(this)
       !remove numerical error from reading file
       do iatom=1, this%natom
          if(abs(this%charge_sp(iatom)) < 1d-6)&
               this%charge_sp(iatom)=0
       end do

#ifdef MPI
    if(mpisize /=1)&
         call rism3d_solv_mpi_clone(this, mpirank,mpicomm)
#endif /*MPI*/    
       !check if delhv0_dT has been set.  If not free xvv_dT
       if(any(this%delhv0_dT==huge(1d0)))then
          this%delhv0_dT=huge(1d0)
          if(safemem_dealloc(this%xvv_dT)/=0)&
               call rism_report_error("Failed to deallocate xvv_dT")
       end if
    end if

    if(sum(abs(this%charge_sp)) > 0d0)then
       this%ionic = .true.
    else
       this%ionic=.false.
    end if
  end subroutine rism3d_solv_new_readxvv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Clone constructor.  Creates new solvent object that is identical to the 
!!!original.
!!!IN:
!!!   this :: object to be copied
!!!   clone :: clone of the object.  No memory space is shared.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_solv_clone(this,clone)
    implicit none
    type(rism3d_solv),intent(in) :: this
    type(rism3d_solv),intent(inout) :: clone
    call rism3d_solv_new(clone, this%dr, this%dt, this%temperature, this%dielconst, &
         this%xappa, this%xikt, this%xikt_dT, this%smear, &
         this%natom, this%nspecies, this%natomspecies, this%nr, this%mult, &
         this%atomname, this%fourier_tbl, this%xvv, this%charge, this%charge_sp, &
         this%rho, this%rho_sp, &
         this%sig_s, this%eps, this%delhv0, this%delhv0_dT, this%xvv_dT)
  end subroutine rism3d_solv_clone

#ifdef MPI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Allocates memory on non-master nodes and then distributes information out 
!!!from the master.  It is assumed that the object on the master already exists.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_solv_mpi_clone(this,rank,comm)
    implicit none
    type(rism3d_solv),intent(inout) :: this
    integer, intent(in) :: rank,comm
    integer :: err
    include 'mpif.h'
    !first distribute the pieces of information needed to allocate memory
    call mpi_bcast(this%nr,1,mpi_integer,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast NR")
    call mpi_bcast(this%natom,1,mpi_integer,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast NATOM")
    call mpi_bcast(this%nspecies,1,mpi_integer,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast NSPECIES")

    !non-master processes should now allocate memory
    if(rank /= 0) then
       call allocate_solv(this)
    end if

    !now distribute the arrays to the non-master processes
    call mpi_bcast(this%temperature,1,mpi_double_precision,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast TEMPERATURE")
    call mpi_bcast(this%dielconst,1,mpi_double_precision,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast DIELCONST")
    call mpi_bcast(this%xappa,1,mpi_double_precision,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast XAPPA")
    call mpi_bcast(this%xikt,1,mpi_double_precision,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast XIKT")
    call mpi_bcast(this%xikt_dT,1,mpi_double_precision,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast XIKT_DT")
    call mpi_bcast(this%smear,1,mpi_double_precision,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast SMEAR")
    call mpi_bcast(this%dr,1,mpi_double_precision,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast DR")
    call mpi_bcast(this%dt,1,mpi_double_precision,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast DT")
    call mpi_bcast(this%xvv,size(this%xvv,1)*size(this%xvv,3)&
         *size(this%xvv,3),mpi_double_precision,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast XVV")
    call mpi_bcast(this%fourier_tbl,size(this%fourier_tbl),&
         mpi_double_precision,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast FOURIER_TBL")
    call mpi_bcast(this%charge,size(this%charge),mpi_double_precision,0,&
         comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast CHARGE")
    call mpi_bcast(this%charge_sp,size(this%charge_sp),mpi_double_precision,0,&
         comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast CHARGE_SP")
    call mpi_bcast(this%delhv0,size(this%delhv0),mpi_double_precision,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast DELHV0")
    call mpi_bcast(this%rho,size(this%rho),mpi_double_precision,0,&
         comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast RHO")
    call mpi_bcast(this%rho_sp,size(this%rho_sp),mpi_double_precision,0,&
         comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast RHO_SP")
    call mpi_bcast(this%sig_s,size(this%sig_s),mpi_double_precision,0,&
         comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast SIG_S")
    call mpi_bcast(this%eps,size(this%eps),mpi_double_precision,0,&
         comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast EPS")

    call mpi_bcast(this%mult,size(this%mult),mpi_integer,0,&
         comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast MULT")
    call mpi_bcast(this%natomspecies,size(this%natomspecies),mpi_integer,0,&
         comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast NATOMSPECIES")
    call mpi_bcast(this%atomname,size(this%atomname)*4,mpi_character,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast ATOMNAME")

    !optional derivative data
    call mpi_bcast(this%delhv0_dT,size(this%delhv0_dT),mpi_double_precision,0,comm,err)
    if(err /=0) call rism_report_error&
         ("RISM3D_SOLV: could not broadcast DELHV0_DT")
    if(this%delhv0_dT(1) /=huge(1d0))then
       call mpi_bcast(this%xvv_dT,product(ubound(this%xvv_dT)),&
            mpi_double_precision,0,comm,err)
       if(err /=0) call rism_report_error&
            ("RISM3D_SOLV: could not broadcast XVVDT")
    end if
    !non-master nodes finish initialization
    if(rank /=0) then
       call calc_fourier_tbl(this)
    end if

    !this is cheaper to recalculate than to broadcast
    if(sum(abs(this%charge_sp)) > 0d0)then
       this%ionic = .true.
    else
       this%ionic=.false.
    end if
  end subroutine rism3d_solv_mpi_clone
#endif /*MPI*/

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!check if we can perform a temperature derivative calculation
!!!(i.e. all the necessary information is available)
!!!IN:
!!!   this : rism3d_solv object
!!!OUT:
!!!    .true. if we can, .false. if we can't
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism3d_solv_canCalc_DT(this) result(can_dT)
    implicit none
    type(rism3d_solv), intent(in) :: this
    logical :: can_dT

    can_dT = associated(this%xvv_dT)
  end function rism3d_solv_canCalc_DT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!destroyer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_solv_destroy(this)
    use safemem
    implicit none
    type(rism3d_solv) :: this
    integer :: err
    if(safemem_dealloc(this%mult)/=0)&
         call rism_report_error("RISM3D_SOLV: failed to deallocate MULT")
    if(safemem_dealloc(this%natomspecies)/=0)&
         call rism_report_error("RISM3D_SOLV: failed to deallocate NATOMSPECIES")
    if(safemem_dealloc(this%atomname)/=0)&
         call rism_report_error("RISM3D_SOLV: failed to deallocate MULT")
    if(safemem_dealloc(this%fourier_tbl)/=0)&
         call rism_report_error("RISM3D_SOLV: failed to deallocate MULT")
    if(safemem_dealloc(this%xvv)/=0)&
         call rism_report_error("RISM3D_SOLV: failed to deallocate MULT")
    if(safemem_dealloc(this%charge)/=0)&
         call rism_report_error("RISM3D_SOLV: failed to deallocate MULT")
    if(safemem_dealloc(this%charge_sp)/=0)&
         call rism_report_error("RISM3D_SOLV: failed to deallocate MULT")
    if(safemem_dealloc(this%rho)/=0)&
         call rism_report_error("RISM3D_SOLV: failed to deallocate MULT")
    if(safemem_dealloc(this%sig_s)/=0)&
         call rism_report_error("RISM3D_SOLV: failed to deallocate MULT")
    if(safemem_dealloc(this%eps)/=0)&
         call rism_report_error("RISM3D_SOLV: failed to deallocate MULT")
    if(safemem_dealloc(this%delhv0)/=0)&
         call rism_report_error("RISM3D_SOLV: failed to deallocate MULT")
    if(safemem_dealloc(this%delhv0_dT)/=0)&
         call rism_report_error("RISM3D_SOLV: failed to deallocate MULT")
  end subroutine rism3d_solv_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!allocate memory at the begining if the simulation.
!!!IN:
!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine allocate_solv(this) 
    use safemem
    implicit none

    type(rism3d_solv), intent(inout) :: this
    integer :: new_limit    ! new stack limit


    !allocate real memory
    this%fourier_tbl => safemem_realloc(this%fourier_tbl,this%nr,.false.)
    this%xvv => safemem_realloc(this%xvv,this%nr,this%natom,this%natom,.false.)
    this%charge => safemem_realloc(this%charge,this%natom,.false.)
    this%charge_sp => safemem_realloc(this%charge_sp,this%natom,.false.)
    this%eps => safemem_realloc(this%eps,this%natom,.false.)
    this%sig_s => safemem_realloc(this%sig_s,this%natom,.false.)
    this%rho => safemem_realloc(this%rho,this%natom,.false.)
    this%rho_sp => safemem_realloc(this%rho_sp,this%nspecies,.false.)
    this%delhv0 => safemem_realloc(this%delhv0,this%natom,.false.)

    this%delhv0_dT => safemem_realloc(this%delhv0_dT,this%natom,.false.)
    !flag that delhv0_dT not set
    this%delhv0_dT = huge(1d0)
    this%xvv_dT => safemem_realloc(this%xvv_dT,this%nr,this%natom,this%natom,.false.)

    !flag that xikt_dT not set
    this%xikt_dT = huge(1d0)

    !allocate integer memory
    this%mult => safemem_realloc(this%mult,this%natom,.false.)
    this%natomspecies => safemem_realloc(this%natomspecies,this%natom,.false.)

    !allocate character (holrith)memory
    this%atomname => safemem_realloc(this%atomname,4, this%natom,.false.)
  end subroutine allocate_solv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Pre-calculate reciprocal  grid spacing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine calc_fourier_tbl(this)
    use constants, only : PI
    implicit none
    type(rism3d_solv), intent(inout) :: this
    integer :: ir
#ifdef RISM_DEBUG
    write(6,*) "CALC_FOURIER_TBL"; call flush(6)
#endif /*RISM_DEBUG*/
    !compute the 1D Fourier grid spacing used for the RISM1D calculations
    this%dt = PI/(this%nr*this%dr)
    do ir = 1,this%nr
       this%fourier_tbl(ir) = this%dt*(ir-1)
    end do
  end subroutine calc_fourier_tbl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Reads in the header data from a Amber style Xvv file.  Determines the version
!!!and reads in POINTER information that determines the sizes of array in the
!!!rest of the file
!!!IN:
!!!   nf       : Fortran unit number for an open file
!!!   natom    : total number of solvent sites
!!!   nspecies : number of molecular species
!!!   nr       : number of grid points in Xvv
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readxvv1(this,nf,natom,nspecies,nr)
    use rism_parm
    implicit none
    type(rism3d_solv), intent(inout) :: this
    integer, intent(in) :: nf
    integer, intent(out) :: natom, nspecies, nr
    integer :: i, iok
    character(len=80) fmt, filename
    character(len=80) fmtin,ifmt,afmt,rfmt,type

    inquire(unit=nf,name=filename)
    fmt = ''
    ifmt = '(10I8)'
    afmt = '(20A4)'
    rfmt = '(5E16.8)'
    call rism_parm_nxtsec_init()
    !     ----- READ THE POINTERS AND THE VERSION -----

    fmtin = ifmt
    type = 'POINTERS'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok,this%xvv_version)
    !check for a version flag
    if(iok==-1) call rism_report_error('Xvv file '//trim(filename)//&
         ': missing %VERSION line')
    !check the version number.  Acceptable values are 0000.001,
    !0001.000 and 0001.001
    if(this%xvv_version == 0.000d0 .or. &
         (this%xvv_version >= 0.0015d0 .and. this%xvv_version <= 0.9995) .or. &
         this%xvv_version > 1.0015d0)then
       call rism_report_error('Xvv file '//trim(filename)//&
            ': bad version number. Must be 0001.001, 0001.000 or 0000.001')
!!$    elseif(this%xvv_version < 1.0005d0)then
!!$       call rism_report_warn("(a,f7.3)","Xvv version: ",dble(this%xvv_version))
!!$       call rism_report_warn("Unable to calculate UC or PMV temperature derivatives")
!!$       if(this%xvv_version < 0.9995d0)then
!!$          call rism_report_warn("Unable to calculate energy/entropy decomposition")
!!$          call rism_report_warn("UC assumes pure water")
!!$       end if
    end if

    read(nf,fmt) NR,NATOM,NSPECIES

  end subroutine readxvv1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Reads the data from a Amber style Xvv file into a solvent object.  Does 
!!!necessary unit conversion. Number of sites and grid size should be set in the
!!!solvent object based on previous call to readxvv1().
!!!IN:
!!!   solv : solvent object
!!!   nf   : Fortran unit number of an open file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine readxvv2(this,nf)
    use constants, only : COULOMB_CONST_E, KB, BOLTZMANN, AVOGADRO
    use rism_parm
    implicit none
    type(rism3d_solv), intent(inout) :: this
    integer, intent(in) :: nf

    integer :: i, iok
    character(len=80) fmt, filename
    character(len=80) fmtin,ifmt,afmt,rfmt,type

    integer itab, itab1, ig,iv,imv, iv1, iv2, iiga, isp


    inquire(unit=nf,name=filename)
    fmt = ''

    ifmt = '(10I8)'
    afmt = '(20A4)'
    rfmt = '(5E16.8)'
    !     ----- READ THE SYMBOLS AND THE CHARGES AND THE MASSES -----

    fmtin = rfmt
    type = 'THERMO'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
    read(nf,fmt) this%temperature,this%dielconst,this%xappa,this%xikt,this%dr,this%smear

    fmtin = ifmt
    type = 'MTV'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
    read(nf,fmt) (this%mult(i),i = 1,this%natom)

    if(this%xvv_version>=1.)then
       fmtin = ifmt
       type = 'NVSP'
       call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
       read(nf,fmt) (this%natomspecies(i),i = 1,this%natom)
    else
       this%natomspecies=huge(1)
    end if

    fmtin = afmt
    type = 'ATOM_NAME'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
    read(nf,fmt) (this%atomname(i),i = 1,this%natom)

    fmtin = rfmt
    type = 'RHOV'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
    read(nf,fmt) (this%rho(i),i = 1,this%natom)

    fmtin = rfmt
    type = 'QSPV'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
    read(nf,fmt) (this%charge_sp(i),i = 1,this%natom)
    !UNIT CONVERSION 
    !from [e] to [sqrt(kT A)]
    if(this%xvv_version<1.)&
         this%charge_sp = this%charge_sp&
         *sqrt(COULOMB_CONST_E/(KB *this%temperature))

    fmtin = rfmt
    type = 'QV'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
    read(nf,fmt) (this%charge(i),i = 1,this%natom)
    !UNIT CONVERSION 
    !from [e] to [sqrt(kT A)]
    if(this%xvv_version<1.)&
         this%charge = this%charge&
         *sqrt(COULOMB_CONST_E/(KB *this%temperature))

    fmtin = rfmt
    type = 'EPSV'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
    read(nf,fmt) (this%eps(i),i = 1,this%natom)
    !UNIT CONVERSION
    !from [J/MOL] to [kT]
    !  this%eps=this%eps/JPKC
    if(this%xvv_version<1.)&
         this%eps=this%eps/ (BOLTZMANN * AVOGADRO *this%temperature) 


    fmtin = rfmt
    if(this%xvv_version<1.)then
       type = 'SIGV'
    else
       type = 'RMIN2V'
    end if
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
    read(nf,fmt) (this%sig_s(i),i = 1,this%natom)

    fmtin = rfmt
    type = 'DELHV0'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
    read(nf,fmt) (this%delhv0(i),i = 1,this%natom)
    !UNIT CONVERSION
    !from [e] to [sqrt(kT A)]
    if(this%xvv_version<1.)&
         this%delhv0 = this%delhv0/sqrt(COULOMB_CONST_E/KB/this%temperature)

    fmtin = rfmt
    type = 'XVV'
    call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
    read(nf,fmt) (((this%xvv(itab,iv1,iv2),itab = 1,this%nr),iv1=1,this%natom),iv2=1,this%natom)

    !Fields for version 1.000
    if(this%xvv_version>=1.)then
       !Temperature derivative information is optional
       
       fmtin = rfmt
       type = 'DELHV0_DT'
       call rism_parm_nxtsec(nf, rism_report_getEUnit(), 1,fmtin,  type,  fmt,  iok)
       if(iok==0)then
          read(nf,fmt) (this%delhv0_dT(i),i = 1,this%natom)
       
          fmtin = rfmt
          type = 'XVV_DT'
          call rism_parm_nxtsec(nf, rism_report_getEUnit(), 1,fmtin,  type,  fmt,  iok)
          if(iok==0)then
             read(nf,fmt) (((this%xvv_dT(itab,iv1,iv2),itab = 1,this%nr),iv1=1,this%natom),iv2=1,this%natom)
          end if
       end if

       if(this%xvv_version>=1.0005)then
          fmtin = rfmt
          type = 'THERMO_DT'
          call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
          read(nf,fmt) this%xikt_dT

          fmtin = rfmt
          type = 'RHOSP'
          call rism_parm_nxtsec(nf, rism_report_getEUnit(), 0,fmtin,  type,  fmt,  iok)
          read(nf,fmt) this%rho_sp
       elseif(this%xvv_version>=1.)then
          !get RHOSP from other information
          do isp =1, this%nspecies
             this%rho_sp(isp) = this%rho(sum(this%natomspecies(1:isp)))/&
                  this%mult(sum(this%natomspecies(1:isp)))
          end do
       end if
    else
       this%rho_sp = 0d0
       this%rho_sp(1) = this%rho(1)
    end if
  end subroutine readxvv2
end module rism3d_solv_c
