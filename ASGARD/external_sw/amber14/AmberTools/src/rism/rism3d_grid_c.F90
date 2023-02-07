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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Grid class for 3D-RISM.  Handles the size and shape of the grid and
!!!manages precalculated tables of wavevectors.
!!!
!!!The class is not MPI aware but explicitly supports grid decomposition. The 
!!!size of the global grid (NG*) is stored as is the local grid (N*).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "../include/dprec.fh"

module rism3d_grid_c
  use safemem
  implicit none

  type rism3d_grid
       !grdspc   :: grid spacing for all of the grids [A]
       !boxlen   :: box size for 3d-rism.  for PBC calculations, these should 
       !            generally be equal [A] 
       !boxvol   :: volume of the box we are doing the calculation in [A^3]
       !voxel    :: volume of a grid point [A^3]
       _REAL_ :: grdspc(3)=0.5d0,boxlen(3),boxvol,voxel

       !ngr        :: number of global r-space grid points in each dimension
       !ngk        :: number of global k-space grid points in each dimension
       !nr         :: number of r-space grid points in each dimension
       !nk         :: number of k-space grid points in each dimension
       !nrOff      :: r-space offset from (0,0,0)
       !nkOff      :: k-space offset from (0,0,0)
       integer ::  ngr(3), ngk(3), nr(3), nk(3), nrOff(3), nkOff(3)
       
       !ngrTotal :: total number of global r-space grid points
       !nrTotal :: total number of local r-space grid points
       !nkTotal :: total number of local k-space grid points
       integer :: ngrTotal, nrTotal, nkTotal

       !nga        :: number of G indices
       integer :: nga
       !ga         :: absolute value of G (THIS%NGA long)
       !gv         :: G
       !g2         :: G^2
       _REAL_,pointer :: ga(:)=>NULL(),gv(:,:)=>NULL(),g2(:)=>NULL()

       !indga      :: index of |G| in accending order of |G|
       integer,pointer :: indga(:)=>NULL()

       !mpirank :: MPI rank number
       !mpisize :: Number of MPI processes
       !mpicomm :: MPI communicator
       integer :: mpirank, mpisize, mpicomm
  end type rism3d_grid
private setup_wavevector
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Constructor for a new rism3d_grid
!!!IN:
!!!   this : grid object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_grid_new(this,mpicomm)
    implicit none
    type(rism3d_grid), intent(inout) :: this
    integer, optional, intent(in) :: mpicomm
    if(present(mpicomm))then
       call rism3d_grid_setmpi(this,mpicomm)
    else
       call rism3d_grid_setmpi(this,0)
    end if
  end subroutine rism3d_grid_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!sets the MPI communicator.  If compiled without MPI, rank=0 and size=1
!!!IN:
!!!   this  : grid object
!!!   comm  : MPI communicator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_grid_setmpi(this,comm)
    implicit none
#if defined(MPI)
    include 'mpif.h'
#endif
    type(rism3d_grid), intent(inout) :: this
    integer, intent(in) :: comm
    integer :: err
       this%mpicomm = comm
       this%mpirank = 0
       this%mpisize = 1
#if defined(MPI)
       if(this%mpicomm ==MPI_COMM_NULL)&
            call rism_report_error("RISM3D_GRID: received NULL MPI communicator")
       call mpi_comm_rank(this%mpicomm,this%mpirank,err)
       if(err /=0) call rism_report_error&
            ("(a,i8)","RISM3D GRID: could not get MPI rank for communicator ",comm)
       call mpi_comm_size(this%mpicomm,this%mpisize,err)
       if(err /=0) call rism_report_error&
            ("(a,i8)","RISM3D GRID: could not get MPI size for communicator ",comm)
#endif /*defined(MPI)*/

  end subroutine rism3d_grid_setmpi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!sets the grid spacing for future calculations
!!!IN:
!!!   this  : grid object
!!!   grdpsc: grid spacing in each dimension [A]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_grid_setSpacing(this,grdspc)
    implicit none
    type(rism3d_grid), intent(inout) :: this
    _REAL_, intent(in) :: grdspc(3)

    this%grdspc = grdspc
  end subroutine rism3d_grid_setSpacing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Constructor for a new rism3d_grid
!!!IN:
!!!   this  : grid object
!!!   grdpsc: grid spacing in each dimension [A]
!!!   ngr   : number of global r-space grid points in each dimension
!!!   ngk   : number of global k-space grid points in each dimension
!!!   nr    : number of r-space grid points in each dimension
!!!   nk    : number of k-space grid points in each dimension
!!!   nrOff : r-space offset from (0,0,0)
!!!   nkOff : k-space offset from (0,0,0)
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_grid_resize(this,grdspc,ngr,ngk,nr,nk,nrOff,nkOff)
    implicit none
    type(rism3d_grid), intent(inout) :: this
    integer, intent(in) ::  ngr(3), ngk(3), nr(3), nk(3), nrOff(3), nkOff(3)
    _REAL_, optional, intent(in) :: grdspc(3)

    this%grdspc = grdspc
   
    this%ngr=ngr
    this%ngk=ngk
    this%nr =nr
    this%nk =nk
    this%nrOff=nrOff
    this%nkOff=nkOff

    this%ngrTotal=product(this%ngr)
    this%nrTotal=product(this%nr)
    this%nkTotal=product(this%nk)

    this%boxlen = this%ngr*this%grdspc
    this%boxvol = product(this%boxlen)
    this%voxel = product(this%grdspc)

    this%indga => safemem_realloc(this%indga,this%nkTotal/2,.false.)
    this%g2 => safemem_realloc(this%g2,this%nkTotal/2,.false.)
    this%gv => safemem_realloc(this%gv,3,this%nkTotal/2,.false.)
    
    call setup_wavevector(this)
  end subroutine rism3d_grid_resize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Destroys a rism3d_grid object
!!!IN:
!!!   this : grid object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_grid_destroy(this)
    implicit none
    type(rism3d_grid), intent(inout) :: this

    this%mpirank = 0
    this%mpisize = 0
    this%mpicomm = 0
    if(safemem_dealloc(this%ga) /=0) then
       call rism_report_error("RISM3D_GRID: failed to deallocate ga")
    end if
    if(safemem_dealloc(this%gv) /=0) then
       call rism_report_error("RISM3D_GRID: failed to deallocate gv")
    end if
    if(safemem_dealloc(this%g2) /=0) then
       call rism_report_error("RISM3D_GRID: failed to deallocate g2")
    end if
    if(safemem_dealloc(this%indga) /=0) then
       call rism_report_error("RISM3D_GRID: failed to deallocate indga")
    end if
  end subroutine rism3d_grid_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                             PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Does first round of set up after a new solvation box size has been set.
!!!Specifically the code
!!! -Calculates the Wavevector G, G^2,                   
!!! -Calculates, Sorts and Indexes Ga=|G|,
!!! -allocates memory for ga and this%xvva
!!!IN:
!!!   this     :: rism3d_grid object    
!!!MODIFIED:
!!!   gv, g2, ga, this%indga, this%nga, this%xvva
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine setup_wavevector (this)
      use constants, only : PI
      use safemem
      use rism_util, only: index_array,checksum
      implicit none
      type(rism3d_grid), intent(inout) :: this

      integer ::  ngx,ngy,ngz, lgx,lgy,lgz, igx,igy,igz, igk,igs,&
           ndimx,ndimy,ndimz
      _REAL_ ::  gx,gy,gz, ga2o,ga2
      integer,pointer :: indg2(:)=>NULL()

      logical :: opened
      integer :: irank, ierr

      indg2 => safemem_realloc(indg2,this%nkTotal/2,.false.)

      !................ geting and checking box arrays sizes .................
      ngx = this%ngr(1)
      ngy = this%ngr(2)
      ngz = this%ngr(3)
#if defined(MPI)
      ndimx = this%nk(1)
      ndimy = this%nk(2)
      ndimz = this%nk(3)
      !.......... initializing G for Real FFT in wrap-around order ...........
      do igz=0,ndimz-1
         do igy=0,ndimy-1
            do igx=0,ndimx/2-1
               igk = 1 + igx + igy*(ndimx/2) + igz*ndimy*(ndimx/2)

               lgx = mod( igx+              (ngx/2-1), ngx) - (ngx/2-1)
               lgy = mod( igz+this%nkOff(3)+(ngy/2-1), ngy) - (ngy/2-1)
               lgz = mod( igy+              (ngz/2-1), ngz) - (ngz/2-1)

               !.................... getting and loading G and G^2 ....................
               gx = 2d0*PI/this%boxlen(1) * lgx
               gy = 2d0*PI/this%boxlen(2) * lgy
               gz = 2d0*PI/this%boxlen(3) * lgz
               this%gv(1,igk) = gx
               this%gv(2,igk) = gy
               this%gv(3,igk) = gz
               this%g2(igk) = gx**2 + gy**2 + gz**2
            enddo
         enddo
      enddo
#else
      !.......... initializing G for Real FFT in wrap-around order ...........
      do igz=0,ngz-1
         do igy=0,ngy-1
            do igx=0,ngx/2-1
               igk = 1 + igx + igy*ngx/2 + igz*ngy*ngx/2

               lgx = mod( igx+(ngx/2-1), ngx) - (ngx/2-1)
               lgy = mod( igy+(ngy/2-1), ngy) - (ngy/2-1)
               lgz = mod( igz+(ngz/2-1), ngz) - (ngz/2-1)

               !.................... getting and loading G and G^2 ....................
               gx = 2d0*PI/this%boxlen(1) * lgx
               gy = 2d0*PI/this%boxlen(2) * lgy
               gz = 2d0*PI/this%boxlen(3) * lgz
               this%gv(1,igk) = gx
               this%gv(2,igk) = gy
               this%gv(3,igk) = gz
               this%g2(igk) = gx**2 + gy**2 + gz**2
            enddo
         enddo
      enddo
      !.......... initializing Nyquist frequency of G for Real FFT ...........
      do igz=0,ngz-1
         do igy=0,ngy-1
            igx = ngx/2
            igk = ngx/2*ngy*ngz + 1 + igy + igz*ngy
            lgx = mod( igx+(ngx/2-1), ngx) - (ngx/2-1)
            lgy = mod( igy+(ngy/2-1), ngy) - (ngy/2-1)
            lgz = mod( igz+(ngz/2-1), ngz) - (ngz/2-1)

            !.................... getting and loading G and G^2 ....................
            gx = 2d0*PI/this%boxlen(1) * lgx
            gy = 2d0*PI/this%boxlen(2) * lgy
            gz = 2d0*PI/this%boxlen(3) * lgz
            this%gv(1,igk) = gx
            this%gv(2,igk) = gy
            this%gv(3,igk) = gz
            this%g2(igk) = gx**2 + gy**2 + gz**2
         enddo
      enddo
#endif /*defined(MPI)*/

      !............. sorting box G^2 by ascending absolute value .............
      call index_array (this%g2,indg2, this%nkTotal/2)
      !............. selecting distinct values and loading Gabs ..............
      !start ga with a size of 1
      this%ga => safemem_realloc(this%ga,1,.false.)
      
      ga2o = -1d0
      this%nga = 0
      do igk=1,this%nkTotal/2
         igs = indg2(igk)
         ga2 = this%g2(igs)
         if (ga2 > ga2o)  then
            this%nga = this%nga + 1
            !bump up the memory.  We will reduce it later if it is too big
            if (this%nga > size(this%ga))  then
               this%ga => safemem_realloc(this%ga,4*this%nga/3,.true.)
            endif
            this%ga(this%nga) = sqrt( ga2)
            ga2o = ga2
         endif
         this%indga(igs) = this%nga
      enddo

      !fix memory size here
      this%ga => safemem_realloc(this%ga,this%nga,.true.)
      ierr=safemem_dealloc(indg2)
    end subroutine setup_wavevector
end module rism3d_grid_c

