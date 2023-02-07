!<compile=optimized>

!The 3D-RISM-KH software found here is copyright (c) 2010 by Andriy Kovalenko,
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
!
!4) I. Omelyan and A. Kovalenko, J. Chem. Phys., 139, 244106 (2013).

#include "../include/dprec.fh"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Object for force-coordinate-extrapolation (FCE) multiple time step (MTS).
!!!
!!!N total previous frames are held in memory for both position and
!!!solvation force.  When we are extrapolating the force for a given
!!!time step we first express the current position as a linear
!!!combination of the previous steps:
!!!
!!!  R^(N+1) ~ \sum^N_i a_i R^i
!!!
!!!the a_i that best satisfy this equation are then used to calculate
!!!an approximation of the solvation:
!!!
!!!  F^(N+1) ~ \sum^N_i a_i F^i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module fce_c
  implicit none
  type fce
     sequence
     !enumerated values for coordinate type to use as the basis for extrapolation
     integer :: CRDPOS=0, CRDDIST=1, CRDXYZ=2

     !nbasis  :: total number of basic coordinate-force pairs used for interpolation
     !nbase   :: number of leading FCE basis vectors

     ! Note that if nbasis > nbase, then selecting of the best nbase points
     ! among all nbasis coordinates is carried out

     !crd     :: 0-position, 1-distance, 2-xyz distance

     !weigh   ::  0-usual (default) or /=0-weighted coordinate minimization
     !            deviations [expensive but more precise]

     !sort    ::  0-no sorting
     !        ::/=0-permorm sorting even through nbasis=nbase or additionally
     !              to selecting if nbasis>nbase 

     ! It should be pointed out that contrary to selecting, the sorting
     ! is not so important because it should lead to the same results as
     ! for no sorting, provided the round-off errors are neglected

     ! Note that no weighting (weigh=0) are performed for trans=0 but
     ! selecting (if nbasis>nbase) and sorting (if sort/=0) are possible
      
     !enormsw ::  0-no minimization of the norms of the solutions
     !        ::/=0-minimization with enormsw-weight

     !trans ::
     !
     ! 0-no coordinate transformation at force extrapolation (very fast
     !   but not precise, can be used only for small outer time steps
     !   which do not exceed 200 fs) [default, not recommended for
     !   larger spacing]
     !
     ! 1-transformation with respect to the first basic point, while
     !   selecting according to the position of the current coordinate
     !   [recommended for large steps up to order of several picoseconds,
     !    because fast and precise]
     !
     ! 2-transformation with respect to the first basic point, while
     !   selecting according to the position of the current coordinate,
     !   i.e. like trans=1, but now using the normal equations method
     !   for the linear least squares coordinate deviation minimization
     !   instead of the QR decomposition approach as in the cases
     !   trans=0,1, and 3. Moreover, the normal equations method is
     !   is complemeted here by the weight minimization of the norms
     !   of the solutions. Note that if the weight of such additional
     !   minimization enormsw=0 then trans=2 is mathematically equivalent
     !   to trans=1. However, an extra precision can be reached when
     !   enormsw accepts small positive values. Note also that within
     !   the QR decomposition method an additional norm solution
     !   minimization is possible only when fcenbase is larger than
     !   the number of the solute degrees of freedom. The case trans=3
     !   is recommended (in conjunction with weigh=0, see below, to reduce
     !   costs) especcialy for huge outer steps and moderate fcenbase
     !
     ! 3-transformation and selecting with respect to the current coordinate
     !   [most precise but time expensive, not recommended for large nbasis]
     !
     ! 4-no normalization, coordinate transformation, weighting, selecting,
     !   and sorting, but with individual force extrapolation and possible
     !   cutting-off the neighbours relatively to each current atom [default
     !   for the original AMBER11 version]
     !

     !nsample :: number of samples collected
     !natom   :: number of atoms
     integer :: nbasis=-1, crd=-1, nsample=0, natom=-1

     !mpirank - MPI rank
     !mpicomm - MPI mpicomm
     !mpisize - number of MPI processes
     !atom0 - first atom in MPI decomposition
     !atomF - final atom in MPI decomposition
     integer :: mpirank=0, mpicomm=0, mpisize=1
     integer :: atom0, atomf

     !cut :: distance cutoff used for creating the basis set
     _REAL_ :: cut

     !cutlist :: for each atoms, list of atoms to include in force extrapolation (natom+1,natom).
     !           the first row element indicates how many atoms follow in this row.
     integer,pointer :: cutlist(:,:) => NULL()
     !refit :: TBA
     logical,pointer :: refit(:) => NULL()

     !force :: previous forces (dimensions:natom:nbasis)
     !coord :: previous coordinates (dimensions:natom:nbasis)
     _REAL_,pointer :: force(:,:,:) => NULL(), coord(:,:,:) => NULL()

     !sforce :: transformed forces (dimensions:natom:nbasis)
     !scoord :: transformed coordinates (dimensions:natom:nbasis)
     _REAL_,pointer :: sforce(:,:,:) => NULL(), scoord(:,:,:) => NULL()

     !! coeff :: the a_i coefficients in the above equation
     _REAL_, pointer :: coeff(:,:) => NULL()

     integer :: nbase=-1, weigh=-1, trans=-1, sort=-1
    _REAL_   :: enormsw=-1.d0

  end type fce

  interface fce_new
     module procedure fce_new_std
  end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!public subroutines and functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
public fce_new, fce_destroy, fce_update, fce_forcea, fce_forceb, fce_forcebm, &
       fce_forcec, fce_force, fce_estimate, fce_transformi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!private subroutines and functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
private nlist

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!public subroutines and functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Constructor. If the this is MPI, then only the rank 0 process parameters will
!!!be used.
!!!IN:
!!!   this :: FCE object
!!!   natom :: number of atoms
!!!   nbasis :: number of all basis vectors
!!!   nbase  :: number of leading basis vectors
!!!   weigh  :: weighted minimization flag
!!!   trans  :: transforming coordinate flag
!!!   sort   :: sorting flag
!!!   enormsw:: solution norm minimization weight
!!!   crd :: coordinate orientation method. 0-position, 1-distance, 2-xyz distance
!!!   cut :: cut off distance for dependence
!!!   o_mpicomm :: (optional) MPI communicator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fce_new_std(this,natom,nbasis,nbase,crd,weigh,trans,sort,enormsw,cut,o_mpicomm)
    use safemem
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif /*MPI*/
    type(fce),intent(inout) :: this
    integer, intent(in) :: natom, nbasis, crd
    integer, intent(in) :: nbase, weigh, trans, sort
    _REAL_, intent(in) :: cut, enormsw
    integer, optional, intent(in) :: o_mpicomm
    integer :: err
#ifdef MPI
    this%mpicomm = 0
    if(present(o_mpicomm)) this%mpicomm = o_mpicomm
    if(this%mpicomm == MPI_COMM_NULL)&
       call rism_report_error("FCE: received NULL MPI communicator")
    call mpi_comm_rank(this%mpicomm,this%mpirank,err)
    if(err /=0) call rism_report_error&
         ("(a,i8)","FCE: could not get MPI rank for communicator ",this%mpicomm)
    call mpi_comm_size(this%mpicomm,this%mpisize,err)
    if(err /=0) call rism_report_error&
         ("(a,i8)","FCE: could not get MPI size for communicator ",this%mpicomm)
#endif /*MPI*/
    if(this%mpirank==0)then
       this%natom = natom
       this%nbasis = nbasis

       this%nbase = nbase
       this%weigh = weigh
       this%trans = trans
       this%sort  = sort
       this%enormsw = enormsw

       this%crd = crd
       this%cut = cut**2

       this%cutlist => safemem_realloc(this%cutlist,natom,natom,.false.)
       this%force => safemem_realloc(this%force,3,natom,nbasis,.false.)
       this%coord => safemem_realloc(this%coord,3,natom,nbasis,.false.)

       this%sforce => safemem_realloc(this%sforce,3,natom,nbasis,.false.)
       this%scoord => safemem_realloc(this%scoord,3,natom,nbasis,.false.)

       this%coeff => safemem_realloc(this%coeff,3,natom,.false.)
    end if

#ifdef MPI
    !first distribute the pieces of information needed to allocate memory
    call mpi_bcast(this%natom,1,mpi_integer,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast NATOM")
    call mpi_bcast(this%nbasis,1,mpi_integer,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast NBASIS")

    call mpi_bcast(this%nbase,1,mpi_integer,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast NBASE")

    call mpi_bcast(this%weigh,1,mpi_integer,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast WEIGH")

    call mpi_bcast(this%trans,1,mpi_integer,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast TRANS")

    call mpi_bcast(this%sort,1,mpi_integer,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast SORT")

    call mpi_bcast(this%enormsw,1,mpi_double_precision,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast ENORMSW")

    call mpi_bcast(this%crd,1,mpi_integer,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast CRD")
    call mpi_bcast(this%cut,1,mpi_double_precision,0,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not broadcast CUT")

    !non-master processes should now allocate memory
    if(this%mpirank /= 0) then
       this%cutlist => safemem_realloc(this%cutlist,this%natom,this%natom,.false.)
       this%force => safemem_realloc(this%force,3,this%natom,this%nbasis,.false.)
       this%coord => safemem_realloc(this%coord,3,this%natom,this%nbasis,.false.)

       this%sforce => safemem_realloc(this%sforce,3,this%natom,this%nbasis,.false.)
       this%scoord => safemem_realloc(this%scoord,3,this%natom,this%nbasis,.false.)

       this%coeff => safemem_realloc(this%coeff,3,this%natom,.false.)
    end if

    !Arrays contain no data yet so there is nothing to transfer

#endif /*MPI*/

    !set local atom range for this process
    this%atom0 = this%natom/this%mpisize*this%mpirank+1
    this%atomF = min(int(this%natom*dble(this%mpirank+1)/dble(this%mpisize)),this%natom)

  end subroutine fce_new_std

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Destroyer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fce_destroy(this)
    use safemem
    implicit none
    type(fce),intent(inout) :: this
    integer :: err
    this%natom = -1
    this%nbasis = -1

    this%nbase = -1
    this%weigh = -1
    this%trans = -1
    this%sort  = -1
    this%enormsw = -1.d0

    this%crd = -1

    err = safemem_dealloc(this%cutlist)
    err = safemem_dealloc(this%force)
    err = safemem_dealloc(this%coord)
    err = safemem_dealloc(this%coeff)

    err = safemem_dealloc(this%sforce)
    err = safemem_dealloc(this%scoord)

  end subroutine fce_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!update the force and coordinate basis vectors
!!!IN:
!!!   this  : FCE object
!!!   force : forces on atoms.  In the case of MPI, it is assumed that the forces
!!!           are distributed across the processes and must be reduced internally
!!!   coord : coordinates of the atoms, assume to the the same for all processes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fce_update(this, force, coord)
    use safemem
    implicit none
#if defined(MPI)
      include 'mpif.h'
#endif /*defined(MPI)*/
    type(fce), intent(inout) :: this
    _REAL_, intent(in) :: force(3,this%natom), coord(3,this%natom)

    !workspace for the case of MPI 1.1
#ifdef MPI
    integer :: err
#  ifndef USE_MPI_IN_PLACE
    _REAL_, pointer :: tforce(:,:)=>NULL()
#  endif /*ifdef USE_MPI_IN_PLACE*/
#endif /*MPI*/

#ifdef RISM_DEBUG
    write(6,*) "FCE_UPDATE"; call flush(6)
    write(6,*) "FCE", this%atom0, this%atomF, this%mpisize, this%mpirank, this%mpicomm, this%cut
#endif /*RISM_DEBUG*/

    !
    !For now we will shift the data through the storage arrays as more data is
    !added.  Later, this should be modified to just update a pointer to the most
    !recent entry
    !
    if(this%nsample > 0)then
       this%coord(:,:,2:min(this%nsample+1,this%nbasis)) = this%coord(:,:,1:min(this%nsample, this%nbasis-1))
       this%force(:,:,2:min(this%nsample+1,this%nbasis)) = this%force(:,:,1:min(this%nsample, this%nbasis-1))
    end if
    this%coord(:,:,1) = coord
    this%force(:,:,1) = force
    !
    !reduce the atoms locally
    !
#ifdef MPI
#  ifdef USE_MPI_IN_PLACE
    call mpi_allreduce(MPI_IN_PLACE,this%force(:,:,1),3*this%natom,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not reduce force")
#  else /*ifdef USE_MPI_IN_PLACE*/
    tforce => safemem_realloc(tforce,3,this%natom,.false.)
    call mpi_allreduce(this%force(:,:,1),tforce,3*this%natom,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
    if(err /=0) call rism_report_error&
         ("FCE: could not reduce force")
    this%force(:,:,1) = tforce
    if(safemem_dealloc(tforce) /=0) call rism_report_error("FCE_UPDATE: deallocate TFORCE failed")
#  endif /*ifdef USE_MPI_IN_PLACE*/
#endif /*MPI*/

    this%nsample = min(this%nsample+1,this%nbasis)

  end subroutine fce_update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Transform the force and coordinate basis vectors
!!!IN:
!!!   this  : FCE object
!!!   force : forces on atoms
!!!   coord : coordinates of the atoms
!!!OUT:
!!!  sforce : transformed forces
!!!  scoord : transformed coordinates
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fce_transformi(this,force,coord,sforce,scoord)

    use safemem
    implicit none

#ifdef MPI
      include 'mpif.h'
#endif

    type(fce), intent(inout) :: this
    _REAL_, intent(in) :: force(3,this%natom), coord(3,this%natom)
    _REAL_ :: sforce(3,this%natom), scoord(3,this%natom)

    integer :: lswork, infos
   _REAL_, pointer :: swork(:)=>NULL()

    DOUBLE PRECISION hmu(4,4),ddd,q0,q1,q2,q3,fis(3,3),wjjj(4)
    DOUBLE PRECISION xiii1,xiii2,xiii3,fiii1,fiii2,fiii3

     integer :: iatm, isa

#ifdef MPI
     integer :: err
     DOUBLE PRECISION hmus(this%nbasis,4,4)
#endif

! The transformation is carried out with respect to the first pair

! The coordinates are known for all the processors
       do iatm = 1, this%natom
    this%scoord(1,iatm,1)=this%coord(1,iatm,1)
    this%scoord(2,iatm,1)=this%coord(2,iatm,1)
    this%scoord(3,iatm,1)=this%coord(3,iatm,1)
       end do

! The forces are distributed across the processors
       do iatm = this%atom0, this%atomF
    this%sforce(1,iatm,1)=this%force(1,iatm,1)
    this%sforce(2,iatm,1)=this%force(2,iatm,1)
    this%sforce(3,iatm,1)=this%force(3,iatm,1)
       end do

! Find optimal rotation of a solute as a whole for the next pairs to
! provide the minimum for the distances between the basic coordinates

#ifdef MPI

! Summing up the quaternion matrix over atoms across the processes

       hmus=0.d0
       do iatm = this%atom0, this%atomF
       fiii1=this%coord(1,iatm,1)
       fiii2=this%coord(2,iatm,1)
       fiii3=this%coord(3,iatm,1)
       do isa=2,this%nsample
       xiii1=this%coord(1,iatm,isa)
       xiii2=this%coord(2,iatm,isa)
       xiii3=this%coord(3,iatm,isa)
       ddd=(xiii1+fiii1)**2+(xiii2+fiii2)**2+(xiii3+fiii3)**2
       hmus(isa,1,1)=hmus(isa,1,1)+((xiii1-fiii1)**2+(xiii2-fiii2)**2+ &
                          (xiii3-fiii3)**2)
       hmus(isa,1,2)=hmus(isa,1,2)+2.d0*(xiii2*fiii3-xiii3*fiii2)
       hmus(isa,1,3)=hmus(isa,1,3)+2.d0*(xiii3*fiii1-xiii1*fiii3)
       hmus(isa,1,4)=hmus(isa,1,4)+2.d0*(xiii1*fiii2-xiii2*fiii1)
       hmus(isa,2,2)=hmus(isa,2,2)+(ddd-2.d0*(xiii1*fiii1+fiii1*xiii1))
       hmus(isa,2,3)=hmus(isa,2,3)     -2.d0*(xiii1*fiii2+fiii1*xiii2)
       hmus(isa,2,4)=hmus(isa,2,4)     -2.d0*(xiii1*fiii3+fiii1*xiii3)
       hmus(isa,3,3)=hmus(isa,3,3)+(ddd-2.d0*(xiii2*fiii2+fiii2*xiii2))
       hmus(isa,3,4)=hmus(isa,3,4)     -2.d0*(xiii2*fiii3+fiii2*xiii3)
       hmus(isa,4,4)=hmus(isa,4,4)+(ddd-2.d0*(xiii3*fiii3+fiii3*xiii3))
       end do
       end do

! Collecting the partial sums from different processes and placing the
! result into the same quaternion matrix known for all processors

    call mpi_allreduce(MPI_IN_PLACE,hmus,this%nbasis*16,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)

       do isa=2,this%nsample

       hmu(:,:)=hmus(isa,:,:)

#else

! Summing up the quaternion matrix over atoms in the single process mode

       do isa=2,this%nsample
       hmu=0.d0
       do iatm = 1, this%natom
       fiii1=this%coord(1,iatm,1)
       fiii2=this%coord(2,iatm,1)
       fiii3=this%coord(3,iatm,1)
       xiii1=this%coord(1,iatm,isa)
       xiii2=this%coord(2,iatm,isa)
       xiii3=this%coord(3,iatm,isa)
       ddd=(xiii1+fiii1)**2+(xiii2+fiii2)**2+(xiii3+fiii3)**2
       hmu(1,1)=hmu(1,1)+((xiii1-fiii1)**2+(xiii2-fiii2)**2+ &
                          (xiii3-fiii3)**2)
       hmu(1,2)=hmu(1,2)+2.d0*(xiii2*fiii3-xiii3*fiii2)
       hmu(1,3)=hmu(1,3)+2.d0*(xiii3*fiii1-xiii1*fiii3)
       hmu(1,4)=hmu(1,4)+2.d0*(xiii1*fiii2-xiii2*fiii1)
       hmu(2,2)=hmu(2,2)+(ddd-2.d0*(xiii1*fiii1+fiii1*xiii1))
       hmu(2,3)=hmu(2,3)     -2.d0*(xiii1*fiii2+fiii1*xiii2)
       hmu(2,4)=hmu(2,4)     -2.d0*(xiii1*fiii3+fiii1*xiii3)
       hmu(3,3)=hmu(3,3)+(ddd-2.d0*(xiii2*fiii2+fiii2*xiii2))
       hmu(3,4)=hmu(3,4)     -2.d0*(xiii2*fiii3+fiii2*xiii3)
       hmu(4,4)=hmu(4,4)+(ddd-2.d0*(xiii3*fiii3+fiii3*xiii3))
       end do

#endif

! The problem is reduced to find eigenvalues and eigenvectors
! of the corresponding quaternion symmetric matrices

! First query the optimal workspace

       lswork = -1
       swork => safemem_realloc(swork,1,.false.)
       CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )
       lswork = INT(swork(1))
       swork => safemem_realloc(swork,lswork,.false.)

! Now solve the eigenproblem

       CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )

       if(infos.ne.0) then
       print *,infos
       stop 'Message from DSYEV'
       end if

! Define the 4-dimensional quaternion corresponding to the smallest eigenvalue

       q0=hmu(1,1)
       q1=hmu(2,1)
       q2=hmu(3,1)
       q3=hmu(4,1)

! Determine the rotational matrix in 3-dimensional space

       fis(1,1)=q0**2+q1**2-q2**2-q3**2
       fis(1,2)=2.d0*(-q0*q3+q1*q2)
       fis(1,3)=2.d0*( q0*q2+q1*q3)
       fis(2,1)=2.d0*( q0*q3+q1*q2)
       fis(2,2)=q0**2+q2**2-q1**2-q3**2
       fis(2,3)=2.d0*(-q0*q1+q2*q3)
       fis(3,1)=2.d0*(-q0*q2+q1*q3)
       fis(3,2)=2.d0*( q0*q1+q2*q3)
       fis(3,3)=q0**2+q3**2-q1**2-q2**2

! Perform the coordinate and force transformations

       do iatm = 1, this%natom
    this%scoord(1,iatm,isa)=fis(1,1)*this%coord(1,iatm,isa)+fis(2,1)*this%coord(2,iatm,isa)+fis(3,1)*this%coord(3,iatm,isa)
    this%scoord(2,iatm,isa)=fis(1,2)*this%coord(1,iatm,isa)+fis(2,2)*this%coord(2,iatm,isa)+fis(3,2)*this%coord(3,iatm,isa)
    this%scoord(3,iatm,isa)=fis(1,3)*this%coord(1,iatm,isa)+fis(2,3)*this%coord(2,iatm,isa)+fis(3,3)*this%coord(3,iatm,isa)
       end do

       do iatm = this%atom0, this%atomF
    this%sforce(1,iatm,isa)=fis(1,1)*this%force(1,iatm,isa)+fis(2,1)*this%force(2,iatm,isa)+fis(3,1)*this%force(3,iatm,isa)
    this%sforce(2,iatm,isa)=fis(1,2)*this%force(1,iatm,isa)+fis(2,2)*this%force(2,iatm,isa)+fis(3,2)*this%force(3,iatm,isa)
    this%sforce(3,iatm,isa)=fis(1,3)*this%force(1,iatm,isa)+fis(2,3)*this%force(2,iatm,isa)+fis(3,3)*this%force(3,iatm,isa)
       end do

       end do

  end subroutine fce_transformi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Using previous forces and coordinates, predicts solvation forces
!!! for the current set of coordinates.  Like  linprojpredict, except
!!! we are finding the forces on each atom individually and rotate the
!!! solute for each to optimize the prediction.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! 0-no transformation [not recommended] and weighting (weigh=0) but
     ! with possible selecting (if nbasis>nbase) and sorting (if sort/=0) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fce_forcea(this,force,coord)

    use safemem
    implicit none

#ifdef MPI
      include 'mpif.h'
#endif

    type(fce), intent(inout) :: this
    _REAL_, intent(out) :: force(3,this%natom)
    _REAL_, intent(in) :: coord(3,this%natom)

    integer :: iatm,id,err

    !LAPACK variables
    integer :: M, N, NN, NRHS=1, LDA, LDB, LWORK, RANK, INFO
    integer, pointer :: JPVT(:)=>NULL()

    !uncertainties of input coordinates in A and B are assumed
    !to be negligible small
    _REAL_ :: RCOND=0.d0
    _REAL_, pointer :: A(:,:)=>NULL(), B(:,:)=>NULL(), work(:)=>NULL()
    _REAL_ :: rrr(this%nbasis),rar(this%nbase)

    integer :: irr(this%nbasis),ird(this%nbase),irs(this%nbase)

    integer :: isa,jsa,is1,iisap1

     N = this%nsample

! No selecting and sorting if nbasis=nbase and sort=0

          if(this%nbasis.eq.this%nbase.and.this%sort.eq.0) then

          do isa=1,this%nbase
          ird(isa)=isa
          end do

          end if

! Calculate the distances between the current solute position and its
! counterparts from the extended basis set if nbasis > nbase or sort/=0 

          if(this%nbasis.gt.this%nbase.or.this%sort.ne.0) then

          do isa=1,N
          rrr(isa)=0.d0
          do iatm = this%atom0, this%atomF
          rrr(isa)=rrr(isa)+(coord(1,iatm)-this%coord(1,iatm,isa))**2+ &
                            (coord(2,iatm)-this%coord(2,iatm,isa))**2+ &
                            (coord(3,iatm)-this%coord(3,iatm,isa))**2
          end do
          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,rrr,N,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#endif

! Normalizing the distances

          do isa=1,N
          rrr(isa)=dsqrt(rrr(isa)/this%natom)
          end do

          end if

! Select the first nbase minimal values among all N=nsample elements
! if N > nbase (note that nbase <= nsample <= nbasis). In other words,
! choose the best force-coordinates pairs.

          if(this%nbasis.gt.this%nbase.and.this%sort.ne.0) then

          do isa=1,N
          irr(isa)=isa
          end do

          if(N.gt.this%nbase) then

          call selects(this%nbase,N,rrr,irr)

          if(this%nbase.ne.1) then

          do isa=1,this%nbase-1
          rar(isa)=rrr(isa)
          end do

! Now perform sorting of the selected minimal values in the descending order

          call hsort(this%nbase-1,rar,irs)

          end if

          irs(this%nbase)=this%nbase

          do isa=1,this%nbase
          ird(isa)=irr(irs(isa))
          end do

          end if

! On the beginning when N=nsample=nbase perform only sorting because
! then selecting is not necessary

          if(N.eq.this%nbase) then

          call hsort(N,rrr,ird)

          end if

          end if

! Only select the first nbase minimal values among all N=nsample elements
! if N > nbase without sorting (sort=0)

          if(this%nbasis.gt.this%nbase.and.this%sort.eq.0) then

          do isa=1,N
          irr(isa)=isa
          end do

          if(N.gt.this%nbase) call selects(this%nbase,N,rrr,irr)

          do isa=1,this%nbase
          ird(isa)=irr(isa)
          end do

          end if

! If sort/=0, perform sorting even for nbasis=nbase, i.e. when selecting
! is not necessary

          if(this%nbasis.eq.this%nbase.and.this%sort.ne.0) then

          call hsort(this%nbase,rrr,ird)

          end if

! The nearest pair

          is1=ird(1)

! Zero-order force approximation

          do iatm = this%atom0, this%atomF
          do id=1,3
          force(id,iatm) = this%force(id,iatm,is1)
          end do
          end do

          NN=this%nbase-1

          if(NN.gt.0) then

! High-order force extrapolation

          JPVT=>safemem_realloc(JPVT,NN,.false.)
          JPVT=0

! Solve the linear least-square problem

          M = 3 * this%natom
          LDA = M
          LDB = max(M,NN)

          A=>safemem_realloc(A,M,NN,.false.)
          B=>safemem_realloc(B,max(NN,M),1,.false.)

! Using multiprocessor technique to fill in the matrix elements

#ifdef MPI
          A=0.d0
          B=0.d0

          jsa=3*(this%atom0-1)
#else
          jsa=0
#endif

          do iatm = this%atom0, this%atomF
          B(jsa+1:jsa+3,1)=coord(:,iatm)-this%coord(:,iatm,is1)
          do isa=1,NN
          iisap1=ird(isa+1)
          A(jsa+1:jsa+3,isa)=this%coord(:,iatm,iisap1)-this%coord(:,iatm,is1)
          end do
          jsa=jsa+3
          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,A,M*NN,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
    call mpi_allreduce(MPI_IN_PLACE,B,LDB,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#endif

          !first determine the amount of memory required
          LWORK=-1
          WORK => safemem_realloc(work,1,.false.)
          call dgelsy(M,NN,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,LWORK,INFO)

          !perform the calculation
          LWORK=WORK(1)
          WORK => safemem_realloc(WORK,LWORK,.false.)
          call dgelsy(M,NN,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,LWORK,INFO)

          if(INFO.ne.0) then
          print *,INFO
          stop 'Message from DGELSY'
          end if

! Extrapolate the forces

          do iatm = this%atom0, this%atomF
          do isa=1,NN
          iisap1=ird(isa+1)
          do id=1,3
          force(id,iatm) = force(id,iatm)+&
          B(isa,1)*(this%force(id,iatm,iisap1)-this%force(id,iatm,is1))
          end do
          end do
          end do

          end if

    err = safemem_dealloc(A)
    err = safemem_dealloc(B)
    err = safemem_dealloc(WORK)
    err = safemem_dealloc(JPVT)

  end subroutine fce_forcea

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! 1-transformation with respect to the first basic point, while
     !   selecting according to the position of the current coordinate
     !   with sorting (if sort/=0) (recommended, because fast and precise)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fce_forceb(this,force,coord,sforce,scoord)

    use safemem
    implicit none

#ifdef MPI
      include 'mpif.h'
#endif

    type(fce), intent(inout) :: this
    _REAL_, intent(out) :: force(3,this%natom)
    _REAL_, intent(in) :: coord(3,this%natom)

    _REAL_, intent(out) :: sforce(3,this%natom)
    _REAL_, intent(out) :: scoord(3,this%natom)

    integer :: iatm,jatm,id,err

    !LAPACK variables
    integer :: M, N, NN, NRHS=1, LDA, LDB, LWORK, RANK, INFO
    integer, pointer :: JPVT(:)=>NULL()

    integer :: lswork, infos
   _REAL_, pointer :: swork(:)=>NULL()

    !uncertainties of input coordinates in A and B are assumed
    !to be negligible small
    _REAL_ :: RCOND=0.d0
    _REAL_, pointer :: A(:,:)=>NULL(), B(:,:)=>NULL(), work(:)=>NULL()
    _REAL_ :: rrr(this%nbasis),rar(this%nbase)

    integer :: irr(this%nbasis),ird(this%nbase),irs(this%nbase)

    integer :: isa,jsa,is1,iisap1

     DOUBLE PRECISION hmu(4,4),ddd,q0,q1,q2,q3,fis(3,3),wjjj(4),rweight
     DOUBLE PRECISION xiii1,xiii2,xiii3,fiii1,fiii2,fiii3,rwrwr2

#ifdef MPI
     integer inum,jsad(this%mpisize)
#endif

     N = this%nsample

! The transformation of the curent coordinate with respect to the first
! basic point is performed by minimization of the distance between them

       hmu=0.d0

! Summing up the quaternion matrix over atoms across the processes

       do iatm = this%atom0, this%atomF
       fiii1=this%coord(1,iatm,1)
       fiii2=this%coord(2,iatm,1)
       fiii3=this%coord(3,iatm,1)
       xiii1=coord(1,iatm)
       xiii2=coord(2,iatm)
       xiii3=coord(3,iatm)
       ddd=(xiii1+fiii1)**2+(xiii2+fiii2)**2+(xiii3+fiii3)**2
       hmu(1,1)=hmu(1,1)+((xiii1-fiii1)**2+(xiii2-fiii2)**2+ &
                          (xiii3-fiii3)**2)
       hmu(1,2)=hmu(1,2)+2.d0*(xiii2*fiii3-xiii3*fiii2)
       hmu(1,3)=hmu(1,3)+2.d0*(xiii3*fiii1-xiii1*fiii3)
       hmu(1,4)=hmu(1,4)+2.d0*(xiii1*fiii2-xiii2*fiii1)
       hmu(2,2)=hmu(2,2)+(ddd-2.d0*(xiii1*fiii1+fiii1*xiii1))
       hmu(2,3)=hmu(2,3)     -2.d0*(xiii1*fiii2+fiii1*xiii2)
       hmu(2,4)=hmu(2,4)     -2.d0*(xiii1*fiii3+fiii1*xiii3)
       hmu(3,3)=hmu(3,3)+(ddd-2.d0*(xiii2*fiii2+fiii2*xiii2))
       hmu(3,4)=hmu(3,4)     -2.d0*(xiii2*fiii3+fiii2*xiii3)
       hmu(4,4)=hmu(4,4)+(ddd-2.d0*(xiii3*fiii3+fiii3*xiii3))
       end do

#ifdef MPI

! Collecting the partial sums from different processes and placing the
! result into the same quaternion matrix known for all processors

    call mpi_allreduce(MPI_IN_PLACE,hmu,16,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)

#endif

! The problem is reduced to find eigenvalues and eigenvectors 
! of the corresponding quaternion symmetric matrix

! First query the optimal workspace

       lswork = -1
       swork => safemem_realloc(swork,1,.false.)
       CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )
       lswork = INT(swork(1))
       swork => safemem_realloc(swork,lswork,.false.)

! Now solve the eigenproblem

       CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )

       if(infos.ne.0) then
       print *,infos
       stop 'Message from DSYEV'
       end if

! Define the 4-dimensional quaternion corresponding to the smallest eigenvalue

       q0=hmu(1,1)
       q1=hmu(2,1)
       q2=hmu(3,1)
       q3=hmu(4,1)

! Determine the rotational matrix in 3-dimensional space

       fis(1,1)=q0**2+q1**2-q2**2-q3**2
       fis(1,2)=2.d0*(-q0*q3+q1*q2)
       fis(1,3)=2.d0*( q0*q2+q1*q3)
       fis(2,1)=2.d0*( q0*q3+q1*q2)
       fis(2,2)=q0**2+q2**2-q1**2-q3**2
       fis(2,3)=2.d0*(-q0*q1+q2*q3)
       fis(3,1)=2.d0*(-q0*q2+q1*q3)
       fis(3,2)=2.d0*( q0*q1+q2*q3)
       fis(3,3)=q0**2+q3**2-q1**2-q2**2

! Perform the transformations

       do iatm = 1, this%natom
   scoord(1,iatm)=fis(1,1)*coord(1,iatm)+fis(2,1)*coord(2,iatm)+fis(3,1)*coord(3,iatm)
   scoord(2,iatm)=fis(1,2)*coord(1,iatm)+fis(2,2)*coord(2,iatm)+fis(3,2)*coord(3,iatm)
   scoord(3,iatm)=fis(1,3)*coord(1,iatm)+fis(2,3)*coord(2,iatm)+fis(3,3)*coord(3,iatm)
       end do

! No selecting and sorting if nbasis=nbase and sort=0

          if(this%nbasis.eq.this%nbase.and.this%sort.eq.0) then

          do isa=1,this%nbase
          ird(isa)=isa
          end do

          end if

! Calculate the distances between the current solute position and its
! counterparts from the extended basis set if nbasis > nbase or sort/=0

          if(this%nbasis.gt.this%nbase.or.this%sort.ne.0) then

          do isa=1,N
          rrr(isa)=0.d0
          do iatm = this%atom0, this%atomF
          rrr(isa)=rrr(isa)+(scoord(1,iatm)-this%scoord(1,iatm,isa))**2+ &
                            (scoord(2,iatm)-this%scoord(2,iatm,isa))**2+ &
                            (scoord(3,iatm)-this%scoord(3,iatm,isa))**2
          end do
          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,rrr,N,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#endif

! Normalizing the distances

          do isa=1,N
          rrr(isa)=dsqrt(rrr(isa)/this%natom)
          end do

          end if

! Select the first nbase minimal values among all N=nsample elements
! if N > nbase (note that nbase <= nsample <= nbasis). In other words,
! choose the best transformed force-coordinates pairs.

          if(this%nbasis.gt.this%nbase.and.this%sort.ne.0) then

          do isa=1,N
          irr(isa)=isa
          end do

          if(N.gt.this%nbase) then

          call selects(this%nbase,N,rrr,irr)

          if(this%nbase.ne.1) then

          do isa=1,this%nbase-1
          rar(isa)=rrr(isa)
          end do

! Now perform sorting of the selected minimal values in the descending order

          call hsort(this%nbase-1,rar,irs)

          end if

          irs(this%nbase)=this%nbase

          do isa=1,this%nbase
          ird(isa)=irr(irs(isa))
          end do

          end if

! On the beginning when N=nsample=nbase perform only sorting because
! then selecting is not necessary

          if(N.eq.this%nbase) then

          call hsort(N,rrr,ird)

          end if

          end if

! Only select the first nbase minimal values among all N=nsample elements
! if N > nbase without sorting (sort=0)

          if(this%nbasis.gt.this%nbase.and.this%sort.eq.0) then

          do isa=1,N
          irr(isa)=isa
          end do

          if(N.gt.this%nbase) call selects(this%nbase,N,rrr,irr)

          do isa=1,this%nbase
          ird(isa)=irr(isa)
          end do

          end if

! If sort/=0, perform sorting even for nbasis=nbase, i.e. when selecting
! is not necessary

          if(this%nbasis.eq.this%nbase.and.this%sort.ne.0) then

          call hsort(this%nbase,rrr,ird)

          end if

! The nearest pair

          is1=ird(1)

! Zero-order force approximation

          do iatm = this%atom0, this%atomF
          do id=1,3
          sforce(id,iatm) = this%sforce(id,iatm,is1)
          end do
          end do

          NN=this%nbase-1

          if(NN.gt.0) then

! High-order force extrapolation

          JPVT=>safemem_realloc(JPVT,NN,.false.)
          JPVT=0

! Solve the linear least-square problem with or without weighting

          if(this%weigh.eq.0) then

          M = 3 * this%natom
          LDA = M
          LDB = max(M,NN)

          A=>safemem_realloc(A,M,NN,.false.)
          B=>safemem_realloc(B,max(NN,M),1,.false.)

! Using multiprocessor technique to fill in the matrix elements

#ifdef MPI
          A=0.d0
          B=0.d0

          jsa=3*(this%atom0-1)
#else
          jsa=0
#endif

          do iatm = this%atom0, this%atomF
          B(jsa+1:jsa+3,1)=scoord(:,iatm)-this%scoord(:,iatm,is1)
          do isa=1,NN
          iisap1=ird(isa+1)
          A(jsa+1:jsa+3,isa)=this%scoord(:,iatm,iisap1)-this%scoord(:,iatm,is1)
          end do
          jsa=jsa+3
          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,A,M*NN,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
    call mpi_allreduce(MPI_IN_PLACE,B,LDB,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#endif

          !first determine the amount of memory required
          LWORK=-1
          WORK => safemem_realloc(work,1,.false.)
          call dgelsy(M,NN,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,LWORK,INFO)

          !perform the calculation
          LWORK=WORK(1)
          WORK => safemem_realloc(WORK,LWORK,.false.)
          call dgelsy(M,NN,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,LWORK,INFO)

          if(INFO.ne.0) then
          print *,INFO
          stop 'Message from DGELSY'
          end if

          else

          jsa=0

#ifdef MPI

          jsad=0
          do iatm = this%atom0, this%atomF
          do jatm = iatm+1, this%natom
          rwrwr2=sum((scoord(1:3,iatm)-scoord(1:3,jatm))**2)
          if(rwrwr2 < this%cut) then
          jsa=jsa+3
          end if
          end do
          end do

    call MPI_Gather(jsa,1,MPI_INTEGER,jsad,1,MPI_INTEGER,0,this%mpicomm,err);

      if(this%mpirank == 0) then
      do inum=2,this%mpisize
      jsad(inum)=jsad(inum-1)+jsad(inum)
      end do
      end if

    call mpi_allreduce(MPI_IN_PLACE,jsad,this%mpisize,MPI_INTEGER,MPI_SUM,this%mpicomm,err)

      if(this%mpirank == 0) jsa=0
      do inum=1,this%mpisize-1
      if(this%mpirank == inum) jsa=jsad(inum)
      end do

          M=jsad(this%mpisize)
#else
          M = 3 * ((this%natom*(this%natom-1))/2)
#endif

          LDA = M
          LDB = max(M,NN)

          A=>safemem_realloc(A,M,NN,.false.)
          B=>safemem_realloc(B,max(NN,M),1,.false.)

          A=0.d0
          B=0.d0

          do iatm = this%atom0, this%atomF
          do jatm = iatm+1, this%natom

          rwrwr2=sum((scoord(1:3,iatm)-scoord(1:3,jatm))**2)

! Use cutting for weighting

          if(rwrwr2 < this%cut) then

          rweight=1.d0/dsqrt(rwrwr2)

          B(jsa+1:jsa+3,1) = rweight*(&
                             scoord(:,iatm)-this%scoord(:,iatm,is1) &
                            -scoord(:,jatm)+this%scoord(:,jatm,is1))
          do isa=1,NN
          iisap1=ird(isa+1)
          A(jsa+1:jsa+3,isa) = rweight*(&
                          this%scoord(:,iatm,iisap1)-this%scoord(:,iatm,is1) &
                         -this%scoord(:,jatm,iisap1)+this%scoord(:,jatm,is1))
          end do
          jsa=jsa+3
          end if
          end do
          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,A,M*NN,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
    call mpi_allreduce(MPI_IN_PLACE,B,LDB,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#else
          M=jsa
#endif

          !first determine the amount of memory required
          LWORK=-1
          WORK => safemem_realloc(work,1,.false.)
          call dgelsy(M,NN,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,LWORK,INFO)

          !perform the calculation
          LWORK=WORK(1)
          WORK => safemem_realloc(WORK,LWORK,.false.)
          call dgelsy(M,NN,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,LWORK,INFO)

          if(INFO.ne.0) then
          print *,INFO
          stop 'Message from DGELSY'
          end if

          end if

! Extrapolate the forces in the transformed space

          do iatm = this%atom0, this%atomF
          do isa=1,NN
          iisap1=ird(isa+1)
          do id=1,3
          sforce(id,iatm) = sforce(id,iatm)+&
          B(isa,1)*(this%sforce(id,iatm,iisap1)-this%sforce(id,iatm,is1))
          end do
          end do
          end do

          end if

! Performing the inverse transformation to obtain the extrapolated forces
! in the usual coordinates

          do iatm = this%atom0, this%atomF
    force(1,iatm)=fis(1,1)*sforce(1,iatm)+fis(1,2)*sforce(2,iatm)+fis(1,3)*sforce(3,iatm)
    force(2,iatm)=fis(2,1)*sforce(1,iatm)+fis(2,2)*sforce(2,iatm)+fis(2,3)*sforce(3,iatm)
    force(3,iatm)=fis(3,1)*sforce(1,iatm)+fis(3,2)*sforce(2,iatm)+fis(3,3)*sforce(3,iatm)
          end do

    err = safemem_dealloc(A)
    err = safemem_dealloc(B)
    err = safemem_dealloc(WORK)
    err = safemem_dealloc(JPVT)
    err = safemem_dealloc(swork)

  end subroutine fce_forceb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! 2-like trans=1 but within the normal equations method complemeted
     !   by the weight minimization of the norms of the solutions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fce_forcebm(this,force,coord,sforce,scoord)

    use safemem
    implicit none

#ifdef MPI
      include 'mpif.h'
#endif

    type(fce), intent(inout) :: this
    _REAL_, intent(out) :: force(3,this%natom)
    _REAL_, intent(in) :: coord(3,this%natom)

    _REAL_, intent(out) :: sforce(3,this%natom)
    _REAL_, intent(out) :: scoord(3,this%natom)

    integer :: iatm,jatm,id,err

    !LAPACK variables
    integer :: N, NN, NRHS=1, LDA, LDB, LWORK, INFO
    integer :: IPIV(this%nbase+1)
   _REAL_, pointer :: WORK(:)=>NULL()

    integer :: lswork, infos
   _REAL_, pointer :: swork(:)=>NULL()

    !uncertainties of input coordinates in A and B are assumed
    !to be negligible small
    _REAL_ :: RCOND=0.d0
    _REAL_ :: A(this%nbase+1,this%nbase+1),B(this%nbase+1,1)
    _REAL_ :: rrr(this%nbasis),rar(this%nbase)

    integer :: irr(this%nbasis),ird(this%nbase),irs(this%nbase)

    integer :: isa,jsa,iisa,jjsa

     DOUBLE PRECISION hmu(4,4),ddd,q0,q1,q2,q3,fis(3,3),wjjj(4)
     DOUBLE PRECISION xiii1,xiii2,xiii3,fiii1,fiii2,fiii3,rwrwr2

     N = this%nsample

! The transformation of the curent coordinate with respect to the first
! basic point is performed by minimization of the distance between them

       hmu=0.d0

! Summing up the quaternion matrix over atoms across the processes

       do iatm = this%atom0, this%atomF
       fiii1=this%coord(1,iatm,1)
       fiii2=this%coord(2,iatm,1)
       fiii3=this%coord(3,iatm,1)
       xiii1=coord(1,iatm)
       xiii2=coord(2,iatm)
       xiii3=coord(3,iatm)
       ddd=(xiii1+fiii1)**2+(xiii2+fiii2)**2+(xiii3+fiii3)**2
       hmu(1,1)=hmu(1,1)+((xiii1-fiii1)**2+(xiii2-fiii2)**2+ &
                          (xiii3-fiii3)**2)
       hmu(1,2)=hmu(1,2)+2.d0*(xiii2*fiii3-xiii3*fiii2)
       hmu(1,3)=hmu(1,3)+2.d0*(xiii3*fiii1-xiii1*fiii3)
       hmu(1,4)=hmu(1,4)+2.d0*(xiii1*fiii2-xiii2*fiii1)
       hmu(2,2)=hmu(2,2)+(ddd-2.d0*(xiii1*fiii1+fiii1*xiii1))
       hmu(2,3)=hmu(2,3)     -2.d0*(xiii1*fiii2+fiii1*xiii2)
       hmu(2,4)=hmu(2,4)     -2.d0*(xiii1*fiii3+fiii1*xiii3)
       hmu(3,3)=hmu(3,3)+(ddd-2.d0*(xiii2*fiii2+fiii2*xiii2))
       hmu(3,4)=hmu(3,4)     -2.d0*(xiii2*fiii3+fiii2*xiii3)
       hmu(4,4)=hmu(4,4)+(ddd-2.d0*(xiii3*fiii3+fiii3*xiii3))
       end do

#ifdef MPI

! Collecting the partial sums from different processes and placing the
! result into the same quaternion matrix known for all processors

    call mpi_allreduce(MPI_IN_PLACE,hmu,16,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)

#endif

! The problem is reduced to find eigenvalues and eigenvectors
! of the corresponding quaternion symmetric matrix

! First query the optimal workspace

       lswork = -1
       swork => safemem_realloc(swork,1,.false.)
       CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )
       lswork = INT(swork(1))
       swork => safemem_realloc(swork,lswork,.false.)

! Now solve the eigenproblem

       CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )

       if(infos.ne.0) then
       print *,infos
       stop 'Message from DSYEV'
       end if

! Define the 4-dimensional quaternion corresponding to the smallest eigenvalue

       q0=hmu(1,1)
       q1=hmu(2,1)
       q2=hmu(3,1)
       q3=hmu(4,1)

! Determine the rotational matrix in 3-dimensional space

       fis(1,1)=q0**2+q1**2-q2**2-q3**2
       fis(1,2)=2.d0*(-q0*q3+q1*q2)
       fis(1,3)=2.d0*( q0*q2+q1*q3)
       fis(2,1)=2.d0*( q0*q3+q1*q2)
       fis(2,2)=q0**2+q2**2-q1**2-q3**2
       fis(2,3)=2.d0*(-q0*q1+q2*q3)
       fis(3,1)=2.d0*(-q0*q2+q1*q3)
       fis(3,2)=2.d0*( q0*q1+q2*q3)
       fis(3,3)=q0**2+q3**2-q1**2-q2**2

! Perform the transformations

       do iatm = 1, this%natom
   scoord(1,iatm)=fis(1,1)*coord(1,iatm)+fis(2,1)*coord(2,iatm)+fis(3,1)*coord(3,iatm)
   scoord(2,iatm)=fis(1,2)*coord(1,iatm)+fis(2,2)*coord(2,iatm)+fis(3,2)*coord(3,iatm)
   scoord(3,iatm)=fis(1,3)*coord(1,iatm)+fis(2,3)*coord(2,iatm)+fis(3,3)*coord(3,iatm)
       end do

! No selecting and sorting if nbasis=nbase and sort=0

          if(this%nbasis.eq.this%nbase.and.this%sort.eq.0) then

          do isa=1,this%nbase
          ird(isa)=isa
          end do

          end if

! Calculate the distances between the current solute position and its
! counterparts from the extended basis set if nbasis > nbase or sort/=0 

          if(this%nbasis.gt.this%nbase.or.this%sort.ne.0) then

          do isa=1,N
          rrr(isa)=0.d0
          do iatm = this%atom0, this%atomF
          rrr(isa)=rrr(isa)+(scoord(1,iatm)-this%scoord(1,iatm,isa))**2+ &
                            (scoord(2,iatm)-this%scoord(2,iatm,isa))**2+ &
                            (scoord(3,iatm)-this%scoord(3,iatm,isa))**2
          end do
          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,rrr,N,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#endif

! Normalizing the distances

          do isa=1,N
          rrr(isa)=dsqrt(rrr(isa)/this%natom)
          end do

          end if

! Select the first nbase minimal values among all N=nsample elements
! if N > nbase (note that nbase <= nsample <= nbasis). In other words,
! choose the best transformed force-coordinates pairs.

          if(this%nbasis.gt.this%nbase.and.this%sort.ne.0) then

          do isa=1,N
          irr(isa)=isa
          end do

          if(N.gt.this%nbase) then

          call selects(this%nbase,N,rrr,irr)

          if(this%nbase.ne.1) then

          do isa=1,this%nbase-1
          rar(isa)=rrr(isa)
          end do

! Now perform sorting of the selected minimal values in the descending order

          call hsort(this%nbase-1,rar,irs)

          end if

          irs(this%nbase)=this%nbase

          do isa=1,this%nbase
          ird(isa)=irr(irs(isa))
          end do

          end if

! On the beginning when N=nsample=nbase perform only sorting because
! then selecting is not necessary

          if(N.eq.this%nbase) then

          call hsort(N,rrr,ird)

          end if

          end if

! Only select the first nbase minimal values among all N=nsample elements
! if N > nbase without sorting (sort=0)

          if(this%nbasis.gt.this%nbase.and.this%sort.eq.0) then

          do isa=1,N
          irr(isa)=isa
          end do

          if(N.gt.this%nbase) call selects(this%nbase,N,rrr,irr)

          do isa=1,this%nbase
          ird(isa)=irr(isa)
          end do

          end if

! If sort/=0, perform sorting even for nbasis=nbase, i.e. when selecting
! is not necessary

          if(this%nbasis.eq.this%nbase.and.this%sort.ne.0) then

          call hsort(this%nbase,rrr,ird)

          end if

! Zero-order force approximation

          do iatm = this%atom0, this%atomF
          do id=1,3
          sforce(id,iatm) = 0.d0
          end do
          end do

          NN=this%nbase+1

          if(NN.gt.1) then

! High-order force extrapolation

          LDA = NN
          LDB = NN
          NRHS = 1

! Solve the linear least-square problem with or without weighting

! Using multiprocessor technique to fill in the symmetric matrix elements

          if(this%weigh.eq.0) then

          do isa=1,NN-1
          iisa=ird(isa)
          do jsa=isa,NN-1
          jjsa=ird(jsa)
          A(isa,jsa)=0.d0

       do iatm = this%atom0, this%atomF
       A(isa,jsa)=A(isa,jsa)+this%scoord(1,iatm,iisa)*this%scoord(1,iatm,jjsa) &
                            +this%scoord(2,iatm,iisa)*this%scoord(2,iatm,jjsa) &
                            +this%scoord(3,iatm,iisa)*this%scoord(3,iatm,jjsa)
          end do
       A(isa,jsa)=A(isa,jsa)/this%natom
          end do
          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,A,NN*NN,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#endif

          do isa=1,NN-1
          iisa=ird(isa)
          B(isa,1)=0.d0
          do iatm = this%atom0, this%atomF
          B(isa,1)=B(isa,1)+scoord(1,iatm)*this%scoord(1,iatm,iisa) &
                           +scoord(2,iatm)*this%scoord(2,iatm,iisa) &
                           +scoord(3,iatm)*this%scoord(3,iatm,iisa)
          end do
          B(isa,1)=B(isa,1)/this%natom
          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,B,NN,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#endif

          else

          A=0.d0

          do isa=1,NN-1
          iisa=ird(isa)
          do jsa=isa,NN-1
          jjsa=ird(jsa)

       do iatm = this%atom0, this%atomF
       do jatm = iatm+1, this%natom

          rwrwr2=sum((scoord(1:3,iatm)-scoord(1:3,jatm))**2)

          if(rwrwr2 < this%cut) then

          A(isa,jsa)=A(isa,jsa)+( &
                      (this%scoord(1,iatm,iisa)-this%scoord(1,jatm,iisa))* &
                      (this%scoord(1,iatm,jjsa)-this%scoord(1,jatm,jjsa))+ &
                      (this%scoord(2,iatm,iisa)-this%scoord(2,jatm,iisa))* &
                      (this%scoord(2,iatm,jjsa)-this%scoord(2,jatm,jjsa))+ &
                      (this%scoord(3,iatm,iisa)-this%scoord(3,jatm,iisa))* &
                      (this%scoord(3,iatm,jjsa)-this%scoord(3,jatm,jjsa)))/rwrwr2

          end if

          end do
          end do

          A(isa,jsa)=A(isa,jsa)/((this%natom*(this%natom-1))/2)

          end do
          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,A,NN*NN,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#endif

          B=0.d0

          do isa=1,NN-1
          iisa=ird(isa)
       do iatm = this%atom0, this%atomF
       do jatm = iatm+1, this%natom

          rwrwr2=sum((scoord(1:3,iatm)-scoord(1:3,jatm))**2)

          if(rwrwr2 < this%cut) then

          B(isa,1)=B(isa,1)+( &
                    (scoord(1,iatm)-scoord(1,jatm))* &
                    (this%scoord(1,iatm,iisa)-this%scoord(1,jatm,iisa))+ &
                    (scoord(2,iatm)-scoord(2,jatm))* &
                    (this%scoord(2,iatm,iisa)-this%scoord(2,jatm,iisa))+ &
                    (scoord(3,iatm)-scoord(3,jatm))* &
                    (this%scoord(3,iatm,iisa)-this%scoord(3,jatm,iisa)))/rwrwr2

          end if

          end do
          end do

          B(isa,1)=B(isa,1)/((this%natom*(this%natom-1))/2)

          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,B,NN,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#endif

          end if

! Complemented minimization of the norms of the solutions with enormsw-weight

          do isa=1,NN-1
          A(isa,isa)=A(isa,isa)+this%enormsw
          end do

! Using Lagrange method to normalize the solutions

          do isa=1,NN-1
          A(isa,NN)=-1.d0
          end do
          A(NN,NN)=0.d0

          B(NN,1)=-1.d0

! Solving the extended system of linear equations

! First query the optimal workspace

          LWORK = -1
          WORK => safemem_realloc(WORK,1,.false.)
          call DSYSV( 'U', NN, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )
          LWORK = INT(WORK(1))
          WORK => safemem_realloc(WORK,LWORK,.false.)

! Perform the actual calculations

          call DSYSV( 'U', NN, NRHS, A, LDA, IPIV, B, LDB, WORK, LWORK, INFO )

          if(INFO.ne.0) then
          print *,INFO
          stop 'Message from DSYSV'
          end if

! Extrapolate the forces in the transformed space

          do iatm = this%atom0, this%atomF
          do isa=1,NN-1
          iisa=ird(isa)
          do id=1,3
          sforce(id,iatm) = sforce(id,iatm)+B(isa,1)*this%sforce(id,iatm,iisa)
          end do
          end do
          end do

          end if

! Performing the inverse transformation to obtain the extrapolated forces
! in the usual coordinates

          do iatm = this%atom0, this%atomF
    force(1,iatm)=fis(1,1)*sforce(1,iatm)+fis(1,2)*sforce(2,iatm)+fis(1,3)*sforce(3,iatm)
    force(2,iatm)=fis(2,1)*sforce(1,iatm)+fis(2,2)*sforce(2,iatm)+fis(2,3)*sforce(3,iatm)
    force(3,iatm)=fis(3,1)*sforce(1,iatm)+fis(3,2)*sforce(2,iatm)+fis(3,3)*sforce(3,iatm)
          end do

    err = safemem_dealloc(swork)
    err = safemem_dealloc(WORK)

  end subroutine fce_forcebm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! 3-transformation and selecting with respect to the current coordinate
     !   (most precise but time-expensive, not recommended for large nbasis)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fce_forcec(this,force,coord,sforce,scoord)

    use safemem
    implicit none

#ifdef MPI
      include 'mpif.h'
#endif

    type(fce), intent(inout) :: this
    _REAL_, intent(out) :: force(3,this%natom)
    _REAL_, intent(in) :: coord(3,this%natom)

    _REAL_, intent(out) :: sforce(3,this%natom)
    _REAL_, intent(out) :: scoord(3,this%natom)

    integer :: iatm,jatm,id,err

    !LAPACK variables
    integer :: M, N, NN, NRHS=1, LDA, LDB, LWORK, RANK, INFO
    integer, pointer :: JPVT(:)=>NULL()

    integer :: lswork, infos
   _REAL_, pointer :: swork(:)=>NULL()

    !uncertainties of input coordinates in A and B are assumed
    !to be negligible small
    _REAL_ :: RCOND=0.d0
    _REAL_, pointer :: A(:,:)=>NULL(), B(:,:)=>NULL(), work(:)=>NULL()
    _REAL_ :: rrr(this%nbasis),rar(this%nbase)

    integer :: irr(this%nbasis),ird(this%nbase),irs(this%nbase)

    integer :: isa,jsa,is1,iisap1

     DOUBLE PRECISION hmu(4,4),ddd,q0,q1,q2,q3,fis(3,3),wjjj(4),rweight
     DOUBLE PRECISION xiii1,xiii2,xiii3,fiii1,fiii2,fiii3,rwrwr2

#ifdef MPI
     DOUBLE PRECISION hmus(this%nbasis,4,4)
     integer inum,jsad(this%mpisize)
#endif

     N = this%nsample

! Transform each basic force-coordinate pairs with respect to the current
! point by minimization of the distance between them and that point

#ifdef MPI

! Summing up the quaternion matrix over atoms across the processes

       hmus=0.d0
       do iatm = this%atom0, this%atomF
       fiii1=coord(1,iatm)
       fiii2=coord(2,iatm)
       fiii3=coord(3,iatm)
       do isa=1,N
       xiii1=this%coord(1,iatm,isa)
       xiii2=this%coord(2,iatm,isa)
       xiii3=this%coord(3,iatm,isa)
       ddd=(xiii1+fiii1)**2+(xiii2+fiii2)**2+(xiii3+fiii3)**2
       hmus(isa,1,1)=hmus(isa,1,1)+((xiii1-fiii1)**2+(xiii2-fiii2)**2+ &
                          (xiii3-fiii3)**2)
       hmus(isa,1,2)=hmus(isa,1,2)+2.d0*(xiii2*fiii3-xiii3*fiii2)
       hmus(isa,1,3)=hmus(isa,1,3)+2.d0*(xiii3*fiii1-xiii1*fiii3)
       hmus(isa,1,4)=hmus(isa,1,4)+2.d0*(xiii1*fiii2-xiii2*fiii1)
       hmus(isa,2,2)=hmus(isa,2,2)+(ddd-2.d0*(xiii1*fiii1+fiii1*xiii1))
       hmus(isa,2,3)=hmus(isa,2,3)     -2.d0*(xiii1*fiii2+fiii1*xiii2)
       hmus(isa,2,4)=hmus(isa,2,4)     -2.d0*(xiii1*fiii3+fiii1*xiii3)
       hmus(isa,3,3)=hmus(isa,3,3)+(ddd-2.d0*(xiii2*fiii2+fiii2*xiii2))
       hmus(isa,3,4)=hmus(isa,3,4)     -2.d0*(xiii2*fiii3+fiii2*xiii3)
       hmus(isa,4,4)=hmus(isa,4,4)+(ddd-2.d0*(xiii3*fiii3+fiii3*xiii3))
       end do
       end do

! Collecting the partial sums from different processes and placing the
! result into the same quaternion matrix known for all processors

    call mpi_allreduce(MPI_IN_PLACE,hmus,this%nbasis*16,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)

       do isa=1,N

       hmu(:,:)=hmus(isa,:,:)

#else

! Summing up the quaternion matrix over atoms in the single process mode

       do isa=1,N
       hmu=0.d0
       do iatm = 1, this%natom
       fiii1=coord(1,iatm)
       fiii2=coord(2,iatm)
       fiii3=coord(3,iatm)
       xiii1=this%coord(1,iatm,isa)
       xiii2=this%coord(2,iatm,isa)
       xiii3=this%coord(3,iatm,isa)
       ddd=(xiii1+fiii1)**2+(xiii2+fiii2)**2+(xiii3+fiii3)**2
       hmu(1,1)=hmu(1,1)+((xiii1-fiii1)**2+(xiii2-fiii2)**2+ &
                          (xiii3-fiii3)**2)
       hmu(1,2)=hmu(1,2)+2.d0*(xiii2*fiii3-xiii3*fiii2)
       hmu(1,3)=hmu(1,3)+2.d0*(xiii3*fiii1-xiii1*fiii3)
       hmu(1,4)=hmu(1,4)+2.d0*(xiii1*fiii2-xiii2*fiii1)
       hmu(2,2)=hmu(2,2)+(ddd-2.d0*(xiii1*fiii1+fiii1*xiii1))
       hmu(2,3)=hmu(2,3)     -2.d0*(xiii1*fiii2+fiii1*xiii2)
       hmu(2,4)=hmu(2,4)     -2.d0*(xiii1*fiii3+fiii1*xiii3)
       hmu(3,3)=hmu(3,3)+(ddd-2.d0*(xiii2*fiii2+fiii2*xiii2))
       hmu(3,4)=hmu(3,4)     -2.d0*(xiii2*fiii3+fiii2*xiii3)
       hmu(4,4)=hmu(4,4)+(ddd-2.d0*(xiii3*fiii3+fiii3*xiii3))
       end do

#endif

! The problem is reduced to find eigenvalues and eigenvectors
! of the corresponding quaternion symmetric matrices

! First query the optimal workspace

       lswork = -1
       swork => safemem_realloc(swork,1,.false.)
       CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )
       lswork = INT(swork(1))
       swork => safemem_realloc(swork,lswork,.false.)

! Now solve the eigenproblem

       CALL  DSYEV( 'V', 'U', 4, hmu, 4, wjjj, swork, lswork, infos )

       if(infos.ne.0) then
       print *,infos
       stop 'Message from DSYEV'
       end if

! Define the 4-dimensional quaternion corresponding to the smallest eigenvalue

       q0=hmu(1,1)
       q1=hmu(2,1)
       q2=hmu(3,1)
       q3=hmu(4,1)

! Determine the rotational matrix in 3-dimensional space

       fis(1,1)=q0**2+q1**2-q2**2-q3**2
       fis(1,2)=2.d0*(-q0*q3+q1*q2)
       fis(1,3)=2.d0*( q0*q2+q1*q3)
       fis(2,1)=2.d0*( q0*q3+q1*q2)
       fis(2,2)=q0**2+q2**2-q1**2-q3**2
       fis(2,3)=2.d0*(-q0*q1+q2*q3)
       fis(3,1)=2.d0*(-q0*q2+q1*q3)
       fis(3,2)=2.d0*( q0*q1+q2*q3)
       fis(3,3)=q0**2+q3**2-q1**2-q2**2

! Perform the coordinate and force transformations

       do iatm = 1, this%natom
    this%scoord(1,iatm,isa)=fis(1,1)*this%coord(1,iatm,isa)+fis(2,1)*this%coord(2,iatm,isa)+fis(3,1)*this%coord(3,iatm,isa)
    this%scoord(2,iatm,isa)=fis(1,2)*this%coord(1,iatm,isa)+fis(2,2)*this%coord(2,iatm,isa)+fis(3,2)*this%coord(3,iatm,isa)
    this%scoord(3,iatm,isa)=fis(1,3)*this%coord(1,iatm,isa)+fis(2,3)*this%coord(2,iatm,isa)+fis(3,3)*this%coord(3,iatm,isa)
       end do

       do iatm = this%atom0, this%atomF
    this%sforce(1,iatm,isa)=fis(1,1)*this%force(1,iatm,isa)+fis(2,1)*this%force(2,iatm,isa)+fis(3,1)*this%force(3,iatm,isa)
    this%sforce(2,iatm,isa)=fis(1,2)*this%force(1,iatm,isa)+fis(2,2)*this%force(2,iatm,isa)+fis(3,2)*this%force(3,iatm,isa)
    this%sforce(3,iatm,isa)=fis(1,3)*this%force(1,iatm,isa)+fis(2,3)*this%force(2,iatm,isa)+fis(3,3)*this%force(3,iatm,isa)
       end do

       if(this%nbasis.gt.this%nbase.or.this%sort.ne.0) rrr(isa)=wjjj(1)

       end do

       do iatm = 1, this%natom
       scoord(1,iatm)=coord(1,iatm)
       scoord(2,iatm)=coord(2,iatm)
       scoord(3,iatm)=coord(3,iatm)
       end do

! No selecting and sorting if nbasis=nbase and sort=0

          if(this%nbasis.eq.this%nbase.and.this%sort.eq.0) then

          do isa=1,this%nbase
          ird(isa)=isa
          end do

          end if

! Select the first nbase minimal values among all N=nsample elements
! if N > nbase (note that nbase <= nsample <= nbasis). In other words,
! choose the best transformed force-coordinates pairs.

          if(this%nbasis.gt.this%nbase.and.this%sort.ne.0) then

          do isa=1,N
          irr(isa)=isa
          end do

          if(N.gt.this%nbase) then

          call selects(this%nbase,N,rrr,irr)

          if(this%nbase.ne.1) then

          do isa=1,this%nbase-1
          rar(isa)=rrr(isa)
          end do

! Now perform sorting of the selected minimal values in the descending order

          call hsort(this%nbase-1,rar,irs)

          end if

          irs(this%nbase)=this%nbase

          do isa=1,this%nbase
          ird(isa)=irr(irs(isa))
          end do

          end if

! On the beginning when N=nsample=nbase perform only sorting because
! then selecting is not necessary

          if(N.eq.this%nbase) then

          call hsort(N,rrr,ird)

          end if

          end if

! Only select the first nbase minimal values among all N=nsample elements
! if N > nbase without sorting (sort=0)

          if(this%nbasis.gt.this%nbase.and.this%sort.eq.0) then

          do isa=1,N
          irr(isa)=isa
          end do

          if(N.gt.this%nbase) call selects(this%nbase,N,rrr,irr)

          do isa=1,this%nbase
          ird(isa)=irr(isa)
          end do

          end if

! If sort/=0, perform sorting even for nbasis=nbase, i.e. when selecting
! is not necessary

          if(this%nbasis.eq.this%nbase.and.this%sort.ne.0) then

          call hsort(this%nbase,rrr,ird)

          end if

! The nearest pair

          is1=ird(1)

! Zero-order force approximation

          do iatm = this%atom0, this%atomF
          do id=1,3
          sforce(id,iatm) = this%sforce(id,iatm,is1)
          end do
          end do

          NN=this%nbase-1

          if(NN.gt.0) then

! High-order force extrapolation

          JPVT=>safemem_realloc(JPVT,NN,.false.)
          JPVT=0

! Solve the linear least-square problem with or without weighting

          if(this%weigh.eq.0) then

          M = 3 * this%natom
          LDA = M
          LDB = max(M,NN)

          A=>safemem_realloc(A,M,NN,.false.)
          B=>safemem_realloc(B,max(NN,M),1,.false.)

! Using multiprocessor technique to fill in the matrix elements

#ifdef MPI
          A=0.d0
          B=0.d0

          jsa=3*(this%atom0-1)
#else
          jsa=0
#endif

          do iatm = this%atom0, this%atomF
          B(jsa+1:jsa+3,1)=scoord(:,iatm)-this%scoord(:,iatm,is1)
          do isa=1,NN
          iisap1=ird(isa+1)
          A(jsa+1:jsa+3,isa)=this%scoord(:,iatm,iisap1)-this%scoord(:,iatm,is1)
          end do
          jsa=jsa+3
          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,A,M*NN,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
    call mpi_allreduce(MPI_IN_PLACE,B,LDB,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#endif

          !first determine the amount of memory required
          LWORK=-1
          WORK => safemem_realloc(work,1,.false.)
          call dgelsy(M,NN,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,LWORK,INFO)

          !perform the calculation
          LWORK=WORK(1)
          WORK => safemem_realloc(WORK,LWORK,.false.)
          call dgelsy(M,NN,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,LWORK,INFO)

          if(INFO.ne.0) then
          print *,INFO
          stop 'Message from DGELSY'
          end if

          else

          jsa=0

#ifdef MPI

          jsad=0
          do iatm = this%atom0, this%atomF
          do jatm = iatm+1, this%natom
          rwrwr2=sum((scoord(1:3,iatm)-scoord(1:3,jatm))**2)
          if(rwrwr2 < this%cut) then
          jsa=jsa+3
          end if
          end do
          end do

    call MPI_Gather(jsa,1,MPI_INTEGER,jsad,1,MPI_INTEGER,0,this%mpicomm,err);

      if(this%mpirank == 0) then
      do inum=2,this%mpisize
      jsad(inum)=jsad(inum-1)+jsad(inum)
      end do
      end if

    call mpi_allreduce(MPI_IN_PLACE,jsad,this%mpisize,MPI_INTEGER,MPI_SUM,this%mpicomm,err)

      if(this%mpirank == 0) jsa=0
      do inum=1,this%mpisize-1
      if(this%mpirank == inum) jsa=jsad(inum)
      end do

          M=jsad(this%mpisize)
#else
          M = 3 * ((this%natom*(this%natom-1))/2)
#endif

          LDA = M
          LDB = max(M,NN)

          A=>safemem_realloc(A,M,NN,.false.)
          B=>safemem_realloc(B,max(NN,M),1,.false.)

          A=0.d0
          B=0.d0

          do iatm = this%atom0, this%atomF
          do jatm = iatm+1, this%natom

          rwrwr2=sum((scoord(1:3,iatm)-scoord(1:3,jatm))**2)

! Use cutting for weighting

          if(rwrwr2 < this%cut) then

          rweight=1.d0/dsqrt(rwrwr2)

          B(jsa+1:jsa+3,1) = rweight*(&
                             scoord(:,iatm)-this%scoord(:,iatm,is1) &
                            -scoord(:,jatm)+this%scoord(:,jatm,is1))
          do isa=1,NN
          iisap1=ird(isa+1)
          A(jsa+1:jsa+3,isa) = rweight*(&
                          this%scoord(:,iatm,iisap1)-this%scoord(:,iatm,is1) &
                         -this%scoord(:,jatm,iisap1)+this%scoord(:,jatm,is1))
          end do
          jsa=jsa+3
          end if
          end do
          end do

#ifdef MPI
    call mpi_allreduce(MPI_IN_PLACE,A,M*NN,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
    call mpi_allreduce(MPI_IN_PLACE,B,LDB,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#else
          M=jsa
#endif

          !first determine the amount of memory required
          LWORK=-1
          WORK => safemem_realloc(work,1,.false.)
          call dgelsy(M,NN,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,LWORK,INFO)

          !perform the calculation
          LWORK=WORK(1)
          WORK => safemem_realloc(WORK,LWORK,.false.)
          call dgelsy(M,NN,NRHS,A,LDA,B,LDB,JPVT,RCOND,RANK,WORK,LWORK,INFO)

          if(INFO.ne.0) then
          print *,INFO
          stop 'Message from DGELSY'
          end if

          end if

! Extrapolate the forces in the transformed space

          do iatm = this%atom0, this%atomF
          do isa=1,NN
          iisap1=ird(isa+1)
          do id=1,3
          sforce(id,iatm) = sforce(id,iatm)+&
          B(isa,1)*(this%sforce(id,iatm,iisap1)-this%sforce(id,iatm,is1))
          end do
          end do
          end do

          end if

          do iatm = this%atom0, this%atomF
          force(1,iatm)=sforce(1,iatm)
          force(2,iatm)=sforce(2,iatm)
          force(3,iatm)=sforce(3,iatm)
          end do

    err = safemem_dealloc(A)
    err = safemem_dealloc(B)
    err = safemem_dealloc(WORK)
    err = safemem_dealloc(JPVT)
    err = safemem_dealloc(swork)

  end subroutine fce_forcec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! 4-no normalization, coordinate transformation, weighting, selecting,
     !   and sorting, but with individual force extrapolation and possible
     !   cutting-off the neighbours [original AMBER11 version]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fce_force(this,force,coord)
    use safemem
    implicit none
    type(fce), intent(inout) :: this
    _REAL_, intent(out) :: force(3,this%natom)
    _REAL_, intent(in) :: coord(3,this%natom)

    integer :: iatm,jatm,id,err
    
    !LAPACK variables
    integer :: M,N,NRHS, LDA,LDB, LWORK,RANK=0, INFO=0
    integer, pointer :: JPVT(:)=>NULL()
    !not really sure what value this should be...
    _REAL_ :: RCOND=0d0
    _REAL_, pointer :: A(:,:)=>NULL(),B(:,:)=>NULL(),work(:)=>NULL()

    _REAL_,external :: ddot

!!$    write(6,*) "FCE_FORCE", this%atom0, this%atomF, this%cut

    force=0
    N = this%nbasis
    JPVT=>safemem_realloc(JPVT,N,.false.)
    JPVT=0
    !the cutoff and non-cutoff situations are handled differently
    if(this%cut>0) then
       !update cutoff list if requested
       call nlist(this)

       !iterate over each atom and calculate the extrapolated force
       do iatm = this%atom0, this%atomF
!          call orient(this)

          M = 3*(this%cutlist(1,iatm)+1)
          LDA = M
          LDB = max(M,N)
          NRHS = 1

          A=>safemem_realloc(A,M,N,.false.)
          B=>safemem_realloc(B,max(N,M),1,.false.)

          !the target atom is never in the cutlist so put it in the first index
          B(1:3,1) = coord(:,iatm)
          A(1:3,:) = this%coord(:,iatm,:)

          !now add the atoms from within the cutoff
          do jatm=2,this%cutlist(1,iatm)+1
             B(jatm*3-2: jatm*3,1) = coord(:,this%cutlist(jatm,iatm))
             A(jatm*3-2: jatm*3,:) = this%coord(:,this%cutlist(jatm,iatm),:)
          end do

          !first determine the amount of memory required
          LWORK=-1
          WORK => safemem_realloc(work,1,.false.)
          call dgelsy(M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, INFO )
          
          !perform the calculation
          LWORK=WORK(1)
          WORK => safemem_realloc(WORK,LWORK,.false.)
          call dgelsy(M, N, NRHS, A, LDA, B, LDB, JPVT, RCOND, RANK, WORK, LWORK, INFO )

          !the force in each dimension is now a dot product
          do id=1,3
             force(id,iatm) = ddot(N,this%force(id,iatm,:),1,B(1:N,1),1)
          end do
!!$          write(6,*) "FCE FORCE",iatm, force(:,iatm)
 !         call unorient(this)
       end do
    else !cutoff
       
    end if !cutoff
    err = safemem_dealloc(A)
    err = safemem_dealloc(B)
    err = safemem_dealloc(WORK)
    err = safemem_dealloc(JPVT)
  end subroutine fce_force

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!private subroutines and functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Very simple method to generate a cutoff list.  The resulting matrix is 2D is
!!!lists the number of neighbour atoms for each atom (first row) and then the
!!!atom numbers themselves.
!!!E.g.
!!!   1 2 1
!!!   2 1 2
!!!   0 3 0
!!!   0 0 0
!!!For this three atom system, atoms 1 and 2, 3 and 2 are neighbours.  I.e.
!!!atoms 1 and 3 each have one neighbour and 2 has two.  The number of
!!!neighbours is list in the first row.  The subsequent rows list the ids of the
!!!neighbours.
!!!IN:
!!!   this : FCE object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine nlist(this)
  implicit none
  type(fce),intent(inout) :: this
  integer :: id,iatm1,iatm2
#ifdef RISM_DEBUG
  write(6,*) "NLIST"; call flush(6)
#endif /*RISM_DEBUG*/

  this%cutlist=0
  do iatm1 = 1,this%natom
     do iatm2 = iatm1+1,this%natom
        if( sum((this%coord(1:3,iatm1,1) - this%coord(1:3,iatm2,1))**2) < this%cut)then
           this%cutlist(1,iatm1) = this%cutlist(1,iatm1) +1
           this%cutlist(1,iatm2) = this%cutlist(1,iatm2) +1
           this%cutlist(this%cutlist(1,iatm1)+1,iatm1) = iatm2
           this%cutlist(this%cutlist(1,iatm2)+1,iatm2) = iatm1
        end if
     end do
  end do
end subroutine nlist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      Estimate the accuracy of the extrapolation by measuring the     !
!      difference between the exact and extrapolated forces at the     !
!      end of the outer time interval                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fce_estimate(this,force,coord,sforce,scoord,forcem,deviat)
    implicit none

#ifdef MPI
      include 'mpif.h'
#endif

    type(fce), intent(inout) :: this
    _REAL_, intent(in) :: force(3,this%natom), coord(3,this%natom)

    _REAL_, intent(out) :: sforce(3,this%natom), scoord(3,this%natom)

    _REAL_ forcem(3,this%natom),stforce(3,this%natom),sas,sasa

    DOUBLE PRECISION deviat

     integer :: iaa

#ifdef MPI
    _REAL_ ssas,ssasa
     integer :: err
#endif

    save sas,sasa

        if(this%nsample.eq.0) then
        sas=0.d0
        sasa=0.d0
        end if

! Calculating the extrapolated 3D-RISM forces

       if(this%nsample.ge.this%nbase) then

       if(this%trans.eq.0) call fce_forcea(this,forcem,coord)

       if(this%trans.eq.1) call fce_forceb(this,forcem,coord,sforce,scoord)

       if(this%trans.eq.2) call fce_forcebm(this,forcem,coord,sforce,scoord)

       if(this%trans.eq.3) call fce_forcec(this,forcem,coord,sforce,scoord)

       if(this%trans.eq.4) call fce_force(this,forcem,coord)

! Summing up the exact 3D-RISM forces from different processes and
! making the result to be known for all processors

#ifdef MPI
    call mpi_allreduce(force,stforce,3*this%natom,MPI_DOUBLE_PRECISION,MPI_SUM,this%mpicomm,err)
#else
    stforce=force
#endif

! Calculating the relative force deviations using multiprocessing

       do iaa = this%atom0, this%atomF

          sas=sas+(forcem(1,iaa)-stforce(1,iaa))**2 &
                 +(forcem(2,iaa)-stforce(2,iaa))**2 &
                 +(forcem(3,iaa)-stforce(3,iaa))**2

          sasa=sasa+stforce(1,iaa)**2+stforce(2,iaa)**2+stforce(3,iaa)**2

       end do

#ifdef MPI

! Summing up the deviations and collecting them within the root processes

    call mpi_reduce(sas,ssas,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,this%mpicomm,err)
    call mpi_reduce(sasa,ssasa,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,this%mpicomm,err)

       if(this%mpirank==0) then

       deviat=dsqrt(ssas/ssasa)

       end if
#else
       deviat=dsqrt(sas/sasa)
#endif

       end if

end subroutine fce_estimate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module fce_c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!************************************************************************
! Modified Quick Select Algorithm
!************************************************************************

! The input array arr will be rearranged to have the kth smallest value
! in location arr(k), with all smaller elements moved to arr(1:k-1)
! (in arbitrary order) and all larger elements in arr(k+1:) 

      subroutine selects(k,n,arr,irr)
      INTEGER k,n
      double precision arr(*)
      INTEGER i,ir,j,l,mid,irr(*)
      double precision a,temp
      INTEGER ia,itemp
      l=1
      ir=n
1     if(ir-l.le.1)then
        if(ir-l.eq.1)then
          if(arr(ir).lt.arr(l))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp

            itemp=irr(l)
            irr(l)=irr(ir)
            irr(ir)=itemp

          endif
        endif
        return
      else
        mid=(l+ir)/2
        temp=arr(mid)
        arr(mid)=arr(l+1)
        arr(l+1)=temp

        itemp=irr(mid)
        irr(mid)=irr(l+1)
        irr(l+1)=itemp

        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp

         itemp=irr(l+1)
         irr(l+1)=irr(ir)
         irr(ir)=itemp

        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp

            itemp=irr(l)
            irr(l)=irr(ir)
            irr(ir)=itemp

        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp

          itemp=irr(l+1)
          irr(l+1)=irr(l)
          irr(l)=itemp

        endif
        i=l+1
        j=ir
        a=arr(l)

        ia=irr(l)

3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp

        itemp=irr(i)
        irr(i)=irr(j)
        irr(j)=itemp

        goto 3
5       arr(l)=arr(j)

        irr(l)=irr(j)

        arr(j)=a

        irr(j)=ia

        if(j.ge.k)ir=j-1
        if(j.le.k)l=i
      endif
      goto 1
      END
!  (C) Copr. 1986-92 Numerical Recipes Software

!************************************************************************
! Super Efficient Sorting (Adv. Eng. Software, 1984, Vol.6, No. 4, p.198
!************************************************************************
      subroutine hsort(n,data,list)
      integer n
      integer list(*)
      double precision data(*)
      integer maxstk, ncut
      parameter (maxstk=32, ncut=12)
      integer ll, lr, lm, nl, nr, ltemp, ist, i, j, k
      integer lstack(maxstk), rstack(maxstk)
      double precision guess, value
      do i=1,n
      list(i)=i
      end do
      ll=1
      lr=n
      ist=0
   10 continue
      if((lr-ll).ge.ncut) then
      nl=ll
      nr=lr
      lm=(ll+lr)/2
      guess=data(list(lm))
   20 continue
      if(data(list(nl)).lt.guess) then
      nl=nl+1
      goto 20
      end if
   30 continue
      if(guess.lt.data(list(nr))) then
      nr=nr-1
      goto 30
      end if
      if(nl.lt.(nr-1)) then
      ltemp=list(nl)
      list(nl)=list(nr)
      list(nr)=ltemp
      nl=nl+1
      nr=nr-1
      goto 20
      end if
      if(nl.le.nr) then
      if(nl.lt.nr) then
      ltemp=list(nl)
      list(nl)=list(nr)
      list(nr)=ltemp
      end if
      nl=nl+1
      nr=nr-1
      end if
      ist=ist+1
      if(nr.lt.lm) then
      lstack(ist)=nl
      rstack(ist)=lr
      lr=nr
      else
      lstack(ist)=ll
      rstack(ist)=nr
      ll=nl
      end if
      goto 10
      end if
      if(ist.ne.0) then
      ll=lstack(ist)
      lr=rstack(ist)
      ist=ist-1
      goto 10
      end if
      j=1
      k=list(1)
      value=data(k)
      do 40 i=2,min(n,ncut)
      if(data(list(i)).lt.value) then
      j=i
      value=data(list(i))
      end if
   40 continue
      list(1)=list(j)
      list(j)=k
      do 60 i=2,n
      j=i
      k=list(i)
      value=data(k)
   50 continue
      if(value.lt.data(list(j-1))) then
      list(j)=list(j-1)
      j=j-1
      goto 50
      end if
      list(j)=k
   60 continue
      return
      end
