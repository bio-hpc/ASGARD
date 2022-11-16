!<compile=optimized>

!The 3D-RISM-KH software found here is copyright (c) 2010-2012 by Andriy Kovalenko,
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

#include "../include/dprec.fh"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Original MDIIS implementation by Andriy Kovalenko.
!!Parallelization by Sergey Gusarov.
!!Further modifications and objectification by Tyler Luchko.
!!New memory model and BLASification Tyler Luchko December, 2009
!!
!!The memory model has been modified such that we are always dealing
!!with contiguous blocks of memory.  As vectors are added to working
!!from, the memory (alredy allocated) is accessed in increasing column
!!numbers.  This helps with cache prediction and loading for BLAS
!!routines.  Also, there is no longer any need to map between vector
!!numbers and their memory locations.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mdiis_blas_c
  use rism_report_c
  use rism_timer_c
  implicit none
  type mdiis_blas
     private
     !istep  :: number of steps since start/restart
     !np     :: number of data points (first dimension length)
     !nVec0   :: current number of vectors
     !nVec    :: max number of vectors
     integer ::  istep = 0, nVec0 = 0, np, nVec
     !overlap :: overlap matrix
     _REAL_,dimension(:,:),pointer :: overlap=>NULL()
     !ratio of the current residual to the minimum residual found that causes a restart
     _REAL_ ::  restart
     integer :: mpirank=0, mpisize=1, mpicomm=0
     !delta :: coefficient for residual gradient.  Should be between 0 and 1.
     !tol   :: target residual
     _REAL_ :: delta, tol
     !pointers to working data
     !xi    :: array of vector data.  np X nVec
     !ri    :: array of residual data.  np X nVec
     _REAL_,pointer :: xi(:,:)=>NULL(), ri(:,:)=>NULL()

     type(rism_timer) :: overlapTimer, restartTimer, &
          lapackTimer, projectTimer
  end type mdiis_blas

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Create new MDIIS object.    Provide convergence parameters and
!!!working memory.  Note that the working memory must me nVec times
!!!larger than the data array.  Further, the number of data points
!!!cannot change between iterations.  xi(:,1) and ri(:,1) are the active vectors.
!!!IN:
!!!   this  :: mdiis object
!!!   delta :: coefficient for residual gradient.  Should be between 0 and 1.
!!!   tol   :: target residual
!!!   xi    :: array of vector data.  np X nVec
!!!   ri    :: array of residual data.  np X nVec
!!!   np    :: number of data points (first dimension length)
!!!   nVec   :: length of DIIS vectors
!!!   restart :: restart threshold factor. Ratio of the current residual to the 
!!!   timer :: parent rism_timer
!!!              minimum residual in the basis that causes a restart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_blas_new_serial (this,delta,tol,restart,timer)
    implicit none
    type(mdiis_blas),intent(inout) :: this
    _REAL_,intent(in) :: delta, tol, restart
    type(rism_timer), intent(inout) :: timer

    this%delta = delta
    this%tol = tol
    this%restart = restart
    this%istep=0
    this%nVec0=0
    call rism_timer_new(this%overlapTimer,"Overlap")
    call rism_timer_setParent(this%overlapTimer,timer)
    call rism_timer_new(this%restartTimer,"Restart")
    call rism_timer_setParent(this%restartTimer,timer)
    call rism_timer_new(this%lapackTimer,"Lapack")
    call rism_timer_setParent(this%lapackTimer,timer)
    call rism_timer_new(this%projectTimer,"Project")
    call rism_timer_setParent(this%projectTimer,timer)
  end subroutine mdiis_blas_new_serial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Create new MDIIS object.  Provide convergence parameters, working
!!!memory and MPI parameters.  Note that the working memory must me
!!!nVec times larger than the data array.  Further, the number of data
!!!points cannot change between iterations.  xi(:,1) and ri(:,1) are
!!!the active vectors.
!!!IN:
!!!   this  :: mdiis object
!!!   delta :: coefficient for residual gradient.  Should be between 0 and 1.
!!!   tol   :: target residual
!!!   xi    :: array of vector data.  np X nVec
!!!   ri    :: array of residual data.  np X nVec
!!!   np    :: number of data points (first dimension length)
!!!   nVec   :: length of DIIS vectors
!!!   restart :: restart threshold factor. Ratio of the current residual to the 
!!!              minimum residual in the basis that causes a restart
!!!   rank  :: MPI process rank
!!!   size  :: Number of processes
!!!   comm  :: MPI communicator
!!!   timer :: parent rism_timer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_blas_new_mpi(this,delta,tol,restart,timer,&
       rank,size,comm)
    implicit none
    type(mdiis_blas),intent(out) :: this
    _REAL_,intent(in) :: delta, tol, restart
    integer, intent(in) :: rank, size, comm
    type(rism_timer), intent(inout) :: timer
    call mdiis_blas_new_serial(this,delta,tol,restart,timer)
    this%mpirank = rank
    this%mpisize = size
    this%mpicomm = comm
  end subroutine mdiis_blas_new_mpi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Destroys the MDIIS object.  All variables are reset except for working memory.
!!!IN:
!!!   this :: mdiis object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_blas_destroy(this)
    use safemem
    implicit none
    type(mdiis_blas),intent(inout) :: this
    integer :: err
    call mdiis_blas_reset(this)
    if(safemem_dealloc(this%overlap)/=0) call rism_report_error("MDIIS_BLAS: failed to deallocate OVERLAP")
    this%delta=0
    this%tol=0
    this%np=0
    this%nVec0=0
    nullify(this%xi)
    nullify(this%ri)
    call rism_timer_destroy(this%overlapTimer)
    call rism_timer_destroy(this%restartTimer)
    call rism_timer_destroy(this%lapackTimer)
    call rism_timer_destroy(this%projectTimer)
  end subroutine mdiis_blas_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Reset the MDIIS object. Progress is reset to a single basis
!!!vector.  All other variables are untouched.
!!!IN:
!!!   this :: mdiis object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_blas_reset(this)
    use safemem
    implicit none
    type(mdiis_blas),intent(inout) :: this
    integer :: err
    this%overlap => safemem_realloc(this%overlap,this%nVec,this%nVec,.false.)
    this%overlap = 0
    this%istep=0
    this%nVec0=0
    this%overlap=0
  end subroutine mdiis_blas_reset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Resizes and resets MDIIS object. Progress is reset to a single basis
!!!vector. Working memory remains intact memory and is
!!!swapped such that xi(:,1) and ri(:,1) become the active vectors.
!!!IN:
!!!   this :: mdiis object
!!!   xi    :: array of vector data.  np X nvec
!!!   ri    :: array of residual data.  np X nvec
!!!   np    :: number of data points (first dimension length)
!!!   nvec   :: length of DIIS vectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_blas_resize(this,xi,ri,np,nvec)
    implicit none
    type(mdiis_blas),intent(inout) :: this
    _REAL_,target, intent(in) :: xi(np,nvec), ri(np,nvec)
    integer,intent(in) ::  np
    integer, intent(in) :: nvec
    !transfer the current working vector to the first index
    this%np = np
    this%nvec = nvec
    this%xi => xi
    this%ri => ri
    call mdiis_blas_reset(this)
  end subroutine mdiis_blas_resize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the number of vector currently used as a basis by the MDIIS object
!!!IN:
!!!   this :: mdiis object
!!!OUT:
!!!   The number of vectors currently used
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function mdiis_blas_getCurrentNVec (this) result(nVec)
    implicit none
    type(mdiis_blas),intent(in) :: this
    integer :: nVec
    nVec = this%nVec0
  end function mdiis_blas_getCurrentNVec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the index number of the working vector
!!!IN: 
!!!   this :: mdiis object
!!!OUT:
!!!    the number of the vector to read/write
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function mdiis_blas_getWorkVector (this) result(index)
    implicit none
    type(mdiis_blas),intent(in) :: this
    integer :: index
    index = 1
  end function mdiis_blas_getWorkVector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!One MDIIS iteration
!!!IN:
!!!   this  :: mdiis object
!!!   rms1  :: residual at the end of the step
!!!   conver:: have we converged
!!!   tolerance :: tolerance for this calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  mdiis_blas_advance (this, rms1,conver,tolerance_o)
    implicit none 
#ifdef MPI
    include 'mpif.h'
#endif /*MPI*/
#include "def_time.h" 
    
    type(mdiis_blas),intent(inout) :: this
    _REAL_,intent(out) ::  rms1
    logical,intent(out) ::  conver
    _REAL_,optional,intent(in) ::  tolerance_o

    _REAL_ :: tolerance

    !isiupd :: vector to be updated
    integer ::  isiupd
    !indx :: for LAPACK routines
    integer ::  indx(0:this%nVec)
    _REAL_ :: toverlap(1:this%nVec)

    !aij :: overlap matrix for LAPACK
    !bi  :: linear coefficients from LAPACK
    _REAL_ :: aij(0:this%nVec,0:this%nVec), bi(0:this%nVec,0:this%nVec)

    !temp :: swap variable
    _REAL_ :: temp
    !nVec1 :: number of vectors, not including the update vector
    integer :: nVec1

    !counter variables
    integer :: is,is1,is2
    integer :: err=0, ierr

#ifdef RISM_DEBUG
    write(6,*)"MDIIS start",this%delta,rms1,this%tol,this%np,this%nVec,conver
#endif
    tolerance = this%tol
    if(present(tolerance_o)) tolerance = tolerance_o
    !.................. initialize counters and switches ...................

    !.................. increment step and MDIIS counters ..................
    this%istep = this%istep + 1
    this%nVec0 = min( this%nVec0+1, this%nVec)

    !........................ specify storage list .........................
    isiupd = mod(this%istep-1,this%nVec-1)+2

    !............. calculate diagonal overlap of new residual ..............
    !                              and
    !........... calculate nondiagonal overlaps of new residual ............
    call rism_timer_start(this%overlapTimer)
    call timer_start(TIME_MDIIS_DATA)
#if defined(MPI)      
    call DGEMV ('T',this%np,this%nVec0,1d0,this%ri,this%np,&
         this%ri,1,0d0,toverlap,1)
    CALL MPI_AllREDUCE(toverlap(1:this%nVec0),this%overlap(1:this%nVec0,1),this%nVec0,&
         MPI_DOUBLE_PRECISION,MPI_SUM, this%mpicomm,ierr)
    if(err /=0) call rism_report_error&
         ("MDIIS_BLAS_UPDATE: could not reduce OVERLAP")
#else
    call DGEMV ('T',this%np,this%nVec0,1d0,this%ri,this%np,&
         this%ri,1,0d0,this%overlap,1)
#endif /*defined(MPI)*/
    call DCOPY(this%nVec0-1,this%overlap(2:this%nVec0,1),1,this%overlap(1,2),this%nVec)
    call timer_stop(TIME_MDIIS_DATA)
    call rism_timer_stop(this%overlapTimer)
 
    !................ get mean square value of new residual ................
    rms1 = sqrt( this%overlap(1,1)/(this%np*this%mpisize))
    !.............. check mean square residual for tolerance ...............
    if (rms1 <= tolerance)  then
       conver = .true.
       return
    else
       conver = .false.
    endif

    call rism_timer_start(this%restartTimer)
    call checkrestart(this,rms1,isiupd)
    call rism_timer_stop(this%restartTimer)
    nVec1 = min(this%nVec0,this%nVec-1)


    !......................... load DIIS matrices .........................
    call rism_timer_start(this%lapackTimer)
    call timer_start(TIME_MDIIS_LAPACK)
    aij(0,0) = 0d0
    bi=0
    bi(0,0) = -1d0
    aij(1:this%nVec0,0) = -1d0 
    aij(0,1:this%nVec0) = -1d0
    do is2=1,this%nVec0
       call DCOPY(this%nVec0,this%overlap(1:this%nVec0,is2),1,aij(1:this%nVec0,is2),1)
    enddo

    !....................... calculate DIIS estimate .......................
    call DGETRF(this%nVec0+1,this%nVec0+1,aij,this%nVec+1,indx,err)
    if(err > 0)then
       call rism_report_error("LU-factorization failed.  U = 0")
    elseif(err<0)then
       err = err*(-1)
       call rism_report_error("LU-factorization failed.")
    endif
    call DGETRS('N',this%nVec0+1,this%nVec0+1,aij,this%nVec+1,indx,bi,this%nVec+1,err)
    if(err < 0)then
       err = err*(-1)
       call rism_report_error("Linear equation solver failed.")
    endif
    call timer_stop(TIME_MDIIS_LAPACK)
    call rism_timer_stop(this%lapackTimer)

    !......... get DIIS minimum, MDIIS correction, and next point ..........
    call rism_timer_start(this%projectTimer)
    !swap bi for the index that will be discarded with the first index
    temp = bi(isiupd,0)
    bi(isiupd,0) = bi(1,0)
    bi(1,0)=temp

    call DSWAP(this%np,this%xi,1,this%xi(1,isiupd),1)
    call DGEMV ('N',this%np,nVec1,1d0,this%xi(1,2),this%np,&
         bi(2,0),1,bi(1,0),this%xi(1,1),1)

    call DSWAP(this%np,this%ri,1,this%ri(1,isiupd),1)
    call DGEMV ('N',this%np,nVec1+1,this%delta,this%ri,this%np,&
         bi(1,0),1,1d0,this%xi,1)

    !.................. reload overlaps of current point ...................
    call DCOPY(this%nVec,this%overlap(:,1),1,this%overlap(:,isiupd),1)
    call DCOPY(this%nVec,this%overlap(1,1),this%nVec,this%overlap(isiupd,1),this%nVec)
    call rism_timer_stop(this%projectTimer)
  end subroutine mdiis_blas_advance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Checks if a restart is require (rms1 > this%restart * smallest rms in memory).
!!!If so, the lowest rms vector is chosen and all others are discarded.
!!!this  :: mdiis object
!!!rms1  :: residual at the end of the step
!!!isiupd:: vector that will be updated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine checkRestart (this,rms1,isiupd)
    implicit none
    type(mdiis_blas),intent(inout) :: this
    _REAL_,intent(in) :: rms1
    integer, intent(out):: isiupd

    _REAL_ ::  rmsmin,rmsmi2
    integer :: isirst, isimin,isimi2, is

    !............... get two residuals, minimal in DIIS set ................
    rmsmin = this%overlap(1,1)
    isimin = 1
    rmsmi2 = rmsmin
    isimi2 = isimin
    do is=2,this%nVec0
       if (this%overlap(is,is) < rmsmin)  then
          rmsmi2 = rmsmin
          isimi2 = isimin
          rmsmin = this%overlap(is,is)
          isimin = is
       elseif (this%overlap(is,is) < rmsmi2)  then
          rmsmi2 = this%overlap(is,is)
          isimi2 = is
       endif
    enddo
    rmsmin = sqrt( rmsmin/this%np/this%mpisize)
    rmsmi2 = sqrt( rmsmi2/this%np/this%mpisize)

    !................ if convergence is poor, restart MDIIS ................
    if (this%nVec0 > 1 .AND. rms1 > this%restart*rmsmin)  then

       !.......... choose restarting vector so as to prevent cycling .........

       !if this%nVec0 = 1 then we have problems
       !if this%nVec0 = 2 ditch the old vector and just use the new vector
       !otherwise we use the minimal residue vector
       
       if(this%nVec0 == 1)then
       elseif(this%nVec0 == 2)then
          isirst = 1
       else
          isirst = isimin
       end if

       !................... restore vector to restart from ...................
       if (isirst /= 1)  then
          call timer_start(TIME_MDIIS_DATA)
          call DCOPY(this%np,this%ri(1:this%np,isirst),1,this%ri(1:this%np,1),1)
          call DCOPY(this%np,this%xi(1:this%np,isirst),1,this%xi(1:this%np,1),1)
          call timer_stop(TIME_MDIIS_DATA)
          this%overlap(1,1) = this%overlap(isirst,isirst)
       endif

       !.................... reset MDIIS vectors counters .....................
       this%nVec0 = 1
       this%istep = 1
       isiupd = mod(this%istep-1,this%nVec-1)+2
    endif
  end subroutine checkRestart
end module mdiis_blas_c
