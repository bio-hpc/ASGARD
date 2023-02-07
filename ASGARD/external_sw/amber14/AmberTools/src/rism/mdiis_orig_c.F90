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
!!
!!This is the reference MDIIS implementation.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mdiis_orig_c
  use safemem
  use rism_timer_c
  implicit none
  type mdiis_orig
     private
     integer ::  istep = 0, nis0 = 0, isimio=0
     integer :: nis, np
     _REAL_,dimension(:,:),pointer :: ovlij=>NULL(), tovlij=>NULL()
     integer :: mpirank=0, mpisize=1, mpicomm=0
     !delta :: coefficient for residual gradient.  Should be between 0 and 1.
     !tol   :: target residual
     _REAL_ :: delta, tol
     !pointers to working data
     !xi    :: array of vector data.  np X nis
     !ri    :: array of residual data.  np X nis
     _REAL_,pointer :: xi(:,:)=>NULL(), ri(:,:)=>NULL()
     !restart restarting comparison factor .....................
     _REAL_ ::  restart

     type(rism_timer) :: overlapTimer, restartTimer, &
          lapackTimer, projectTimer

  end type mdiis_orig

  interface mdiis_orig_new
     module procedure mdiis_orig_new_serial, mdiis_orig_new_mpi
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Create new MDIIS object.    Provide convergence parameters and
!!!working memory.  Note that the working memory must me nis times
!!!larger than the data array.  Further, the number of data points
!!!cannot change between iterations.  xi(:,1) and ri(:,1) are the active vectors.
!!!IN:
!!!   this  :: mdiis object
!!!   delta :: coefficient for residual gradient.  Should be between 0 and 1.
!!!   tol   :: target residual
!!!   xi    :: array of vector data.  np X nis
!!!   ri    :: array of residual data.  np X nis
!!!   np    :: number of data points (first dimension length)
!!!   nis   :: length of DIIS vectors
!!!   restart :: restart threshold factor. Ratio of the current residual to the 
!!!              minimum residual in the basis that causes a restart
!!!   timer :: parent rism_timer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_orig_new_serial (this,delta,tol,restart,timer)
    implicit none
    type(mdiis_orig),intent(out) :: this
    _REAL_,intent(in) :: delta, tol, restart
    type(rism_timer), intent(inout) :: timer

    this%delta = delta
    this%tol = tol
    this%restart = restart
    this%istep = 0
    call rism_timer_new(this%overlapTimer,"Overlap")
    call rism_timer_setParent(this%overlapTimer,timer)
    call rism_timer_new(this%restartTimer,"Restart")
    call rism_timer_setParent(this%restartTimer,timer)
    call rism_timer_new(this%lapackTimer,"Lapack")
    call rism_timer_setParent(this%lapackTimer,timer)
    call rism_timer_new(this%projectTimer,"Project")
    call rism_timer_setParent(this%projectTimer,timer)
  end subroutine mdiis_orig_new_serial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Create new MDIIS object.  Provide convergence parameters, working
!!!memory and MPI parameters.  Note that the working memory must me
!!!nis times larger than the data array.  Further, the number of data
!!!points cannot change between iterations.  xi(:,1) and ri(:,1) are
!!!the active vectors.
!!!IN:
!!!   this  :: mdiis object
!!!   delta :: coefficient for residual gradient.  Should be between 0 and 1.
!!!   tol   :: target residual
!!!   xi    :: array of vector data.  np X nis
!!!   ri    :: array of residual data.  np X nis
!!!   np    :: number of data points (first dimension length)
!!!   nis   :: length of DIIS vectors
!!!   restart :: restart threshold factor. Ratio of the current residual to the 
!!!              minimum residual in the basis that causes a restart
!!!   rank  :: MPI process rank
!!!   size  :: Number of processes
!!!   comm  :: MPI communicator
!!!   timer :: parent rism_timer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_orig_new_mpi(this,delta,tol,restart,timer,&
       rank,size,comm)
    implicit none
    type(mdiis_orig),intent(out) :: this
    _REAL_,intent(in) :: delta, tol, restart
    integer, intent(in) :: rank, size, comm
    type(rism_timer), intent(inout) :: timer
#ifdef RISM_DEBUG
    write(6,*) "MDIIS_ORIG_NEW_MPI",rank,size,comm
#endif /*RISM_DEBUG*/    
    call mdiis_orig_new_serial(this,delta,tol,restart,timer)
    this%mpirank = rank
    this%mpisize = size
    this%mpicomm = comm
#ifdef RISM_DEBUG
    write(6,*) "MDIIS_ORIG_NEW_MPI",this%mpirank,this%mpisize,this%mpicomm
#endif /*RISM_DEBUG*/    
  end subroutine mdiis_orig_new_mpi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Destroys the MDIIS object.  All variables are reset except for
!!!working memory.
!!!IN:
!!!   this :: mdiis object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_orig_destroy(this)
    implicit none
    type(mdiis_orig),intent(inout) :: this
    integer :: err
    call mdiis_orig_reset(this)
    this%delta=0
    this%tol=0
    this%np=0
    this%nis=0
    if(safemem_dealloc(this%ovlij) /= 0) call rism_report_error("MDIIS_ORIG: failed to deallocate OVLIJ")
    if(safemem_dealloc(this%tovlij) /= 0) call rism_report_error("MDIIS_ORIG: failed to deallocate TOVLIJ")
    nullify(this%xi)
    nullify(this%ri)
    call rism_timer_destroy(this%overlapTimer)
    call rism_timer_destroy(this%restartTimer)
    call rism_timer_destroy(this%lapackTimer)
    call rism_timer_destroy(this%projectTimer)
  end subroutine mdiis_orig_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Reset the MDIIS object. Progress is reset to a single basis
!!!vector.  All other variables are untouched.
!!!IN:
!!!   this :: mdiis object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_orig_reset(this)
    use safemem
    implicit none
    type(mdiis_orig),intent(inout) :: this
    integer :: err

    this%ovlij=>safemem_realloc(this%ovlij,this%nis,this%nis,.false.)
    this%tovlij=>safemem_realloc(this%tovlij,this%nis,this%nis,.false.)
    this%nis0 = 0
    this%isimio=0
    this%ovlij=0
    this%tovlij=0
    this%istep=0
  end subroutine mdiis_orig_reset

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
  subroutine mdiis_orig_resize(this,xi,ri,np,nvec)
    implicit none
    type(mdiis_orig),intent(inout) :: this
    _REAL_,target, intent(in) :: xi(np,nvec), ri(np,nvec)
    integer,intent(in) ::  np
    integer, intent(in) :: nvec
    !transfer the current working vector to the first index
    this%np = np
    this%nis = nvec
    this%xi => xi
    this%ri => ri
    call mdiis_orig_reset(this)
  end subroutine mdiis_orig_resize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the number of vector currently used as a basis by the MDIIS object
!!!IN:
!!!   this :: mdiis object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function mdiis_orig_getCurrentNVec (this) result(nVec)
    implicit none
    type(mdiis_orig),intent(in) :: this
    integer :: nVec
    nVec = this%nis0
  end function mdiis_orig_getCurrentNVec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the index number of the working vector
!!!IN: 
!!!   this :: mdiis object
!!!OUT:
!!!    the number of the vector to read/write
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function mdiis_orig_getWorkVector (this) result(index)
    implicit none
    type(mdiis_orig),intent(in) :: this
    integer :: index
    index = 1
  end function mdiis_orig_getWorkVector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!One MDIIS iteration
!!!IN:
!!!   this  :: mdiis object
!!!   rms1  :: residual at the end of the step
!!!   conver:: have we converged
!!!   tolerance :: tolerance for this calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  mdiis_orig_advance (this, rms1,conver,tolerance_o)
    use rism_util, only : checksum
    implicit none 
#ifdef MPI
    include 'mpif.h'
#endif /*MPI*/
#include "def_time.h" 
    
    type(mdiis_orig),intent(inout) :: this
    logical ::  conver
    _REAL_ ::  rms1
    _REAL_,optional,intent(in) ::  tolerance_o

    _REAL_ :: tolerance

    integer ::  isiupd, isirst, isimin,isimi2, &
         ip, is,is1,is2, isi,is1i,is2i
    integer ::  isindx(this%nis), indx(0:this%nis)
    _REAL_ ::  rmsmin,rmsmi2, xs,rs
    _REAL_ :: aij(0:this%nis,0:this%nis)
    _REAL_ :: taij(0:this%nis,0:this%nis)
    _REAL_ :: bi(0:this%nis,0:this%nis), detsgn
    _REAL_ :: tbi(0:this%nis,0:this%nis)
    integer :: err=0, ierr

    tolerance = this%tol
    if(present(tolerance_o)) tolerance = tolerance_o
#ifdef RISM_DEBUG
    write(6,*)"MDIIS start",this%delta,rms1,tolerance,this%np,this%nis,conver
#endif

    !.................. initialize counters and switches ...................

    !.................. increment step and MDIIS counters ..................
    this%istep = this%istep + 1
    this%nis0 = min( this%nis0+1, this%nis)

    !........................ specify storage list .........................
    isindx(1) = 1
    do is=2,this%nis
       isindx(is) = 2 + mod( this%istep-is+this%nis-1, this%nis-1)
    enddo
    isiupd = isindx(this%nis)

    !............. calculate diagonal overlap of new residual ..............
    this%ovlij(1,1) = 0d0
    this%tovlij(1,1) = 0d0
    call rism_timer_start(this%overlapTimer)
    call timer_start(TIME_MDIIS_DATA)
#if defined(MPI)      
    do ip=1,this%np
       this%tovlij(1,1) = this%tovlij(1,1) + this%ri(ip,1)**2
    enddo
    CALL MPI_AllREDUCE(this%tovlij(1,1),this%ovlij(1,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,     &
         this%mpicomm,ierr)
    if(err /=0) call rism_report_error&
         ("MDIIS_BLAS2_UPDATE: could not reduce OVERLAP")
#else
    do ip=1,this%np
       this%ovlij(1,1) = this%ovlij(1,1) + this%ri(ip,1)**2
    enddo
#endif /*defined(MPI)*/
    call timer_stop(TIME_MDIIS_DATA)
    call rism_timer_stop(this%overlapTimer)
    !................ get mean square value of new residual ................
    rms1 = sqrt( this%ovlij(1,1)/(this%np*this%mpisize))
#ifdef RISM_DEBUG
    write(6,*) "RMS1",this%mpirank,rms1,checksum(this%ri,product(shape(this%ri)),this%mpicomm),&
         sqrt(checksum(this%ri**2,product(shape(this%ri)),this%mpicomm)/this%np/this%mpisize),&
         checksum(this%ri**2,product(shape(this%ri)),this%mpicomm),this%np,this%mpisize
    write(6,*) "shape(this%ri)",shape(this%ri)
#endif /*RISM_DEBUG*/
    !.............. check mean square residual for tolerance ...............
    if (rms1 <= tolerance)  then
       conver = .true.
       return
    else
       conver = .false.
    endif

    !............... get two residuals, minimal in DIIS set ................
    rmsmin = this%ovlij(1,1)
    isimin = 1
    rmsmi2 = rmsmin
    isimi2 = isimin
    do is=2,this%nis0
       isi = isindx(is)
       if (this%ovlij(isi,isi) < rmsmin)  then
          rmsmi2 = rmsmin
          isimi2 = isimin
          rmsmin = this%ovlij(isi,isi)
          isimin = isi
       elseif (this%ovlij(isi,isi) < rmsmi2)  then
          rmsmi2 = this%ovlij(isi,isi)
          isimi2 = isi
       endif
    enddo
      rmsmin = sqrt( rmsmin/this%np/this%mpisize)
      rmsmi2 = sqrt( rmsmi2/this%np/this%mpisize)

    !................ if convergence is poor, restart MDIIS ................
    call rism_timer_start(this%restartTimer)
    if (this%nis0 > 1 .AND. rms1 > this%restart*rmsmin)  then

       !.......... choose restarting vector so as to prevent cycling .........
       if (isimin == this%isimio .AND. this%nis0 < this%nis)  then
          isirst = isimi2
       else
          isirst = isimin
       endif

       !................... restore vector to restart from ...................
       if (isirst /= 1)  then
          call timer_start(TIME_MDIIS_DATA)
          do ip=1,this%np
             this%ri(ip,1) = this%ri(ip,isirst)
             this%xi(ip,1) = this%xi(ip,isirst)
          enddo
          call timer_stop(TIME_MDIIS_DATA)
          this%ovlij(1,1) = this%ovlij(isirst,isirst)
       endif

       !.................... reset MDIIS vectors counters .....................
       this%nis0 = 1
       this%isimio = isiupd

    endif
    call rism_timer_stop(this%restartTimer)

    !........... calculate nondiagonal overlaps of new residual ............
    call rism_timer_start(this%overlapTimer)
    do is=2,this%nis0
       isi = isindx(is)
       this%ovlij(isi,1) = 0d0
    call timer_start(TIME_MDIIS_DATA)
#if defined(MPI)      
       this%tovlij(isi,1) = 0d0
       do ip=1,this%np
          this%tovlij(isi,1) = this%tovlij(isi,1) + this%ri(ip,isi)*this%ri(ip,1)
       enddo
       CALL MPI_AllREDUCE(this%tovlij(isi,1),this%ovlij(isi,1),1,MPI_DOUBLE_PRECISION,MPI_SUM,     &
            this%mpicomm,ierr)
       if(err /=0) call rism_report_error&
            ("MDIIS_ORIG_UPDATE: could not reduce OVERLAP")
#else
       do ip=1,this%np
          this%ovlij(isi,1) = this%ovlij(isi,1) + this%ri(ip,isi)*this%ri(ip,1)
       enddo
#endif /*defined(MPI)*/
    call timer_stop(TIME_MDIIS_DATA)
       this%ovlij(1,isi) = this%ovlij(isi,1)
    enddo
    call rism_timer_stop(this%overlapTimer)
    !......................... load DIIS matrices ..........................
    call rism_timer_start(this%lapackTimer)
    call timer_start(TIME_MDIIS_LAPACK)
    aij(0,0) = 0d0
    bi=0
    bi(0,0) = -1d0
    do is=1,this%nis0
       aij(is,0) = -1d0 
       aij(0,is) = -1d0
    enddo
    do is2=1,this%nis0
       do is1=1,this%nis0
          is1i = isindx(is1)
          is2i = isindx(is2)
          aij(is1,is2) = this%ovlij(is1i,is2i)
       enddo
    enddo

    !....................... calculate DIIS estimate .......................
    call DGETRF(this%nis0+1,this%nis0+1,aij,this%nis+1,indx,err)
    if(err > 0)then
       call rism_report_error("LU-factorization failed.  U = 0")
    elseif(err<0)then
       err = err*(-1)
       call rism_report_error("LU-factorization failed.")
    endif
    call DGETRS('N',this%nis0+1,this%nis0+1,aij,this%nis+1,indx,bi,this%nis+1,err)
    if(err < 0)then
       err = err*(-1)
       call rism_report_error("Linear equation solver failed.")
    endif
    call timer_stop(TIME_MDIIS_LAPACK)
    call rism_timer_stop(this%lapackTimer)

!    call timer_start(TIME_MDIIS_DATA)
    !......... get DIIS minimum, MDIIS correction, and next point ..........
    call rism_timer_start(this%projectTimer)
    do ip=1,this%np
       xs = 0d0
       rs = 0d0
       do is=1,this%nis0
          isi = isindx(is)
          xs = xs + bi(is,0)*this%xi(ip,isi)
          rs = rs + bi(is,0)*this%ri(ip,isi)
       enddo
       this%xi(ip,isiupd) = this%xi(ip,1)
       this%ri(ip,isiupd) = this%ri(ip,1)
       this%xi(ip,1) = xs + this%delta*rs
    enddo
!    call timer_stop(TIME_MDIIS_DATA)

    !.................. reload overlaps of current point ...................
    do is=1,this%nis
       isi = isindx(is)
       this%ovlij(isi,isiupd) = this%ovlij(isi,1)
    enddo

    do is=1,this%nis
       isi = isindx(is)
       this%ovlij(isiupd,isi) = this%ovlij(1,isi)
    enddo
    call rism_timer_stop(this%projectTimer)
  end subroutine mdiis_orig_advance
end module mdiis_orig_c
