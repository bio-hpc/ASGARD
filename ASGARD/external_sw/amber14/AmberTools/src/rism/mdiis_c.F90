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
!!Generic object interface for various MDIIS implementations.  A
!!specific implementation is chosen by the user at run time via the
!!constructor.  All subroutines and functions should then be called in
!!the same manner with only details of execution being different.
!!
!!MDIIS uses solutions and residual (whatever needs to be minimized
!!for each data point) from previous iterations to predict a new
!!solution with a lower residual.  For all methods, working memory for
!!the data and the residual must be allocated and passed to the
!!constructor.  This memory must not be reallocated/deallocated before
!!the MDIIS object is destroyed.  Data for solutions and residuals are
!!vectors and the memory for a data vector may be of any rank (in
!!practise > 2).  The total memory size should be the vector length *
!!the number of past solutions you wish the MDIIS routine to use.

!!Depending on the implementation chosen, the active vector to
!!read/write data from may change location in working memory after
!!calling MDIIS_ADVANCE(). Use MDIIS_GETVECTORINDEX() to get the
!!active vector number.
!!
!!To add a new MDIIS implementation a new class must be written and
!!registered here and in the Makefile.  Within this class, the new
!!implementation must be registered with a 'USE' statement, added as a
!!type (with enumeration constant) and appropriate calls made from
!!each subroutine/function.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mdiis_c
  !add new module with renaming here
  use mdiis_orig_c
  use mdiis_blas_c
  use mdiis_blas2_c
  use rism_timer_c
  implicit none

  type mdiis
     !enumerated values, add yours here
     integer  :: KOVALENKO=0, KOVALENKO_OPT=1, KOVALENKO_OPT2=2
     !pointer to MDIIS objects, add yours here
     type(mdiis_orig),pointer :: orig=>null()
     type(mdiis_blas),pointer :: blas=>null()
     type(mdiis_blas2),pointer :: blas2=>null()
     type(rism_timer) :: advanceTimer, initTimer
  end type mdiis

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
!!!              minimum residual in the basis that causes a restart
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_new (this,method,delta,tol, restart)
    implicit none
    type(mdiis),intent(inout) :: this
    integer, intent(in) :: method
    _REAL_,intent(in) :: delta, tol, restart

    _REAL_ ::rms
    logical :: con
    call rism_timer_new(this%initTimer,"MDIIS Setup")
    call rism_timer_start(this%initTimer)
    call rism_timer_new(this%advanceTimer,"MDIIS Advance")

    if(method == this%KOVALENKO)then
       allocate(this%orig)
       call mdiis_orig_new_serial(this%orig,delta,tol,restart,&
           this%advanceTimer)
    elseif(method == this%KOVALENKO_OPT)then
       allocate(this%blas)
       call mdiis_blas_new_serial(this%blas,delta,tol,restart,&
           this%advanceTimer)
    elseif(method == this%KOVALENKO_OPT2)then
       allocate(this%blas2)
      call mdiis_blas2_new_serial(this%blas2,delta,tol,restart,&
           this%advanceTimer)
    end if
    !add elsif for your method here
    call rism_timer_stop(this%initTimer)
  end subroutine mdiis_new

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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_new_mpi (this,method,delta,tol, restart, &
       rank, size, comm)
    implicit none
    type(mdiis),intent(out) :: this
    integer, intent(in) :: method
    _REAL_,intent(in) :: delta, tol, restart
    integer, intent(in) :: rank, size, comm
    call rism_timer_new(this%initTimer,"MDIIS Setup")
    call rism_timer_start(this%initTimer)
    call rism_timer_new(this%advanceTimer,"MDIIS Advance")

    if(method == this%KOVALENKO)then
       allocate(this%orig)
       call mdiis_orig_new_mpi(this%orig,delta,tol,restart,&
            this%advanceTimer, rank,size,comm)
    elseif(method == this%KOVALENKO_OPT)then
       allocate(this%blas)
       call mdiis_blas_new_mpi(this%blas,delta,tol,restart,&
           this%advanceTimer,rank,size,comm)
    elseif(method == this%KOVALENKO_OPT2)then
       allocate(this%blas2)
       call mdiis_blas2_new_mpi(this%blas2,delta,tol,restart,&
            this%advanceTimer,rank,size,comm)
    end if
    !add elsif for your method here
    call rism_timer_stop(this%initTimer)
  end subroutine mdiis_new_mpi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Destroy MDIIS object. All internal variables and pointers are
!!!reset. Working memory remains intact memory is swapped (if
!!!necessary) such that xi(:,1) and ri(:,1) become the active vectors.
!!!IN:
!!!   this :: mdiis object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_destroy(this)
    implicit none
    type(mdiis),intent(inout) :: this
    if(associated(this%orig))then
       call mdiis_orig_destroy(this%orig)
       deallocate(this%orig)
    end if
    if(associated(this%blas))then
       call mdiis_blas_destroy(this%blas)
       deallocate(this%blas)
    end if
    if(associated(this%blas2))then
       call mdiis_blas2_destroy(this%blas2)
       deallocate(this%blas2)
    end if
    !add elsif for your method here
    call rism_timer_destroy(this%initTimer)
    call rism_timer_destroy(this%advanceTimer)
  end subroutine mdiis_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Reset MDIIS object. Progress is reset to a single basis
!!!vector.Working memory remains intact memory is swapped (if
!!!necessary) such that xi(:,1) and ri(:,1) become the active vectors.
!!!All other variables are untouched.
!!!IN:
!!!   this :: mdiis object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_reset(this)
    implicit none
    type(mdiis),intent(inout) :: this
    if(associated(this%orig))then
       call mdiis_orig_reset(this%orig)
    end if
    if(associated(this%blas))then
       call mdiis_blas_reset(this%blas)
    end if
    if(associated(this%blas2))then
       call mdiis_blas2_reset(this%blas2)
    end if
    !add elsif for your method here
  end subroutine mdiis_reset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets data vectors and resets MDIIS object. Progress is reset to a single basis
!!!vector. Working memory remains intact memory is swapped (if
!!!necessary) such that xi(:,1) and ri(:,1) become the active vectors.
!!!All other variables are untouched.
!!!IN:
!!!   this :: mdiis object
!!!   xi    :: array of vector data.  np X nvec
!!!   ri    :: array of residual data.  np X nvec
!!!   np    :: number of data points (first dimension length)
!!!   nvec   :: length of DIIS vectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_setData(this,xi,ri,np,nvec)
    implicit none
    type(mdiis),intent(inout) :: this
    _REAL_,target, intent(in) :: xi(np,nvec), ri(np,nvec)
    integer,intent(in) ::  np
    integer, intent(in) :: nvec
    call mdiis_resize(this,xi,ri,np,nvec)
  end subroutine mdiis_setData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Resizes and resets MDIIS object. Progress is reset to a single basis
!!!vector. Working memory remains intact memory is swapped (if
!!!necessary) such that xi(:,1) and ri(:,1) become the active vectors.
!!!All other variables are untouched.
!!!IN:
!!!   this :: mdiis object
!!!   xi    :: array of vector data.  np X nvec
!!!   ri    :: array of residual data.  np X nvec
!!!   np    :: number of data points (first dimension length)
!!!   nvec   :: length of DIIS vectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_resize(this,xi,ri,np,nvec)
    implicit none
    type(mdiis),intent(inout) :: this
    _REAL_,target, intent(in) :: xi(np,nvec), ri(np,nvec)
    integer,intent(in) ::  np
    integer, intent(in) :: nvec
    if(associated(this%orig))then
       call mdiis_orig_resize(this%orig,xi,ri,np,nvec)
    end if
    if(associated(this%blas))then
       call mdiis_blas_resize(this%blas,xi,ri,np,nvec)
    end if
    if(associated(this%blas2))then
       call mdiis_blas2_resize(this%blas2,xi,ri,np,nvec)
    end if
    !add elsif for your method here
  end subroutine mdiis_resize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Set parent(s) for this timer
!!!IN:
!!!   this : rism1d object
!!!   parentAdvance : parent for the advance subroutine timer
!!!   parentOther   : (optional) parent for everything else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mdiis_setTimerParent(this, parentAdvance, parentOther)
    implicit none
    type(mdiis), intent(inout) :: this
    type(rism_timer), intent(inout) :: parentAdvance
    type(rism_timer), optional, intent(inout) :: parentOther
    call rism_timer_start(this%initTimer)
    call rism_timer_setParent(this%advanceTimer,parentAdvance)
    if(present(parentOther))then
       call rism_timer_setParent(this%initTimer,parentOther)
    end if
    call rism_timer_stop(this%initTimer)
  end subroutine mdiis_setTimerParent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the number of vector currently used as a basis by the MDIIS object
!!!IN:
!!!   this :: mdiis object
!!!OUT:
!!!   The number of vectors currently used
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function getCurrentNVec (this) result(nVec)
    implicit none
    type(mdiis),intent(inout) :: this
    integer :: nVec
    call rism_timer_start(this%initTimer)
    if(associated(this%orig))then
       nVec=mdiis_orig_getCurrentNVec(this%orig)
    elseif(associated(this%blas))then
       nVec=mdiis_blas_getCurrentNVec(this%blas)
    elseif(associated(this%blas2))then
       nVec=mdiis_blas2_getCurrentNVec(this%blas2)
    end if
    !add elsif for your method here
    call rism_timer_stop(this%initTimer)
  end function getCurrentNVec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the index number of the working vector
!!!IN: 
!!!   this :: mdiis object
!!!OUT:
!!!    the number of the vector to read/write
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function mdiis_getWorkVector (this) result(index)
    implicit none
    type(mdiis),intent(inout) :: this
    integer :: index
    call rism_timer_start(this%initTimer)
    if(associated(this%orig))then
       index=mdiis_orig_getWorkVector(this%orig)
    elseif(associated(this%blas))then
       index=mdiis_blas_getWorkVector(this%blas)
    elseif(associated(this%blas2))then
       index=mdiis_blas2_getWorkVector(this%blas2)
    end if
    !add elsif for your method here
    call rism_timer_stop(this%initTimer)
  end function mdiis_getWorkVector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!One MDIIS iteration. Query mdiis_getVectorNumer() to get the new active vector.
!!!IN:
!!!   this  :: mdiis object
!!!   rms1  :: residual at the end of the step
!!!   conver:: have we converged
!!!   tolerance :: tolerance for this calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  mdiis_advance (this, rms1,conver,tolerance_o)
    implicit none 
    type(mdiis),intent(inout) :: this
    _REAL_, intent(out) ::  rms1
    logical,intent(inout) ::  conver
    _REAL_,optional,intent(in) ::  tolerance_o
    call rism_timer_start(this%advanceTimer)
#ifdef RISM_DEBUG
    write(6,*)"MDIIS_advance generic",rms1,conver
#endif

    if(associated(this%orig))then
       if(present(tolerance_o))then
          call mdiis_orig_advance(this%orig,rms1,conver,tolerance_o)
       else
          call mdiis_orig_advance(this%orig,rms1,conver)
       end if
    elseif(associated(this%blas))then
       if(present(tolerance_o))then
          call mdiis_blas_advance(this%blas,rms1,conver,tolerance_o)
       else
          call mdiis_blas_advance(this%blas,rms1,conver)
       end if
    elseif(associated(this%blas2))then
       if(present(tolerance_o))then
          call mdiis_blas2_advance(this%blas2,rms1,conver,tolerance_o)
       else
          call mdiis_blas2_advance(this%blas2,rms1,conver)
       end if
    end if
    !add elsif for your method here
    call rism_timer_stop(this%advanceTimer)
  end subroutine mdiis_advance

end module mdiis_c
