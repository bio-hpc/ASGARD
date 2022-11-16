!The 3D-RISM-KH software found here is copyright (c) 2012 by 
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
!!!Debug module for 3D-RISM.  By using MPI, grid and solvent data, it
!!!allows printing out of arrays and other data from any part of the
!!!program.  It is also useful for comparing MPI and serial as arrays
!!!can be printed in the correct order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  module rism3d_debug_c
    use rism3d_grid_c
    use rism3d_solv_c
    implicit none
    private
    
    type(rism3d_grid),pointer :: grid=>NULL()
    type(rism3d_solv), pointer :: solv=>NULL()
    integer :: mpirank, mpisize,mpicomm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!prints data w/ or w/o FFT padding.  See subroutines below for details.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    interface rism3d_debug_print
       module procedure rism3d_debug_print_with_k
       module procedure rism3d_debug_print_r
       module procedure rism3d_debug_print_r1
       module procedure rism3d_debug_print_r2
    end interface
    public rism3d_debug_new, rism3d_debug_print
  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!create a new instance.  There is only one global instance of this module.
!!!IN:
!!!   grid3d :: grid object for 3D-RISM
!!!   solvent :: solvent object for 3D-RISM
!!!   rank :: mpirank
!!!   size :: mpisize
!!!   comm :: mpicomm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_debug_new(grid3d,solvent,rank, size, comm)
      implicit none
      type(rism3d_grid),target :: grid3d
      type(rism3d_solv),target :: solvent
      integer, intent(in) :: rank, size, comm
      grid => grid3d
      solv => solvent
      mpirank = rank
      mpisize = size
      mpicomm = comm
    end subroutine rism3d_debug_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Print an array with k-space padding (e.g. guv or huv).  To get the
!!!correct output order, specify if it is in r-space or not and if it
!!!is in Numerical Recipes FFTlayout or not
!!!IN:
!!!   data :: nk X nsolv
!!!   rspace :: .true. for real data
!!!   label :: text printed on each line of output
!!!   nrformat :: .true. if Numerical Recipes FFT memory layout
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_debug_print_with_k(data,rspace,label, nrformat)
      implicit none
      _REAL_, intent(in) :: data(:,:)
      logical, intent(in) :: rspace, nrformat
      character(len=*), intent(in) :: label
      !padding :: print FFT padding indices for r-space data
      logical :: padding
      integer :: iv, ix, iy, iz, ir, dim(3), irank
      integer :: offset(3)
      integer :: ierr
      character(len=40) :: fmt
      fmt = '(a,i3,i3,i5,i4,i4,i4,1p,e24.16)'
      dim = grid%ngr
      padding = .false.

      if(.not.nrformat)&
           dim(1) = dim(1)+2

      if(rspace)then
         dim(3) = dim(3)/mpisize
         offset = grid%nroff
      else
         dim(2) = dim(2)/mpisize
         offset = grid%nkoff
      endif

      if(rspace)then
         !case of r-space data
         do iv =1, solv%natom
            do irank=0, mpisize-1
#ifdef MPI
               !this method is not guaranteed to work and should be
               !changed to sending data to the master node for
               !printing
               call mpi_barrier(mpicomm,ierr)
#endif
               if(irank == mpirank)then
                  do iz=1, dim(3)
                     do iy = 1, dim(2)
                        do ix =1, dim(1)
                           ir = ix + (iy-1)*(dim(1)) + (iz-1)*(dim(1)*dim(2))
                           if(.not. padding .and. .not.nrformat .and. &
                                (ix == dim(1) .or. ix == dim(1)-1)) cycle
                           write(0,fmt) trim(label),irank,iv,ir,&
                                ix+irank*offset(1), iy+irank*offset(2), iz+irank*offset(3), &
                                data(ir,iv)
                           call flush(0)
                        end do
                        if(nrformat .and. padding)then
                           ir=product(dim) + 1 + (iy-1)*2 + (iz-1)*dim(2)*2
                           write(0,fmt) trim(label),irank,iv,ir,&
                                dim(1)+1+irank*offset(1), iy+irank*offset(2), iz+irank*offset(3), &
                                data(ir,iv)                              
                           write(0,fmt) trim(label),irank,iv,ir+1,&
                                dim(1)+2+irank*offset(1), iy+irank*offset(2), iz+irank*offset(3), &
                                data(ir+1,iv)                              
                        end if
                        call flush(0)
                     end do
                  end do
               end if
#ifdef MPI
               call mpi_barrier(mpicomm,ierr)
#endif
            end do
         end do
      else
         !case of k-space data
         do iv =1, solv%natom
            do iz=1, dim(3)
               do irank=0, mpisize-1
#ifdef MPI
                  call mpi_barrier(mpicomm,ierr)
#endif
                  if(irank == mpirank)then
                     do iy = 1, dim(2)
                        do ix =1, dim(1)
                           ir = ix + (iz-1)*(dim(1)) + (iy-1)*(dim(1)*dim(3))
                           write(0,fmt) trim(label),irank, iv,ir,&
                                ix+irank*offset(1), iy+irank*offset(3), iz+irank*offset(2), &
                                data(ir,iv)
                           call flush(0)
                        end do
                     end do
                  end if
#ifdef MPI
                  call mpi_barrier(mpicomm,ierr)
#endif
               end do
            end do
         end do

      end if

    end subroutine rism3d_debug_print_with_k
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Print real data array
!!!IN:
!!!   data :: nx X ny X nz X nsolv
!!!   label :: text printed on each line of output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_debug_print_r(data,label)
      use safemem
      implicit none
#ifdef MPI
      include 'mpif.h'
      integer :: ierr, status(MPI_STATUS_SIZE)
#endif
      _REAL_, intent(in) :: data(:,:,:,:)
      character(len=*), intent(in) :: label
      integer :: iv, ix, iy, iz, ir, dim(3), irank
      integer :: offset(3)
      character(len=40) :: fmt
      _REAL_,pointer :: temp(:,:,:)=> NULL()
      fmt = '(a,i3,i3,i4,i4,i4,1p,e24.16)'
      dim = grid%ngr
      dim(3) = dim(3)/mpisize
      offset = (/0,0,dim(3)/)
      temp => safemem_realloc(temp,dim(1),dim(2),dim(3),.false.)
      do iv =1, solv%natom
#ifdef MPI
         !non-head node process send data to the master
         if(mpirank /= 0)then
            call MPI_Send(data(:,:,:,iv), product(dim), MPI_DOUBLE_PRECISION, 0, 0,&
                 mpicomm, ierr)
         endif
#endif
         if(mpirank == 0)then
            !master node either uses local data or receives data to print
            do irank =0, mpisize-1
               if(irank == 0) then
                  temp=data(:,:,:,iv)
               else
#ifdef MPI
                  call MPI_Recv(temp, product(dim), MPI_DOUBLE_PRECISION, irank, 0, &
                       mpicomm, status, ierr)
#endif
               end if
               do iz=1, dim(3)
                  do iy = 1, dim(2)
                     do ix =1, dim(1)
                        write(0,fmt) trim(label),irank,iv,&
                             ix+irank*offset(1), iy+irank*offset(2), iz+irank*offset(3), &
                             temp(ix,iy,iz)
                        call flush(0)
                     end do
                  end do
               end do
            end do
         end if
      end do
      if(safemem_dealloc(temp) /= 0)&
           call rism_report_error("RISM3D_DEBUG_PRINT_R: Deallocating temp failed")
    end subroutine rism3d_debug_print_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Print real data array
!!!IN:
!!!   data :: (nx X ny X nz X nsolv)
!!!   label :: text printed on each line of output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_debug_print_r1(data,label)
      use safemem
      implicit none
#ifdef MPI
      include 'mpif.h'
      integer :: ierr, status(MPI_STATUS_SIZE)
#endif
      _REAL_, target, intent(in) :: data(:)
      character(len=*), intent(in) :: label
      _REAL_, pointer :: temp(:,:,:,:) => NULL()

      temp(1:grid%nr(1),1:grid%nr(2),1:grid%nr(3),1:solv%natom) => data
      call rism3d_debug_print(temp,label)
    end subroutine rism3d_debug_print_r1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Print real data array
!!!IN:
!!!   data :: (nx X ny X nz) X nsolv
!!!   label :: text printed on each line of output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_debug_print_r2(data,label)
      use safemem
      implicit none
#ifdef MPI
      include 'mpif.h'
      integer :: ierr, status(MPI_STATUS_SIZE)
#endif
      _REAL_, target, contiguous, intent(in) :: data(:,:)
      character(len=*), intent(in) :: label
      _REAL_, pointer :: temp(:,:,:,:) => NULL()

      temp(1:grid%nr(1),1:grid%nr(2),1:grid%nr(3),1:solv%natom) => data(:,:)
      call rism3d_debug_print(temp,label)
    end subroutine rism3d_debug_print_r2
  end module rism3d_debug_c
