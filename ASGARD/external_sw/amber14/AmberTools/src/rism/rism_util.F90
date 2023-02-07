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
module rism_util
  use rism_report_c
  contains
    subroutine corr_drift(ff,mass,natu &
#ifdef MPI
         ,rank,size,comm &
#endif /*MPI*/
         )
      implicit none
#if MPI
      include 'mpif.h'
      integer, intent(in) :: rank,size,comm
#endif /*MPI*/
      integer, intent(in) :: natu
      _REAL_,intent(inout) :: ff(3,natu)
      _REAL_,intent(in) :: mass(natu)
      
      integer :: id,iatu
      _REAL_ :: totmass,totfrc(3)
      _REAL_ :: rel_drift(3),drift,magtotfrc,norm(3)
      integer :: err
      _REAL_ :: mpitmp(3)
#ifdef RISM_DEBUG
      write(6,*) "CORR_DRIFT"
#endif /*RISM_DEBUG*/
      do id = 1, 3
         totfrc(id) = sum(ff(id,1:natu))
#ifdef RISM_DEBUG
         write (6,*) "ID=",id
         write (6,*) ff(id,1:natu)
#endif /*RISM_DEBUG*/
      end do
#if defined(MPI) && defined(MPI)
#  ifdef USE_MPI_IN_PLACE
      call MPI_ALLREDUCE(MPI_IN_PLACE,totfrc,3,MPI_DOUBLE_PRECISION,&
           MPI_SUM,comm,err)
#  else
      call MPI_ALLREDUCE(totfrc,mpitmp,3,MPI_DOUBLE_PRECISION,&
           MPI_SUM,comm,err)      
      totfrc = mpitmp
#  endif /*USE_MPI_IN_PLACE*/
      if(err /=0) call rism_report_error&
           ("RISM3D CORR_DRIFT: could not reduce TOTFRC")
#endif /*defined(MPI)*/

#ifdef RISM_DEBUG
      write(6,*)"checking drift..."
      write(6,*)"x",totfrc(1)
      write(6,*)"y",totfrc(2)
      write(6,*)"z",totfrc(3)
      write(6,*)"corrected drift..."
#endif /*RISM_DEBUG*/

!
!The force is distributed across all so the correction must be too.
!I.e., the correction is divided by the number of processors
!

      totmass = sum(mass)
      magtotfrc = sqrt(sum(totfrc**2))
      do iatu=1,natu
#if defined(MPI)
         ff(:,iatu) = ff(:,iatu) - totfrc*mass(iatu)/totmass/size
#else
         ff(:,iatu) = ff(:,iatu) - totfrc*mass(iatu)/totmass
#endif /*defined(MPI)*/
      end do
#ifdef RISM_DEBUG
      do id = 1, 3
         totfrc(id) = sum(ff(id,1:natu))
      end do      
#  if defined(MPI)
#    ifdef USE_MPI_IN_PLACE
      call MPI_ALLREDUCE(MPI_IN_PLACE,totfrc,3,MPI_DOUBLE_PRECISION,&
           MPI_SUM,comm,err)
#    else
      call MPI_ALLREDUCE(totfrc,mpitmp,3,MPI_DOUBLE_PRECISION,&
           MPI_SUM,comm,err)      
      totfrc = mpitmp
#    endif /*USE_MPI_IN_PLACE*/
      if(err /=0) call rism_report_error&
           ("RISM3D CORR_DRIFT: could not reduce TOTFRC")
#  endif /*defined(MPI)*/
      write(6,*)"x",totfrc(1)
      write(6,*)"y",totfrc(2)
      write(6,*)"z",totfrc(3)
      write(6,*)'outputting ff to file:  ala.ff'
      !#endif
#endif /*RISM_DEBUG*/
    end subroutine corr_drift

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate the center of mass
!IN:
!   ratu   :: the x,y,z position of each solute atom.
!   ratucm :: will hold the center of mass
!   mass   :: mass of each atom
!   natu   :: the number of solute atoms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_cm(ratu,ratucm,mass,natu)
  implicit none
  integer,intent(in) :: natu
  _REAL_,intent(in) :: ratu(3,natu),mass(natu)
  _REAL_,intent(out) :: ratucm(3)
  integer :: id
#ifdef RISM_DEBUG
  write(6,*)"CALC CM", natu
  write(6,*)"RATU", ratu
  write(6,*)"MASS", mass
  call flush(6)
#endif /*RISM_DEBUG*/
  ratucm=0
  do id=1,3
     ratucm(id) = sum(ratu(id,1:natu)*mass(1:natu))
  end do
  ratucm = ratucm/sum(mass)
#ifdef RISM_DEBUG
  write(6,*)"ratucm", ratucm
  call flush(6)
#endif /*RISM_DEBUG*/
end subroutine calc_cm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the cross product of A and B and places it in C.  It is assumed
!!!that A, B and C are of length three.
!!!IN:
!!!   a :: three element array
!!!   b :: three element array
!!!   c :: three element array, on return contains the cross product of A and B
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cross(a,b,c)
  implicit none
  _REAL_,intent(in) :: a(3),b(3)
  _REAL_,intent(out) :: c(3)
  c(1) = a(2)*b(3)-a(3)*b(2)
  c(2) =-a(1)*b(3)+a(3)*b(1)
  c(3) = a(1)*b(2)-a(2)*b(1)  
end subroutine cross

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Translate the system by the given translation vector.
!IN:
!   ratu   :: the x,y,z position of each solute atom.  This is modified.
!   natu   :: the number of solute atoms
!   trans  :: translation vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine translate(ratu,natu,trans)
  implicit none
  integer,intent(in) :: natu
  _REAL_,intent(inout) :: ratu(3,natu)
  _REAL_, intent(in) :: trans(3)

  integer :: iatu

  do iatu = 1, natu
     ratu(1:3,iatu) = ratu(1:3,iatu) + trans
  end do
end subroutine translate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculates the rotational velocity and corresponding kinetic energy 
!using finite difference between two frames.
!IN:
!   ratu1   :: the current x,y,z position of each solute atom.
!   ratu0   :: the previous x,y,z position of each solute atom.
!   mass    :: mass of each atom
!   natu    :: the number of solute atoms
!   dt      :: difference in time between the two frames
!OUT:
!    three element array of the rotational velocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function calc_rotation(ratu1,ratu0,mass,natu,dt)
  implicit none
  _REAL_,intent(in) :: ratu1(3,natu),ratu0(3,natu),mass(natu),dt
  integer,intent(in) :: natu
  _REAL_ :: calc_rotation(3)

  _REAL_ :: cm0(3),cm1(3),angvel(3),r0(3),r1(3), erot,rxv(3), v(3,natu), proj(3),dir(3)
  _REAL_ :: moi(6)
  integer :: iatu,id,info
  integer :: ipiv(3)
#ifdef RISM_DEBUG
  write(6,*) "CALC_ROTATION",dt
#endif /*RISM_DEBUG*/
  call calc_cm(ratu0,cm0,mass,natu)
  call calc_cm(ratu1,cm1,mass,natu)
  !calculate the total angular momentum (place in angvel)
  angvel=0d0
  do iatu=1,natu
     r0 = ratu0(:,iatu)-cm0
     r1 = ratu1(:,iatu)-cm1
!!$     v(:,iatu)=(r1-r0)/(dt*20.455d0)
     v(:,iatu)=(r1-r0)/(dt)
     
     call cross(r1,v(:,iatu),rxv)
!!$     angvel = angvel + mass(iatu)*rxv/sum(r1**2)
     angvel = angvel + mass(iatu)*rxv
  end do
#ifdef RISM_DEBUG
  write(6,*) "ANGULAR MOMENTUM", angvel
  write(6,'(3(g16.3))') v
#endif /*RISM_DEBUG*/
  moi = calc_moi(ratu1,mass)
#ifdef RISM_DEBUG
  write(6,*) "MOI", moi
#endif /*RISM_DEBUG*/
  call DSPSV('U',3,1,moi,ipiv,angvel,3,info)
#ifdef RISM_DEBUG
  write(6,*) "ANGULAR VELOCITY", angvel
#endif /*RISM_DEBUG*/

  dir = angvel/sqrt(sum(angvel**2))
  
  !calculate moment of inertia w.r.t. the axis of rotation
  moi(1) = 0d0
  do iatu=1,natu
     r1 = ratu1(:,iatu)-cm1
     proj = sum(r1*dir)/sum(dir**2)*dir
#ifdef RISM_DEBUG
     write(6,*) "rxv",(r1-proj)*dir
#endif /*RISM_DEBUG*/
     moi(1) = moi(1) + sum((r1-proj)**2)*mass(iatu)
  end do
#ifdef RISM_DEBUG
  write(6,*) "ANGULAR MOMENTUM", moi(1)*angvel
#endif /*RISM_DEBUG*/
  calc_rotation = angvel
!!$  moi=0
!!$  erot=0
!!$  do iatu=1,natu
!!$     r0 = ratu0(:,iatu)-cm0
!!$     call cross(r0,v(:,iatu),rxv)
!!$     proj = sum(r0*angvel)/sum(angvel**2)*angvel
!!$     moi=moi+mass(iatu)*sum((r0-proj)**2)
!!$     erot = erot + .5*mass(iatu)*sum((r0-proj)**2)*sum((rxv/sum(r0**2))**2)
!!$  end do
!!$  write(6,*) moi
!!$  do id=1,3
!!$     write(6,*) id,sum(mass*v(id,:))
!!$     write(6,*) id,angvel(id),sum(angvel**2)
!!$  end do
!!$  write(6,*) "EROT", 0.5*moi*sum(angvel**2), erot
!!$  write(6,*) "EROT", erot

end function calc_rotation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate the principal axes of the atom distribution.  
!This follows, in part, mofi() in nmode/thermo.f 
!IN:
!   ratu   :: the x,y,z position of each solute atom. (3,natom)
!   mass   :: mass of each atom
!   pa     :: the three prinicpal axes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calc_pa(ratu,mass,pa)
  implicit none
  _REAL_,intent(in) :: ratu(:,:),mass(:)
  _REAL_, intent(out) :: pa(3,3)
  integer :: id,ier

  !t        : moment of inertia tensor in 1-d, in upper triangular, 
  !           column-major format (xx, xy, yy, xz, yz, zz)
  !eigenval : This will be the moment of interia in each direction
  !work     : temp space for the algorithm
  _REAL_ :: t(6),eigenval(3),eigenvec(3,3),work(3*3)
#ifdef RISM_DEBUG
  write(6,*)"CALC PA"
  call flush(6)
#endif /*RISM_DEBUG*/
  
  call dspev('V','U',3,calc_moi(ratu,mass),eigenval,pa,3,work,ier)

end subroutine calc_pa


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate the moment of inertia in upper triangular form.
!           column-major format (xx, xy, yy, xz, yz, zz)
!This follows, in part, mofi() in nmode/thermo.f 
!IN:
!   ratu   :: the x,y,z position of each solute atom. (3,natu)
!   mass   :: mass of each atom
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function calc_moi(ratu,mass)
  implicit none
  _REAL_,intent(in) :: ratu(:,:),mass(:)

  !calc_moi        : moment of inertia tensor in 1-d, in upper triangular, 
  _REAL_ :: calc_moi(6)
#ifdef RISM_DEBUG
  write(6,*)"CALC PA"
  call flush(6)
#endif /*RISM_DEBUG*/
  calc_moi(1) = sum(mass*(ratu(2,:)**2+ratu(3,:)**2))
  calc_moi(3) = sum(mass*(ratu(1,:)**2+ratu(3,:)**2))
  calc_moi(6) = sum(mass*(ratu(1,:)**2+ratu(2,:)**2))
  calc_moi(2) = -sum(mass*(ratu(1,:)*ratu(2,:)))
  calc_moi(4) = -sum(mass*(ratu(1,:)*ratu(3,:)))
  calc_moi(5) = -sum(mass*(ratu(2,:)*ratu(3,:)))

end function calc_moi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Assumes that the principal axes have been aligned to coincide with
!the x-, y- and z- axes for both structures. The strctures are then
!compared to eachother and the first rotated to best fit the second
!while maintaining the orientation of the PA

!IN:
!   ratu  :: the x,y,z position of each solute atom.
!   ratu2 :: the x,y,z position of each solute atom.
!   natu  :: the number of solute atoms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine alignorient(ratu,ratu2,natu,backquat)
  use constants, only : PI
  use quaternion, only : rotate_quat, quat_mult
  use rism_report_c
  implicit none
  integer,intent(in) :: natu
  _REAL_,intent(inout) :: ratu(3,natu)
  _REAL_,intent(in) :: ratu2(3,natu)
  _REAL_, intent(out) :: backquat(4)
  _REAL_ :: rotaxis(3), quat(4)
  integer :: axis(3)
  integer :: id,iatu

  !here we calculate the product of each coordinate of each atom in
  !the two structures.  If both have had their PA aligned to the
  !coordinate system there are zero to three 180 degree rotations that
  !will align the two molecules.  The rotation axes are those with
  !over all positive coordinate products (some atom movement may
  !distort the molecule).

#ifdef RISM_DEBUG
  write(6,*) "ALIGNORIENT"
#endif /* RISM_DEBUG*/
  do id = 1, 3
     rotaxis(id) = sum(ratu(id,:)*ratu2(id,:))
     if(rotaxis(id) <0) then
        axis(id) = 0
     else
        axis(id)=1
     end if
  end do
#ifdef RISM_DEBUG
  write(6,*) "rotate axis", axis, rotaxis
#endif /*RISM_DEBUG*/
  if(sum(axis)==0)then
     call rism_report_error("ALIGNORIENT failed")
  elseif(sum(axis)==3)then
     !already aligned
     backquat=(/1d0,0d0,0d0,0d0/)
     return
  end if
  backquat=0d0
  do id = 1, 3
     if(axis(id) == 1)then
        rotaxis =0d0
        rotaxis(id) = 1d0

!!$        quat(2:4) = rotaxis*sin(PI/2d0)
!!$        quat(1) = cos(PI/2d0)
        quat=0d0
        quat(id+1) = 1d0

        do iatu = 1 ,natu
           call rotate_quat(ratu(1:3,iatu),quat)
        end do
#ifdef RISM_DEBUG
        write(6,*) "QUAT", quat, backquat
#endif /*RISM_DEBUG*/
        if(sum(backquat) == 0d0 ) then
           backquat(1) = quat(1)
           backquat(2:4) = -1d0*quat(2:4)
!           backquat(2:4) = quat(2:4)
        else
           quat(2:4) = -1d0*quat(2:4)
           call quat_mult(quat,backquat,backquat)
        end if
#ifdef RISM_DEBUG
        write(6,*) "QUAT", quat, backquat
#endif /*RISM_DEBUG*/
     end if
  end do
#ifdef RISM_DEBUG
  write(6,*) ratu(:,1)
#endif /*RISM_DEBUG*/
end subroutine alignorient

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Rotates the system such that the previously calculated pricipal axes
!coincide with the x-,y- and z-axes.
!IN:
!   ratu   :: the x,y,z position of each solute atom.  (3,natom) This is modified.
!   pa     :: the three prinicpal axes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine orient_pa(ratu,pa,backquat)
  use constants, only : PI
  use quaternion, only : new_quat,rotate_quat,quat_mult
  implicit none
  _REAL_,intent(inout) :: ratu(:,:)
  _REAL_, intent(inout) :: pa(3,3)
  _REAL_, intent(out) :: backquat(4)
  _REAL_ :: angle, quat(4),dir(3),xaxis(3)=(/1d0,0d0,0d0/),yaxis(3)=(/0d0,1d0,0d0/),&
       checkv(3)
  integer :: natu

  integer :: iatu,ipa
#ifdef RISM_DEBUG
  write(6,*) "ORIENT_PA"
#endif /*RISM_DEBUG*/
  natu=ubound(ratu,2)

  !
  !rotate first principal axis the x-axis
  !

  !get the angle.  This is always positive
  angle = acos(min(1d0,max(-1d0,dot_product(pa(1:3,1),xaxis))))
  if(angle < PI - 1d-6 .and. angle > -PI+1d-6)then
     call cross(pa(1:3,1),xaxis,dir)
  else
     dir = yaxis
  end if
  
  !get the cross product between pa(:,1) and dir.  This will determine
  !if we should rotate using +ive or -ive angle
  call cross(dir,pa(1:3,1),checkv)
  if(dot_product(checkv,xaxis) < 0)then
     angle = -angle
  end if
  quat = new_quat(angle,dir)
  do iatu = 1 ,natu
     call rotate_quat(ratu(1:3,iatu),quat)
  end do
  do ipa = 1,3
     call rotate_quat(pa(1:3,ipa),quat)
  end do
  backquat = quat

  !Next the second PA (this also places the third one as well)
  !max/min protects against round-off errors in normalized vectors 
  angle = acos(min(1d0,max(-1d0,dot_product(pa(1:3,2),yaxis))))
  dir = xaxis
  call cross(dir,pa(1:3,2),checkv)
  
  if(dot_product(checkv,yaxis) < 0)then
     angle = -angle
  end if

  quat = new_quat(angle,dir)
  do iatu = 1 ,natu
     call rotate_quat(ratu(1:3,iatu),quat)
  end do
  do ipa = 1,3
     call rotate_quat(pa(1:3,ipa),quat)
  end do

  call quat_mult(quat,backquat,backquat)
  backquat(2:4) = -1d0*backquat(2:4)

#ifdef RISM_DEBUG
  write(6,*) ratu(:,1)
#endif /*RISM_DEBUG*/

end subroutine orient_pa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Tests if number is prime. Note: this is a slow algorithm.  It checks if 2 or
!!!any odd number (except 1) less than the square root of number is a factor.
!!!IN:
!!!   number : number to test
!!!OUT:
!!!    .true. if prime, otherwise .false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function isprime(number)
  implicit none
  logical :: isprime
  integer :: i,number,j
  isprime = .false.
  i = abs(number)
  if(mod(i,2)==0) return 
  j=3
  do while (j**2 < i)
     if(mod(i,j)==0) return 
     j = j+2
  end do
  isprime = .true.
end function isprime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Checks if the the number is factorizable by the list of numbers given. 
!!!Works best for prime numbers.  It is also best to pre-sort the list from 
!!!highest to lowest.
!!!IN:
!!!   number : number to test
!!!   factors: potential factors
!!!OUT:
!!!    .true. if factorizable, otherwise .false.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function isfactorable(number,factor)
  implicit none
  integer, intent(in) :: number, factor(:)
  logical :: isfactorable
  integer :: i, num, numold
  isfactorable=.false.
  num = number
  do i=1, ubound(factor,1)
     numold=num+1
     do while (numold>num)
        numold=num
        if(mod(num,factor(i)) == 0)then
           num=num/factor(i)
        end if
     end do
  end do
  if(num == 1) isfactorable=.true.
end function isfactorable


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Fills PTR with indices of VAL to give the values of VAL in accending order.
!!!I.e.  VAL(PTR) will be in accending order.  The reproduces the functionality
!!!of INDEXX from Numerical Recipes.
!!!IN:
!!!   val :: array of values
!!!   ptr :: will be an array of VAL indices.  I.e. it is a 'pointer' to the VAL
!!!          elements
!!!   n   :: length of the arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine index_array(val,ptr,n)
  implicit none
  _REAL_, intent(in) :: val(n)
  integer, intent(out) :: ptr(n)
  integer, intent(in) :: n
  integer :: i, temp

  !Following Numerical Recipes, we are using heapsort to sort the array in place.
  !For the actual heapsort algorithm we have follow Wikipedia heapsort article
  !from 2010/02/13.

  !initialize ptr
  ptr =  (/(i, i=1,n)/)

  !build the heap
  do i =  (n+1)/2, 1, -1
     call siftdown(val,ptr,i,n,n)
  end do

  !Traverse the heap to produce the sorted order
  do i = n, 1, -1
     !move the largest element (found in the first index) to the last index
     temp = ptr(i)
     ptr(i) = ptr(1)
     ptr(1) = temp
     !bring up the next largest element to the first index
     call siftdown(val(1:i-1),ptr(1:i-1),1,i-1,n)
  end do
end subroutine index_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Used in the index_array heapsort.  Brings the largest value between I and M
!!!to index I while ensuring that values at 2*I and 2*I+1 are less than the value
!!!at I.
!!!IN:
!!!   val :: array of values
!!!   ptr :: will be an array of VAL indices.  I.e. it is a 'pointer' to the VAL
!!!          elements
!!!   i   :: the start element to use for the arrays
!!!   m   :: the end element to use for the arrays
!!!   n   :: length of the arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine siftdown(val,ptr,i,n,m)
  implicit none
  _REAL_, intent(in) :: val(m)
  integer, intent(inout) :: ptr(m)
  integer, intent(in) :: i,n,m
  integer:: root,child
  integer :: temp
  root=i
  do while (root*2 <= n)
     child = root*2
     if(child+1 <= n) then
        if(val(ptr(child+1))> val(ptr(child))) child = child+1
     end if
     if(child <= n) then
        if(val(ptr(child))> val(ptr(root)))then
           temp = ptr(child) 
           ptr(child) = ptr(root)
           ptr(root) = temp
           root = child
        else 
           exit
        end if
     end if
  end do
end subroutine siftdown

subroutine linLeastSqFit(basis, xa, ya, coef, o_sig, o_chisq, o_varcoef)
  use safemem
  implicit none
  interface 
     subroutine basis(x,xn,n)
       implicit none
       _REAL_, intent(in) :: x
       _REAL_, intent(out) :: xn(n)
       integer, intent(in) ::n
     end subroutine basis
  end interface
  _REAL_,intent(in) :: xa(:), ya(:)
  _REAL_,intent(out) :: coef(:)
  _REAL_,optional,intent(in) :: o_sig(:)
  _REAL_,optional,intent(out) :: o_chisq, o_varcoef(:)
  
  _REAL_ :: sig(size(ya)), A(size(xa),size(coef))
  !B : working vector for DGELSD.  On input it contains ya(), on
  !    output it contains coef
  !SV : singular values of A from DGELSD
  _REAL_ :: B(size(ya),1), SV(size(coef))
  !work : work space for DGELSD
  _REAL_, pointer :: work(:)=>NULL()
  !info : return value of DGELSD
  !rank : effective rank of A determined by DGELSD
  integer :: info, rank
  integer, pointer :: iwork(:)=>NULL()
  integer :: i, j, nrow,ncol
  sig = 1d0
  if(present(o_sig)) sig = o_sig

  !prepare input array
  ncol=size(coef)
  nrow=size(xa)
!  write(0,*) "nrow, ncol", nrow, ncol
!  A(:,1) = 1d0
!  A(:,1) = xa**2
!  do j = 2, ncol
!     do i = 1, nrow
!        A(i,j) = A(i,j-1)*xa(i)
!     end do
!  end do
  do i = 1, nrow
     call basis(xa(i),A(i,:),ncol)
  end do
  do i = 1,nrow
!     write(0,*) A(i,:)
  end do
  B(:,1) = ya
  !get work memory requirement
!  write(0,*) "B",B
  iwork=>safemem_realloc(iwork,1)
  work=>safemem_realloc(work,1)
  call DGELSD(nrow,ncol,1,A,nrow,B,nrow,SV,-1d0,rank,work,-1,iwork,info)
  work=>safemem_realloc(work,int(work(1)))
  iwork=>safemem_realloc(iwork,iwork(1))
  call DGELSD(nrow,ncol,1,A,nrow,B,nrow,SV,-1d0,rank,work,size(work),iwork,info)
  if(info/=0)&
       call rism_report_error("(a,i4)","POLYFIT: DGELSD failed: ",info)
  if(safemem_dealloc(work) /=0) &
       call rism_report_error("POLYFIT: Failed to deallocate workspace")
  if(safemem_dealloc(iwork) /=0) &
       call rism_report_error("POLYFIT: Failed to deallocate workspace")
  coef = B(1:ncol,1)

  !get chisq
  if(present(o_chisq))then
     o_chisq = 0d0
     !recompute A
!!$      A(:,1) = xa*2
!!$     do j = 2, ncol
!!$        do i = 1, nrow
!!$           A(i,j) = A(i,j-1)*xa(i)
!!$        end do
!!$     end do
     do i = 1, nrow
        call basis(xa(i),A(i,:),ncol)
     end do
     !reuse B to get ya - A*coef
     B(:,1) = ya
     call DGEMV('N',nrow,ncol,-1d0,A,nrow,coef,1,1d0,B,1)
     o_chisq = sum(B(:,1)**2/sig)
  end if
end subroutine linLeastSqFit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!'Progressive' polynomial interpolation using Neville's algorithm. Adds
!!!one point at a time, starting from the lowest indicies, to the
!!!interpolation until either the data is exhausted or the relative
!!!difference begins to increase.  Can be used as a drop in
!!!replacement for poly_interp.
!!!IN:
!!!   xa :: array of x values
!!!   ya :: array of y values
!!!   n  :: number of array elements, must be >2
!!!   x  :: argument to the polynomial
!!!   y  :: (out) polynomial value at x
!!!   dy :: (out) error estimate in y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine poly_interp_progressive(xa, ya,n, x, y, dy)
  implicit none
  _REAL_, intent(in) :: xa(n),ya(n),x
  _REAL_, intent(out) :: y, dy
  integer, intent(in) :: n
  _REAL_ :: P(n),rel_diff,rel_diff0, y0,dy0
  integer :: i,m
  rel_diff=huge(1d0)
  rel_diff0=huge(1d0)
  y0=huge(1d0)
  P=ya
  do m = 1,n-2
     !Neville's algorithm
     do i = 1, n-m
        p(i) = (x-xa(i+m))*p(i) - (x-xa(i))*p(i+1)
        p(i) = p(i)/(xa(i)-xa(i+m))
     end do
     y=((x-xa(n))*p(1) - (x-xa(1))*p(2))/(xa(1)-xa(n))

     !compute relative difference between this and the previous iteration
     if(y0 /= huge(1d0))then
        rel_diff = abs((y0-y)/y)
     end if
     !estimated error
     dy = (p(1)+p(2)-2d0*y)/2d0
     !if the relative difference is not decreasing, we're done
     if(rel_diff0 < rel_diff)then
        y=y0
        return
     end if
     !update old values
     y0=y
     dy0=dy
     rel_diff0=rel_diff
  end do
  !if we get here, we are just returning the y from using all points
end subroutine poly_interp_progressive

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Polynomial interpolation using Neville's algorithm.
!!!Serves as a drop in replacement for Numerical Recipes POLINT subroutine; 
!!!however, the error estimate is done differently but is of the same order.
!!!IN:
!!!   xa :: array of x values
!!!   ya :: array of y values
!!!   n  :: number of array elements
!!!   x  :: argument to the polynomial
!!!   y  :: polynomial value at x
!!!   dy :: error estimate in y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine poly_interp(xa, ya,n, x, y, dy)
  implicit none
  _REAL_, intent(in) :: xa(n),ya(n),x
  _REAL_, intent(out) :: y, dy
  integer, intent(in) :: n
  _REAL_ :: P(n)
  integer :: i,m
  P=ya
  do m = 1,n-2
     do i = 1, n-m
        p(i) = (x-xa(i+m))*p(i) - (x-xa(i))*p(i+1)
        p(i) = p(i)/(xa(i)-xa(i+m))
     end do
  end do
  y=((x-xa(n))*p(1) - (x-xa(1))*p(2))/(xa(1)-xa(n))
  dy = (p(1)+p(2)-2d0*y)/2d0
end subroutine poly_interp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Computes support abscissas (nodes) and weights for Gaussian quadrature with
!!!Legendre polynomials.  This is a drop in replacement for Numerical Recipes
!!!GAULEG.
!!!IN:
!!!   a :: lower integration bound
!!!   b :: upper integration bound
!!!   x :: abscissas
!!!   weight :: weights
!!!   n :: number of points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gaussquad_legendre (a,b,x,weight,n)
  use constants, only : PI
  implicit none
  _REAL_,intent(in) :: a, b
  _REAL_,intent(out) :: x(n), weight(n)
  integer, intent(in) :: n
  
  !nroot :: roots are symmetric about zero we only need to find one half
  integer :: nroot

  integer :: iroot
  
  _REAL_ :: root0, root,p,dp
  !eps :: convergence criterion for roots
  _REAL_ :: eps=1d-14

  nroot = (n+1)/2
  do iroot = 1, nroot
     !compute roots
     !initial guess
     root0 = cos(PI*(iroot-0.25d0)/(n+0.5d0))
     do while(.true.) !Newton-Rhapson.  Very well behaved so we risk the infinite loop
        call legendre(root0,p,dp,n)
        root = root0-p/dp
        if(abs(root-root0) < eps) exit
        root0=root
     end do
     !scale roots
     x(iroot) = a+(1d0-root)*(b-a)/2d0
     x(n-iroot+1) = b-(1d0-root)*(b-a)/2d0
     !compute scaled integral
     weight(iroot) = (b-a)/ ((1d0-root*root)*dp*dp)
     weight(n-iroot+1) = weight(iroot)
  end do
end subroutine gaussquad_legendre

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the value of the Nth order Lengendre polynomial at X and its first
!!!derivative.
!!!IN:
!!!   x :: -1 < x < 1
!!!   y :: P_n(x)
!!!   dy :: P'_n(x)
!!!   n :: order of the polynomial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine legendre(x,y,dy,n)
  implicit none
  _REAL_, intent(in) :: x
  _REAL_, intent(out) :: y,dy
  integer, intent(in) :: n
  _REAL_ :: y0,y1
  integer :: i
  y0 = 1d0
  y = 1d0
  dy = 0d0
  if(n==0) return
  y1 = x
  y = x
  dy = 1d0
  if(n==1) return
  do i = 2,n
     y = ((2d0*i-1d0)*x*y1 -(i-1d0)*y0)/(i)
     dy = (x*y-y1)*i/(x*x-1d0)
     y0=y1
     y1=y
  end do
end subroutine legendre

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Computes support abscissas (nodes) and weights for Gaussian quadrature with
!!!Laguerre polynomials.  The domain of integration is always [0,infinity).  For
!!!double precission n=180 seems to be the largest order that is reasonable. 
!!!Stroud and Secrest suggest n=42 as the largest to avoid floating point overflow.
!!!This may be due single precision.
!!!
!!!This also differs from traditional Gauss-Laguerre in that weight(i) include
!!!a factor of exp(x(i)) so this may be applied much like Gauss-Legendre.  However,
!!!to achieve suitable results the integrated function should converge faster exp(-x).
!!!IN:
!!!   x :: abscissas
!!!   weight :: weights
!!!   n :: number of points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gaussquad_laguerre (x,weight,n)
  use constants, only : PI
  implicit none
  _REAL_,intent(out) :: x(n), weight(n)
  integer, intent(in) :: n

  !alpha : used for the associate Laguerre polynomial.  0 gives the standard polynomial
  integer,parameter ::  alpha = 0

  !Root guesses
  !root0   : current guess
  !root00  : guess from last iteration
  !root000 : guess from two iterations ago
  !root    : solved root
  !p       : value of polynomial
  !dp      : derivative of polynomial
  _REAL_ :: root0=0, root00=0, root000=0, root,p,dp

  !eps :: convergence criterion for roots
  _REAL_ :: eps=1d-12

  integer :: iroot, icount

  !Root guesses are from A. H. Stroud and D. Secrest. Gaussian Quadrature Formulas. Prentice-Hall. (1966)
  !these are slightly modified in that we use the previous real roots to seed our next guess and not
  !the previous guesses.
  do iroot = 1, n
     !compute roots
     !initial guess
     if(iroot==1)then
        root0=(1d0+alpha)*(3d0+0.92d0*alpha)/(1d0+2.4d0*n+1.8d0*alpha)
     elseif(iroot==2)then
        root0 = root00 + (15d0+6.25d0*alpha)/(1d0+0.9*alpha + 2.5d0*n)
     else
        root0 = root00 +(root00-root000)/(1d0+0.3*alpha) * ((1+2.55d0*iroot)/(1.9d0*iroot)&
             +1.26d0*iroot*alpha/(1d0+3.5d0*iroot))
     end if
     icount = 0
     do while(.true.) !Newton-Rhapson.
        call laguerre(root0,p,dp,n)
        root = root0-p/dp
        if(abs(root-root0)/root < eps) exit
        icount = icount+1
        if(icount >1000)then
           write(0,'(a,i4,e24.16,e24.16)') "ERROR: gaussquad_laguerre: Newton-Rhapson failed.  Try a lower order.",&
                iroot, root0, root
           stop
        end if
        root0=root
     end do
     root000 = root00
     root00 = root0
     root0 = root
     x(iroot) = root
     call laguerre(root,p,dp,n+1)
     weight(iroot) = x(iroot)/(dble(n+1)*p)**2*exp(x(iroot))
     if(iroot > 1 .and. x(iroot) <= x(iroot-1))then
        write(0,'(a,i4)') "ERROR: gaussquad_laguerre failed to find unique roots for order", n
        stop
     end if
  end do
end subroutine gaussquad_laguerre

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the value of the Nth order Laguerre polynomial at X and its first
!!!derivative.
!!!IN:
!!!   x :: 0 < x < infinity
!!!   y :: P_n(x)
!!!   dy :: P'_n(x)
!!!   n :: order of the polynomial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine laguerre(x,y,dy,n)
  implicit none
  _REAL_, intent(in) :: x
  _REAL_, intent(out) :: y,dy
  integer, intent(in) :: n
  _REAL_ :: y0,y1
  integer :: i
  y0 = 1d0
  y = 1d0
  dy = 0d0
  if(n==0) return
  y1 = -x+1d0
  y = -x+1d0
  dy = -1d0
  if(n==1) return
  do i = 2,n
     y = ((2d0*i-1d0 -x)*y1-dble(i-1)*y0)/dble(i)
     dy = i*(y-y1)/x
     y0=y1
     y1=y
  end do
end subroutine laguerre

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Computes zeroth order spherical Bessel function.  Adapted from the GNU 
!!!Scientific library.
!!!IN:
!!!   x :: calculate the spherical Bessel function a x
!!!   o_err :: (optional) absolute error in the result.  For abs(x)>0.5d0 this is
!!!            an approximation of the GSL version since we don't know the error
!!!            for the intrinsic sin(x)
!!!OUT:
!!!    the value of the zeroth spherical Bessel function at x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function spherical_bessel_j0(x,o_err) result(val)
  implicit none
  _REAL_, intent(in) :: x
  _REAL_, optional, intent(out) :: o_err
  _REAL_ :: val
  _REAL_ :: ax
  _REAL_, parameter :: c1 = -1d0/6d0, &
       c2 =  1d0/120d0,&
       c3 = -1d0/5040d0,&
       c4 =  1d0/362880d0,&
       c5 = -1d0/39916800d0,&
       c6 =  1d0/6227020800d0
  ax = abs(x)
  if(ax < 0.5d0) then
     ax = ax*ax
     val = 1d0 + ax*(c1 + ax*(c2 + ax*(c3 + ax*(c4 + ax*(c5 + ax*c6)))))
     if(present(o_err)) o_err = 2d0 * epsilon(1d0) * abs(val)
  else
     val  = sin(x)/x;
     if(present(o_err))then
        !we don't know what the error is for intrinsic sin(x)!
        o_err  = abs(epsilon(1d0)/x)  + 2.0 * epsilon(1d0) * abs(val)
     end if
  end if
end function spherical_bessel_j0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Computes first order spherical Bessel function.  Adapted from the GNU 
!!!Scientific library.
!!!IN:
!!!   x :: calculate the spherical Bessel function a x
!!!   o_err :: (optional) absolute error in the result.  For abs(x)>0.25d0 this is
!!!            an approximation of the GSL version since we don't know the error
!!!            for the intrinsics sin(x) and cos(x)
!!!OUT:
!!!    the value of the first spherical Bessel function at x
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function spherical_bessel_j1(x,o_err) result(val)
  implicit none
  _REAL_, intent(in) :: x
  _REAL_, optional, intent(out) :: o_err
  _REAL_ :: val
  _REAL_ :: ax
  _REAL_, parameter :: c1 = -1d0/10d0,&
       c2 =  1d0/280d0,&
       c3 = -1d0/15120d0,&
       c4 =  1d0/1330560d0,&
       c5 = -1d0/172972800d0

  ax = abs(x)
  if(x==0)then
     val=0d0
     if(present(o_err)) o_err=0d0
  elseif(ax < 3.1d0*tiny(1d0)) then
     !GSL  gives an underflow error here.  We'll just approximate it as zero
     val=0d0
     if(present(o_err)) o_err=0d0
  elseif(ax < 0.25d0) then
     ax = ax*ax
     val = x/3d0*(1d0 + ax*(c1 + ax*(c2 + ax*(c3 + ax*(c4 + ax*c5)))))
     if(present(o_err)) o_err = 2d0 * epsilon(1d0) * abs(val)
  else
     val  = (sin(x)/x - cos(x))/x;
     if(present(o_err))then
        !we don't know what the error is for intrinsic sin(x)!
        o_err  = (abs(epsilon(1d0)/x)+epsilon(1d0))/abs(x)  + 2.0 * epsilon(1d0) * abs(val)
        o_err = o_err+ 2.0 * epsilon(1d0)*(abs(sin(x)/(x*x)) + abs(cos(x)/x))
     end if
  end if
end function spherical_bessel_j1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Does a MPI sum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function checksum(a,n,comm)
  implicit none
#ifdef MPI
  include 'mpif.h'
#endif /*MPI*/
  _REAL_, intent(in) :: a(n)
  integer, intent(in) :: n,comm
  _REAL_ :: checksum, temp
  integer :: err
  checksum = sum(a)
#ifdef MPI
#ifdef USE_MPI_IN_PLACE
    call mpi_allreduce(MPI_IN_PLACE,checksum,1,MPI_DOUBLE_PRECISION,mpi_sum,comm,err)
#else
    call mpi_allreduce(checksum,temp,1,MPI_DOUBLE_PRECISION,mpi_sum,comm,err)
    checksum = temp
#endif /*USE_MPI_IN_PLACE*/
    if(err /=0) call rism_report_error&
         ("RISM3D CHECKSUM: could not reduce CHECKSUM")
#endif /*MPI*/
  
end function checksum


function lcm(a,b)
  implicit none
  integer :: lcm,a,b
  lcm = a*b/gcd(a,b)
end function lcm

recursive function gcd(s,t) result(res)
  implicit none
  integer :: a,b,c,s,t,res
  a=s
  b=t
!  write(6,*)"GCD",a,b
  if(b>a)then
    c=a
    a=b
    b=c
  end if
  if(b==0) then
    res = a
    return
  end if
  a=mod(a,b)
  res= gcd(b,a)
!  write(6,*)"GCD",gcd
end function gcd
!!endfix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Convert string to upper case
!!!IN:
!!!   cc : string will be converted to upper case in place
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  caseup (cc)
  implicit none
  character(len=*)  cc
  character(len=1) ::  bc
  integer*1  bi
  equivalence (bc,bi)
  integer ::  i
  do i=1,len(cc)
     bc = cc(i:i)
     if (bi >= ichar('a') .AND. bi <= ichar('z')) &
          cc(i:i) = char(ibclr(bi,5))
  enddo
  return
end subroutine caseup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Convert string to lower case
!!!IN:
!!!   cc : string will be converted to lower case in place
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  caselow (cc)
  implicit none
  character(len=*)  cc
  character(len=1) ::  bc
  integer*1  bi
  equivalence (bc,bi)
  integer ::  i
  do i=1,len(cc)
     bc = cc(i:i)
     if (bi >= ichar('A') .AND. bi <= ichar('Z')) &
          cc(i:i) = char(ibset(bi,5))
  enddo
  return
end subroutine caselow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Find a free fortran unit for temporary use.
!!!IN:
!!!   start : (optional) minimum unit number that is acceptible 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function freeUnit(start) result(unit)
  implicit none
  integer, optional, intent(in) :: start
  integer :: unit
  logical :: opened
  opened = .true.
  !skip commonly used unit numbers
  unit = 10
  !user defined minimum unit number
  if(present(start)) unit = start
  !search for free unit
  unit = unit-1
  do while(opened)
     unit = unit+1
     inquire(unit,opened=opened)
  end do
end function freeUnit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Intel CPUs typically use extended (80-bit) precision when
!!!possible.  However, different compilers may or may not use the
!!!extended precision value when printing.  To ensure consistency,
!!!we test each value and values less than TINY() are printed as
!!!zero.  These very small values are actually meaningless and
!!!break our testing proceedures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
elemental function rmExPrec(ep) result(rp)
  implicit none
  _REAL_, intent(in) :: ep
  _REAL_ :: rp
  if(abs(ep) < tiny(ep))then
     rp=0d0
  else
     rp=ep
  end if
end function rmExPrec
end module rism_util
