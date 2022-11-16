#include "../include/dprec.fh"
      module interactmask
        _REAL_,  allocatable :: rk0   (:) ! Initial BOND_FORCE_CONSTANT
        _REAL_,  allocatable :: tk0   (:) ! Initial ANGLE_FORCE_CONSTANT
        _REAL_,  allocatable :: pk0   (:) ! Initial DIHEDRAL_FORCE_CONSTANT
        _REAL_,  allocatable :: cn1_0 (:) ! Initial LENNARD_JONES_ACOEF
        _REAL_,  allocatable :: cn2_0 (:) ! Initial LENNARD_JONES_BCOEF

        integer mskerr,mski
      end module

! ***********************************************************************
! ***********************************************************************
!

! Modify the RK Bond Force Constant depending of the mask vector
! n   Number of atoms on the mask
! mask  Integer vector of dimension n with 0 or 1.
!       mask(i) = 0 -> All the bonds with atom number i involved won't be
!                      considered
!       mask(i) = 1 -> All the bonds with atom number i involved will be
!                      considered
!
subroutine bondmask(nbin,ib,jb,icb,mask)
 
   use parms, only:rk
   use interactmask
   
   integer nbin
   integer ib(*),jb(*),icb(*),mask(*)
    
    integer i,ata,atb,n

    do i = 1,nbin
      ata   = ib(i)
      atb   = jb(i)
      n     = icb(i)
      rk(n) = rk0(n)*mask(ata)*mask(atb)
    enddo
    
    return
end subroutine bondmask

! ***********************************************************************
! ***********************************************************************
!

subroutine anglmask(nbain,it,jt,kt,ict,mask)

   use parms, only:tk
   use interactmask

   integer nbain
   integer it(*),jt(*),kt(*),ict(*),mask(*)
    
   integer i,ata,atb,atc,n

   do i = 1,nbin
     ata   = it(i)
     atb   = jt(i)
     atc   = kt(i)
     n     = ict(i)
     tk(n) = tk0(n)*mask(ata)*mask(atb)*mask(atc)
   enddo
    
   return
end subroutine anglmask

! ***********************************************************************
! ***********************************************************************
!

subroutine ephimask(nphiin,ip,jp,kp,lp,icp,mask)

   use parms, only:pk
   use interactmask
   
   integer nphiin
   integer ip(*),jp(*),kp(*),lp(*),icp(*),mask(*)
    
   integer i,ata,atb,atc,atd,n

   do i = 1,nphiin
     ata   = ip(i)
     atb   = jp(i)
     atc   = kp(i)
     atd   = lp(i)
     n     = icp(i)
     pk(n) = pk0(n)*mask(ata)*mask(atb)*mask(atc)*mask(atd)
   enddo
    
   return
end subroutine ephimask
