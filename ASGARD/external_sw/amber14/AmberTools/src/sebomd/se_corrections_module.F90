module se_corrections_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!     Author   :: Antoine MARION 
!!     Date     :: 2013-12-10
!!     Function :: Perform corrections to SE Hamiltonians
!!                 - PM3-PIF [1,2,3]
!!                 - PM3-MAIS [1,2]
!!            
!!                 - MM torsion correction for amide 
!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 use sebomd_arrays, only : put_Z, &
                           put_Zeff, &
                           get_interaction
 use se_corrections_tools, only : LookFor_PEP, GetRij, &
                                  putConstants, &
                                  initPIF, computePIF, &
                                  initMAIS, computeMAIS, &
                                  computeGauss, compPEP
                            
 implicit none
 private
 
  
 integer :: natoms, Npep
 !!
 logical :: if_pbc
 logical :: doPIF, doMAIS, doPEP


 public :: SEcorr_init, Apply_se_corrections

 contains
  subroutine SEcorr_init(n,ctype,mm_pep,k_pep,let_pm3,zatms,zchg,symbol2,agauss,bgauss,cgauss, &
                          pbc,xyz,cst1,cst2,cst3,ierror)
   !!! Initialize options and data for SE correction
   !!! ctype gives the type of correction to perform
   !!! ctype = 0 :: No correction
   !!!         1 :: PIF2
   !!!         2 :: PIF3
   !!!         3 :: MAIS1
   !!!         4 :: MAIS2
   !!! PIF1 is automatically included in PIF2 and PIF3 
   !!! for water-water interactions

   implicit none
   integer, intent(in) :: n
   integer, intent(in), dimension(:) :: zatms
   double precision, intent(in), dimension(0:) :: zchg
   double precision, intent(in) :: xyz(3,n)
   double precision, dimension (:,:), intent(in) :: agauss, bgauss, cgauss
   logical, intent(in) :: mm_pep, pbc, let_pm3
   double precision, intent(in) :: cst1,cst2,cst3
   character(len=4), intent(in), dimension(:) :: symbol2
   integer, intent(in) :: ctype
   double precision, intent(in) :: k_pep
   integer, intent(inout) :: ierror
   ! local
   integer :: i

   natoms = n
   if_pbc = pbc

   call putConstants(cst1,cst2,cst3)

   doPIF  = .false.
   doMAIS = .false.
   doPEP  = .false.

   !! Store atomic numbers and core charges in sebomd_arrays
   do i=1,natoms
    call put_Z(i,zatms(i))
    call put_Zeff(i,zchg(zatms(i)))
   enddo

   if (mm_pep) then
     call LookFor_PEP(natoms,xyz,k_pep,if_pbc,symbol2,Npep)
     write(6,'(" ----------------------------------------------------------- ")')
     write(6,'("                   peptide correction ")')
     write(6,'("")')
     if (Npep.ne.0) then
      doPEP = .true.
      write(6,'(" Number of amide group(s) found for MM corection:",i10)') Npep/2
     else
      write(6,'(" NO amide group found for MM corection")')
     endif
     write(6,'("")')
     write(6,'(" ----------------------------------------------------------- ")')
   endif

   if ((ctype.ne.0)) then
     if ((ctype.eq.1).or.(ctype.eq.2)) then ! Means PIF1,2 or 3
      doPIF = .true.
      call initPIF(natoms,ctype,let_pm3,symbol2,agauss,bgauss,cgauss,ierror)
     elseif ((ctype.eq.3).or.(ctype.eq.4)) then ! Means MAIS1 or MAIS2
      doMAIS = .true.
      call initMAIS(natoms,ctype,let_pm3,agauss,bgauss,cgauss,symbol2,ierror)
     else
      write(6,'(a,i5)') "Problem with SE correction option, ctype = ",ctype
     endif
   endif
   
   return

  end subroutine SEcorr_init

  subroutine Apply_se_corrections(xyz,e_se,e_pep,grad,vir_se)
   !!! Performs SE corrections
   implicit none
   double precision, intent(in) :: xyz(3,natoms)
   double precision, intent(out) :: grad(3,natoms)
   double precision, intent(out) :: vir_se(3)
   double precision, intent(out) :: e_se, e_pep 

   e_se = 0.d0
   e_pep = 0.d0
   grad(:,:) = 0.d0 
   vir_se(:) = 0.d0
   
   if (doPIF) then
    call ApplyPIF(xyz,e_se,grad,vir_se)
   elseif (doMAIS) then
    call ApplyMAIS(xyz,e_se,grad,vir_se)
   endif

   if (doPEP) call ApplyPEP(xyz,e_pep,grad)

   return

   end subroutine Apply_se_corrections

   subroutine ApplyPIF(xyz,e_se,grad,vir_se)
    !!! Compute PIF correction for all the intermolecular
    !!! interactions.
    !!! Returns the corrections to energy (e_se), gradient (grad)
    !!! and virial (vir_se)
    implicit none
    double precision, intent(inout) :: e_se
    double precision, intent(inout) :: grad(3,natoms)
    double precision, intent(inout) :: vir_se(3)
    double precision, intent(in) :: xyz(3,natoms)
    ! local
    integer :: i,j,k
    integer :: int_type
    double precision :: rij
    double precision :: uij(3)
    double precision :: e_gauss, de_gauss
    double precision :: e_pif, de_pif
    double precision :: de_se
#ifdef MPI
#include "mpif.h"
    integer :: nproc, myid, commsebomd
#endif

#ifdef MPI
    call se_mpi_vars(nproc, myid, commsebomd)
    do i = myid+2,natoms,nproc
#else
    do i = 2,natoms
#endif
     do j = 1,i-1
      call get_interaction(i,j,int_type)
      
      if (int_type.ne.0) then ! Modify the interaction

       call GetRij(i,j,xyz,if_pbc,rij,uij)
       call computeGauss(i,j,rij,e_gauss,de_gauss)
       call computePIF(i,j,rij,int_type,e_pif,de_pif)

       e_se = e_se - e_gauss + e_pif
       de_se = de_pif - de_gauss

       do k = 1,3
        grad(k,i) = grad(k,i) + de_se * uij(k)
        grad(k,j) = grad(k,j) - de_se * uij(k)
        vir_se(k) = vir_se(k) + de_se*rij*uij(k)**2
       enddo

      endif

     enddo
    enddo

#ifdef MPI
    call se_mpi_allreduce_dble(e_se,1)
    call se_mpi_allreduce_dble(grad,natoms*3)
    call se_mpi_allreduce_dble(vir_se,3)
#endif

    return

   end subroutine ApplyPIF

   subroutine ApplyMAIS(xyz,e_se,grad,vir_se)
    !!! Compute MAIS correction for both inter and intra
    !!! molecular interactions.
    !!! Returns the corrections to energy (e_se), gradient (grad)
    !!! and virial (vir_se)
    implicit none
    double precision, intent(inout) :: e_se
    double precision, intent(inout) :: grad(3,natoms)
    double precision, intent(inout) :: vir_se(3)
    double precision, intent(in) :: xyz(3,natoms)
    ! local
    integer :: i,j,k
    integer :: int_type
    double precision :: rij
    double precision :: uij(3)
    double precision :: e_gauss, de_gauss
    double precision :: e_mais, de_mais
    double precision :: de_se
#ifdef MPI
#include "mpif.h"
    integer :: nproc, myid, commsebomd
#endif

#ifdef MPI
    call se_mpi_vars(nproc, myid, commsebomd)
    do i = myid+2,natoms,nproc
#else
    do i = 2,natoms
#endif
     do j = 1,i-1
      call get_interaction(i,j,int_type)
      
      if (int_type.ne.0) then ! Modify the interaction

       call GetRij(i,j,xyz,if_pbc,rij,uij)
       call computeGauss(i,j,rij,e_gauss,de_gauss)
       call computeMAIS(i,j,rij,int_type,e_mais,de_mais)

       e_se = e_se - e_gauss + e_mais
       de_se = de_mais - de_gauss

       do k = 1,3
        grad(k,i) = grad(k,i) + de_se * uij(k)
        grad(k,j) = grad(k,j) - de_se * uij(k)
        vir_se(k) = vir_se(k) + de_se*rij*uij(k)**2
       enddo

      endif

     enddo
    enddo

#ifdef MPI
    call se_mpi_allreduce_dble(e_se,1)
    call se_mpi_allreduce_dble(grad,natoms*3)
    call se_mpi_allreduce_dble(vir_se,3)
#endif

    return

   end subroutine ApplyMAIS


   subroutine ApplyPEP(xyz,e_pep,grad)
    !!! Compute the MM peptidic correction for the 
    !!! dihedral identified at the init step
    !!! Returns the corrections to energy (e_pep) and gradient (grad)
    implicit none
    double precision :: e_pep
    double precision :: grad(3,natoms)
    double precision :: xyz(3,natoms)
    ! local
    integer :: ipep

    e_pep = 0.d0

    do ipep = 1,Npep
     call compPEP(ipep,if_pbc,xyz,e_pep,grad)
    enddo

    ! e_pep in eV
    ! grad in kcal/mol/Ang

    return

  end subroutine ApplyPEP

end module se_corrections_module
