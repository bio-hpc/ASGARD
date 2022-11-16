module sebomd_arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!     Author   :: Antoine MARION 
!!     Date     :: 2013-12-10
!!     Function :: Stores arrays needed to perform
!!                 SemiEmpirical corrections in sebomd
!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 implicit none
 private

 integer, allocatable, dimension(:) :: Zatms
 integer, allocatable, dimension(:,:) :: PEPlist
 integer, allocatable, dimension(:,:) :: interactions
 double precision, allocatable, dimension(:) :: Zatms_eff
 
 public :: init_sebomd_arrays
 public :: cleanup_sebomd_arrays
 public :: put_Z, get_Z
 public :: put_Zeff, get_Zeff
 public :: put_PEP, get_PEP
 public :: put_all_interactions,put_interaction, get_interaction
 
 contains
  
  subroutine init_sebomd_arrays(n)
   ! allocates SEBOMD arrays
   ! use sebomd_module, only : sebomd_obj
   implicit none
   integer, intent(in) :: n

   allocate(Zatms(n))
   allocate(Zatms_eff(n))
   ! PEPlist
   ! The MM correction of petide bonds applies
   ! a constraint on two dihedrals per amide function.
   ! We need to store two lists of four atoms per amide 
   ! group.
   ! One need at least 6 atoms to have an amide in the molecule.
   ! The maximum number of peptide bond present in a molecule of 
   ! n atoms is thus n/6. Since we need to store two dihedral per
   ! amide, the array has a maximum size of 2n/6 = n/3.
   allocate(PEPlist(n/3,4)) 
   allocate(interactions(n,n))

   return
  end subroutine init_sebomd_arrays

  subroutine put_Z(i,z)
   implicit none
   integer, intent(in) :: z
   integer, intent(in) :: i

   Zatms(i) = z

   return

  end subroutine put_Z

  subroutine put_Zeff(i,zeff)
   implicit none
   double precision, intent(in) :: zeff
   integer, intent(in) :: i

   Zatms_eff(i) = zeff

   return

  end subroutine put_Zeff

  subroutine get_Z(i,z)
   implicit none
   integer, intent(in)  :: i
   integer, intent(out) :: z

   z = Zatms(i)

   return

  end subroutine get_Z

  subroutine get_Zeff(i,zeff)
   implicit none
   integer, intent(in)  :: i
   double precision, intent(out) :: zeff

   zeff = Zatms_eff(i)

   return

  end subroutine get_Zeff

  subroutine put_PEP(ipep,natm,iatm)
   implicit none
   integer, intent(in) :: ipep,natm,iatm

   PEPlist(ipep,natm) = iatm

   return
  
  end subroutine put_PEP

  subroutine get_PEP(ipep,pep_b)
   implicit none
   integer, intent(in) :: ipep
   integer, intent(out) :: pep_b(4)

   pep_b(1) = PEPlist(ipep,1)
   pep_b(2) = PEPlist(ipep,2)
   pep_b(3) = PEPlist(ipep,3)
   pep_b(4) = PEPlist(ipep,4)
   
   return

  end subroutine get_PEP

  subroutine put_all_interactions(int_type)
   implicit none
   integer, intent(in) :: int_type

   interactions(:,:) = int_type

   return

  end subroutine

  subroutine put_interaction(i,j,int_type)
   implicit none
   integer, intent(in) :: i,j,int_type

   interactions(i,j) = int_type

   if (i.ne.j) interactions(j,i) = interactions(i,j)

   return

  end subroutine put_interaction

  subroutine get_interaction(i,j,int_type)
   implicit none
   integer, intent(in) :: i,j
   integer, intent(out) :: int_type

   int_type = interactions(i,j)

   return
  
  end subroutine get_interaction

  subroutine cleanup_sebomd_arrays()
   ! deallocates SEBOMD arrays
   implicit none

   if(allocated(Zatms)) deallocate(Zatms)
   if(allocated(Zatms_eff)) deallocate(Zatms_eff)
   if(allocated(PEPlist)) deallocate(PEPlist)
   if(allocated(interactions)) deallocate(interactions)

   return
  end subroutine cleanup_sebomd_arrays

end module sebomd_arrays
