#include "copyright.h"
#include "../include/assert.fh"
#include "../include/dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++
!This module contains data for the softcore
!potentials usable in TI calculations
!+++++++++++++++++++++++++++++++++++++++++++++

module softcore

  use findmask

  implicit none

#include "parallel.h"
#ifdef MPI
   include 'mpif.h'

  integer ist(MPI_STATUS_SIZE), partner
#endif
  _REAL_,save :: scalpha, scbeta                             ! parameter alpha in the SC-potential or beta for electrostatics
  integer, save :: ifsc, sceeorder, tishake
  character(len=256), save :: scmask
  integer,save :: nsoftcore, natom_partner, nsoftcore_partner, extra_atoms, nmixed      ! numbers of softcore atoms
  integer, dimension(:),allocatable,save :: nsc, nsc_partner                            ! keeps track of which atoms are softcore atoms
  integer, dimension(:),allocatable,save :: molecule_type, molecule_type_partner        ! keeps track of which molecules are partially softcore
  integer, dimension(:),allocatable,save :: frcmask_sc, frcmask, frcmask_partner        ! for efficient force mixing

  integer :: ier_alloc, ierr
  integer,save :: logdvdl, dvdl_pnt, dvdl_norest
  _REAL_, dimension(:),allocatable,save :: dvdl_values       ! used to accumulate all dvdl values
  _REAL_, dimension(:),allocatable,save :: foureps, sigma6   ! the van der waals parameters extracted from cn1 and cn2
  _REAL_, save :: weight0, weight1                           ! for linear mixing w0 is 1-l, w1 is l
  logical, save :: isProcessV1
  _REAL_,save :: sc_dvdl, sc_tot_dvdl,sc_tot_dvdl_partner    ! this is used to accumulate the dvdl contribution due to the 
                                                             ! softcore vdw potential which is added later in mix_frcti
  _REAL_,save :: sc_dvdl_ee, sc_tot_dvdl_ee,sc_tot_dvdl_partner_ee ! dvdl contribution due to softcore electrostatics

  integer, parameter :: ti_ene_cnt = 19           ! update as needed

  _REAL_,save :: oneweight, sc_ener(ti_ene_cnt)                      ! inverse of the lambda-weight, energies of softcore parts
  _REAL_,save :: sc_ener_ave(ti_ene_cnt), sc_ener_rms(ti_ene_cnt)
  _REAL_,save :: sc_ener_tmp(ti_ene_cnt), sc_ener_tmp2(ti_ene_cnt), sc_ener_old(ti_ene_cnt), sc_ener_old2(ti_ene_cnt) ! for running averages
  ! sc_energy array contains:
  ! (1) BOND
  ! (2) ANGLE
  ! (3) DIHEDRAL
  ! (4) 1-4 NB 
  ! (5) 1-4 EE
  ! (6) KINETIC ENERGY
  ! (7) INTERNAL VDW
  ! (8) VDW SOFTCORE DERIVATIVE TERM
  ! (9) INTERNAL ELECTROSTATICS
  ! (10) EE SOFTCORE DERIVATIVE TERM
  ! (11) ADDITIONAL DERIVATIVE TERM
  ! (12) POT ENERGY
  ! (13) TOT ENERGY
  ! (14-19) RESTRAINT ENERGY (type 1-6), following enmr

  _REAL_,save :: sc_temp, boltzmann_factor          ! T and Ekin of only the SC atoms factor is 1/2*k(cal)*deg
  _REAL_,save :: dynlmb                             ! for dynamically rising lambda, the amount of increase
  integer,save :: sc_deg, sc_dof_shaked             ! degrees of freedom due to the softcore atoms
  integer,save :: res_ti_region                     ! used to define/decouple restraints (see ti_check_res)
  integer,save :: emil_sc                           ! for EMIL calculations
  contains

!===================================================================================================

#ifdef MPI
    subroutine setup_sc(natom,nres,igraph,isymbl,ipres, &
         lbres,crd, ntypes, clambda, nstlim)

      integer, intent(in) :: natom, nres, ipres(*),ntypes, nstlim
      character(len=4), intent(in) :: igraph(*), isymbl(*), lbres(*)
      _REAL_, intent(in) :: crd(*), clambda

      call sc_allocate_arrays(natom, ntypes, nstlim)

      call atommask(natom,nres,0,igraph,isymbl,ipres,  &
           lbres,crd,scmask,nsc)
      nsoftcore=sum(nsc)

      if (sanderrank==0) then
         write(6,*) '      '
         write(6,'(a,a,a,i5,a)')  '     Softcore Mask ', scmask(1:len_trim(scmask)), &
              ' matches ', nsoftcore,' atoms'
      end if

      if (sanderrank==0 .and. ifsc /=2 ) then
         call mpi_barrier( commmaster, ierr)
         partner = ieor(masterrank,1)
         call mpi_sendrecv( natom, 1, MPI_INTEGER, partner, 5, &
                            natom_partner, 1, MPI_INTEGER, partner, 5, &
                            commmaster, ist, ierr )

         write(6,'(a,i1,a,i1)') '     this run corresponds to V',masterrank, &
                  ', its softcore atoms interact fully for lambda=',masterrank
         write(6,'(a,i6,a,i6,a)') '     this process: ',natom, &
                  ' atoms, partner process: ',natom_partner,' atoms'
         extra_atoms = natom_partner - natom
         isProcessV1 = (masterrank == 1)

         allocate (nsc_partner(natom+extra_atoms),stat=ier_alloc)

         if (ier_alloc /= 0) call sander_bomb('allocate_arrays[softcore.f]', &
              'cant allocate nsc_partner','')

         call mpi_barrier( commmaster, ierr)

         call mpi_sendrecv( nsc, natom, MPI_INTEGER, partner, 5, &
                    nsc_partner, natom+extra_atoms, MPI_INTEGER, partner, 5, &
                    commmaster, ist, ierr )

         nsoftcore_partner=sum(nsc_partner)
         nmixed=natom-nsoftcore
         if (nmixed /= natom_partner-nsoftcore_partner) then
            call sander_bomb('setup_sc', &
                 'All non-softcore atoms must be identical in both systems','')
         end if

         call build_force_masks(natom)

      end if

      if ( ifsc /= 2 ) then 
         call mpi_bcast(extra_atoms,1,MPI_INTEGER,0,commsander,ierr)
         call mpi_bcast(isProcessV1,1,MPI_LOGICAL,0,commsander,ierr)
      end if

      if (sanderrank==0 .and. ifsc /= 2) then
         call sc_check_and_adjust_crd(natom, crd)
      end if

      call calc_softcore_parameters(ntypes)

      if (isProcessV1) then
         oneweight = 1.0d0 / clambda
      else
         oneweight = 1.0d0 / ( 1.0d0 - clambda )
      end if

      weight0 = 1.0d0 - clambda
      weight1 = 1.0d0 - weight0

      sc_dvdl=0.0d0
      sc_tot_dvdl=0.0d0
      sc_tot_dvdl_partner=0.0d0

      sc_dvdl_ee=0.0d0
      sc_tot_dvdl_ee=0.0d0
      sc_tot_dvdl_partner_ee=0.0d0

      ! Zero energy arrays
      sc_ener(1:ti_ene_cnt) = 0.0d0
      sc_ener_ave(1:ti_ene_cnt) = 0.0d0
      sc_ener_rms(1:ti_ene_cnt) = 0.0d0
      sc_ener_tmp(1:ti_ene_cnt) = 0.0d0
      sc_ener_tmp2(1:ti_ene_cnt) = 0.0d0
      sc_ener_old(1:ti_ene_cnt) = 0.0d0
      sc_ener_old2(1:ti_ene_cnt) = 0.0d0

     RETURN 
    end subroutine setup_sc

!===================================================================================================

    subroutine sc_degrees_o_freedom(ndfmin)

      integer, intent(in) :: ndfmin

      if ( nsoftcore > 0 ) then
         if (ifsc==2) then
            ! System is completely SC atoms
            sc_deg = ( 3 * nsoftcore ) - ndfmin
         else
            ! For now just assume that each SC atoms has its full 3 degrees
            ! of freedom w no additional internal constraints
            sc_deg = 3 * nsoftcore
         end if
         ! Adjust for SHAKE constraints
         sc_deg = sc_deg - sc_dof_shaked
         if (sc_dof_shaked==0) then
            if (sanderrank==0) &
               write(6,'(a,i4)')  '   DOF for the SC part of the system: ', sc_deg
         else
            if (sanderrank==0) &
               write(6,'(a,i4,a,i4)')  '   DOF for the SC part of the system: ', &
                 sc_deg, ' SHAKE constraints in the SC region: ', sc_dof_shaked
         end if
         boltzmann_factor = 2 * 4.184 / ( 8.31441d-3 * sc_deg )
      else
         boltzmann_factor = 0.0d0
      end if

      RETURN
    end subroutine sc_degrees_o_freedom

#endif

!===================================================================================================

    subroutine sc_allocate_arrays(natom, ntypes, nstlim)

      integer, intent(in) :: natom, ntypes, nstlim

      allocate (nsc(natom),stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('allocate_arrays[softcore.f]','cant allocate nsc','')

      allocate (foureps(ntypes*ntypes),stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('allocate_arrays[softcore.f]','cant allocate foureps','')

      allocate (sigma6(ntypes*ntypes),stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('allocate_arrays[softcore.f]','cant allocate sigma6','')

      if (logdvdl ==1) then
         dvdl_pnt=0
         allocate (dvdl_values(nstlim),stat=ier_alloc)
         if (ier_alloc /= 0) call sander_bomb('allocate_arrays[softcore.f]', &
              'cant allocate dvdl_values','')
      end if

      RETURN
    end subroutine sc_allocate_arrays

!===================================================================================================

    subroutine calc_softcore_parameters(ntypes)

      use parms, only: cn1, cn2

      integer, intent(in) :: ntypes
      integer i

      do i=1,((ntypes+1)*ntypes)/2
         if (cn1(i) .ne. 0.0d0) then ! catch zero vdW hydrogens
            sigma6(i)=cn2(i)/cn1(i)
            foureps(i)=cn2(i)*sigma6(i)
         else
            sigma6(i)=1.0d0          ! will always be multiplied by zero
            foureps(i)=0.0d0
         end if
      end do

      RETURN
    end subroutine calc_softcore_parameters

!===================================================================================================

#ifdef MPI
    subroutine sc_change_clambda(clambda)
      
      _REAL_,intent(in) :: clambda

      if (isProcessV1) then
         oneweight = 1.0d0 / clambda
      else
         oneweight = 1.0d0 / ( 1.0d0 - clambda )
      end if

      weight0 = 1.0d0 - clambda
      weight1 = 1.0d0 - weight0

      if (numtasks > 1) then
         call mpi_bcast(clambda,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)
         call mpi_bcast(oneweight,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)
         call mpi_bcast(weight0,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)
         call mpi_bcast(weight1,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)
      end if

    end subroutine sc_change_clambda
#endif

!===================================================================================================

    subroutine build_force_masks(natom)

      integer natom,i,j,k

      allocate (frcmask_sc(3*nsoftcore),stat=ier_alloc)
      allocate (frcmask(3*nmixed),stat=ier_alloc)
      allocate (frcmask_partner(3*nmixed),stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('abuild_force_masks[softcore.f]', &
           'cant allocate masks','')

      j=1
      k=1
      do i=1,natom
         if (nsc(i)==1) then
            frcmask_sc(j)   = 3*(i-1)+1
            frcmask_sc(j+1) = 3*(i-1)+2
            frcmask_sc(j+2) = 3*(i-1)+3
            j=j+3
         else
            frcmask(k)   = 3*(i-1)+1
            frcmask(k+1) = 3*(i-1)+2
            frcmask(k+2) = 3*(i-1)+3
            k=k+3
         end if
      end do
      k=1
      do i=1,natom_partner
         if (nsc_partner(i)==1) then
         else
            frcmask_partner(k)   = 3*(i-1)+1
            frcmask_partner(k+1) = 3*(i-1)+2
            frcmask_partner(k+2) = 3*(i-1)+3
            k=k+3
         end if
      end do

      RETURN
    end subroutine build_force_masks

!===================================================================================================

    subroutine mix_frc_sc(dvdl,w0,w1,f,fcopy)

      _REAL_ dvdl,w0,w1
      _REAL_ f(*),fcopy(*)
      integer i

      ! First, adjust dvdl for the softcore contributions

      if (isProcessV1) then
         sc_ener(8) = w1 * sc_tot_dvdl + w0 * sc_tot_dvdl_partner
         sc_ener(10) = w1 * sc_tot_dvdl_ee + w0 * sc_tot_dvdl_partner_ee
      else
         sc_ener(8) = w1 * sc_tot_dvdl_partner + w0 * sc_tot_dvdl
         sc_ener(10) = w1 * sc_tot_dvdl_partner_ee + w0 * sc_tot_dvdl_ee
      end if

      ! SC_DER = SC_EEL_DER + SC_VDW_DER
      sc_ener(11) = sc_ener(8) + sc_ener(10)

      ! Add the pot energy terms here - needed for minimization 
      sc_ener(12) = sc_ener(1) + sc_ener(2) + sc_ener(3) + sc_ener(4) + &
                    sc_ener(5) + sc_ener(7) + sc_ener(9)      

      dvdl = dvdl + sc_ener(11)

      sc_tot_dvdl = 0.0d0
      sc_tot_dvdl_partner = 0.0d0

      sc_tot_dvdl_ee = 0.0d0
      sc_tot_dvdl_partner_ee = 0.0d0

      ! Then, mix the forces using the dual-topology approach
      if ( isProcessV1 ) then
         do i=1,3*nmixed
            f(frcmask(i)) = w1 * f(frcmask(i)) + w0 * fcopy(frcmask_partner(i))
         end do
         do i=1,3*nsoftcore
            ! The following scales down ALL forces on appearing softcore atoms
            ! non-lambda depended SC forces are scaled UP by the same amount in ene.f
            f(frcmask_sc(i)) = w1 * f(frcmask_sc(i))
         end do
      else
         do i=1,3*nmixed
            fcopy(frcmask(i)) = w1 * f(frcmask_partner(i)) + w0 * fcopy(frcmask(i))
         end do
         do i=1,3*nsoftcore
            ! The following scales down ALL forces on disappearing softcore atoms
            ! non-lambda depended SC forces are scaled UP by the same amount in ene.f
            fcopy(frcmask_sc(i)) = w0 * fcopy(frcmask_sc(i))
         end do
      end if

      RETURN
    end subroutine mix_frc_sc

!===================================================================================================

#ifdef MPI
    subroutine sc_check_and_adjust_crd(natom, crd)

      ! Input variables
      integer, intent(in) :: natom
      _REAL_ crd(3*natom)

      ! Local
      _REAL_ crd_partner(3*(natom+extra_atoms))
      integer i, nadj

      ! Get the coordinates from the partner process
      call mpi_sendrecv( crd, 3*natom, MPI_DOUBLE_PRECISION, partner, 5, &
            crd_partner, 3*(natom+extra_atoms), MPI_DOUBLE_PRECISION, partner, &
            5, commmaster, ist, ierr )

      ! Go over all common coordinates and check for deviations
      ! deviating by more than 0.1 A leads to termination

      write(6,'(a)') '     Checking for mismatched coordinates.'
      nadj=0

      do i=1,3*nmixed
         if ( crd(frcmask(i)) /= crd_partner(frcmask_partner(i)) ) then
            if (abs(crd(frcmask(i))-crd_partner(frcmask_partner(i))) > 0.1 ) then
               write (6,'(a,i5,a,i5,a)') '     WARNING: Local coordinate ', &
                   frcmask(i),' differs from partner coordinate ', &
                   frcmask_partner(i),' !'
               call sander_bomb('sc_check_and_adjust', &
                   'Atom coordinate disagreement','Check input files.')
            else
               nadj = nadj + 1
               if (isProcessV1) then
                  crd(frcmask(i)) = crd_partner(frcmask_partner(i))
               end if
               if (nadj<11) then
                  if (nadj<10) then
                     write (6,'(a,i5,a,i5,a)') &
                         '     WARNING: Local coordinate ',frcmask(i), &
                         ' differs from partner coordinate ', &
                         frcmask_partner(i),' !'
                     if (isProcessV1) then
                        write (6,'(a)') '     Deviation is small, using partner coordinate.'
                     else
                        write (6,'(a)') '     Deviation is small, changing partner coordinate.'
                     end if
                  else
                     write (6,'(a)') '     ... making more adjustments ...'
                  end if
               end if
            end if
         else
            ! Everything ok, do nothing
         end if
      end do
      if (nadj>9) then
         write (6,'(a,i5,a)') '     A total of ',nadj, &
              ' small coordinate adjustments were made, check results carefully.'
      end if

    end subroutine sc_check_and_adjust_crd

!===================================================================================================

    subroutine sc_check_perturbed_molecules(nummols, molsiz)

      integer nummols,molsiz(*)

      integer nummols_partner
      integer i, iatom, num, imol, temp

      ! Check for every molecule if it is *partially* softcore
      ! on either this or the partner side

      ! First communicate the number of molecules

      !write (6,*) 'Debugging info for pscale:'
      
      call mpi_sendrecv( nummols, 1, MPI_INTEGER, partner, 5, &
           nummols_partner, 1, MPI_INTEGER, partner, 5, &
           commmaster, ist, ierr )

      ! Prepare the molecule type array
      allocate (molecule_type(nummols),stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('sc_check_perturbed_molecules[softcore.f]', &
           'cant allocate molecule_type','')

      allocate (molecule_type_partner(nummols_partner),stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('sc_check_perturbed_molecules[softcore.f]', &
           'cant allocate molecule_type_partner','')

      !write (6,*) nummols, 'Molecules on this side, ', nummols_partner, ' on the partner side'

      ! check which molecules are partially softcore
      i=0
      do imol = 1, nummols
         temp=0
         num = molsiz(imol)
         do iatom = 1, num
            i = i + 1
            temp = temp + nsc(i)
         end do
         if (temp == 0) then
            ! no atoms in this molecule are softcore
            ! its C.O.M. needs not be communicated later
            ! unless the partner molecule is partially softcore
            molecule_type(imol)=0
         else 
            if (temp == num) then
               ! all atoms in this molecule are softcore
               ! its C.O.M. needs not be communicated later
               molecule_type(imol) = 2
            else
               ! some atoms in this molecule are softcore
               ! its C.O.M. need to be communicated later
               molecule_type(imol) = 1
            end if
         end if
      end do

      ! exchange the molecule type arrays

      call mpi_sendrecv( molecule_type, nummols, MPI_INTEGER, partner, 5, &
           molecule_type_partner, nummols_partner, MPI_INTEGER, partner, 5, &
           commmaster, ist, ierr )

      ! now go through the arrays and check which non-soft core 
      ! molecule on this side has a partially soft core molecule 
      ! on the partner side

      i=0
      do imol = 1, nummols
         if (molecule_type(imol) == 2) then
            ! this molecule is completely softcore on this side
            ! nothing corresponds to it on the partner -> do nothing
         else
            i=i+1
            if (molecule_type_partner(i) == 2) then
               ! skip all complete softcore partner molecules
               do
                  i = i + 1
                  if (molecule_type_partner(i) /= 2) then
                     exit
                  end if
               end do
            end if
            if (molecule_type(imol) == 0 .and. molecule_type_partner(i) == 1) then
               ! This molecule has no softcore atoms but its partner molecule has
               ! therefore prepare for sending its C.O.M. later
               molecule_type(imol) = 1
               write (6,'(a,i5,a,a)') '     Molecule ',imol, ' is nonsoftcore, but ', &
                    'its partner is partially softcore, therefore its C.O.M. will be exchanged.' 
            end if
         end if
      end do

      !write (6,'(a)') 'Final molecule types:'
      
      imol=0
      do i = 1, nummols
         if (molecule_type(i) == 2) then
            write (6,'(a,i5,a)') '     Molecule ',i, &
                 ' is completely softcore and skipped for C.O.M..'
         end if
         if (molecule_type(i) == 1) then
            write (6,'(a,i5,a)') '     Molecule ',i, &
                 ' is partially softcore on this side or the corresponding partner molecule is.'
         end if
         if (molecule_type(i) == 0) then
            imol = imol + 1
         end if
      end do
      !write (6, '(i5,a)') imol,'Molecules are nonsoftcore on both sides'

      return

    end subroutine sc_check_perturbed_molecules
!===================================================================================================

    subroutine sc_pscale(x,amass,nummols,molsiz,oldrecip,ucell)
      ! This subroutine is copied from ew_pscale, but includes master-to-master
      ! communication to get average center of masses for perturbed molecules

      _REAL_ x(3,*),amass(*)
      _REAL_ oldrecip(3,3),ucell(3,3)
      integer nummols,molsiz(*)
      integer imol,iatom,i,num,isave
      _REAL_ mass,massmol,xmol,ymol,zmol,fm1,fm2,fm3, &
           xmolnu,ymolnu,zmolnu
      _REAL_ xmol_partner, ymol_partner, zmol_partner

      partner = ieor(masterrank,1)

      ! apply C.O.M. based pressure scaling, atom based is not available here

      i = 0
      do imol = 1,nummols
         massmol = 0.d0
         xmol = 0.d0
         ymol = 0.d0
         zmol = 0.d0
         num = molsiz(imol)
         isave = i
         
         ! get c.o.m. of molecule

         do iatom = 1,num
            i = i + 1
            mass = amass(i)
            massmol = massmol + mass
            xmol = xmol + mass*x(1,i)
            ymol = ymol + mass*x(2,i)
            zmol = zmol + mass*x(3,i)
         end do

         xmol = xmol / massmol
         ymol = ymol / massmol
         zmol = zmol / massmol

         if (molecule_type(imol) == 1) then 
            ! calculate the average C.O.M. of the V0 and V1 versions of the molecule
            call mpi_sendrecv( xmol, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                 xmol_partner, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                 commmaster, ist, ierr )
            call mpi_sendrecv( ymol, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                 ymol_partner, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                 commmaster, ist, ierr )
            call mpi_sendrecv( zmol, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                 zmol_partner, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                 commmaster, ist, ierr )
            xmol = (xmol + xmol_partner) / 2.0d0
            ymol = (ymol + ymol_partner) / 2.0d0
            zmol = (zmol + zmol_partner) / 2.0d0
         end if

         !   ...now get fracs for c.o.m. using old cell params
         
         fm1 = xmol*oldrecip(1,1)+ymol*oldrecip(2,1)+ &
               zmol*oldrecip(3,1)
         fm2 = xmol*oldrecip(1,2)+ymol*oldrecip(2,2)+ &
               zmol*oldrecip(3,2)
         fm3 = xmol*oldrecip(1,3)+ymol*oldrecip(2,3)+ &
               zmol*oldrecip(3,3)
         !   ...use these with new cell params to get new c.o.m. cartesians
         
         xmolnu = fm1*ucell(1,1)+fm2*ucell(1,2)+fm3*ucell(1,3)
         ymolnu = fm1*ucell(2,1)+fm2*ucell(2,2)+fm3*ucell(2,3)
         zmolnu = fm1*ucell(3,1)+fm2*ucell(3,2)+fm3*ucell(3,3)
         i = isave
         
         !   ...now rigidly translate molecule
         do iatom = 1,num
            i = i + 1
            x(1,i) = x(1,i) + xmolnu - xmol
            x(2,i) = x(2,i) + ymolnu - ymol
            x(3,i) = x(3,i) + zmolnu - zmol
         end do
      end do  !  imol = 1,nummols

   return
end subroutine sc_pscale 

!===================================================================================================

    subroutine mix_temp_scaling(factor, v)

      _REAL_ factor, v(*), factor_sc
      _REAL_ factor_partner
      integer i

      partner = ieor(masterrank,1)
      call mpi_sendrecv( factor, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                         factor_partner, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                         commmaster, ist, ierr )

      factor_sc = factor

      if (isProcessV1) then
         factor = weight1 * factor + weight0 * factor_partner
      else
         factor = weight1 * factor_partner + weight0 * factor
      end if

      factor_sc = factor_sc / factor

      do i=1,3*nsoftcore
         ! The following results in softcore atom velocities effectively being scaled only by
         ! the temp-factor of their own process. This prevents instabilities at high/low lambdas
         v(frcmask_sc(i)) = factor_sc * v(frcmask_sc(i))
      end do

      RETURN
    end subroutine mix_temp_scaling

!===================================================================================================

    subroutine sc_mix_sum(sum)

      _REAL_ sum, sum_partner

      partner = ieor(masterrank,1)
      call mpi_sendrecv( sum, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                         sum_partner, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                         commmaster, ist, ierr )

      if (isProcessV1) then
         sum = weight1 * sum + weight0 * sum_partner
      else
         sum = weight1 * sum_partner + weight0 * sum
      end if

      RETURN
    end subroutine sc_mix_sum
#endif

!===================================================================================================

    subroutine sc_lngdyn(winv,amass,v,f,sdfac, c_explic, c_implic, istart, iend, nr, dtx)

      use random, only: gauss
      ! Passed in
      integer istart, iend, nr
      _REAL_ winv(*),amass(*),v(*),f(*),sdfac,c_explic,c_implic, dtx

      ! Local
      integer j, i3, k
      _REAL_ wfac, aamass, rsd, fln, fln_buffer(3*nsoftcore)

      k=1
      do j=1, nr
         if (nsc(j)==1) then
            aamass = amass(j)
            rsd = sdfac*sqrt(aamass)
            call gauss( 0.d0, rsd, fln )
            fln_buffer(k)=fln
            call gauss( 0.d0, rsd, fln )
            fln_buffer(k+1)=fln
            call gauss( 0.d0, rsd, fln )
            fln_buffer(k+2)=fln
            k=k+3
         end if
      end do

      do j=1, extra_atoms
         call gauss( 0.d0, 1.d0, fln )
         call gauss( 0.d0, 1.d0, fln )
         call gauss( 0.d0, 1.d0, fln )
      end do

      k=1

      do j=1, istart-1
         if (nsc(j)==0) then
            call gauss( 0.d0, 1.d0, fln )
            call gauss( 0.d0, 1.d0, fln )
            call gauss( 0.d0, 1.d0, fln )
         else
            k=k+3
         end if
      end do

      i3 = 3*(istart-1)
      do j=istart, iend

         wfac = winv(j) * dtx
         aamass = amass(j)
         
         rsd = sdfac*sqrt(aamass)

         if (nsc(j)==0) then
            call gauss( 0.d0, rsd, fln )
            v(i3+1) = (v(i3+1)*c_explic + (f(i3+1)+fln)*wfac) * c_implic
            call gauss( 0.d0, rsd, fln )
            v(i3+2) = (v(i3+2)*c_explic + (f(i3+2)+fln)*wfac) * c_implic
            call gauss( 0.d0, rsd, fln )
            v(i3+3) = (v(i3+3)*c_explic + (f(i3+3)+fln)*wfac) * c_implic      
         else
            v(i3+1) = (v(i3+1)*c_explic + (f(i3+1)+fln_buffer(k))*wfac) * c_implic
            v(i3+2) = (v(i3+2)*c_explic + (f(i3+2)+fln_buffer(k+1))*wfac) * c_implic
            v(i3+3) = (v(i3+3)*c_explic + (f(i3+3)+fln_buffer(k+2))*wfac) * c_implic
            k=k+3
         end if
         i3 = i3+3
      end do

      do j=iend+1, nr
         if (nsc(j)==0) then
            call gauss( 0.d0, 1.d0, fln )
            call gauss( 0.d0, 1.d0, fln )
            call gauss( 0.d0, 1.d0, fln )
         end if
      end do

      RETURN
    end subroutine sc_lngdyn

!===================================================================================================

#ifdef MPI
    subroutine sc_mix_position(vcmx,vcmy,vcmz,clambda)

      _REAL_ vcmx,vcmy,vcmz,clambda
      _REAL_ vcm(3), vcm_partner(3)

      vcm(1)=vcmx
      vcm(2)=vcmy
      vcm(3)=vcmz

      partner = ieor(masterrank,1)
      call mpi_sendrecv( vcm, 3, MPI_DOUBLE_PRECISION, partner, 5, &
                         vcm_partner, 3, MPI_DOUBLE_PRECISION, partner, 5, &
                         commmaster, ist, ierr )

      if (isProcessV1) then
         vcm(1:3) = clambda * vcm(1:3) + (1.0d0-clambda) * vcm_partner(1:3)
      else
         vcm(1:3) = clambda * vcm_partner(1:3) + (1.0d0-clambda) * vcm(1:3)
      end if

      vcmx=vcm(1)
      vcmy=vcm(2)
      vcmz=vcm(3)

      RETURN
    end subroutine sc_mix_position
#endif

!===================================================================================================

#ifdef MPI
    subroutine sc_mix_velocities(v, nr3)

      _REAL_ v(nr3),vcopy(nr3+3*extra_atoms)
      integer nr3, i

      partner = ieor(masterrank,1)
      call mpi_sendrecv( v, nr3, MPI_DOUBLE_PRECISION, partner, 5, &
                         vcopy, nr3+3*extra_atoms, MPI_DOUBLE_PRECISION, &
                         partner, 5, commmaster, ist, ierr )

      ! Mix velocities of shared atoms that have become out of sync
      ! because of COM removal.

      if ( isProcessV1 ) then
         do i=1,3*nmixed
            v(frcmask(i)) = weight1 * v(frcmask(i)) + weight0 * vcopy(frcmask_partner(i))
         end do
         do i=1,3*nsoftcore
            ! The following would scale down velocities on non-shared atoms
            ! v(frcmask_sc(i)) = weight1 * v(frcmask_sc(i))
         end do
      else
         do i=1,3*nmixed
            v(frcmask(i)) = weight1 * vcopy(frcmask_partner(i)) + weight0 * v(frcmask(i))
         end do
         do i=1,3*nsoftcore
            ! The following would scale down velocities on non-shared atoms
            ! v(frcmask_sc(i)) = weight0 * v(frcmask_sc(i))
         end do
      end if

      RETURN
    end subroutine sc_mix_velocities
#endif

!===================================================================================================

#ifdef MPI
    subroutine sc_sync_x(x,nr3)
      
      _REAL_ x(nr3), xpartner(nr3+3*extra_atoms)
      integer nr3, i

      partner = ieor(masterrank,1)
      call mpi_sendrecv( x, nr3, MPI_DOUBLE_PRECISION, partner, 5, &
                         xpartner, nr3+3*extra_atoms, MPI_DOUBLE_PRECISION, &
                         partner, 5, commmaster, ist, ierr )

      ! Synchronize coordinates to those of Process V0
      ! softcore atoms can be ignored here
      if (isProcessV1) then
         do i=1,3*nmixed
            x(frcmask(i)) = xpartner(frcmask_partner(i))
         end do
      end if

      RETURN
    end subroutine sc_sync_x

#endif

!===================================================================================================

#ifdef MPI
    subroutine sc_compare(x,nr3, tag)
      
      _REAL_ x(nr3), xpartner(nr3+3*extra_atoms)
      integer nr3, i
      _REAL_ delta, delta2
      character(len=3) tag

      partner = ieor(masterrank,1)
      call mpi_sendrecv( x, nr3, MPI_DOUBLE_PRECISION, partner, 5, &
                         xpartner, nr3+3*extra_atoms, MPI_DOUBLE_PRECISION, &
                         partner, 5, commmaster, ist, ierr )

      ! Check for differences between coords or velocities
      ! softcore atoms can be ignored here
      delta2 = 0.0d0
      do i=1,3*nmixed
         delta =  x(frcmask(i)) - xpartner(frcmask_partner(i))
         delta2 = delta2 + delta*delta
         if ( delta > 0.001 .or.  delta < -0.001 ) then
            ! stop the run if major desynchronisation has occured
            write (6,'(a,a,a,i6,a,f18.14)') 'ERROR: Large sync deviation of ', &
                 tag, ' coordinate ', i, ' delta ', delta
            call sander_bomb('sc_compare[softcore.F90]','MPI processes are out of sync','')
         end if
      end do
      ! Print a warning if coordinates have desynchronized
      if (delta2 > 1.0d-20) then
         write (6,'(a,a,a)') 'WARNING: MPI processes of V0 and V1 have slightly desynchronized.', &
              'This can be caused by compiler loss of precision errors, ', &
              'check your results very carefully if this warning occurs repeatedly'
      end if
      
      RETURN
    end subroutine sc_compare
#endif

!===================================================================================================

    subroutine adj_dvdl_stat(edvdl, edvdl_r)
      use state
      type(state_rec)  :: edvdl, edvdl_r
      _REAL_ vdw_val, ee_val, tot_val
   
      if ( isProcessV1 ) then
         vdw_val = weight1 * sc_tot_dvdl + weight0 * sc_tot_dvdl_partner
         ee_val  = weight1 * sc_tot_dvdl_ee + weight0 * sc_tot_dvdl_partner_ee
      else
         vdw_val = weight1 * sc_tot_dvdl_partner + weight0 * sc_tot_dvdl
         ee_val  = weight1 * sc_tot_dvdl_partner_ee + weight0 * sc_tot_dvdl_ee
      end if

      edvdl%pot%vdw = edvdl%pot%vdw + vdw_val
      edvdl_r%pot%vdw = edvdl_r%pot%vdw + vdw_val

      edvdl%pot%elec = edvdl%pot%elec + ee_val
      edvdl_r%pot%elec = edvdl_r%pot%elec + ee_val

      tot_val = vdw_val + ee_val
      edvdl%pot%tot = edvdl%pot%tot + tot_val
      edvdl_r%pot%tot = edvdl_r%pot%tot + tot_val

      RETURN
    end subroutine adj_dvdl_stat

!===================================================================================================

    subroutine log_dvdl(dvdl)

      _REAL_ dvdl

      if (logdvdl /= 0) then
         dvdl_pnt = dvdl_pnt+1
         dvdl_values(dvdl_pnt)=dvdl;
      end if

      RETURN
    end subroutine log_dvdl

!===================================================================================================

    subroutine sc_print_dvdl_values()

      integer i

      if (logdvdl /= 0) then
         write (6,'(a,i8,a)') 'Summary of dvdl values over ',dvdl_pnt,' steps:'
         do i=1, dvdl_pnt
            write (6,'(f11.4)') dvdl_values(i)
         end do
         write (6,'(a)') 'End of dvdl summary'
      end if

      RETURN
    end subroutine sc_print_dvdl_values

!===================================================================================================

    subroutine sc_nomix_frc(f,nr3,ener)
      use state

      _REAL_ f(*)
      type(state_rec) :: ener
      _REAL_ dvdl
      integer nr3

      ! 'partner' etot is zero, so dvdl is minus epot plus sc-contributions
      dvdl= -ener%pot%tot + weight0 * sc_tot_dvdl
      sc_ener(8) = weight0 * sc_tot_dvdl
      sc_ener(10) = weight0 * sc_tot_dvdl_ee
      sc_ener(11) = sc_ener(8) + sc_ener(10)
      ! include sc_dvdl into vdw and potential energies:
      ener%pot%tot = ener%pot%tot + weight0 * sc_tot_dvdl
      ener%pot%vdw = ener%pot%vdw + weight0 * sc_tot_dvdl


      ! scale down energies and forces
      ener = ener * weight0
      f(1:nr3) = weight0 * f(1:nr3)
      ener%pot%dvdl = dvdl    

      RETURN
    end subroutine sc_nomix_frc

!===================================================================================================

#ifdef MPI
    subroutine calc_softcore_ekin(amass,v,vold,istart,iend)

      integer istart,iend
      _REAL_ amass(*), v(*), vold(*), sc_ekin

      !local
      integer j,m,i3,ierr
      _REAL_ aamass, temp

      i3= 3 * (istart-1)
      sc_ekin=0.0d0

      do j = istart,iend
         aamass = amass(j)
         do m = 1,3
            i3 = i3+1
            if (nsc(j) == 1) then
               temp = aamass * 0.25d0 * ( v(i3) + vold(i3) )**2
               sc_ekin = sc_ekin + temp
            end if
         end do
      end do
      sc_ekin = sc_ekin * 0.5d0

      ! collect ekin from the nodes
      temp=0.0d0
      call mpi_reduce(sc_ekin, temp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, commsander, ierr)
      sc_ekin = temp

      sc_ener(6) = sc_ekin

      RETURN
    end subroutine calc_softcore_ekin
#endif

!===================================================================================================

subroutine ti_check_res(nmrat, nmrcom, resttype, i)

  implicit none

! Formal arguments:
  ! The nmr code will read past what is defined here as the bounds for the array
  integer, intent(in) :: nmrat(16, *)
  integer, intent(in) :: nmrcom(2, *)
  integer, intent(in) :: resttype(*)
  integer, intent(in) :: i  !index into nmrat array

! Local variables:
  integer             :: iat
  integer             :: jat
  integer             :: ip
  integer             :: j
  integer             :: max_iat
  integer             :: iat_num

  res_ti_region = 0

  if (resttype(i) .le. 4) then
    max_iat = resttype(i) + 1
  else 
    max_iat = 8
  end if

  do iat = 1, max_iat
     jat = nmrat(iat, i)
     if (jat .lt. 0) then
        do ip = -nmrat(iat, i), -nmrat(8 + iat, i)
           do j = nmrcom(1, ip), nmrcom(2, ip)
              iat_num = j
              if(nsc(iat_num) .ne. 0) then
                 res_ti_region = 1
                 return
              end if
           end do
        end do
     else
        iat_num = (jat/3) + 1
        if(nsc(iat_num) .ne. 0) then
           res_ti_region = 1
           return
        end if      
     end if
  end do

  return

end subroutine ti_check_res

!===================================================================================================

    subroutine sc_print_energies(channel, energies)

      integer channel
      _REAL_ energies(ti_ene_cnt)

      sc_temp = energies(6) * boltzmann_factor

      if (nsoftcore>0) then
         write (channel,900) nsoftcore, sc_temp
         write (channel,950) energies(13),energies(6),energies(12)
         write (channel,1000) energies(1),energies(2),energies(3)
         write (channel,1100) energies(4),energies(5),energies(7)
         write (channel,1200) energies(9)
         write (channel,1300) energies(14),energies(15),energies(16)
         write (channel,1400) energies(17),energies(18),energies(19)
         write (channel,1600) energies(10),energies(8),energies(11)
         write (channel,2000)
      end if

       900 format (2x,'Softcore part of the system: ',i5,' atoms,',9x,'TEMP(K)    = ',f14.2)
       950 format (1x,'SC_Etot=    ',f11.4,2x,'SC_EKtot=    ',f11.4,2x,'SC_EPtot   =    ',f11.4)
      1000 format (1x,'SC_BOND=    ',f11.4,2x,'SC_ANGLE=    ',f11.4,2x,'SC_DIHED   =    ',f11.4)
      1100 format (1x,'SC_14NB=    ',f11.4,2x,'SC_14EEL=    ',f11.4,2x,'SC_VDW     =    ',f11.4)
      1200 format (1x,'SC_EEL =    ',f11.4)
      1300 format (1x,'SC_RES_DIST=',f11.4,2x,'SC_RES_ANG=  ',f11.4,2x,'SC_RES_TORS=    ',f11.4)
      1400 format (1x,'SC_RES_PLPT=',f11.4,2x,'SC_RES_PLPL= ',f11.4,2x,'SC_RES_GEN =    ',f11.4)
      1600 format (1x,'SC_EEL_DER= ',f11.4,2x,'SC_VDW_DER=  ',f11.4,2x,'SC_DERIV   =    ',f11.4)
      2000 format (t2,78('-'),/)
      RETURN
    end subroutine sc_print_energies

!===================================================================================================

    subroutine summarize_ti_changes(natom, resat)

      integer natom, i
      character(len=14) resat(natom)

      write (6,'(a)') '      TI atoms summary'
      write (6,*)

      do i=1,natom
         if ( nsc(i) == 1 ) then
            write(6,'(a,i6,a,a)') ' Atom: ', i, ' - ',resat(i)(1:13)
         end if
      end do

      write (6,2000)

      2000 format (t2,78('-'),/)
      RETURN
    end subroutine summarize_ti_changes

!===================================================================================================

    subroutine cleanup_sc()
      if(allocated(nsc)) deallocate (nsc,stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('cleanup_sc[softcore.f]','cant deallocate nsc','')

      if(allocated(foureps)) deallocate (foureps,stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('cleanup_sc[softcore.f]','cant deallocate foureps','')

      if(allocated(sigma6)) deallocate (sigma6,stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('cleanup_sc[softcore.f]','cant deallocate sigma6','')

      if(allocated(nsc_partner)) deallocate (nsc_partner, stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('cleanup_sc[softcore.f]', &
           'cant deallocate nsc_partner','')

      if(allocated(frcmask)) deallocate (frcmask, stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('cleanup_sc[softcore.f]','cant deallocate frcmask','')

      if(allocated(frcmask_sc)) deallocate (frcmask_sc, stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('cleanup_sc[softcore.f]', &
           'cant deallocate frcmask_sc','')

      if(allocated(frcmask_partner)) deallocate (frcmask_partner, stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('cleanup_sc[softcore.f]', &
           'cant deallocate frcmask_partner','')

      if(allocated(dvdl_values)) deallocate (dvdl_values, stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('cleanup_sc[softcore.f]', &
           'cant deallocate dvdl_values','')

      if(allocated(molecule_type)) deallocate (molecule_type, stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('cleanup_sc[softcore.f]', &
           'cant deallocate molecule_type','')

      if(allocated(molecule_type_partner)) deallocate (molecule_type_partner, stat=ier_alloc)
      if (ier_alloc /= 0) call sander_bomb('cleanup_sc[softcore.f]', &
           'cant deallocate molecule_type_partner','')

      RETURN
    end subroutine cleanup_sc

!===================================================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Compute and emit net charge, neutralizing that from roundoff error.
!-----------------------------------------------------------------------
!     --- TI_CHECK_NEUTRAL ---
!-----------------------------------------------------------------------
!
!     This routine will remove any net charge resulting from
!     conversion of the low precision parm topology charges.
!     It has been modified for TI, so that atoms that are common to V0/V1
!     have the same charges.
!
#ifdef MPI
subroutine ti_check_neutral(charge,natom)

   use constants, only : INV_AMBER_ELECTROSTATIC, TEN_TO_MINUS2, zero  
   implicit none
   _REAL_  charge(*)
   integer natom
   _REAL_             :: charge_partner(natom + extra_atoms)
   ! charge_lst(i): 0 - ith charge must be the same, 1 - charge can vary
   integer            :: charge_lst(natom)
   integer            :: i, j
   integer            :: nsum(3), nsum_partner
   integer            :: nsum_net(2)
   integer            :: nchange(2)
   integer            :: ichange
   _REAL_             :: sum_val(3), sum_val_partner
   _REAL_             :: sum_val_net(3)
   integer            :: myproc, myproc_partner
   integer            :: iidx, jidx
#  include "box.h"
#  include "ew_cntrl.h"
#  include "extra.h"
#  include "parallel.h"
   include 'mpif.h'

   nsum(:) = 0
   sum_val(:) = 0.d0
   nchange(:) = 0

   if (master) then 
     partner = ieor(masterrank,1)

     if (masterrank .eq. 1) then
       myproc = 2
       myproc_partner = 1
     else
       myproc = 1
       myproc_partner = 2
     end if

     ! Get the charge array from the partner
     call mpi_barrier( commmaster, ierr )
     call mpi_sendrecv( charge, natom, MPI_DOUBLE_PRECISION, partner, 5, &
        charge_partner, natom+extra_atoms, MPI_DOUBLE_PRECISION, partner, 5, &
        commmaster, ist, ierr )
     charge_lst(:) = 0
     if (ifsc .eq. 0) then
       do i = 1, natom
         if (abs((charge(i) - charge_partner(i))*INV_AMBER_ELECTROSTATIC) > 1.d-10) then
           charge_lst(i) = 1
         end if
       end do
     else
       do i = 1, nmixed
         iidx = (frcmask(3*(i-1)+1) + 2) / 3
         jidx = (frcmask_partner(3*(i-1)+1) + 2) / 3
         if (abs((charge(iidx) - charge_partner(jidx))*INV_AMBER_ELECTROSTATIC) > 1.d-10) then
           charge_lst(i) = 1
         end if
       end do
       charge_lst(:) = charge_lst(:) + nsc(:)
     end if

     do j = 1, natom
       if (charge_lst(j) .eq. 1) then
         sum_val(myproc) = sum_val(myproc) + charge(j)
         nsum(myproc) = nsum(myproc) + 1
       else
         sum_val(3) = sum_val(3) + charge(j)
         nsum(3) = nsum(3) + 1
       end if
     end do

     call mpi_sendrecv( sum_val(myproc), 1, MPI_DOUBLE_PRECISION, partner, 5, &
                        sum_val_partner, 1, MPI_DOUBLE_PRECISION, partner, 5, &
                        commmaster, ist, ierr )
     sum_val(myproc_partner) = sum_val_partner
     
     call mpi_sendrecv( nsum(myproc), 1, MPI_INTEGER, partner, 5, &
                        nsum_partner, 1, MPI_INTEGER, partner, 5, &
                        commmaster, ist, ierr )
     nsum(myproc_partner) = nsum_partner

     do i = 1, 2
       sum_val_net(i) = sum_val(i) + sum_val(3)
       nsum_net(i) = nsum(i) + nsum(3)
       if (master) write(6, '(/,5x,a,i2,a,f12.8)') &
         'Sum of charges for TI region ', i, ' = ', sum_val_net(i) * INV_AMBER_ELECTROSTATIC      

       if (abs(sum_val_net(i)*INV_AMBER_ELECTROSTATIC) .gt. TEN_TO_MINUS2) then
         if (master) &
           write(6, '(5x,a,/)') 'Assuming uniform neutralizing plasma'
       else
         if (master) &
           write(6, '(5x,a,/)') 'Forcing neutrality...'
           nchange(i) = 1
           ichange = i
       end if
     end do

     if (nchange(1)+nchange(2) .eq. 1) then
         if (nsum_net(ichange) .gt. 0) then
           sum_val_net(ichange) = sum_val_net(ichange) / nsum_net(ichange)
         else
           sum_val_net(ichange) = 0.d0
         end if

         do j = 1, natom
           if (ichange .eq. myproc .or. charge_lst(j) .eq. 0) then
             charge(j) = charge(j) - sum_val_net(ichange)
           end if
         end do 
     else if (nchange(1)+nchange(2) .eq. 2) then
       ! We need to remove the net charge on two effective simulation regions
       ! We first remove the average charge, then remove the remaining amount 
       ! from the softcore atoms. 
     
       ! Check for regions that have no charge
       if (nsum_net(1) .ne. 0 .and. nsum_net(2) .ne. 0) then
         sum_val_net(3) = 0.5d0 * (sum_val_net(1) / nsum_net(1) + &
                          sum_val_net(2) / nsum_net(2))
       else if (nsum_net(1) .ne. 0) then
         sum_val_net(3) = sum_val_net(1) / nsum_net(1)
       else if (nsum_net(2) .ne. 0) then
         sum_val_net(3) = sum_val_net(2) / nsum_net(2)
       else
         sum_val_net(3) = 0.d0
       end if

       ! For completely VDW runs/perturb to nothing runs these can be zero
       if (nsum(1) .gt. 0) then
         sum_val_net(1) = (sum_val_net(1) - sum_val_net(3) * nsum(3)) / nsum(1)
       else
         sum_val_net(1) = 0.d0
       end if

       if (nsum(2) .gt. 0) then
         sum_val_net(2) = (sum_val_net(2) - sum_val_net(3) * nsum(3)) / nsum(2)
       else
         sum_val_net(2) = 0.d0
       end if

       do j = 1, natom
         if (charge_lst(j) .eq. 0) then
           charge(j) = charge(j) - sum_val_net(3)
         else
           charge(j) = charge(j) - sum_val_net(myproc)
         end if
       end do         
     end if

   end if

   call mpi_bcast(charge,natom,MPI_DOUBLE_PRECISION,0,commsander,ierr)

   return
end subroutine ti_check_neutral 
#endif

end module softcore
                                                            ! === TODO === !

! softcore and 10-12 are not compatible

!===================================================================================================
