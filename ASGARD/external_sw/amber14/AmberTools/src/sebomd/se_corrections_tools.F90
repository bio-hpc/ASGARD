module se_corrections_tools
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!     Author   :: Antoine MARION 
!!     Date     :: 2013-12-10
!!     Function :: All needed tools to perform 
!!                 SE corrections.
!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 use sebomd_arrays, only : get_Z, get_Zeff, &
                           put_PEP, get_PEP, &
                           put_all_interactions, put_interaction, get_interaction
 use se_corrections_params, only : initPIFparams, checkPIF, getPIFparams, &
                                   initMAISparams, checkMAIS, getMAISparams
                            
 implicit none
 private
 
 double precision, dimension (4,83) :: agauss, bgauss, cgauss
 integer :: n_gauss
 !! Frc cst for MM peptidic correction
 double precision :: k_pep
 !! Cosntants
 double precision :: hartree2eV, Bohr2Ang, ev2kcal
  

 public :: LookFor_PEP
 public :: initPIF
 public :: initMAIS
 public :: computePIF
 public :: computeMAIS
 public :: computeGauss
 public :: GetRij
 public :: putConstants
 public :: compPEP

 contains
  
  subroutine putConstants(c1,c2,c3)
   !! Store the constants comming from sebomd
   !! to remain consistent with the rest of the code.
   !!
   !! hartree2eV, Bohr2Ang and ev2kcal
   implicit none
   double precision, intent(in) :: c1,c2,c3

   hartree2eV = c1
   Bohr2Ang   = c2
   ev2kcal    = c3

   return

  end subroutine putConstants
  
  subroutine computeGauss(i,j,rij,e_gauss,de_gauss)
   !! Compute g_pm3(A,B) or g_am1(A,B) term and
   !! first derivative for a couple of atoms A and B 
   !! g(A,B) = sum_{k=1}^{ng} (
   !!             a_{k,A}e^{-b_{k,A}(R_{AB}-c{{k,A})^2}
   !!           + a_{k,B}e^{-b_{k,B}(R_{AB}-c{{k,B})^2} )
   !!
   !! ng = 2,4 for AM1 and 2 for PM3
   implicit none
   integer, intent(in) :: i,j
   double precision, intent(in) :: rij
   double precision, intent(out) :: e_gauss, de_gauss
   ! local
   integer :: zi, zj
   double precision :: znuci, znucj
   double precision :: ZZ
   integer :: ng
   double precision :: agi, bgi, cgi, agj, bgj, cgj
   double precision :: argi, argj
   double precision :: r, r2

   e_gauss = 0.d0
   de_gauss = 0.d0

   
   call Get_Z(i,zi)
   call Get_Z(j,zj)
   call Get_Zeff(i,znuci)
   call Get_Zeff(j,znucj)

    
   do ng = 1,n_gauss
    agi = agauss(ng,zi)
    bgi = bgauss(ng,zi)
    cgi = cgauss(ng,zi)

    agj = agauss(ng,zj)
    bgj = bgauss(ng,zj)
    cgj = cgauss(ng,zj)

    argi = bgi*(rij-cgi)**2
    argj = bgj*(rij-cgj)**2

    if (argi.lt.25.0d0) then
     e_gauss = e_gauss + agi*exp(-argi)
     de_gauss = de_gauss - agi*2.0d0*bgi*(rij-cgi)*exp(-argi)
    endif

    if (argj.lt.25.0d0) then
     e_gauss = e_gauss + agj*exp(-argj)
     de_gauss = de_gauss - agj*2.0d0*bgj*(rij-cgj)*exp(-argj)
    endif
   enddo

   ZZ = znuci*znucj

   r = 1.0d0/rij
   r2 = r*r

   ! Gradient in eV.Ang-1
   de_gauss = (de_gauss*ZZ*r - e_gauss*ZZ*r2)
   ! Gradient in kcal.mol-1.Ang-1
   de_gauss = de_gauss * ev2kcal

   ! Energy in eV
   e_gauss = e_gauss*ZZ*r

   return

  end subroutine computeGauss

  subroutine computeMAIS(i,j,rij,int_type,emais,demais)
   !! Compute g_mais(A,B) term and first derivative 
   !! for a couple of atoms A and B.
   !! g_{mais}(A,B) =  sum_{k=1}^{3} ( 
   !!    \alpha_{k,AB} e^{-\beta_{k,AB}(\gamma_{k,AB} - R_{AB})^2 )
   !!
   !! refs ::
   !!   - MAIS1 (Water):
   !!     *M.I. Bernal-Uruchurtu, M.F. Ruiz-Lopez Chem. Phys. Lett. 330 (2000) 118-124
   !!   - MAIS2 (H-Cl in water):
   !!     *O. Arillo-Flores, M.F. Ruiz-Lopez, M.I. Bernal-Uruchurtu Theor. Chem. Acc. 118 (2007) 425-435
   implicit none
   integer, intent(in) :: i,j,int_type
   double precision, intent(in) :: rij
   double precision, intent(out) :: emais, demais
   ! local
   integer :: zi, zj, n
   double precision :: r, gau
   double precision, dimension(3) :: a,b,c

   call get_Z(i,zi)
   call get_Z(j,zj)

   call getMAISparams(zi,zj,int_type,a,b,c)

   r = rij / Bohr2Ang
   
   emais = 0.d0
   demais = 0.d0

   do n=1,3

    gau = a(n)*dexp(-b(n)*(c(n)-r)**2)

    emais = emais + gau
    demais = demais + 2.0d0*b(n)*(c(n)-r)*gau

   enddo

   ! Energy in eV
   emais = emais * Hartree2eV
   ! Gradient in kcal.mol-1.Ang-1
   demais = demais * Hartree2eV * ev2kcal / Bohr2Ang

   return

  end subroutine computeMAIS

  subroutine computePIF(i,j,rij,int_type,epif,depif)
   !! Compute g_pif(A,B) term and first derivative 
   !! for a couple of atoms A and B located on two different molecules.
   !! g_{pif}(A,B) =   \alpha_{AB} e^{-\beta_{AB}R_{AB} 
   !!               +  \chi_{AB}     / R_{AB}^6
   !!               +  \delta_{AB}   / R_{AB}^8
   !!               +  \epsilon_{AB} / R_{AB}^10
   !!
   !! refs ::
   !!   - PIF1 (water-water):
   !!     *M.I. Bernal-Uruchurtu, M.T.C. Matins-Costa, C. Millot, M.F. Ruiz-Lopez J. Comput. Chem. 21(2000) 572-581
   !!   - PIF2 (solute-water + HCl-water):
   !!     *W. Harb, M.I. Bernal-Uruchurtu, M.F. Ruiz-Lopez Theor. Chem. Acc. 112 (2004) 204-216
   !!     *O. Arillo-Flores, M.F. Ruiz-Lopez, M.I. Bernal-Uruchurtu Theor. Chem. Acc. 118 (2007) 425-435
   !!   - PIF3 (Htype dependent PIF):
   !!     *A. Marion, G. Monard, M.F. Ruiz-Lopez, F. Ingrosso ????????
   implicit none
   integer, intent(in) :: i,j,int_type
   double precision, intent(in) :: rij
   double precision, intent(out) :: epif, depif
   ! local
   integer :: zi, zj
   double precision :: r
   double precision :: alpha, beta, chi, delta, epsilon
   double precision :: r1,r2,r6,r8,r10
   double precision :: e1,e2,e3,e4,de1,de2,de3,de4


   call get_Z(i,zi)
   call get_Z(j,zj)

   call getPIFparams(zi,zj,int_type,alpha,beta,chi,delta,epsilon) 

   r = rij / Bohr2Ang

   r1 = 1.0d0/r
   r2 = r1*r1
   r6 = r2*r2*r2
   r8 = r6*r2
   r10 = r8*r2

   e1 = alpha*dexp(-beta*r)
   e2 = chi*r6
   e3 = delta*r8
   e4 = epsilon*r10

   epif = e1 + e2 + e3 + e4

   de1 =   -beta*e1
   de2 =  -6.0d0*e2
   de3 =  -8.0d0*e3
   de4 = -10.0d0*e4


   depif = de1 + (de2 + de3 + de4)*r1


   ! Energy in eV
   epif = epif * Hartree2eV
   ! Gradient in kcal.mol-1.Ang-1
   depif = depif * Hartree2eV * ev2kcal / Bohr2Ang

   return

  end subroutine computePIF
  
  subroutine LookFor_PEP(natoms,xyz,kp,if_pbc,symbol2,Npep)
  !! Find amide groups
  !!  X1    X2
  !!   \   /
  !!     N
  !!     |
  !!     C
  !!   //  \
  !!  O     R
  !!
  !! 1 :: Look for N type atoms
  !!      (in Amber, 'N ' type is for nitrogen of amide)
  !!      (in Gaff, 'n ' type is for nitrogen of amide)
  !! 2 :: Look for the atoms connected to this nitrogen (using cutoff=1.5)
  !!      if the atom is C/c type (carbonyl in Amber) then store the C atoms
  !!      else store X1 and X2 atoms
  !! 3 :: Look for atoms connected to the C atom (cutoff = 1.5)
  !!      if this atom is type O/o (carbonyl in Amber/Gaff) store Oatm
  !!
  !! Npep_bonds = number of peptidic bonds found
  !! PEPlist(ipep,k) = contains all the dihedrals on which we
  !!                   should add MM potential
  !!                   ipep = index of the dihedral
  !!                   k = 1 :: X(1,2)
  !!                     = 2 :: N
  !!                     = 3 :: C
  !!                     = 4 :: O 
    implicit none
    integer, intent(in) :: natoms
    character (len=4), intent(in) :: symbol2(natoms)
    double precision, intent(in) :: xyz(3,natoms)
    double precision, intent(in) :: kp
    logical, intent(in) :: if_pbc
    integer, intent(out) :: Npep
    ! local
    integer :: i,j,k
    integer :: zj
    double precision :: rij, uij(3)
    double precision :: cutoff, cutoff_NH
    integer :: Xatms(2)
    integer :: Natm
    integer :: Catm
    integer :: Oatm
    integer :: cnt_X

    cutoff = 1.8d0 ! CutOff for covalent bonds
    cutoff_NH = 1.3d0 ! CutOff for NH covalent bonds

    k_pep = kp / ev2kcal

    Npep = 0

    do i=1,natoms
     Natm = 0
     Catm = 0
     Oatm = 0
     Xatms(:) = 0

     if (symbol2(i).eq.'N'.or.symbol2(i).eq.'n') then
      Natm = i
      cnt_X = 0
      do j=1,natoms
       if ((j.ne.i) & 
        .and.(symbol2(j).ne.'OW') &
        .and.(symbol2(j).ne.'HW')) then ! No need to check for water atoms

        call GetRij(i,j,xyz,if_pbc,rij,uij)
        call Get_Z(j,zj)

        if (rij.le.cutoff) then ! Check only the atoms close to N
         if (.not.(symbol2(j).eq.'C'.or.symbol2(j).eq.'c')) then ! Means we found X1 or X2
          if (zj.eq.1) then ! means we found an Hydrogen atom
           if (rij.le.cutoff_NH) then ! check if this H is bonded to N or just H-bonded
             cnt_X = cnt_X + 1
             Xatms(cnt_X) = j
           endif
          else
           cnt_X = cnt_X + 1
           Xatms(cnt_X) = j
          endif
         else ! Means we found C
          Catm = j
          do k = 1,natoms

           if ((k.ne.j).and.(symbol2(k).eq.'O'.or.symbol2(k).eq.'o')) then ! Look only for O type atoms
            call GetRij(j,k,xyz,if_pbc,rij,uij)
            if (rij.le.cutoff) Oatm = k
           endif ! O

          enddo ! k
         endif ! X or C
        endif ! cutoff
       endif ! j .not. water atom
      enddo ! j

      if    ((Natm.ne.0) & ! Check if the pep bond is complete
        .and.(Catm.ne.0) &
        .and.(Oatm.ne.0) &
        .and.(Xatms(1).ne.0) &
        .and.(Xatms(2).ne.0)) then
       Npep = Npep + 1
        call put_PEP(Npep,1,Xatms(1))
        call put_PEP(Npep,2,Natm)
        call put_PEP(Npep,3,Catm)
        call put_PEP(Npep,4,Oatm)
       Npep = Npep + 1
        call put_PEP(Npep,1,Xatms(2))
        call put_PEP(Npep,2,Natm)
        call put_PEP(Npep,3,Catm)
        call put_PEP(Npep,4,Oatm)
      endif ! if complete pep bond

     endif ! if N
    enddo ! i

    return
  end subroutine LookFor_PEP 

  subroutine initMAIS(natoms,ctype,let_pm3,a,b,c,symbol2,ierror)
    !! init MAIS interactions types (int_type) and check if the 
    !! parameters are all available.
    !! If they are missing parameters, let_pm3 controls
    !! rather the calculation continues or not
    !!
    !! int_type = 0 :: PM3
    !!            1 :: MAIS1
    !!            2 :: MAIS2
    !!
    !! Interaction type tells which parameters to use (MAIS 1 or 2)
    !! 
    !! refs ::
    !!   - MAIS1 (Water):
    !!     *M.I. Bernal-Uruchurtu, M.F. Ruiz-Lopez Chem. Phys. Lett. 330 (2000) 118-124
    !!   - MAIS2 (H-Cl in water):
    !!     *O. Arillo-Flores, M.F. Ruiz-Lopez, M.I. Bernal-Uruchurtu Theor. Chem. Acc. 118 (2007) 425-435
    implicit none
    integer, intent(in) :: natoms, ctype
    character (len=4), intent(in) :: symbol2(natoms)
    integer, intent(inout) :: ierror
    double precision, dimension (:,:), intent(in) :: a,b,c
    logical, intent(in) :: let_pm3
    ! local
    integer :: i,j,int_type,maistype
    integer :: zi, zj
    logical :: if_param


    if (ctype == 3) then
     maistype = 1
    elseif (ctype == 4) then
     maistype = 2
    endif

    agauss(:,:) = a(:,:)
    bgauss(:,:) = b(:,:)
    cgauss(:,:) = c(:,:)

    call put_all_interactions(0)

    n_gauss = 2

    call initMAISparams

    do i=1,natoms-1
     call get_Z(i,zi)
     do j=i+1,natoms
      call get_Z(j,zj)

      Call checkMAIS(zi,zj,maistype,if_param)

      if ((.not.if_param).and.(let_pm3)) then
       int_type = 0
       write(6,'("PM3 parameters used for interaction ", a4,"-",a4)') symbol2(i),symbol2(j)
      elseif (if_param) then
       int_type = maistype
      else
       ierror = 1
       write(6,'("Error :: No MAIS parameters for interaction ", a4,"-",a4)') symbol2(i),symbol2(j)
       write(6,'("         Please set let_pm3 = 1 to use PM3 parameters")')
       exit
      endif

      call put_interaction(i,j,int_type)

     enddo
     if (ierror.ne.0) exit
    enddo

    return

  end subroutine initMAIS

  subroutine initPIF(natoms,ctype,let_pm3,symbol2,a,b,c,ierror)
    !! init PIF interactions types (int_type) and check if the 
    !! parameters are all available.
    !! If they are missing parameters, let_pm3 controls
    !! rather the calculation continues or not.
    !!
    !! int_type = 0 :: PM3
    !!            1 :: PIF1
    !!            2 :: PIF2
    !!            3 :: PIF3
    !! 
    !! Interaction type tells which parameters to use (PIF 1,2 or 3)
    !!
    !! PIF1 for water-water interactions
    !! PIF2 for solute water interactions
    !! PIF3 for Htype dependent solute-water interactions  
    !!
    !! refs ::
    !!   - PIF1 (water-water):
    !!     *M.I. Bernal-Uruchurtu, M.T.C. Matins-Costa, C. Millot, M.F. Ruiz-Lopez J. Comput. Chem. 21(2000) 572-581
    !!   - PIF2 (solute-water + HCl-water):
    !!     *W. Harb, M.I. Bernal-Uruchurtu, M.F. Ruiz-Lopez Theor. Chem. Acc. 112 (2004) 204-216
    !!     *O. Arillo-Flores, M.F. Ruiz-Lopez, M.I. Bernal-Uruchurtu Theor. Chem. Acc. 118 (2007) 425-435
    !!   - PIF3 (Htype dependent PIF):
    !!     *A. Marion, G. Monard, M.F. Ruiz-Lopez, F. Ingrosso ????????
    use se_inter, only : se_interres
    implicit none
    integer, intent(in) :: natoms, ctype
    character (len=4), intent(in) :: symbol2(natoms)
    integer, intent(inout) :: ierror
    double precision, dimension (:,:), intent(in) :: a,b,c
    logical, intent(in) :: let_pm3
    ! local
    integer :: i,j,int_type,resinter, cntH
    integer :: zi, zj, pif1, pif2
    character (len=4) :: sym_i, sym_j
    logical :: H_hydrophobic(natoms)
    logical :: if_param
    logical :: PIF3

    agauss(:,:) = a(:,:)
    bgauss(:,:) = b(:,:)
    cgauss(:,:) = c(:,:)

    n_gauss = 2

    pif1 = 1
    pif2 = 2
    cntH = 0

    call initPIFparams

    call put_all_interactions(0)

    if (ctype.eq.2) then
     PIF3 = .true.
    else
     PIF3 = .false.
    endif

!!!! 2013-12-11
!!!! List of H types in Amber
!!!! Source http://ambermd.org/doc/AtomTypesTableWorkspace.xhtml
!------------------------------------------------------
! Symbol |  Description           | Hydrophobic? (PIF3)
!------------------------------------------------------
!  H     |  H bonded to nitrogen  |      false
!
!  hn    |  Idem in gaff          |      false
!
!  hN    |  Idem in lipid11       |      false
!
!  HW    |  H in water            |      false
!
!  hw    |  Idem in gaff          |      false
!
!  HO    |  H bonded to oxygen    |      false
!
!  Ho    |  Idem in GLYCAM        |      false
!
!  hO    |  Idem in lipid11       |      false
!
!  hR    |  Idem in lipid11??     |      false
!
!  ho    |  Idem in gaff          |      false
!
!  HS    |  H bonded to sulfure   |      false
!
!  hs    |  Idem in gaff          |      false
!
!  hp    |  H bonded to phostphate|      false
!
!  H1    |  Aliphatic H bonded    |      true 
!        |  to C with 1 electron  |      
!        |  withdrawing group     |      
!
!  h1    |  Idem in gaff          |      true 
!
!  hE    |  Idem in lipid11       |      true 
!
!  hS    |  Idem in lipid11??     |      false (not sure)
!
!  H2    |  Aliphatic H bonded    |      true 
!        |  to C with 2 electron  |      
!        |  withdrawing group     |      
!
!  h2    |  Idem in gaff          |      true 
!
!  H3    |  Aliphatic H bonded    |      true 
!        |  to C with 3 electron  |      
!        |  withdrawing group     |      
!
!  h3    |  Idem in gaff          |      true 
!
!  HA    |  H aromatic bonded to  |      true 
!        |  C without electron    |      
!        |  withdrawing group     |      
!
!  ha    |  Idem in gaff          |      true 
!
!  Ha    |  Idem in GLYCAM        |      true 
!
!  hB    |  Idem in lipid11       |      true 
!
!  H4    |  H aromatic bonded to  |      true 
!        |  C with 1 electron     |      
!        |  withdrawing group     |      
!
!  h4    |  Idem in gaff          |      true 
!
!  H5    |  H aromatic bonded to  |      true 
!        |  C with 2 electron     |      
!        |  withdrawing group     |      
!        |  + HCOO group          |      
!
!  h5    |  Idem in gaff          |      true 
!
!  HC    |  H coneted to          |      true 
!        |  aliphatic C           |           
!
!  hc    |  H coneted to          |      true 
!        |  aliphatic C (gaff)    |           
!
!  Hc    |  H coneted to          |      true 
!        |  aliphatic C (GLYCAM)  |           
!
!  hA    |  H coneted to          |      true 
!        |  aliphatic C (lipid11) |           
!
!  HP    |  H coneted to C next   |      true 
!        |  to positively         | 
!        |  charged group         | 
!
!  Hp    |  Idem in GLYCAM        |      true 
!
!  hX    |  Idem in lipid11       |      true 
!
!  hx    |  Idem in gaff          |      true 
!
!  HZ    |  H coneted to sp C     |      true 
!--------------------------------------------

    ! Check if we have hydrophobic H for PIF3
    if (PIF3) then
     do i=1,natoms
      sym_i = symbol2(i)
      if     ((sym_i.eq.'HC') &
        .or.  (sym_i.eq.'hc') &
        .or.  (sym_i.eq.'Hc') &
        .or.  (sym_i.eq.'hA') &
        .or.  (sym_i.eq.'HP') &
        .or.  (sym_i.eq.'Hp') &
        .or.  (sym_i.eq.'hX') &
        .or.  (sym_i.eq.'hx') &
        .or.  (sym_i.eq.'HZ') &
        .or.  (sym_i.eq.'HA') &
        .or.  (sym_i.eq.'ha') &
        .or.  (sym_i.eq.'Ha') &
        .or.  (sym_i.eq.'hB') &
        .or.  (sym_i.eq.'H1') &
        .or.  (sym_i.eq.'h1') &
        .or.  (sym_i.eq.'hE') &
        .or.  (sym_i.eq.'H2') &
        .or.  (sym_i.eq.'h2') &
        .or.  (sym_i.eq.'H3') &
        .or.  (sym_i.eq.'h3') &
        .or.  (sym_i.eq.'H4') &
        .or.  (sym_i.eq.'h4') &
        .or.  (sym_i.eq.'H5') &
        .or.  (sym_i.eq.'h5')) then
        H_hydrophobic(i) = .true.
        cntH = cntH + 1
      else
        H_hydrophobic(i) = .false.
      endif
     enddo
     write(6,'("")')
     write(6,'(" -------------------------------------------------------------- ")')
     write(6,'("      Applying PIF3 intermolecular potential ")')
     write(6,'("")')
     write(6,'(" Number of H considered as hydrophobic for PIF3:", i12)') cntH
     write(6,'(" -------------------------------------------------------------- ")')
    endif

       

    do i=1,natoms-1
     call get_Z(i,zi)
     sym_i = symbol2(i)
     do j=i+1,natoms
      call get_Z(j,zj)
      sym_j = symbol2(j)

      call se_interres(i,j,resinter)

      if (resinter.eq.0) then ! intramolecular interaction (keep PM3)

       int_type = 0

      elseif (resinter.eq.1) then ! Water-Water interaction

       Call checkPIF(zi,zj,pif1,if_param)
       if (if_param) then
        int_type = 1
       else
        write(6,'("Error while checking parameters for Water-Water interaction")')
        ierror = 1
        exit
       endif

      elseif (resinter.eq.2) then ! Solute-Water interactions

       if ((PIF3) .and. (H_hydrophobic(i).or.H_hydrophobic(j))) then  ! PIF3
        int_type = 3
       else ! Standard PIF2
        int_type = 2
        Call checkPIF(zi,zj,pif2,if_param)
        if (.not.if_param) then
         if (let_pm3) then
          int_type = 0
         else
          write(6,'("Interaction ", a4,"-",a4," is not defined with PIF")') symbol2(i),symbol2(j)
          write(6,'("Use let_pm3 = 1 to use PM3 for missing interactions")')
          ierror = 1
          exit
         endif
        endif

       endif

      endif ! resinter

      if (ierror.ne.0) exit

      call put_interaction(i,j,int_type)

     enddo ! j
     if (ierror.ne.0) exit
    enddo ! i
    
    return

  end subroutine initPIF

  subroutine GetRij(i,j,xyz,if_pbc,rij,uij)
    !! Returns Rij and the unitary vector (Uij)
    !! depending on if pbc is turned on or off
    implicit none
    integer, intent(in) :: i,j
    double precision, intent(in) :: xyz(3,*)
    logical, intent(in) :: if_pbc
    double precision, intent(out) :: rij
    double precision, intent(out) :: uij(3)
    ! local
    integer :: k
    double precision :: xi,yi,zi
    double precision :: xj,yj,zj
     
    xi = XYZ(1,i)
    yi = XYZ(2,i)
    zi = XYZ(3,i)


    if (if_pbc) then
      call se_pbcxyz(i,j,xj,yj,zj)
    else
     xj = XYZ(1,j)
     yj = XYZ(2,j)
     zj = XYZ(3,j)
    endif

    uij(1) = xi - xj
    uij(2) = yi - yj
    uij(3) = zi - zj

    rij = 0.d0
    do k=1,3
     rij = rij + uij(k)**2
    enddo

    rij = dsqrt(rij)

    uij(:) = uij(:)/rij
    
    return

  end subroutine GetRij

  subroutine compPEP(ipep,pbc,xyz,e_pep,grad)
   !! Computes peptidic correction for all the identified dihedrals.
   !! Returns energy and gradient (numerical)
   !!
   !! E_{pep} = k_{pep} * sin(\phi)^2
   implicit none
   integer, intent(in) :: ipep
   logical, intent(in) :: pbc
   double precision, intent(inout) :: e_pep
   double precision, intent(inout) :: grad(3,*)
   double precision, intent(inout) :: xyz(3,*)
   ! local
   integer :: j, k
   double precision :: dih
   integer, dimension(4) :: pep_b
   double precision :: delX, eb, ef

   delX = 1.0d-8

   call get_PEP(ipep,pep_b)
   call se_dihedrpbc(pbc,xyz,pep_b(1), &
                             pep_b(2), &
                             pep_b(3), &
                             pep_b(4), &
                             dih)
   e_pep = e_pep + k_pep*dsin(dih)**2

   do j = 1,4
    do k=1,3
     xyz(k,pep_b(j)) = &
     xyz(k,pep_b(j)) - delX

     call se_dihedrpbc(pbc,xyz,pep_b(1), &
                               pep_b(2), &
                               pep_b(3), &
                               pep_b(4), &
                               dih)

     eb = k_pep*dsin(dih)**2

     xyz(k,pep_b(j)) = &
     xyz(k,pep_b(j)) + 2.d0*delX

     call se_dihedrpbc(pbc,xyz,pep_b(1), &
                               pep_b(2), &
                               pep_b(3), &
                               pep_b(4), &
                               dih)

     ef = k_pep*dsin(dih)**2

     xyz(k,pep_b(j)) = &
     xyz(k,pep_b(j)) - delX

     grad(k,pep_b(j)) = &
     grad(k,pep_b(j)) - &
     (eb -  ef)/(2.d0*delX)*eV2kcal
    enddo
   enddo

   return

  end subroutine compPEP


end module se_corrections_tools

