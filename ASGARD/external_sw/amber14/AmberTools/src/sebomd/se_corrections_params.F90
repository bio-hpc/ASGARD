module se_corrections_params
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!     Author   :: Antoine MARION 
!!     Date     :: 2013-12-10
!!     Function :: Stores and outputs SemiEmpirical parameters
!!                 - PIF1
!!                   *M.I. Bernal-Uruchurtu, M.T.C. Matins-Costa, 
!!                    C. Millot, M.F. Ruiz-Lopez J. Comput. Chem.
!!                    21(2000) 572-581
!!                 - PIF2
!!                   *W. Harb, M.I. Bernal-Uruchurtu, 
!!                    M.F. Ruiz-Lopez Theor. Chem. Acc. 
!!                    112 (2004) 204-216
!!                   *O. Arillo-Flores, M.F. Ruiz-Lopez,
!!                    M.I. Bernal-Uruchurtu Theor. Chem. Acc. 
!!                    118 (2007) 425-435 
!!                 - PIF3
!!                   *A. Marion, G. Monard, M.F. Ruiz-Lopez,
!!                    F. Ingrosso  ????????
!!                 - MAIS1
!!                   *M.I. Bernal-Uruchurtu, M.F. Ruiz-Lopez 
!!                    Chem. Phys. Lett. 330 (2000) 118-124
!!                 - MAIS2
!!                   *O. Arillo-Flores, M.F. Ruiz-Lopez,
!!                    M.I. Bernal-Uruchurtu Theor. Chem. Acc. 
!!                    118 (2007) 425-435 
!!     Details  :: Parameters are stored in small tables.
!!                 Pointers are used to make correspondance
!!                 between atomic number and position in 
!!                 the table.
!!      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 implicit none
 private

 !!!! Pointers for parameters matrices
 !! dimension is 83 to fit with PM3 
 integer, dimension(83) :: pointersPIF1, pointersPIF2, pointersPIF3
 integer, dimension(83) :: pointersMAIS1, pointersMAIS2
 !!!! PIF1 parameters
 !! Water-Water interaction
 integer, parameter :: nPIF1atm = 2
 double precision, dimension(nPIF1atm,nPIF1atm,5) :: PIF1
 logical, dimension(nPIF1atm,nPIF1atm) :: PIF1test
 !!!! PIF2 parameters
 !! Solute-Water interaction
 integer, parameter :: nPIF2atm = 5
 double precision, dimension(nPIF2atm,nPIF2atm,5) :: PIF2
 logical, dimension(nPIF2atm,nPIF2atm) :: PIF2test
 !!!! PIF3 parameters
 !! Solute-Water interaction (Hydrophobic)
 integer, parameter :: nPIF3atm = 2
 double precision, dimension(nPIF3atm,nPIF3atm,5) :: PIF3
 logical, dimension(nPIF3atm,nPIF3atm) :: PIF3test
 !!!! MAIS1 parameters
 !! Water-Water
 integer, parameter :: nMAIS1atm = 2
 double precision, dimension(nMAIS1atm,nMAIS1atm,9) :: MAIS1
 logical, dimension(nMAIS1atm,nMAIS1atm) :: MAIS1test
 !!!! MAIS2 parameters
 !! Water-Water + HCl
 integer, parameter :: nMAIS2atm = 3
 double precision, dimension(nMAIS2atm,nMAIS2atm,9) :: MAIS2
 logical, dimension(nMAIS2atm,nMAIS2atm) :: MAIS2test
  

 !public :: GetPIF
 public :: CheckPIF
 public :: CheckMAIS
 public :: initPIFparams
 public :: initMAISparams
 public :: getPIFparams
 public :: getMAISparams

 contains

 subroutine getMAISparams(zi,zj,int_type,a,b,c)
  !! Returns Mais parameters (a,b and c)
  !! for a couple of atoms, depending on the type
  !! of interaction selected (MAIS1 or MAIS2)
  !!
  !! a,b and c are dimension 3 vetors 
  implicit none
  integer, intent(in) :: zi, zj, int_type
  double precision, intent(out), dimension(3) :: a,b,c
  ! local
  integer :: p1, p2

  if (int_type.eq.1) then

   p1 = pointersMAIS1(zi)
   p2 = pointersMAIS1(zj)
   a(1) = MAIS1(p1,p2,1)
   a(2) = MAIS1(p1,p2,2)
   a(3) = MAIS1(p1,p2,3)
   b(1) = MAIS1(p1,p2,4)
   b(2) = MAIS1(p1,p2,5)
   b(3) = MAIS1(p1,p2,6)
   c(1) = MAIS1(p1,p2,7)
   c(2) = MAIS1(p1,p2,8)
   c(3) = MAIS1(p1,p2,9)

  elseif (int_type.eq.1) then

   p1 = pointersMAIS2(zi)
   p2 = pointersMAIS2(zj)
   a(1) = MAIS2(p1,p2,1)
   a(2) = MAIS2(p1,p2,2)
   a(3) = MAIS2(p1,p2,3)
   b(1) = MAIS2(p1,p2,4)
   b(2) = MAIS2(p1,p2,5)
   b(3) = MAIS2(p1,p2,6)
   c(1) = MAIS2(p1,p2,7)
   c(2) = MAIS2(p1,p2,8)
   c(3) = MAIS2(p1,p2,9)
   
  endif

  return

 end subroutine getMAISparams

 subroutine getPIFparams(zi,zj,int_type,a,b,c,d,e)
  !! Returns PIF parameters (a,b and c)
  !! for a couple of atoms, depending on the type
  !! of interaction selected (PIF1,2 or 3)
  implicit none
  integer, intent(in) :: zi, zj, int_type
  double precision, intent(out) :: a,b,c,d,e
  ! local
  integer :: p1, p2

  if (int_type.eq.1) then

   p1 = pointersPIF1(zi)
   p2 = pointersPIF1(zj)
   a = PIF1(p1,p2,1)
   b = PIF1(p1,p2,2)
   c = PIF1(p1,p2,3)
   d = PIF1(p1,p2,4)
   e = PIF1(p1,p2,5)

  elseif (int_type.eq.2) then

   p1 = pointersPIF2(zi)
   p2 = pointersPIF2(zj)
   a = PIF2(p1,p2,1)
   b = PIF2(p1,p2,2)
   c = PIF2(p1,p2,3)
   d = PIF2(p1,p2,4)
   e = PIF2(p1,p2,5)
   
  elseif (int_type.eq.3) then

   p1 = pointersPIF3(zi)
   p2 = pointersPIF3(zj)
   a = PIF3(p1,p2,1)
   b = PIF3(p1,p2,2)
   c = PIF3(p1,p2,3)
   d = PIF3(p1,p2,4)
   e = PIF3(p1,p2,5)
   
  endif

  return

 end subroutine getPIFparams

 subroutine initPIFparams()
  implicit none

  call initPIF1
  call initPIF2
  call initPIF3

  return
 
 end subroutine initPIFparams

 subroutine initMAISparams()
  implicit none

  call initMAIS1
  call initMAIS2

  return
 
 end subroutine initMAISparams

 subroutine checkMAIS(zi,zj,int_type,if_param)
  implicit none
  integer, intent(in) :: zi, zj, int_type
  logical, intent(out) :: if_param

  ! local
  integer :: p1,p2

  if_param = .false.

  if (int_type.eq.1) then
   p1 = pointersMAIS1(zi)
   p2 = pointersMAIS1(zj)
   if ((p1.ne.0) .and. (p2.ne.0)) if_param = MAIS1test(p1,p2)
  elseif (int_type.eq.2) then
   p1 = pointersMAIS2(zi)
   p2 = pointersMAIS2(zj)
   if ((p1.ne.0) .and. (p2.ne.0)) if_param = MAIS2test(p1,p2)
  else
   print*, 'Problem checkMAIS'
  endif

  return

 end subroutine checkMAIS

 subroutine checkPIF(zi,zj,int_type,if_param)
  implicit none
  integer, intent(in) :: zi, zj, int_type
  logical, intent(out) :: if_param

  ! local
  integer :: p1,p2

  if_param = .false.

  if (int_type.eq.1) then
   p1 = pointersPIF1(zi)
   p2 = pointersPIF1(zj)
   if ((p1.ne.0) .and. (p2.ne.0)) if_param = PIF1test(p1,p2)
  elseif (int_type.eq.2) then
   p1 = pointersPIF2(zi)
   p2 = pointersPIF2(zj)
   if ((p1.ne.0) .and. (p2.ne.0)) if_param = PIF2test(p1,p2)
  elseif (int_type.eq.3) then
   p1 = pointersPIF3(zi)
   p2 = pointersPIF3(zj)
   if ((p1.ne.0) .and. (p2.ne.0)) if_param = PIF3test(p1,p2)
  else
   print*, 'Problem checkPIF'
  endif

  return

 end subroutine checkPIF

 subroutine initMAIS1()
  implicit none
  integer :: i,j

  MAIS1(:,:,:) = 0.d0
  MAIS1test(:,:) = .false.
  pointersMAIS1(:) = 0


  !! Atom pointers

  pointersMAIS1(1) = 1 ! H
  pointersMAIS1(8) = 2 ! O

  !! H-H
  MAIS1test(1,1) = .true.
  MAIS1(1,1,1) =  0.000362d0  ! alpha_1 
  MAIS1(1,1,2) =  0.009138d0  ! alpha_2
  MAIS1(1,1,3) =  0.007175d0  ! alpha_3
  MAIS1(1,1,4) =  0.385490d0  ! beta_1
  MAIS1(1,1,5) =  0.227530d0  ! beta_2
  MAIS1(1,1,6) =  3.013020d0  ! beta_3
  MAIS1(1,1,7) =  6.230900d0  ! gamma_1
  MAIS1(1,1,8) =  2.211680d0  ! gamma_2
  MAIS1(1,1,9) =  2.225720d0  ! gamma_3
   
  !! H-O
  MAIS1test(1,2) = .true.
  MAIS1(1,2,1) = -0.071576d0
  MAIS1(1,2,2) =  0.015255d0
  MAIS1(1,2,3) =  0.029575d0
  MAIS1(1,2,4) =  0.192055d0
  MAIS1(1,2,5) =  1.284460d0
  MAIS1(1,2,6) =  0.345024d0
  MAIS1(1,2,7) =  1.441740d0
  MAIS1(1,2,8) =  2.193810d0
  MAIS1(1,2,9) =  2.586020d0

  !! O-O
  MAIS1test(2,2) = .true.
  MAIS1(2,2,1) = -6.797745d0
  MAIS1(2,2,2) =  6.912401d0
  MAIS1(2,2,3) =  0.074600d0
  MAIS1(2,2,4) =  0.326159d0
  MAIS1(2,2,5) =  0.320157d0
  MAIS1(2,2,6) =  1.268150d0
  MAIS1(2,2,7) =  1.211370d0
  MAIS1(2,2,8) =  1.209450d0
  MAIS1(2,2,9) =  2.579530d0

  ! Store the upper triangle of the parameters matrix

  do j=1,nMAIS1atm
   do i=j+1,nMAIS1atm
    MAIS1test(i,j) = MAIS1test(j,i)
    MAIS1(i,j,:) = MAIS1(j,i,:)
   enddo
  enddo

 return
 
 end subroutine initMAIS1

 subroutine initMAIS2()
  implicit none
  integer :: i,j

  MAIS2(:,:,:) = 0.d0
  MAIS2test(:,:) = .false.
  pointersMAIS2(:) = 0


  !! Atom pointers

  pointersMAIS2(1)  = 1 ! H
  pointersMAIS2(8)  = 2 ! O
  pointersMAIS2(17) = 3 ! Cl

  !! H-H
  MAIS2test(1,1) = .true.
  MAIS2(1,1,1) =  0.000362d0  ! alpha_1 
  MAIS2(1,1,2) =  0.009138d0  ! alpha_2
  MAIS2(1,1,3) =  0.007175d0  ! alpha_3
  MAIS2(1,1,4) =  0.385490d0  ! beta_1
  MAIS2(1,1,5) =  0.227530d0  ! beta_2
  MAIS2(1,1,6) =  3.013020d0  ! beta_3
  MAIS2(1,1,7) =  6.230900d0  ! gamma_1
  MAIS2(1,1,8) =  2.211680d0  ! gamma_2
  MAIS2(1,1,9) =  2.225720d0  ! gamma_3
   
  !! H-O
  MAIS2test(1,2) = .true.
  MAIS2(1,2,1) =  0.038146d0
  MAIS2(1,2,2) = -0.042715d0 
  MAIS2(1,2,3) = -0.039872d0 
  MAIS2(1,2,4) =  0.564949d0
  MAIS2(1,2,5) =  0.518414d0
  MAIS2(1,2,6) =  0.342029d0
  MAIS2(1,2,7) =  4.128030d0
  MAIS2(1,2,8) =  4.058190d0
  MAIS2(1,2,9) =  0.529303d0

  !! H-Cl
  MAIS2test(1,3) = .true.
  MAIS2(1,3,1) = -0.247371d0
  MAIS2(1,3,2) =  0.245420d0
  MAIS2(1,3,3) = -0.009203d0 
  MAIS2(1,3,4) =  1.520440d0
  MAIS2(1,3,5) =  1.649960d0
  MAIS2(1,3,6) =  0.121566d0
  MAIS2(1,3,7) =  2.176850d0
  MAIS2(1,3,8) =  2.223590d0
  MAIS2(1,3,9) =  0.841882d0

  !! O-O
  MAIS2test(2,2) = .true.
  MAIS2(2,2,1) = -6.797745d0
  MAIS2(2,2,2) =  6.912401d0
  MAIS2(2,2,3) =  0.074600d0
  MAIS2(2,2,4) =  0.326159d0
  MAIS2(2,2,5) =  0.320157d0
  MAIS2(2,2,6) =  1.268150d0
  MAIS2(2,2,7) =  1.211370d0
  MAIS2(2,2,8) =  1.209450d0
  MAIS2(2,2,9) =  2.579530d0

  !! O-Cl
  MAIS2test(2,3) = .true.
  MAIS2(2,3,1) =  0.017759d0
  MAIS2(2,3,2) = -0.001574d0 
  MAIS2(2,3,3) = -0.065971d0 
  MAIS2(2,3,4) =  1.022950d0
  MAIS2(2,3,5) =  0.369802d0
  MAIS2(2,3,6) =  0.591341d0
  MAIS2(2,3,7) =  4.543160d0
  MAIS2(2,3,8) =  7.538830d0
  MAIS2(2,3,9) =  1.989830d0

  ! Store the upper triangle of the parameters matrix

  do j=1,nMAIS2atm
   do i=j+1,nMAIS2atm
    MAIS2test(i,j) = MAIS2test(j,i)
    MAIS2(i,j,:) = MAIS2(j,i,:)
   enddo
  enddo

 return
 
 end subroutine initMAIS2

 subroutine initPIF1()
  implicit none
  integer :: i,j

  PIF1(:,:,:) = 0.d0
  PIF1test(:,:) = .false.
  pointersPIF1(:) = 0


  !! Atom pointers

  pointersPIF1(1) = 1 ! H
  pointersPIF1(8) = 2 ! O

  !! H-H
  PIF1test(1,1) = .true.
  PIF1(1,1,1) =      0.5197876d0 ! alpha 
  PIF1(1,1,2) =      2.4899800d0 ! beta
  PIF1(1,1,3) =     50.1010700d0 ! chi
  PIF1(1,1,4) =   -676.3223000d0 ! delta
  PIF1(1,1,5) =   2711.7680000d0 ! epsilon
   
  !! H-O
  PIF1test(1,2) = .true.
  PIF1(1,2,1) =     49.3655400d0
  PIF1(1,2,2) =      2.2764670d0
  PIF1(1,2,3) =    -38.2875000d0
  PIF1(1,2,4) =   -145.0515000d0
  PIF1(1,2,5) =    651.8225000d0

  !! O-O
  PIF1test(2,2) = .true.
  PIF1(2,2,1) =     12.4595800d0
  PIF1(2,2,2) =      2.5027340d0
  PIF1(2,2,3) =   -266.2244000d0
  PIF1(2,2,4) =  12587.5300000d0
  PIF1(2,2,5) = -89486.9100000d0

  ! Store the upper triangle of the parameters matrix

  do j=1,nPIF1atm
   do i=j+1,nPIF1atm
    PIF1test(i,j) = PIF1test(j,i)
    PIF1(i,j,:) = PIF1(j,i,:)
   enddo
  enddo

 return
 
 end subroutine initPIF1
  
 subroutine initPIF2()
  implicit none
  integer :: i,j

  PIF2(:,:,:) = 0.d0
  PIF2test(:,:) = .false.
  pointersPIF2(:) = 0


  !! Atom pointers

  pointersPIF2(1) = 1 ! H
  pointersPIF2(6) = 2 ! C
  pointersPIF2(7) = 3 ! N
  pointersPIF2(8) = 4 ! O
  pointersPIF2(17)= 5 ! Cl

  !! H-H
  PIF2test(1,1) = .true.
  PIF2(1,1,1) =      2.4794900d0 ! alpha
  PIF2(1,1,2) =      2.8961200d0 ! beta
  PIF2(1,1,3) =     47.2620000d0 ! chi
  PIF2(1,1,4) =   -505.7300000d0 ! delta
  PIF2(1,1,5) =   1549.8100000d0 ! epsilon
  
  !! H-C
  PIF2test(1,2) = .true.
  PIF2(1,2,1) =      0.5312600d0
  PIF2(1,2,2) =      0.8305610d0
  PIF2(1,2,3) =   -353.4700000d0
  PIF2(1,2,4) =   4643.4000000d0
  PIF2(1,2,5) =  -7083.6300000d0

  !! H-N
  PIF2test(1,3) = .true.
  PIF2(1,3,1) =     33.1138100d0
  PIF2(1,3,2) =      1.7829930d0
  PIF2(1,3,3) =   -172.4870000d0
  PIF2(1,3,4) =    158.4300000d0
  PIF2(1,3,5) =   2423.9700000d0

  !! H-O
  PIF2test(1,4) = .true.
  PIF2(1,4,1) =     29.3251700d0
  PIF2(1,4,2) =      2.0926480d0
  PIF2(1,4,3) =    -55.7790000d0
  PIF2(1,4,4) =     44.5300000d0
  PIF2(1,4,5) =    313.4400000d0

  !! H-Cl
  PIF2test(1,5) = .true.
  PIF2(1,5,1) =     27.4253000d0
  PIF2(1,5,2) =      2.5211050d0
  PIF2(1,5,3) =    -19.9352200d0
  PIF2(1,5,4) =     29.0613800d0
  PIF2(1,5,5) =    280.3394000d0

  !! C-N
  PIF2test(2,3) = .true.
  PIF2(2,3,1) =      0.9469700d0
  PIF2(2,3,2) =      0.9813460d0
  PIF2(2,3,3) =      2.0850000d0
  PIF2(2,3,4) =      6.1200000d0
  PIF2(2,3,5) =    102.3000000d0

  !! C-O
  PIF2test(2,4) = .true.
  PIF2(2,4,1) =      0.8788400d0
  PIF2(2,4,2) =      0.8860320d0
  PIF2(2,4,3) =      1.8790000d0
  PIF2(2,4,4) =     33.9500000d0
  PIF2(2,4,5) =    843.6700000d0

  !! N-N
  PIF2test(3,3) = .true.
  PIF2(3,3,1) =      1.4353300d0
  PIF2(3,3,2) =      0.8415580d0
  PIF2(3,3,3) =      1.9410000d0
  PIF2(3,3,4) =     35.4500000d0
  PIF2(3,3,5) =    887.7400000d0

  !! N-O
  PIF2test(3,4) = .true.
  PIF2(3,4,1) =      1.9470200d0
  PIF2(3,4,2) =      1.0161940d0
  PIF2(3,4,3) =      1.9410000d0
  PIF2(3,4,4) =     35.4500000d0
  PIF2(3,4,5) =    887.7400000d0

  !! O-O
  PIF2test(4,4) = .true.
  PIF2(4,4,1) =     75.1546600d0
  PIF2(4,4,2) =      1.5120630d0
  PIF2(4,4,3) =    338.6560000d0
  PIF2(4,4,4) = -32185.6800000d0
  PIF2(4,4,5) = 274634.1000000d0

  !! O-Cl
  PIF2test(4,5) = .true.
  PIF2(4,5,1) =     35.2773200d0
  PIF2(4,5,2) =      1.2929060d0
  PIF2(4,5,3) =   -680.9553000d0
  PIF2(4,5,4) =     19.1353500d0
  PIF2(4,5,5) =  57189.9700000d0



  ! Store the upper triangle of the parameters matrix

  do j=1,nPIF2atm
   do i=j+1,nPIF2atm
    PIF2test(i,j) = PIF2test(j,i)
    PIF2(i,j,:) = PIF2(j,i,:)
   enddo
  enddo

 return
 
 end subroutine initPIF2

 subroutine initPIF3()
  implicit none
  integer :: i,j

  PIF3(:,:,:) = 0.d0
  PIF3test(:,:) = .false.
  pointersPIF3(:) = 0


  !! Atom pointers

  pointersPIF3(1) = 1 ! H
  pointersPIF3(8) = 2 ! O

  !! H-H
  PIF3test(1,1) = .true.
  PIF3(1,1,1) =    370.0761307d0
  PIF3(1,1,2) =      3.5857863d0
  PIF3(1,1,3) =     29.8960036d0
  PIF3(1,1,4) =   -315.7477361d0
  PIF3(1,1,5) =    513.5394769d0
   
  !! H-O
  PIF3test(1,2) = .true.
  PIF3(1,2,1) =     18.5576298d0
  PIF3(1,2,2) =      1.8252006d0
  PIF3(1,2,3) =    -78.1098760d0
  PIF3(1,2,4) =    277.1073977d0
  PIF3(1,2,5) =   -305.5240028d0

  ! Store the upper triangle of the parameters matrix

  do j=1,nPIF3atm
   do i=j+1,nPIF3atm
    PIF3test(i,j) = PIF3test(j,i)
    PIF3(i,j,:) = PIF3(j,i,:)
   enddo
  enddo

 return
 
 end subroutine initPIF3

end module se_corrections_params




















