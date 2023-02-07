! <compile=optimized>
#include "../include/assert.fh"
module xray_reciprocal_space_module
   ! Reflection and spacegroup symmetry operations
   use xray_globals_module
   implicit none

   !---------------------------------------------------------
   ! Typedefs:

   type hkl_data_scale_type
      real(real_kind), pointer :: f(:) => NULL()
      real(real_kind), pointer :: sigma(:) => NULL()
      real(real_kind) :: scale=1
      logical :: refine_scale = .TRUE.
      real(real_kind) :: bfactor(6)=0
      logical :: refine_bfactor = .TRUE.
      logical :: anisotropic = .TRUE. ! FALSE means constrain B11=B22=B33
      !                   and B12=B13=B23==0
      logical :: isotropic = .TRUE.   ! FALSE means constrain B11+B22+B33==0
      character(len=16) :: name = ""  ! Name for output messages
   end type hkl_data_scale_type

   !---------------------------------------------------------
   ! Constants:

   ! Lattice types: P, A, B, C, I, F, or R
   character(len=1), parameter :: lattice_types(230) = (/ &
         'P','P','P','P','C',  'P','P','C','C','P', & ! 10
   'P','C','P','P','C',  'P','P','P','P','C', & ! 20
   'C','F','I','I','P',  'P','P','P','P','P', & ! 30
   'P','P','P','P','C',  'C','C','A','A','A', & ! 40
   'A','F','F','I','I',  'I','P','P','P','P', & ! 50
   'P','P','P','P','P',  'P','P','P','P','P', & ! 60
   'P','P','C','C','C',  'C','C','C','F','F', & ! 70
   'I','I','I','I','P',  'P','P','P','I','I', & ! 80
   'P','I','P','P','P',  'P','I','I','P','P', & ! 90
   'P','P','P','P','P',  'P','I','I','P','P', & !100
   'P','P','P','P','P',  'P','I','I','I','I', & !110
   'P','P','P','P','P',  'P','P','P','I','I', & !120
   'I','I','P','P','P',  'P','P','P','P','P', & !130
   'P','P','P','P','P',  'P','P','P','I','I', & !140
   'I','I','P','P','P',  'R','P','R','P','P', & !150
   'P','P','P','P','R',  'P','P','P','P','R', & !160
   'R','P','P','P','P',  'R','R','P','P','P', & !170
   'P','P','P','P','P',  'P','P','P','P','P', & !180
   'P','P','P','P','P',  'P','P','P','P','P', & !190
   'P','P','P','P','P',  'F','I','P','I','P', & !200
   'P','F','F','I','P',  'I','P','P','F','F', & !210
   'I','P','P','I','P',  'F','I','P','F','I', & !220
   'P','P','P','P','F',  'F','F','F','I','F' /) !230

   ! Laue type code, but with an extra group to distinguish P321 from P312.
   integer, parameter :: recip_au_types(230) = (/ &
         1,  1,  2,  2,  2,    2,  2,  2,  2,  2, & ! 10
   2,  2,  2,  2,  2,    3,  3,  3,  3,  3, & ! 20
   3,  3,  3,  3,  3,    3,  3,  3,  3,  3, & ! 30
   3,  3,  3,  3,  3,    3,  3,  3,  3,  3, & ! 40
   3,  3,  3,  3,  3,    3,  3,  3,  3,  3, & ! 50
   3,  3,  3,  3,  3,    3,  3,  3,  3,  3, & ! 60
   3,  3,  3,  3,  3,    3,  3,  3,  3,  3, & ! 70
   3,  3,  3,  3,  4,    4,  4,  4,  4,  4, & ! 80
   4,  4,  4,  4,  4,    4,  4,  4,  5,  5, & ! 90
   5,  5,  5,  5,  5,    5,  5,  5,  5,  5, & !100
   5,  5,  5,  5,  5,    5,  5,  5,  5,  5, & !110
   5,  5,  5,  5,  5,    5,  5,  5,  5,  5, & !120
   5,  5,  5,  5,  5,    5,  5,  5,  5,  5, & !130
   5,  5,  5,  5,  5,    5,  5,  5,  5,  5, & !140
   5,  5,  6,  6,  6,    6,  6,  6,  8,  7, & !150
   8,  7,  8,  7,  8,    7,  8,  7,  8,  7, & !160
   8,  7,  8,  7,  8,    7,  8,  9,  9,  9, & !170
   9,  9,  9,  9,  9,    9, 10, 10, 10, 10, & !180
   10, 10, 10, 10, 10,   10, 10, 10, 10, 10, & !190
   10, 10, 10, 10, 11,   11, 11, 11, 11, 11, & !200
   11, 11, 11, 11, 11,   11, 11, 11, 12, 12, & !210
   12, 12, 12, 12, 12,   12, 12, 12, 12, 12, & !220
   12, 12, 12, 12, 12,   12, 12, 12, 12, 12 /) !230

   ! This is the au_code index:
   !------------------------------------------------------------------
   !code, Laue group, CCP4 recipriocal asymmetric-unit limits
   !====  ==========  =======================================
   ! 1    PG1         h k l : l >= 0
   !                  h k 0 : h >= 0
   !                  0 k 0 : k >= 0
   !
   ! 2    P2/m        h k l : h >= 0, l >= 0
   !                  h k 0 : h >= 0
   !
   ! 3    Pmmm        h k l : h >= 0, k >= 0, l >= 0
   !
   ! 4    P4/m        h k l : h >= 0, k >= 0, l >= 0
   !                  0 k l : k >  0
   !                  0 0 l : l >  0
   !
   ! 5    P4/mmm      h k l : h >= 0, k >= 0, l >= 0, h >= k
   !                  0 0 l : l >  0
   !
   ! 6    P3 (R3)     h k l : h >= 0, k > 0
   !                  0 0 l : l >  0
   !
   ! 7    P312        h k l : h >= 0, k >= 0, k <= h   (all l)
   !                  h 0 l : l >= 0
   !
   ! 8    P321        h k l : h >= 0, k >= 0, k <= h   (all l)
   !                  h h l : l >= 0
   !
   ! 9    P6/m        h k l : h >= 0, k >= 0, l >= 0
   !                  0 k l : k >  0
   !                  0 0 l : l >  0
   !
   ! 10   P6/mmm      h k l : h >= 0, k >= 0, l >= 0, h >= k
   !                  0 0 l : l >  0
   !
   ! 11   P32         h k l : h >= 0, k >= 0, l >= 0, l>=h,
   !             and  h k h : k >= h
   !                  h k l : k >  h if l >  h
   !
   ! 12   P432        h k l : h >= 0, k >= 0, l >= 0, k >= l and l>= h
   !------------------------------------------------------------------


   integer, parameter :: reciprocal_octants_start(3,12) = reshape ( (/ &
         -1,-1, 1,  & ! 1)  -1     Triclinic
         -1, 1, 1,  & ! 2)  2/m    Monoclinic
          1, 1, 1,  & ! 3)  mmm    Orthorhombic
          1,-1, 1,  & ! 4)  4/m    Tetragonal
          1, 1, 1,  & ! 5)  4/mmm  Tetragonal
          1, 1,-1,  & ! 6)  -3     Trigonal (including Rhombohedral)
          1, 1,-1,  & ! 7)  -3m    Trigonal (312 group)
          1, 1,-1,  & ! 8)  -3m    Trigonal (321 group)
          1, 1, 1,  & ! 9)  6/m    Hexagonal
          1, 1, 1,  & ! 10) 6/mmm  Hexagonal
          1, 1, 1,  & ! 11) m-3    Cubic
          1, 1, 1   & ! 12) m-3m   Cubic
          /), (/3,12/))

   ! This is the symmetry_system index:
   integer, parameter :: &
         SYMM_TRICLINIC = 1, &
         SYMM_MONOCLINIC = 2, &
         SYMM_ORTHORHOMBIC = 3, &
         SYMM_TETRAGONAL = 4, &
         SYMM_TRIGONAL = 5, &
         SYMM_HEXAGONAL = 6, &
         SYMM_CUBIC = 7

   ! Get the symmetery system for a given au_code
   integer, parameter :: symmetry_system(12) = (/ &
         SYMM_TRICLINIC, &
         SYMM_MONOCLINIC, &
         SYMM_ORTHORHOMBIC, &
         SYMM_TETRAGONAL, &
         SYMM_TETRAGONAL, &
         SYMM_TRIGONAL, &
         SYMM_TRIGONAL, &
         SYMM_TRIGONAL, &
         SYMM_HEXAGONAL, &
         SYMM_HEXAGONAL, &
         SYMM_CUBIC, &
         SYMM_CUBIC &
         /)

   ! Scale refinement constants:
   integer, parameter :: B11=1, B22=2, B33=3, B12=4, B13=5, B23=6
   integer, parameter :: params_per_set = 7
   integer, parameter :: i_scale=7 ! index to scale member of a parameter group

   ! Status return: negative values are error indicators
   integer, parameter :: &
         STATUS_BAD_MATRIX = -3, &
         STATUS_ALL_FIXED  = -2, &
         STATUS_DIVERGENT  = -1, &
         STATUS_CONVERGED  =  1, &
         STATUS_MAX_CYCLES =  2

   real(real_kind), parameter :: REAL_EPSILON = 1E-20_rk_
   real(real_kind), parameter :: COS_EPSILON = 0.001_rk_

!==============================================================================
contains

   subroutine scale_data( &
            num_hkl,hkl_index,hkl_selected, &
            num_sets,data, &
            scale_min,scale_max, &
            bfactor_min,bfactor_max, &
            max_cycles,tolerance, &
            reference_set,status, print )
      use xray_utils_module, only: pack_index
      implicit none
      integer, intent(in) :: num_hkl
      integer, intent(in) :: hkl_index(3,num_hkl)
      logical, intent(in), optional :: hkl_selected(num_hkl)

      integer, intent(in) :: num_sets
      type(hkl_data_scale_type), intent(inout) :: data(num_sets)

      real(real_kind), intent(in) :: scale_min, scale_max
      real(real_kind), intent(in) :: bfactor_min, bfactor_max
      integer, intent(in) :: max_cycles
      real(real_kind), intent(in) :: tolerance
      integer, intent(in), optional :: reference_set
      integer, intent(out) :: status
      logical, intent(in) :: print

      !-----------------------------------------------------------------
      integer :: num_params, num_refined
      real(real_kind), allocatable :: sig2i(:)
      real(real_kind) :: rcell(3)

      logical :: constrain_rhombohedral ! B12=B13=B23
      integer :: constrain_anisotropy
      integer, parameter :: constrain_B11_B22_B33 = 7 ! b'111'
      integer, parameter :: constrain_B11_B22     = 6 ! b'110'
      integer, parameter :: constrain_B11_B33     = 5 ! b'101'
      integer, parameter :: constrain_B22_B33     = 3 ! b'011'

      real(real_kind) :: reference_sum
      integer :: ref_set

      !-----------------------------------------------------------------
      integer :: i,Icycle,i_set, alloc_status
      real(real_kind) :: cos_rangle(3)
      logical :: anisoB_refine(6), anisoB_constrain(6)
      real(real_kind) :: lambda, Btrace
      real(real_kind) :: chi_squared,chi_squared_try
      real(real_kind) :: rfactor,rfactor_try

      !-------------------------------------------------------------------------
      ! Automatic arrays:
      integer :: refindex(num_sets*params_per_set)
      logical :: param_refine(num_sets*params_per_set)
      real(real_kind), dimension(num_sets*params_per_set) :: parameters
      real(real_kind), dimension(num_sets*params_per_set) :: param_try
      ! Allocatable arrays:
      real(real_kind), allocatable :: alpha(:,:), alpha_try(:,:), &
            beta(:), beta_try(:,:)
      !-------------------------------------------------------------------------

      num_params = params_per_set * num_sets
      rcell = recip_cell(1:3)

      !=============================================================
      ! First step:
      ! Determine which parameters are refined, and
      ! remove redundant parameters from aniso-B refinement.
      !
      !-------------------------------------------------------------
      ! Impose restrictions according to the crystal symmetry
      ! Triclinc         [none]
      ! Monoclinic                       B13=B23=0 when beta=alpha=90
      !                                  B12=B23=0 when gamma=alpha=90
      !                                  B12=B13=0 when gamma=beta=90
      ! Orthorhombic                     B12=B13=B23=0
      ! Tetragonal       B11=B22     and B12=B13=B23=0
      ! Rhombohedral     B11=B22=B33 and B12=B13=B23
      ! Hexagonal        B11=B22     and B13=B23=0
      ! Cubic            B11=B22=B33 and B12=B13=B23=0 (=isotropic)

      if (any(data%anisotropic)) then
         cos_rangle = cos(recip_cell(4:6))

         anisoB_refine(1:3) = .TRUE.
         anisoB_refine(4:6) = abs(cos_rangle(3:1:-1)) > COS_EPSILON
         anisoB_constrain = .FALSE.

         ! The following test works for any unit_cell choice, but assumes
         ! that values are never coincidentally exactly equal.
         constrain_rhombohedral = .FALSE.
         constrain_anisotropy = 0
         if (rcell(1)==rcell(2)) then
            if (rcell(2)==rcell(3)) then
               constrain_anisotropy = constrain_B11_B22_B33
               anisoB_refine(1:3) = (/.TRUE.,.FALSE.,.FALSE./)
               if (cos_rangle(1)==cos_rangle(2) &
                        .and. cos_rangle(2)==cos_rangle(3) &
                        .and. abs(cos_rangle(1)) > COS_EPSILON ) then
                  constrain_rhombohedral = .TRUE.
                  anisoB_refine(4:6) = (/.TRUE.,.FALSE.,.FALSE./)
               end if
            else ! a==b/=c
               constrain_anisotropy = constrain_B11_B22
               anisoB_refine(1:3) = (/.TRUE.,.FALSE.,.TRUE./)
            end if
         else if (rcell(3)==rcell(1)) then
            constrain_anisotropy = constrain_B11_B33
            anisoB_refine(1:3) = (/.TRUE.,.TRUE.,.FALSE./)
         else if (rcell(2)==rcell(3)) then
            constrain_anisotropy = constrain_B22_B33
            anisoB_refine(1:3) = (/.TRUE.,.TRUE.,.FALSE./)
         end if
      end if

      !-------------------------------------------------------------
      do i=1,num_sets
         param_refine((i-1)*params_per_set+i_scale) = data(i)%refine_scale
         if (.not.data(i)%refine_bfactor) then
            param_refine((i-1)*params_per_set+B11 : &
                  (i-1)*params_per_set+B23) = .FALSE.
         else if (data(i)%anisotropic) then
            param_refine((i-1)*params_per_set+B11 : &
                  (i-1)*params_per_set+B23) = anisoB_refine
         else
            param_refine((i-1)*params_per_set+B11 : &
                  (i-1)*params_per_set+B23) = &
                  (/.TRUE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE./)
         end if
      end do
      call pack_index(param_refine,refindex,num_refined)

      if (num_refined==0) then
         status = STATUS_ALL_FIXED
         return
      end if

      allocate(alpha(num_refined,num_refined), &
            alpha_try(num_refined,num_refined), &
            beta(num_refined), &
            beta_try(num_refined,1), &
            stat=alloc_status)
      REQUIRE(alloc_status==0)

      do i_set = 1,num_sets
         if (associated(data(i_set)%sigma)) then
            if (.not. allocated(sig2i)) then
               allocate(sig2i(num_hkl),stat=alloc_status)
               REQUIRE(alloc_status==0)
               sig2i = (data(i_set)%sigma)**2
            else
               sig2i = sig2i + (data(i_set)%sigma)**2
            end if
         end if
      end do
      if (allocated(sig2i)) sig2i=1.0_rk_ / sig2i

      !--------------------------------------------------------------------
      do i = 1,num_sets
         parameters((i-1)*params_per_set+B11: &
               (i-1)*params_per_set+B23) = data(i)%bfactor
         parameters((i-1)*params_per_set+i_scale) = data(i)%scale
      end do

      !--------------------------------------------------------------------
      status = STATUS_MAX_CYCLES

      lambda = 1.0E-4
      if (present(reference_set)) then
         ref_set = reference_set
      else
         ref_set=1
      end if
      reference_sum = sum(data(ref_set)%f)
      call evaluate_params(parameters,alpha,beta,chi_squared,rfactor,.false.)

      if( print ) then
         write(stdout,'(A)') repeat('-',80)
         write(stdout,'(A)') "Initial state:"
         call write_status()
      endif

      LSQ_CYCLE: do Icycle = 1,max_cycles

         beta_try(:,1) = beta
         alpha_try = alpha
         !forall (i=1:num_refined) alpha_try(i,i)=alpha_try(i,i)+lambda
         forall (i=1:num_refined)
            alpha_try(i,i)=alpha_try(i,i) * (1.0_rk_+lambda)
         end forall

         call gauss_jordan_matrix_invert(alpha_try,beta_try)
         if (maxval(beta_try) < tolerance) then
            status = STATUS_CONVERGED
            exit LSQ_CYCLE
         end if

         param_try = parameters+unpack(beta_try(:,1),param_refine,0.0_rk_)
         call evaluate_params(param_try,alpha_try, &
                beta_try(:,1),chi_squared_try,rfactor_try,.false.)

         if (chi_squared_try <= chi_squared) then
         end if

         if (chi_squared_try <= chi_squared) then
            if( print ) then
               write(stdout,'(A)') repeat('-',80)
               write(stdout,'(A,I4,A,ES10.3)')  &
                  'Cycle:',Icycle,'  d(chi^2)=',chi_squared_try-chi_squared
            endif
            lambda = 0.1_rk_*lambda
            chi_squared = chi_squared_try
            rfactor = rfactor_try
            alpha = alpha_try
            beta = beta_try(:,1)
            parameters = param_try
            if( print ) call write_status()
         else if (chi_squared_try == chi_squared) then
            status = STATUS_CONVERGED
         else
            lambda = 10.0_rk_*lambda
         end if

         !---------------------------------------------------------------
         if (lambda>1E+6) exit

         if (status == STATUS_CONVERGED) exit LSQ_CYCLE

      end do LSQ_CYCLE

      !--------------------------------------------------------------------
      do i = 1,num_sets
         if (data(i)%anisotropic) then
            data(i)%bfactor = &
                  parameters((i-1)*params_per_set+B11:(i-1)*params_per_set+B23)
         else
            data(i)%bfactor(1:3) = parameters((i-1)*params_per_set+B11)
            data(i)%bfactor(4:6) = 0.0
         end if
         data(i)%scale = parameters((i-1)*params_per_set+i_scale)
      end do

      if( print ) then
         write(stdout,'(A)') repeat('-',80)
         write(stdout,'(A)') "Final state:"
         call write_status()
         write(stdout,'(A)') repeat('-',80)
         call evaluate_params(parameters,alpha,beta,chi_squared,rfactor,print)
      endif

      return

      !--------------------------------------------------------------------
   contains
      subroutine write_status()
         implicit none
         write(stdout,'(A,ES6.0,3(A,F9.5))') &
               'lambda=',lambda,', chi^2=',chi_squared,', R-factor=',rfactor
         write(stdout,'(A5,14X,A5,3X,A4,6(5X,A3))') &
               'name:','scale','Biso','B11','B22','B33','B12','B13','B23'
         do i = 1,num_sets
            Btrace = sum(parameters((i-1)*params_per_set+B11: &
                                    (i-1)*params_per_set+B33))/3.0
            write(stdout,'(A16,F8.4,7(1X,F7.3))') &
                  data(i)%name, &
                  parameters((i-1)*params_per_set+i_scale), Btrace, &
                  parameters((i-1)*params_per_set+B11: &
                                    (i-1)*params_per_set+B23) - Btrace
         end do
      end subroutine write_status

      subroutine evaluate_params(params,alpha,beta,chi_squared,rfactor,print)
         implicit none
         real(real_kind), intent(inout), target :: params(:)
         real(real_kind), intent(out) :: alpha(:,:), beta(:)
         real(real_kind), intent(out) :: chi_squared, rfactor
         logical, intent(in) :: print

         real(real_kind), dimension(num_sets*params_per_set), target :: deriv
         integer :: i,j,i_hkl, i_set
         real(real_kind) :: residual, weight, residual_sum, reference_sum
         real(real_kind) :: f, mqSS(6)
         real(real_kind), pointer :: p_scale, p_bfactor(:)
         real(real_kind), pointer :: d_scale, d_bfactor(:)
         real(real_kind) :: foo(3)

         integer, parameter :: Bi_index(6) = (/ 1,2,3,1,1,2 /)
         integer, parameter :: Bj_index(6) = (/ 1,2,3,2,3,3 /)
         real(real_kind), parameter :: Bij_const(6) = &
               (/-0.25,-0.25,-0.25,-0.50,-0.50,-0.50/)
         real(real_kind) :: S(3)

         !--------------------------------------------------------------------
         ! Apply parameter limit constraints:
         do i_set = 1,num_sets
            p_bfactor => params((i_set-1)*params_per_set+1: &
                                (i_set-1)*params_per_set+6)
            p_scale => params((i_set-1)*params_per_set+i_scale)

            if (.not. data(i_set)%anisotropic) then
               p_bfactor(B11:B33) = max(min(sum(p_bfactor(B11:B33)) / &
                                3.0,bfactor_max),bfactor_min)
            end if
            if (data(i_set)%scale > 0 ) then
               p_scale = max(min(p_scale,scale_max),scale_min)
            else
               p_scale = max(min(p_scale,-scale_min),-scale_max)
            end if

         end do
         !--------------------------------------------------------------------

         alpha = 0.0
         beta = 0.0
         chi_squared = 0.0
         rfactor = 0.0
         residual_sum = 0.0
         reference_sum = 0.0

         HKL: do i_hkl = 1,num_hkl

            if (present(hkl_selected)) then
               if (.not. hkl_selected(i_hkl)) cycle
            end if

            S = real(hkl_index(:,i_hkl),real_kind) * rcell(:)
            mqSS = Bij_const * S(Bi_index)*S(Bj_index)

            residual = 0.0
            SET: do i_set = 1,num_sets
               p_scale   => params((i_set-1)*params_per_set+i_scale)
               p_bfactor => params((i_set-1)*params_per_set+1: &
                                   (i_set-1)*params_per_set+6)
               d_scale   =>  deriv((i_set-1)*params_per_set+i_scale)
               d_bfactor =>  deriv((i_set-1)*params_per_set+1: &
                                   (i_set-1)*params_per_set+6)

               f = data(i_set)%f(i_hkl) * exp(sum(mqSS*p_bfactor))
               d_scale = f
               f = f * p_scale
               foo(i_set) = f
               residual = residual - f
               d_bfactor = mqSS * f

               if (p_scale>0.0) reference_sum = reference_sum+f

               !----------------------------------------------------------------
               ! Apply derivative constraints:
               if ( .false. .and.  &
                        data(i_set)%refine_bfactor) then
                  if (.not. data(i_set)%anisotropic) then
                     d_bfactor(1:3) = sum(d_bfactor(1:3)) !!!!/3.0  !!FIXME!!
                     d_bfactor(4:6) = 0.0
                  else
                     if (constrain_rhombohedral) then
                        d_bfactor(4:6) = sum(d_bfactor(4:6))!!/3.0
                     end if
                     select case(constrain_anisotropy)
                     case(constrain_B11_B22_B33)
                     d_bfactor(1:3) = sum(d_bfactor(1:3))!!/3.0
                     case(constrain_B11_B22)
                     d_bfactor(1:2) = sum(d_bfactor(1:2))!!/2.0
                     case(constrain_B11_B33)
                     d_bfactor(1:3:2) = sum(d_bfactor(1:3:2))!!/2.0
                     case(constrain_B22_B33)
                     d_bfactor(2:3) = sum(d_bfactor(2:3))!!/2.0
                     end select
                  end if
               end if
               !----------------------------------------------------------------
            end do SET

            if (allocated(sig2i)) then
               do i = 1,num_refined
                  weight = deriv(refindex(i)) * sig2i(i_hkl)
                  forall (j=1:i)
                  alpha(i,j)=alpha(i,j) + weight * deriv(refindex(j))
                  end forall
                  beta(i) = beta(i)+residual*weight
               end do
               chi_squared = chi_squared + residual**2 * sig2i(i_hkl)
            else
               do i = 1,num_refined
                  weight = deriv(refindex(i))
                  forall (j=1:i)
                  alpha(i,j)=alpha(i,j) + weight * deriv(refindex(j))
                  end forall
                  beta(i) = beta(i)+residual*weight
               end do
               chi_squared = chi_squared + residual**2
            end if
            residual_sum = residual_sum + abs(residual)
            if( print ) &
               write(6,'(3i5,2f12.2)') hkl_index(1,i_hkl),hkl_index(2,i_hkl), &
               hkl_index(3,i_hkl), foo(1),-foo(2)

         end do HKL

         do i = 2,num_refined
            forall (j=1:i-1) alpha(j,i)=alpha(i,j)
         end do

         ! partially normalized chi-squared: approximately an R^2-factor
         chi_squared = chi_squared / reference_sum
         rfactor = residual_sum/reference_sum
         !!rfactor = rfactor / reference_sum * &
         !!    abs(params((ref_set-1)*params_per_set+i_scale))
         !write(stdout,*) 'Rsum=',residual,'; ref_sum=',ref_sum,'; R=',rfactor

         ! From HKL Package:
         !     R linear = SUM ( ABS(I - <I>)) / SUM (I)
         !     R square = SUM ( (I - <I>) ** 2) / SUM (I ** 2)
         !     Chi**2   = SUM ( (I - <I>) ** 2) / (Error ** 2 * N / (N-1) ) )

      end subroutine evaluate_params

      subroutine gauss_jordan_matrix_invert(a,b)
         implicit none
         real(real_kind), intent(inout) :: a(:,:), b(:,:)
         external :: dgesv
         integer :: lda, ldb, n, nrhs, info
         integer :: ipiv(size(a,2))
         !_ASSERT(size(a,1)>=size(a,2))
         !_ASSERT(size(a,2)==size(b,2))
         lda = size(a,1)
         n = size(a,2)
         ldb = size(b,1)
         nrhs = size(b,2)
         call dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      end subroutine gauss_jordan_matrix_invert

   end subroutine scale_data

   !+ Derive matrices, lookup symmetry operators, etc., from unit cell.
   !+ Mostly equivalent to Ewald get_ucell().
   subroutine derive_cell_info(a,b,c,alpha,beta,gamma,ncode)
      use xray_utils_module, only: M_D2R, M_R2D, inverse
      implicit none
      real(real_kind), intent(in) :: a,b,c,alpha,beta,gamma
      integer, intent(in), optional :: ncode ! 1 is standard for PDB and AMBER.
      ! GLOBALS:
      ! input: cell -> unit_cell
      ! output: recip_cell, frac_to_orth, orth_to_frac, cell_volume

      real(real_kind) :: sin_alpha, sin_beta, sin_gamma
      real(real_kind) :: cos_alpha, cos_beta, cos_gamma
      real(real_kind) :: cos_ralpha, cos_rbeta, cos_rgamma
      real(real_kind) :: sin_ralpha, sin_rbeta, sin_rgamma
      integer :: ncode_

      unit_cell = (/a,b,c,alpha,beta,gamma/)

      cos_alpha = cos(alpha * M_D2R)
      sin_alpha = sin(alpha * M_D2R)
      cos_beta  = cos(beta  * M_D2R)
      sin_beta  = sin(beta  * M_D2R)
      cos_gamma = cos(gamma * M_D2R)
      sin_gamma = sin(gamma * M_D2R)

      cell_volume = a * b * c &
            * sqrt(1.0D0 + 2.0D0*cos_alpha*cos_beta*cos_gamma &
            - cos_alpha**2 - cos_beta**2 - cos_gamma**2 )

      cos_ralpha=(cos_beta *cos_gamma-cos_alpha)/(sin_beta *sin_gamma)
      sin_ralpha=sqrt(1.0D0 - cos_ralpha**2)
      cos_rbeta =(cos_gamma*cos_alpha-cos_beta )/(sin_gamma*sin_alpha)
      sin_rbeta=sqrt(1.0D0 - cos_rbeta**2)
      cos_rgamma=(cos_alpha*cos_beta -cos_gamma)/(sin_alpha*sin_beta )
      sin_rgamma=sqrt(1.0D0 - cos_rgamma**2)

      recip_cell(1) = b * c * sin_alpha / cell_volume
      recip_cell(2) = a * c * sin_beta  / cell_volume
      recip_cell(3) = a * b * sin_gamma / cell_volume
      recip_cell(4) = acos(cos_ralpha) * M_R2D
      recip_cell(5) = acos(cos_rbeta ) * M_R2D
      recip_cell(6) = acos(cos_rgamma) * M_R2D

      ! NOTE:  recip_volume = 1.0 / cell_volume

      ! Fractional to orthogonal matrix is equivalent to
      ! the 3 real-space unit-cell axis vectors.
      ! Here, we use a fixed orthogonalization choice.
      ! (Note: there is no translation; origins are always the same for both.)

      !   ncode = 1 - A along X   C* along Z
      !   ncode = 2 - B along X   A* along Z
      !   ncode = 3 - C along X   B* along Z
      !   ncode = 4 - A+B along X   C* along Z
      !   ncode = 5 - C along Z   A* along X
      !   ncode = 6 - A along X   B* along Y

      if (present(ncode)) then
         ncode_ = ncode
      else
         ncode_ = 1
      end if
      select case (ncode_)
      case (1)
      frac_to_orth=reshape( &
            (/ A,           0.0_rk_,                0.0_rk_, &
               B*cos_gamma, B*sin_gamma,           0.0_rk_, &
               C*cos_beta, -C*sin_beta*cos_ralpha, C*sin_beta*sin_ralpha /), &
            (/3,3/));
      case(2)
      frac_to_orth=reshape( &
            (/ A*cos_gamma, -A*sin_gamma*cos_rbeta, A*sin_gamma*sin_rbeta, &
               B,            0.0_rk_,                0.0_rk_, &
               C*cos_alpha,  C*sin_alpha,           0.0_rk_ /), &
            (/3,3/));
      case(3)
      frac_to_orth=reshape( &
            (/ A*cos_beta,   A*sin_beta,             0.0_rk_, &
               B*cos_alpha, -B*sin_alpha*cos_rgamma, B*sin_alpha*sin_rgamma, &
               C,            0.0_rk_,                 0.0_rk_ /), &
            (/3,3/));
      case(4)
      frac_to_orth=reshape( &
            (/ A/2.0, -A*sin_gamma, 0.0_rk_, &
               A/2.0,  A*sin_gamma, 0.0_rk_, &
               0.0_rk_, 0.0_rk_,      C /), &
            (/3,3/));
      case(5)
      frac_to_orth=reshape( &
            (/ A*sin_beta*sin_rgamma, -A*sin_beta*cos_rgamma, A*cos_beta, &
               0.0_rk_,                 B*sin_alpha,           B*cos_alpha, &
               0.0_rk_,                 0.0_rk_,                C /), &
            (/3,3/));
      case(6)
      frac_to_orth=reshape( &
            (/ A,           0.0_rk_,                  0.0_rk_, &
               B*cos_gamma, B*sin_gamma*sin_ralpha, -B*sin_gamma*cos_ralpha, &
               C*cos_beta,  0.0_rk_,                  C*sin_beta /), &
            (/3,3/));
      end select

      orth_to_frac = transpose(inverse(frac_to_orth))

   end subroutine derive_cell_info

   subroutine hkl_generate(resolution_range)
      use xray_globals_module, only: hkl_index
      implicit none
      real(real_kind), intent(in) :: resolution_range(2)
      ! GLOBAL: integer, allocatable :: hkl_index(:,:)

      real(real_kind) :: S_min, S_max, S_previous, S
      integer :: hstep, kstep, lstep, h, k, l, alloc_status
      logical :: hfound, kfound
      integer :: nrefl ! reflection array size
      integer :: nhkl ! number of HKLs added to array
      integer :: hkl_start(3)

      if (resolution_range(1)>resolution_range(2)) then
         S_min = 1.0_rk_ / resolution_range(1)
         S_max = 1.0_rk_ / resolution_range(2)
      else
         S_min = 1.0_rk_ / resolution_range(2)
         S_max = 1.0_rk_ / resolution_range(1)
      end if

      hkl_start = reciprocal_octants_start(:,au_type)

      nhkl=0
      nrefl=1000
      allocate(hkl_index(3,nrefl),stat=alloc_status)
      REQUIRE(alloc_status==0)
      H_STEP: do hstep = hkl_start(1),1,2
         ! Start negatives at one to avoid double-counting zero plane.
         H = merge(0,-1,(hstep>0))
         hfound=.true.
         H_SEARCH: do ! Keep stepping in h until nothing found in k steps
            !write(6,'(A,I5)') "H=",h; call flush(6)
            hfound=.false.
            K_STEP: do kstep = hkl_start(2),1,2
               K = merge(0,-1,(kstep>0))
               K_SEARCH: do ! Keep stepping in k until nothing found in l steps
                  kfound=.false.
                  L_STEP: do lstep = hkl_start(3),1,2
                     L = merge(0,-1,(lstep>0))
                     ! Step in L if: S is below max res. (S_max),
                     !            OR the resolution is getting lower.
                     S_previous=0.0_rk_
                     L_SEARCH: do
                        S=hkl_S((/H,K,L/),orth_to_frac)
                        if (S>S_min .and. S<S_max) then
                           ! Reflection in resolution range
                           kfound=.true.
                           if (hkl_octant_au((/h,k,l/)) &
                                    .and. .not. hkl_is_sysabs(h,k,l) ) then
                              ! Reflection in recip. au
                              if (nhkl>=size(hkl_index,1)) then
                                 if (nhkl>1000000) then
                                    write(stderr,'(A)') &
            'In hkl_generate(): Too many reflections generated (>1000000) BUG?'
                                    call mexit(stdout,1)
                                 end if
                                 ! Reallocate hkl_index:
                                 nrefl = nrefl + min(1000,nrefl/4)
                                 call realloc(nrefl)
                              end if

                              nhkl=nhkl+1
                              hkl_index(:,nhkl) = (/h,k,l/)
                           end if
                        end if
                        ! Exit the loop if we have gone outside of
                        ! the resolution limit.
                        ! Also ensure that subsequent steps are
                        ! toward higher resolution (or the first step).
                        if (S>S_max .and. S>S_previous &
                              .and. S_previous>0.0) then
                           exit L_SEARCH
                        end if
                        S_previous=S
                        L=L+lstep
                     end do L_SEARCH
                  end do L_STEP
                  if(.not.kfound) exit K_SEARCH
                  hfound=hfound.or.kfound
                  K=K+kstep
               end do K_SEARCH
            end do K_STEP
            H=H+hstep
            if(.not.hfound) exit H_SEARCH
         end do H_SEARCH
      end do H_STEP

      ! Free excess array storage;
      ! shrink HKL allocation back to the actual size.
      nrefl = nhkl
      call realloc(nrefl)

      call sort_hkl(hkl_index,nrefl)

   contains
      subroutine realloc(nrefl)
         integer, intent(in) :: nrefl
         integer, allocatable :: hkl_index_save(:,:)
         integer :: ncopy
#ifdef HAVE_MOVE_ALLOC
         ! move_alloc() removes an extra copy step, but is F2003.
         call move_alloc(hkl_index,hkl_index_save)
#else
         allocate(hkl_index_save(3,size(hkl_index,2)),stat=alloc_status)
         REQUIRE(alloc_status==0)
         hkl_index_save = hkl_index
         deallocate(hkl_index)
#endif
         allocate(hkl_index(3,nrefl),stat=alloc_status)
         REQUIRE(alloc_status==0)
         ncopy=min(size(hkl_index_save,2),nrefl)
         hkl_index(:,1:nrefl)=hkl_index_save(:,1:nrefl)
         deallocate(hkl_index_save)
      end subroutine realloc
   end subroutine hkl_generate

   ! Sort reflections by H,K,L
   subroutine sort_hkl(hkl,hkl_size)
      implicit none
      integer, intent(in) :: hkl_size
      integer, intent(inout) :: hkl(3,hkl_size)
      integer, parameter :: QSORT_THRESHOLD = 64
      integer, parameter :: SORT_DIMENSION = 2

      ! The sort array must be named 'array'.
#define array hkl
#include "sort_inline.inc"
#undef array

   contains
      pure function sort_compare(a,b) result(less_than)
         implicit none
         logical :: less_than
         integer, intent(in) :: a, b
         less_than = .false.
         if (hkl(1,a)<hkl(1,b)) then
            less_than=.true.
         else if (hkl(1,a)==hkl(1,b)) then
            if (hkl(2,a)<hkl(2,b)) then
               less_than = .true.
            else if (hkl(2,a)==hkl(2,b)) then
               if (hkl(3,a)<hkl(3,b)) then
                  less_than = .true.
               end if
            end if
         end if
      end function sort_compare
      subroutine sort_swap(a,b)
         implicit none
         integer, intent(in) :: a,b
         integer :: hold(3)
         hold=hkl(:,a)
         hkl(:,a)=hkl(:,b)
         hkl(:,b)=hold
      end subroutine sort_swap
      subroutine sort_shift(a,b)
         implicit none
         integer, intent(in) :: a,b
         integer :: hold(3)
         hold=hkl(:,b)
         hkl(:,a+1:b)=hkl(:,a:b-1)
         hkl(:,a)=hold
      end subroutine sort_shift
   end subroutine sort_hkl

   function hkl_in_au(hkl) result(result)
      implicit none
      logical :: result
      integer, intent(in) :: hkl(3)
      integer :: i
      integer :: oct(3)
      oct = reciprocal_octants_start(:,au_type)
      result = .false.
      do i=1,3
        if (hkl(i)<0 .and. oct(i)>0) return
      end do
      result = hkl_octant_au(hkl)
   end function hkl_in_au
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !+ Check if HKL is in recip AU; must already be in valid octant
   function hkl_octant_au(hkl) result(result)
      implicit none
      logical :: result
      integer, intent(in) :: hkl(3)
      integer h,k,l
      h=hkl(1); k=hkl(2); l=hkl(3)
      select case(au_type)
      case (1);  result = ((l/=0 .or. h>=0) .and. (h/=0 .or. l/=0 .or. k>=0))
      case (2);  result = (l/=0 .or. h>=0)
      case (3);  result = .true. ! Only requires octant constraint
      case (4);  result = ((h/=0 .or. k>=0) .and. (h<=0 .or. k>0) )
      case (5);  result = (h>=k)
      case (6);  result = (k/=0)
      case (7);  result = (h>=k .and. (h/=0 .or. l>=0))  ! h0l requires l>=0
      case (8);  result = (h>=k .and. (h/=k .or. l>=0) ) ! hhl requires l>=0
      case (9);  result = ((h/=0 .or. k>=0) .and. (h<=0 .or. k>0))
      case (10); result = (h>=k)
      case (11); result = (l>=h .and. (l/=h .or. k>=h) .and. (l<=h .or. k>h))
      case (12); result = (k>=l .and. l>=h )
      end select
   end function hkl_octant_au

   subroutine hkl_to_au(hkl,symmop,friedel)
      implicit none
      integer, intent(inout) :: hkl(3)
      integer, intent(out), optional :: symmop
      logical, intent(out), optional :: friedel
      integer, save :: isym = 1
      integer :: i, hkl2(3)
      if (present(friedel)) friedel = .false.
      do i=1,num_symmops
         hkl2=int(matmul(transpose(symmop_inv(1:3,1:3,isym)),hkl))
         if (hkl_in_au(hkl2)) then
            hkl=hkl2
            if (present(symmop)) symmop=i
            return
         end if
         isym = modulo(isym,num_symmops) + 1
      end do
      if (present(friedel)) friedel = .true.
      hkl = -hkl
      do i=1,num_symmops
         hkl2=int(matmul(transpose(symmop_inv(1:3,1:3,isym)),hkl))
         if (hkl_in_au(hkl2)) then
            hkl=hkl2
            if (present(symmop)) symmop=i
            return
         end if
      end do
      stop 'BUG: cannot map reflection to eh recip AU!'

!      phase_shift = M_TWOPI*sum(hkl*symmop(:,4,isym))/RTH
!      if (friedel) then
!        refl = refl*dconjg(cmplx(cos(-phase_shift),sin(-phase_shift)))
!      else
!        refl = refl*cmplx(cos(phase_shift),sin(phase_shift))
!      end if
!      ! More work for HL coeffs!
   end subroutine hkl_to_au

! Leave HKL sorting as an external utility?
   subroutine hkl_reduce()
     implicit none
     !integer :: i
     !do i=1,num_hkl
     !   call hkl_to_au(hkl_index(:,i))
     !end do
     stop 'FIXME'
   end subroutine hkl_reduce

   function hkl_S(hkl,orth_to_frac) result(S)
      implicit none
      real(kind=real_kind) :: S
      integer, intent(in) :: hkl(3)
      real(kind=real_kind), intent(in) :: orth_to_frac(3,3)
      S = sqrt(sum(matmul(hkl,orth_to_frac)**2))
   end function hkl_S

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  1 TRICLINIC -1
   !  2 MONOCLINIC, 2/m
   !  3 ORTHOROMBIC, mmm
   !  4 TETRAGONAL, 4/m
   !  5 TETRAGONAL, 4/mmm
   !  6 TRIGONAL,-3
   !  7 TRIGONAL, 312
   !  8 TRIGONAL, 321
   !  9 HEXAGONAL 6/m
   !! 10 HEXAGONAL 6/mmm
   !! 11 CUBIC m3
   !! 12 CUBIC m3m
   !
   !subroutine hkl_au_oct(sgnum, h,k,l) ! Define octants for recip. sp. au
   !  integer sgnum, h, k, l
   !  integer htab(12),ktab(12),ltab(12)
   !  integer asymmetric_unit_octants(3,12)
   !
   !  type crystal_lattice_type
   !    character(len=12) :: system
   !    character(len=7) :: point_group
   !    character(len=1) :: lattice_type
   !    integer :: laue_type
   !    integer :: octant(3)
   !  end type crystal_lattice_type
   !  type(crystal_lattice_type) :: crystal_lattice(12) = (/ &
   !       crystal_lattice_type( "Triclinic",     "-1", 'P', 1, [-1,-1, 1]) &
   !       crystal_lattice_type( "Monoclinic",   "2/m", 'P', 1, [-1,-1, 1]) &
   !       crystal_lattice_type( "Monoclinic",   "2/m", 'P', 1, [-1,-1, 1]) &
   !       crystal_lattice_type( "Orthorhombic", "mmm", 'P', 1, [-1,-1, 1]) &
   !       crystal_lattice_type( "Orthorhombic", "mmm", 'P', 1, [-1,-1, 1]) &
   !       crystal_lattice_type( "Orthorhombic", "mmm", 'P', 1, [-1,-1, 1]) &
   !       crystal_lattice_type( "Orthorhombic", "mmm", 'P', 1, [-1,-1, 1]) &
   !       crystal_lattice_type( "Tetragonal",   "4/m", 'P', 1, [-1,-1, 1]) &
   !       crystal_lattice_type( "Tetragonal", "4/mmm", 'P', 1, [-1,-1, 1]) &
   !       crystal_lattice_type( "Trigonal",      "-3", 'P', 1, [-1,-1, 1]) &
   !       crystal_lattice_type( "Trigonal",     "312", 'P', 1, [-1,-1, 1]) &
   !       crystal_lattice_type( "Trigonal",     "321", 'P', 1, [-1,-1, 1]) &
   !       crystal_lattice_type( "Hexagonal",    "6/m", 'P', 1, [-1,-1, 1]) &
   !       crystal_lattice_type( "Hexagonal",  "6/mmm", 'P', 1, [-1,-1, 1]) &
   !       crystal_lattice_type( "Cubic",         "m3", 'P', 1, [-1,-1, 1]) &
   !       crystal_lattice_type( "Cubic",        "m3m", 'P', 1, [-1,-1, 1]) &
   !       crystal_lattice_type( "Cubic",        "m3m", 'P', 1, [-1,-1, 1]) &
   !  /)
   !
   !  data htab  / -1, -1,  1,  1,   1,   1,   1,   1,   1,   1,   1,   1 /
   !  data ktab  / -1,  1,  1, -1,   1,   1,   1,   1,   1,   1,   1,   1 /
   !  data ltab  /  1,  1,  1,  1,   1,  -1,  -1,  -1,   1,   1,   1,   1 /
   !
   !  h=htab(au_type_index(sgnum))
   !  k=ktab(au_type_index(sgnum))
   !  l=ltab(au_type_index(sgnum))
   !!end subroutine
   !!
   !subroutine sg_name_to_num(sgname,sgnum)
   !  character(len=*)sgname
   !  integer sgnum
   !
   !  do sgnum=1,230
   !    if (sgname_index(sgnum) == sgname) return
   !  end do
   !  sgnum=0
   !
   !end subroutine sg_name_to_num
   !
   !
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !+ determine if h,k,l is a systematic absence
   !
   ! Uses logic from CCP4, so reciprocal AU is CCP4 convention.
   function hkl_is_sysabs(h,k,l) result(absent)
      implicit none
      logical :: absent
      integer, intent(in) :: h,k,l

      integer :: isym
      real(real_kind) :: phase
      ! max. rounding error for phase, in degrees:
      real(real_kind), parameter :: PHASE_EPSILON = 0.05
      real(real_kind) :: hkl(3)

      if (num_symmops<2) then
         absent=.FALSE.
         return
      end if

      ! Check for general absence conditions:
      absent = .TRUE.
      select case(lattice_type)
      case ('P'); continue
      case ('A'); if (mod(k+l,2) /= 0) return
      case ('B'); if (mod(h+l,2) /= 0) return
      case ('C'); if (mod(h+k,2) /= 0) return
      case ('I'); if (mod(k+k+l,2) /= 0) return
      case ('R'); if (mod(k+l-h,3) /= 0) return
      case ('F'); if (iand(h+k,1) /= 0 .or. iand(k+l,1) /= 0) return
      case default; STOP 'ERROR: bad symmetry group ID'
      end select

      ! All other systematic absences must lie in a principle zone.
      if ((h/=0).and.(k/=0).and.(l/=0).and.((h/=k).or.(h/=l))) then
         absent = .FALSE.
         return
      end if

      ! Systematic absence test, for reflections in a principle zone:
      ! For any pair of symm. operations, if (h' k' l')==(h k l) AND the phases
      ! differ, it is a systematic absence. (Note: If phases match, it is centric)
      hkl = (/h,k,l/)
      do isym = 2, num_symmops
         if ( all( matmul(hkl, symmop(1:3,1:3,isym)) == hkl ) ) then
            phase = sum(hkl*symmop(1:3,4,isym))
            if (abs(phase-nint(phase)) > PHASE_EPSILON) then
               absent = .TRUE.
               return
            end if
         end if
      end do
      absent = .FALSE.
   end function hkl_is_sysabs

end module xray_reciprocal_space_module
