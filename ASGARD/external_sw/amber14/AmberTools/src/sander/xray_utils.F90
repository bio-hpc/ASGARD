! <compile=optimized>
#include "../include/assert.fh"
! General utility functions and constants for X-ray code, mostly math.
! It may be good to make some of the common code and/or use existing utilities.
!
module xray_utils_module
   use xray_globals_module, only: real_kind, rk_, stdout
   use constants, only: M_D2R => DEG_TO_RAD, M_R2D => RAD_TO_DEG

   !----------------------------------------------------------------------------
   ! Generic function interfaces:
   interface polar
      module procedure polar_c4, polar_c8
   end interface polar

   interface complex_from_polar
      module procedure complex_from_polar_c4, complex_from_polar_c8
   end interface complex_from_polar

   interface amplitude
      module procedure amplitude_c4, amplitude_c8
   end interface amplitude

   interface phase
      module procedure phase_c4, phase_c8
   end interface phase

   interface operator(.cross.)
      module procedure vector3_cross
   end interface operator(.cross.)

   interface normalize
      module procedure vector3_normalize
   end interface normalize

   interface inverse
      module procedure matrix3_inverse_pure
   end interface inverse
   interface inverse
      module procedure matrix3_inverse_det
   end interface inverse
   interface determinant
      module procedure matrix3_determinant
   end interface determinant
   interface det
      module procedure matrix3_determinant
   end interface det

   interface operator(.mult.)
      module procedure coordinate_transform_multiply
   end interface operator(.mult.)

   interface operator(.mult.)
      module procedure vector3_matrix3_product
   end interface operator(.mult.)

contains

   pure function coordinate_transform_multiply(a,b) result(product)
      real(real_kind), dimension(3,4) :: product
      real(real_kind), dimension(3,4), intent(in) :: a,b
      integer i,j
      do i=1,3
         do j=1,3
            product(i,j) = sum(a(1:3,j) * b(i,1:3))
         end do
         product(i,4) = sum(a(1:3,i) * b(:,4)) + a(i,4)
      end do
   end function coordinate_transform_multiply

   pure function coordinate_transform_inverse(xform) result(inverse)
      real(real_kind), dimension(3,4) :: inverse
      real(real_kind), dimension(3,4), intent(in) :: xform
      inverse(1:3,1:3) = matrix3_inverse_pure(xform(1:3,1:3))
      inverse(1:3,4) = matmul(-1.0*xform(1:3,4),inverse(1:3,1:3))
   end function coordinate_transform_inverse

   pure function vector3_normalize(v)
      real(real_kind), intent(in) :: v(3)
      real(real_kind) :: vector3_normalize(3)
      ! Caution: this does not check for divide-by-zero, because it is pure.
      vector3_normalize(1:3)=v(1:3) / sqrt(sum(v**2))
   end function vector3_normalize

   pure function vector3_cross(b,c)
      real(real_kind), intent(in) :: b(3), c(3)
      real(real_kind) :: vector3_cross(3)
      vector3_cross(1) = b(2)*c(3) - c(2)*b(3)
      vector3_cross(2) = b(3)*c(1) - c(3)*b(1)
      vector3_cross(3) = b(1)*c(2) - c(1)*b(2)
   end function vector3_cross

   pure function matrix3_determinant(matrix) result(determinant)
      real(real_kind), intent(in) :: matrix(3,3)
      real(real_kind) :: determinant
      determinant = &
            matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2)) &
            - matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1)) &
            + matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
   end function matrix3_determinant

   pure function matrix3_inverse_pure(matrix) result(inverse)
      real(real_kind) :: inverse(3,3)
      real(real_kind), intent(in) :: matrix(3,3)
      real(real_kind) :: determinant
      real(real_kind), parameter :: epsilon = tiny(determinant)*100.0_rk_
      determinant = det(matrix)
      if ( abs(determinant) > epsilon ) then
         inverse(1,:)=vector3_cross(matrix(:,2),matrix(:,3))
         inverse(2,:)=vector3_cross(matrix(:,3),matrix(:,1))
         inverse(3,:)=vector3_cross(matrix(:,1),matrix(:,2))
         inverse = inverse / determinant
      else
         inverse = 0.0
      end if
   end function matrix3_inverse_pure

   function matrix3_inverse_det(matrix,determinant) result(inverse)
      real(real_kind)             :: inverse(3,3)
      real(real_kind), intent(in) :: matrix(3,3)
      real(real_kind), intent(out) :: determinant
      real(real_kind), parameter :: epsilon = tiny(determinant)*100.0_rk_
      determinant = det(matrix)
      if ( abs(determinant) > epsilon ) then
         inverse(:,1)=vector3_cross(matrix(:,2),matrix(:,3))
         inverse(:,2)=vector3_cross(matrix(:,3),matrix(:,1))
         inverse(:,3)=vector3_cross(matrix(:,1),matrix(:,2))
         inverse = inverse / determinant
      else
         inverse = 0.0
         determinant = 0.0
      end if
   end function matrix3_inverse_det

   pure function vector3_matrix3_product(vector,matrix) result(product)
      real(real_kind) :: product(3)
      real(real_kind), intent(in) :: vector(3)
      real(real_kind), intent(in) :: matrix(3,3)
      product(1) = sum(vector(1:3) * matrix(1:3,1))
      product(2) = sum(vector(1:3) * matrix(1:3,2))
      product(3) = sum(vector(1:3) * matrix(1:3,3))
   end function vector3_matrix3_product

   complex(4) elemental function polar_c4(c)
      complex(4), intent(in) :: c
      real(4) :: phase
      if (abs(c) < 1E-20) then
         phase=0
      else
         phase = real(atan2(aimag(c),real(c)) * M_R2D,kind=4)
         if (phase<0) phase=phase+360.0
      end if
      polar_c4 = cmplx(abs(c),phase)
   end function polar_c4

   complex(8) elemental function polar_c8(c)
      complex(8), intent(in) :: c
      real(8) :: phase
      if (abs(c) < 1D-20) then
         phase=0
      else
         phase = atan2(aimag(c),real(c)) * M_R2D
         if (phase<0) phase=phase+360.0
      end if
      polar_c8 = cmplx(abs(c),phase)
   end function polar_c8

   complex(4) elemental function complex_from_polar_c4(p) result(c)
      complex(4), intent(in) :: p
      real(4) :: phase, amplitude
      phase = real(aimag(p) * M_D2R,kind=4)
      amplitude = real(p, kind=4)
      c = amplitude*cmplx(cos(phase),sin(phase))
   end function complex_from_polar_c4

   complex(8) elemental function complex_from_polar_c8(p) result(c)
      complex(8), intent(in) :: p
      real(8) :: phase, amplitude
      phase = aimag(p) * M_D2R
      amplitude = real(p)
      c = amplitude*cmplx(cos(phase),sin(phase))
   end function complex_from_polar_c8

   real(4) elemental function amplitude_c4(c) result(amplitude)
      complex(4), intent(in) :: c
      amplitude = abs(c)
   end function amplitude_c4

   real(8) elemental function amplitude_c8(c) result(amplitude)
      complex(8), intent(in) :: c
      amplitude = abs(c)
   end function amplitude_c8

   real(4) elemental function phase_c4(c) result(phase)
      complex(4), intent(in) :: c
      if (abs(c)<1E-20) then
         phase=0.0
      else
         phase = real(atan2(aimag(c),real(c)) * M_R2D,kind=4)
         if (phase<0) phase=phase+360.0
      end if
   end function phase_c4

   real(8) elemental function phase_c8(c) result(phase)
      complex(8), intent(in) :: c
      if (abs(c)<1D-20) then
         phase=0.0
      else
         phase = atan2(aimag(c),real(c)) * M_R2D
         if (phase<0) phase=phase+360.0
      end if
   end function phase_c8

#if 0  /* not used at present  */
#ifdef MKL
   function erf_inv(x)
      real(real_kind), intent(in) :: x
      real(real_kind) :: erf_inv
      INTERFACE verfinv
        SUBROUTINE vserfinv(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
        SUBROUTINE vderfinv(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE verfinv
      real(real_kind) :: a(1)
      call verfinv(1,(/x/),a)
      erf_inv = a(1)
   end function erf_inv
   function erfc_inv(x)
      real(real_kind), intent(in) :: x
      real(real_kind) :: erfc_inv
      INTERFACE verfcinv
        SUBROUTINE vserfcinv(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
        SUBROUTINE vderfcinv(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE verfcinv
      real(real_kind) :: a(1)
      call verfcinv(1,(/x/),a)
      erfc_inv = a(1)
   end function erfc_inv
#else
   !----------------------------------------------------------
   !  Approximations to the inverse of the error function, from:
   !    J. M. Blair, C. A. Edwards, and J. H. Johnson,
   !    "Rational Chebyshev Approximations for the Inverse
   !    of the Error Function", Mathematics of Computation,
   !    30 (1976) 827-830 + microfiche appendix.
   pure function erf_inv(x)
      real(real_kind), intent(in) :: x
      real(real_kind) :: erf_inv
      real(real_kind) :: t
      real(real_kind), parameter :: p1(0:2) = (/ -13.0959967422_rk_, 26.785225760_rk_, &
            -9.289057635_rk_ /)
      real(real_kind), parameter :: q1(0:3) = (/ -12.0749426297_rk_, 30.960614529_rk_, &
            -17.149977991_rk_, 1.00000000_rk_ /)
      real(real_kind), parameter :: p2(0:3) = (/ -0.12402565221_rk_, 1.0688059574_rk_, &
            -1.9594556078_rk_, 0.4230581357_rk_ /)
      real(real_kind), parameter :: q2(0:3) = (/ -0.08827697997_rk_, 0.8900743359_rk_, &
            -2.1757031196_rk_, 1.0000000000_rk_ /)
      real(real_kind), parameter :: p3(0:5) = &
            (/ 0.1550470003116_rk_, 1.382719649631_rk_, &
            0.690969348887_rk_, -1.128081391617_rk_, &
            0.680544246825_rk_, -0.16444156791_rk_ /)
      real(real_kind), parameter :: q3(0:2) = &
            (/ 0.155024849822_rk_, 1.385228141995_rk_, &
            1.000000000000_rk_ /)

      ! This approximation, taken from Table 10 of Blair et al., is valid
      ! for |x|<=0.75 and has a maximum relative error of 4.47 x 10^-8.
      if(abs(x) <= 0.75_rk_) then
         t = x**2 - 0.75_rk_**2
         erf_inv = x*(p1(0)+t*(p1(1)+t*p1(2)))/(q1(0)+t*(q1(1)+t*(q1(2)+t*q1(3))))

         ! This approximation, taken from Table 29 of Blair et al., is valid
         ! for .75<=|x|<=.9375 and has a maximum relative error of 4.17 x 10^-8.
      else if(abs(x) <= 0.9375_rk_) then
         t = x**2 - 0.9375_rk_**2
         erf_inv = x*(p2(0)+t*(p2(1)+t*(p2(2)+t*p2(3)))) &
               / (q2(0)+t*(q2(1)+t*(q2(2)+t*q2(3))))

         ! This approximation, taken from Table 50 of Blair et al.
         ! It is valid for .9375 <= |x| <= 1-10^-100 and has a
         ! maximum relative error of 2.45 x 10^-8.
      else if(abs(x) < 1.0) then
         t=1.0_rk_/sqrt(-log(1.0_rk_-abs(x)))
         erf_inv=sign((p3(0)/t+p3(1)+t*(p3(2)+t*(p3(3)+t*(p3(4)+t*p3(5))))) &
               / (q3(0)+t*(q3(1)+t*(q3(2)))),x)
      else
         erf_inv = HUGE(1.0_rk_)
      end if

   end function erf_inv

   pure function erfc_inv(x)
      real(real_kind), intent(in) :: x
      real(real_kind) :: erfc_inv
      real(real_kind) :: t, x1
      logical :: negative
      real(real_kind), parameter :: p1(0:5) = &
            (/ 0.1550470003116_rk_, 1.382719649631_rk_, &
            0.690969348887_rk_, -1.128081391617_rk_, &
            0.680544246825_rk_, -0.16444156791_rk_ /)
      real(real_kind), parameter :: q1(0:2) = &
            (/ 0.155024849822_rk_, 1.385228141995_rk_, &
            1.000000000000_rk_ /)
      real(real_kind), parameter :: p2(0:3) = &
            (/ 0.00980456202915_rk_, 0.363667889171_rk_, &
            0.97302949837_rk_,-0.5374947401_rk_ /)
      real(real_kind), parameter :: q2(0:2) = &
            (/ 0.00980451277802_rk_,0.363699971544_rk_, &
            1.000000000000_rk_ /)

      negative = x > 1.0_rk_
      if (negative) then
         x1 = 1.0_rk_-x
      else
         x1 = x
      end if

      if ( x >= 0.0625_rk_) then
         erfc_inv = erf_inv(1.0_rk_ - x)

         ! Table 50 of Blair et al.
      else if(x >= 1.0e-100_rk_) then
         t = 1.0_rk_/sqrt(-log(x))
         erfc_inv = (p1(0)/t+p1(1)+t*(p1(2)+t*(p1(3) &
               +t*(p1(4)+t*p1(5)))))/(q1(0)+t*(q1(1)+t*q1(2)))

         ! This approximation, taken from Table 70 of Blair et al.
         ! It is valid for 10^-1000 <= |x| <= 10^-100 and has a
         ! maximum relative error of 2.45 x 10^-8.
      else if(x > tiny(1.0_rk_)) then
         t = 1.0_rk_/sqrt(-log(x))
         erfc_inv = (p2(0)/t+p2(1)+t*(p2(2)+t*p2(3)))/(q2(0)+t*(q2(1)+t*q2(2)))
      else if (x<0.0_rk_) then
         !erfc_inv = nan8 !0.0_rk_/0.0_rk_ ! Bad input; return Not-a-Number
         erfc_inv = 0 !FIXME?
      else
         erfc_inv = HUGE(1.0_rk_)
      end if
      if (negative) erfc_inv = -erfc_inv
   end function erfc_inv
#endif /* MKL */
#endif /* 0 */

   ! Given a logical-slection array, return an array of indices into the
   ! selection's target array, and the number of selected indices.
   subroutine pack_index(selection,index,count)
      implicit none
      logical, intent(in) :: selection(:)
      integer, intent(out) :: index(size(selection,1))
      integer, intent(out), optional :: count
      ! locals
      integer :: i, n
      ! begin
      n=0
      do i = 1,size(selection)
         if (selection(i)) then
            n = n + 1
            index(n) = i
         end if
      end do
      if (present(count)) count = n
   end subroutine pack_index

   subroutine write_map_ccp4(unit, extent, rho, mask, scale)
      use xray_globals_module, only: &
            unit_cell, spacegroup_number, &
            num_symmops, symmop, grid_size, stdout
      implicit none
      integer, intent(in) :: unit
      integer :: extent(2,3)
      real(real_kind) :: rho(0:grid_size(1), &
                             0:grid_size(2), &
                             0:grid_size(3))
      integer, intent(in), optional :: mask( &
                             0:grid_size(1), &
                             0:grid_size(2), &
                             0:grid_size(3))
      logical, intent(in), optional :: scale
      ! local
      logical :: qscale
      real(real_kind) :: map_min,map_max,map_avg,map_sigma
      integer :: a, b, c, nr, isym, u, v, w
      integer :: i, j
      integer :: abc(3), as, bs, cs
      real(real_kind) :: sect_avg, sect_sigma, ed
      equivalence(as,abc(1))
      equivalence(bs,abc(2))
      equivalence(cs,abc(3))
      integer :: order_uvw(3)
      integer :: alloc_stat
      real(real_kind), allocatable :: section(:,:)
      real(real_kind), parameter :: r_epsilon = 1e-20
      if (present(scale)) then
         qscale = scale
      else
         qscale = .false.
      end if
      !-----------------------------------------------------
      ! get statistics for the whole map:
      map_min=0
      map_max=0
      map_avg=0
      map_sigma=0
      nr=0
      do c=lbound(rho,3),ubound(rho,3)
         do b=lbound(rho,2),ubound(rho,2)
            do a=lbound(rho,1),ubound(rho,1)
               if (present(mask)) then
                  if (mask(a,b,c) == 0) cycle
               end if
               ed=rho(a,b,c)
               map_avg=map_avg+ed
               map_sigma=map_sigma+ed**2
               if (nr==0) then
                  map_min=ed
                  map_max=ed
                  nr=1
               else
                  map_min=min(map_min,ed)
                  map_max=max(map_max,ed)
                  nr=nr+1
               end if
            end do
         end do
      end do

      ! print overall statistics
      map_avg=map_avg/nr
      map_sigma=sqrt(max(0.0_rk_,map_sigma/nr-map_avg**2))
      write(stdout,'(A,2(F14.5,A))') &
            ' MAP: avg. density in unit cell=', map_avg, &
            ' sigma=',map_sigma,' e/A^3'
      write(stdout,'(2(A,F14.5))') &
            ' MAP: minimum =', map_min, '  maximum =',map_max

      if (qscale.and.map_sigma < r_epsilon) then
         write(stdout,'(A)') &
               ' XMAPX: sigma very small for map.  No scaling performed.'
         map_sigma=1
      end if

      ! This also sets order_uvw()
      call write_header(unit)
      allocate(section(extent(1,order_uvw(1)):extent(2,order_uvw(1)), &
                       extent(1,order_uvw(2)):extent(2,order_uvw(2))), &
                       stat=alloc_stat)
      REQUIRE(alloc_stat==0)
      isym=1
      SLOW_AXIS: do w=extent(1,order_uvw(3)),extent(2,order_uvw(3))
         abc(order_uvw(3))=w
         sect_avg=0
         sect_sigma=0
         nr=0
         section=0
         MEDIUM_AXIS: do v=extent(1,order_uvw(2)),extent(2,order_uvw(2))
            abc(order_uvw(2))=v
            FAST_AXIS: do u=extent(1,order_uvw(1)),extent(2,order_uvw(1))
               abc(order_uvw(1))=u
               ! Iterate isym only if the current one does not match,
               ! because the same symmop will be used in the majority of cases.
               DO_SYMMOP: do i = 1,num_symmops
                  do j=1,3
                     if (modulo(nint(sum(symmop(1:3,j,isym)*abc(:)) &
                           + symmop(j,4,isym)*grid_size(j)),grid_size(j)) &
                           >= grid_size(j)) then
                        isym=modulo(isym,num_symmops)+1
                        cycle DO_SYMMOP
                     end if
                  end do
                  exit DO_SYMMOP
               end do DO_SYMMOP
               if (present(mask)) then
                  if (mask(as,bs,cs) == 0) cycle FAST_AXIS
               end if
               ed=rho(as,bs,cs)
               section(u,v)=ed
               sect_avg=sect_avg+ed
               sect_sigma=sect_sigma+ed**2
               nr=nr+1
            end do FAST_AXIS
         end do MEDIUM_AXIS

         ! print info about the section
         sect_avg=sect_avg/nr
         sect_sigma=sqrt(max(0.0_rk_,sect_sigma/nr-sect_avg**2))
         write(stdout,'(A,I4,A,F14.5,A,F14.5,A)') &
               ' MAP: section #',w,' avg. dens.=',sect_avg, &
               ' sigma=',sect_sigma,' e/A^3'

         ! scale the section
         if (qscale) section=section-map_avg/map_sigma

         write(unit) section

      end do SLOW_AXIS

   contains

      ! NOTE: (U,V,W) === (Fast,Medium,Slow)
      ! UNIT MUST be opened for RAW BINARY I/O !
      subroutine write_header(unit)
         integer, intent(in) :: unit
         !-------------------------------------------------------------------------
         ! MAP header structure:
         integer(4) :: map_header(256)

         integer(4) :: h_map_uvw_extent(3) ! Map Extent
         integer(4) :: h_map_mode
         integer(4) :: h_map_uvw_origin(3) ! Map Origin
         integer(4) :: h_map_xyz_size(3) ! Grid size
         real(4)    :: h_map_unit_cell(6)
         integer(4) :: h_map_order_uvw(3)
         real(4)    :: h_map_min,h_map_max,h_map_mean
         integer(4) :: h_map_spacegroup_number
         integer(4) :: h_map_symmop_num_chars
         integer(4) :: h_map_skew_flag
         real(4)    :: h_map_skew_matrix(3,3), h_map_skew_translation(3)
         integer(4) :: h_map_pad(15)
         integer(4) :: h_map_tag, h_map_arch_stamp
         real(4)    :: h_map_rms
         integer(4) :: h_map_num_labels
         integer(4) :: h_map_int_labels(20,10) ! == 10 * char(len=80)

         equivalence(map_header( 1),h_map_uvw_extent)
         equivalence(map_header( 4),h_map_mode)
         equivalence(map_header( 5),h_map_uvw_origin)
         equivalence(map_header( 8),h_map_xyz_size)
         equivalence(map_header(11),h_map_unit_cell)
         equivalence(map_header(17),h_map_order_uvw)
         equivalence(map_header(20),h_map_min)
         equivalence(map_header(21),h_map_max)
         equivalence(map_header(22),h_map_mean)
         equivalence(map_header(23),h_map_spacegroup_number)
         equivalence(map_header(24),h_map_symmop_num_chars)
         equivalence(map_header(25),h_map_skew_flag)
         equivalence(map_header(26),h_map_skew_matrix)
         equivalence(map_header(35),h_map_skew_translation)
         equivalence(map_header(36),h_map_pad)
         equivalence(map_header(53),h_map_tag)
         equivalence(map_header(54),h_map_arch_stamp)
         equivalence(map_header(55),h_map_rms)
         equivalence(map_header(56),h_map_num_labels)
         equivalence(map_header(57),h_map_int_labels)

         !-------------------------------------------------------------------------
         !  spacegroups axis ordering: 1=YXZ, 2=ZXY
         integer, parameter :: axis_order_index(23) = &
               (/2,2,2,2,1,1,1,1,1,2,1,1,1,1,1,2,2,2,1,2,2,1,2/)
         integer :: axis_order_uvw(3,2) = reshape((/  2,1,3,  3,1,2  /),(/3,2/))
         character(len=80) :: label
         integer(4) :: native_IT, native_FT

         ! integer :: j, l, ld, lu, lt
         integer :: i, sg
         character(len=60) :: filename
         ! character(len=12) :: username
         ! character(len=11) :: date
         ! character(len=8) :: time
         !-------------------------------------------------------------------------
         ! Allow for CCP4-style alternate origin numbers:
         sg = mod(spacegroup_number,1000)

         if (sg>0 .and. sg<size(axis_order_index)) then
            order_uvw = axis_order_uvw(:,axis_order_index(sg))
         else
            order_uvw = axis_order_uvw(:,1)
         end if

         map_header = 0
         h_map_int_labels=ICHAR(' ')
         h_map_order_uvw = order_uvw
         h_map_uvw_origin(:) = extent(1,order_uvw(:))
         h_map_uvw_extent(:) = extent(2,order_uvw(:)) - h_map_uvw_origin(:) + 1
         h_map_mode=2
         h_map_xyz_size(:) = grid_size(:)
         h_map_unit_cell(:) = real(unit_cell,kind=4)
         h_map_spacegroup_number = spacegroup_number
         h_map_max = real(map_max, kind=4)
         h_map_min = real(map_min, kind=4)
         h_map_mean = real(map_avg, kind=4)
         h_map_rms = real(map_sigma,kind=4)

         !----------------------------------------------------------
         ! Generate the endianness flag:
         native_IT = IAND(15,transfer(char(4)//char(3)//char(2)//char(1),0_4))
         native_FT = native_IT
         h_map_arch_stamp = transfer( CHAR(native_FT*16+native_FT) &
               // CHAR(native_IT*16+1) &
               // CHAR(0) // CHAR(0), &
               h_map_arch_stamp)
         h_map_tag = transfer('MAP ',0_4)

         !----------------------------------------------------------
         inquire(unit=unit,name=filename)
         label='FILENAME="'//trim(filename)//'"'
         h_map_int_labels(:,1) = transfer(label,h_map_int_labels(:,1))

         !label='DATE:'//date(1:ld)//'  '// &
               !      time(1:lt)//'       created by user: '//trim(user_name)
         !h_map_int_labels(:,2) = transfer(label,h_map_int_labels(:,2))

         !h_map_int_labels(:,3) = transfer(label,h_map_int_labels(:,3))

         ! Add up to 7 lines of remark titles:
         !i=3
         !do while(i-3<ntitle .and. i<10)
         !   i=i+1
         !   label = title(i-3)
         !   h_map_int_labels(:,i) = transfer(label,h_map_int_labels(:,i))
         !end do
         i=1
         h_map_num_labels = i

         label=" "
         i=i+1
         do while(i<10)
            i=i+1
            h_map_int_labels(:,i) = transfer(label,h_map_int_labels(:,i))
         end do

         write(unit) map_header

      end subroutine write_header
   end subroutine write_map_ccp4

   subroutine write_map_ezd(unit, extent, rho) ! mask, scale not used
      use xray_globals_module, only: unit_cell, grid_size, stdout
      implicit none
      integer, intent(in) :: unit
      integer :: extent(2,3)
      real(real_kind) :: rho(0:grid_size(1), &
                             0:grid_size(2), &
                             0:grid_size(3))
      ! local
      real(real_kind) :: map_scale
      ! begin
      map_scale = 1
      write(unit,'(A)')        'EZD_MAP'
      write(unit,'(A,6F10.5)') 'CELL ',unit_cell
      write(unit,'(A,3I6)')    'ORIGIN ',extent(1,:)
      write(unit,'(A,3I6)')    'EXTENT ',extent(2,:)-extent(1,:)+1
      write(unit,'(A,3I6)')    'GRID ',grid_size
      write(unit,'(A,F10.5)')  'SCALE ',map_scale
      write(unit,'(A)')        'MAP'
      call write_data(rho,product(extent(2,:)-extent(1,:)+1))
      write(unit,'(A)')        'END'
   contains
      subroutine write_data(data,ndata)
         integer, intent(in) :: ndata
         real(real_kind), intent(in) :: data(ndata)
         ! local
         integer :: i, len
         character(len=(7*12)) :: line
         ! begin
         do i=1,ndata,7
            write(line,'(7f12.5)') data(i:min(i+6,ndata))
            call pack_string(line,len)
            write(unit,'(A)') line(1:len)
         end do
      end subroutine write_data
   end subroutine write_map_ezd

   subroutine pack_string(string,len)
      implicit none
      character(len=*), intent(inout) :: string
      integer, intent(out) :: len
      ! local
      integer :: i, n
      ! begin
      len=len_trim(string)
      i=1
      do while (i <= len)
         if (string(i:i) == ' ') then
            n = blank_len(string(i:))
            if (n > 1) then
               n=n-1
               string(i+1:(len-n))=string((i+1+n):len)
               len=len-n
            end if
         end if
         i=i+1
      end do
      return
   end subroutine pack_string

   function blank_len(string) result(length)
      integer :: length
      character(len=*) :: string
      ! local
      integer :: i
      ! begin
      i=1
      do while (string(i:i)==' ')
         i=i+1
         if (i>len(string)) exit
      end do
      length = i-1
   end function blank_len

   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! Allocate the next unused logial UNIT number from a predefined range.
   function allocate_lun(lun_return) result(lun)
      integer :: lun
      integer, intent(out), optional :: lun_return
      integer, parameter :: FILE_UNIT_FIRST=201, FILE_UNIT_LAST=250
      integer, save :: unit = FILE_UNIT_FIRST-1
      ! locals
      integer :: i
      logical :: opened
      ! begin
      lun = FILE_UNIT_FIRST
      do i = FILE_UNIT_FIRST, FILE_UNIT_LAST
         unit = unit + 1
         if (unit > FILE_UNIT_LAST) unit = FILE_UNIT_FIRST
         inquire(unit=unit,opened=opened)
         if (.not. opened) then
            lun = unit
            if (present(lun_return)) lun_return = unit
            return
         end if
      end do
      write(stdout,'(A)') &
            'ERROR in ALLOCATE_LUN(): ran out of available LUNs!'
      call mexit(stdout,2)
   end function allocate_lun

end module xray_utils_module
