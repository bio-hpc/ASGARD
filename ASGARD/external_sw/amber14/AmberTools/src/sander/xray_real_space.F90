! <compile=optimized>
#include "../include/assert.fh"
module xray_real_space_module
   !  grid_atomic_density     --  Compute the unit-cell electron density from coordinates.
   !                              Used to set up the real-space grid for FFT.
   !
   use xray_globals_module, only: real_kind, rk_, stdin, stdout, stderr
   implicit none
   public

contains
   !----------------------------------------------------------------------------
   subroutine grid_atomic_density( &
            num_atoms, au_coor, tempFactor, occupancy, scatter_type, &
            num_scatter_types, num_scatter_coeff, scatter_coeffs, &
            d_tempFactor, d_occupancy, d_XYZ)

      !-------------------------------------------------------------------------
      use xray_globals_module, only: density_map, fft_bfactor_sharpen, &
            fft_density_tolerance, fft_radius_min, fft_radius_max, &
            angstrom_to_grid, frac_to_grid, grid_to_orth, grid_size
      use constants, only: M_PI => PI, M_FOURPI => FOURPI
      implicit none
      type real_ptr
         real(real_kind), pointer :: p
      end type real_ptr

      integer,         intent(in) :: num_atoms
      real(real_kind), intent(in) :: au_coor(3,num_atoms)
      real(real_kind), intent(in) :: tempFactor(num_atoms), occupancy(num_atoms)
      integer,         intent(in) :: scatter_type(num_atoms)
      integer,         intent(in) :: num_scatter_coeff
      integer,         intent(in) :: num_scatter_types
      real(real_kind), intent(in) :: &
                           scatter_coeffs(2,num_scatter_coeff,num_scatter_types)
      real(real_kind), intent(inout), optional :: &
            d_tempFactor(num_atoms), d_occupancy(num_atoms), d_XYZ(3,num_atoms)

      integer :: i, j, iatom, num_points
      integer :: dx, dy, dz, ix, iy, iz
      integer :: grid_point(3)
      real(real_kind) :: rcut, rcut_squared, rcut_max
      real(real_kind) :: grid_xyz(3), delta_xyz(3), dist_squared
      real(real_kind), dimension(3) :: rvec_z, rvec_yz, rvec_xyz, offset
      logical :: density_calc

      integer :: dir, max_points, alloc_status

      real(real_kind) :: usage, min_usage, max_usage

      real(real_kind), dimension(num_scatter_coeff) :: acoeffs, nbcoeffs
      real(real_kind) :: c0,c1,c2
      integer :: sphere_extent(3)

      ! ------- Allocatable/Pointer arrays:
      real(real_kind), allocatable :: radius_squared(:), radius_vector(:,:)
      real(real_kind), allocatable :: d(:), a(:), e(:)
      type(real_ptr), allocatable :: points(:)

      real(real_kind), parameter :: REAL_EPSILON = 1e-30_rk_
      real(real_kind), parameter :: ACOEFF_EPSILON = 1e-30_rk_
      real(real_kind), parameter :: M_FOURPI2 = M_FOURPI * M_PI
      ! ----------------------------------------------------------------------

      write(*,*)'num scatter types=',num_scatter_types
      write(*,*)'num_atoms = ',num_atoms
      write(*,'(A,3F8.3)') 'FFT grid spacing = ',angstrom_to_grid**(-1)

      ! ----------------------------------------------------------------------
      density_calc = .not. ( present(d_tempFactor) &
                     .or. present(d_occupancy) .or. present(d_XYZ) )

      ! ----------------------------------------------------------------------
      ! Determine the largest cutoff radius, and allocate arrays accordingly.
      rcut_max = 0
      do iatom = 1,num_atoms
         call get_real_space_coeffs(iatom,acoeffs,nbcoeffs,rcut)
         rcut_max = max(rcut_max,rcut)
      end do
      write(*,*) 'Rcut_max=',rcut_max,', FFT_radius_min=',fft_radius_min
      if (rcut_max<fft_radius_min) then
         write(stderr,'(A)') &
          'ERROR in xray_real_space.f:grid_atomic_density(): Rcut_max is too small!'
         !call mexit(stderr,1)
      end if

      !FIXME: update for non-orthognal cells
      sphere_extent(:) = int(rcut_max * angstrom_to_grid(:)) + 1
      max_points = product(sphere_extent(:) * 2 + 1)
#ifdef _XRAY_DEBUG
      write(*,*)'sphere_extent=',sphere_extent
      write(*,*)'max_points=',max_points
#endif

      allocate( &
            points(max_points), &          !  grid pointers (for data scattering)
            radius_squared(max_points), &  !  R^2, for vectorized calculations
            radius_vector(3,max_points), &  !  R-vector, for vectorized calculations

            ! temporary arrays for vectorized loops:
            d(max_points), &
            a(max_points), &
            e(max_points), &
            stat=alloc_status)
      REQUIRE(alloc_status==0)

#ifdef _XRAY_DEBUG
      write(*,*) "RCUT_MAX = ",rcut_max
      write(*,*) "SPHERE_MAX = ",sphere_extent(:)
      write(*,*) 'NUM_ATOMS = ',num_atoms
#endif

      min_usage = 1
      max_usage = 0

      ! loop over all atoms.
      do iatom = 1,num_atoms
         ! Adjust atom scatter coefficients and calulate the cutoff radius:
         
         call get_real_space_coeffs(iatom,acoeffs,nbcoeffs,rcut)

         !write(*,*)'tempFactor=',tempFactor(iatom)
         !write(*,*) 'sharpen=', fft_bfactor_sharpen
         !write(*,*)'acoeffs=',acoeffs
         !write(*,*)'nbcoeffs=',nbcoeffs

         !FIXME: there should be a B-factor limit to avoid this truncation.
         rcut = min(rcut,fft_radius_max) ! Limit to defined maximum
         rcut_squared = rcut**2

         ! Calculate the exact grid coordinate for this atom:
         grid_xyz(:) = frac_to_grid * au_coor(:,iatom)

         ! Find the nearest lattice point:
         grid_point(:) = nint(grid_xyz(:))

         ! Set delta_xyz(:) to the grid - space offset from the lattice point.
         delta_xyz(:) = grid_point(:) - grid_xyz(:)

         ! Calculate the real-space (orthogonal Angstrom) grid point offset.
         offset = matmul(grid_to_orth,delta_xyz)

         ! Determine the minimum grid box to encompass the atomic density sphere:
         sphere_extent(:) = min(grid_size(:), &
               int(rcut * angstrom_to_grid(:) + 0.5_rk_ + REAL_EPSILON))
         num_points = product(sphere_extent(:)*2+1)

         if (num_points == 0) then
            write(stderr,'(A)') 'SPHERE SIZE IS ZERO!'
            write(stderr,*) 'grid_size = ',grid_size(:)
            write(stderr,*) 'rcut = ',rcut
            write(stderr,*) 'angstrom_to_grid = ',angstrom_to_grid(:)
            call mexit(stderr,1)
         end if

         ! --------------------------------------------------------------------
         ! Gather data for vectorized calculations:
         ! Note: an orthogonal cell can utilize some additional shortcuts
         ! in the sum-squared calculation.
         num_points = 0
         do dz = -sphere_extent(3),sphere_extent(3)
            iz = modulo(grid_point(3)+dz,grid_size(3))
            ! This is the Z component of the radius vector, plus the
            ! offset from the atom center to the nearest grid point:
            rvec_z = grid_to_orth(:,3) * real(dz,real_kind) + offset

            do dy = -sphere_extent(2),sphere_extent(2)
               rvec_yz = grid_to_orth(:,2) * real(dy,real_kind) + rvec_z
               ! This check eliminates almost 25% of inner loops
               if (sum(rvec_yz**2) > rcut_squared) cycle
               iy = modulo(grid_point(2)+dy,grid_size(2))

               ! Loop both directions from the center to avoid
               ! looping in areas outside of the cutoff:
               do dir = -1,1,2
                  do dx = max(dir,0), dir*sphere_extent(1), dir
                     rvec_xyz = grid_to_orth(:,1) * real(dx,real_kind) + rvec_yz
                     dist_squared = sum(rvec_xyz**2)
                     if (dist_squared > rcut_squared) exit
                     ix = modulo(grid_point(1)+dx,grid_size(1))
                     num_points = num_points + 1
                     d(num_points) = density_map(ix,iy,iz)
                     if (density_calc) then
                        points(num_points)%p => density_map(ix,iy,iz)
                     end if
                     radius_squared(num_points) = dist_squared
                     if (present(d_XYZ)) then
                        radius_vector(:,num_points) = rvec_xyz
                     end if
                  end do ! dx
               end do ! dir
            end do ! dy
         end do ! dz

         usage = real(num_points)/real(max_points)
         min_usage = min(usage,min_usage)
         max_usage = max(usage,max_usage)

         ! --------------------------------------------------------------------
         ! Vectorized calculations (one pass per coefficient pair):
         do i = 1,num_scatter_coeff
            if (nbcoeffs(i) /= 0.0_rk_) then
               forall (j=1:num_points) e(j) = exp(nbcoeffs(i) * radius_squared(j))
            else
               e(1:num_points) = 1.0_rk_
            end if
            !write(*,*) 'exp=',e(1:num_points:32)

            if (density_calc) then
               ! density + = Acoeff * exp(Bcoeff * radius_squared)
               ! where Bcoeff is negative of standard
               d(1:num_points) = d(1:num_points) + acoeffs(i) * e(1:num_points)
            else
               if (present(d_XYZ)) then
                  do j = 1,3
                     d_XYZ(j,iatom) = d_XYZ(j,iatom) - 2.0_rk_ * &
                           sum(d(1:num_points) * radius_vector(j,1:num_points) &
                           * nbcoeffs(i) * acoeffs(i) * e(1:num_points))
                  end do
               end if
               if (present(d_tempFactor)) then
                  c0 =  1.0_rk_ / (scatter_coeffs(1,i,scatter_type(iatom)) &
                        + occupancy(iatom) - fft_bfactor_sharpen)
                  c1 = -1.5_rk_ * acoeffs(i) * c0
                  c2 =  M_FOURPI2 * acoeffs(i) * c0**2
                  d_tempFactor(iatom) = d_tempFactor(iatom) + &
                        sum( d(1:num_points) &
                       * (c1 + c2 * radius_squared(1:num_points)) * e(1:num_points) )
               end if
               if (present(d_occupancy)) then
                  d_occupancy(iatom) = d_occupancy(iatom) + &
                        sum( d(1:num_points) * acoeffs(i) * e(1:num_points) )
               end if

            ! d(f')
            ! dFprime + = map(xyz) * AEEC * EXP(nbcoeffs(5) * radius_squared)
            ! AEEC = bfactor(iatom) * (SQRT(4.0_rk_ * M_PI/(occupancy(iatom) + BSCAL)))**3
            end if

         end do

         ! --------------------------------------------------------------------
         ! Scatter data for vectorized calculations (only for density_calc):
         if (density_calc) then
            ! (Note: vector notation not allowed with pointer members)
            do i = 1,num_points
               !write(*,*) 'i=',i,', d=',d(i)
               points(i)%p = d(i)
            end do
         end if

      end do !iatom

!FIXME: 
      write(stdout,'(2(A,F8.3))') &
            'Box usage: min = ',min_usage,', max = ',max_usage
      ! Actual sphere/box volume ~ = 0.52

   contains
      subroutine get_real_space_coeffs(iatom,acoeffs,nbcoeffs,rcut)
         integer, intent(in) :: iatom
         real(real_kind), intent(out), dimension(num_scatter_coeff) :: acoeffs, nbcoeffs
         real(real_kind), intent(out) :: rcut
         ! local
         real(real_kind) :: rcut_squared, b
         integer :: i
         ! Calculate the real-space scatter coefficients and the cutoff radius:
         rcut_squared = 0
         do i = 1,num_scatter_coeff
            b = scatter_coeffs(2,i,scatter_type(iatom)) &
                  + tempFactor(iatom) - fft_bfactor_sharpen
            if (b<0.1) then
               stop 'ERROR: adjusted B scatter coefficient too small or negative!' !FIXME
            end if
            ! This is the real-space A coefficient:
            acoeffs(i) = occupancy(iatom) &
                  * scatter_coeffs(1,i,scatter_type(iatom)) &
                  * (sqrt(M_FOURPI/b))**3
            b = M_FOURPI2/b
            ! This is the negative, real-space B coefficient:
            nbcoeffs(i) = -b
            rcut_squared = max(rcut_squared, (fft_density_tolerance + &
                        log(max(abs(acoeffs(i)),ACOEFF_EPSILON)))/b)
         end do
         rcut = sqrt(rcut_squared)

         !FIXME: there should be a B-factor limit to avoid this truncation.
         rcut = min(rcut,fft_radius_max) ! Limit to defined maximum

      end subroutine get_real_space_coeffs

   end subroutine grid_atomic_density

end module xray_real_space_module
