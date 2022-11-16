! <compile=optimized>
#include "../include/assert.fh"
module xray_FFT_module
! This module has the FFT versions of routines for calculating X-ray forces.
! Fourier (reciprocal) space values are stored as list of H,K,L values rather
! than a 3D array. Normally, only H,K,L indices with an observed value are
! saved.
!
! These routines all pass F90-style dimensioned arrays.
! This version has several dependent modules commented out using CPP
! conditionals.
!
!  SUBROUTINES:
!
!  FFT_dXYZBQ_dF  --  Calculate the derivative of the coordinate,
!                              B, or occupancy versus the structure factor.
!
!  FFT_Fcalc      --  Calculate structure factors from coordinates.
!
!  CNS_FFTab      --  This reduces the FFT grid to the observerved HKL list.
!                     and therefore defines the choice of reciprocal asymmetric
!                     unit. It was derived from CNS. This is actually
!                     somewhat complicated when symmetry is present, due to
!                     many special-case tests for reflections tha lie
!                     on symmetry planes or axes.
!
   use xray_globals_module
   use constants, only: M_TWOPI => TWOPI
   use xray_reciprocal_space_module
   use xray_real_space_module
   use xray_fftpack_module
   implicit none

   !-------------------------------------------------------------------
   ! Descriptor for the FFT implementation:
   type(fft3c_descriptor), save :: fft3_desc

   !-------------------------------------------------------------------
contains

   !-------------------------------------------------------------------
   ! Thus is a generic front-end to utilize a specific FFT library.
   subroutine FFT_setup()
      integer :: i
      if (spacegroup_number/=1) stop 'ERROR: only P1 symmetry is supported.'
      write(stdout,'(A)') '---------------- XRAY FFT SETUP -------------------'
      write(stdout,'(A,3F9.3,3F7.2)') 'Unit cell: ',unit_cell
      write(stdout,'(A,F9.5)') 'FFT: Requested grid spacing: ',fft_grid_spacing
      grid_size=ceiling(unit_cell(1:3)/fft_grid_spacing)
      write(stdout,'(A,3I8)') 'FFT: Initial grid size: ',grid_size
      do i=1,3
         call fft_adjust_size(grid_size(i),2)
      end do
      write(stdout,'(A,3I8)') 'FFT: Adjusted grid size: ',grid_size
      frac_to_grid = grid_size
      write(stdout,*) frac_to_grid
      angstrom_to_grid = frac_to_grid / unit_cell(1:3)
      write(stdout,*) angstrom_to_grid
      do i=1,3
        orth_to_grid(:,i) = orth_to_frac(:,i) * frac_to_grid(i)
        grid_to_orth(:,i) = frac_to_orth(:,i) / frac_to_grid(i)
      end do
      write(stdout,'(A,3F8.3)') 'FFT: Actual grid spacing: ',angstrom_to_grid**(-1)
      call fft3c_new(fft3_desc,grid_size)
      allocate(density_map(0:grid_size(1)-1, &
                           0:grid_size(2)-1, &
                           0:grid_size(3)))
   end subroutine FFT_setup

   subroutine FFT_forward()
   end subroutine FFT_forward

   subroutine FFT_backward()
   end subroutine FFT_backward

   ! ==========================================================================

   subroutine FFT_Fcalc( &
            !num_hkl,hkl,Fcalc,mSS4,hkl_selected, & ! reflections
            num_hkl,Fcalc,hkl_selected, & ! reflections
      num_atoms,xyz,tempFactor,occupancy,scatter_type_index, & ! atoms
      num_types,num_coeff,coeffs) ! scatter coeff.
      use xray_utils_module, only: allocate_lun, write_map_ezd
      implicit none

      integer, intent(in) :: num_hkl
      ! integer, intent(in) :: hkl(3,num_hkl)                ! H,K,L indices
      complex(real_kind), intent(out) :: Fcalc(num_hkl)    ! Fcalc(h,k,l)
      ! real(real_kind), intent(in) :: mSS4(num_hkl)         ! -(S^2/4)
      integer, intent(in), optional :: hkl_selected(num_hkl)

      integer, intent(in) :: num_atoms
      real(real_kind), intent(in) :: xyz(3,num_atoms) !fractional
      real(real_kind), intent(in) :: tempFactor(num_atoms)
      real(real_kind), intent(in), optional :: occupancy(num_atoms)
      integer, intent(in) :: scatter_type_index(num_atoms)

      integer, intent(in) :: num_coeff, num_types
      real(real_kind), intent(in) :: coeffs(2,num_coeff,num_types)

      ! integer :: i, j, ih,ik,il, unit
      integer :: unit

      ! Calculate the high-res limit
      !FIXME: this calculation is wrong.
      !if (present(hkl_selected)) then
      !  high_res_limit = 1.0_rk_/sqrt( -4.0 * minval(mSS4,mask = hkl_selected))
      !else
      !  high_res_limit = 1.0_rk_/sqrt( -4.0 * minval(mSS4))
      !end if

      density_map = 0

      call grid_atomic_density( &
            num_atoms, xyz, tempFactor, occupancy, scatter_type_index, &
            num_types, num_coeff, coeffs)

      call amopen(allocate_lun(unit),'xray_atom.map','U','F','R')
      call write_map_ezd(unit, &
          reshape((/0,grid_size(1)-1,0,grid_size(2)-1,0,grid_size(3)/),(/2,3/)), &
          density_map)
      stop

      !FIXME: call FFT_forward(FFT)

      call FFT_grid_to_hkl(density_map,num_hkl,hkl_index,Fcalc,hkl_selected)
#if 0
      ! Extract reflections from the FFT that are present in our selected HKL indices.
      do i = 1,num_hkl
         if (present(hkl_selected)) then
            if (hkl_selected(i)==0) cycle
         end if

         ! The FFT calculates the zero - plane and the half where H > 0
         ! Use point symmetry to ensure H >= 0
         ! Point inversion requires converting to the complex conjugate.
         ! NOTE: From testing, the need for complex-conjugate is reversed from
         ! what I expected. (TODO: Look over the math again.)
         ! Negative K,L values wrap at the origin.
         if (hkl(1,i)<0) then
            ih = -hkl(1,i)
            ik = -hkl(2,i); if (ik<0) ik = ik + recip_size(2)
            il = -hkl(3,i); if (il<0) il = il + recip_size(3)
            Fcalc(i) = fft_normalization * cmplx(recip_grid(1,ih,ik,il),recip_grid(1,ih,ik,il)) &
                  * exp( fft_bfactor_sharpen * mSS4(i))
         else
            ih = hkl(1,i)
            ik = hkl(2,i); if (ik<0) ik = ik + recip_size(2)
            il = hkl(3,i); if (il<0) il = il + recip_size(3)
            ! complex-conjugate:
            Fcalc(i) = fft_normalization * cmplx(recip_grid(1,ih,ik,il),-recip_grid(1,ih,ik,il)) &
                  * exp( fft_bfactor_sharpen * mSS4(i))
         end if
      end do
#endif
   end subroutine FFT_Fcalc

   subroutine FFT_dXYZBQ_dF(&
            num_hkl,dF,hkl_selected, & ! reflections
      num_atoms,xyz,tempFactor,occupancy,scatter_type_index, & ! atoms
      num_coeff, num_types, coeffs, & ! scatter coeff.
      dxyz,d_occupancy,d_tempFactor & ! output derivatives
      )   ! TODO: d_aniso_Bij

      implicit none

      ! reciprocal space arrays:
      integer, intent(in) :: num_hkl
      ! integer, intent(in) :: hkl(3,num_hkl)
      complex(real_kind), intent(in) :: dF(num_hkl)
      ! real(real_kind), intent(in) :: mSS4(num_hkl)
      integer, intent(in), optional :: hkl_selected(num_hkl)

      ! coordinate arrays:
      integer, intent(in) :: num_atoms
      real(real_kind), intent(in) :: xyz(3,num_atoms)
      real(real_kind), intent(in) :: tempFactor(num_atoms)
      real(real_kind), intent(in), optional :: occupancy(num_atoms)
      integer, intent(in) :: scatter_type_index(num_atoms)
      real(real_kind), intent(out), optional :: dxyz(3,num_atoms)
      real(real_kind), intent(out), optional :: d_occupancy(num_atoms)
      real(real_kind), intent(out), optional :: d_tempFactor(num_atoms)
      !real(real_kind), intent(out), optional :: d_aniso_Bij(6,num_atoms)

      ! scatter coefficient array:
      integer, intent(in) :: num_coeff, num_types
      real(real_kind), intent(in) :: coeffs(num_coeff,2,num_types)

      ! integer :: ihkl, isym, i, j, ih,ik,il
      ! real(real_kind) :: scale

      ! This array holds the scatter factor for each atom at a given HKL index.
      ! Alternatively, atoms could be sorted by type, but that is likely inefficient,
      ! because sorting is slow and that order would be inefficient for non - bond forces.
      ! -----------------------------------------------------------------------------
      ! Map reflections to the grid, applying symmetry operations:
      call FFT_hkl_to_grid(density_map,num_hkl,hkl_index,dF,hkl_selected)
#if 0
      scale = 1.0_rk_/fft_normalization
      recip_grid = (0.0_rk_,0.0_rk_)
      do isym = 1,num_symmops
         do ihkl = 1,num_hkl
            if (present(hkl_selected)) then
               if (hkl_selected(ihkl)==0) cycle
            end if
            ! Note: H==0 plane requires both +/- symmetry points.
            if (hkl(1,ihkl) <= 0) then
               ih = -hkl(1,ihkl)
               ik = -hkl(2,ihkl); if (ik<0) ik = ik + recip_size(2)
               il = -hkl(3,ihkl); if (il<0) il = il + recip_size(3)
               recip_grid(ih,ik,il) = recip_grid(ih,ik,il) &
                     + dF(ihkl) * scale / exp( fft_bfactor_sharpen * mSS4(ihkl))
            end if
            if (hkl(1,ihkl) >= 0) then
               ik = hkl(1,i)
               ik = hkl(2,i); if (ik<0) ik = ik + recip_size(2)
               il = hkl(3,i); if (il<0) il = il + recip_size(3)
               recip_grid(ih,ik,il) = recip_grid(ih,ik,il) + conjg(dF(ihkl) &
                     * scale / exp( fft_bfactor_sharpen * mSS4(ihkl)))
            end if
         end do
      end do
#endif
      !FIXME: call FFT_backward(FFT)

      call grid_atomic_density( &
            num_atoms, xyz, tempFactor, occupancy, scatter_type_index, &
            num_types, num_coeff, coeffs, &
            d_tempFactor = d_tempFactor, &
            d_occupancy = d_occupancy, &
            d_XYZ = dxyz)

   end subroutine FFT_dXYZBQ_dF

   ! -------------------------------------------------
   ! Converts from HKL reflection list to FFT grid.
   ! Applies symmetry operators.
   subroutine FFT_hkl_to_grid(recip_grid, &
            num_hkl,hkl_index,refl,hkl_selected)

      real(fft_kind) :: &
            recip_grid(2,0:recip_size(1)-1,0:recip_size(2)-1,0:recip_size(3-1))
      integer, intent(in) :: num_hkl
      integer, intent(in) :: hkl_index(3,num_hkl)
      complex(real_kind), intent(in) :: refl(num_hkl)
      integer, intent(in), optional :: hkl_selected(num_hkl)

      real(real_kind) :: scale
      real(real_kind) :: grid_to_phase(3), phase
      integer :: isym, ihkl, h,k,l,hm,km,lm
      real(real_kind) :: rsymm_hkl(3), s(3,3), t(3)
      integer :: symm_hkl(3)
      real(real_kind) :: r
      logical :: Friedel
      complex(real_kind) :: f

      grid_to_phase = M_TWOPI/real(grid_size,real_kind)
      scale = 1.0_rk_/fft_normalization

      recip_grid = 0.0_rk_

      ! loop over all symmetry operators
      do isym = 1,num_symmops
         ! fractional translation converted to phase-shift:
         t = symmop(1:3,4,isym) * M_TWOPI
         ! Reciprocal symmop is the transpose of the real-space symmop-op
         s = transpose(symmop(1:3,1:3,isym))

         do ihkl = 1,num_hkl
            if(present(hkl_selected)) then
               if (hkl_selected(ihkl)==0) cycle
            end if

            rsymm_hkl = matmul(real(hkl_index(1:3,ihkl)),s)
            phase = sum(rsymm_hkl * t) ! FIXME
            f = cmplx(cos(phase),sin(phase))
            ! Convert to array indices (when needed):
            symm_hkl = modulo(nint(rsymm_hkl),grid_size)

            ! Map h,k,l to the H>=0 hemisphere:
            Friedel = (symm_hkl(1) > grid_size(1)/2)
            if (Friedel) symm_hkl = modulo(-symm_hkl,grid_size)

            ! Note: the derivative of a constant factor is conjugate,
            ! as are the factors for the subsequent fourier transform.
            f = f * conjg(refl(ihkl)) &
                  * scale / exp( fft_bfactor_sharpen * mSS4(ihkl))
            if (Friedel) f = conjg(f)
            recip_grid(:,symm_hkl(1),symm_hkl(2),symm_hkl(3)) = &
                  recip_grid(:,symm_hkl(1),symm_hkl(2),symm_hkl(3)) &
                  + (/real(f),aimag(f)/)

         end do !ihkl
      end do !isym

      ! The H == 0 plane contains part of the symmetric hemisphere.
      ! The data must be symmetrically summed and distributed.
      ! From CNS: Note the factor 1/2 in the calling routine.
      do l = 0,grid_size(3)/2
         lm = modulo(l,grid_size(3))
         if (l /= lm) cycle
         do k = 0,grid_size(2)/2
            km = modulo(k,grid_size(2))

            ! add complementary conjugates of 0,k,l and 0,km,lm
            r = recip_grid(1,0,k,l) + recip_grid(1,0,km,lm)
            recip_grid(1,0,k,l) = r
            recip_grid(1,0,km,lm) = r

            r = recip_grid(2,0,k,l) - recip_grid(2,0,km,lm)
            recip_grid(1,0,k,l) = r
            recip_grid(1,0,km,lm) = -r

            do h = 1,merge(grid_size(1)/2, grid_size(1)-1, k == km)
               !recip_grid(h,k,l) = recip_grid(h,k,l) &
               !                  + conjg(recip_grid(grid_size(1)-h,km,lm))
               !recip_grid(grid_size(1)-h,km,lm) = dconjg(recip_grid(h,k,l))
               hm = grid_size(1)-h

               r = recip_grid(1,h,k,l) + recip_grid(1,hm,km,lm)
               recip_grid(1,h,k,l) = r
               recip_grid(1,hm,km,lm) = r

               r = recip_grid(2,h,k,l) - recip_grid(2,hm,km,lm)
               recip_grid(1,h,k,l) = r
               recip_grid(1,hm,km,lm) = -r

            end do
         end do
      end do
   end subroutine FFT_hkl_to_grid

   ! -------------------------------------------------
   ! Converts between FFT grid format and reflection list format.
   ! Applies symmetry operators.
   subroutine FFT_grid_to_hkl(recip_grid, &
            num_hkl,hkl_index,refl,hkl_selected)
      implicit none
      real(fft_kind) :: &
            recip_grid(2,0:recip_size(1)-1,0:recip_size(2)-1,0:recip_size(3-1))
      integer, intent(in) :: num_hkl
      integer, intent(in) :: hkl_index(3,num_hkl)
      complex(real_kind), intent(out) :: refl(num_hkl)
      integer, intent(in), optional :: hkl_selected(num_hkl)

      real(real_kind) :: grid_to_phase(3), phase
      integer :: isym, ihkl
      real(real_kind) :: rsymm_hkl(3), s(3,3), t(3)
      integer :: symm_hkl(3)
      logical :: Friedel
      complex(real_kind) :: f
      real(real_kind) :: r(2)

      grid_to_phase = M_TWOPI/real(grid_size,real_kind)

      ! fill the reflection array with initial indentity-symmop values:
      do ihkl = 1,num_hkl
         if(present(hkl_selected)) then
            if (hkl_selected(ihkl)==0) cycle ! set unselected to zero?
         end if

         rsymm_hkl = matmul(real(hkl_index(1:3,ihkl)),s)
         symm_hkl = modulo(nint(rsymm_hkl),grid_size)

         ! Map h,k,l to the H >= 0 hemisphere:
         Friedel = (symm_hkl(1) > grid_size(1)/2)
         if (Friedel) symm_hkl = modulo(-symm_hkl,grid_size)

         r = recip_grid(:,symm_hkl(1),symm_hkl(2),symm_hkl(3))
         if (Friedel) r(2) = -r(2) !conjg(r)
         refl(ihkl) = fft_normalization * cmplx(r(1),r(2))

      end do !ihkl

      ! loop over all other symmetry operators
      do isym = 2,num_symmops
         t = symmop(1:3,4,isym) * M_TWOPI
         s = transpose(symmop(1:3,1:3,isym))

         do ihkl = 1,num_hkl
            if(present(hkl_selected)) then
               if (hkl_selected(ihkl)==0) cycle
            end if

            ! compute the phase shift from the reciprocal
            ! symmetry index (always zero for P1)
            rsymm_hkl = matmul(real(hkl_index(1:3,ihkl)),s)
            phase = sum(rsymm_hkl * t) ! FIXME!
            f = cmplx(cos(phase),sin(phase))

            ! Convert to array indices (when needed):
            symm_hkl = modulo(nint(rsymm_hkl),grid_size)

            ! Map h,k,l to the H >= 0 hemisphere:
            Friedel = (symm_hkl(1) > grid_size(1)/2)
            if (Friedel) symm_hkl = modulo(-symm_hkl,grid_size)

            r = recip_grid(:,symm_hkl(1),symm_hkl(2),symm_hkl(3))
            if (Friedel) r(2) = -r(2) !conjg(r)
            refl(ihkl) = refl(ihkl) + f * fft_normalization * cmplx(r(1),r(2))

         end do !ihkl
      end do !isym

   end subroutine FFT_grid_to_hkl

   !=====================================================================
   ! Increment NUM to the next value with no factors greater than
   ! FFT_MAX_PRIME, that also contains the factor FACTOR.
   pure subroutine fft_adjust_size (num,factor)
      implicit none
      integer, intent(inout) :: num
      integer, intent(in), optional :: factor ! required factor
      do
         if (mod(num,factor)==0) then
            if (max_factor(num,fft_max_prime)) return
         end if
         num = num + 1
      end do
   end subroutine fft_adjust_size

   !=====================================================================
   ! Return FALSE if any factors of NUM are greater than FFT_MAX_PRIME
   pure function max_factor(num, max) result(ok)
      implicit none
      logical :: ok
      integer, intent(in) :: num, max
      ! local
      integer :: n, p
      ! begin
      ok = .true.
      if (num <= 1) return ! OK (assumes no negative numbers)
      if (max < 2) then
         ok = .false.
         return
      end if
      n = num
      do while (.not.btest(n,0))
         n = n / 2
      end do
      if (n == 1) return
      do p = 3, max, 2
         do while (mod(n, p) == 0)
            n = n / p
         end do
         if (n == 1) return
      end do
      ok = .false.
      return
   end function max_factor
   !=======================================================================

   subroutine test_fft()
      implicit none
      ! Test FFT library to ensure it actually works.
      ! local
      integer :: i, u, v, w, n, alloc_status
      logical :: error
      integer :: grid_size(3) !, rm(1)
      type(fft1c_descriptor) :: cfft
      type(fft1cr_descriptor) :: crfft
      type(fft3c_descriptor) :: cfft3
      logical, parameter :: verbose=.true.
      real(fft_kind) :: err
      ! allocatable
      complex(fft_kind), allocatable, target :: cmap(:,:,:), cmap_orig(:,:,:)
      ! complex/real equivalence
      integer, parameter :: MAXSEQ=256
      complex(fft_kind) :: cseq(0:MAXSEQ-1)
      real(fft_kind) :: rseq(0:MAXSEQ*2-1)
      equivalence(cseq,rseq)
      real(fft_kind) :: rseq_orig(0:MAXSEQ-1)

      ! parameter
      real(fft_kind), parameter :: max_err = 0.01
      real(fft_kind), parameter :: rzero(1) = (/0.0/)
      complex(fft_kind), parameter :: czero(1) = rzero

      ! begin

      ! Set up complex-to-complex transform
      grid_size = (/ 4*3, 3*5, 5*2 /)
      do i = 1, 3
         call fft_adjust_size(grid_size(i))
      end do
      allocate(cmap(grid_size(1),grid_size(2),grid_size(3)), &
            cmap_orig(grid_size(1),grid_size(2),grid_size(3)), &
            stat=alloc_status)
      REQUIRE(alloc_status==0)

      ! FFT3C -> conjugate complex -> FFT3C -> VFYCMAP
      ! Fill CMAP with a pattern:
      do u = 1, grid_size(1)
         do v = 1, grid_size(2)
            do w = 1, grid_size(3)
               cmap(u,v,w) = cmplx( &
                     real(mod(u+v+w,11),fft_kind), &
                     real(mod(u+v+w,13),fft_kind))
            end do
         end do
      end do
      cmap_orig = cmap

      call fft3c_new(cfft3,grid_size,error)
      if (error) then
         write(stderr,'(A)') 'Fatal Error: Call FFT3C_NEW failed.'
         call mexit(stderr,1)
      end if
      call fft3c(cfft3,cmap, error)
      if (error) then
         write(stderr,'(A)') 'Fatal Error: Call FFT3C 1 failed.'
         call mexit(stderr,1)
      end if
      cmap=conjg(cmap)/real(product(grid_size),fft_kind)
      call fft3c(cfft3,cmap, error)
      if (error) then
         write(stderr,'(A)') 'Fatal Error: Call FFT3C 2 failed.'
         call mexit(stderr,1)
      end if

      ! See if we got the right numbers back
      err = max(maxval(abs(real(cmap)-real(cmap_orig))), &
            maxval(abs(aimag(cmap)-aimag(cmap_orig))))
      if (err > max_err) then
         write(stderr,'(A)') 'Fatal Error: Test of FFT3C failed.'
         call mexit(stderr,1)
      end if
      deallocate(cmap,cmap_orig)

      ! Test single precision complex-to-complex-to-real transform
      n = 1
      do
         n = n + 1
         call fft_adjust_size(n,2)
         if (n>MAXSEQ) exit

         !call mfftcr(1, nn, rm)
         ! SFFT1C -> SFFT1CR (with unique half) -> VFYSSEQ
         !!!call isseq(n,cseq)
         do i = 0,n-1
            rseq_orig(i) = real(mod(i+4,7),fft_kind)
         end do
         cseq = rseq_orig

         call fft1c_new(cfft,n,error)
         if (error) then
            write(stderr,'(A)') 'Fatal Error: CALL SFFT1C 1 failed.'
            call mexit(6,1)
         end if
         call fft1c(cfft, cseq)
         call fft1c_delete(cfft)

         call fft1cr_new(crfft,n,error)
         if (error) then
            write(stderr,'(A)') 'Fatal Error: CALL SFFT1CR 1 failed.'
            call mexit(6,1)
         end if
         call fft1cr(crfft, rseq)
         call fft1cr_delete(crfft)

         ! transformed real is scaled and reversed.
         error = any(nint(rseq(n-1:0)/real(n,fft_kind))/=nint(rseq_orig))
         if (error) then
            write(stderr, '(1X,A,I3)') &
                  'Fatal Error: Test of SFFT1C and SFFT1CR failed for N = ',n
            exit
         end if
         if (verbose) then
            write(stdout, '(1X,A,I3)') &
                  'Test of SFFT1C and SFFT1CR succeeded for N = ',n
         end if
      end do
      return
   end subroutine test_fft

end module xray_FFT_module
