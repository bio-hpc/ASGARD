! <compile=optimized>
#include "../include/assert.fh"
module xray_fourier_module
! This module has the fourier (non-FFT) routines for calculating X-ray forces.
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
!  fourier_dXYZBQ_dF  --  Calculate the derivative of the coordinate,
!                              B, or occupancy versus the structure factor.
!                              This one uses the direct (exact) Fourier.
!
!  dTarget_dF         --  Calculate a structure-factor restraint force
!                              from the difference of Fobs and Fcalc. The
!                              residual restraint target-function is a
!                              harmonic restraint in reciprocal space.
!
! FUNCTIONS:
!
! atom_scatter_factor_mss4 --  Calculate the atomic scatter (f) at a specific 
!                              resolution, given the Gaussian coefficients, 
!                              and a modified resolution, defined as -S*S/4.0, 
!                              where S is the reciprocal resolution. (mss4 
!                              stands for "Minus S Squared over 4") See comments
!                              at the beginning of xray_fourier_Fcalc()

   use xray_globals_module
   use constants, only: M_TWOPI => TWOPI
   use xray_reciprocal_space_module
   use xray_real_space_module
   implicit none

   !-------------------------------------------------------------------
contains

   !-------------------------------------------------------------------
   ! Caution: Some literature uses S to represent S^2
   !
   ! (for d*, see p. 93 of Glusker, Lewis, Rossi)
   ! S == d* = sqrt(sum(HKL * orth_to_frac)^2) = sqrt(-4*mSS4)
   ! mSS4 = -S*S/4
   ! mSS4 is more relevant to the formulas used than S.

   subroutine fourier_Fcalc( &
      num_hkl,hkl,Fcalc,mSS4,hkl_selected, & ! reflections
      num_atoms,xyz,tempFactor,occupancy,scatter_type_index)

      implicit none
      integer, intent(in) :: num_hkl
      integer, intent(in) :: hkl(3,num_hkl)
      complex(real_kind), intent(out) :: Fcalc(num_hkl)
      real(real_kind), intent(in) :: mSS4(num_hkl)
      integer, intent(in), optional :: hkl_selected(num_hkl)

      integer, intent(in) :: num_atoms
      real(real_kind), intent(in), target :: xyz(3,num_atoms)
      real(real_kind), intent(in) :: tempFactor(num_atoms)
      real(real_kind), intent(in), target, optional :: occupancy(num_atoms)
      ! index into the scatter - type coeffs table for each atom:
      integer, intent(in), target :: scatter_type_index(num_atoms) 

      ! locals
      integer :: ihkl, i
      ! automatic
      real(real_kind) :: atomic_scatter_factor(num_scatter_types)
      real(real_kind) :: f(num_atoms), angle(num_atoms)

      !  use the following if you expect num_atoms to change during a run:
      ! integer :: alloc_status
      ! allocate(angle(num_atoms), f(num_atoms), stat=alloc_status)
      ! REQUIRE(alloc_status==0)

      do ihkl = 1, num_hkl
         if (present(hkl_selected)) then
            if (hkl_selected(ihkl)==0) cycle
         end if

         ! NOTE: the atomic scatter factors do not change for a given atom type
         ! and hkl index as long as the unit cell is unchanged. If the number 
         ! of atom types is small, it may be worth saving a pre-calculated 
         ! array of num_type by num_reflections.
         ! (this is *not* what is being done here, however)

         do i = 1, num_scatter_types
             atomic_scatter_factor(i) = &
               atom_scatter_factor_mss4(scatter_coefficients(:,:,i),mSS4(ihkl))
         end do

         ! Fhkl = SUM( fj * exp(2 * M_PI * i * (h * xj + k * yj + l * zj)) ),
         !      j = 1,num_selected_atoms
         ! where:
         !    The sum is versus j, over all selected atoms
         !    fj is the atomic scatter for atom j:   atomic_scatter_factor(j)
         !    h,k,l are from the hkl list for index ihkl:   hkl(1:3,ihkl)
         !    x,y,z are coordinates for atom j:   xyz(1:3,j)
         !        xyz(:) may be a reduced list.
         !
         ! Rather than using a complex exponential where the real part is 
         ! always zero, this is optimized to calculate sin and cosine parts,
         ! then convert to a complex number
         ! after the A and B components are summed over all selected atoms.
         ! This can be written as:
         !
         ! Ahkl = SUM( fj * cos(2 * M_PI * (h * xj + k * yj + l * zj)) ), 
         ! Bhkl = SUM( fj * sin(2 * M_PI * (h * xj + k * yj + l * zj)) ),
         !    j = 1,num_selected_atoms

         f(:) = exp( mSS4(ihkl) * tempFactor(:) ) &
              * atomic_scatter_factor(scatter_type_index(:))
         if (present(occupancy)) then
            f(:) = f(:)*occupancy(:)
         endif
         angle(:) = matmul(M_TWOPI * hkl(1:3,ihkl),xyz(1:3,:))
         Fcalc(ihkl) = cmplx( sum(f(:) * sin(angle(:))), &
              sum(f(:) * cos(angle(:))), rk_ )
      end do

   end subroutine fourier_Fcalc

   ! -------------------------------------------------------------------------
   ! Calculate dXYZ, dTempFactor and/or dOccupancy with respect to dF,
   ! using the direct Fourier sum.

   ! Note Q is short for occupancy; B is short for tempFactor (aka B - factor).
   ! Small molecule software uses a U parameter in place of B - factor.
   ! U has a physical meaning, but B - factor is a more natural fit to Fouriers.

   !     isotropic B - factor = 8 * pi ** 2 * isotropic - U
   !     isotropic U = [U(1,1) + U(2,2) + U(3,3)]/3.0

   subroutine fourier_dXYZBQ_dF( &
            num_hkl,hkl,dF,mSS4,hkl_selected, & ! reflections
            num_atoms,xyz,tempFactor,scatter_type_index, & ! atoms
            occupancy, &  ! (put optional args after all non-optional ones)
            dxyz,d_occupancy,d_tempFactor & ! output derivatives
      )   ! TODO: d_aniso_Bij  (anisotropic B refinement)
      implicit none

      ! reciprocal space arrays:
      integer, intent(in) :: num_hkl
      integer, intent(in) :: hkl(3,num_hkl)
      complex(real_kind), intent(in) :: dF(num_hkl)
      real(real_kind), intent(in) :: mSS4(num_hkl)
      integer, intent(in) :: hkl_selected(num_hkl)

      ! coordinate arrays:
      integer, intent(in) :: num_atoms
      real(real_kind), intent(in) :: xyz(3,num_atoms)
      real(real_kind), intent(in) :: tempFactor(num_atoms)
      integer, intent(in) :: scatter_type_index(num_atoms)
      real(real_kind), intent(in), optional :: occupancy(num_atoms)

      ! output derivatives:
      real(real_kind), intent(out), optional :: dxyz(3,num_atoms)
      real(real_kind), intent(out), optional :: d_occupancy(num_atoms)
      real(real_kind), intent(out), optional :: d_tempFactor(num_atoms)
      !real(real_kind), intent(out), optional :: d_aniso_Bij(6,num_atoms)

      ! locals
      integer :: ihkl, iatom, i
      real(real_kind) :: atomic_scatter_factor(num_scatter_types)
      real(real_kind) :: dhkl(3)
      complex(real_kind) :: f
      real(real_kind) :: phase

      dxyz = 0._rk_
      REFLECTION: do ihkl = 1,num_hkl
         !if (present(hkl_selected)) then
            if (hkl_selected(ihkl)==0) cycle REFLECTION
         !end if

         do i = 1,num_scatter_types
            atomic_scatter_factor(i) = &
               atom_scatter_factor_mss4(scatter_coefficients(:,:,i),mSS4(ihkl))
         end do

         ! FIXME: symmetry operations are included here, and require an 
         ! additional loop.
         ! This code is currently limited to P1.
         dhkl = hkl(:,ihkl) * M_TWOPI ! * symmop...

         ! --------------------------------------------------------------------
         ATOM: do iatom = 1,num_atoms
            phase = sum( dhkl * xyz(:,iatom) )
            f = atomic_scatter_factor(scatter_type_index(iatom)) &
                  * exp(mSS4(ihkl) * tempFactor(iatom)) 
   
            f = f * cmplx(sin(phase),cos(phase), rk_)

            if (present(d_occupancy)) then
               d_occupancy(iatom) = d_occupancy(iatom) + &
                 real(f) * real(dF(ihkl)) + aimag(f) * aimag(dF(ihkl))
            end if

            if (present(occupancy)) then
               f = f * occupancy(iatom)
            end if

            if (present(d_tempFactor)) then
               d_tempFactor(iatom) = d_tempFactor(iatom) &
                  + ( real(f) * real(dF(ihkl)) + aimag(f) * aimag(dF(ihkl)) ) &
                  * mSS4(ihkl)
            end if

            if (present(dxyz)) then
               dxyz(:,iatom) = dxyz(:,iatom) + dhkl(:) * &
                   ( aimag(f) * real(dF(ihkl)) - real(f) * aimag(dF(ihkl)) )
            end if

         end do ATOM

      end do REFLECTION

   end subroutine fourier_dXYZBQ_dF

   ! -------------------------------------------------------------------------
   ! R - factor = sum(abs(Fobs - Fcalc)) / sum(Fobs)
   ! This routine computes the force gradient on Fcalc as a harmonic
   ! restraint

   subroutine dTarget_dF(num_hkl,fobs,Fcalc,weight,selected,deriv, &
         residual,xray_energy)
      implicit none
      integer, intent(in) :: num_hkl
      real(real_kind), intent(in) :: fobs(:)
      complex(real_kind), intent(in) :: Fcalc(:)
      integer, intent(in), optional :: selected(:)
      real(real_kind), intent(in), optional :: weight (:)
      complex(real_kind), intent(out), optional :: deriv(:)
      real(real_kind), intent(out), optional :: residual
      real(real_kind), intent(out), optional :: xray_energy

      real(real_kind) :: sum_fo_fc, sum_fo_fo, sum_fc_fc
      real(real_kind) :: abs_Fcalc(num_hkl)
      real(real_kind) :: Fcalc_scale, norm_scale
      real(real_kind), parameter :: F_EPSILON = 1.0e-20_rk_
      integer :: num_selected, i

      abs_Fcalc(:) = abs(Fcalc(:))

#if 0   /* no weights for now!  */
      if (present(weight)) then
         if (present(selected)) then
            sum_fo_fc = sum(weight * fobs * abs_Fcalc,selected/=0)
            sum_fo_fo = sum(weight * fobs ** 2,selected/=0)
            sum_fc_fc = sum(weight * abs_Fcalc ** 2,selected/=0)
         else
            sum_fo_fc = sum(weight * fobs * abs_Fcalc)
            sum_fo_fo = sum(weight * fobs ** 2)
            sum_fc_fc = sum(weight * abs_Fcalc ** 2)
         end if
      else
#endif
         if (present(selected)) then
            sum_fo_fc = sum(fobs * abs_Fcalc,selected/=0)
            sum_fo_fo = sum(fobs ** 2,selected/=0)
            sum_fc_fc = sum(abs_Fcalc ** 2,selected/=0)
            num_selected = 0
            do i=1,num_hkl
              if( selected(i) /= 0 ) num_selected = num_selected + 1
            end do
         else
            sum_fo_fc = sum(fobs * abs_Fcalc)
            sum_fo_fo = sum(fobs ** 2)
            sum_fc_fc = sum(abs_Fcalc ** 2)
         end if
#if 0
      end if
#endif

      if (sum_fc_fc<F_EPSILON .or. sum_fo_fo<F_EPSILON) then
         if (present(deriv)) then
            if (present(selected)) then
               where (selected/=0) deriv = (0.0_rk_,0.0_rk_)
            else
               deriv = (0.0_rk_,0.0_rk_)
            end if
         end if
         if (present(residual)) residual = -1.0_rk_
         return
      end if

      Fcalc_scale = sum_fo_fc / sum_fc_fc
      norm_scale = 1.0_rk_ / sum_fo_fo

      ! Note: when Fcalc is approximately zero the phase is undefined, 
      ! so no force can be determined even if the energy is high. (Similar 
      ! to a linear bond angle.)

      ! We get the phase from Fcalc, so it needs to be divided by abs(Fcalc)
      ! to become a unit vector.
      !
      ! deriv = Fcalc/abs(Fcalc) * K * ( Fobs - Fcalc_scale * abs(Fcalc) )
      !       = Fcalc * K * ( Fobs/abs(Fcalc) - Fcalc_scale )
      !      ....plus the terms arising from differentiating Fcalc_scale with
      !      respect to Fcalc

      if (present(deriv)) then
#if 0  /* no weights for now!  */
         if (present(weight)) then
            if (present(selected)) then
               where (selected/=0)
               where (abs_Fcalc > F_EPSILON)
               deriv = Fcalc * ( - 2.0_rk_ * Fcalc_scale * &
                     weight * norm_scale * ( fobs / abs_Fcalc - Fcalc_scale ) )
            else where
               deriv = (0.0_rk_,0.0_rk_)
               end where
               end where
            else ! no selected
               where (abs_Fcalc > F_EPSILON)
               deriv = Fcalc * ( - 2.0_rk_ * Fcalc_scale * &
                     weight * norm_scale * ( fobs / abs_Fcalc - Fcalc_scale ) )
            else where
               deriv = (0.0_rk_,0.0_rk_)
               end where
            end if
         else ! no weight
#endif
            if (present(selected)) then
               where (selected/=0)
                  where (abs_Fcalc > F_EPSILON)
                     deriv(:) = - 2.0_rk_ * Fcalc(:) * norm_scale * &
                        ( fobs(:) - Fcalc_scale*abs_Fcalc(:) ) *  &
#if 1
                        ( Fcalc_scale/abs_Fcalc(:) )  !Joe's version
#else
                        ! add in terms from derivative of Fcalc_scale:
                        ( Fcalc_scale/abs_Fcalc(:) + &
                          fobs(:)/sum_fc_fc - 2._rk_*abs_Fcalc(:) * &
                          sum_fo_fc/(sum_fc_fc**2) )
#endif
                  else where
                     deriv = (0.0_rk_,0.0_rk_)
                  end where
                  else where
                     deriv = (0.0_rk_,0.0_rk_)
               end where
            else ! no selected (and no weight)
               where (abs_Fcalc > F_EPSILON)
                  deriv = (0.0_rk_,0.0_rk_)  !FIXME
               else where
                  deriv = (0.0_rk_,0.0_rk_)
               end where
            end if
         end if
#if 0
      end if
#endif

      ! Do the appropriate R - factor calculation, depending on weights 
      ! and/or a selection mask.

      if (present(residual)) then
         if (present(weight)) then
            if (present(selected)) then
               residual = sum(weight * abs(fobs - Fcalc_scale * abs_Fcalc), &
                 selected/=0) / sum(fobs,selected/=0)
            else
               residual = sum (weight * abs( fobs - Fcalc_scale * abs_Fcalc )) &
                / sum(fobs)
            end if
         else
            if (present(selected)) then
               residual = sum(abs(fobs - Fcalc_scale * abs_Fcalc), &
                  selected/=0) / sum(fobs, selected/=0)
               xray_energy = norm_scale * &
                 sum((fobs - Fcalc_scale * abs_Fcalc)**2, selected/=0)
            else
               residual = sum (abs(fobs - Fcalc_scale * abs_Fcalc)) / sum(fobs)
            end if
         end if
      end if

   end subroutine dTarget_dF

   function atom_scatter_factor_mss4(coeffs,mss4) result(sfac)
      real(real_kind) :: sfac
      real(real_kind), intent(in) :: coeffs(2,scatter_ncoeffs), mss4
      sfac = coeffs(1,scatter_ncoeffs) + &
             sum( coeffs(1,1:scatter_ncoeffs-1) &
             * exp(mss4*coeffs(2,1:scatter_ncoeffs-1)))
   end function atom_scatter_factor_mss4

   !dac addition:
   subroutine get_mss4(num_hkl,hkl_index,mSS4)
      integer, intent(in) :: num_hkl
      integer, intent(in) :: hkl_index(3,num_hkl)
      real(real_kind), intent(inout) :: mss4(:)
      integer :: ihkl, h,k,l
      real(real_kind) :: a,b,c,alpha,beta,gamma,V,S2
      real(real_kind) :: sina,cosa,sinb,cosb,sing,cosg
      real(real_kind) :: astar,bstar,cstar,cosas,cosbs,cosgs

      a = unit_cell(1)
      b = unit_cell(2)
      c = unit_cell(3)
      alpha = 3.14159*unit_cell(4)/180.
      beta = 3.14159*unit_cell(5)/180.
      gamma = 3.14159*unit_cell(6)/180.
      sina = sin(alpha)
      cosa = cos(alpha)
      sinb = sin(beta)
      cosb = cos(beta)
      sing = sin(gamma)
      cosg = cos(gamma)

      V = a*b*c*sqrt( 1.d0 - cosa**2 - cosb**2 -cosg**2 &
          + 2.d0*cosa*cosb*cosg )
      astar = b*c*sina/V
      bstar = a*c*sinb/V
      cstar = b*a*sing/V
      cosas = (cosb*cosg - cosa)/(sinb*sing)
      cosbs = (cosa*cosg - cosb)/(sina*sing)
      cosgs = (cosb*cosa - cosg)/(sinb*sina)

      ! work from p. 93 of Glusker, Lewis, Rossi:

      do ihkl=1,num_hkl
        h = hkl_index(1,ihkl)
        k = hkl_index(2,ihkl)
        l = hkl_index(3,ihkl)
        S2 = (h*astar)**2 + (k*bstar)**2 + (l*cstar)**2  &
           + 2*k*l*bstar*cstar*cosas + 2*l*h*cstar*astar*cosbs  &
           + 2*h*k*astar*bstar*cosgs 
        mSS4(ihkl) = -S2/4.
      end do
   end subroutine get_mss4

end module xray_fourier_module
