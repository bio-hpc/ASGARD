!+ Poisson Equation potential calculation using fast fourier transforms
!+ See "Outline_Of_pb_fft_and_subroutines.txt" for more information
!+ Used to find coulombic potential field.
!+ Currently only handles poisson equation on uniform dielectric constant
!+ will scale values by epsin.
!+ If needed, can easily be modified to compute helmholtz eqn as well
!+   would only need to add kappa variable to be passed in arguments of pb_fft
!+   and feed comment out the line:
!+              fft_kappa = 0
!+
!+ Usage:
!+               Description                                    |   Arg Name    |   Arg #(s)
!+              --------------                                     ---------        --------
!+      Input:          
!+              -Charge grid                                    | cggrid        |   1
!+              -grid dimensions (# of nodes in each dimension) | nx,ny,nz      |   2-4
!+              -grid spacing                                   | gridh         |   6
!+               (unimplemented place holder, set to 1)
!+              -dielectric constant. Assumed uniform;          | epsin         |   7
!+               set to 1 for vacuum case
!+               usage.
!+      Output:
!+              -electrostatic potential grid                   | phi_grid      |   5



#ifdef FFTW
subroutine pb_fft(cggrid,nx,ny,nz,phi_grid,gridh,epsin,fft_kappa)

use FFTW3
   implicit none

#  include "../include/md.h"
#  include "pb_constants.h"

   ! Passed variables:
   _REAL_ gridh !grid spacing, alias of h 
   _REAL_ cggrid(nx,ny,nz) !alias for charge grid passed in.
    !cggrid is passed directly to fft solver so internal construction will no longer
    !be needed.
   integer nx,ny,nz !Number of grid nodes in x, y, and z direction
   _REAL_ phi_grid(nx,ny,nz) !Holds phi grid data - output
   _REAL_ epsin  !Dielectric constant. Assumed spatially uniform

   ! Local variables:
   _REAL_ epsout !dummy variable both set to epsin, but may be employed in later versions
   _REAL_ fd_phi_grid(nx,ny,nz) !holds values for fd phi grid calculations
   _REAL_ fft_kappa !dummy variable for kappa (salt grid coefficient) unused for now
   integer omega(3) !used to pass harmonic vector to green's funtion

   complex(C_DOUBLE_COMPLEX), pointer :: ft_in(:,:,:) !input buffer for fft transformation of charge / phi grids
   type(C_PTR) :: inptr

   complex(C_DOUBLE_COMPLEX), pointer :: ft_out(:,:,:) !output buffer for fft transformtion of charge / phi grids
   type(C_PTR) :: outptr
   
   type(C_PTR) :: planf !forward fftw transform plan
   type(C_PTR) :: planr !reverse fftw transorm plan

   real*8 gf_con !dummy variable for green's function output
   integer ii,ki,li,mi,xi,yi,zi,xj,yj,zj,cgi,iatom,iatm,jatm !counters

        write(6,'(a)') '-debug: pb_fft: initializing';flush(6)

   inptr = fftw_alloc_complex(int(nx*ny*nz,C_SIZE_T))
   call c_f_pointer(inptr,ft_in,[nx,ny,nz])
   
   outptr = fftw_alloc_complex(int(nx*ny*nz,C_SIZE_T))
   call c_f_pointer(outptr,ft_out,[nx,ny,nz])

   !epsin = 1 !epsin can be used to set the uniform dielectric value now but
              !non-uniform dielectric still not supported
   epsout = epsin
   !fft_kappa = 0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                     BEGIN PHI GRID CALCULATION BLOCKS                       !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
   !++++++++++++++++++++Begin FFT Phi Grid Calculation Block+++++++++++++++++++++!
   !+ generate the fourier transform "plans" for the forward, planf,            +!
   !+ and inverse, planr, fourier transforms, optimized for cggrid              +!
   !+ for now, full complex transforms are implemented.                         +!
   !+ imaginary compenents of potential grid will be discarded by recasting     +!
   !+ as real.                                                                  +!
   !+ Can likely be revised using rfftw3d to perform real to complex            +!
   !+ transform. Will yield some memory savings and possibly boost              +!
   !+ performance slightly. Unclear as to effect this will have on accuracy     +!
   !+ Currently this process is done all in the main block. Depending on how    +!
   !+ the solving process is ultimately organized, it may be necessary to       +!
   !+ encapsulate it in its own subroutine.                                     +!
  ! write(6,*) 'pb_fft.f: building fft plans'; flush(6)
   write(6,'(a)') '-debug: pb_fft: building fft plans';flush(6) 
   planf = fftw_plan_dft_3d(nz,ny,nx,ft_in,ft_out,FFTW_FORWARD,FFTW_ESTIMATE)
   planr = fftw_plan_dft_3d(nz,ny,nx,ft_out,ft_in,FFTW_BACKWARD,FFTW_ESTIMATE)
!#ifdef FFTW2
!   call fftw3d_f77_create_plan(planf,nx,ny,nz,FFTW_FORWARD,FFTW_ESTIMATE)
!   call fftw3d_f77_create_plan(planr,nx,ny,nz,FFTW_BACKWARD,FFTW_ESTIMATE)
!#endif
   !--------------------------Forward Transform Block----------------------------!
   !+ copy the charge grid cggrid into the fourier transform input buffer       +!
   !+ ft_in_buffer = cggrid. Since ft_in_buffer needs to be of type double      +!
   !+ and cggrid is type _REAL_ some typecasting may be needed first to         +!
   !+ get this to run correctly, but fortran appears to handle this itself      +!
   write(6,'(a)') '-debug: pb_fft: copying cggrid to ft_in';flush(6)
   ft_in= cggrid;
   ! apply the fourier transform using planf
   write(6,'(a)') '-debug: pb_fft.f: performing forward transform'; flush(6)
!#ifdef FFTW2
!   call fftwnd_f77_one(planf,ft_in,ft_out);
!#endif
   call fftw_execute_dft(planf,ft_in,ft_out)
   !ft_out = ft_out
   !-----------------------end of forward transform block------------------------!
   !---------------Fourier Space Green Function Transform Block------------------!
   !+ apply (multiply) the greens function to each node in transformed          +!
   !+ charge grid ft_out. We will then have fourier space version               +!
   !+ of potential grid                                                         +!
   write(6,'(a)') '-debug: pb_fft.f: applying green function'; flush(6)
   do ki = 1,nx
      do li = 1,ny
         do mi = 1,nz
            omega(1) = ki
            omega(2) = li
            omega(3) = mi
            call green_fun(gf_con,omega,nx,ny,nz,gridh,epsin,fft_kappa)
            ft_out(ki,li,mi) = ft_out(ki,li,mi) * gf_con;
         end do
      end do
   end do
   !-----------end of fourier space green function transform block---------------!
   !---------------------Reverse Fourier Transform Block-------------------------!
   !+ take the inverse transform using planr to get real-space                  +!
   !+ version of potential grid.                                                +!
   write(6,'(a)') '-debug: pb_fft.f: performing reverse transform'; flush(6)
!#ifdef FFTW2
!   call fftwnd_f77_one(planr,ft_out,ft_in);
!#endif
   call fftw_execute_dft(planr,ft_out,ft_in)
   !+ recast ft_in as real and divide by node volume the former may be          +!
   !+ unnecessary if later revision uses rfftw instead to ensure real           +!
   !+ values. The latter is needed because the fftw package does not            +!
   !+ normalize automatically so this must be done to ensure proper             +!
   !+ scaling                                                                   +!
   write(6,'(a)') '-debug: pb_fft: scaling output'; flush(6)
   phi_grid = real(ft_in) / (nx * ny * nz)*INV_FOURPI/gridh**3; 

  ! write(6,*) 'pb_fft: finished'; flush(6)
   !+ using phi_grid we can now calculate                                       +!
   !+ forces and sum potentials to get pbfrc and eel respectively               +!
   !+ potentials summed by add_atom_potential                                   +!
   !+ forces calculated by calc_force                                           +!
   !+ both energy and force calculations are unimplemented at the moment pending+!
   !+ verification of proper photential grid generation.                        +!
   !------------------end of reverse fourier transform block---------------------!
   !++++++++++++++++++End of FFT Phi Grid Calculation Block++++++++++++++++++++++!
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
   call fftw_free(inptr)
   call fftw_free(outptr)

contains
!+++++++++++++++++Harmonic Green Function Subroutine++++++++++++++++++++++++++++++++
!+ function that calculates the greens function for our differential opperator
!+ in this case we are solving a more general helmholtz equation
!+ which can be used in computing a linearized poisson boltzmann eqn solution
!+ delta phi - (fft_kappa^2)*phi = -4*pi*rho
subroutine green_fun(out,omega,nx,ny,nz,gridh,epsin,fft_kappa)
   
    !returns type double precision to interface with ft_in_buffer better

    !passed variables
    integer omega(3) !k-space harmonics vector, input for green's function
    integer nx,ny,nz
    !_REAL_ meps
    _REAL_ gridh !grid spacing. uniform for now but could later be used for
                 !anisotropicly spaced rectangular grids
    _REAL_ out !output
    _REAL_ fft_kappa !Kappa value for computing helmholtz like equation

    !internal variables
    _REAL_ gm_h !scaling / conversion factor
    _REAL_ gridhx,gridhy,gridhz !place holders for future use
    _REAL_ denomx, denomy, denomz, denom, epsout, epsin
    !in current configuration this subroutine computes k-space version of poisson
    !equation for a periodic system. However, it can easily be modified for use
    !in solving a helmholtz equation with uniform dielectric.
    !epsin = 1    !allows dielectric to be set to any spatialy uniform value 
    !fft_kappa = 0 !this is set to zero now since only poisson equation is considered
                  !if set to nonzero, this solves helmholtz equation; equivalent to
                  !linearized poisson boltzmann equation
    epsout = epsin  !only uniform dielectric constant is considered for now.

    gridhx = gridh;
    gridhy = gridh;
    gridhz = gridh;
    !only uniformly spaced grid is assumed for now. Can easily be modified to
    !support nonuniform grid by passing gridh as an array of 3 reals and updating
    !gridhx, gridhy and gridhz assignment appropriately.

    gm_h = 8.0d0*ACOS(0.0d0)*(gridhx * gridhy * gridhz) ** (2.0d0)

    !meps = 2.0d0 ** (-53d0);

    denomx = (2.0d0 * cos(2.0d0 * PI *(omega(1) - 1)/nx) &
            -(2.0d0 + (fft_kappa * gridhx / epsout) ** 2.0d0))
    denomy = (2.0d0 * cos(2.0d0 * PI *(omega(2) - 1)/ny) &
            -(2.0d0 + (fft_kappa * gridhy / epsout) ** 2.0d0))
    denomz = (2.0d0 * cos(2.0d0 * PI *(omega(3) - 1)/nz) &
            -(2.0d0 + (fft_kappa * gridhz / epsout) ** 2.0d0))
    denom = epsin * (denomx * (gridhy * gridhz) ** 2.0d0 &
            + denomy * (gridhx * gridhz) ** 2.0d0 &
            + denomz * (gridhx * gridhy) ** 2.0d0)
    if (denom == 0.0d0) then !need to avoid divide by zero error.
       out = 0.0d0
    else
       out = -gm_h / denom
    end if

end subroutine green_fun

end subroutine pb_fft
#endif
