#include "copyright.h"
#include "../include/dprec.fh"

module pb_lsolver

   implicit none

#  include "pb_constants.h"

   integer l_xm, l_ym, l_zm, l_xmym, l_xmymzm
   integer l_itn, l_maxitn, l_bcopt
   integer mg_nlevel, ncyc_before, ncyc_after
   _REAL_ l_fmiccg, l_wsor, l_fmiccg2
   _REAL_ l_norm, l_inorm, l_accept, l_epsout, l_pbkappa, l_h

   ! All

   integer,allocatable :: pl_ind(:,:)
   integer,allocatable :: pl_ind3d(:,:,:)
   _REAL_,allocatable :: l_am1(:), l_am2(:), l_am3(:)
   _REAL_,allocatable :: l_am4(:), l_am5(:), l_am6(:)
   _REAL_,allocatable :: l_bv(:), l_rd(:), l_ad(:)
   _REAL_,allocatable :: l_pv(:), l_tv(:), l_zv(:)

   ! PICCG

   _REAL_,allocatable :: pl_am1(:,:,:), pl_am2(:,:,:), pl_am3(:,:,:)
   _REAL_,allocatable :: pl_am4(:,:,:), pl_am5(:,:,:), pl_am6(:,:,:)
   _REAL_,allocatable :: pl_bv(:,:,:), pl_ad(:,:,:), pl_rd(:,:,:)
   _REAL_,allocatable :: pl_pv(:,:,:), pl_tv(:,:,:), pl_zv(:,:,:)

   ! MG

   integer gid
   integer,allocatable :: mg_index(:), mg_index_ext(:),mg_x_idx(:),mg_size(:,:)
   _REAL_,allocatable :: mg_onorm(:)
   _REAL_,allocatable :: l_rv(:), l_iv(:), l_bz(:), l_xv(:)
   _REAL_,allocatable :: l_scratch_am1(:), l_scratch_am2(:), l_scratch_am3(:)
   _REAL_,allocatable :: l_scratch_vf(:), l_scratch_vc(:), l_scratch_bz(:)

contains

!===========================================================================

subroutine init_param( nx,ny,nz,nxny,nxnynz,&
                       p_maxitn,p_bcopt,p_fmiccg,p_fmiccg2,p_accept,p_pbkappa,&
                       p_epsout,p_h,p_wsor )

   implicit none

   integer nx,ny,nz,nxny,nxnynz,p_maxitn,p_bcopt
   _REAL_ p_fmiccg,p_accept,p_epsout,p_pbkappa,p_wsor,p_h,p_fmiccg2

   l_xm = nx
   l_ym = ny
   l_zm = nz
   l_xmym = nxny
   l_xmymzm = nxnynz
   l_maxitn = p_maxitn
   l_bcopt = p_bcopt
   l_accept = p_accept

   ! ICCG

   l_fmiccg = p_fmiccg
   l_fmiccg2= p_fmiccg2

   ! MG

   mg_nlevel = 4
   ncyc_before = 10
   ncyc_after = 10
   l_pbkappa = p_pbkappa
   l_epsout = p_epsout
   l_h = p_h

   ! SOR

   l_wsor = p_wsor


end subroutine

!===========================================================================

subroutine allocate_array( solvopt )

   implicit none

   integer solvopt

   integer l,m,n,i

   select case (solvopt)
   case (1)
      if ( l_bcopt /= 10 ) then
         allocate( l_ad(1:l_xmymzm+l_xmym),l_am1(1-l_xmym:l_xmymzm+l_xmym))
         allocate(l_am2(1-l_xmym:l_xmymzm+l_xmym),l_am3(1-l_xmym:l_xmymzm+l_xmym))
         allocate( l_rd(1-l_xmym:l_xmymzm),l_bv(1-l_xmym:l_xmymzm))
         allocate( l_tv(1-l_xmym:l_xmymzm+l_xmym),l_zv(1-l_xmym:l_xmymzm+l_xmym))
         allocate( l_pv(1-l_xmym:l_xmymzm+l_xmym))
      else
         allocate( pl_bv(0:l_xm,0:l_ym,0:l_zm),pl_tv(1:l_xm+1,1:l_ym+1,1:l_zm+1) )
         allocate( pl_zv(0:l_xm,0:l_ym,0:l_zm),pl_pv(0:l_xm,0:l_ym,0:l_zm) )
         allocate( pl_ad(0:l_xm,0:l_ym,0:l_zm),pl_rd(0:l_xm,0:l_ym,0:l_zm) )
         allocate( pl_am1(0:l_xm,0:l_ym,0:l_zm),pl_am2(0:l_xm,0:l_ym,0:l_zm) )
         allocate( pl_am3(0:l_xm,0:l_ym,0:l_zm),pl_am4(0:l_xm,0:l_ym,0:l_zm) )
         allocate( pl_am5(0:l_xm,0:l_ym,0:l_zm),pl_am6(0:l_xm,0:l_ym,0:l_zm) )
         allocate( pl_ind(1:l_xmymzm,3) )
      end if
   case (2)
      allocate ( mg_index(1:mg_nlevel+1),mg_index_ext(1:mg_nlevel+1))
      allocate ( mg_x_idx(1:mg_nlevel+1),mg_size(1:3,1:mg_nlevel) )
      allocate ( mg_onorm(1:mg_nlevel) )

      mg_index_ext(1) = 1
      mg_index(1) = 1
      mg_x_idx(1) = 1
      mg_size(1,1) = l_xm
      mg_size(2,1) = l_ym
      mg_size(3,1) = l_zm
      m = l_xmymzm
      l = m + l_xmym
      n = l + l_xmym

      allocate( l_zv(1:m) )
      if ( l_bcopt == 10 ) then
         allocate( l_scratch_vc( 1-(l_xm+3)*(l_ym+3):(l_xm+3)*(l_ym+3)*(l_zm+3)+(l_xm+3)*(l_ym+3)) )
         allocate( l_scratch_vf( 1-(l_xm+3)*(l_ym+3):(l_xm+3)*(l_ym+3)*(l_zm+3)+(l_xm+3)*(l_ym+3)) )
         allocate( l_scratch_bz( 1-(l_xm+3)*(l_ym+3):(l_xm+3)*(l_ym+3)*(l_zm+3)+(l_xm+3)*(l_ym+3)) )
         allocate( l_scratch_am1(1-(l_xm+3)*(l_ym+3):(l_xm+3)*(l_ym+3)*(l_zm+3)) )
         allocate( l_scratch_am2(1-(l_xm+3)*(l_ym+3):(l_xm+3)*(l_ym+3)*(l_zm+3)) )
         allocate( l_scratch_am3(1-(l_xm+3)*(l_ym+3):(l_xm+3)*(l_ym+3)*(l_zm+3)) )
      end if

      do i = 2, mg_nlevel
         mg_index_ext(i) = 1 + l
         mg_index(i) = 1 + m
         mg_x_idx(i) = 1 + n
         if ( l_bcopt == 10 ) then
            mg_size(1:3,i) = ( mg_size(1:3,i-1) ) / 2 ! even no. grid points
         else
            mg_size(1:3,i) = ( mg_size(1:3,i-1) - 1 ) / 2 ! odd no. grid points
         end if
         m = m + mg_size(1,i) * mg_size(2,i) * mg_size(3,i)
         l = l + mg_size(1,i) * mg_size(2,i) * mg_size(3,i) + mg_size(1,i) * mg_size(2,i)
         n = n + mg_size(1,i) * mg_size(2,i) * mg_size(3,i) + 2 * mg_size(1,i) * mg_size(2,i)
      end do
      mg_index_ext(i) = 1 + l
      mg_index(i) = 1 + m
      mg_x_idx(i) = 1 + n

      allocate ( l_ad(1:m), l_bv(1:m), l_rv(1:m), l_iv(1:m), l_bz(1:m) )
      allocate ( l_am1(1:l), l_am2(1:l), l_am3(1:l) )
      allocate ( l_xv(1:n) )
      if ( l_bcopt == 10 ) then
         allocate( pl_ind3d(0:(sum((2+mg_size(1,1:mg_nlevel)))),0:(sum((2+mg_size(2,1:mg_nlevel)))),&
                  0:(sum((2+mg_size(3,1:mg_nlevel)))) ) )
      end if
   case (3)
      allocate(l_ad(1:l_xmymzm),l_am1(1-l_xmym:l_xmymzm),l_am2(1-l_xmym:l_xmymzm),l_am3(1-l_xmym:l_xmymzm))
      allocate(l_bv(1:l_xmymzm),l_pv(1-l_xmym:l_xmymzm+l_xmym),l_zv(1:l_xmymzm))
      if (l_bcopt==10) then 
         allocate(l_am4(1-l_xmym:l_xmymzm),l_am5(1-l_xmym:l_xmymzm),l_am6(1-l_xmym:l_xmymzm))
      end if
   case (4)
      allocate(l_ad(1:l_xmymzm),l_am1(1-l_xmym:l_xmymzm),l_am2(1-l_xmym:l_xmymzm),l_am3(1-l_xmym:l_xmymzm))
      allocate(l_bv(1:l_xmymzm),l_zv(1:l_xmymzm))
      if (l_bcopt==10) then 
         allocate(pl_ind3d(0:(l_xm+1),0:(l_ym+1),0:(l_zm+1)))
         allocate(pl_ind(1:l_xmymzm,1:3))
      end if
   end select

end subroutine allocate_array

!===========================================================================

subroutine deallocate_array(solvopt)

   implicit none
   integer solvopt

   select case (solvopt)
   case (1)
      if ( l_bcopt /= 10 ) then
         deallocate( l_ad,l_am1 )
         deallocate( l_am2,l_am3 )
         deallocate( l_rd,l_bv )
         deallocate( l_tv,l_zv )
         deallocate( l_pv )
      else
         deallocate( pl_bv,pl_tv )
         deallocate( pl_zv,pl_pv )
         deallocate( pl_ad,pl_rd )
         deallocate( pl_am1,pl_am2 )
         deallocate( pl_am3,pl_am4 )
         deallocate( pl_am5,pl_am6 )
         deallocate( pl_ind )
      end if
   case (2)
      deallocate( mg_index, mg_index_ext )
      deallocate( mg_x_idx,mg_size )
      deallocate( mg_onorm )
      deallocate( l_zv )
      deallocate( l_ad, l_bv, l_rv, l_iv, l_bz )
      deallocate( l_am1, l_am2, l_am3 )
      deallocate( l_xv )
      if (l_bcopt == 10) then
         deallocate( l_scratch_vc )
         deallocate( l_scratch_vf )
         deallocate( l_scratch_bz )
         deallocate( l_scratch_am1 )
         deallocate( l_scratch_am2 )
         deallocate( l_scratch_am3 )
         deallocate( pl_ind3d )
      end if
   case (3)
      deallocate( l_ad,l_am1,l_am2,l_am3 )
      if (l_bcopt == 10) then 
         deallocate( l_am4,l_am5,l_am6 )
      end if
      deallocate( l_bv,l_pv,l_zv )
   case (4)
      deallocate( l_ad,l_am1,l_am2,l_am3)
      deallocate( l_bv,l_zv )
      if ( l_bcopt == 10 ) then
         deallocate( pl_ind )
         deallocate( pl_ind3d )
      end if
   end select

end subroutine deallocate_array

!===========================================================================

subroutine init_array( solvopt,epsx,epsy,epsz,p_bv,p_iv,p_xs )
   
   implicit none

   integer solvopt
   _REAL_ epsx(*),epsy(*),epsz(*)
   _REAL_ p_bv(1:l_xmymzm),p_iv(1:l_xmymzm)
   _REAL_ p_xs(1-l_xmym:l_xmymzm+l_xmym)

   integer lxmym,l,m,n,i,j,k
   integer ii, xi, yi, zi
   integer xsi, ysi, zsi, xn, yn, zn
   _REAL_,allocatable :: lepsx(:), lepsy(:), lepsz(:) 
   _REAL_ lfactor

   select case ( solvopt )

   case (1) ! ICCG

      ! non-periodic systems

      if ( l_bcopt /= 10 ) then

      ! set up A bands

         call feedepsintoam( l_xm, l_ym, l_zm, &
                             l_am1(1:l_xmymzm), l_am2(1:l_xmymzm), l_am3(1:l_xmymzm), &
                             epsx, epsy, epsz )
         l_am1(l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
         l_am2(l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
         l_am3(l_xmymzm+1:l_xmymzm+l_xmym) = ZERO
         call feedepsintoad( l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz )
         if ( l_pbkappa /= ZERO ) then
            lfactor = l_epsout*(l_h*l_pbkappa)**2
            l_ad(1:l_xmymzm) = lfactor*p_iv(1:l_xmymzm) + l_ad(1:l_xmymzm)
         end if

      ! for non-periodic grids, set up padding elements of A bands

         l_am1(1-l_xmym:0) = ZERO
         l_am2(1-l_xmym:0) = ZERO
         l_am3(1-l_xmym:0) = ZERO
         call setupper(l_am1(1),l_am2(1),l_am3(1))

      ! more padding elements for other arrays and initialization

         l_ad (l_xmymzm+1:l_xmymzm+l_xmym) = ZERO

         l_bv (1-l_xmym:0) = ZERO
         l_bv (1:l_xmymzm) = p_bv(1:l_xmymzm)
         l_pv (1-l_xmym:0) = ZERO
         l_pv (l_xmymzm+1:l_xmymzm+l_xmym) = ZERO

         l_tv = ZERO
         l_zv = ZERO
         l_rd = ONE

      ! periodic systems

      else

         ! initialize all to zero first

         pl_ad = ZERO
         pl_am1 = ZERO; pl_am2 = ZERO; pl_am3 = ZERO 
         pl_am4 = ZERO; pl_am5 = ZERO; pl_am6 = ZERO
         pl_bv = ZERO
         pl_tv = ZERO
         pl_zv = ZERO
         pl_pv = ZERO
         pl_rd = ZERO

         ! set the 1D to 3D indexing array, pl_ind

         do k = 1, l_zm; do j = 1, l_ym; do i = 1, l_xm
            ii = i+l_xm*(j-1+l_ym*(k-1))
            pl_ind(ii,1) = i; pl_ind(ii,2) = j; pl_ind(ii,3) = k
         end do; end do; end do

         ! set b vector

         call loadbv_piccg( pl_bv, p_bv(1:l_xmymzm) )

         ! set A bands

         call feedepsintoam_piccg( l_xm, l_ym, l_zm, pl_am1, pl_am2, pl_am3, &
                                   epsx, epsy, epsz )
         call feedepsintoad_piccg( l_xm, l_ym, l_zm, pl_ad, &
                                   epsx, epsy, epsz, &
                                   p_iv )
         call setoutter_piccg( pl_am1, pl_am2, pl_am3, pl_am4, pl_am5, pl_am6 )
      end if

   case (3) ! CG

      ! set up A bands

      call feedepsintoam( l_xm, l_ym, l_zm, &
                          l_am1(1:l_xmymzm), l_am2(1:l_xmymzm), l_am3(1:l_xmymzm), &
                          epsx, epsy, epsz )
      if ( l_pbkappa == ZERO ) then
         call feedepsintoad( l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz )
      else
         lfactor = l_epsout*(l_h*l_pbkappa)**2
         call feedepsintoad( l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz )
         l_ad(1:l_xmymzm) = l_ad(1:l_xmymzm) + lfactor*p_iv(1:l_xmymzm)
      end if
       
      ! padding the lower faces.

      l_am1(1-l_xmym:0) = ZERO
      l_am2(1-l_xmym:0) = ZERO
      l_am3(1-l_xmym:0) = ZERO

      ! padding the upper faces.

      if ( l_bcopt == 10 ) then
         call setupper2( l_am1(1),l_am2(1),l_am3(1), &
                            l_am4(1),l_am5(1),l_am6(1) )
      else
         call setupper ( l_am1(1),l_am2(1),l_am3(1) )
      end if

      l_bv(1:l_xmymzm) = p_bv(1:l_xmymzm)
      l_pv(1-l_xmym:0) = ZERO
      l_pv(l_xmymzm+1:l_xmymzm+l_xmym) = ZERO

   case (4) ! SOR

      ! set up A bands

      call feedepsintoam( l_xm, l_ym, l_zm, &
                          l_am1(1:l_xmymzm), l_am2(1:l_xmymzm), l_am3(1:l_xmymzm), &
                          epsx, epsy, epsz )
      if ( l_pbkappa == ZERO ) then
         call feedepsintoad( l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz )
      else
         lfactor = l_epsout*(l_h*l_pbkappa)**2
         call feedepsintoad( l_xm, l_ym, l_zm, l_ad, epsx, epsy, epsz )
         l_ad(1:l_xmymzm) = lfactor*p_iv(1:l_xmymzm) + l_ad(1:l_xmymzm)
      end if

      ! for periodic grids, use indexing arrays from 1d to 3d and 3d to 1d
      ! pl_ind is for interior points only

      if ( l_bcopt == 10 ) then
         do zi = 0, (l_zm+1); do yi = 0, (l_ym+1); do xi = 0, (l_xm+1)
            ii = 1 + mod(xi-1+l_xm,l_xm) + l_xm*mod(yi-1+l_ym,l_ym) + l_xm*l_ym*mod(zi-1+l_zm,l_zm)
            pl_ind3d(xi,yi,zi) = ii
            if ( xi > 0 .and. xi <= l_xm .and. yi > 0 .and. yi <= l_ym .and. zi >0 .and. zi<=l_zm ) then
               pl_ind(ii,1:3) = (/xi,yi,zi/)
            end if
         end do; end do; end do

      ! for non-periodic grids, set up padding elements

      else
         ! the lower faces.

         l_am1(1-l_xmym:0) = ZERO
         l_am2(1-l_xmym:0) = ZERO
         l_am3(1-l_xmym:0) = ZERO

         ! the upper faces.

         call setupper(l_am1(1), l_am2(1), l_am3(1))
      end if

      l_bv(1:l_xmymzm) = p_bv(1:l_xmymzm)

   case (2) ! MG

      l_ad  = ZERO
      l_am1 = ZERO
      l_am2 = ZERO
      l_am3 = ZERO
      l_iv  = ZERO
      l_bv  = ZERO
      l_rv  = ZERO
      l_zv  = ZERO
      l_xv  = ZERO
      l_xv(1+l_xmym:l_xmymzm+l_xmym) = p_xs(1:l_xmymzm)
      l_bv(1:l_xmymzm) = p_bv(1:l_xmymzm)

      m = 0
      do i = 1, mg_nlevel
         m = m + (mg_size(1,i)+1) * (mg_size(2,i)+1) * (mg_size(3,i)+1)
      end do
      allocate( lepsx(1:m), lepsy(1:m), lepsz(1:m) )
      call feedepsintoam( l_xm, l_ym, l_zm, lepsx(1:l_xmymzm),  &
                          lepsy(1:l_xmymzm), lepsz(1:l_xmymzm), &
                                              epsx, epsy, epsz )

      lfactor = l_epsout*(l_h*l_pbkappa)**2
      l_iv(1:l_xmymzm) = p_iv(1:l_xmymzm)

      i = 1
      m = mg_index(i)
      n = mg_index_ext(i)
      lxmym = mg_size(1,i)*mg_size(2,i)
      call set_am_ad(lepsx(m),lepsy(m),lepsz(m),l_iv(m), &
                     l_am1(n+lxmym), l_am2(n+lxmym), l_am3(n+lxmym), &
                     l_ad(m), l_bz(m), &
                     mg_size(1,i),mg_size(2,i),mg_size(3,i),lfactor,l_epsout)

      do i = 2, mg_nlevel
         l = mg_index(i-1)
         m = mg_index(i)
         n = mg_index_ext(i)
         call restrict_eps_map(lepsx(l),lepsy(l),lepsz(l),mg_size(1,i-1),mg_size(2,i-1),mg_size(3,i-1), &
                               lepsx(m),lepsy(m),lepsz(m),mg_size(1,i),mg_size(2,i),mg_size(3,i))
         call restrict_iv(l_iv(l),mg_size(1,i-1),mg_size(2,i-1),mg_size(3,i-1),&
                          l_iv(m),mg_size(1,i),mg_size(2,i),mg_size(3,i) )
         lfactor = lfactor * 4
         lxmym = mg_size(1,i)*mg_size(2,i)

         call set_am_ad(lepsx(m),lepsy(m),lepsz(m),l_iv(m), &
                        l_am1(n+lxmym), l_am2(n+lxmym), l_am3(n+lxmym), &
                        l_ad(m), l_bz(m), &
                        mg_size(1,i),mg_size(2,i),mg_size(3,i),lfactor,l_epsout)
      end do
      deallocate(lepsx,lepsy,lepsz)

      ! set up indexing array for the periodic version

      if ( l_bcopt == 10 ) then
         i = 1
         xn = mg_size(1,1); yn = mg_size(2,1); zn = mg_size(3,1)
         xsi = 0; ysi = 0; zsi = 0
         do zi = 0, (zn+1); do yi = 0, (yn+1); do xi = 0, (xn+1)
            ii = 1+mod(xi+xn-1,xn)+xn*mod(yi-1+yn,yn)+xn*yn*mod(zi-1+zn,zn)
            pl_ind3d(xi+xsi,yi+ysi,zi+zsi) = ii
         end do; end do; end do

         do i = 2, mg_nlevel
            xn = mg_size(1,i); yn=mg_size(2,i); zn=mg_size(3,i)
            xsi = sum((2+mg_size(1,1:(i-1)))); ysi=sum((2+mg_size(2,1:(i-1)))); zsi=sum((2+mg_size(3,1:(i-1))))
            do zi=0,(zn+1); do yi=0,(yn+1); do xi=0,(xn+1)
               ii=1+mod(xi+xn-1,xn)+xn*mod(yi-1+yn,yn)+xn*yn*mod(zi-1+zn,zn)
               pl_ind3d(xi+xsi,yi+ysi,zi+zsi) = ii
            end do; end do; end do
         end do
      end if

   end select

contains

subroutine feedepsintoam(xm, ym, zm, am1, am2, am3, eps1, eps2, eps3)

   implicit none

   integer xm, ym, zm
   _REAL_ am1(1:xm,1:ym,1:zm)
   _REAL_ am2(1:xm,1:ym,1:zm)
   _REAL_ am3(1:xm,1:ym,1:zm)
   _REAL_ eps1(0:xm,1:ym,1:zm)
   _REAL_ eps2(1:xm,0:ym,1:zm)
   _REAL_ eps3(1:xm,1:ym,0:zm)

   am1(1:xm,1:ym,1:zm) = eps1(1:xm,1:ym,1:zm)
   am2(1:xm,1:ym,1:zm) = eps2(1:xm,1:ym,1:zm)
   am3(1:xm,1:ym,1:zm) = eps3(1:xm,1:ym,1:zm)


end subroutine feedepsintoam

subroutine feedepsintoad(xm, ym, zm, ad, eps1, eps2, eps3)

   implicit none

   integer xm, ym, zm
   _REAL_ ad(1:xm,1:ym,1:zm)
   _REAL_ eps1(0:xm,1:ym,1:zm)
   _REAL_ eps2(1:xm,0:ym,1:zm)
   _REAL_ eps3(1:xm,1:ym,0:zm)

   ad(1:xm,1:ym,1:zm) =                    eps1(1:xm,  1:ym,  1:zm)
   ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps1(0:xm-1,1:ym,  1:zm)
   ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps2(1:xm,  1:ym,  1:zm)
   ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps2(1:xm,  0:ym-1,1:zm)
   ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps3(1:xm,  1:ym,  1:zm)
   ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps3(1:xm,  1:ym,  0:zm-1)


end subroutine feedepsintoad

subroutine feedepsintoam_piccg(xm, ym, zm, am1, am2, am3, eps1, eps2, eps3)

   implicit none

   integer xm,ym,zm
   _REAL_ am1(0:xm,0:ym,0:zm)
   _REAL_ am2(0:xm,0:ym,0:zm)
   _REAL_ am3(0:xm,0:ym,0:zm)
   _REAL_ eps1(0:xm,1:ym,1:zm)
   _REAL_ eps2(1:xm,0:ym,1:zm)
   _REAL_ eps3(1:xm,1:ym,0:zm)

   am1(1:xm,1:ym,1:zm) = eps1(1:xm,1:ym,1:zm)
   am2(1:xm,1:ym,1:zm) = eps2(1:xm,1:ym,1:zm)
   am3(1:xm,1:ym,1:zm) = eps3(1:xm,1:ym,1:zm)


end subroutine feedepsintoam_piccg

subroutine feedepsintoad_piccg(xm, ym, zm, ad, eps1, eps2, eps3, piv)

   implicit none

   integer xm,ym,zm
   _REAL_ ad(0:xm,0:ym,0:zm)
   _REAL_ eps1(0:xm,1:ym,1:zm)
   _REAL_ eps2(1:xm,0:ym,1:zm)
   _REAL_ eps3(1:xm,1:ym,0:zm)
   _REAL_ piv(1:xm,1:ym,1:zm)

   _REAL_ lfactor

    ad(1:xm,1:ym,1:zm) =                    eps1(1:xm,  1:ym,  1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps1(0:xm-1,1:ym,  1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps2(1:xm,  1:ym,  1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps2(1:xm,  0:ym-1,1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps3(1:xm,  1:ym,  1:zm)
    ad(1:xm,1:ym,1:zm) = ad(1:xm,1:ym,1:zm)+eps3(1:xm,  1:ym,  0:zm-1)
    if ( l_pbkappa /= ZERO ) then
       lfactor = l_epsout*(l_h*l_pbkappa)**2
       ad(1:xm,1:ym,1:zm) = piv(1:xm,1:ym,1:zm)*lfactor + ad(1:xm,1:ym,1:zm)
    end if


end subroutine feedepsintoad_piccg

subroutine loadbv_piccg(lbv,pbv)

   implicit none

   _REAL_ lbv(0:l_xm,0:l_ym,0:l_zm)
   _REAL_ pbv(1:l_xm,1:l_ym,1:l_zm)

   integer i,j,k,ii

   lbv(1:l_xm,1:l_ym,1:l_zm) = pbv(1:l_xm,1:l_ym,1:l_zm)


end subroutine loadbv_piccg

subroutine setoutter_piccg(l_am1,l_am2,l_am3,l_am4,l_am5,l_am6)

   implicit none

   _REAL_  l_am1(0:l_xm,0:l_ym,0:l_zm)
   _REAL_  l_am2(0:l_xm,0:l_ym,0:l_zm)
   _REAL_  l_am3(0:l_xm,0:l_ym,0:l_zm)
   _REAL_  l_am4(0:l_xm,0:l_ym,0:l_zm)
   _REAL_  l_am5(0:l_xm,0:l_ym,0:l_zm)
   _REAL_  l_am6(0:l_xm,0:l_ym,0:l_zm)

   integer xi,yi,zi

   do yi = 1, l_ym; do zi = 1, l_zm
      l_am4(1   ,yi,zi) = l_am1(l_xm,yi,zi)
      l_am1(l_xm,yi,zi) = 0
   end do; end do

   do xi = 1, l_xm; do zi = 1, l_zm
      l_am5(xi,1   ,zi) = l_am2(xi,l_ym,zi)
      l_am2(xi,l_ym,zi) = 0
   end do; end do

   do xi = 1, l_xm; do yi = 1, l_ym
      l_am6(xi,yi,1   ) = l_am3(xi,yi,l_zm)
      l_am3(xi,yi,l_zm) = 0
   end do; end do

end subroutine setoutter_piccg

subroutine setupper( l_am1, l_am2, l_am3 )

   implicit none

   _REAL_ l_am1(l_xm,l_ym,l_zm), l_am2(l_xm,l_ym,l_zm), l_am3(l_xm,l_ym,l_zm)

   integer i,j,k

   do j = 1, l_ym; do k = 1, l_zm
      l_am1(l_xm,j,k) = ZERO
   end do; end do
   do i = 1, l_xm; do k = 1, l_zm
      l_am2(i,l_ym,k) = ZERO
   end do; end do
   do i = 1, l_xm; do j = 1, l_ym
      l_am3(i,j,l_zm) = ZERO
   end do; end do

end subroutine setupper

subroutine setupper2( l_am1, l_am2, l_am3, l_am4, l_am5, l_am6 )

   implicit none

   _REAL_ l_am1(l_xm,l_ym,l_zm), l_am2(l_xm,l_ym,l_zm), l_am3(l_xm,l_ym,l_zm)
   _REAL_ l_am4(l_xm,l_ym,l_zm), l_am5(l_xm,l_ym,l_zm), l_am6(l_xm,l_ym,l_zm)
   integer i, j, k

   ! setting x faces';flush(6)

   do j = 1, l_ym; do k = 1, l_zm
      l_am4(1,j,k) = l_am1(l_xm,j,k)
      l_am1(l_xm,j,k) = ZERO
   end do; end do

   ! setting y faces';flush(6)

   do i = 1, l_xm; do k = 1, l_zm
      l_am5(i,1,k) = l_am2(i,l_ym,k)
      l_am2(i,l_ym,k) = ZERO
   end do; end do

   ! setting z faces';flush(6)

   do i = 1, l_xm; do j = 1, l_ym
      l_am6(i,j,1) = l_am3(i,j,l_zm)
      l_am3(i,j,l_zm) = ZERO
   end do; end do

end subroutine setupper2


end subroutine init_array
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! pb_iccg
!
! ICGG core routine for linearized FDPB equation, non-periodic version.  
! This version uses 1-d padded arrays for performance.
!
! Authors:
! Ray Luo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine pb_iccg( phi, xs )

   implicit none

   ! Passed variables

   ! phi is solution array for future uses in energy and force calculation.
   ! xs is the internal scaled solution that also provides an initial guess.

   _REAL_ phi(l_xmymzm), xs(1-l_xmym:l_xmymzm+l_xmym)

   ! Local variables

   logical uconvg
   integer i, j
   _REAL_ alpha, beta, pdotz, bdotb1, bdotb2

   ! initialization

   ! WMBS - This modifies the diagonal row (d(i)^-1) using
   !        the formula from Luo R., David L., and Gilson M. K. J.C.P. 2002
   !        see page 1246
   !        l_rd stores the reciprical of the diagonal matrix D
   !        from the M=(L+D)D^-1(D+U) preconditioner splitting

   do i = 1, l_xmymzm
      l_rd(i) = ONE/( l_ad(i) - &
         l_am1(i-1     )*(         l_am1(i-1     )+l_fmiccg*l_am2(i-1     )+l_fmiccg*l_am3(i-1     ))*l_rd(i-1     ) - &
         l_am2(i-l_xm  )*(l_fmiccg*l_am1(i-l_xm  )+         l_am2(i-l_xm  )+l_fmiccg*l_am3(i-l_xm  ))*l_rd(i-l_xm  ) - &
         l_am3(i-l_xmym)*(l_fmiccg*l_am1(i-l_xmym)+l_fmiccg*l_am2(i-l_xmym)+         l_am3(i-l_xmym))*l_rd(i-l_xmym) )
   end do

   do i = 1, l_xmymzm
      ! use l_rd to modify the l_ad terms. In the unmodified this would have been
      ! skipped since the diagonal would be unity.
      l_ad(i) = l_ad(i)*l_rd(i)
      l_rd(i) = sqrt(l_rd(i))
      l_bv(i) = l_bv(i)*l_rd(i) ! first step to precondition initial guess in l_bv (diagonal mult.)

      ! WMBS - L matrix sub-diagonal term building
      l_am1(i-1     ) = l_am1(i-1     )*l_rd(i)*l_rd(i-1     ) ! L sub-diagonal of previous x terms
      l_am2(i-l_xm  ) = l_am2(i-l_xm  )*l_rd(i)*l_rd(i-l_xm  ) ! L sub-diagonal of previous y terms
      l_am3(i-l_xmym) = l_am3(i-l_xmym)*l_rd(i)*l_rd(i-l_xmym) ! L sub-diagonal of previous z terms
   end do


   l_inorm = ZERO
   do i = 1, l_xmymzm 
      l_ad(i) = l_ad(i) - TWO ! this would yield -1 in unmodified version

      ! WMBS - apply L matrix to l_bv.
      l_bv(i) = l_bv(i) + l_am1(i-1     )*l_bv(i-1     ) &
                        + l_am2(i-l_xm  )*l_bv(i-l_xm  ) &
                        + l_am3(i-l_xmym)*l_bv(i-l_xmym)
      l_inorm = l_inorm + abs(l_bv(i))
   end do

   do i = l_xmymzm, 1, -1
      ! WMBS - l_tv is the transpose of the lower triangular
      !        this loop is applying (D+U) part of M to initial guess xs
      !        l_tv(i) = 0 when i > l_xmymzm or i < 1 here 
      l_tv(i) = xs(i) + l_am1(i     )*l_tv(i+1     ) &
                      + l_am2(i     )*l_tv(i+l_xm  ) &
                      + l_am3(i     )*l_tv(i+l_xmym)
   end do

   do i = 1, l_xmymzm
      l_zv(i) = xs(i) + l_ad (i       )*l_tv(i       ) &
                      + l_am1(i-1     )*l_zv(i-1     ) &
                      + l_am2(i-l_xm  )*l_zv(i-l_xm  ) &
                      + l_am3(i-l_xmym)*l_zv(i-l_xmym)
   end do

   bdotb1 = ZERO
   do i = l_xmymzm, 1, -1
      l_zv(i) = l_zv(i) + l_tv(i)
      l_bv(i) = l_bv(i) - l_zv(i)

      ! iteration 0.

      bdotb1 = bdotb1 + l_bv(i)*l_bv(i)
      l_pv(i)  = l_bv(i)

      ! first step of the matrix vector multiplication, see below
      ! WMBS - apply (D+U) part of M 
      l_tv(i) = l_pv(i) + l_am1(i     )*l_tv(i+1   ) &
                    + l_am2(i     )*l_tv(i+l_xm  ) &
                    + l_am3(i     )*l_tv(i+l_xmym)
   end do

   l_itn = 0
   uconvg = .true.

   ! main loop

   do while ( uconvg )

      ! second and third steps of the matrix vector multiplication

      pdotz = ZERO
      do i = 1, l_xmymzm+l_xmym
      ! WMBS - Apply (L+D)D^1 part of M to goes to l_xmymzm+l_xmym to accomodate
      !        (D+U)*p matrix vector multiplication
         l_zv(i) = l_pv(i) + l_ad (i       )*l_tv(i       ) &
                           + l_am1(i-1     )*l_zv(i-1     ) &
                           + l_am2(i-l_xm  )*l_zv(i-l_xm  ) &
                           + l_am3(i-l_xmym)*l_zv(i-l_xmym)

         j = i - l_xmym
         l_zv(j) = l_zv(j) + l_tv(j)

         pdotz = pdotz + l_pv(j)*l_zv(j)
      end do
      alpha = bdotb1/pdotz

      l_norm = ZERO
      bdotb2 = ZERO
      l_itn = l_itn + 1
      do i = 1, l_xmymzm
         xs(i) = xs(i) + alpha*l_pv(i)
         l_bv(i) = l_bv(i) - alpha*l_zv(i)
         l_norm = l_norm  + abs(l_bv(i))

         bdotb2 = bdotb2 + l_bv(i)*l_bv(i)
      end do

      ! check convergence

      if ( l_itn >= l_maxitn .or. l_norm <= l_accept*l_inorm ) then

         uconvg = .false.
         if ( l_itn >= l_maxitn ) then
            write(6,'(a)') 'PB Warning in pb_iccg(): maxitn exceeded!'
         end if

      else

         beta = bdotb2/bdotb1
         bdotb1 = bdotb2

         ! first step of the matrix vector multiplication

         do i = l_xmymzm, 1, -1
            l_pv(i) = l_bv(i) + beta*l_pv(i)

            l_tv(i) = l_pv(i) + l_am1(i)*l_tv(i+1     ) &
                              + l_am2(i)*l_tv(i+l_xm  ) &
                              + l_am3(i)*l_tv(i+l_xmym)
         end do
      end if
   end do  !  while ( uconvg )

   ! back scaling of the solution

   do i = l_xmymzm, 1, -1
      l_tv(i)  = xs(i) + l_am1(i)*l_tv(i+1     ) &
                       + l_am2(i)*l_tv(i+l_xm  ) &
                       + l_am3(i)*l_tv(i+l_xmym)

      phi(i) = l_tv(i)*l_rd(i)
   end do


end subroutine pb_iccg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! pb_piccg
!
! ICGG core routine for linearized FDPB equation, periodic version.  
! This version uses 3-d arrays for ease of implementation.
!
! Authors:
! Wesley Smith
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine pb_piccg( phi, xs )

   implicit none

   ! Passed variables

   _REAL_ phi(1:l_xm,1:l_ym,1:l_zm), xs(0:l_xm,0:l_ym,0:l_zm)

   ! Local variables

   logical uconvg
   integer ii
   integer xi, yi, zi, xp, yp, zp
   _REAL_ alpha, beta, pdotz, bdotb1, bdotb2

   ! initialization
 
   call pb_piccg_initialize(pl_ind,pl_rd,pl_ad,pl_am1,pl_am2,pl_am3,&
           pl_am4,pl_am5,pl_am6,pl_bv,pl_tv,pl_pv,pl_zv,            &
           alpha,bdotb1,bdotb2,pdotz,beta,l_inorm,uconvg,           &
           xs)

   l_itn=0
   uconvg = .true.

   ! main loop

   do while ( uconvg )

      ! second and third steps of the matrix vector multiplication

      pdotz = ZERO
      do ii = 1, l_xmymzm-l_xmym
         xi = pl_ind(ii,1); yi = pl_ind(ii,2); zi = pl_ind(ii,3)
         xp = mod(xi,l_xm)+1; yp = mod(yi,l_ym)+1; zp = mod(zi,l_zm)+1
         pl_zv(xi,yi,zi) = pl_pv (xi,yi,zi  ) + pl_ad (xi,yi,zi  )*pl_tv(xi,yi,zi  )&
                                              + pl_am1(xi-1,yi,zi)*pl_zv(xi-1,yi,zi)&
                                              + pl_am2(xi,yi-1,zi)*pl_zv(xi,yi-1,zi)&
                                              + pl_am3(xi,yi,zi-1)*pl_zv(xi,yi,zi-1)&
                                              + pl_am4(xp,yi,zi  )*pl_zv(xp,yi,zi  )&
                                              + pl_am5(xi,yp,zi  )*pl_zv(xi,yp,zi  )&
                                              + pl_am6(xi,yi,zp  )*pl_zv(xi,yi,zp  )    
      end do

      ! we are not done yet ...

      zi = l_zm
      do xi = 1, l_xm; do yi = 1, l_ym
         xp = mod(xi,l_xm)+1; yp = mod(yi,l_ym)+1; zp = mod(zi,l_zm)+1
         pl_zv(xi,yi,zi) = pl_pv (xi,yi,zi  ) + pl_ad (xi,yi,zi  )*pl_tv(xi,yi,zi  )&
                                              + pl_am1(xi-1,yi,zi)*pl_zv(xi-1,yi,zi)&
                                              + pl_am2(xi,yi-1,zi)*pl_zv(xi,yi-1,zi)&
                                              + pl_am3(xi,yi,zi-1)*pl_zv(xi,yi,zi-1)&
                                              + pl_am4(xp,yi,zi  )*pl_zv(xp,yi,zi  )&
                                              + pl_am5(xi,yp,zi  )*pl_zv(xi,yp,zi  )&
                                              + pl_am6(xi,yi,zp  )*pl_zv(xi,yi,zp  )

         pl_zv(xi,yi,1) = pl_zv(xi,yi,1) + pl_tv(xi,yi,1)
         pdotz = pdotz + pl_pv(xi,yi,1) * pl_zv(xi,yi,1)
      end do; end do

      pl_zv(1:l_xm,1:l_ym,2:l_zm) = pl_zv(1:l_xm,1:l_ym,2:l_zm) + pl_tv(1:l_xm,1:l_ym,2:l_zm)
      pdotz = pdotz + sum(pl_pv(1:l_xm,1:l_ym,2:l_zm)*pl_zv(1:l_xm,1:l_ym,2:l_zm))

      alpha = bdotb1/pdotz
      l_norm = ZERO
      bdotb2 = ZERO
      l_itn = l_itn + 1

      do ii = 1, l_xmymzm
         xi = pl_ind(ii,1); yi = pl_ind(ii,2); zi = pl_ind(ii,3)
         xs(xi,yi,zi) = xs(xi,yi,zi) + alpha*pl_pv(xi,yi,zi)
         pl_bv(xi,yi,zi) = pl_bv(xi,yi,zi) - alpha*pl_zv(xi,yi,zi)
         l_norm = l_norm + abs(pl_bv(xi,yi,zi))

         bdotb2 = bdotb2 + pl_bv(xi,yi,zi)**2
      end do

      ! check convergence

      if ( l_itn >= l_maxitn .or. l_norm <= l_accept*l_inorm ) then
         uconvg = .false.
         if ( l_itn >= l_maxitn ) then
            write(6,'(a)') 'PB Warning in pb_piccg(): maxitn exceeded!'
         end if
      else
         beta = bdotb2/bdotb1
         bdotb1 = bdotb2

         ! first step of the matrix vector multiplication
         ! WMBS-  need to update tv for next loop as in initialization.
         !        periodicity can't be done ahead since tv needs an updated pv which is tangled in
         !        the same loop. We will soak extra zero addition steps instead of building bv ahead
         !        of time in extra loop over all grids

         pl_pv(1:l_xm,1:l_ym,1:l_zm) = pl_bv(1:l_xm,1:l_ym,1:l_zm) + beta*pl_pv(1:l_xm,1:l_ym,1:l_zm)
         do ii = l_xmymzm, 1, -1
            xi = pl_ind(ii,1); yi = pl_ind(ii,2); zi = pl_ind(ii,3)
            xp = mod(xi+l_xm-2,l_xm)+1; yp = mod(yi+l_ym-2,l_ym)+1; zp = mod(zi+l_zm-2,l_zm)+1 
            pl_tv(xi,yi,zi) = pl_pv (xi,yi,zi) + pl_am1(xi,yi,zi)*pl_tv(xi+1,yi,zi) &
                                               + pl_am2(xi,yi,zi)*pl_tv(xi,yi+1,zi) &
                                               + pl_am3(xi,yi,zi)*pl_tv(xi,yi,zi+1) &
                                               + pl_am4(xi,yi,zi)*pl_tv(xp,yi,zi  ) &
                                               + pl_am5(xi,yi,zi)*pl_tv(xi,yp,zi  ) &
                                               + pl_am6(xi,yi,zi)*pl_tv(xi,yi,zp  )
         end do
      end if

   end do  !  while ( uconvg )

   ! back scaling of the solution

   do ii = l_xmymzm, 1, -1
      xi = pl_ind(ii,1); yi = pl_ind(ii,2); zi = pl_ind(ii,3)   
      xp = mod(xi+l_xm-2,l_xm)+1; yp = mod(yi+l_ym-2,l_ym)+1; zp = mod(zi+l_zm-2,l_zm)+1 
      pl_tv(xi,yi,zi) = xs(xi,yi,zi) + pl_am1(xi,yi,zi)*pl_tv(xi+1,yi,zi) &
                                     + pl_am2(xi,yi,zi)*pl_tv(xi,yi+1,zi) &
                                     + pl_am3(xi,yi,zi)*pl_tv(xi,yi,zi+1) &
                                     + pl_am4(xi,yi,zi)*pl_tv(xp,yi,zi  ) &
                                     + pl_am5(xi,yi,zi)*pl_tv(xi,yp,zi  ) &
                                     + pl_am6(xi,yi,zi)*pl_tv(xi,yi,zp  )

      phi(xi,yi,zi) = pl_tv(xi,yi,zi)*pl_rd(xi,yi,zi)
   end do

   ! periodic systems have a non-trivial null space. Most notably, the null space will
   ! include the uniform distribution. In order to attain consistent results, the solution
   ! must be purged of its content by subratcting the projection of the soultion distribution
   ! onto the uniform distribution. This is particularly important for non-charge neutral
   ! systems to which a uniform neutralizing charge was added in order to solve the system.

   phi(1:l_xm,1:l_ym, 1:l_zm) = phi(1:l_xm, 1:l_ym, 1:l_zm) - sum(phi(1:l_xm,1:l_ym,1:l_zm)) / l_xmymzm

contains

subroutine pb_piccg_initialize(lind,lrd,lad,lam1,lam2,lam3,lam4,lam5,lam6,&
            lbv,ltv,lpv,lzv,alpha,bdotb1,bdotb2,pdotz,beta,inorm,unconv,&
            lxs)

   implicit none

   integer lind(1:l_xmymzm,3)
   _REAL_ alpha, bdotb1, bdotb2, pdotz, beta, inorm
   _REAL_ lrd (0:l_xm, 0:l_ym, 0:l_zm), lad (0:l_xm  ,0:l_ym  ,0:l_zm  )
   _REAL_ lam1(0:l_xm, 0:l_ym, 0:l_zm), lam2(0:l_xm  ,0:l_ym  ,0:l_zm  )
   _REAL_ lam3(0:l_xm, 0:l_ym, 0:l_zm), lam4(0:l_xm  ,0:l_ym  ,0:l_zm  )
   _REAL_ lam5(0:l_xm, 0:l_ym, 0:l_zm), lam6(0:l_xm  ,0:l_ym  ,0:l_zm  )
   _REAL_ lbv (0:l_xm, 0:l_ym, 0:l_zm), lpv (0:l_xm  ,0:l_ym  ,0:l_zm  )
   _REAL_ lzv (0:l_xm, 0:l_ym, 0:l_zm), ltv (1:l_xm+1,1:l_ym+1,1:l_zm+1)
   _REAL_ lxs (0:l_xm, 0:l_ym, 0:l_zm)
 
   logical unconv
   integer xi,yi,zi,xp,yp,zp,i,j,k,ii
   _REAL_ netcrg

   ! neutralize grid charge in periodic Poisson systems
   ! salt would automatically do it in the periodic PB systems.

   if ( l_pbkappa == ZERO ) then
      netcrg = -sum(lbv(1:l_xm,1:l_ym,1:l_zm))/l_xmymzm
      if ( abs(netcrg) > 1.0d-9 ) then
         write(6,'(a)') 'PB Warning: System net charge detected and removed for bcopt=10'
         lbv(1:l_xm,1:l_ym,1:l_zm) = lbv(1:l_xm,1:l_ym,1:l_zm) + netcrg           
      end if 
   end if

   lrd(1:l_xm,1:l_ym,1:l_zm) = ONE

   ! first step in setup of l_rd working array
   ! WMBS - This seems to modify the diagonal row (d(i)^-1) using:
   !        the formula from Luo R., David L., and Gilson M. K. J.C.P. 2002
   !        see page 1246
   !        l_rd seems to store the reciprical of the diagonal matrix D
   !        from the M=(L+D)D^-1(D+U) preconditioner splitting

   do ii = 1, l_xmymzm
      xi = lind(ii,1); yi = lind(ii,2); zi = lind(ii,3)
      xp = mod(xi,l_xm)+1; yp = mod(yi,l_ym)+1; zp=mod(zi,l_zm)+1 
      lrd(xi,yi,zi) = ONE/( lad(xi,yi,zi)-                                    &  
                     lam1(xi-1,yi,zi)*(lam1(xi-1,yi,zi)     + l_fmiccg   *(   & 
                        lam2(xi-1,yi,zi)+ lam3(xi-1,yi,zi)) + l_fmiccg2  *(   &
                        lam4(xi-1,yi,zi)+ lam5(xi-1,yi,zi)  + lam6(xi-1,yi,zi)&
                        ))*lrd(xi-1,yi,zi) -                                  &
                     lam2(xi,yi-1,zi)*(lam2(xi,yi-1,zi)     + l_fmiccg   *(   &
                        lam1(xi,yi-1,zi)+ lam3(xi,yi-1,zi)) + l_fmiccg2  *(   &
                        lam4(xi,yi-1,zi)+ lam5(xi,yi-1,zi)  + lam6(xi,yi-1,zi)&
                        ))*lrd(xi,yi-1,zi) -                                  &
                     lam3(xi,yi,zi-1)*(lam3(xi,yi,zi-1)     + l_fmiccg   *(   &
                        lam1(xi,yi,zi-1)+ lam2(xi,yi,zi-1)) + l_fmiccg2  *(   &
                        lam4(xi,yi,zi-1)+ lam5(xi,yi,zi-1)  + lam6(xi,yi,zi-1)&
                        ))*lrd(xi,yi,zi-1) -                                  &
                     lam4(xp,yi,zi)*(lam4(xp,yi,zi)         + l_fmiccg2  *(   &
                        lam5(xp,yi,zi)+ lam6(xp,yi,zi))     + l_fmiccg   *(   &
                        lam1(xp,yi,zi)+ lam2(xp,yi,zi)      + lam3(xp,yi,zi)  &
                        ))*lrd(xp,yi,zi) -                                    &
                     lam5(xi,yp,zi)*(lam5(xi,yp,zi)         + l_fmiccg2  *(   &
                        lam4(xi,yp,zi)+ lam6(xi,yp,zi))     + l_fmiccg   *(   &
                        lam1(xi,yp,zi)+ lam2(xi,yp,zi)      + lam3(xi,yp,zi)  &
                        ))*lrd(xi,yp,zi) -                                    &
                     lam6(xi,yi,zp)*(lam6(xi,yi,zp)         + l_fmiccg2  *(   &
                        lam4(xi,yi,zp)+ lam5(xi,yi,zp))     + l_fmiccg   *(   &
                        lam1(xi,yi,zp)+ lam2(xi,yi,zp)      + lam3(xi,yi,zp)  &
                        ))*lrd(xi,yi,zp)                                      &
                     )
   end do

   ! WMBS - L matrix sub-diagonal term building

   do ii = 1, l_xmymzm
      xi = lind(ii,1); yi = lind(ii,2); zi = lind(ii,3)
      lad(xi,yi,zi) = lad(xi,yi,zi)*lrd(xi,yi,zi)
      lrd(xi,yi,zi) = sqrt(lrd(xi,yi,zi))
      lbv(xi,yi,zi) = lbv(xi,yi,zi)*lrd(xi,yi,zi)
      lam1(xi-1,yi,zi) = lam1(xi-1,yi,zi)*lrd(xi,yi,zi)*lrd(xi-1,yi,zi)
      lam2(xi,yi-1,zi) = lam2(xi,yi-1,zi)*lrd(xi,yi,zi)*lrd(xi,yi-1,zi)
      lam3(xi,yi,zi-1) = lam3(xi,yi,zi-1)*lrd(xi,yi,zi)*lrd(xi,yi,zi-1)
   end do

   do yi = 1, l_ym; do zi = 1, l_zm
      lam4(1,yi,zi) = lam4(1,yi,zi)*lrd(1,yi,zi)*lrd(l_xm,yi,zi)
   end do; end do
   do xi = 1, l_xm; do zi = 1, l_zm
      lam5(xi,1,zi) = lam5(xi,1,zi)*lrd(xi,1,zi)*lrd(xi,l_ym,zi)
   end do; end do
   do xi = 1, l_xm; do yi = 1, l_ym
      lam6(xi,yi,1) = lam6(xi,yi,1)*lrd(xi,yi,1)*lrd(xi,yi,l_zm)
   end do; end do

   ! Continue preconditioning of initial guess - l_bv

   inorm = ZERO

   do ii = 1, l_xmymzm
      xi = lind(ii,1); yi = lind(ii,2); zi = lind(ii,3)
      xp = mod(xi,l_xm)+1; yp = mod(yi,l_ym)+1; zp = mod(zi,l_zm)+1
      lad(xi,yi,zi) = lad(xi,yi,zi) - TWO
      lbv(xi,yi,zi) = lbv(xi,yi,zi)                         &
                        +  lam1(xi-1,yi,zi)*lbv(xi-1,yi,zi) &
                        +  lam2(xi,yi-1,zi)*lbv(xi,yi-1,zi) &
                        +  lam3(xi,yi,zi-1)*lbv(xi,yi,zi-1) &
                        +  lam4(xp,yi,zi  )*lbv(xp,yi,zi  ) &
                        +  lam5(xi,yp,zi  )*lbv(xi,yp,zi  ) &
                        +  lam6(xi,yi,zp  )*lbv(xi,yi,zp  ) 
      inorm = inorm + ABS(lbv(xi,yi,zi))
   end do

   ! U*p explicit conditioning

   do ii=l_xmymzm, 1, -1
      xi = lind(ii,1); yi = lind(ii,2); zi = lind(ii,3)
      xp = mod(xi+l_xm-2,l_xm)+1
      yp = mod(yi+l_ym-2,l_ym)+1
      zp = mod(zi+l_zm-2,l_zm)+1
      ltv(xi,yi,zi) = lxs(xi,yi,zi) + lam1(xi,yi,zi)*ltv(xi+1,yi,zi) &
                                    + lam2(xi,yi,zi)*ltv(xi,yi+1,zi) &
                                    + lam3(xi,yi,zi)*ltv(xi,yi,zi+1) &
                                    + lam4(xi,yi,zi)*ltv(xp,yi,zi  ) &
                                    + lam5(xi,yi,zi)*ltv(xi,yp,zi  ) &
                                    + lam6(xi,yi,zi)*ltv(xi,yi,zp  ) 
   end do   

   ! WMBS - loops for applying the L and D matrices
   ! L*p explicit conditioning

   do ii = 1, l_xmymzm
      xi = lind(ii,1); yi = lind(ii,2); zi = lind(ii,3)
      xp = mod(xi,l_xm)+1; yp = mod(yi,l_ym)+1; zp = mod(zi,l_zm)+1
      lzv(xi,yi,zi) = lxs(xi,yi,zi) + lad(xi,yi,zi   )*ltv(xi,yi,zi  )  &
                                    + lam1(xi-1,yi,zi)*lzv(xi-1,yi,zi)  &
                                    + lam2(xi,yi-1,zi)*lzv(xi,yi-1,zi)  &
                                    + lam3(xi,yi,zi-1)*lzv(xi,yi,zi-1)  &
                                    + lam4(xp,yi,zi  )*lzv(xp,yi,zi  )  &
                                    + lam5(xi,yp,zi  )*lzv(xi,yp,zi  )  &
                                    + lam6(xi,yi,zp  )*lzv(xi,yi,zp  )      
   end do

   ! WMBS - need to add periodicity for l_tv. This could be a bit
   !        if we use additional l_am4 - l_am6 with nonzero terms
   !        only on boundaries, we can simply extend l_tv setting
   !        line

   bdotb1 = ZERO
   do ii = l_xmymzm, 1, -1
      xi = lind(ii,1); yi = lind(ii,2); zi = lind(ii,3)
      xp = mod(xi+l_xm-2,l_xm)+1 !wrap x=1 to x=l_xm
      yp = mod(yi+l_ym-2,l_ym)+1 !wrap y=1 to y=l_ym
      zp = mod(zi+l_zm-2,l_zm)+1 !wrap z=1 to z=l_zm

      lzv(xi,yi,zi) = lzv(xi,yi,zi) + ltv(xi,yi,zi)
      lbv(xi,yi,zi) = lbv(xi,yi,zi) - lzv(xi,yi,zi)

      bdotb1 = bdotb1 + lbv(xi,yi,zi)**2
      lpv(xi,yi,zi) = lbv(xi,yi,zi)

      ! WMBS - apply (D+U) part of M
      ! am4 through 6 store periodic dielectric on 0 indices
      ! ltv needs to wrap to element 1 for periodicity

      ltv(xi,yi,zi) = lpv(xi,yi,zi) + lam1(xi,yi,zi)*ltv(xi+1,yi,zi) &
                                    + lam2(xi,yi,zi)*ltv(xi,yi+1,zi) &
                                    + lam3(xi,yi,zi)*ltv(xi,yi,zi+1) &
                                    + lam4(xi,yi,zi)*ltv(xp,yi,zi  ) &
                                    + lam5(xi,yi,zi)*ltv(xi,yp,zi  ) &
                                    + lam6(xi,yi,zi)*ltv(xi,yi,zp  )
                      
   end do
   l_itn = 0
   unconv = .true.


end subroutine pb_piccg_initialize


end subroutine pb_piccg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! pb_sor
!
! SOR core routine for linearized FDPB equation, non-periodic version.  
!
! Authors:
! Jun Wang
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine pb_sor(phi,xs)

   implicit none

   ! Passed variables

   _REAL_ phi(l_xmymzm), xs(1-l_xmym:l_xmymzm+l_xmym)

   ! Local variables

   logical uconvg
   integer ii
   _REAL_ wsor1

   ! calculate initial norm

   l_inorm = sum( abs( l_bv(1:l_xmymzm) ) )
   wsor1 = ONE - l_wsor
   uconvg = .true.
   l_itn = 0

   l_zv(1:l_xmymzm) = ONE/l_ad(1:l_xmymzm)

   ! main loop

   do while ( uconvg )

      do ii = 1, l_xmymzm
         xs(ii) = wsor1*xs(ii) + l_wsor * (l_am1(ii-1     ) * xs(ii-1     ) + &
                                           l_am1(ii       ) * xs(ii+1     ) + &
                                           l_am2(ii-l_xm  ) * xs(ii-l_xm  ) + &
                                           l_am2(ii       ) * xs(ii+l_xm  ) + &
                                           l_am3(ii-l_xmym) * xs(ii-l_xmym) + &
                                           l_am3(ii       ) * xs(ii+l_xmym) + &
                                           l_bv(ii        ))* l_zv(ii     )
      end do

      l_itn = l_itn + 1

      do ii = 1, l_xmymzm
         phi(ii) = l_am1(ii-1     ) * xs(ii-1     ) + l_am1(ii)*xs(ii+1     ) + &
                   l_am2(ii-l_xm  ) * xs(ii-l_xm  ) + l_am2(ii)*xs(ii+l_xm  ) + &
                   l_am3(ii-l_xmym) * xs(ii-l_xmym) + l_am3(ii)*xs(ii+l_xmym) + &
                   l_bv(ii) - l_ad(ii)* xs(ii)
      end do
      l_norm = sum(abs(phi(1:l_xmymzm)))

      if ( l_itn >= l_maxitn .or. l_norm <= l_accept*l_inorm ) then
         uconvg = .false.
         if ( l_itn >= l_maxitn ) then
            write(6,'(a)') 'PBMD WARNING: SOR maxitn exceeded!'
         end if
      end if

   end do

   phi(1:l_xmymzm) = xs(1:l_xmymzm)


end subroutine pb_sor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! pb_psor
!
! SOR core routine for linearized FDPB equation, periodic version.  
!
! Authors:
! Wesley Smith
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine pb_psor(phi,xs)

   implicit none

   ! Passed variables

   _REAL_ phi(l_xmymzm), xs(1-l_xmym:l_xmymzm+l_xmym)

   ! Local variables

   logical uconvg
   integer ii,xi,yi,zi
   _REAL_ wsor1
   _REAL_ netcrg

   ! neutralize grid charge in periodic Poisson systems
   ! salt would automatically do it in the periodic PB systems.

   if ( l_pbkappa == ZERO ) then
      netcrg = -sum(l_bv(1:l_xmymzm))/l_xmymzm
      if ( abs(netcrg) > 1.0d-9 ) then
         write(6,'(a)') 'PB Warning: System net charge detected and removed for bcopt=10'
         l_bv(1:l_xmymzm) = l_bv(1:l_xmymzm) + netcrg
      end if
   end if

   ! calculate initial norm

   l_inorm = sum( ABS( l_bv(1:l_xmymzm) ) )
   wsor1 = ONE - l_wsor
   uconvg = .true.
   l_itn = 0

   l_zv(1:l_xmymzm) = ONE/l_ad(1:l_xmymzm)

   ! main loop

   do while ( uconvg )

      do ii = 1, l_xmymzm
         xi = pl_ind(ii,1); yi = pl_ind(ii,2); zi = pl_ind(ii,3)
         xs(pl_ind3d(xi,yi,zi)) = wsor1*xs(pl_ind3d(xi,yi,zi)) + &
                               l_wsor*( l_am1(pl_ind3d(xi-1,yi,zi)) * xs(pl_ind3d(xi-1,yi,zi)) + &
                                        l_am1(pl_ind3d(xi,yi,zi  )) * xs(pl_ind3d(xi+1,yi,zi)) + &
                                        l_am2(pl_ind3d(xi,yi-1,zi)) * xs(pl_ind3d(xi,yi-1,zi)) + &
                                        l_am2(pl_ind3d(xi,yi,zi  )) * xs(pl_ind3d(xi,yi+1,zi)) + &
                                        l_am3(pl_ind3d(xi,yi,zi-1)) * xs(pl_ind3d(xi,yi,zi-1)) + &
                                        l_am3(pl_ind3d(xi,yi,zi  )) * xs(pl_ind3d(xi,yi,zi+1)) + &
                                        l_bv(pl_ind3d(xi,yi,zi  )) ) * l_zv(pl_ind3d(xi,yi,zi))
      end do

      l_itn = l_itn + 1

      forall ( xi = 1:l_xm, yi = 1:l_ym, zi = 1:l_zm )
         phi(pl_ind3d(xi,yi,zi)) =  l_am1(pl_ind3d(xi-1,yi,zi)) * xs(pl_ind3d(xi-1,yi,zi)) + & 
                                    l_am1(pl_ind3d(xi,yi,zi  )) * xs(pl_ind3d(xi+1,yi,zi)) + &
                                    l_am2(pl_ind3d(xi,yi-1,zi)) * xs(pl_ind3d(xi,yi-1,zi)) + &
                                    l_am2(pl_ind3d(xi,yi,zi  )) * xs(pl_ind3d(xi,yi+1,zi)) + &
                                    l_am3(pl_ind3d(xi,yi,zi-1)) * xs(pl_ind3d(xi,yi,zi-1)) + &
                                    l_am3(pl_ind3d(xi,yi,zi  )) * xs(pl_ind3d(xi,yi,zi+1)) + &
                                    l_bv(pl_ind3d(xi,yi,zi  )) - l_ad(pl_ind3d(xi,yi,zi  ))* &
                                    xs(pl_ind3d(xi,yi,zi  ))
      end forall

      l_norm = sum(ABS(phi(1:l_xmymzm)))
      if ( l_itn >= l_maxitn .or. l_norm <= l_accept*l_inorm ) then
         uconvg = .false.
         if ( l_itn >= l_maxitn ) then
            write(6,'(a)') 'PBMD WARNING: SOR maxitn exceeded!'
         end if
      end if

   end do

   ! periodic systems have a non-trivial null space. Most notably, the null space will
   ! include the uniform distribution. In order to attain consistent results, the solution
   ! must be purged of its content by subratcting the projection of the soultion distribution
   ! onto the uniform distribution. This is particularly important for non-charge neutral
   ! systems to which a uniform neutralizing charge was added in order to solve the system.

   phi(1:l_xmymzm) = xs(1:l_xmymzm) - sum(xs(1:l_xmymzm)) / l_xmymzm


end subroutine pb_psor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! pb_pcg
!
! CG core routine for linearized FDPB equation, periodic version.  
!
! Authors:
! Jun Wang
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine pb_pcg ( phi, xs )

   ! Passed variables

   _REAL_ phi(l_xmymzm), xs(1-l_xmym:l_xmymzm+l_xmym)

   ! Local variables

   logical uconvg
   integer i, j, k, ii, jj
   _REAL_ alpha, beta, pdotz, bdotb1, bdotb2
   _REAL_ netcrg

   ! neutralize grid charge in periodic Poisson systems
   ! salt would automatically do it in the periodic PB systems.

   if ( l_pbkappa == ZERO ) then
      netcrg = -sum(l_bv(1:l_xmymzm))/l_xmymzm
      if ( ABS(netcrg) > 1.0d-9 ) then
         write(6,'(a)') 'PB Warning: System net charge detected and removed for bcopt=10'
         l_bv(1:l_xmymzm) = l_bv(1:l_xmymzm) + netcrg
      end if
   end if

   ! compute ||b||

   l_inorm = sum(ABS(l_bv(1:l_xmymzm)))

   !write(6,'(a,e10.3)')  'itn & norm ', 0, l_inorm
   
   ! compute b - A * x(0) and save it in r(0)
   ! p(0) = r(0)
   !
   ! iteration 0:
   !   compute <r(0),r(0)>

   l_itn = 0
   bdotb1 = ZERO

   ! running edge periodicity loops

   do j = 1, l_ym; do k = 1, l_zm
      ii = 1+(j-1)*l_xm+(k-1)*l_xmym
      jj = ii + l_xm - 1
      l_bv(ii)=l_bv(ii)+l_am4(ii)*xs(jj)
      l_bv(jj)=l_bv(jj)+l_am4(ii)*xs(ii)
   end do; end do
   do i = 1, l_xm; do k = 1, l_zm
      ii = i+(k-1)*l_xmym
      jj = ii + l_xmym - l_xm
      l_bv(ii)=l_bv(ii)+l_am5(ii)*xs(jj)
      l_bv(jj)=l_bv(jj)+l_am5(ii)*xs(ii)
   end do; end do
   do i = 1, l_xm; do j = 1, l_ym
      ii = i+(j-1)*l_xm
      jj = ii + l_xmymzm - l_xmym
      l_bv(ii)=l_bv(ii)+l_am6(ii)*xs(jj)
      l_bv(jj)=l_bv(jj)+l_am6(ii)*xs(ii)
   end do; end do

   ! running central update loop

   do i = 1, l_xmymzm
      l_bv(i) = l_bv(i) + l_am3(i-l_xmym)*xs(i-l_xmym) &
                        + l_am2(i-l_xm  )*xs(i-l_xm  ) &
                        + l_am1(i-1     )*xs(i-1     ) &
                        - l_ad (i       )*xs(i       ) &
                        + l_am1(i       )*xs(i+1     ) &
                        + l_am2(i       )*xs(i+l_xm  ) &
                        + l_am3(i       )*xs(i+l_xmym)
      bdotb1 = bdotb1 + l_bv(i)*l_bv(i)
      l_pv(i) = l_bv(i)
   end do

   ! main loop

   uconvg = .true.
   do while ( uconvg )

      ! iteration i:
      ! compute Ap(i) = A * p(i)
      ! compute alpha(i) = <r(i),r(i)>/<p(i),Ap(i)>

      pdotz = ZERO
      l_zv = ZERO
      do j = 1, l_ym; do k = 1, l_zm
         ii = 1+(j-1)*l_xm+(k-1)*l_xmym
         jj = ii + l_xm - 1
         l_zv(ii)=l_zv(ii)-l_am4(ii)*l_pv(jj)
         l_zv(jj)=l_zv(jj)-l_am4(ii)*l_pv(ii)
      end do; end do
      do i = 1, l_xm; do k = 1, l_zm
         ii = i+(k-1)*l_xmym
         jj = ii + l_xmym - l_xm
         l_zv(ii)=l_zv(ii)-l_am5(ii)*l_pv(jj)
         l_zv(jj)=l_zv(jj)-l_am5(ii)*l_pv(ii)
      end do; end do
      do i = 1, l_xm; do j = 1, l_ym
         ii = i+(j-1)*l_xm
         jj = ii + l_xmymzm - l_xmym
         l_zv(ii)=l_zv(ii)-l_am6(ii)*l_pv(jj)
         l_zv(jj)=l_zv(jj)-l_am6(ii)*l_pv(ii)
      end do; end do

      do i = 1, l_xmymzm
         l_zv(i) = l_ad (i       )*l_pv(i       ) &
                 - l_am3(i-l_xmym)*l_pv(i-l_xmym) &
                 - l_am2(i-l_xm  )*l_pv(i-l_xm  ) &
                 - l_am1(i-1     )*l_pv(i-1     ) &
                 - l_am1(i       )*l_pv(i+1     ) &
                 - l_am2(i       )*l_pv(i+l_xm  ) &
                 - l_am3(i       )*l_pv(i+l_xmym) &
                 + l_zv(i)

         pdotz = pdotz + l_pv(i)*l_zv(i)
      end do

      l_itn = l_itn + 1

      ! update x(i+1) = x(i) + alpha(i) p(i)
      !        r(i+1) = r(i) - alpha(i) Ap(i)

      alpha = bdotb1/pdotz
      l_norm = ZERO
      bdotb2 = ZERO
      do i = 1, l_xmymzm
         xs(i)   = xs(i)   + alpha*l_pv(i)
         l_bv(i) = l_bv(i) - alpha*l_zv(i)
         l_norm  = l_norm  +   ABS(l_bv(i))

      ! compute beta(i) = <r(i+1),r(i+1)>/<r(i),r(i)>, part one

         bdotb2  = bdotb2  + l_bv(i)*l_bv(i)
      end do

      !write(6,'(a,i5,3e10.3)')  'itn & norm ',l_itn, l_norm, l_inorm, l_accept

      ! check convergence

      if ( l_itn >= l_maxitn .or. l_norm <= l_accept*l_inorm ) then

         uconvg = .false.
         if ( l_itn >= l_maxitn ) then
            write(6,'(a)') 'PBMD WARNING: CG maxitn exceeded!'
         end if

      else

      ! compute beta(i) = <r(i+1),r(i+1)>/<r(i),r(i)>, part two

         beta   = bdotb2/bdotb1
         bdotb1 = bdotb2

      ! update p(i+1) = r(i+1) + beta(i) p(i)

         l_pv(1:l_xmymzm) = l_bv(1:l_xmymzm) + beta*l_pv(1:l_xmymzm)
      end if
   end do

   ! periodic systems have a non-trivial null space. Most notably, the null space will
   ! include the uniform distribution. In order to attain consistent results, the solution
   ! must be purged of its content by subratcting the projection of the soultion distribution
   ! onto the uniform distribution. This is particularly important for non-charge neutral
   ! systems to which a uniform neutralizing charge was added in order to solve the system.

   phi(1:l_xmymzm) = xs(1:l_xmymzm) - sum(xs(1:l_xmymzm)) / l_xmymzm


end subroutine pb_pcg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! pb_cg
!
! CG core routine for linearized FDPB equation, non-periodic version.  
!
! Authors:
! Ray Luo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine pb_cg(phi,xs)

   implicit none

   ! Passed variables

   _REAL_ phi(l_xmymzm), xs(1-l_xmym:l_xmymzm+l_xmym)

   ! Local variables

   logical uconvg
   integer i, j
   _REAL_ alpha, beta, pdotz, bdotb1, bdotb2

   l_inorm = sum(abs(l_bv(1:l_xmymzm)))

   !  compute b - A * x(0) and save it in r(0)
   !  p(0) = r(0)

   !  iteration 0:
   !  compute <r(0),r(0)>

   l_itn = 0
   bdotb1 = ZERO
   do i = 1, l_xmymzm
      l_bv(i)  = l_bv(i) + l_am3(i-l_xmym)*xs(i-l_xmym) &
                         + l_am2(i-l_xm  )*xs(i-l_xm  ) &
                         + l_am1(i-1     )*xs(i-1     ) &
                         - l_ad (i       )*xs(i       ) &
                         + l_am1(i       )*xs(i+1     ) &
                         + l_am2(i       )*xs(i+l_xm  ) &
                         + l_am3(i       )*xs(i+l_xmym)
      bdotb1 = bdotb1 + l_bv(i)*l_bv(i)
      l_pv(i) = l_bv(i)
   end do

   ! main loop

   uconvg = .true.
   do while ( uconvg )

      ! iteration i:

      ! compute Ap(i) = A * p(i)
      ! compute alpha(i) = <r(i),r(i)>/<p(i),Ap(i)>

      pdotz = ZERO
      do i = 1, l_xmymzm
         l_zv(i) = - l_am3(i-l_xmym)*l_pv(i-l_xmym) &
                   - l_am2(i-l_xm  )*l_pv(i-l_xm  ) &
                   - l_am1(i-1     )*l_pv(i-1     ) &
                   - l_am1(i       )*l_pv(i+1     ) &
                   - l_am2(i       )*l_pv(i+l_xm  ) &
                   - l_am3(i       )*l_pv(i+l_xmym) &
                   + l_ad (i       )*l_pv(i       )  
         pdotz = pdotz + l_pv(i)*l_zv(i)
      end do

      ! iteration i+1:

      l_itn = l_itn + 1

      ! update x(i+1) = x(i) + alpha(i) p(i)
      !        r(i+1) = r(i) - alpha(i) Ap(i)

      alpha = bdotb1/pdotz
      l_norm = ZERO
      bdotb2 = ZERO
      do i = 1,l_xmymzm
         xs(i) = xs(i) + alpha*l_pv(i)
         l_bv(i) = l_bv(i) - alpha*l_zv(i)
         l_norm  = l_norm  +  abs(l_bv(i))

      ! compute beta(i) = <r(i+1),r(i+1)>/<r(i),r(i)>, part one

         bdotb2 = bdotb2 + l_bv(i)*l_bv(i)
      end do

      !write(6,'(a,i5,3e10.3)')  'itn & norm ',l_itn, l_norm, l_inorm, l_accept

      ! check convergence

      if ( l_itn >= l_maxitn .or. l_norm <= l_accept*l_inorm ) then

         uconvg = .false.
         if ( l_itn >= l_maxitn ) then
            write(6,'(a)') 'PBMD WARNING: CG l_maxitn exceeded!'
         end if

      else

      ! compute beta(i) = <r(i+1),r(i+1)>/<r(i),r(i)>, part two

         beta   = bdotb2/bdotb1
         bdotb1 = bdotb2

      ! update p(i+1) = r(i+1) + beta(i) p(i)

         l_pv(1:l_xmymzm) = l_bv(1:l_xmymzm) + beta*l_pv(1:l_xmymzm)
      end if
   end do

   phi(1:l_xmymzm) = xs(1:l_xmymzm)


end subroutine pb_cg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! pb_mg
!
! MG core routine for linearized FDPB equation, non-periodic and periodic versions.  
!
! Authors:
! Jun Wang, Wesley Smith
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine pb_mg(phi,xs)

   implicit none

   ! Passed variables

   _REAL_ phi(l_xmymzm), xs(1-l_xmym:l_xmymzm+l_xmym)

   ! Local variables

   logical uconvg
   integer j,lxly,lxlylz,p1,p2,gid_v
   _REAL_  netcrg

   ! neutralize grid charge in periodic Poisson systems
   ! salt would automatically do it in the periodic PB systems.

   if ( l_bcopt == 10 .and. l_pbkappa == ZERO ) then
      netcrg = -sum(l_bv(1:l_xmymzm))/l_xmymzm
      if ( ABS(netcrg) > 1.0d-9 ) then
         write(6,'(a)') 'PB Warning: System net charge detected and removed for bcopt=10'
         l_bv(1:l_xmymzm) = l_bv(1:l_xmymzm) + netcrg
      end if
   end if 
 
   l_inorm = sum(abs(l_bv(1:l_xmymzm)))
   l_itn = 0
   uconvg = .true.

   ! main loop

   do while ( uconvg )

      mg_onorm = 9.9d99

      ! this goes up (initially 1) 

      do j = 1, mg_nlevel-1
         gid = j
         lxly = mg_size(1,j)*mg_size(2,j)
         lxlylz = lxly*mg_size(3,j)
         call relax(l_xv(mg_x_idx(j)),mg_size(1,j),mg_size(2,j),mg_size(3,j), &
                    l_am1(mg_index_ext(j)),l_am2(mg_index_ext(j)),l_am3(mg_index_ext(j)), &
                    l_ad(mg_index(j)),l_bv(mg_index(j)),l_rv(mg_index(j)), &
                    lxly,lxlylz,ncyc_before,l_accept,mg_onorm(j) )
         call restrict_bv(l_rv(mg_index(j)), mg_size(1,j), mg_size(2,j), mg_size(3,j), &
                          l_bv(mg_index(j+1)), mg_size(1,j+1),mg_size(2,j+1),mg_size(3,j+1) )
         l_xv(mg_x_idx(j+1):mg_x_idx(j+2)-1) = ZERO
      end do

      ! this is the top j = 4

         lxly = mg_size(1,j)*mg_size(2,j)
         lxlylz = lxly*mg_size(3,j)
         gid = j
         call relax(l_xv(mg_x_idx(j)),mg_size(1,j),mg_size(2,j),mg_size(3,j), &
                    l_am1(mg_index_ext(j)),l_am2(mg_index_ext(j)),l_am3(mg_index_ext(j)), &
                    l_ad(mg_index(j)),l_bv(mg_index(j)),l_rv(mg_index(j)), &
                    lxly,lxlylz,-1,l_accept,mg_onorm(j))

      ! this goes down (finally 1)

      do j = mg_nlevel-1, 1, -1    
         gid = j
         lxly = mg_size(1,j)*mg_size(2,j)
         lxlylz = lxly*mg_size(3,j)
         p1 = mg_x_idx(j+1)+mg_size(1,j+1)*mg_size(2,j+1)
         p2 = mg_x_idx(j)+lxly
         if ( l_bcopt == 10 ) &
         call pinterpolate(l_xv(p1),mg_size(1,j+1), mg_size(2,j+1), mg_size(3,j+1),&
                           l_xv(p2),mg_size(1,j),   mg_size(2,j),   mg_size(3,j) , &
                           l_am1(mg_index_ext(j)),l_am2(mg_index_ext(j)), &
                           l_am3(mg_index_ext(j)),l_bz(mg_index(j)),l_epsout)
         if ( l_bcopt /= 10 ) &
         call  interpolate(l_xv(p1), mg_size(1,j+1), mg_size(2,j+1), mg_size(3,j+1),&
                           l_xv(p2), mg_size(1,j),   mg_size(2,j),   mg_size(3,j) , &
                           l_am1(mg_index_ext(j)),l_am2(mg_index_ext(j)), &
                           l_am3(mg_index_ext(j)),l_bz(mg_index(j)),l_epsout)

         call relax(l_xv(mg_x_idx(j)),mg_size(1,j),mg_size(2,j),mg_size(3,j), &
                    l_am1(mg_index_ext(j)),l_am2(mg_index_ext(j)),l_am3(mg_index_ext(j)), &
                    l_ad(mg_index(j)),l_bv(mg_index(j)), l_rv(mg_index(j)), &
                    lxly,lxlylz,ncyc_after,l_accept,mg_onorm(j) )
      end do

      l_itn = l_itn + 1
      l_norm = sum(abs(l_rv(1:l_xmymzm)))
      !write(6,'(a)') 'pb_mg itn & norm', litn, lnorm
  
      ! check convergence
  
      if ( l_itn >= l_maxitn .or. l_norm <= l_inorm*l_accept ) then
         uconvg = .false.
         if ( l_itn >= l_maxitn ) then
            write(6,'(a)') 'PB_MG WARNING: maxitn exceeded!'
         endif
      end if

   end do

   xs (1:l_xmymzm) = l_xv(l_xmym+1:l_xmym+l_xmymzm)

   if ( l_bcopt /= 10 ) phi(1:l_xmymzm) = xs(1:l_xmymzm)
   if ( l_bcopt == 10 ) phi(1:l_xmymzm) = xs(1:l_xmymzm) - sum(xs(1:l_xmymzm))/l_xmymzm
  
end subroutine pb_mg

!===========================================================================

subroutine relax(xs,nx,ny,nz,lam1,lam2,lam3,lad,lbv,lrv,nxny,nxnynz,ncyc,accept,onorm)
 
   implicit none

   integer nx,ny,nz,nxny,nxnynz,ncyc
   _REAL_ xs(1-nxny:nxnynz+nxny), lam1(1-nxny:nxnynz), lam2(1-nxny:nxnynz), lam3(1-nxny:nxnynz)
   _REAL_ lad(1:nxnynz), lbv(1:nxnynz), lrv(1:nxnynz)
   _REAL_ accept, onorm

   logical luconvg
   integer ii,litn,itmax, xi,yi,zi,xsi,ysi,zsi
   integer itn_checknorm
   _REAL_ wsor, wsor1, linorm, lnorm, lxmymzm
   !_REAL_  netcrg

   if ( ncyc > 0 ) then
      itn_checknorm = ncyc
      wsor = 1.0d0
      wsor1 = ONE - wsor
   else
      itn_checknorm = 10
      wsor = 1.9d0
      wsor1 = ONE - wsor
   end if

   ! neutralize the charge grid for periodic systems

   !if ( l_bcopt == 10 ) then
   !   netcrg = -sum(lbv(1:nxnynz))/nxnynz
   !   if ( ABS(netcrg) > 1.0d-9 ) then
   !      write(6,'(a)') 'PB Warning: System net charge detected and removed for bcopt=10'
   !      lbv(1:nxnynz) = l_bv(1:nxnynz) + netcrg
   !   end if
   !end if   

   linorm = sum(abs(lbv(1:nxnynz)))
   l_zv(1:nxnynz) = ONE/lad(1:nxnynz)

   litn = 0
   itmax = 1000
   luconvg = .true.
   !write(6,'(a,2i5,2e10.3)')  ' relax gid itn & norm ', gid, litn, linorm, onorm

   ! start the sor iteration ...

   do while ( luconvg )

      ! periodic systems

      if ( l_bcopt == 10 ) then 
         if ( gid == 1 ) then
            xsi = 0; ysi = 0; zsi = 0
         else
            xsi = sum(2+mg_size(1,1:(gid-1))); ysi = sum(2+mg_size(2,1:(gid-1))); zsi = sum(2+mg_size(3,1:(gid-1)))
         end if

         do zi =1, nz; do yi = 1, ny; do xi = 1, nx ! apparently this cannot be done out of order (i.e. no forall)
            xs(pl_ind3d(xsi+xi,ysi+yi,zsi+zi)) = wsor1 * xs(pl_ind3d(xsi+xi,ysi+yi,zsi+zi)) + wsor * (      &  
                  lam1( pl_ind3d(xsi+xi-1,ysi+yi  ,zsi+zi  ) )*xs( pl_ind3d(xsi+xi-1,ysi+yi  ,zsi+zi  ) ) + &
                  lam1( pl_ind3d(xsi+xi  ,ysi+yi  ,zsi+zi  ) )*xs( pl_ind3d(xsi+xi+1,ysi+yi  ,zsi+zi  ) ) + &
                  lam2( pl_ind3d(xsi+xi  ,ysi+yi-1,zsi+zi  ) )*xs( pl_ind3d(xsi+xi  ,ysi+yi-1,zsi+zi  ) ) + &
                  lam2( pl_ind3d(xsi+xi  ,ysi+yi  ,zsi+zi  ) )*xs( pl_ind3d(xsi+xi  ,ysi+yi+1,zsi+zi  ) ) + &
                  lam3( pl_ind3d(xsi+xi  ,ysi+yi  ,zsi+zi-1) )*xs( pl_ind3d(xsi+xi  ,ysi+yi  ,zsi+zi-1) ) + &
                  lam3( pl_ind3d(xsi+xi  ,ysi+yi  ,zsi+zi  ) )*xs( pl_ind3d(xsi+xi  ,ysi+yi  ,zsi+zi+1) ) + &
                lbv(pl_ind3d(xsi+xi,ysi+yi,zsi+zi)) ) * l_zv(pl_ind3d(xsi+xi,ysi+yi,zsi+zi))
         end do; end do; end do

         litn = litn + 1

         if ( mod(litn,itn_checknorm) == 0 ) then
            forall ( zi = 1:nz, yi = 1:ny, xi = 1:nx ) ! this can be done out of order
               lrv(pl_ind3d(xsi+xi,ysi+yi,zsi+zi)) = lbv(pl_ind3d(xsi+xi,ysi+yi,zsi+zi))- &
                  lad ( pl_ind3d(xsi+xi  ,ysi+yi  ,zsi+zi  ) )*xs( pl_ind3d(xsi+xi  ,ysi+yi  ,zsi+zi  ) ) + &
                  lam1( pl_ind3d(xsi+xi-1,ysi+yi  ,zsi+zi  ) )*xs( pl_ind3d(xsi+xi-1,ysi+yi  ,zsi+zi  ) ) + &
                  lam1( pl_ind3d(xsi+xi  ,ysi+yi  ,zsi+zi  ) )*xs( pl_ind3d(xsi+xi+1,ysi+yi  ,zsi+zi  ) ) + &
                  lam2( pl_ind3d(xsi+xi  ,ysi+yi-1,zsi+zi  ) )*xs( pl_ind3d(xsi+xi  ,ysi+yi-1,zsi+zi  ) ) + &
                  lam2( pl_ind3d(xsi+xi  ,ysi+yi  ,zsi+zi  ) )*xs( pl_ind3d(xsi+xi  ,ysi+yi+1,zsi+zi  ) ) + &
                  lam3( pl_ind3d(xsi+xi  ,ysi+yi  ,zsi+zi-1) )*xs( pl_ind3d(xsi+xi  ,ysi+yi  ,zsi+zi-1) ) + &
                  lam3( pl_ind3d(xsi+xi  ,ysi+yi  ,zsi+zi  ) )*xs( pl_ind3d(xsi+xi  ,ysi+yi  ,zsi+zi+1) )  
            end forall
            
            lnorm = sum(abs(lrv(1:nxnynz)))
            !write(6,*)  ' relax itn & norm ', litn, lnorm, onorm

            ! check convergence
  
            if ( litn .ge. itmax .or. ( ncyc .gt. 0 .and. (litn .ge. ncyc .and. lnorm < onorm ) ) &
                 .or. lnorm .le. accept*linorm ) then
               luconvg = .false.
               if ( ncyc .gt. 0 .and. litn .ge. ncyc .and. lnorm > onorm ) then
                  write(6,'(a,3i6,2f12.4)') "PB_MG FAILED: ncyc, itn, gid, norm, onorm", ncyc, litn, gid, lnorm, onorm
                  stop
               end if

               if ( ncyc .gt. 0 ) onorm = lnorm 
               if ( litn .ge. itmax ) then
                  write(6,'(a)') 'PB_MG WARNING: SOR maxitn exceeded!'
               end if
            end if
         end if

      ! non-periodic systems

      else
 
         do ii = 1, nxnynz, 1
            xs(ii) = wsor1*xs(ii) + wsor * (lam1(ii-1   ) * xs(ii-1   ) + &
                                            lam1(ii     ) * xs(ii+1   ) + &
                                            lam2(ii-nx  ) * xs(ii-nx  ) + &
                                            lam2(ii     ) * xs(ii+nx  ) + &
                                            lam3(ii-nxny) * xs(ii-nxny) + &
                                            lam3(ii     ) * xs(ii+nxny) + &
                                             lbv(ii     )             ) * l_zv(ii)
         end do

         litn = litn + 1
         
         if ( mod(litn,itn_checknorm) == 0 ) then
            forall (ii = 1:nxnynz)
               lrv(ii) = lam1(ii-1) * xs(ii-1) + lam1(ii)* xs(ii+1) + &
                        lam2(ii-nx) * xs(ii-nx) + lam2(ii)*xs(ii+nx) + &
                        lam3(ii-nxny) * xs(ii-nxny) + lam3(ii)*xs(ii+nxny) + &
                        lbv(ii) - lad(ii)* xs(ii)
            end forall

            lnorm = sum(abs(lrv(1:nxnynz)))
  
            ! check convergence
 
            if ( litn .ge. itmax .or. ( ncyc .gt. 0 .and. (litn .ge. ncyc .and. lnorm < onorm ) ) &
                 .or. lnorm .le. accept*linorm ) then
               luconvg = .false.
               if ( ncyc .gt. 0 .and. litn .ge. ncyc .and. lnorm > onorm ) then
                  write(6,'(a,3i6,2f12.4)') "PB_MG FAILED: ncyc, itn, gid, norm, onorm", ncyc, litn, gid, lnorm, onorm
                  stop
               end if

               if ( ncyc .gt. 0 ) onorm = lnorm 
               if ( litn .ge. itmax ) then
                  write(6,'(a)') 'PB_MG WARNING: SOR maxitn exceeded!'
               end if
            end if
         end if
      end if

   end do


end subroutine relax

!===========================================================================
      
subroutine set_am_ad( epsx,epsy,epsz,iv,lam1,lam2,lam3,lad,lbz, &
                      xn,yn,zn,lfactor,epsout )

    implicit none

   integer xn,yn,zn
   _REAL_  lfactor, epsout
   _REAL_  epsx(xn,yn,zn), epsy(xn,yn,zn), epsz(xn,yn,zn), iv(xn,yn,zn)
   _REAL_  lam1(xn,yn,zn),lam2(xn,yn,zn),lam3(xn,yn,zn)
   _REAL_  lad(xn,yn,zn),lbz(xn,yn,zn)

   integer  i,j,k,i1,j1,k1
   _REAL_ lam1t,lam2t,lam3t

   lam1(1:xn,1:yn,1:zn) = epsx(1:xn,1:yn,1:zn)
   lam2(1:xn,1:yn,1:zn) = epsy(1:xn,1:yn,1:zn)
   lam3(1:xn,1:yn,1:zn) = epsz(1:xn,1:yn,1:zn)

   do k = 1, zn
      k1 = k-1
      do j = 1, yn
         j1 = j-1
         do i = 1, xn
            i1 = i-1
            lam1t = epsout
            if ( i1 /= 0 ) lam1t = lam1(i1,j,k)
            lam2t = epsout
            if ( j1 /= 0 ) lam2t = lam2(i,j1,k)
            lam3t = epsout
            if ( k1 /= 0 ) lam3t = lam3(i,j,k1)
            lad(i,j,k) = lam1t + lam1(i,j,k) + lam2t + lam2(i,j,k) &
                       + lam3t + lam3(i,j,k) 
         end do
      end do
   end do

   lbz = lfactor*iv
   lad = lad + lbz

   if ( l_bcopt /= 10 ) then
      do k = 1, zn; do j = 1, yn
         lam1(xn,j,k) = ZERO
      end do; end do
      do k = 1, zn; do i = 1, xn
         lam2(i,yn,k) = ZERO
      end do; end do
      do j = 1, yn; do i = 1, xn
         lam3(i,j,zn) = ZERO
      end do; end do
   end if

end subroutine set_am_ad

!===========================================================================

subroutine restrict_eps_map( epsx, epsy, epsz, xn, yn, zn, epsxr, epsyr, epszr, xnr, ynr, znr )

   implicit none
 
   integer xn, yn, zn, xnr, ynr, znr
   _REAL_ epsx(xn,yn,zn),epsxr(xnr,ynr,znr)
   _REAL_ epsy(xn,yn,zn),epsyr(xnr,ynr,znr)
   _REAL_ epsz(xn,yn,zn),epszr(xnr,ynr,znr)

   integer i,j,k,i2,j2,k2
   integer k2m,k2p,j2m,j2p,i2m,i2p

   if ( l_bcopt == 10 ) then

      do k = 1, znr
         k2 = 2*k-1; k2m = modulo(k2-2,zn)+1; k2p = modulo(k2,zn)+1
         do j = 1, ynr
            j2 = 2*j-1; j2m = modulo(j2-2,yn)+1; j2p = modulo(j2,yn)+1
            do i = 1, xnr
               i2 = 2*i-1; i2m = modulo(i2-2,xn)+1; i2p = modulo(i2,xn)+1

               epsxr(i,j,k) = hmav(epsx(i2  ,j2  ,k2  ),epsx(i2p ,j2  ,k2  ))/4.d0 + &
                             (hmav(epsx(i2  ,j2m ,k2  ),epsx(i2p ,j2m ,k2  )) + &
                              hmav(epsx(i2  ,j2p ,k2  ),epsx(i2p ,j2p ,k2  )) + &
                              hmav(epsx(i2  ,j2  ,k2m ),epsx(i2p ,j2  ,k2m )) + &
                              hmav(epsx(i2  ,j2  ,k2p ),epsx(i2p ,j2  ,k2p )))/8.d0 + &
                             (hmav(epsx(i2  ,j2m ,k2m ),epsx(i2p ,j2m ,k2m )) + &
                              hmav(epsx(i2  ,j2p ,k2m ),epsx(i2p ,j2p ,k2m )) + &
                              hmav(epsx(i2  ,j2m ,k2p ),epsx(i2p ,j2m ,k2p )) + &
                              hmav(epsx(i2  ,j2p ,k2p ),epsx(i2p ,j2p ,k2p )))/16.d0 

               epsyr(i,j,k) = hmav(epsy(i2  ,j2  ,k2  ),epsy(i2  ,j2p ,k2  ))/4.d0 + &
                             (hmav(epsy(i2m ,j2  ,k2  ),epsy(i2m ,j2p ,k2  )) + &
                              hmav(epsy(i2p ,j2  ,k2  ),epsy(i2p ,j2p ,k2  )) + &
                              hmav(epsy(i2  ,j2  ,k2m ),epsy(i2  ,j2p ,k2m )) + &
                              hmav(epsy(i2  ,j2  ,k2p ),epsy(i2  ,j2p ,k2p )))/8.d0 + &
                             (hmav(epsy(i2m ,j2  ,k2m ),epsy(i2m ,j2p ,k2m )) + &
                              hmav(epsy(i2p ,j2  ,k2m ),epsy(i2p ,j2p ,k2m )) + &
                              hmav(epsy(i2m ,j2  ,k2p ),epsy(i2m ,j2p ,k2p )) + &
                              hmav(epsy(i2p ,j2  ,k2p ),epsy(i2p ,j2p ,k2p )))/16.d0 

               epszr(i,j,k) = hmav(epsz(i2  ,j2  ,k2  ),epsz(i2  ,j2  ,k2p ))/4.d0 + &
                             (hmav(epsz(i2  ,j2m ,k2  ),epsz(i2  ,j2m ,k2p )) + &
                              hmav(epsz(i2  ,j2p ,k2  ),epsz(i2  ,j2p ,k2p )) + &
                              hmav(epsz(i2m ,j2  ,k2  ),epsz(i2m ,j2  ,k2p )) + &
                              hmav(epsz(i2p ,j2  ,k2  ),epsz(i2p ,j2  ,k2p )))/8.d0 + &
                             (hmav(epsz(i2m ,j2m ,k2  ),epsz(i2m ,j2m ,k2p )) + &
                              hmav(epsz(i2m ,j2p ,k2  ),epsz(i2m ,j2p ,k2p )) + &
                              hmav(epsz(i2p ,j2m ,k2  ),epsz(i2p ,j2m ,k2p )) + &
                              hmav(epsz(i2p ,j2p ,k2  ),epsz(i2p ,j2p ,k2p )))/16.d0 
            end do
         end do
      end do         

   else

      do k = 1, znr
         k2 = 2*k
         do j = 1, ynr
            j2 = 2*j
            do i = 1, xnr
               i2 = 2*i
               epsxr(i,j,k) = hmav(epsx(i2  ,j2  ,k2  ),epsx(i2+1,j2  ,k2  ))/4.d0 + &
                             (hmav(epsx(i2  ,j2-1,k2  ),epsx(i2+1,j2-1,k2  )) + &
                              hmav(epsx(i2  ,j2+1,k2  ),epsx(i2+1,j2+1,k2  )) + &
                              hmav(epsx(i2  ,j2  ,k2-1),epsx(i2+1,j2  ,k2-1)) + &
                              hmav(epsx(i2  ,j2  ,k2+1),epsx(i2+1,j2  ,k2+1)))/8.d0 + &
                             (hmav(epsx(i2  ,j2-1,k2-1),epsx(i2+1,j2-1,k2-1)) + &
                              hmav(epsx(i2  ,j2+1,k2-1),epsx(i2+1,j2+1,k2-1)) + &
                              hmav(epsx(i2  ,j2-1,k2+1),epsx(i2+1,j2-1,k2+1)) + &
                              hmav(epsx(i2  ,j2+1,k2+1),epsx(i2+1,j2+1,k2+1)))/16.d0 

               epsyr(i,j,k) = hmav(epsy(i2  ,j2  ,k2  ),epsy(i2  ,j2+1,k2  ))/4.d0 + &
                             (hmav(epsy(i2-1,j2  ,k2  ),epsy(i2-1,j2+1,k2  )) + & 
                              hmav(epsy(i2+1,j2  ,k2  ),epsy(i2+1,j2+1,k2  )) + &
                              hmav(epsy(i2  ,j2  ,k2-1),epsy(i2  ,j2+1,k2-1)) + &
                              hmav(epsy(i2  ,j2  ,k2+1),epsy(i2  ,j2+1,k2+1)))/8.d0 + &
                             (hmav(epsy(i2-1,j2  ,k2-1),epsy(i2-1,j2+1,k2-1)) + &
                              hmav(epsy(i2+1,j2  ,k2-1),epsy(i2+1,j2+1,k2-1)) + &
                              hmav(epsy(i2-1,j2  ,k2+1),epsy(i2-1,j2+1,k2+1)) + &
                              hmav(epsy(i2+1,j2  ,k2+1),epsy(i2+1,j2+1,k2+1)))/16.d0 

               epszr(i,j,k) = hmav(epsz(i2  ,j2  ,k2  ),epsz(i2  ,j2  ,k2+1))/4.d0 + &
                             (hmav(epsz(i2  ,j2-1,k2  ),epsz(i2  ,j2-1,k2+1)) + &
                              hmav(epsz(i2  ,j2+1,k2  ),epsz(i2  ,j2+1,k2+1)) + &
                              hmav(epsz(i2-1,j2  ,k2  ),epsz(i2-1,j2  ,k2+1)) + &
                              hmav(epsz(i2+1,j2  ,k2  ),epsz(i2+1,j2  ,k2+1)))/8.d0 + &
                             (hmav(epsz(i2-1,j2-1,k2  ),epsz(i2-1,j2-1,k2+1)) + &
                              hmav(epsz(i2-1,j2+1,k2  ),epsz(i2-1,j2+1,k2+1)) + &
                              hmav(epsz(i2+1,j2-1,k2  ),epsz(i2+1,j2-1,k2+1)) + &
                              hmav(epsz(i2+1,j2+1,k2  ),epsz(i2+1,j2+1,k2+1)))/16.d0 
            end do
         end do
      end do
   end if

contains 

function hmav(a,b)

   _REAL_ hmav, a, b

   hmav = 2.d0*a*b/(a+b)

end function hmav

function hmav4(a,b,c,d)

   _REAL_ hmav4, a,b,c,d

   hmav4 = 4.d0/(1.d0/a+1.d0/b+1.d0/c+1.d0/d)

end function hmav4

end subroutine restrict_eps_map

!===========================================================================

subroutine restrict_iv(ivf, xn, yn, zn, ivr, xnr, ynr, znr)

   implicit none

   integer xn, yn, zn, xnr, ynr, znr
   _REAL_ ivf(xn,yn,zn),ivr(xnr,ynr,znr)

   integer i,j,k,i2,j2,k2
   integer i2m,i2p,j2m,j2p,k2m,k2p

   if ( l_bcopt == 10 ) then
      do k = 1, znr
         k2 = 2*k-1; k2m = modulo(k2-2,zn)+1; k2p = modulo(k2,zn)+1
         do j = 1, ynr
            j2 = 2*j-1; j2m = modulo(j2-2,yn)+1; j2p = modulo(j2,yn)+1
            do i = 1, xnr
               i2 = 2*i-1; i2m = modulo(i2-2,xn)+1; i2p = modulo(i2,xn)+1

               ivr(i,j,k) =       ( ivf(i2m ,j2m ,k2m ) + TWO * ivf(i2  ,j2m ,k2m ) + ivf(i2p ,j2m ,k2m ) ) +  &
                           TWO  * ( ivf(i2m ,j2  ,k2m ) + TWO * ivf(i2  ,j2  ,k2m ) + ivf(i2p ,j2  ,k2m ) ) +  &
                                  ( ivf(i2m ,j2p ,k2m ) + TWO * ivf(i2  ,j2p ,k2m ) + ivf(i2p ,j2p ,k2m ) ) +  &
                           TWO  * ( ivf(i2m ,j2m ,k2  ) + TWO * ivf(i2  ,j2m ,k2  ) + ivf(i2p ,j2m ,k2  ) ) +  &
                           FOUR * ( ivf(i2m ,j2  ,k2  ) + TWO * ivf(i2  ,j2  ,k2  ) + ivf(i2p ,j2  ,k2  ) ) +  &
                           TWO  * ( ivf(i2m ,j2p ,k2  ) + TWO * ivf(i2  ,j2p ,k2  ) + ivf(i2p ,j2p ,k2  ) ) +  &
                                  ( ivf(i2m ,j2m ,k2p ) + TWO * ivf(i2  ,j2m ,k2p ) + ivf(i2p ,j2m ,k2p ) ) +  &
                           TWO  * ( ivf(i2m ,j2  ,k2p ) + TWO * ivf(i2  ,j2  ,k2p ) + ivf(i2p ,j2  ,k2p ) ) +  &
                                  ( ivf(i2m ,j2p ,k2p ) + TWO * ivf(i2  ,j2p ,k2p ) + ivf(i2p ,j2p ,k2p ) )
               ivr(i,j,k) = ivr(i,j,k) / 64.d0
            end do
         end do
      end do

   else

      do k = 1, znr
         k2 = k*2
         do j = 1, ynr
            j2 = j*2
            do i = 1, xnr
               i2 = i*2
               ivr(i,j,k) =       ( ivf(i2-1,j2-1,k2-1) + TWO * ivf(i2  ,j2-1,k2-1) + ivf(i2+1,j2-1,k2-1) ) +  &
                           TWO  * ( ivf(i2-1,j2  ,k2-1) + TWO * ivf(i2  ,j2  ,k2-1) + ivf(i2+1,j2  ,k2-1) ) +  &
                                  ( ivf(i2-1,j2+1,k2-1) + TWO * ivf(i2  ,j2+1,k2-1) + ivf(i2+1,j2+1,k2-1) ) +  &
                           TWO  * ( ivf(i2-1,j2-1,k2  ) + TWO * ivf(i2  ,j2-1,k2  ) + ivf(i2+1,j2-1,k2  ) ) +  &
                           FOUR * ( ivf(i2-1,j2  ,k2  ) + TWO * ivf(i2  ,j2  ,k2  ) + ivf(i2+1,j2  ,k2  ) ) +  &
                           TWO  * ( ivf(i2-1,j2+1,k2  ) + TWO * ivf(i2  ,j2+1,k2  ) + ivf(i2+1,j2+1,k2  ) ) +  &
                                  ( ivf(i2-1,j2-1,k2+1) + TWO * ivf(i2  ,j2-1,k2+1) + ivf(i2+1,j2-1,k2+1) ) +  &
                           TWO  * ( ivf(i2-1,j2  ,k2+1) + TWO * ivf(i2  ,j2  ,k2+1) + ivf(i2+1,j2  ,k2+1) ) +  &
                                  ( ivf(i2-1,j2+1,k2+1) + TWO * ivf(i2  ,j2+1,k2+1) + ivf(i2+1,j2+1,k2+1) )
               ivr(i,j,k) = ivr(i,j,k) / 64.d0
            end do
         end do
      end do

   end if

end subroutine restrict_iv

!===========================================================================

subroutine restrict_bv(bvf, xn, yn, zn, bvr, xnr, ynr, znr)
 
   implicit none 

   integer xn, yn, zn, xnr, ynr, znr
   _REAL_ bvf(xn,yn,zn),bvr(xnr,ynr,znr)

   integer i,j,k,i2,j2,k2
   integer i2m,i2p,j2m,j2p,k2m,k2p

   if ( l_bcopt == 10 ) then

      do k = 1, znr
         k2 = 2*k-1; k2m = modulo(k2-2,zn)+1; k2p = modulo(k2,zn)+1
         do j = 1, ynr
            j2 = 2*j-1; j2m = modulo(j2-2,yn)+1; j2p = modulo(j2,yn)+1
            do i = 1, xnr
               i2 = 2*i-1; i2m = modulo(i2-2,xn)+1; i2p = modulo(i2,xn)+1

               bvr(i,j,k) =       ( bvf(i2m ,j2m ,k2m ) + TWO * bvf(i2  ,j2m ,k2m ) + bvf(i2p ,j2m ,k2m ) ) +  &
                           TWO  * ( bvf(i2m ,j2  ,k2m ) + TWO * bvf(i2  ,j2  ,k2m ) + bvf(i2p ,j2  ,k2m ) ) +  &
                                  ( bvf(i2m ,j2p ,k2m ) + TWO * bvf(i2  ,j2p ,k2m ) + bvf(i2p ,j2p ,k2m ) ) +  &
                           TWO  * ( bvf(i2m ,j2m ,k2  ) + TWO * bvf(i2  ,j2m ,k2  ) + bvf(i2p ,j2m ,k2  ) ) +  &
                           FOUR * ( bvf(i2m ,j2  ,k2  ) + TWO * bvf(i2  ,j2  ,k2  ) + bvf(i2p ,j2  ,k2  ) ) +  &
                           TWO  * ( bvf(i2m ,j2p ,k2  ) + TWO * bvf(i2  ,j2p ,k2  ) + bvf(i2p ,j2p ,k2  ) ) +  &
                                  ( bvf(i2m ,j2m ,k2p ) + TWO * bvf(i2  ,j2m ,k2p ) + bvf(i2p ,j2m ,k2p ) ) +  &
                           TWO  * ( bvf(i2m ,j2  ,k2p ) + TWO * bvf(i2  ,j2  ,k2p ) + bvf(i2p ,j2  ,k2p ) ) +  &
                                  ( bvf(i2m ,j2p ,k2p ) + TWO * bvf(i2  ,j2p ,k2p ) + bvf(i2p ,j2p ,k2p ) )
               bvr(i,j,k) = bvr(i,j,k) / 16.d0                 
            end do
         end do
      end do

   else

      do k = 1, znr
         k2 = k*2
         do j = 1, ynr
            j2 = j*2
            do i = 1, xnr
               i2 = i*2
               bvr(i,j,k) =       ( bvf(i2-1,j2-1,k2-1) + TWO * bvf(i2  ,j2-1,k2-1) + bvf(i2+1,j2-1,k2-1) ) +  &
                           TWO  * ( bvf(i2-1,j2  ,k2-1) + TWO * bvf(i2  ,j2  ,k2-1) + bvf(i2+1,j2  ,k2-1) ) +  &
                                  ( bvf(i2-1,j2+1,k2-1) + TWO * bvf(i2  ,j2+1,k2-1) + bvf(i2+1,j2+1,k2-1) ) +  &
                           TWO  * ( bvf(i2-1,j2-1,k2  ) + TWO * bvf(i2  ,j2-1,k2  ) + bvf(i2+1,j2-1,k2  ) ) +  &
                           FOUR * ( bvf(i2-1,j2  ,k2  ) + TWO * bvf(i2  ,j2  ,k2  ) + bvf(i2+1,j2  ,k2  ) ) +  &
                           TWO  * ( bvf(i2-1,j2+1,k2  ) + TWO * bvf(i2  ,j2+1,k2  ) + bvf(i2+1,j2+1,k2  ) ) +  &
                                  ( bvf(i2-1,j2-1,k2+1) + TWO * bvf(i2  ,j2-1,k2+1) + bvf(i2+1,j2-1,k2+1) ) +  &
                           TWO  * ( bvf(i2-1,j2  ,k2+1) + TWO * bvf(i2  ,j2  ,k2+1) + bvf(i2+1,j2  ,k2+1) ) +  &
                                  ( bvf(i2-1,j2+1,k2+1) + TWO * bvf(i2  ,j2+1,k2+1) + bvf(i2+1,j2+1,k2+1) )
               bvr(i,j,k) = bvr(i,j,k) / 16.d0
            end do
         end do
      end do

   end if


end subroutine restrict_bv

!===========================================================================
! In the interest of not reinventing the wheel, the existing prolongation algorithm
! is used. This subroutine is simply a wrapper to jury-rig it using a "virtual grid".
! The idea here is to use scratch grids for which the upper faces have been padded.
! The lower edges of the coarse and dielectric scratch grid are assigned to the upper padding
! faces. These scratch grids can then be interpolated just as if this where dirichlet boundaries
! After that, the inner portion of the fine scratch grid is transfered back to the working fine
! grid. I.e. effictively scratch(2:xnf+1,2:ynf+1,2:znf+1) -> fineV(1:xnf,1:ynf,1:znf)...
subroutine pinterpolate( lxvc, xnc, ync, znc, lxvf, xnf, ynf, znf, lam1, lam2, lam3,lbzf,epsout)

   implicit none

   integer xnc,ync,znc,xnf,ynf,znf
   _REAL_ lxvc(1:xnc*ync*znc),lxvf(1:xnf*ynf*znf), epsout
   _REAL_ lam1(1-xnf*ynf:xnf*ynf*znf),lam2(1-xnf*ynf:xnf*ynf*znf),lam3(1-xnf*ynf:xnf*ynf*znf)
   _REAL_ lbzf(1:xnf*ynf*znf)

   integer i,j,k,xniyni,xniynizni,ii,ii2
   integer sic,sjc,skc,sif,sjf,skf,sii

   l_scratch_am1 = ZERO
   l_scratch_am2 = ZERO
   l_scratch_am3 = ZERO
   l_scratch_vc  = ZERO
   l_scratch_vf  = ZERO

   ! transfer coarse grid onto scratch, wrapping lower faces to upper faces

   do skc = 1, znc+1
      k = modulo(skc-1,znc)+1
      do sjc=1,ync+1
         j = modulo(sjc-1,ync)+1
         do sic=1,xnc+1
            i = modulo(sic-1,xnc)+1
            ii=i+(j-1)*xnc+(k-1)*xnc*ync
            sii=sic+(sjc-1)*(xnc+1)+(skc-1)*(xnc+1)*(ync+1)
            l_scratch_vc(sii)=lxvc(ii)
   end do; end do; end do

   ! Transfer dielectric tensor arrays to scratch grids with periodic wrapping

   do skf = 1, znf+3
      k = modulo(skf-2,znf)+1
      do sjf = 1, ynf+3
         j = modulo(sjf-2,ynf)+1
         do sif = 1, xnf+3
            i = modulo(sif-2,xnf)+1
            ii = i+(j-1)*xnf+(k-1)*xnf*ynf
            sii = sif+(sjf-1)*(xnf+3)+(skf-1)*(xnf+3)*(ynf+3)
            l_scratch_am1(sii) = lam1(ii)
            l_scratch_am2(sii) = lam2(ii)
            l_scratch_am3(sii) = lam3(ii)
            l_scratch_bz(sii)  = lbzf(ii)
   end do; end do; end do   

   ! Call interpolate(), feeding in scratch grids and appropriate dimensions

   call interpolate( l_scratch_vc(1:(xnc+1)*(ync+1)*(znc+1)), xnc+1, ync+1, znc+1,&
                     l_scratch_vf(1:(xnf+3)*(ynf+3)*(znf+3)), xnf+3, ynf+3, znf+3,&
                     l_scratch_am1(1-(xnf+3)*(ynf+3):(xnf+3)*(ynf+3)*(znf+3)),&
                     l_scratch_am2(1-(xnf+3)*(ynf+3):(xnf+3)*(ynf+3)*(znf+3)),& 
                     l_scratch_am3(1-(xnf+3)*(ynf+3):(xnf+3)*(ynf+3)*(znf+3)),&
                     l_scratch_bz(1:(xnf+3)*(ynf+3)*(znf+3)), epsout )

   ! Transfer scratch fine grid back to working fine grid.
   ! The inner 2:xn+1 by 2:yn+1 by 2:zn+1 is transfered to 1:xn by 1:yn by 1:zn of the
   ! working fine grid

   do skf = 2, znf+1
      k = skf-1;
      do sjf = 2, ynf+1
         j = sjf-1
         do sif = 2, xnf+1
            i = sif-1 
            ii = i+(j-1)*xnf+(k-1)*xnf*ynf
            sii = sif+(sjf-1)*(xnf+3)+(skf-1)*(xnf+3)*(ynf+3)
            lxvf(ii) = lxvf(ii)+l_scratch_vf(sii)
   end do; end do; end do

end subroutine

subroutine interpolate(v, xn, yn ,zn, vi, xni, yni, zni, lam1, lam2, lam3, lbz, epsout )

   integer xn,yn,zn,xni,yni,zni
   _REAL_ v(1:xn*yn*zn),vi(1:xni*yni*zni), epsout
   _REAL_ lam1(1-xni*yni:xni*yni*zni),lam2(1-xni*yni:xni*yni*zni),lam3(1-xni*yni:xni*yni*zni)
   _REAL_ lbz(1:xni*yni*zni)

   integer i,j,k,xniyni,xniynizni,ii,ii2

   if ( xn*2+1 /= xni .or. yn*2+1 /= yni .or. zn*2+1 /= zni ) then
      write (6,*) "Interpolation failed because of incorrect dimension"
      write (6,*) xn, yn, zn
      write (6,*) xni, yni, zni
      stop
   end if

   xniyni = xni*yni
   xniynizni = xni*yni*zni

   lam1(1-xniyni:0) = epsout
   lam2(1-xniyni:0) = epsout
   lam3(1-xniyni:0) = epsout
   do k = 1, zni; do j = 1, yni
      lam1(j*xni+(k-1)*xniyni) = epsout
   end do; end do
   do k = 1, zni; do i = 1, xni
      lam2(i-xni+k*xniyni) = epsout
   end do; end do
   do j = 1, yni; do i = 1, xni
      lam3(i+(j-1)*xni+xniynizni-xniyni) = epsout
   end do; end do

   ii = 0
   do k = 1, zn; do j=1 ,yn; do i=1, xn !iterate over each coarse grid
      ii = ii + 1 !coarse grid index
      ii2 = (k*2-1)*xniyni + (j*2-1)*xni + i*2 !fine grid index
      vi(ii2) = vi(ii2) + v(ii) !injection onto coincident fine grid
      !shift_1 defines direction to each of six adjacent type 2 points
      !shift_2 defines direction to the type 3 points adjacent to the type 2 points
      !shift_3 defines direction to the type 4 points adjacent to the type 3 points
                    !vi xnyn   xnynzn    l   v     lbz am_1 shift1   am2 shift2 am3 shift3
      call ipl_chain(vi,xniyni,xniynizni,ii2,v(ii),lbz,lam1,     -1,lam2,xni,  lam3,xniyni,xn,yn,zn)
      call ipl_chain(vi,xniyni,xniynizni,ii2,v(ii),lbz,lam1,     +1,lam2,xni,  lam3,xniyni,xn,yn,zn)
      call ipl_chain(vi,xniyni,xniynizni,ii2,v(ii),lbz,lam2,   -xni,lam1,  1,  lam3,xniyni,xn,yn,zn)
      call ipl_chain(vi,xniyni,xniynizni,ii2,v(ii),lbz,lam2,   +xni,lam1,  1,  lam3,xniyni,xn,yn,zn)
      call ipl_chain(vi,xniyni,xniynizni,ii2,v(ii),lbz,lam3,-xniyni,lam2,xni,  lam1,     1,xn,yn,zn)
      call ipl_chain(vi,xniyni,xniynizni,ii2,v(ii),lbz,lam3,+xniyni,lam2,xni,  lam1,     1,xn,yn,zn)
   end do; end do; end do

   lam1(1-xniyni:0) = ZERO
   lam2(1-xniyni:0) = ZERO
   lam3(1-xniyni:0) = ZERO
   do k = 1, zni; do j = 1, yni
      lam1(j*xni+(k-1)*xniyni) = ZERO
   end do; end do
   do k = 1, zni; do i = 1, xni
      lam2(i-xni+k*xniyni) = ZERO
   end do; end do
   do j = 1, yni; do i = 1, xni
      lam3(i+(j-1)*xni+xniynizni-xniyni) = ZERO
   end do; end do

contains

subroutine ipl_chain(vi,xnyn,xnynzn,l,v,lbz,am_1,shift_1,am_2,shift_2,am_3,shift_3,xn,yn,zn)

   implicit none

   integer xnyn,xnynzn,l,shift_1,shift_2,shift_3,xn,yn,zn
   _REAL_ vi(1:xnynzn),v,am_1(1-xnyn:xnynzn),am_2(1-xnyn:xnynzn),am_3(1-xnyn:xnynzn)
   _REAL_ lbz(1:xnynzn)

   integer l1
   _REAL_ v1
   
   ! v: value of coarse grid, vi: entire fine grid array

   v1 = ipl_comp1(v,l,lbz,am_1,xnyn,xnynzn,shift_1) !calculates contribution to type 2 point
   l1 = l + shift_1 !location of fine grid defined by coarsegrid at l and shift_1
   vi(l1) = vi(l1) + v1 !add contribution of coarse grid at l to fine grid at l1
      
   call ipl_chain2(vi,xnyn,xnynzn,l1,v1,lbz,am_1,shift_1,am_2,-shift_2,am_3,shift_3,xn,yn,zn)
   call ipl_chain2(vi,xnyn,xnynzn,l1,v1,lbz,am_1,shift_1,am_2,shift_2,am_3,shift_3,xn,yn,zn)
   call ipl_chain2(vi,xnyn,xnynzn,l1,v1,lbz,am_1,shift_1,am_3,-shift_3,am_2,shift_2,xn,yn,zn)
   call ipl_chain2(vi,xnyn,xnynzn,l1,v1,lbz,am_1,shift_1,am_3,shift_3,am_2,shift_2,xn,yn,zn)

end subroutine ipl_chain

subroutine ipl_chain2(vi,xnyn,xnynzn,l1,v1,lbz,am_1,shift_1,am_2,shift_2,am_3,shift_3,xn,yn,zn)

   implicit none

   integer xnyn,xnynzn,l1,shift_1,shift_2,shift_3,xn,yn,zn
   _REAL_ vi(1:xnynzn),v1,am_1(1-xnyn:xnynzn),am_2(1-xnyn:xnynzn),am_3(1-xnyn:xnynzn)
   _REAL_ lbz(1:xnynzn)

   integer l2
   _REAL_ v2

   v2 = ipl_comp2(v1,l1,lbz,am_1,am_2,xnyn,xnynzn,shift_1,shift_2)
   l2 = l1 + shift_2
   vi(l2) = vi(l2) + v2
   vi(l2-shift_3)= vi(l2-shift_3) + ipl_comp3(v2,l2,lbz,am_1,am_2,am_3,xnyn,xnynzn,shift_1,shift_2,-shift_3) 
   vi(l2+shift_3)= vi(l2+shift_3) + ipl_comp3(v2,l2,lbz,am_1,am_2,am_3,xnyn,xnynzn,shift_1,shift_2,+shift_3) 

end subroutine ipl_chain2

!compute contribution of coarse grid at l to a type 2 fine grid at l+shift_1
function ipl_comp1(v,l,lbz,am_1,xnyn,xnynzn,shift_1) 

   implicit none

   integer l,xnyn,xnynzn,shift_1
   _REAL_ ipl_comp1,v,am_1(1-xnyn:xnynzn),lbz(1:xnynzn)
   if ( shift_1 < 0 ) then 
      ipl_comp1 = v * am_1(l+shift_1) / ( lbz(l+shift_1) + am_1(l+2*shift_1) + am_1(l+shift_1) )
   else
      ipl_comp1 = v * am_1(l) / ( lbz(l+shift_1) + am_1(l) + am_1(l+shift_1) )
   end if
   return

end function ipl_comp1

function ipl_comp2(v,l,lbz,am_1,am_2,xnyn,xnynzn,shift_1,shift_2) 

   implicit none

   integer l,xnyn,xnynzn,shift_1,shift_2
   _REAL_ ipl_comp2,v,am_1(1-xnyn:xnynzn),am_2(1-xnyn:xnynzn)
   _REAL_ lbz(1:xnynzn)

   _REAL_ lad

   lad = am_1(l+shift_2) + am_1(l+shift_2-abs(shift_1)) + lbz(l+shift_2)
   if ( shift_2 < 0 ) then
      ipl_comp2 = v * am_2(l+shift_2) / ( am_2(l+2*shift_2) + am_2(l+shift_2) + lad)
   else
      ipl_comp2 = v * am_2(l) / ( am_2(l) + am_2(l+shift_2) + lad )
   end if
   return

end function ipl_comp2
            
function ipl_comp3(v,l,lbz,am_1,am_2,am_3,xnyn,xnynzn,shift_1,shift_2,shift_3) 

   implicit none

   integer l,xnyn,xnynzn,shift_1,shift_2,shift_3
   _REAL_ ipl_comp3,v,am_1(1-xnyn:xnynzn),am_2(1-xnyn:xnynzn),am_3(1-xnyn:xnynzn)
   _REAL_ lbz(1:xnynzn)

   _REAL_ lad

   lad = am_1(l+shift_3) + am_1(l+shift_3-abs(shift_1)) + am_2(l+shift_3) + &
         am_2(l+shift_3-abs(shift_2)) + lbz(l+shift_3)
   if ( shift_3 < 0 ) then
      ipl_comp3 = v * am_3(l+shift_3) / ( am_3(l+2*shift_3) + am_3(l+shift_3) + lad)
   else
      ipl_comp3 = v * am_3(l) / ( am_3(l) + am_3(l+shift_3) + lad )
   end if
   return

end function ipl_comp3

end subroutine interpolate

end module pb_lsolver
