! <compile=optimized>

!The 3D-RISM-KH software found here is copyright (c) 2011-2012 by 
!Andriy Kovalenko, Tyler Luchko and David A. Case.
!
!This program is free software: you can redistribute it and/or modify it
!under the terms of the GNU General Public License as published by the Free
!Software Foundation, either version 3 of the License, or (at your option)
!any later version.
!
!This program is distributed in the hope that it will be useful, but
!WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
!or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!for more details.
!
!You should have received a copy of the GNU General Public License in the
!../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
!
!Users of the 3D-RISM capability found here are requested to acknowledge
!use of the software in reports and publications.  Such acknowledgement
!should include the following citations:
!
!1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999); 
!ibid. 112:10391-10417 (2000).   
!
!2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by  
!F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.  
!
!3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
!and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010). 

#include "../include/dprec.fh"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Hypernetted chain (HNC) equation closure class for 3D-RISM.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  module rism3d_hnc_c
    use rism3d_potential_c
    use rism3d_grid_c
    !the HNC type
    type rism3d_hnc
       type(rism3d_potential),pointer :: pot => NULL()
       !grid : points to grid in potential object
       type(rism3d_grid),pointer :: grid => NULL()
    end type rism3d_hnc

    public rism3d_hnc_new, rism3d_hnc_destroy!, rism3d_hnc_guv
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Initializes the HNC closure
!!!IN:
!!!   this : HNC object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_hnc_new(this,pot)
      implicit none
      type(rism3d_hnc), intent(inout) :: this
      type(rism3d_potential), target, intent(in) :: pot
      this%pot => pot
      this%grid => this%pot%grid
    end subroutine rism3d_hnc_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Guv from Uuv, Huv, and Cuv using the HNC closure
!!!IN:
!!!   this : the HNC closure object
!!!   guv  : site-site pair correlation function
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_hnc_guv(this,guv, huv, cuv)
      implicit none
      type(rism3d_hnc), intent(in) :: this
      _REAL_, intent(out) :: guv(:,:)
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      integer :: iv, ir, ix, iy, iz, ig, ig1
      _REAL_ :: tuv

      do iv = 1,this%pot%solv%natom
         do iz = 1, this%grid%nr(3)
            do iy = 1, this%grid%nr(2)
               do ix = 1, this%grid%nr(1)
                  ig1 = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(2)*this%grid%nr(1)
#if defined(MPI)
                  ig = ix + (iy-1)*(this%grid%nr(1)+2) + (iz-1)*(this%grid%nr(1)+2)*this%grid%nr(2)
#else
                  ig = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(1)*this%grid%nr(2)
#endif /*defined(MPI)*/
                  tuv = -this%pot%uuv(ix,iy,iz,iv) + huv(ig,iv) - cuv(ix,iy,iz,iv)
                  guv(ig,iv) = exp(tuv)
               end do
            end do
         end do
      end do
    end subroutine rism3d_hnc_guv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Guv_dT from Uuv, Huv_dT, Cuv_dT and Guv using the HNC closure
!!!IN:
!!!   this    : the HNC closure object
!!!   guv_dT  : temperature derivative site-site pair correlation function
!!!   huv_dT  : temperature derivative site-site total correlation function
!!!   cuv_dT  : temperature derivative site-site direct correlation function
!!!   guv     : site-site pair correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_hnc_guv_dT(this,guv_dT, huv_dT, cuv_dT, guv)
      implicit none
      type(rism3d_hnc), intent(in) :: this
      _REAL_, intent(out) :: guv_dT(:,:)
      _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:), guv(:,:)
      integer :: iv, ir, ix, iy, iz, ig, ig1

      do iv = 1,this%pot%solv%natom
         do iz = 1, this%grid%nr(3)
            do iy = 1, this%grid%nr(2)
               do ix = 1, this%grid%nr(1)
                  ig1 = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(2)*this%grid%nr(1)
#if defined(MPI)
                  ig = ix + (iy-1)*(this%grid%nr(1)+2) + (iz-1)*(this%grid%nr(1)+2)*this%grid%nr(2)
#else
                  ig = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(1)*this%grid%nr(2)
#endif /*defined(MPI)*/

                  guv_dT(ig,iv) = (this%pot%uuv(ix,iy,iz,iv) + huv_dT(ig,iv) - cuv_dT(ix,iy,iz,iv))&
                       * guv(ig,iv)
               end do
            end do
         end do
      end do
    end subroutine rism3d_hnc_guv_dT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the excess chemical potential in kT for each site
!!!IN:
!!!   this : the closure object
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_hnc_exchem(this, huv, cuv) result(exchem)
      implicit none
      type(rism3d_hnc), intent(inout) :: this
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      _REAL_ :: exchem(this%pot%solv%natom)
      _REAL_ :: tuv
      integer :: ix, iy, iz, iv, ig, igk
      exchem = 0.d0
      do iv=1,this%pot%solv%natom
         do iz=1,this%grid%nr(3)
            do iy=1,this%grid%nr(2)
               do ix=1,this%grid%nr(1)
                  exchem(iv) = exchem(iv)&
                       +rism3d_hnc_exchem_IJK(this,huv,cuv,(/ix,iy,iz/),iv)
               end do
            end do
         end do
         exchem(iv) =  exchem(iv)*this%grid%voxel
      enddo

    end function rism3d_hnc_exchem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the excess chemical potential in kT at the requested grid point
!!!IN:
!!!   this : the closure object
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!   ijk  : 3d-grid index
!!!   iv   : solvent site
!!!OUT:
!!!   The contribution to the excess chemical potential from grid point i, j, k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_hnc_exchem_IJK(this, huv, cuv,ijk,iv) result(exchem)
      implicit none
      type(rism3d_hnc), intent(in) :: this
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      integer, intent(in) :: ijk(3), iv
      _REAL_ :: exchem
      _REAL_ :: tuv
      integer :: ix, iy, iz, igk
      ix=ijk(1)
      iy=ijk(2)
      iz=ijk(3)
#if defined(MPI)
      igk = ix + (iy-1)*(this%grid%nr(1)+2) + (iz-1)*this%grid%nr(2)*(this%grid%nr(1)+2)
#else
      igk = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(2)*this%grid%nr(1)
#endif /*defined(MPI)*/
      tuv = huv(igk,iv) - cuv(ix,iy,iz,iv)
      exchem = 0.5d0*huv(igk,iv)*tuv - cuv(ix,iy,iz,iv)
      exchem = exchem*this%pot%solv%rho(iv)
    end function rism3d_hnc_exchem_IJK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the short range part of the excess chemical potential w/ 
!!!asymptotic correction in kT for each site
!!!IN:
!!!   this : the closure object
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_hnc_aexchem(this, huv, cuv) result(exchem)
      implicit none
      type(rism3d_hnc), intent(inout) :: this
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      _REAL_ :: exchem(this%pot%solv%natom)
      _REAL_ :: exchemh2lr(this%pot%solv%natom),exchemhclr(this%pot%solv%natom)
      _REAL_ :: tuv, huvlr, cuvlr, tuvlr
      integer :: ix, iy, iz, iv, ig, igk

      if(.not.this%pot%solu%charged)then
         exchem = rism3d_hnc_exchem(this,huv,cuv)
         return
      end if

      exchem = 0.d0
      !
      !Long range part
      !
      call rism3d_potential_int_h2_hc(this%pot,exchemh2lr,exchemhclr)

      !
      !Short range part
      !
      huvlr=0d0
      do iv=1,this%pot%solv%natom
         do iz=1,this%grid%nr(3)
            do iy=1,this%grid%nr(2)
               do ix=1,this%grid%nr(1)
                  ig = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(2)*this%grid%nr(1)
#if defined(MPI)
                  igk = ix + (iy-1)*(this%grid%nr(1)+2) + (iz-1)*this%grid%nr(2)*(this%grid%nr(1)+2)
#else
                  igk = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(2)*this%grid%nr(1)
#endif /*defined(MPI)*/
                  tuv = huv(igk,iv) - cuv(ix,iy,iz,iv)
                  if(this%pot%solv%ionic)&
                       huvlr = this%pot%solv%charge_sp(iv)*this%pot%asymhr(ig)
                  cuvlr = this%pot%solv%charge(iv)*this%pot%asymcr(ig)
                  tuvlr = huvlr - cuvlr
                  exchem(iv) = exchem(iv) + 0.5d0*huv(igk,iv)*tuv - cuv(ix,iy,iz,iv) &
                       - (0.5d0*huvlr*tuvlr)
               end do
            end do
         end do
         exchem(iv) =  this%pot%solv%rho(iv)*(exchem(iv)*this%grid%voxel &
              +  (exchemh2lr(iv) - exchemhclr(iv))/2d0)
      enddo

    end function rism3d_hnc_aexchem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the solvation energy, dE = dmu + dTS, from the temperature derivative w/o
!!!asymptotic correction in kT for each site
!!!IN:
!!!   this    : the HNC closure object
!!!   huv_dT  : temperature derivative site-site total correlation function
!!!   cuv_dT  : temperature derivative site-site direct correlation function
!!!   huv     : site-site total correlation function
!!!   cuv     : site-site direct correlation function
!!!OUT:
!!!    solvation energy without long-range correction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_hnc_solvEne(this, huv_dT, cuv_dT, huv, cuv) result(solvEne)
      implicit none
      type(rism3d_hnc), intent(inout) :: this
      _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      _REAL_ :: solvEne(this%pot%solv%natom)
      integer :: ix, iy, iz, iv, ig, igk

      solvEne = 0.d0
      do iv=1,this%pot%solv%natom
         do iz=1,this%grid%nr(3)
            do iy=1,this%grid%nr(2)
               do ix=1,this%grid%nr(1)
!!$#if defined(MPI)
!!$                  igk = ix + (iy-1)*(this%grid%nr(1)+2) + (iz-1)*this%grid%nr(2)*(this%grid%nr(1)+2)
!!$#else
!!$                  igk = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(2)*this%grid%nr(1)
!!$#endif /*defined(MPI)*/
!!$                  solvEne(iv) = solvEne(iv) + huv(igk,iv)*huv_dT(igk,iv) &
!!$                       - (1.d0+0.5d0*huv(igk,iv))*cuv_dT(ix,iy,iz,iv) &
!!$                       - 0.5d0*cuv(ix,iy,iz,iv)*huv_dT(igk,iv)
                  solvene(iv)= solvene(iv)+&
                       rism3d_hnc_solvEne_IJK(this, huv_dT, cuv_dT, huv, cuv, (/ix,iy,iz/), iv)
               end do
            end do
         end do
         solvEne(iv) =  solvEne(iv)*this%grid%voxel
      enddo
    end function rism3d_hnc_solvEne

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the solvation energy, dE = dmu + dTS, from the
!!!temperature derivative w/o asymptotic correction in kT for the grid
!!!point.
!!!IN:
!!!   this    : the HNC closure object
!!!   huv_dT  : temperature derivative site-site total correlation function
!!!   cuv_dT  : temperature derivative site-site direct correlation function
!!!   huv     : site-site total correlation function
!!!   cuv     : site-site direct correlation function
!!!   ijk  : 3d-grid index
!!!OUT:
!!!    Contribution to the solvation energy from grid point i, j, k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_hnc_solvEne_IJK(this, huv_dT, cuv_dT, huv, cuv, ijk, iv) result(solvEne)
      implicit none
      type(rism3d_hnc), intent(in) :: this
      _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      integer, intent(in) :: ijk(3), iv
      _REAL_ :: solvEne
      integer :: ix, iy, iz, igk
      ix=ijk(1)
      iy=ijk(2)
      iz=ijk(3)
#if defined(MPI)
      igk = ix + (iy-1)*(this%grid%nr(1)+2) + (iz-1)*this%grid%nr(2)*(this%grid%nr(1)+2)
#else
      igk = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(2)*this%grid%nr(1)
#endif /*defined(MPI)*/
      solvEne = - huv(igk,iv)*huv_dT(igk,iv) &
           + (1.d0+0.5d0*huv(igk,iv))*cuv_dT(ix,iy,iz,iv) &
           + 0.5d0*cuv(ix,iy,iz,iv)*huv_dT(igk,iv)
      solvEne = solvEne *this%pot%solv%rho(iv)
    end function rism3d_hnc_solvEne_IJK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the solvation energy, dE = dmu + dTS, from the temperature derivative w/
!!!asymptotic correction in kT for each site
!!!IN:
!!!   this    : the HNC closure object
!!!   huv_dT  : temperature derivative site-site total correlation function
!!!   cuv_dT  : temperature derivative site-site direct correlation function
!!!   huv     : site-site total correlation function
!!!   cuv     : site-site direct correlation function
!!!OUT:
!!!    excess chemical potential with long-range correction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_hnc_asolvEne(this, huv_dT, cuv_dT, huv, cuv) result(solvEne)
      implicit none
      type(rism3d_hnc), intent(inout) :: this
      _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      _REAL_ :: solvEne(this%pot%solv%natom)
      _REAL_ :: solvEneh2lr(this%pot%solv%natom),solvEnehclr(this%pot%solv%natom)
      _REAL_ :: tuv, huvlr, cuvlr
      integer :: ix, iy, iz, iv, ig, igk

      if(.not.this%pot%solu%charged)then
         solvene = rism3d_hnc_solvene(this,huv_dT,cuv_dT,huv,cuv)
         return
      end if

      !
      !long range part
      !
      call rism3d_potential_int_h2_hc(this%pot,solvEneh2lr,solvEnehclr)

      !
      !short range part
      !
      solvEne = 0.d0
      huvlr=0d0
      do iv=1,this%pot%solv%natom
         do iz=1,this%grid%nr(3)
            do iy=1,this%grid%nr(2)
               do ix=1,this%grid%nr(1)
                  ig = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(2)*this%grid%nr(1)
#if defined(MPI)
                  igk = ix + (iy-1)*(this%grid%nr(1)+2) + (iz-1)*this%grid%nr(2)*(this%grid%nr(1)+2)
#else
                  igk = ix + (iy-1)*this%grid%nr(1) + (iz-1)*this%grid%nr(2)*this%grid%nr(1)
#endif /*defined(MPI)*/
                  tuv = huv(igk,iv) - cuv(ix,iy,iz,iv)
                  if(this%pot%solv%ionic)&
                       huvlr = this%pot%solv%charge_sp(iv)*this%pot%asymhr(ig)
                  cuvlr = this%pot%solv%charge(iv)*this%pot%asymcr(ig)
                  solvEne(iv) = solvEne(iv) + huv(igk,iv)*huv_dT(igk,iv) &
                                                 - (1.d0+0.5d0*huv(igk,iv))*cuv_dT(ix,iy,iz,iv) &
                                                 - 0.5d0*cuv(ix,iy,iz,iv)*huv_dT(igk,iv)

!!??                     solvEne(iv) = solvEne(iv) + huvlr*huvlr - huvlr*cuvlr - cuvlr
                  solvEne(iv) = solvEne(iv) + huvlr*huvlr - huvlr*cuvlr
               end do
            end do
         end do
         solvEne(iv) =  -1d0*(this%pot%solv%rho(iv)*(solvEne(iv)*this%grid%voxel &
              -  (solvEneh2lr(iv) - solvEnehclr(iv))))
      enddo

    end function rism3d_hnc_asolvEne

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Frees memory and resets the HNC closure
!!!IN:
!!!   this : HNC object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_hnc_destroy(this)
      implicit none
      type(rism3d_hnc), intent(inout) :: this
      nullify(this%pot)
      nullify(this%grid)
    end subroutine rism3d_hnc_destroy
  end module rism3d_hnc_c
