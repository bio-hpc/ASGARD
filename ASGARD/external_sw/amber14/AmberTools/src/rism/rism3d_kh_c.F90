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
!!!Kovalenko-Hirata closure class for 1D-RISM.
!!!Kovalenko, A.; Hirata, F. J. Chem. Phys. 1999, 110, 10095–10112.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  module rism3d_kh_c
    use rism3d_potential_c
    use rism3d_grid_c
    !the KH type
    type rism3d_kh
       type(rism3d_potential),pointer :: pot => NULL()
       !grid : points to grid in potential object
       type(rism3d_grid),pointer :: grid => NULL()
    end type rism3d_kh

    public rism3d_kh_new, rism3d_kh_destroy!, rism3d_kh_guv
  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Initializes the KH closure
!!!IN:
!!!   this : KH object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_kh_new(this,pot)
      implicit none
      type(rism3d_kh), intent(inout) :: this
      type(rism3d_potential), target, intent(in) :: pot
      this%pot => pot
      this%grid => this%pot%grid
    end subroutine rism3d_kh_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Guv from Uuv, Huv, and Cuv using the KH closure
!!!IN:
!!!   this : the KH closure object
!!!   guv  : site-site pair correlation function
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_kh_guv(this,guv, huv, cuv)
      implicit none
      type(rism3d_kh), intent(in) :: this
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
                  if(tuv >= 0d0)then
                     guv(ig,iv) = 1d0 + tuv
                  else
                     guv(ig,iv) = exp(tuv)
                  endif
               end do
            end do
         end do
      end do
    end subroutine rism3d_kh_guv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Guv_dT from Uuv, Huv_dT, Cuv_dT and Guv using the KH closure
!!!IN:
!!!   this    : the KH closure object
!!!   guv_dT  : temperature derivative site-site pair correlation function
!!!   huv_dT  : temperature derivative site-site total correlation function
!!!   cuv_dT  : temperature derivative site-site direct correlation function
!!!   guv     : site-site pair correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_kh_guv_dT(this,guv_dT, huv_dT, cuv_dT, guv)
      implicit none
      type(rism3d_kh), intent(in) :: this
      _REAL_, intent(out) :: guv_dT(:,:)
      _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:), guv(:,:)
      integer :: iv, ir, ix, iy, iz, ig, ig1
      _REAL_ :: tuv_dT

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

                  tuv_dT = huv_dT(ig,iv) - cuv_dT(ix,iy,iz,iv)

                  if(guv(ig,iv) >= 1d0)then
                     guv_dT(ig,iv) = this%pot%uuv(ix,iy,iz,iv) + tuv_dT
                  else
                     guv_dT(ig,iv) = (this%pot%uuv(ix,iy,iz,iv) + tuv_dT) * guv(ig,iv)
                  endif
               end do
            end do
         end do
      end do
    end subroutine rism3d_kh_guv_dT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the excess chemical potential in kT for each site
!!!IN:
!!!   this : the closure object
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_kh_exchem(this, huv, cuv) result(exchem)
      implicit none
      type(rism3d_kh), intent(inout) :: this
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      _REAL_ :: exchem(this%pot%solv%natom)
      _REAL_ :: tuv
      integer :: ix, iy, iz, iv, ig, igk
      exchem = 0.d0
      do iv=1,this%pot%solv%natom
         do iz=1,this%grid%nr(3)
            do iy=1,this%grid%nr(2)
               do ix=1,this%grid%nr(1)
                  exchem(iv) = exchem(iv)+rism3d_kh_exchem_IJK(this,huv,cuv,(/ix,iy,iz/),iv)
               end do
            end do
         end do
         exchem(iv) =  exchem(iv)*this%grid%voxel
      enddo
    end function rism3d_kh_exchem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the excess chemical potential in kT at the requested grid point
!!!IN:
!!!   this : the closure object
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!   ijk  : 3d-grid index
!!!OUT:
!!!   The contribution to the excess chemical potential from grid point i, j, k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_kh_exchem_IJK(this, huv, cuv,ijk,iv) result(exchem)
      implicit none
      type(rism3d_kh), intent(in) :: this
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
      if (huv(igk,iv) > 0d0)  then
         exchem = -(1.d0+0.5d0*huv(igk,iv))*cuv(ix,iy,iz,iv)
      else
         exchem = (0.5d0*huv(igk,iv)*tuv - cuv(ix,iy,iz,iv))
      endif
      exchem = exchem*this%pot%solv%rho(iv)
    end function rism3d_kh_exchem_IJK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the short range part of the excess chemical potential w/ 
!!!asymptotic correction in kT for each site
!!!IN:
!!!   this : the closure object
!!!   huv  : site-site total correlation function
!!!   cuv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_kh_aexchem(this, huv, cuv) result(exchem)
      implicit none
      type(rism3d_kh), intent(inout) :: this
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      _REAL_ :: exchem(this%pot%solv%natom)
      _REAL_ :: exchemh2lr(this%pot%solv%natom),exchemhclr(this%pot%solv%natom)
      _REAL_ :: tuv, huvlr, cuvlr
      integer :: ix, iy, iz, iv, ig, igk

      if(.not.this%pot%solu%charged)then
         exchem = rism3d_kh_exchem(this,huv,cuv)
         return
      end if

      !
      !long range part
      !
      call rism3d_potential_int_h2_hc(this%pot,exchemh2lr,exchemhclr)
      do iv=1,this%pot%solv%natom
         if(this%pot%solu%totcharge*this%pot%solv%charge_sp(iv) > 0.d0) &
!!         if(this%pot%solu%totcharge*this%pot%solv%charge_sp(iv) <= 0.d0) &
              exchemh2lr(iv)=0
      end do

      !
      !short range part
      !
      exchem = 0.d0
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
                  if (huv(igk,iv) > 0d0)  then
                     exchem(iv) = exchem(iv) - (1.d0+0.5d0*huv(igk,iv))*cuv(ix,iy,iz,iv)
                  else
                     exchem(iv) = exchem(iv) + 0.5d0*huv(igk,iv)*tuv - cuv(ix,iy,iz,iv)
                  endif

                  if (this%pot%solu%totcharge*this%pot%solv%charge_sp(iv) <= 0.d0) then
!!$                     exchem(iv) = exchem(iv) + 0.5d0*huvlr*cuvlr + cuvlr
                     exchem(iv) = exchem(iv) + 0.5d0*huvlr*cuvlr
                  else
!!$                     exchem(iv) = exchem(iv) - 0.5d0*huvlr*huvlr + 0.5d0*huvlr*cuvlr + cuvlr
                     exchem(iv) = exchem(iv) - 0.5d0*huvlr*huvlr + 0.5d0*huvlr*cuvlr
                  endif
               end do
            end do
         end do
         exchem(iv) =  this%pot%solv%rho(iv)*(exchem(iv)*this%grid%voxel &
              +  (exchemh2lr(iv) - exchemhclr(iv))/2d0)
      enddo

    end function rism3d_kh_aexchem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the solvation energy, dE = dmu + dTS, from the temperature derivative w/o
!!!asymptotic correction in kT for each site
!!!IN:
!!!   this    : the KH closure object
!!!   huv_dT  : temperature derivative site-site total correlation function
!!!   cuv_dT  : temperature derivative site-site direct correlation function
!!!   huv     : site-site total correlation function
!!!   cuv     : site-site direct correlation function
!!!OUT:
!!!    solvation energy without long-range correction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_kh_solvEne(this, huv_dT, cuv_dT, huv, cuv) result(solvEne)
      implicit none
      type(rism3d_kh), intent(inout) :: this
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
!!$                  if (huv(igk,iv) > 0d0)  then
!!$                     solvEne(iv) = solvEne(iv) + (1.d0+0.5d0*huv(igk,iv))*cuv_dT(ix,iy,iz,iv) &
!!$                                                 + 0.5d0*cuv(ix,iy,iz,iv)*huv_dT(igk,iv)
!!$                  else
!!$                     solvEne(iv) = solvEne(iv) - huv(igk,iv)*huv_dT(igk,iv) &
!!$                                                 + (1.d0+0.5d0*huv(igk,iv))*cuv_dT(ix,iy,iz,iv) &
!!$                                                 + 0.5d0*cuv(ix,iy,iz,iv)*huv_dT(igk,iv)
!!$                  endif
                  solvene(iv)= solvene(iv)+&
                       rism3d_kh_solvEne_IJK(this, huv_dT, cuv_dT, huv, cuv, (/ix,iy,iz/), iv)
               end do
            end do
         end do
         solvEne(iv) =  solvEne(iv)*this%grid%voxel
      enddo
    end function rism3d_kh_solvEne

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the solvation energy, dE = dmu + dTS, from the
!!!temperature derivative w/o asymptotic correction in kT for the grid
!!!point.
!!!IN:
!!!   this    : the KH closure object
!!!   huv_dT  : temperature derivative site-site total correlation function
!!!   cuv_dT  : temperature derivative site-site direct correlation function
!!!   huv     : site-site total correlation function
!!!   cuv     : site-site direct correlation function
!!!   ijk  : 3d-grid index
!!!OUT:
!!!    Contribution to the solvation energy from grid point i, j, k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_kh_solvEne_IJK(this, huv_dT, cuv_dT, huv, cuv, ijk, iv) result(solvEne)
      implicit none
      type(rism3d_kh), intent(in) :: this
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
      if (huv(igk,iv) > 0d0)  then
         solvEne =   (1.d0+0.5d0*huv(igk,iv))*cuv_dT(ix,iy,iz,iv) &
              + 0.5d0*cuv(ix,iy,iz,iv)*huv_dT(igk,iv)
      else
         solvEne =  -huv(igk,iv)*huv_dT(igk,iv) &
              + (1.d0+0.5d0*huv(igk,iv))*cuv_dT(ix,iy,iz,iv) &
              + 0.5d0*cuv(ix,iy,iz,iv)*huv_dT(igk,iv)
      endif
      solvEne = solvEne*this%pot%solv%rho(iv)
    end function rism3d_kh_solvEne_IJK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the solvation energy, dE = dmu + dTS, from the temperature derivative w/
!!!asymptotic correction in kT for each site
!!!IN:
!!!   this    : the KH closure object
!!!   huv_dT  : temperature derivative site-site total correlation function
!!!   cuv_dT  : temperature derivative site-site direct correlation function
!!!   huv     : site-site total correlation function
!!!   cuv     : site-site direct correlation function
!!!OUT:
!!!    excess chemical potential with long-range correction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_kh_asolvEne(this, huv_dT, cuv_dT, huv, cuv) result(solvEne)
      implicit none
      type(rism3d_kh), intent(inout) :: this
      _REAL_, intent(in) :: huv_dT(:,:),cuv_dT(:,:,:,:)
      _REAL_, intent(in) :: huv(:,:),cuv(:,:,:,:)
      _REAL_ :: solvEne(this%pot%solv%natom)
      _REAL_ :: solvEneh2lr(this%pot%solv%natom),solvEnehclr(this%pot%solv%natom)
      _REAL_ :: tuv, huvlr, cuvlr
      integer :: ix, iy, iz, iv, ig, igk

      if(.not.this%pot%solu%charged)then
         solvene = rism3d_kh_solvene(this,huv_dT,cuv_dT,huv,cuv)
         return
      end if

      !
      !long range part
      !
      call rism3d_potential_int_h2_hc(this%pot,solvEneh2lr,solvEnehclr)
!      print *,'solvEneh2lr',solvEneh2lr
!      print *,'solvEnehclr',solvEnehclr
      do iv=1,this%pot%solv%natom
!!??         if(this%pot%solu%totcharge*this%pot%solv%charge_sp(iv) > 0.d0) &
         if(this%pot%solu%totcharge*this%pot%solv%charge_sp(iv) <= 0.d0) &
              solvEneh2lr(iv)=0
      end do

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
                  if (huv(igk,iv) > 0d0)  then
                     solvEne(iv) = solvEne(iv) - (1.d0+0.5d0*huv(igk,iv))*cuv_dT(ix,iy,iz,iv) &
                                                 - 0.5d0*cuv(ix,iy,iz,iv)*huv_dT(igk,iv)
                  else
                     solvEne(iv) = solvEne(iv) + huv(igk,iv)*huv_dT(igk,iv) &
                                                 - (1.d0+0.5d0*huv(igk,iv))*cuv_dT(ix,iy,iz,iv) &
                                                 - 0.5d0*cuv(ix,iy,iz,iv)*huv_dT(igk,iv)
                  endif

                  if (this%pot%solu%totcharge*this%pot%solv%charge_sp(iv) <= 0.d0) then
!!??                     solvEne(iv) = solvEne(iv) - huvlr*cuvlr - cuvlr
                     solvEne(iv) = solvEne(iv) - huvlr*cuvlr
                  else
!!??                     solvEne(iv) = solvEne(iv) + huvlr*huvlr - huvlr*cuvlr - cuvlr
                     solvEne(iv) = solvEne(iv) + huvlr*huvlr - huvlr*cuvlr
                  endif
               end do
            end do
         end do
         solvEne(iv) =  -1d0*(this%pot%solv%rho(iv)*(solvEne(iv)*this%grid%voxel &
              -  (solvEneh2lr(iv) - solvEnehclr(iv))))
      enddo

    end function rism3d_kh_asolvEne

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Frees memory and resets the KH closure
!!!IN:
!!!   this : KH object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_kh_destroy(this)
      implicit none
      type(rism3d_kh), intent(inout) :: this
      nullify(this%pot)
      nullify(this%grid)
    end subroutine rism3d_kh_destroy
  end module rism3d_kh_c
