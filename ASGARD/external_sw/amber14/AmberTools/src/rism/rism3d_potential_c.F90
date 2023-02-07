!<compile=optimized>

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Potential class for 3D-RISM. Used to calculate/store quantities that are
!potential dependent and do not change while converging a solution calculation.
!
!Pointers to solute and solvent objects are maintained.  So, if values change in
!these objects, they are automatically used in the potential calculation.
!
!This class is generally MPI agnostic.  That is, MPI is only an issue for setting
!the grid size where both information about the size of the local slab and the
!total grid must be supplied.  There is no MPI communication within the class.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "../include/dprec.fh"

  module rism3d_potential_c
    use safemem
    use rism3d_solu_c
    use rism3d_solv_c
    use rism3d_grid_c
    use rism_timer_c
#include "def_time.h"
    type rism3d_potential
       !cut2 :: cutoff**2 for rism potential calculations
       !cut_hlr :: cutoff for long range asymptotics of h
       !cut_clr :: cutoff for long range asymptotics of c
       _REAL_ :: cut2, cut_hlr, cut_clr

       !grid :: pointer to grid object
       type(rism3d_grid), pointer :: grid=>NULL()

       !solu :: pointer to solute object
       !solv :: pointer to solvent object
       type(rism3d_solu), pointer :: solu=>NULL()
       type(rism3d_solv), pointer :: solv=>NULL()

       !timer        : total time for potential and asymptotics
       !ljTimer      : time for LJ potential
       !coulombTimer : timer for Coulomb potential
       !asympTimer   : timer for asymptotics
       type(rism_timer) :: timer, ljTimer, coulombTimer, asympTimer

       !uuv        :: potential energy of the solvent about the solute.  This is recalculated for each 
       !              solution but  we want to reserve the memory and may want to change it if the box 
       !              dimensions change [kT]
       _REAL_, pointer :: uuv(:,:,:,:)=>NULL()
       !siguv :: solute-solvent sigma interaction matrix.  Calculated once and 
       !         used in later calls [A]
       _REAL_, pointer :: siguv(:,:) => NULL()
       !epsuv :: solute-solvent epsilon interaction matrix.  Calculated once and
       !          used in later calls [kT]
       _REAL_, pointer :: epsuv(:,:) => NULL()

       !Asymptotics are part of the supercell formalism and, as this
       !may not be used, asymptotics grids are only allocated if/when
       !the asymptotics calcluation is requested. In particular, if
       !the solvent is not ionic, asymh[rk] is not allocated.  This
       !should be transparent to the user
       !asymcr     :: supercell real space asymptotics for C
       !asymck     :: supercell k space asymptotics for C
       !asymhr     :: supercell real space asymptotics for H
       !asymhk     :: supercell k space asymptotics for H
       _REAL_, pointer :: asymcr(:)=>NULL(),asymck(:)=>NULL(),asymhr(:)=>NULL(),&
            asymhk(:)=>NULL()
       
       !huvk0 :: (2,solv%natom) long-range part of Huv(k) at k=0.
       !huvk0_dT :: (2,solv%natom) long-range part of temperature derivative Huv(k) at k=0.
       _REAL_, pointer :: huvk0(:,:)=>NULL(), huvk0_dT(:,:)=>NULL()

   end type rism3d_potential

    private prep_UV, uljuv, ucoulu

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Constructor.
!!!IN:
!!!   this   : potential object
!!!   grid   : grid object.  A pointer to this will be retained
!!!   solv   : solvent object.  A pointer to this will be retained
!!!   solu   : solute object.  A pointer to this will be retained
!!!   ng3    : number of grid points in each dimension.  If this is an MPI 
!!!            calculation, pass the local grid size
!!!   cut   : cutoff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_potential_new(this,grid,solv,solu,cut)
      implicit none
      type(rism3d_potential), intent(inout) :: this
      type(rism3d_grid), target, intent(in) :: grid
      type(rism3d_solu), target, intent(in) :: solu
      type(rism3d_solv), target, intent(in) :: solv
      _REAL_, intent(in):: cut
      this%grid=>grid
      this%solv=>solv
      this%solu=>solu
      call rism3d_potential_setCut(this,cut)
      this%siguv => safemem_realloc(this%siguv,this%solu%natom,this%solv%natom,.false.)
      this%epsuv => safemem_realloc(this%epsuv,this%solu%natom,this%solv%natom,.false.)
      call prep_UV(this)
      call rism_timer_new(this%timer,"Potential")
      call rism_timer_new(this%ljTimer,"Lennard-Jones")
      call rism_timer_setParent(this%ljTimer,this%timer)
      call rism_timer_new(this%coulombTimer,"Coulomb")
      call rism_timer_setParent(this%coulombTimer,this%timer)
      call rism_timer_new(this%asympTimer,"Asymptotics")
      call rism_timer_setParent(this%asympTimer,this%timer)
    end subroutine rism3d_potential_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Set parent for this timer
!!!IN:
!!!   this : rism3d_potential object
!!!   parent : parent timer object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_potential_setTimerParent(this, parent)
    implicit none
    type(rism3d_potential), intent(inout) :: this
    type(rism_timer), intent(inout) :: parent
    call rism_timer_start(this%timer)
    call rism_timer_setParent(this%timer,parent)
    call rism_timer_stop(this%timer)
  end subroutine rism3d_potential_setTimerParent

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets the cut off distance for potential, force and long range asymptotics 
!!!calculations.  Note that the long range asymptotics cut off is calculated 
!!!from the solvent parameters    
!!!IN:
!!!   this :: potential object
!!!   cut  :: distance cutoff for potential and force calculations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_potential_setcut(this,cut)
      implicit none
      type(rism3d_potential),intent(inout) :: this
      _REAL_, intent(in) :: cut

      !exp_quarter_xappa2_smear2 : exp(-1/4*(kappa^2)*(smear^2))
      !half_xappa_smear      : 1/2* kappa   * smear
      !asymhr_coeff_at_cut   : long range h without charge or 1/r dependence (asymhr/(Qq)*r)
      !asymcr_coeff_at_cut   : long range c without charge or 1/r dependence (asymcr/(Qq)*r)
      !tcut                  : trial cutoff
      _REAL_ :: exp_quarter_xappa2_smear2, asymhr_coeff_at_cut, asymcr_coeff_at_cut,&
           half_xappa_smear, tcut

      _REAL_,external ::      erfc

      !assign potential cutoff
      !ensure that the square won't overflow
      this%cut2 = min(sqrt(huge(1d0)),cut)
      this%cut2 = this%cut2**2

      half_xappa_smear = 0.5d0*this%solv%xappa*this%solv%smear
      exp_quarter_xappa2_smear2 = exp(half_xappa_smear**2)


      !This is really a very simpled minded approach to finding the cut of distance.
      !Basically, it is a linear search for a value where the non-1/r part is less
      !than an error value.  This really should be done with Newton-Rhapson instead.
      tcut=4d0
      this%cut_hlr=-1d0
      this%cut_clr=-1d0
      do while((this%cut_hlr <0 .or. this%cut_clr<0) .and. tcut < 1000d0)
         asymhr_coeff_at_cut=exp_quarter_xappa2_smear2 &
              *(exp(-this%solv%xappa*tcut)*erfc(half_xappa_smear - tcut/this%solv%smear) &
              - exp(this%solv%xappa*tcut)*erfc(half_xappa_smear+tcut/this%solv%smear))/2d0
         asymcr_coeff_at_cut=(1d0-erfc(tcut/this%solv%smear))
         if(exp_quarter_xappa2_smear2-asymhr_coeff_at_cut < 1d-7 .and. this%cut_hlr<0) &
              this%cut_hlr = tcut
         if(1d0-asymcr_coeff_at_cut < 1d-7 .and. this%cut_clr<0) &
              this%cut_clr = tcut
         tcut = tcut+.1d0
      end do

!!$      if(this%cut_hlr<0) this%cut_hlr = HUGE(1d0)
!!$      if(this%cut_clr<0) this%cut_clr = HUGE(1d0)
      this%cut_hlr = HUGE(1d0)
      this%cut_clr = HUGE(1d0)
    end subroutine rism3d_potential_setcut

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the potential on to the grid.
!!!IN:
!!!   this   : potential object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_potential_calc(this)
      use rism_util, only : checksum
      implicit none
      type(rism3d_potential), intent(inout) :: this
      call rism_timer_start(this%timer)
      !ensure a grid size has been set
      if(.not.associated(this%grid%gv))then
         call rism_report_error("rism3d_potential_calc: grid size not set")
         stop
      end if
      !check if the grid size has changed
      if(ubound(this%uuv,1) /= this%grid%nr(1) .or.&
         ubound(this%uuv,2) /= this%grid%nr(2) .or.&
         ubound(this%uuv,3) /= this%grid%nr(3) .or.&
         .not. associated(this%uuv))then
         this%uuv => safemem_realloc(this%uuv,this%grid%nr(1),this%grid%nr(2), this%grid%nr(3),&
              this%solv%natom,.false.)
      end if

      call timer_start(TIME_UCOULU)
      call rism_timer_start(this%coulombTimer)
      this%uuv=0
      if(this%solu%charged)&
           call  ucoulu (this, this%uuv)
      call rism_timer_stop(this%coulombTimer)
      call timer_stop(TIME_UCOULU)

      call timer_start(TIME_ULJUV)
      call rism_timer_start(this%ljTimer)
      call  uljuv (this, this%uuv)
      call rism_timer_stop(this%ljTimer)
      call timer_stop(TIME_ULJUV)
      call rism_timer_stop(this%timer)
    end subroutine rism3d_potential_calc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Set r-space supercell asymptotic function H
!!!IN:
!!!  this  :: rism3d potential object
!!!  ahr   :: pointer to pre-calculated long range asymptotics.  This class will
!!!           try to deallocate this memory when destroyed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  rism3d_potential_setAsymph(this,ahr)
      implicit none
      type(rism3d_potential), intent(inout) :: this
      _REAL_, target, intent(in) :: ahr(:)
      this%asymhr => ahr
    end subroutine rism3d_potential_setAsymph

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Set r-space supercell asymptotic function C
!!!IN:
!!!  this  :: rism3d potential object
!!!  acr   :: pointer to pre-calculated long range asymptotics.  This class will
!!!           try to deallocate this memory when destroyed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  rism3d_potential_setAsympc(this,acr)
      implicit none
      type(rism3d_potential), intent(inout) :: this
      _REAL_, target, intent(in) :: acr(:)
      this%asymcr => acr
    end subroutine rism3d_potential_setAsympc
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Synthesizes supercell asymptotic function of C and H
!!!IN:
!!!  this  :: rism3d potential object
!!!  huvk0 :: (2,solv%natom) long-range part of Huv(k) at k=0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine  rism3d_potential_asympch(this)
      use constants, only : PI, FOURPI
      use rism_util, only : checksum
      implicit none
#if defined(MPI)
      include 'mpif.h'
#endif /*defined(MPI)*/
      type(rism3d_potential), intent(inout) :: this
      _REAL_ :: uc1gc,uc1gh
      integer ::  igx,igy,igz, ig, iu,iv, ig0
      _REAL_ ::  phase,sumcos,sumsin, rx,ry,rz, delx,dely,delz,delyz2, ra, &
           uc1g,sr,uc,ucs,ucsc
      _REAL_ :: xappa2,smear2_4
      !offset :: z-axis offset (important for spatially distributed MPI)
      _REAL_ :: offset

!      if(.not.this%solu%charged) return

      call rism_timer_start(this%asympTimer)

      !ensure a grid size has been set
      if(.not.associated(this%grid%gv))then
         call rism_report_error("rism3d_potential_asympch: grid size not set")
      end if
      
      !allocate memory if the grid size has changed
      if(this%grid%nkTotal /= ubound(this%asymck,1) .or. &
           product(this%grid%nr) /= ubound(this%asymcr,1) .or. &
           .not. associated(this%asymcr) .or. &
           .not. associated(this%asymhr))then
         this%asymcr => safemem_realloc(this%asymcr,product(this%grid%nr),.false.)
         this%asymck => safemem_realloc(this%asymck,this%grid%nkTotal,.false.)
         if(this%solv%ionic)then
            this%asymhr => safemem_realloc(this%asymhr,product(this%grid%nr),.false.)
            this%asymhk => safemem_realloc(this%asymhk,this%grid%nkTotal,.false.)
         end if
      end if
      !allocate long range part of Huv(k) if not done already
      if(.not.associated(this%huvk0))then
         this%huvk0 => safemem_realloc(this%huvk0, 2, this%solv%natom)
         this%huvk0_dT => safemem_realloc(this%huvk0_dT, 2, this%solv%natom)
      end if

      smear2_4 = this%solv%smear**2/4d0
      xappa2=this%solv%xappa**2

      !.....initialize
      this%asymck=0d0
      if(this%solv%ionic)&
           this%asymhk=0d0

      !.....asymptotic functions for Cuv and Huv in R-space
      !..................... enumerating box grid points .....................
      offset = this%grid%grdspc(3)*(this%grid%nrOff(3))
      do igz=0,this%grid%nr(3)-1
         rz = igz*this%grid%grdspc(3)+offset
         do igy=0,this%grid%nr(2)-1
            ry = igy*this%grid%grdspc(2)
            do igx=0,this%grid%nr(1)-1
               ig = 1 + igx + igy*this%grid%nr(1) + igz*this%grid%nr(2)*this%grid%nr(1)
               rx = igx*this%grid%grdspc(1)
               if(this%solv%ionic)&
                    this%asymhr(ig) =  hlr(this,(/rx,ry,rz/))
               this%asymcr(ig) =  clr(this,(/rx,ry,rz/))
            enddo
         enddo
      enddo

      !.....getting long-range part of Cuv and Huv in K-space.

      !k=0 is on the master node. Set it to zero as it is handled separately
      if(this%grid%nkOff(3)==0)then
         this%asymck(1) = 0d0
         this%asymck(2) = 0d0
         if(this%solv%ionic)then
            this%asymhk(1) = 0d0
            this%asymhk(2) = 0d0
         end if
         ig0=2
      else
         ig0=1
      end if
      do ig=ig0,this%grid%nkTotal/2
         sumcos = 0d0
         sumsin = 0d0
         do iu=1,this%solu%natom
            phase =  dot_product(this%grid%gv(:,ig),this%solu%ratu(:,iu))
            sumcos = sumcos + this%solu%charge(iu)*cos(phase)
            sumsin = sumsin + this%solu%charge(iu)*sin(phase)
         enddo
         !.....c in k-space
         uc1g = FOURPI * exp( -smear2_4*this%grid%g2(ig))
         uc1gc = uc1g/this%grid%g2(ig)
         this%asymck(2*ig-1) = uc1gc * sumcos/this%grid%boxvol
         this%asymck(2*ig)   = uc1gc * sumsin/this%grid%boxvol

         !.....h in k-space
         if(this%solv%ionic)then
            uc1gh = uc1g/(this%grid%g2(ig) + xappa2)
            this%asymhk(2*ig)   = uc1gh * sumsin/this%grid%boxvol
            this%asymhk(2*ig-1) = uc1gh * sumcos/this%grid%boxvol
         end if
      enddo

      if(this%grid%nkOff(3)==0)then
         !.....getting the difference between
         !.....long-range asymptotic functions of huv at k=0.
         sumcos = 0d0
         sumsin = 0d0
         do iu=1,this%solu%natom
            phase =  dot_product(this%grid%gv(:,1),this%solu%ratu(:,iu))
            sumcos = sumcos + this%solu%charge(iu)*cos(phase)
            sumsin = sumsin + this%solu%charge(iu)*sin(phase)
         enddo
         do iv=1,this%solv%natom
            this%huvk0(1,iv) = this%solv%delhv0(iv) * sumcos / product(this%grid%boxlen)
            this%huvk0(2,iv) = this%solv%delhv0(iv) * sumsin / product(this%grid%boxlen)
            !make sure we have read the information to do this
            if(this%solv%delhv0_dT(1) /= huge(1d0))then
               this%huvk0_dT(1,iv) = this%solv%delhv0_dT(iv) * sumcos / product(this%grid%boxlen)
               this%huvk0_dT(2,iv) = this%solv%delhv0_dT(iv) * sumsin / product(this%grid%boxlen)
            end if
         enddo
      end if
      call rism_timer_stop(this%asympTimer)
    end subroutine rism3d_potential_asympch

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates integrals of the long range asymptotics h**2 and h*c for all space
!!!(i.e. not just the grid).  These terms are used by the HNC term in the 
!!!chemical potential so we calculate them together.
!!!IN:
!!!   this : the closure object
!!!   h2   : integral of asymhr**2
!!!   hc   : integral of asymhr*asymcr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine rism3d_potential_int_h2_hc(this,h2,hc)
      use constants, only : PI
      use rism_util, only : gaussquad_legendre, checksum
      implicit none
      type(rism3d_potential), intent(in) :: this
      _REAL_, intent(out) :: h2(this%solv%natom), hc(this%solv%natom)

      _REAL_ :: sumhc,sumh2,sumb,xarg,x2arg,denom
      _REAL_ :: dx, dy , dz, r2, k, k2
      integer :: i, j, ik, iv , n0, nf
      integer, parameter :: Nintmx = 200
      _REAL_ :: argum(Nintmx),weight(Nintmx)

      if(.not.this%solu%charged .or. .not.this%solv%ionic)then
         h2=0d0
         hc=0d0
         return
      end if

      !get the numerical long-range contribution.  This uses Gauss-Legendre
      !quadrature where the integration bounds have been converted from 
      !0 to infinity to 0 to 1. 
      !See Genheden et al. J. Phys. Chem. B 114, 8505 (2010) Eq. 14

      !.....initialize gauss-legendre integration
      call gaussquad_legendre (0d0,1d0,argum,weight,Nintmx)
      sumhc = 0.d0
      sumh2 = 0.d0

      !.....k loop

      !parallelize over k-loop.  This is good for a maximum of Nintmx processes.
      !select range of sum for this process      
      !get the start point for this process
      n0 = Nintmx/this%grid%mpisize*this%grid%mpirank +1
      !the final point is the start point for one process
      !higher.  Ensure we do not exceed the desired range of the
      !sum
      nf = min(Nintmx/this%grid%mpisize*(this%grid%mpirank+1), Nintmx)
      do ik = n0, nf

         k = argum(ik)/(1.d0-argum(ik))

         !.....Bessel part
         sumb = 0.d0

         do i=2,this%solu%natom
            do j=1,i-1

               !......... site separation  ..........
               dx = this%solu%ratu(1,i) - this%solu%ratu(1,j)
               dy = this%solu%ratu(2,i) - this%solu%ratu(2,j)
               dz = this%solu%ratu(3,i) - this%solu%ratu(3,j)

               r2 = dx*dx + dy*dy + dz*dz

               xarg = this%solv%xappa*k*sqrt(r2)
               if(xarg == 0.d0) &
                    sumb = sumb + 1.d0*this%solu%charge(i)*this%solu%charge(j)
               if(xarg /= 0.d0) &
                    sumb = sumb + sin(xarg)/xarg *this%solu%charge(i)*this%solu%charge(j)

            enddo
         enddo

         sumb = sumb*2.d0

         do i=1,this%solu%natom
            sumb = sumb + 1.d0*this%solu%charge(i)*this%solu%charge(i)
         enddo

         !.....end of Bessel part

         x2arg = (k*this%solv%xappa*this%solv%smear)**2
         k2 = k**2

         denom = argum(ik)**2+(1.0d0-argum(ik))**2

         sumhc = sumhc - &
              exp(-x2arg/2.d0)/denom*sumb*weight(ik)
         sumh2 = sumh2 + &
              exp(-x2arg/2.d0)*argum(ik)**2/(denom**2)*sumb*weight(ik)

      enddo
      !.....end of k-loop

      !.....site xappa
      do iv = 1,this%solv%natom
         h2(iv) = 8.d0*PI/(this%solv%dielconst) &
              *this%solv%charge_sp(iv)**2
         hc(iv) = 8.d0*PI/(this%solv%dielconst) &
              *this%solv%charge(iv)*this%solv%charge_sp(iv)
      enddo

      !.....unit for xappa
      do iv = 1,this%solv%natom
         h2(iv) = h2(iv)
         hc(iv) = hc(iv)
      enddo

      !.....divided by total xappa
      if(this%solv%xappa /= 0.d0)then
         do iv = 1,this%solv%natom
            h2(iv) = h2(iv) / this%solv%xappa
            hc(iv) = hc(iv) / this%solv%xappa
         enddo
      else
         do iv = 1,this%solv%natom
            h2(iv) = 0.d0
            hc(iv) = 0.d0
         enddo
      endif

      !.....multiplying xmulr by integral
      do iv = 1,this%solv%natom
         h2(iv) =   h2(iv)*sumh2
         hc(iv) = - hc(iv)*sumhc
      enddo

      !.....J devided by KT
      do iv = 1,this%solv%natom
         h2(iv) = h2(iv) / (PI*this%solv%dielconst)
         hc(iv) = hc(iv) / PI
      enddo
    end subroutine rism3d_potential_int_h2_hc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates integral of the long range asymptotics of h**n for all space
!!!(i.e. not just the grid).
!!!IN:
!!!   this  : the closure object
!!!   power : power to multiply t* to
!!!OUT:
!!!   integral of asymhr**n
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_potential_int_hn(this,power) result(integral)
      use constants, only : PI
      use rism_util, only : gaussquad_laguerre
      implicit none
      type(rism3d_potential), intent(in) :: this
      integer, intent(in) :: power
      _REAL_ :: integral(this%solv%natom)
      _REAL_ :: temp(this%solv%natom)
      _REAL_ :: total, shell, r(3)
      _REAL_, parameter :: err =1d-8
      _REAL_ :: sm, tstar
      integer :: iu, iv, istep, ishell, ix, iy, iz, ig, n0, nf, npt
      !nstep :: steps for Gauss-Laguerre quadrature
      integer, parameter :: nstep = 100, maxshell = 1000
      _REAL_ :: argum(nstep),weight(nstep)

      integral = 0d0
      if(.not.this%solv%ionic .or. .not.this%solu%charged) return

      if(power == 2)then
         call rism3d_potential_int_h2_hc(this,integral,temp)
      elseif(power == 1) then
         integral = rism3d_potential_int_h(this)
      else
         !since this integral converges faster than exp(r), we can use Gauss-Laguerre
         !quadrature
!!$      call gaussquad_laguerre(argum,weight,nstep)

         !Instead use a brute force grid based method for now
         !first sum the precalculated values
         total = sum(this%asymhr**power)

         !add sum over shells around the original grid until the
         !relative error falls below our selected threshold The shells
         !are constructed as faces of the grid expanded by one point
         !at each boundary. The faces are two points larger along one
         !axis, such that the faces 'interlock' to cover the edges.
         !Corners still need to be handled separately.

         !NOTE: parallel version implemented here will not give
         !identical results to the serial version since different
         !processes may converge in different number of iterations.
         !This is not a practical issue for the typical precission
         !needed for 3D-RISM.
         do ishell = 1, maxshell

            shell = 0
            ig = 0
            !xy plane
            !select range of sum for this process
            !total number of points in teh outside loop to evaluate
            npt = (this%grid%ngr(2)+2*ishell)
            !get the start point for this process
            n0 = npt/this%grid%mpisize*this%grid%mpirank - ishell
            !the final point is the start point for one process
            !higher.  Ensure we do not exceed the desired range of the
            !sum
            nf = min(npt/this%grid%mpisize*(this%grid%mpirank+1) - ishell-1&
                 , this%grid%ngr(2)-1+ishell)
            do iy = n0,nf
               r(2) = iy*this%grid%grdspc(2)
               do ix = 1-ishell,this%grid%ngr(1)-2+ishell
                  r(1) = ix*this%grid%grdspc(1)
                  !+z
                  r(3) = (this%grid%ngr(3)-1+ishell)*this%grid%grdspc(3)
                  shell = shell +  hlr(this,r)**power
                  !-z
                  r(3) = -ishell*this%grid%grdspc(3)
                  shell = shell +  hlr(this,r)**power
                  ig = ig + 2
               end do
            end do
            !yz plane
            npt = (this%grid%ngr(3)+2*ishell)
            n0 = npt/this%grid%mpisize*this%grid%mpirank - ishell
            nf = min(npt/this%grid%mpisize*(this%grid%mpirank+1) - ishell-1&
                 , this%grid%ngr(3)-1+ishell)
            do iz = n0,nf
               r(3) = iz*this%grid%grdspc(3)
               do iy = 1-ishell,this%grid%ngr(2)-2+ishell
                  r(2) = iy*this%grid%grdspc(2)
                  !+x
                  r(1) = (this%grid%ngr(1)-1+ishell)*this%grid%grdspc(1)
                  shell = shell +  hlr(this,r)**power
                  !-x
                  r(1) = -ishell*this%grid%grdspc(1)
                  shell = shell +  hlr(this,r)**power
                  ig = ig + 2
               end do
            end do
            !zx plane
            npt = (this%grid%ngr(3)+2*ishell-2)
            n0 = npt/this%grid%mpisize*this%grid%mpirank +1- ishell
            nf = min(npt/this%grid%mpisize*(this%grid%mpirank+1) - ishell&
                 , this%grid%ngr(3)-2+ishell)
            do iz = n0,nf
               r(3) = iz*this%grid%grdspc(3)
               do ix = -ishell,this%grid%ngr(1)-1+ishell
                  r(1) = ix*this%grid%grdspc(1)
                  !+y
                  r(2) = (this%grid%ngr(2)-1+ishell)*this%grid%grdspc(2)
                  shell = shell +  hlr(this,r)**power
                  !-y
                  r(2) = -ishell*this%grid%grdspc(2)
                  shell = shell +  hlr(this,r)**power
                  ig = ig + 2
               end do
            end do
            !corners
            if(this%grid%mpirank == 0)then
               do ix = -ishell, this%grid%ngr(1)-1+ishell, this%grid%ngr(1)-1+2*ishell
                  do iy = -ishell, this%grid%ngr(2)-1+ishell, this%grid%ngr(2)-1+2*ishell
                     do iz = -ishell, this%grid%ngr(3)-1+ishell, this%grid%ngr(3)-1+2*ishell
                        r(1) = ix*this%grid%grdspc(1)
                        r(2) = iy*this%grid%grdspc(2)
                        r(3) = iz*this%grid%grdspc(3)
                        shell = shell +  hlr(this,r)**power
                        ig = ig + 1
                     end do
                  end do
               end do
            end if
            total = total+shell
            if(abs(shell/total) < err)then
               exit
            end if
         end do
         if(abs(shell/total) > err)&
              call rism_report_warn("Long range asymptotics integration failed to converge")
         integral = total*(this%solv%charge_sp)**power*this%grid%voxel
      end if
    end function rism3d_potential_int_hn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates integral of the long range asymptotics h for all space
!!!(i.e. not just the grid).
!!!IN:
!!!   this  : the closure object
!!!OUT:
!!!   integral of asymhr for each solvent site
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    function rism3d_potential_int_h(this) result(h)
      use constants, only : PI
      implicit none
      type(rism3d_potential), intent(in) :: this
      _REAL_ :: h(this%solv%natom)
      _REAL_ :: qut

      if(this%grid%mpirank /=0)then
         h=0
         return
      end if

      !we have an analytic expression for this integral so we only need to sum this
      !over the solute and solvent sites
      qut = sum(this%solu%charge)
      h=0d0
      where(this%solv%charge_sp/=0d0) &
           h=-4d0*pi/(this%solv%dielconst*this%solv%xappa**2)&
           *this%solv%charge_sp*qut
    end function rism3d_potential_int_h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Frees all memory and resets values.
!!!IN:
!!!   this   : potential object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism3d_potential_destroy(this)
    implicit none
    type(rism3d_potential), intent(inout) :: this
    call rism_timer_destroy(this%timer)
    call rism_timer_destroy(this%ljTimer)
    call rism_timer_destroy(this%coulombTimer)
    call rism_timer_destroy(this%asympTimer)
    call rism_timer_destroy(this%timer)
    nullify(this%grid)
    nullify(this%solv)
    nullify(this%solu)
    this%cut2 = 0
    if(safemem_dealloc(this%uuv) /=0)&
         call rism_report_error("Uuv deallocation failed")
    if(safemem_dealloc(this%siguv) /=0)&
         call rism_report_error("SIGuv deallocation failed")
    if(safemem_dealloc(this%epsuv) /=0)&
         call rism_report_error("EPSuv deallocation failed")
    if(safemem_dealloc(this%asymcr) /=0)&
         call rism_report_error("asympcr deallocation failed")
    if(safemem_dealloc(this%asymck) /=0)&
         call rism_report_error("asymck deallocation failed")
    if(safemem_dealloc(this%asymhr) /=0)&
         call rism_report_error("asymhr deallocation failed")
    if(safemem_dealloc(this%asymhk) /=0)&
         call rism_report_error("asymhk deallocation failed")
    if(safemem_dealloc(this%huvk0) /=0)&
         call rism_report_error("huvk0 deallocation failed")
    if(safemem_dealloc(this%huvk0_dT) /=0)&
         call rism_report_error("huvk0_dT deallocation failed")
  end subroutine rism3d_potential_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                              PRIVATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Sets up precalculated solute-solvent interactions.  I.e., LJ mixing.
!!!IN:
!!!   this   : potential object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine prep_UV(this)
      use safemem
      implicit none
      type(rism3d_potential), intent(inout) :: this
      integer :: iv, iu

      !........... setting solute-solvent LJ potential parameters ............
      do iv=1,this%solv%natom
         do iu=1,this%solu%natom
            this%siguv(iu,iv) = (this%solu%sig_s(iu)+this%solv%sig_s(iv))
            this%epsuv(iu,iv) = sqrt( this%solu%eps(iu)*this%solv%eps(iv))
         enddo
      enddo
    end subroutine prep_UV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Tabulating the 12-6 Lennard-Jones potential
!!!in the box subject to the minimal image condition
!!!IN:
!!!   this   : potential object
!!!   ulj    : grid to add potential to
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  uljuv (this, ulj)
    use rism_util, only : checksum
    implicit none
#if defined(MPI)
      include 'mpif.h'
#endif /*defined(MPI)*/
      type(rism3d_potential) :: this
    _REAL_, intent(inout) :: ulj(:,:,:,:)
    _REAL_ :: offset
    !number of gridpoints in each direction to use
    !closest point to solute atom
    integer :: grdpnts(3),cp(3),first(3),last(3)
    integer ::  igx,igy,igz, ig, iu, iv, id
    _REAL_ ::  rx,ry,rz, dx2,dy2,dz2, rs2,rs6i,r2,dyz2
    _REAL_ ::  rcor, rcor2
    parameter (rcor=0.2d0, rcor2=rcor**2)

    !..................... enumerating box grid points .....................
    !calculate the number of grid point necessary to cover this cutoff range.
    !In case of a large cutoff, we need to protect against integer overflow
    grdpnts = nint(min(dble(huge(1)-1),sqrt(this%cut2)/this%grid%grdspc))+1
    offset = this%grid%grdspc(3)*(this%grid%nrOff(3))
    do iu=1,this%solu%natom       
       do id=1,3
          cp(id) = nint(this%solu%ratu(id,iu)/this%grid%grdspc(id))
       end do
       cp(3) = cp(3) - this%grid%nrOff(3)
       do id =1,3
          !Note: we have to protect against cp(id)+grdpnts(id) overflowing
          if(cp(id) > 0)then
             first(id) = min(max(1,cp(id)-grdpnts(id)), this%grid%nr(id))
             last(id) = max(min(this%grid%nr(id),&
                  cp(id)+min(huge(1)-cp(id),grdpnts(id))), 1)
          else
             first(id) = min(max(1,cp(id)-min(huge(1)+cp(id),grdpnts(id))), this%grid%nr(id))
             last(id) = max(min(this%grid%nr(id),&
                  cp(id)+grdpnts(id)), 1)
          end if
       end do
       do igz=first(3),last(3)
          rz = (igz-1)*this%grid%grdspc(3)+offset
          dz2 = (rz - this%solu%ratu(3,iu))**2
          do igy=first(2),last(2)
             ry = (igy-1)*this%grid%grdspc(2)
             dy2 = (ry - this%solu%ratu(2,iu))**2
             dyz2 = dz2+dy2
             do igx=first(1),last(1)
                rx = (igx-1)*this%grid%grdspc(1)
                
                !......... site separation subject to minimal image condition ..........
                dx2 = (rx - this%solu%ratu(1,iu))**2
                
                r2 = dx2 + dyz2
                if(r2<this%cut2)then
                   !...................... maintaining potential Ulj ......................
                   do iv=1,this%solv%natom
                      rs2 = r2/this%siguv(iu,iv)**2
                      if (rs2 < rcor2)  rs2 = rcor2
                      rs6i = 1d0/rs2**3
                      ulj(igx,igy,igz,iv) = ulj(igx,igy,igz,iv) &
                           + this%epsuv(iu,iv)*rs6i*(rs6i-2.d0)
                   enddo
                endif
             enddo
          enddo
       enddo
    enddo
  end subroutine uljuv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Synthesizing the Coulomb solute potential in a smaller box.
!!!The exact electrostatic potential is calculated for on a sparse grid (1/8
!!!the number of points or 1/2 density) and for all gridpoints within cutoff
!!!of an atom.  To ensure that all interpolated points are enclosed by eight
!!!neighbours, the sparse grid is extended by one index in each dimension.
!!!IN:
!!!   this   : potential object
!!!   ucu    : grid to add potential to
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine  ucoulu (this,ucu)
    use rism_util, only : checksum
    use safemem
    implicit none
    type(rism3d_potential) :: this
    _REAL_,intent(inout) :: ucu(:,:,:,:)
    _REAL_,pointer ::  smucu(:,:,:)=>NULL()
    integer :: smng3(3)
    _REAL_ :: smratio(3)
    !for MPI spatial decomposition
    !offset :: distance in z-axis due to spatial decomp
    _REAL_ :: offset, smoffset
    !smstart :: initial grid point where the sparse and fine grids have the same value
    integer :: smstart
    integer ::  igx,igy,igz, ig, iu, iv
    _REAL_ ::  r2,ra,rx,ry,rz, dx,dy,dz,dz2,dyz2
    _REAL_ ::  rcor, rcor2
    parameter (rcor=0.002, rcor2=rcor**2)
    integer :: gxstep,gxstart, first(3),last(3)
    integer :: iigx, iigy, iigz
    _REAL_ :: tmp(1)
    !linear spacing of the grid
    _REAL_ :: smgrdspc(3)
    !number of gridpoints in each direction to use
    !closest point to solute atom
    integer :: grdpnts(3),cp(3),id,iu2
    !untouched :: flag indicating that this index of the array has not been set
    _REAL_ :: untouched = HUGE(0d0)
    !............................ screen output ............................
#ifdef RISM_DEBUG
    write(6,*)'tabulating fast solute Coulomb potential ...'
    call flush(6)
#endif /*RISM_DEBUG*/
    smucu=>safemem_realloc(smucu,((this%grid%nr(1)+1)/2+1),((this%grid%nr(2)+1)/2+1),((this%grid%nr(3)+1)/2+1),.false.)
    !......................... clearing summators ..........................
    ucu(:,:,:,1)=untouched
    !
    !Calculate exact potential over the entire sparse grid
    !
    smucu=0
    !ng3+1 takes care of values that are odd
    smng3 = (this%grid%nr+1)/2+1
    smgrdspc=this%grid%grdspc*2d0
    smratio = this%grid%grdspc/smgrdspc
    !calculate the number of grid point necessary to cover this cutoff range.
    !In case of a large cutoff, we need to protect against integer overflow
    grdpnts = nint(min(dble(huge(1)-1),sqrt(this%cut2)/this%grid%grdspc))+1
    offset = this%grid%grdspc(3)*(this%grid%nrOff(3))
    smstart = mod(this%grid%nrOff(3),2)+1
    !if we are start on an odd index (nz_start starts from 0) shift the small grid
    !to a lower offset
    smoffset = offset - (smstart-1)*this%grid%grdspc(3)
    !..................... enumerating box grid points .....................

    do iu=1,this%solu%natom
       do igz=1,smng3(3)
          rz = (igz-1)*smgrdspc(3)+smoffset
          dz = rz - this%solu%ratu(3,iu)
          dz2 = dz*dz
          do igy=1,smng3(2)
             ry = (igy-1)*smgrdspc(2)
             dy = ry - this%solu%ratu(2,iu)
             dyz2=dy*dy+dz2
             do igx=1,smng3(1)
!                ig = 1 + igx + igy*smng3(1) + igz*smng3(2)*smng3(1)
                rx = (igx-1)*smgrdspc(1)
                !...................... summing over solute sites ......................
                !......... getting site separation ..........
                dx = rx - this%solu%ratu(1,iu)

                ra = max(sqrt(dx*dx + dyz2),rcor)
                smucu(igx,igy,igz) = smucu(igx,igy,igz) + this%solu%charge(iu)/ra
            enddo
          enddo
       enddo
    enddo
    !
    !transfer data from the sparse grid to the dense grid
    !
    do igz=smstart,this%grid%nr(3),2
       do igy=1,this%grid%nr(2),2
          do igx=1,this%grid%nr(1),2
             ucu(igx,igy,igz,1) &
                  = smucu((igx+1)/2,(igy+1)/2,smstart -1 + (igz+1)/2)
          enddo
       end do
    end do

    !
    !Calculate exact potential for each grid point within the cutoff
    !
    do iu=1,this%solu%natom
       !determine grid points within cutoff
       do id=1,3
          cp(id) = nint(this%solu%ratu(id,iu)/this%grid%grdspc(id)) !central point (grid point nearest atom)
       end do
       cp(3) = cp(3) - this%grid%nrOff(3)
       do id =1,3
          first(id) = min(max(1,cp(id)-grdpnts(id)),this%grid%nr(id)) !smallest grid index within cutoff
          !Note: we have to protect against cp(id)+grdpnts(id) overflowing
          last(id) = max(min(this%grid%nr(id),&
               cp(id)+min(huge(1)-cp(id),grdpnts(id))), 1) !largest grid index within cuto
       end do
       do igz=first(3),last(3)
          rz = (igz-1)*this%grid%grdspc(3)+offset
          dz = rz - this%solu%ratu(3,iu)
          dz2=dz*dz
          do igy=first(2),last(2)
             ry = (igy-1)*this%grid%grdspc(2)
             dy = ry - this%solu%ratu(2,iu)
             dyz2=dy*dy+dz2
             do igx=first(1),last(1)
                if(ucu(igx,igy,igz,1) == untouched)then
                   rx = (igx-1)*this%grid%grdspc(1)
                   
                   !......... site separation subject to minimal image condition ..........
                   dx = rx - this%solu%ratu(1,iu)
                   r2 = dx*dx + dyz2
                   if(r2<this%cut2) then
                      ucu(igx,igy,igz,1) = 0d0
                      do iu2=1,this%solu%natom
                         dx = rx - this%solu%ratu(1,iu2)
                         dy = ry - this%solu%ratu(2,iu2)
                         dz = rz - this%solu%ratu(3,iu2)
                         
                         ra = max(sqrt(dx*dx + dy*dy + dz*dz),rcor)
                         ucu(igx,igy,igz,1) = ucu(igx,igy,igz,1) + this%solu%charge(iu2)/ra
                      enddo
                   end if
                endif
             end do
          enddo
       enddo
    enddo

    !
    !Interpolate missing values for remaining grid points.
    !
    do igz=1,this%grid%nr(3)
       rz = (igz-1)*this%grid%grdspc(3)+offset
          first(3) = (igz+1)/2
          last(3) = first(3)+1
       do igy=1,this%grid%nr(2)
          ry = (igy-1)*this%grid%grdspc(2)
          first(2) = (igy+1)/2
          last(2) = first(2)+1
          if(mod(igy,2)/=0 .and.mod(igz+smstart-1,2)/=0)then
             gxstep=2
             gxstart=2
          else
             gxstep=1
             gxstart=1
          endif
          do igx=gxstart,this%grid%nr(1),gxstep
             rx = (igx-1)*this%grid%grdspc(1)
             first(1) = (igx+1)/2
             last(1) = first(1)+1
             if(ucu(igx,igy,igz,1) /= untouched)then
                cycle
             endif
             call blend_103(&
                  (dble(igx-1)*smratio(1)-first(1)+1)/(last(1)-first(1)),&
                  (dble(igy-1)*smratio(2)-first(2)+1)/(last(2)-first(2)),&
                  (dble(igz-1+smstart-1)*smratio(3)-first(3)+1)/(last(3)-first(3)),&
                  smucu(first(1),first(2),first(3)),&
                  smucu(first(1),first(2),last(3)),&
                  smucu(first(1),last(2),first(3)),&
                  smucu(first(1),last(2),last(3)),&
                  smucu(last(1),first(2),first(3)),&
                  smucu(last(1),first(2),last(3)),&
                  smucu(last(1),last(2),first(3)),&
                  smucu(last(1),last(2),last(3)),&
                  ucu(igx,igy,igz,1))
          enddo
       end do
    end do

    do iv=this%solv%natom,1,-1
          ucu(:,:,:,iv) =  ucu(:,:,:,1)*this%solv%charge(iv) 
    end do
    ucu=ucu
    if(safemem_dealloc(smucu)/=0)&
         call rism_report_error("Failed to deallocate memory in UCOULU")

    return
  end subroutine ucoulu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the values of the long range asymptotics of h(r) divided by
!!!the species charge at the given point.  A cut off is used to
!!!accelerate the calculation, which should be precomputed by
!!!_setcut().
!!!IN:
!!!   r : position
!!!OUT:
!!!   asymhr at r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function hlr(this,r)
    use constants, only : pi
    implicit none
    type(rism3d_potential),intent(in) :: this
    _REAL_, intent(in) :: r(3)
    _REAL_ :: hlr
    !ra                        : distance
    !half_xappa_smear          : 1/2* kappa   * smear
    _REAL_ :: ra, half_xappa_smear
    integer :: iu
    _REAL_, external :: erfc

    half_xappa_smear = 0.5d0*this%solv%xappa*this%solv%smear
    hlr=0
    if(.not.this%solv%ionic) return
    do iu=1,this%solu%natom
       ra = sqrt(sum((r(:) - this%solu%ratu(:,iu))**2))
       if(ra > this%cut_hlr)then
          !                     !.....h in r-space
          hlr=hlr - this%solu%charge(iu)/ra
       elseif(ra == 0d0)then
          !.....h in r-space
!!$           hlr=hlr - this%solu%charge(iu)*&
!!$               (2d0*exp(-(this%solv%smear*this%solv%xappa/2d0)**2)/(sqrt(PI)*this%solv%smear) &
!!$               - this%solv%xappa * erfc(this%solv%xappa*this%solv%smear/2d0))
          hlr=hlr - this%solu%charge(iu)*&
               (2d0/(sqrt(PI)*this%solv%smear) &
               -exp(half_xappa_smear**2) &
               * this%solv%xappa * erfc(half_xappa_smear))&
               /exp(half_xappa_smear**2)
          !outside the loop we multiply by the inverse of the last line.
          !Since the statement is rarely executed, this should save 
          !a bit of time
       else
          !.....h in r-space
          hlr=hlr - this%solu%charge(iu)/ra &
               *(exp(-this%solv%xappa*ra)*erfc(half_xappa_smear - ra/this%solv%smear) &
               - exp(this%solv%xappa*ra)*erfc(half_xappa_smear+ra/this%solv%smear))/2d0
       end if
    end do
    hlr = hlr*exp(half_xappa_smear**2)/this%solv%dielconst
  end function hlr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the values of the long range asymptotics of c(r) divided by
!!!the species charge at the given point.  A cut off is used to
!!!accelerate the calculation, which should be precomputed by
!!!_setcut().
!!!IN:
!!!   r : position
!!!OUT:
!!!   asymcr at r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function clr(this,r)
    use constants, only : pi
    implicit none
    type(rism3d_potential),intent(in) :: this
    _REAL_, intent(in) :: r(3)
    _REAL_ :: clr
    !ra                        : distance
    _REAL_ :: ra
    integer :: iu
    _REAL_, external :: erfc

    clr=0
    if(.not.this%solu%charged) return
    do iu=1,this%solu%natom
       ra = sqrt(sum((r(:) - this%solu%ratu(:,iu))**2))
       if(ra > this%cut_clr)then
          !.....c in r-space
          clr = clr - this%solu%charge(iu)/ra
       elseif(ra == 0d0)then
          !.....c in r-space
          clr = clr - this%solu%charge(iu)/&
               (sqrt(PI)*this%solv%smear)*2d0
       else
          !.....c in r-space
          clr = clr - this%solu%charge(iu)/ra*(1d0-erfc(ra/this%solv%smear))
       endif
    end do
  end function clr
end module rism3d_potential_c
