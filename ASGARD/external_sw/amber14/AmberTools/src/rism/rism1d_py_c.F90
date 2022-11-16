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
!!!Percus-Yevick closure class for 1D-RISM.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module rism1d_py_c
  use safemem
  !the PY type is actually empty
  type rism1d_py
     logical py
  end type rism1d_py

  public rism1d_py_new, rism1d_py_destroy, rism1d_py_gvv, rism1d_py_exchem,&
       rism1d_py_exchemIon, rism1d_py_pressureFE, rism1d_py_freeEnergy
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Initializes the PY closure
!!!IN:
!!!   this : PY object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_py_new(this)
    implicit none
    type(rism1d_py), intent(inout) :: this
  end subroutine rism1d_py_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Gvv from Uvv, Hvv, and Cvv using the PY closure
!!!IN:
!!!   this : the PY closure object
!!!   gvv  : site-site pair correlation function
!!!   uvv  : site-site potential
!!!   hvv  : site-site total correlation function
!!!   cvv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_py_gvv(this,gvv, uvv, hvv, cvv)
    implicit none
    type(rism1d_py), intent(in) :: this
    _REAL_, intent(out) :: gvv(:,:)
    _REAL_, intent(in) :: uvv(:,:),hvv(:,:),cvv(:,:)
    integer :: ivv, ir
    _REAL_ :: tvv
    do ivv = 1, ubound(gvv,2)
       do ir = 1, ubound(gvv,1)
          tvv = (hvv(ir,ivv) - cvv(ir,ivv))
          gvv(ir,ivv) = exp( - uvv(ir,ivv)) * (1.d0 + tvv)
       end do
    end do
  end subroutine rism1d_py_gvv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Bvv (bridge function) from Uvv, Gvv, and Cvv using the PY closure
!!!IN:
!!!   this : the PY closure object
!!!   uvv  : site-site potential
!!!   gvv  : site-site total correlation function
!!!   cvv  : site-site direct correlation function
!!!OUT:
!!!   bvv  : (pointer) site-site bridge function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_py_bvv(this, uvv, gvv, cvv) result(bvv)
    implicit none
    type(rism1d_py), intent(in) :: this
    _REAL_, pointer :: bvv(:,:)
    _REAL_, intent(in) :: uvv(:,:),gvv(:,:),cvv(:,:)
    integer :: ivv, ir
    _REAL_ :: tvv
    nullify(bvv)
    bvv => safemem_realloc(bvv, ubound(gvv,1), ubound(gvv,2))
    do ivv = 1, ubound(bvv,2)
       do ir = 1, ubound(bvv,1)
          tvv = (gvv(ir,ivv) -1d0 - cvv(ir,ivv))
          bvv(ir,ivv) = -tvv + log(1.d0 + tvv)
       end do
    end do
  end function rism1d_py_bvv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the excess chemical potential in kT
!!!IN:
!!!   this : the closure object
!!!   gvv  : site-site pair correlation function
!!!   cvv  : site-site direct correlation function
!!!   mtv  : multiplicity of each site
!!!   rhov : number density of each site
!!!   rhotrgt : (optional) The final, physical density for thermodynamics.  This
!!!             can be used as an effective correction for some closures that
!!!             over estimate the pressure (e.g. PY and PY)
!!!   dr   : r-space grid spacing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_py_exchem(this, gvv, cvv, mtv, jvv, rhov, rhotrgtv, dr) result(exchem)
    use constants, only : pi
    implicit none
    type(rism1d_py), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),cvv(:,:), rhov(:), rhotrgtv(:), dr
    integer, intent(in) :: mtv(:),jvv(:,:)
    _REAL_ :: exchem(ubound(mtv,1))
    _REAL_ :: r, h2c,hvv, exchemv
    integer :: ivv, iv1, iv2, ir
    
    !No analytic expression for chemical potential so return NaN
    exchem = 0d0
    exchem = exchem/exchem
end function rism1d_py_exchem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the ionic(?) excess chemical potential in kT.  This seems to take
!!!the contribution only from the last solvent type.  I'm not sure what the point is
!!!here.
!!!IN:
!!!   this : the closure object
!!!   gvv  : site-site pair correlation function
!!!   cvv  : site-site direct correlation function
!!!   mtv  : multiplicity of each site
!!!   rhov : number density of each site
!!!   rhotrgt : (optional) The final, physical density for thermodynamics.  This
!!!             can be used as an effective correction for some closures that
!!!             over estimate the pressure (e.g. PY and PY)
!!!   dr   : r-space grid spacing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_py_exchemIon(this, gvv, cvv, mtv, jvv, rhov, rhotrgtv, dr) result(exchem)
    use constants, only : pi
    implicit none
    type(rism1d_py), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),cvv(:,:), rhov(:), rhotrgtv(:), dr
    integer, intent(in) :: mtv(:),jvv(:,:)
    _REAL_ :: exchem(ubound(mtv,1))
    _REAL_ :: r, h2c,hvv, exchemv
    integer :: ivv, iv1, iv2, ir
    
    !No analytic expression for chemical potential so return NaN
    exchem = 0d0
    exchem = exchem/exchem
  end function rism1d_py_exchemIon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the r-space contribution to the pressure in internal units
!!!using the free energy route
!!!IN:
!!!   this : the closure object
!!!   gvv  : site-site pair correlation function
!!!   cvv  : site-site direct correlation function
!!!   mtv  : multiplicity of each site
!!!   rhov : number density of each site
!!!   dr   : r-space grid spacing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_py_pressureFE(this, gvv, cvv, mtv, rhov, dr) result(pr)
    use constants, only : pi
    implicit none
    type(rism1d_py), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),cvv(:,:),rhov(:),dr
    integer, intent(in) :: mtv(:)
    _REAL_ :: pr
    _REAL_ :: prv
    _REAL_ :: r, hvv, h2c
    integer :: ir,ivv,iv1,iv2, cnt
    !r-space
    !No analytic expression for pressure along free energy route so return NaN
    pr = 0d0
    pr = pr/pr
  end function rism1d_py_pressureFE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the r-space contribution to the free energy in internal units
!!!IN:
!!!   this : the closure object
!!!   gvv  : site-site pair correlation function
!!!   cvv  : site-site direct correlation function
!!!   mtv  : multiplicity of each site
!!!   rhov : number density of each site
!!!   dr   : r-space grid spacing
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_py_freeEnergy(this, gvv, cvv, mtv, rhov, dr) result(fe)
    use constants, only : pi
    implicit none
    type(rism1d_py), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:),cvv(:,:),rhov(:),dr
    integer, intent(in) :: mtv(:)
    _REAL_ :: fe
    _REAL_ :: fev
    _REAL_ :: r, hvv, h2c
    integer :: ir,ivv,iv1,iv2, cnt
    !r-space
    !No analytic expression for free energy so return NaN
    fe = 0d0
    fe = fe/fe
  end function rism1d_py_freeEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Frees memory and resets the PY closure
!!!IN:
!!!   this : PY object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_py_destroy(this)
    implicit none
    type(rism1d_py), intent(inout) :: this
  end subroutine rism1d_py_destroy
end module rism1d_py_c
