!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Implicit Solvation
!
! Implicit Solvation contains numerical implicit solvation
! treatments of biomolecules for both electrostatic and
! nonelectrostatic interactions.
!
! The electrostatic interaction is modeled by the linearized
! Poisson-Boltzmann equation. The implementation (module
! poisson_boltzmann) uses the finite difference method and the
! particle-particle particle-mesh strategy to compute the total
! electrostatic energy and forces for molecular mechanics and
! dynamics simulations. See J. Chem. Phys. 119:11035-11047, 2003
! and J. Comp. Chem. 23:1244-1253, 2002 for implementation details.
! See Tan, Yang and Luo, In Preparation for the parameters used.
! Portions of the code are modified from the UHBD program (Comp.
! Phys. Comm. 91:57-95, 1995). See code documentation for more
! information.
!
! The nonelectrostatic interaction is modeled by two separate
! terms: dispersion and cavity (hydrophobic), in this implementation
! (module dispersion_cavity). The dispersion term is computed by a
! numerical integration over the solvent accessible surface area.
! The cavity term is modeled by a term proportional to the solvent
! accessible surface area and with a single proportional constant.
! See Tan and Luo, In Preparation for implementation detail and
! parameters used.
!
! Both implementations requires a numerical solvent accessible surface
! area calculation, as implemented in module solvent_accessibility. See
! Luo, In Preparation for implementation detail.
!
! Ray Luo (rluo@uci.edu)
! Department of Molecular Biology and Biochemistry
! University of California, Irvine
!
! Additional contributing authors are listed in the code
! documentation.
!
! Copyright (c) 2004-2006. The Regents of the University of
! California. Portions Copyright (c) 1989-2006. University of
! Houston.
!
! This file is part of Implicit Solvation.
!
! Implicit Solvation is free software; you can redistribute it
! and/or modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2 of
! the License, or (at your option) any later version.
!
! Implicit Solvation is distributed in the hope that it will be
! useful, but WITHOUT ANY WARRANTY; without even the implied
! warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Implicit Solvation; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston,
! MA 02111-1307, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
