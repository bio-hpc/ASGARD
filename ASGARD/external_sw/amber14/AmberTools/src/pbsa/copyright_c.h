/*******************************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! PBSA
!
! An analysis program for solvent-mediated electrostatic and nonelectrostatic
! interactions in biomolecules.
!
! Please acknowledge your use of PBSA by citing:
!
! Luo, David, and Gilson, J. Comp. Chem. 23:1244-1253, 2002.
!
! Major Developers:
!
! Ray Luo, corresponding author
! Jun Wang, linear Poisson-Boltzmann numerical solvers
! Qin Cai, nonlinear Poisson-Boltzmann numerical solvers
! Xiang Ye, electrostatic energy and force numerical algorithms
! Meng-Juei Hsieh, program interface and parallel implementations
! Chuck Tan, nonelectrostatic energy and force numerical algorithms
!
! Additional contributing authors are listed in the code documentation.
!
! Departments of Molecular Biology and Biochemistry and Biomedical Engineering,
! University of California, Irvine, California
!
! Copyright (c) 2004-2009. The Regents of the University of California.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!
! Portion of the routine "pb_exmol" is modified from the UHBD program
! (Comp. Phys. Comm. 91:57-95, 1995), copyrighted by University of Houston,
! 1989-2009. See program documentation for more information.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Overview:
!
! PBSA models the electrostatic solvation interaction by the Poisson-Boltzmann
! equation. The implementation (module poisson_boltzmann) uses the finite-
! difference method to solve the partial differential equation. Both linear
! and full nonlinear numerical solvers are implemented. Please refer to the
! following publications:
!
! 1. Luo, David, and Gilson, J. Comp. Chem. 23:1244-1253, 2002.
! 2. Qin, Wang, Zhao, and Luo, J. Chem. Phys. 130:14501, 2009.
! 3. Wang and Luo, J. Comp. Chem. In Press, 2010.
!
! for implementation details of the linear solvers. The full nonlinear solvers
! are documented in:
!
! 1. Cai, Wang, and Luo, J. Chem. Comp. Theo. In Press, 2009.
!
! The electrostatic energy and forces are computed based on the finite-
! difference grid potentials as discussed in the following publications:
!
! 1. Lu and Luo, J. Chem. Phys. 119:11035-11047, 2003.
! 2. Xiang, Wang, and Luo, J. Chem. Phys. Submited. 2009.
!
! The dielectric models and molecular surface used in the electrostatic 
! solvation model are documented in:
!
! 1. Wang, Ye, and Luo, J. Phys. Chem. Submitted, 2009.
! 2. Ye and Luo, J. Phys. Chem. Submitted, 2009.
!
! The parameters used in the electrostatic solvation model are documented in
! the following publication:
!
! 1. Tan, Yang and Luo, J. Phys. Chem. 110:18680-18687, 2006.
!
! PBSA models the nonelectrostatic solvation interaction by two separate
! terms, dispersion (or van der Waals or attractive) and cavity (or
! hydrophobic or repulsive), in this implementation (module dispersion_cavity).
! The dispersion term is computed by a numerical integration over the solvent
! accessible surface area. The cavity term is modeled by a term proportional to
! the molecular surface or volume with a single proportional constant. See the
! following publication for details:
! 
! Tan, Tan and Luo, J. Phys. Chem. 111:12263-12274, 2007.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file is part of PBSA.
!
! PBSA is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! PBSA is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PBSA; if not, write to the Free Software Foundation, Inc., 59 Temple Place,
! Suite 330, Boston, MA 02111-1307, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*******************************************************************************/
