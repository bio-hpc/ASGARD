#!/bin/csh

### ###########################################
### ### HYBRID SOLVENT MM-PB/SA CALCULATION ###
### ### Alexander Metz, Holger Gohlke       ### 
### ###########################################
###
### System: Solvated ligand of IL2 (PDB-CODE 1M48)
###
### 1. Prior to the calculation of hybrid solvation free energyies, image the 
### solvent back into the simulation cell and center the solvent to the origin.
### The latter makes defining the center of the solvation sphere trivial 
### (XCAP=0.0; YCAP = 0.0; ZCAP=0.0; see mm_pbsa.in file). Note that all of the 
### solvent molecules should be considered during snapshot creation 
### (as in the 01_SNAPS_SOLV directory); the use of just a sphere or shell of water 
### molecules may cause problems if counterions are able to cross the explicit/implicit
### border between different snapshots.
### 2. Calculation of gasphase energies is done in an additional step. Snapshots 
### without solvent molecules need to be generated for this (as in the 
### 02_SNAPS_DESOLV directory). 
### 3. Parameter and topology files for the solvated and nonsolvated systems are in the
### directories 03_PRM_S and 04_PRM_DES, resp.
### 4. Calculation of the hybrid solvation free energy and the gasphase energy is
### carried out in the directories 05_HYB_SOLV_ENE and 06_GASPHASE_ENE, resp.
###

#  setenv AMBERHOME path_to_Amber12_directory

  cd $AMBERHOME/src/mm_pbsa/Examples/07_MMPBSA_Hyb/05_HYB_SOLV_ENE/
  echo "CALCULATING HYBRID SOLVATION ENERGY"
  $AMBERHOME/bin/mm_pbsa.pl mm_pbsa.in >&! mm_pbsa.out

  cd $AMBERHOME/src/mm_pbsa/Examples/07_MMPBSA_Hyb/06_GASPHASE_ENE/
  echo "CALCULATING INTERNAL ENERGY"
  $AMBERHOME/bin/mm_pbsa.pl mm_pbsa.in >&! mm_pbsa.out


### 5. From these two calculations, effective energies are obtained by adding the 
### corresponding values in the columns of the ~_statistics.out files manually.
### That is, add the PBTOT energy of 05_HYB_SOLV_ENE/~_statistics.out to the GAS
### energy of 06_GASPHASE_ENE/~_statistics.out.
###
### Below are some explanations of the output of the hybrid solvation free energy calculation:
###
### EDISPER: Part of the nonpolar solvation free energy due to attractive dispersion interactions,
###          calculated by summing van der Waals interactions between the solute and the solvent 
###          molecules. Interactions between solvent molecules or within the solute are not considered.
### ELRAELE: Electrostatic part of the solvation free energy, calculated using a linear response 
###          approximation taking into account the explicit solvent.
###          Note that the output in the ~.all.out files is in units of kT (see below).
### EPB:     Electrostatic part of the solvation free energy originating from the reaction field
###          outside the region of explicit solvent. This contribution is calculated using the
###          Poisson-Boltzmann approach.
###          Note that the output in the ~.all.out files is in units of kT (see below).
### corrected reaction field energy: Sum of ELRAELE and EPB
### surface area: The solvent-excluded surface area in [A^2].
### ECAVITY: The solvent-excluded surface area in [A^2].
### PBSUR  : Part of the nonpolar solvation free energy due to the formation a cavity for the solute.
###          Taken to be proportional to the solvent-excluded surface area, with the factor
###          gamma = 0.069 [kcal/mol*A^2] and the constant beta = 0.0 [kcal/mol].
### PBCAL  : PBCAL = PBELE + GAS
###          As no electrostatic gasphase contribution is calculated during the hybrid solvation step,
###          PBCAL = PBELE.
### PBSOL  : PBSOL = EDISPER + PBSUR + ELRAELE (converted to [kcal/mol]) + EPB (converted to [kcal/mol]); 
###          the total solvation free energy.
### PBELE  : PBELE = ELRAELE (converted to [kcal/mol]) + EPB (converted to [kcal/mol]);
###          total electrostatic part of the solvation free energy.
### PBTOT  : PBTOT = PBSOL
###          To get to the effective energy for snapshots, the gasphase energy has to be added to this number.
###
### The ELRAELE, EPB, and "corrected reaction field energy" in the ~.all.out files are given in units of kT.
### The conversion factor from kT units to kcal/mol is: 
### (300.0 [K] * 8.314 [J/mol*K]) / 1000.0 * 4.184 [cal] = 0.5961281070745697897 kcal/kT

exit
