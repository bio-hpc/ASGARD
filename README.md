## ASGARD: Automated Scripts for Gromacs to Analyse and Run molecular Dynamics

ASGARD is a tool that allows via an automated MD workflow to analyse performed MD protein or protein-ligand complex simulations and to easily generate the corresponding report 

### Installation
  1. git clone https://github.com/bio-hpc/ASGARD.git
  2. git clone git@github.com:bio-hpc/ASGARD.git
  3. gh repo clone bio-hpc/ASGARD 
  4. Download the .zip and unzip it in the supercomputing centers you are going to use

### Download singularity image 
Needed to secure compatibility with all cluster.

wget "https://drive.google.com/uc?export=download&id=1Q9ifMDEaoxhsI9eealvsf_cTvKCYkJ8H&confirm=t" -O singularity/ASGARD.simg

### Available Analyses

**1. TARGET**: for Protein MD simulations
  - Root Mean Square Desviation (RMSD)  [ All molecules, protein and ligand ] <br />
  - RMSD Fluctuation  [ For steps, protein and ligand ] <br />
  - Radius of gyration  [ Total and around axes ] <br />
  - Distance center of mass  [ Protein and ligand ] <br />
  - Solvent Accesible Surface (SASA)  [ Over time ] <br /> 
  - Kabsch and Sander dictionary of protein secundary structure (DSSP)  [ Num. aminoacid in each ss and evolution ] <br />

**2. TARGET_QUERY**: for Protein-Ligand Complex MD simulations
  - Same as above (also calculated for the ligand), and Protein-Ligand Interactions analyses: <br />
     - MM-PBSA, interaction energy values, number HBonds during simulation and interactions between ligand and protein in last frame  <br />

### Needed packages
1. **Singularity**>=3.6.0
### Aditional Software (in external_sw folder)
1. **Poseview [Required ChemAxon]**:  generates publication-quality 2D structure-diagrams of protein-ligand complexes.
2. **GROMACS**: to use g_mmpbsa (found in ASGARD/analyze_trajectory/extra/) and different type of force fields, and to create topology with topology/generate_topology.
   
### Usage instructions

Firstly, you need a folder with all outputs files from GROMACS Molecular Dynamics 

Parameters needed to launch ASGARD are the folder with GROMACS results files (-f) (xtc, tpr, gro, mpd used, files and, pdb target and mol2 query files) and profile which you want (-p) , TARGET (a unique protein simulations) and TARGET_QUERY (protein-ligand simultions)

For instance, you can launch ASGARD this way:

sh ASGARD.sh -f 6b0t_6b0t_lig_MD_100ns/ -p TARGET_QUERY

If there are problems with launching ASGARD.sh:

singularity exec --bind ${PWD} singularity/ASGARD.simg ./ASGARD.sh -f 6b0t_6b0t_lig_MD_100ns/ -p TARGET_QUERY

