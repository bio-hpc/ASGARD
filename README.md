## ASGARD: Automated Scripts for Gromacs to Analyse and Run molecular Dynamics

[![DOI](https://zenodo.org/badge/doi/10.26434/chemrxiv-2023-3m9mk.svg?style=svg)](https://doi.org/10.26434/chemrxiv-2023-3m9mk)


ASGARD is a tool that allows via an automated MD workflow to analyse performed MD protein or protein-ligand complex simulations and to easily generate the corresponding report 

### Installation

Using git clone
```bash
  git clone https://github.com/bio-hpc/ASGARD.git
```
```bash
  git clone git@github.com:bio-hpc/ASGARD.git
```

Download the .zip and unzip it in the supercomputing centers you are going to use


### Download singularity image 
The singularity image generated for ASGARD and needed to secure compatibility with any cluster.

```bash
wget "https://drive.google.com/uc?export=download&id=1Q9ifMDEaoxhsI9eealvsf_cTvKCYkJ8H&confirm=t" -O singularity/ASGARD.simg
```

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
     - MM-PBSA <br />
     - Interaction energy values <br />
     - HBonds contacts during simulation <br />
     - 2D and 3D Interactions diagrams in last frame  <br />

**3. TARGET_QUERY_NO_HB**: for Protein-Ligand Complex MD simulations
  - Skip Hydrogen bonds analysis for a faster analysis

### Needed packages
The only thing that you need to run ASGARD is the Singularity package installed in your cluster
1. **Singularity**>=3.6.0

### Aditional Software (in external_sw folder)
1. **Poseview [Required License]**:  generates publication-quality 2D structure-diagrams of protein-ligand complexes. To access your free evaluate license, please enter to https://www.biosolveit.de/license/evaluation/ 
2. **GROMACS**: to use g_mmpbsa (found in ASGARD/analyze_trajectory/extra/) and different type of force fields, and to create topology with topology/generate_topology.

### Usage instructions

Firstly, you need a folder with all outputs files from GROMACS Molecular Dynamics 

Parameters needed to launch ASGARD are the folder with GROMACS results files (-f) (xtc, tpr, gro, mpd used) and profile which you want (-p): TARGET (a single protein simulations) and TARGET_QUERY (protein-ligand simultions), and also the group name of the ligand that you want to analyze (-l)

For instance, you can launch ASGARD this way:
```bash
./ASGARD.sh -f StudyCase2/ -p TARGET_QUERY
```
