## ASGARD: Automated Scripts for Gromacs to Analyse and Run molecular Dynamics

[![DOI](https://zenodo.org/badge/doi/10.1080/07391102.2024.2349527.svg?style=svg)](https://doi.org/10.1080/07391102.2024.2349527)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10213139.svg)](https://doi.org/10.5281/zenodo.10213139)

ASGARD is a tool that automates the analysis of molecular dynamics (MD) simulations of proteins or protein-ligand complexes, and generates a corresponding summary report through a streamlined workflow.

### Installation

**1.** Using git clone
```bash
  git clone https://github.com/bio-hpc/ASGARD.git
```
```bash
  git clone git@github.com:bio-hpc/ASGARD.git
```

**2.** Download the .zip and unzip it in the High performance computing (HPC) cluster you are going to use it


### Download singularity image 
The ASGARD.simg singularity image is needed to secure compatibility with any cluster that supports Singularity.

```bash
wget --no-check-certificate 'https://drive.usercontent.google.com/download?id=1wHvmtUVhUz9DAzPnqeZPE7MVX5fRERGb&export=download&authuser=1&confirm=t&uuid=0c83343d-17fe-4282-bf07-9a2321537a9a&at=APZUnTW_78yhd6klINcZBOjxIU6g:1706872870521' -O singularity/ASGARD.simg
```

### Available Analyses

**1. TARGET**: for Protein MD simulations
  - Root Mean Square Deviation (RMSD)  [ All molecules, protein and ligand ] <br />
  - RMSD Fluctuation  [ For steps, protein and ligand ] <br />
  - Radius of gyration  [ Total and around axes ] <br />
  - Distance center of mass  [ Protein and ligand ] <br />
  - Solvent Accessible Surface Area (SASA)  [ Over time ] <br /> 
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
1. **Poseview [Required License]**:  generates publication-quality 2D structure-diagrams of protein-ligand complexes. A free evaluation license can be requested at https://www.biosolveit.de/license/evaluation/ 
2. **GROMACS**: to use g_mmpbsa (found in ASGARD/analyze_trajectory/extra/) and different type of force fields, and to create topology with topology/generate_topology.

### Usage instructions

To use ASGARD, you will need a folder containing all the output files from your GROMACS molecular dynamics simulations.

The required parameters to launch ASGARD are:
- The folder containing the GROMACS results files (-f), which should include files with the following extensions: xtc, tpr, gro, and mpd.
- The profile you want to use (-p), which can be either TARGET (for single protein simulations) or TARGET_QUERY (for protein-ligand simulations).

For instance, you can launch ASGARD this way:
```bash
./ASGARD.sh -f StudyCase2/ -p TARGET_QUERY
```

Optional parameters are available to modify the default configuration:
  - Ligand group (-lig) [ Name of the group for the atoms you want to designate as the ligand ] <br />
  - Reference structure (-ref) [ Structure used for the RMSD and RMSF calculations (.tpr or .gro formats)  ] <br />


We recommend launching some of the case studies to test the tool and understand its functionality. The case studies can be found at the following Zenodo link:
https://doi.org/10.5281/zenodo.10213139
Additionally, you can read the article about this tool, published in the Journal of Biomolecular Structure and Dynamics: "Enhancing MD simulations: ASGARD's automated analysis for GROMACS"
https://doi.org/10.1080/07391102.2024.2349527
