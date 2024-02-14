#!/bin/bash
name_shuttlemol="./ASGARD.sh"
path_analize_results="ASGARD/analyze_results/"
path_external_sw="ASGARD/external_sw/"
printHelpAux(){
	a=`echo $1 | sed 's/\\\t/\ /g'`
	printf "${GREEN} %8s ${CYAN}%-25s${NONE} %-30s \n" "$a" "$2" "$3"
}
printGlobalHelp(){
	printf "${GREEN} %6s ${CYAN}%5s${CYAN}%s ${NONE}%-s \n" "-$1 )" "[ $2 ]" "$3" "$4"
}

printFiles(){
	printf "${GREEN} %6s ${CYAN}%5s${CYAN}%s ${NONE}%-s \n" "$1" "[ $2 ]" "$3" "$4"
}


printAnalysis(){
	printf "${GREEN} %s ${NONE}%s${NONE}%s ${NONE}%-s \n" "-$1 " "[ $2 ]"
}

echo ""
option=`echo "$1" | tr '[:lower:]' '[:upper:]'`

cat ${pathSL}logo_ASGARD.txt

echo ""

# ASGARD analysis

  echo -e "${PURPLE}Welcome to ASGARD (Automated Scripts for Gromacs to Analyze and Run Dynamics)${NONE} "
  echo " "
  echo -e "${CYAN}ASGARD is a group of python scripts to analyze Molecular Dynamics trajectories obtained by GROMACS${NONE} "
  printAnalysis "Root Mean Square Desviation (RMSD)" "All molecules, protein and ligand"
  printAnalysis "RMSD Fluctuation" "For steps, protein and ligand"
  printAnalysis "Radius of gyration" "Total and around axes"
  printAnalysis "Distance center of mass" "Protein and ligand"
  printAnalysis "Solvent Accesible Surface (SASA)" "Over time"
  printAnalysis "Kabsch and Sander dictionary of protein secundary structure (DSSP)" "Num of aminoacid in each ss and evolution" 
  printAnalysis "Protein-Ligand Interactions (HBonds, MM/PBSA, scores for each binding site residue)" 

# inputs 

if [  -z "$option" ];then
  echo ""
  printf "\t${CYAN}%s${NONE} \n"  "[Input options]"
  echo -e "${CYAN}____________________________________________________${CYAN} "
  echo ""
  echo -e "${PURPLE} Input directory${NONE} "
  echo -e "${PURPLE}____________________________________________________${NONE} "
  echo ""
  printGlobalHelp "f" "folder" "" "Absolute path where MD results files are found"
  echo -e "${PURPLE}It must contain the next files inside them${NONE} "
  printFiles "" ".top" "" "Topology file"
	printFiles "" ".xtc" "" "Trayectory file in xtc format"
  printFiles "" ".trr" "" "Trayectory file in trr format"
	printFiles "" ".gro" "" ".gro topology file "
  printFiles "" ".npt_1.gro" "" ".gro topology file used for npt equilibration step "
  printFiles "" ".tpr" "" ".tpr topology file" 
  printFiles "" ".mdp" "" ".mdp configuration file "
  printFiles "" ".mol2" "" "ligand structure mol2 file "
  printFiles "" ".pdb" "" "target structure pdb file"
  printFiles "" ".edr" "" "Energies file "
    echo ""
  echo ""
  echo -e "${PURPLE} Profiles${NONE} "
  echo -e "${PURPLE}____________________________________________________${NONE} "
  echo ""
  printGlobalHelp "p" "str" "" "Analysis Profile"
  echo -e "${PURPLE}This tool has different profiles depending on type of Molecular Dynamics simulation. ${NONE} "
  echo -e "${PURPLE}Select one of the following profile for this option. ${NONE} "
  printFiles "" "TARGET_QUERY" "" "Analysis for Protein-Ligand Simulation"
	printFiles "" "TARGET" "" "Analysis for Protein Simulation"
	printFiles "" "TEST" "" "Testing"
  
  echo ""
  printf "\t${CYAN}%s${NONE} \n"  "[Optional parameters]"
  echo -e "${CYAN}____________________________________________________${CYAN} "
  echo ""
  echo -e "${PURPLE}In case you need to modify the analysis configuration${NONE} "
  printGlobalHelp "ref" "reference_structure" "" "Choose the structure taken as a reference in RMSD and RMSF calculations (advisable to place it within the Input directory)"
  printGlobalHelp "lig" "ligand group" "" "Indicate the name of the ligand group or that group of atoms you want to take as a ligand "
  echo ""
  echo ""
  printf "\t${CYAN}%s${NONE} \n"  "[Output]"
  echo -e "${CYAN}____________________________________________________${CYAN} "
  echo ""
  echo -e "${PURPLE} Output files${NONE} "
  echo -e "${PURPLE}____________________________________________________${NONE} "
  echo ""
  printFiles "" "folder" "" "Folder with analysis results files and images, and files used by ASGARD"
	printFiles "" ".pdf" "" "PDF Report generated with the generated files (in the folder)"
  

	echo ""
else


option=`echo "$1" | tr '[:lower:]' '[:upper:]'`

exit
fi
