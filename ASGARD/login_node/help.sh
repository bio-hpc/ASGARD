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


#printHelpSwexclusivo()
#{
#
#		if [ -z $software ];then
#			echo ""
#			printHelpAux "Opciones extras"
#			path_login_node=$pathSL
#			source ${pathSL}special_params.sh
#			read_template_params $option
#			for i in `seq 1 $contadorFile`;do
#				IFS='::' read -ra ADDR <<< "${file[$i]}"
#				printHelpAux "(${ADDR[0]})\t" "${ADDR[4]}" "${ADDR[6]}"
#			done
#		fi
#
#}
echo ""
option=`echo "$1" | tr '[:lower:]' '[:upper:]'`

cat ${pathSL}logo_ASGARD.txt

echo ""

#analisis de ASGARD

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
  #printf "\t${CYAN}%s${NONE} \n"  "[Profiles]"
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
  echo -e "${PURPLE} Output files${NONE} "
  echo -e "${PURPLE}____________________________________________________${NONE} "
  echo ""
#  printGlobalHelp "p" "str" "" "Analysis Profile"
#  echo -e "${PURPLE}This tool has different profiles depending on type of Molecular Dynamics simulation. Select one of the following profile for this option: ${NONE} "
  printFiles "" "folder" "" "Folder with analysis results files and images, and files used by ASGARD"
	printFiles "" ".pdf" "" "PDF Report generated with the generated files (in the folder)"
#	printGlobalHelp "o" "N" "" "Technique `readTechniquesSW ${pathSL}techniques/ SLTechnique`"
#	printGlobalHelp "s" "N" "" "Available Software `readTechniquesSW ${pathSL}templateParams/ template`"
#	printGlobalHelp "se" "O" " [ N | Y ]" "Secuencial Mode Y, queue manager N. ( Default N )."
#	printGlobalHelp "de" "O" " [ 1-10  ]" "Debug mode range[1?-?10?]. It can only be used in sequencuial mode (def?a?ult 0)."
#	printGlobalHelp "el" "O" " [ N | Y ]" "Insert data in excell bio-hpc sheet. (Defaul y)."
#	printGlobalHelp "prp" "N" "" "Protocol used to prepare recptor or query."
#	printGlobalHelp "prl" "N" "" "Protocol used to prepare ligand/s."

	echo ""
else


option=`echo "$1" | tr '[:lower:]' '[:upper:]'`

#case $option in
#	#
#	#	explicar profiles
#	#
#	GR)
#		printOp
#		
#		
#		printExample "${name_shuttlemol} -p targets/3s5z.pdb -l queries/testDm2/ -o VS -s GR -j 1  [-GP]"
#
#		echo -e "\tError Comun cuando se modifica/genera la proteina con maestro:"
#		echo -e "\tmodificar todos los NMA por NME y cambiar"
#		echo -e "\t\tATOM   5228  CA  NMA A1020A    -72.002   7.253 -39.562  1.00  0.00           C  "
#		echo -e "\tpor"
#		echo -e "\t\tATOM   5228  CH3 NME A1020A    -72.002   7.253 -39.562  1.00  0.00           C  "
#		echo ""
#		echo -e "${GREEN}Ojo la proteina debe ser pdb y el ligando mols"
#
#
#		echo -e "${GREEN}Ojo para BIPHSIC_SYSTEMS el primier ligando que se duplicara como debe tener el nombre menor ordenado alfabeticamente en la carpeta"
#
#		echo -e "python ${path_external_sw}/gromacs/topology/generate_topology.py -t projects/nlrp3/nlrp3.pdb -q projects/nlrp3/lig/ -p TARGET_QUERY -g gmx_mpi ${NONE}"
#
#		echo ""
#;;
#	GRN)
#		printOp
#		printHelpAux "-d" "[    ]" "Carpeta de la prueba realizada previamente"
#		printExample "${name_shuttlemol} -p protein -l folderLiagnd -o VS -s GRN -j X -prl X -prp X -stepSimulacionDM X [-GP [1-2]] \nExample: ${name_shuttlemol} -t projects/testgrn/1le0.pdb  -o VS -s GRN -prp na -prl na  -el n -j 1 -se y -de 100 -step_md 100  -d VS_GR_1le0_target_2019-03-18 -gp 2 -co 24";;
#	GRC)
#		printOp
#		printHelpAux "-d" "[    ]" "Carpeta de la prueba realizada previamente"
#
#		printExample "${name_shuttlemol} -p protein -l folderLiagnd -o VS -s GRN -j X -prl X -prp X [-GP [1-2]] \nExample: ./sm.sh -p proyectos/3s5z/3s5z_b.pdb -l proyectos/3s5z/queries/ -o VS -s GRC -el n -j 1 -prp gromacs -prl gauss  -d test3s5z2/  -se y -de 10";;
#	*)
#		printHelpAux "No se encuentra ayuda"
#esac
exit
fi

#
#printTitle()
#{
#	title=`cat ${pathSL}templateParams/template${option}.txt |grep "#LanzNombre" | cut -d : -f 2`
#	echo -e "${BLUE}\t\t\t\t\t ${title} ${NONE}"
#	echo -e "${CYAN}-------------------------------------------------------------------------------------------------------${NONE}"
#
#}
#printTitletechnique()
#{
#	echo -e "${BLUE}\t\t\t\t\t $1 ${NONE}"
#	echo -e "${BLUE}\t\t\t\t\t ${title} ${NONE}"
#	echo -e "${CYAN}-------------------------------------------------------------------------------------------------------${NONE}"
#
#}
#printExampleTechnique()
#{
#	echo""
#	echo -e "${YELLOW}Example: $1${NONE}"
#	echo ""
#	echo ""
#}
##
##	Muestra en la ayuda de las tecnicas los softwares utilizables
##
#printOpTchnique()
#{
#	cad="["
#	for i in `ls ${pathSL}templateParams/*`; do
#		ext=`cat $i |grep "#lanzOptions" | awk -F: '{print $2}'`
#		if [[ $ext == *"$option"* ]]; then
#			aux=`cat $i |grep "#LanzNombre" | awk -F: '{print $2}' |awk -F- '{print $1}'`
#				cad="${cad} ${aux}|"
#		fi
#
#	done
#
#	cad=${cad%?}
#
#	cad="$cad]"
#	printf "${GREEN} %8s ${CYAN}%-135.135s${NONE}\n" "-s|-S" "$cad"
#
#}
#printManual()
#{
#	echo -e "${GREEN}Manual: ${BLUE}$1${NONE}"
#}
#printExample()
#{
#
#	printHelpAux "-j" "[    ]" "numero de jobs que se dividira la prueba"
#	printHelpSwexclusivo
#	echo ""
#	printHelpAux "Opciones por defecto"
#
#	extensionesProt=`cat ${pathSL}templateParams/template${option}.txt |grep "#LanzExtensionProtQ" | cut -d : -f 2`
#	extensionesLig=`cat ${pathSL}templateParams/template${option}.txt |grep "#LanzExtensionligB" | cut -d : -f 2`
#	cores=`cat ${pathSL}templateParams/template${option}.txt |grep "#LanzCores" | cut -d : -f 2`
#	time_experiment=`cat ${pathSL}templateParams/template${option}.txt |grep "#LanzTimeExecution" | cut -d : -f 2`
#	grid=`cat ${pathSL}templateParams/template${option}.txt |grep "#LanzGrid" | cut -d : -f 2`
#	mem=`cat ${pathSL}templateParams/template${option}.txt |grep "#LanzMem" | cut -d : -f 2`
#	sizeGridX=`cat ${pathSL}templateParams/template${option}.txt |grep "#lanzSizeGridX" | cut -d : -f 2`
#	sizeGridY=`cat ${pathSL}templateParams/template${option}.txt |grep "#lanzSizeGridY" | cut -d : -f 2`
#	sizeGridZ=`cat ${pathSL}templateParams/template${option}.txt |grep "#lanzSizeGridZ" | cut -d : -f 2`
#
#	printHelpAux "Receptor/Query Extension:" ${extensionesProt} ""
#	printHelpAux "ligands Extension:" ${extensionesLig} ""
#	printHelpAux "Grid:" ${grid}
#	printHelpAux "Cores:" ${cores}
#	printHelpAux "RAM:" ${mem}
#	printHelpAux "TimeDocking:" ${time_experiment}" sec"
#	if [ $sizeGridX != 0 ];then #si la variable no esta en 0 queire decir que se le indico un tamaño de grid
#		printHelpAux "SizeGridX: " ${sizeGridX}" Å"
#		printHelpAux "SizeGridY: " ${sizeGridY}" Å"
#		printHelpAux "SizeGridZ: " ${sizeGridZ}" Å"
#	fi
#
#
#   	#	fi
#	#done < ${pathSL}SW.txt
#	echo""
#
#	echo -e "${YELLOW}Example: $1${NONE}"
#	echo ""
#	echo ""
#}
#___________________________________________________________________________________
#
#	Lee las posibles tecnicas BD, VS, VSB...
#___________________________________________________________________________________
#function readTechniquesSW (  )
#{
#	techniques=""
#	for i in $(ls $1);do
#
#		i="${i%.*}"
#		i="${i##*$2}"
#		techniques=${techniques}" | "${i}
#
#	done
#	techniques="[ ${techniques:2:${#techniques}-1} ]"
#	echo $techniques
#}
##
##	recoge las posibles extensiones del fichero
##
#function getOptions( )
#{
#	opt=`cat ${pathSL}templateParams/template${option}.txt |grep $1 | cut -d : -f 2`
#	opt=`echo $opt | sed s/,/" | "/g`
#	opt="[ "$opt" ]"
#	echo $opt
#
#}
#
##tmb imprime el software y en el futuro las posibles optiones
#printOp()
#{
#	printTitle
#	printHelpAux "-p"  "`getOptions "\#LanzExtensionProtQ"`"  "receptor o query"
#	printHelpAux "-l"  "`getOptions "\#LanzExtensionligB"`"  "query o carpeta de queries"
#	printHelpAux "-s" "[ $option ]"   "software a utilizar"
#	printHelpAux "-o" "`getOptions "\#lanzOptions"`"   "software a utilizar"
#	if [ -z `cat ${pathSL}templateParams/template${option}.txt |grep "\LanzCoordX" | cut -d : -f 2` ];then #softwares que requieren x y z
#		printHelpAux "-x|-X" "[  ]"  "Pos X del docking solo VS"
#		printHelpAux "-y|-Y" "[  ]"  "Pos y del docking solo VS"
#		printHelpAux "-z|-Z" "[  ]"  "Pos z del docking solo VS"
#	fi
#
#}



#
#	Lee las optiones de ShuttleMol/node_login/template
##
#printHelpSwexclusivo()
#{
#
#		if [ -z $software ];then
#			echo ""
#			printHelpAux "Opciones extras"
#			path_login_node=$pathSL
#			source ${pathSL}special_params.sh
#			read_template_params $option
#			for i in `seq 1 $contadorFile`;do
#				IFS='::' read -ra ADDR <<< "${file[$i]}"
#				printHelpAux "(${ADDR[0]})\t" "${ADDR[4]}" "${ADDR[6]}"
#			done
#		fi
#
#}
#echo ""
#option=`echo "$1" | tr '[:lower:]' '[:upper:]'`
#
#cat ${pathSL}logo_shuttlemol.txt
#
##entrada de datos desde shuttlemol
#
#  echo "Running ASGARD (Automated Scripts for Gromacs to Analyze and Running Dynamics)"
#  echo " "
#  echo -e "${PURPLE}Welcome to ASGARD${NONE} "
#  echo -e "${PURPLE}ASGARD is a group of python scripts to analyze Molecular Dynamics trajectories obtained by GROMACS${NONE} "
#  echo -e "${PURPLE}This tool is able to perform the subsequent analyses:${NONE} "
#  echo -e "${PURPLE}This tool is able to perform the subsequent analyses:${NONE} "
#  echo -e "${PURPLE}This tool is able to perform the subsequent analyses:${NONE} "
#  echo -e "${PURPLE}This tool is able to perform the subsequent analyses:${NONE} "
#  echo -e "${PURPLE}This tool is able to perform the subsequent analyses:${NONE} "
#
#if [  -z "$option" ];then
#  echo ""
#  printf "\t${CYAN}%s${NONE} \n"  "[Input options]"
#  echo ""
#  echo -e "${PURPLE}-d) Input directory${NONE} "
#  printGlobalHelp "d" ".top" "" "Topology file"
#  echo -e "${PURPLE}____________________________________________________${NONE} "
#  echo -e "${PURPLE}It must contain the next files inside them${NONE} "
#  printGlobalHelp "d" ".top" "" "Topology file"
#	printGlobalHelp "q" ".xtc" "" "Trayectory file in xtc format"
#	printGlobalHelp "t" ".gro" "" "file "
#	printGlobalHelp "o" "N" "" "Technique `readTechniquesSW ${pathSL}techniques/ SLTechnique`"
#	printGlobalHelp "s" "N" "" "Available Software `readTechniquesSW ${pathSL}templateParams/ template`"
#	printGlobalHelp "se" "O" " [ N | Y ]" "Secuencial Mode Y, queue manager N. ( Default N )."
#	printGlobalHelp "de" "O" " [ 1-10  ]" "Debug mode range[1​-​10​]. It can only be used in sequencuial mode (def​a​ult 0)."
#	printGlobalHelp "el" "O" " [ N | Y ]" "Insert data in excell bio-hpc sheet. (Defaul y)."
#	printGlobalHelp "prp" "N" "" "Protocol used to prepare recptor or query."
#	printGlobalHelp "prl" "N" "" "Protocol used to prepare ligand/s."
#
#  echo ""
#  printf "\t${CYAN}%s${NONE} \n"  "[Input options]"
#  echo ""
#  echo -e "${PURPLE}-d) Input directory${NONE} "
#  echo -e "${PURPLE}____________________________________________________${NONE} "
#  echo -e "${PURPLE}It must contain the next files inside them${NONE} "
#  printGlobalHelp "d" ".top" "" "Folder where​ test data​ will​ be​ save​d​."
#	printGlobalHelp "q" ".xtc" "" "R​eceptor​​ or query​'s​​ file​. "
#	printGlobalHelp "t" ".gro" "" "Ligand's​ file​ or directory"
#	printGlobalHelp "o" "N" "" "Technique `readTechniquesSW ${pathSL}techniques/ SLTechnique`"
#	printGlobalHelp "s" "N" "" "Available Software `readTechniquesSW ${pathSL}templateParams/ template`"
#	printGlobalHelp "se" "O" " [ N | Y ]" "Secuencial Mode Y, queue manager N. ( Default N )."
#	printGlobalHelp "de" "O" " [ 1-10  ]" "Debug mode range[1​-​10​]. It can only be used in sequencuial mode (def​a​ult 0)."
#	printGlobalHelp "el" "O" " [ N | Y ]" "Insert data in excell bio-hpc sheet. (Defaul y)."
#	printGlobalHelp "prp" "N" "" "Protocol used to prepare recptor or query."
#	printGlobalHelp "prl" "N" "" "Protocol used to prepare ligand/s."
 
 
#
#	echo ""
#	printf "\t${CYAN}%s${NONE} \n"  "[ O=optional, D=Depends, N=Necesary ]" 
#	echo ""
#	echo -e "${PURPLE}Global Options${NONE} "
#	echo -e "${PURPLE}____________________________________________________${NONE} "
#	printGlobalHelp "d" "O" "" "Folder where​ test data​ will​ be​ save​d​."
#	printGlobalHelp "q" "D" "" "R​eceptor​​ or query​'s​​ file​. "
#	printGlobalHelp "t" "N" "" "Ligand's​ file​ or directory"
#	printGlobalHelp "o" "N" "" "Technique `readTechniquesSW ${pathSL}techniques/ SLTechnique`"
#	printGlobalHelp "s" "N" "" "Available Software `readTechniquesSW ${pathSL}templateParams/ template`"
#	printGlobalHelp "se" "O" " [ N | Y ]" "Secuencial Mode Y, queue manager N. ( Default N )."
#	printGlobalHelp "de" "O" " [ 1-10  ]" "Debug mode range[1​-​10​]. It can only be used in sequencuial mode (def​a​ult 0)."
#	printGlobalHelp "el" "O" " [ N | Y ]" "Insert data in excell bio-hpc sheet. (Defaul y)."
#	printGlobalHelp "prp" "N" "" "Protocol used to prepare recptor or query."
#	printGlobalHelp "prl" "N" "" "Protocol used to prepare ligand/s."
#	printGlobalHelp "hi" "O" " [ Y | N ]" "At the end of calculations generate graphs and pymol session (Default n)."
#	printGlobalHelp "v" "O" "" "Show version of shuttlemol."
#	printGlobalHelp "h" "O" " [ | Opt | SW ]" "Show help or specific help for a technique or a program." 
#	printGlobalHelp "h" "O" " [ allsw ]" "Show all programs​ (no​t​ all software​​ work​)​." 
#	printGlobalHelp "h" "O" " [ alltc ]" "Show all techniques​ (no​t​ all techniques​ work​)​." 
#	printGlobalHelp "h" "O" " [ prtocols ]" "Show all protocols to prepare proteins and ligands." 
#	printGlobalHelp "em" "O" "" "Send email when finish all jobs. (only with -hi y parameter)" 
#
#	
#	echo ""
#	echo -e "${PURPLE}Queue Manager Options${NONE} "
#	echo -e "${PURPLE}____________________________________________________${NONE} "
#	printGlobalHelp "td" "O" "" "Set e​xecution time of a program. Default ./sm.sh -h software (Not all programs have this option)."
#	printGlobalHelp "tj" "O" "" "Time allocated to the jo​​b​.​ ​Default multiply -td  by runs in a jo​b​.​ ​F​o​r​ mat 00:05:00"
#	printGlobalHelp "mm" "O" "" "Memory reserved for a job."
#	printGlobalHelp "j" "O" "" "Number of jobs that will​ be​ sen​t​ to supercomputer "
#	printGlobalHelp "ni" "O" "" "Priority of the job managers, Sbatch [0-10000] being 10000 the lowest. Qsub [-1024- + 1023] default or higher priority the higher the number."
#	printGlobalHelp "co" "O" "" "Number of cores for execution​.​ Importan​t:​ can't reserve more cores than contain​ed in​ the node"
#	printGlobalHelp "nm" "O" "" "Number of nodes for execution."
#	printGlobalHelp "gp" "O" "" "Number of GPUs for executio​.​ Important​: ​can't reserve th​e​ s​e​resources if the supercomputer don​'​t have​ them​."
#	printGlobalHelp "nj" "O" "" "Name of job."
#
#	echo ""
#	echo -e "${PURPLE}Docking Options${NONE} "
#	echo -e "${PURPLE}____________________________________________________${NONE} "
#	printGlobalHelp "nc" "O" "" "Number of conformations for output docking."
#	printGlobalHelp "x" "D" "" "X coordinate of the center of grid."
#	printGlobalHelp "y" "D" "" "Y coordinate of the center of grid."
#	printGlobalHelp "z" "D" "" "Z coordinate of the center of grid."
#	#printGlobalHelp "fx" "O" "" "File of flexibility for vina."
#	printGlobalHelp "an" "O" "" "Size of cube when used BDC."
#	printGlobalHelp "gx" "O" "" "Size X of grid for some programs os docking."
#	printGlobalHelp "gy" "O" "" "Size Y of grid for some programs os docking."
#	printGlobalHelp "gz" "O" "" "Size Z of grid for some programs os docking."
#	printGlobalHelp "be" "O" "" "Exhaustiveness en BD Ddefault(1)"
#	printGlobalHelp "bda" "O" "" "Atomo docking in BD Ddefault(CA)"
#
#		
#
# 	echo ""
#	echo -e "${PURPLE}Opciones Similaridad${NONE} "
#	echo -e "${PURPLE}____________________________________________________${NONE} "
#	printGlobalHelp "sc" "O" "" "Threshold for some similarity programs. Results lower to this score will be deleted."
#	echo ""
#	echo ""

	

#	echo ""
#else
#option=`echo "$1" | tr '[:lower:]' '[:upper:]'`
#
#case $option in
#	#
#	#	Tecnicas de shuttlemol
#	#
#	GR)
#		printOp
#		
#		
#		printExample "${name_shuttlemol} -p targets/3s5z.pdb -l queries/testDm2/ -o VS -s GR -j 1  [-GP]"
#
#		echo -e "\tError Comun cuando se modifica/genera la proteina con maestro:"
#		echo -e "\tmodificar todos los NMA por NME y cambiar"
#		echo -e "\t\tATOM   5228  CA  NMA A1020A    -72.002   7.253 -39.562  1.00  0.00           C  "
#		echo -e "\tpor"
#		echo -e "\t\tATOM   5228  CH3 NME A1020A    -72.002   7.253 -39.562  1.00  0.00           C  "
#		echo ""
#		echo -e "${GREEN}Ojo la proteina debe ser pdb y el ligando mols"
#
#
#		echo -e "${GREEN}Ojo para BIPHSIC_SYSTEMS el primier ligando que se duplicara como debe tener el nombre menor ordenado alfabeticamente en la carpeta"
#
#		echo -e "python ${path_external_sw}/gromacs/topology/generate_topology.py -t projects/nlrp3/nlrp3.pdb -q projects/nlrp3/lig/ -p TARGET_QUERY -g gmx_mpi ${NONE}"
#
#		echo ""
#;;
#	GRN)
#		printOp
#		printHelpAux "-d" "[    ]" "Carpeta de la prueba realizada previamente"
#		printExample "${name_shuttlemol} -p protein -l folderLiagnd -o VS -s GRN -j X -prl X -prp X -stepSimulacionDM X [-GP [1-2]] \nExample: ${name_shuttlemol} -t projects/testgrn/1le0.pdb  -o VS -s GRN -prp na -prl na  -el n -j 1 -se y -de 100 -step_md 100  -d VS_GR_1le0_target_2019-03-18 -gp 2 -co 24";;
#	GRC)
#		printOp
#		printHelpAux "-d" "[    ]" "Carpeta de la prueba realizada previamente"
#
#		printExample "${name_shuttlemol} -p protein -l folderLiagnd -o VS -s GRN -j X -prl X -prp X [-GP [1-2]] \nExample: ./sm.sh -p proyectos/3s5z/3s5z_b.pdb -l proyectos/3s5z/queries/ -o VS -s GRC -el n -j 1 -prp gromacs -prl gauss  -d test3s5z2/  -se y -de 10";;
#	*)
#		printHelpAux "No se encuentra ayuda"
#esac
#exit
#fi
