#!/usr/bin/env bash
#
#	Genera fichero de grid, la grid y el fichero de docking para autodock4
#_______________________________________________________________________________________________________________________

#
#		Fichero Grid
#_____________________________________________________________________________________________

if [ -f ${CWD}${query}  ];then
#busco en la target los tipos de atomos que son
typeAtomsProt=`cat ${CWD}${target} |awk '{print  substr($0,78)}' |grep -v Type |grep -v ^_ |sort -V -u |sed '/^$/d' |tr '\n' ' '`  ##busca los typos de atomos de la prot
#busco en el query los tipos de atomos que son
typeAtomsLig=`cat ${CWD}${query}  |awk '{print  substr($0,78)}' |grep -v Type |grep -v ^_ |sort -V -u |sed '/^$/d' |tr '\n' ' '` 	##busca tipos de atomos del lig



echo "npts $gridSizeX $gridSizeY $gridSizeZ                        # num.grid points in xyz"    > ${out_aux}_grid.gpf
echo "gridfld ${out_aux}.maps.fld           # grid_data_file" 							        >>	${out_aux}_grid.gpf
#echo "spacing 0.375                        # spacing(A)	"								        >>	${out_aux}_grid.gpf
echo "spacing 0.19                        # spacing(A)	"								        >>	${out_aux}_grid.gpf
echo "receptor_types ${typeAtomsProt}       # receptor atom types"						        >>	${out_aux}_grid.gpf
echo "ligand_types ${typeAtomsLig}         # ligand atom types"							        >>	${out_aux}_grid.gpf
echo "receptor ${CWD}${target}                 # macromolecule"						            >>	${out_aux}_grid.gpf
echo "gridcenter ${x} ${y} ${z}            # xyz-coordinates or auto "					        >>	${out_aux}_grid.gpf
echo "smooth 0.5                           # store minimum energy w/in rad(A) #"		        >>	${out_aux}_grid.gpf
arr=$(echo $typeAtomsLig | tr " " "\n")													        #divido los atomosLigs
for sal in $arr																			
do
	echo "map $out_aux.${sal}.map                    # atom-specific affinity map	" 	        >>	${out_aux}_grid.gpf
done
echo "elecmap ${out_aux}.e.map              # electrostatic potential map "				        >>	${out_aux}_grid.gpf
echo "dsolvmap ${out_aux}.d.map             # desolvation potential map "				        >>	${out_aux}_grid.gpf
echo "dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, constant " 	        >>	${out_aux}_grid.gpf


#_________________________________________________________________________________________________________________
#
#	Fichero de docking
#___________________________________________________________________________________________________________________
echo "autodock_parameter_version 4.2       # used by autodock to validate parameter set #	"				>  ${out_aux}_dock.dpf
echo "outlev 1                             # diagnostic output level						"				>> ${out_aux}_dock.dpf
echo "intelec                              # calculate internal electrostatics				"				>> ${out_aux}_dock.dpf
#echo "seed pid time                        # seeds for random generator						"				>> ${out_aux}_dock.dpf
echo "ligand_types ${typeAtomsLig}         # atoms types in ligand 							"				>> ${out_aux}_dock.dpf
echo "fld ${out_aux}.maps.fld           	   # grid_data_file #								"				>> ${out_aux}_dock.dpf
arr=$(echo $typeAtomsLig | tr " " "\n")													#divido los atomosLigs
for sal in $arr																			
do
#	echo "eno"
	echo "map $out_aux.${sal}.map                    # atom-specific affinity map			" 				>>	${out_aux}_dock.dpf
done

echo "elecmap ${out_aux}.e.map              # electrostatics map#							"				>> ${out_aux}_dock.dpf
echo "desolvmap ${out_aux}.d.map            # desolvation map#								"				>> ${out_aux}_dock.dpf
echo "move ${CWD}$query                    # small molecule#								"				>> ${out_aux}_dock.dpf
#echo "about $x $y $z 				       # small molecule centeri #igualk centro grid 	"				>> ${out_aux}_dock.dpf
#echo "tran0 random                         # initial coordinates/A or random				"				>> ${out_aux}_dock.dpf
#echo "quaternion0 random                   # initial orientation							"				>> ${out_aux}_dock.dpf
#echo "dihe0 random                         # initial dihedrals (relative) or random			"				>> ${out_aux}_dock.dpf
#echo "torsdof ${torsiones}                 # torsional degrees of freedom #torsionLigand	"				>> ${out_aux}_dock.dpf
echo "seed 2015 2015					   #semilla											"				>> ${out_aux}_dock.dpf
#echo "rmstol 2.0                           # cluster_tolerance/A 							"				>> ${out_aux}_dock.dpf
#echo "extnrg 1000.0                        # external grid energy							"				>> ${out_aux}_dock.dpf
#echo "e0max 0.0 10000                      # max initial energy; max number of retries		"				>> ${out_aux}_dock.dpf
echo "ga_pop_size 150                      # number of individuals in population			"				>> ${out_aux}_dock.dpf
echo "ga_num_evals 2500000                 # maximum number of energy evaluations			"				>> ${out_aux}_dock.dpf
echo "ga_num_generations 27000             # maximum number of generations 					"				>> ${out_aux}_dock.dpf
#echo "ga_elitism 1                         # number of top individuals to survive to next generation 	"	>> ${out_aux}_dock.dpf
#echo "ga_mutation_rate 0.02                # rate of gene mutation 							"				>> ${out_aux}_dock.dpf
#echo "ga_crossover_rate 0.8                # rate of crossover 								"				>> ${out_aux}_dock.dpf
#echo "ga_window_size 10                    # 												"				>> ${out_aux}_dock.dpf
#echo "ga_cauchy_alpha 0.0                  # Alpha parameter of Cauchy distribution			"				>> ${out_aux}_dock.dpf
#echo "ga_cauchy_beta 1.0                   # Beta parameter Cauchy distribution 			"				>> ${out_aux}_dock.dpf
echo "set_ga                               # set the above parameters for GA or LGA 		"				>> ${out_aux}_dock.dpf
#echo "sw_max_its 300                       # iterations of Solis & Wets local search 		"				>> ${out_aux}_dock.dpf
#echo "sw_max_succ 4                        # consecutive successes before changing rho 		"				>> ${out_aux}_dock.dpf
#echo "sw_max_fail 4                        # consecutive failures before changing rho 		"				>> ${out_aux}_dock.dpf
#echo "sw_rho 1.0                           # size of local search space to sample 			"				>> ${out_aux}_dock.dpf
#echo "sw_lb_rho 0.01                       # lower bound on rho 							"				>> ${out_aux}_dock.dpf
#echo "ls_search_freq 0.06                  # probability of performing local search on individual 		"	>> ${out_aux}_dock.dpf
echo "set_psw1                             # set the above pseudo-Solis & Wets parameters 	"				>> ${out_aux}_dock.dpf
#echo "unbound_model bound                  # state of unbound ligand 						"				>> ${out_aux}_dock.dpf
echo "ga_run 10                            # do this many hybrid GA-LS runs					"				>> ${out_aux}_dock.dpf
echo "rmstol 0.5                           # cluster_tolerance/A 							"				>> ${out_aux}_dock.dpf
echo "analysis                             # perform a ranked cluster analysis				"				>> ${out_aux}_dock.dpf


fi




