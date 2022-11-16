#!/bin/bash
#
#	Script para Continuar Ejecucion de Gromacs Automaticamente
#	informacion sacada de:
#	https://cluster.earlham.edu/wiki/index.php/Checkpoint_and_Restarting
#
BLUE='\033[0;34m'
NONE='\033[00m' #no color
topologias="N/A"

#Fatal error:
#Too many LINCS warnings (5896)
#If you know what you are doing you can adjust the lincs warning threshold in your mdp file
#or set the environment variable GMX_MAXCONSTRWARN to -1,
#but normally it is better to fix the problem
#For more information and tips for troubleshooting, please check the GROMACS
#website at http://www.gromacs.org/Documentation/Errors
export GMX_MAXCONSTRWARN=-1 
execute()
{
	echo -e "${BLUE}$1${NONE}"
	eval "$1"
	error=$?
	echo -e "${BLUE}ERROR: $error ${NONE}"
	if [ "$error" != "0" ];then
			exit;
	fi
}
ejecutar()
{

if [ "${nodos}" != "1" ];then
	cores=`expr ${cores} \* ${nodos}` #En caso de que se usen varios nodos mucltiplico cores por nodos
	thrads="-ntmpi $cores "
else
	thrads="-ntomp $cores "
fi

IFS='-' read -r -a array <<< "$opt_aux"
#
#	Leo los dos posibles parametros de entrada
#
for element in "${array[@]}"
do
	aux=`echo $element | cut -f1 -d' '`
	param=`echo $element | cut -f2 -d' '`
	case "$aux" in
		versionGromacs)			versionG=$param;;
		mode)					mode=$param;;
	esac
done
IFS=$OLDIFS;
#
#	Dependiendo de la version de gromascs
#
if [ "$versionG" == "4" ];then
 	mpi="_mpi"
 	gromacs5=""
 	graph="g_"
elif [ "$versionG" == "5" ];then
 	mpi=""
 	gromacs5="gmx_mpi "
 	graph=""
else
 	echo "Verion ${versionG} Gromacs incorrecta " 
 	exit
fi 
#
#	GPUR o no
#
if [ "${GPU}" != "N/A" ];then
	echo "ENTR A GPU"
	gpu="-nb gpu"
else
	echo "ENTR NO GPU"
	gpu=""
fi

salida=${salida%.*}
salida=`echo "${salida}" | sed -e "s,GRC,GR,g"`
out_grid=${folderGrid}${salida}
outTxt=${folderTxt}${salida}
out_aux=${folderOutUcm}${salida}
out_molec=${folderMolec}${salida}
out_energies=${folderEnergy}${salida}


execute "$gromacs5 mdrun${mpi} -s ${out_molec}_mdrum_simulation.tpr -cpi ${out_molec}_mdrum_simulation.cpt ${gpu} ${thrads} -append -deffnm ${out_molec}_mdrum_simulation -g ${out_aux}_mdrum_simulation.log"

#
#	Centrar simulacion y genera el pdb
#
printf "Protein\nSystem\n" | gmx_mpi trjconv -s ${out_molec}_mdrum_simulation.tpr -f ${out_molec}_mdrum_simulation.xtc -center -ur compact -pbc mol -o ${out_molec}_mdrum_simulationOriginal_1.xtc
printf "Backbone\nSystem\n" | gmx_mpi trjconv -s ${out_molec}_mdrum_simulation.tpr  -f ${out_molec}_mdrum_simulationOriginal_1.xtc -fit rot+trans -o ${out_molec}_mdrum_simulation_center.xtc
execute "rm  ${out_molec}_mdrum_simulationOriginal_1.xtc"
execute "$gromacs5 editconf${mpi} -f ${out_molec}_mdrum_simulation.gro -o ${out_molec}_mdrum_simulation.pdb" #no hace falta pero convierto la simulacion a pdb

#
#	Esto en el futuro hay que ShuttleMol como job
#
if [ "${mode}" != "" ];then
	echo  "python ${pathSW}gromacs/analizarResults/getResults.py ${out_molec} ${mode} $versionG"
	python ${pathSW}gromacs/analizarResults/getResults.py ${out_molec} ${mode} $versionG 
fi
exit


}


#
#	Convierto a pdb y centro las aguas
#
#salida=`echo "${salida}" | sed -e "s,GRC,GR,g"`
#out_molec=`echo "${out_molec}" | sed -e "s,GRC,GR,g"`
#aux=`echo "${folderMolec}" | sed -e "s,GRC,GR,g"`
#out_aux=`echo "${out_aux}" | sed -e "s,GRC,GR,g"`
#out_molec=${aux}/${ini}/${salida} ##TEST
#echo "$gromacs5 mdrun${mpi} -s ${out_molec}_mdrum_simulation.tpr -cpi ${out_molec}_mdrum_simulation.cpt ${gpu} ${thrads} -append -deffnm ${out_molec}_mdrum_simulation -g ${out_aux}_mdrum_simulation.log"
#cd $CWD 					##TEST
#mv $out_molec* $folderMolec 	##TEST
#printf "Protein\nSystem\n" | gmx_mpi trjconv -s ${out_molecAux}_mdrum_simulation.tpr -f ${out_molecAux}_mdrum_simulationOriginal.xtc -center -ur compact -pbc mol -o ${out_molecAux}_mdrum_simulationOriginal_1.xtc
#printf "Backbone\nSystem\n" | gmx_mpi trjconv -s ${out_molecAux}_mdrum_simulation.tpr  -f ${out_molecAux}_mdrum_simulationOriginal_1.xtc -fit rot+trans -o ${out_molecAux}_mdrum_simulation.xtc
#execute "$gromacs5 editconf${mpi} -f ${out_molecAux}_mdrum_simulation.gro -o ${out_molecAux}_mdrum_simulation.pdb" #no hace falta pero convierto la simulacion a pdb