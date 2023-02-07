#!/usr/bin/env bash
#
#   Description: Parametros de ASGARD
# ______________________________________________________________________________________________________________________


#______________________________________________________________________________________________________________________
#
#			Parametros posibles
#_________________________________________________________________________________________________________________________
read_file_template()
{
	source ${pathSL}special_params.sh
	read_template_params $software

	if [ "$error" == "0" ];then

		for i in `seq 1 $contadorFile`;do
			error=12
			IFS='::' read -ra ADDR <<< "${file[$i]}"
			if [ "$1" == "${ADDR[4]}" ];then
				if [ "${ADDR[2]}" == "Y" ];then
					error=0		
					allComand=`echo "$allComand" |sed -e "s,${ADDR[4]} $2,,g"` #se elimina de all comand y se guarda en parametros especiales
					optAdicionals="$optAdicionals ${ADDR[4]} $2"
					shift
					break
				else
					error=0
					if [[ "$1" != \-* ]];then #esto quiere decir que sera un parametro con un delante -
						error=12 #si es parametro  el siguiente debe empezare por guion
					fi
					allComand=`echo "$allComand" |sed -e "s,${ADDR[4]},,g"` #se elimina de all comand y se guarda en parametros especiales
					optAdicionals="$optAdicionals ${ADDR[4]}"
					break
				fi
			fi
		done
	fi

}
function empty_variable()
{
	if [ -z $2 ];then
		val="N/A"
	else
		val=$2
	fi
    if [ -z "$1" ];then
        echo $val
    else
        echo $1
    fi
}
#
#	Parametros por defecto (deberian aparecer todos e incializarlos para que no den problemas)
#


#
#if [ -z $CWD ];then
#	 CWD=${PWD}/ 	                         	#	Si se llama dentro de un job la ruta sera $HOME, por eso se le puede indicar previamente
#fi 						                        #	Si esta vacio la rutala ruta del script para trabajar con rutas completas en los jobs y asi evitar posibles errores
#
#if [ -z "${name_job}" ];then
#    name_job=""			                        	#	nombre del job paara indicarselo al gestor de colas
#fi
#
#folder_experiment=`empty_variable $folder_experiment` # Carpeta de la prueba
#		                        
#	                        	
#
#
#scoreCorte=`empty_variable $scoreCorte 0` 		#	para algunos softwares de similaridad (Se usa poco)
#numPoses=`empty_variable $numPoses 1`			#	numero de conformaciones que se genraran en softwares de coking (vina, LF )
#torsion=`empty_variable $torsion 12`			# 	no se usa, antiguamente era para indcar el maximo nuemro de torsiones a vina
#numAminoacdo=`empty_variable $numAminoacdo 0`
#
#
#if  [ -z "$optAdicionals" ];then
#	optAdicionals=""
#fi
#if [ -z "$allComand" ];then
#    allComand="$0 $@"
#fi
#
#
#nodos=`empty_variable $nodos 1` 				#	por defecto siempre se solicitara un nodo solo
#
#debug=`empty_variable $debug -1` 			    #	solo funciona con secuencial y es para ir mostrando los comandos que va utilizando shuttlemol
#
#name_target=`empty_variable $name_target no_target `       #	nombre de la target sin extension ni ruta
#name_query=`empty_variable $name_query`         #	nombre query sin extension ni patch
#check=`empty_variable $check`
#disease=`empty_variable $check`                 #	enfermedad, No se usa
#localitation=`empty_variable $localitation`		#	localizacion del cluster
#proteinName=`empty_variable $proteinName`		#	nombre de la target
#subset=`empty_variable  $subset`		    	#	datos de excel
#flex=`empty_variable $flex`				        #	flexibilidad vina solo
#flexFile=`empty_variable $flexFile`			    #	si con vina se le introduce un fichero de flexibilidad generado previamente, defecto no
#secuencial=`empty_variable $secuencial `        #	shuttlemol sin enviar jobs al cluster
#dynamic=`empty_variable $dynamic `              #	flexibilidad dinamica necesita conexsion con el servidor 1and1
#outJob=`empty_variable $outJob `
#histograms=`empty_variable $histograms `        #	Generar historgramas o no
#grid=`empty_variable $grid `                    #	parametro para saber si el software usa grids
#cores=`empty_variable $cores `	                #   num de cores a utilizar
#mem=`empty_variable $mem `                      #   num memoria a utilizar
#time_job=`empty_variable $time_job `              # 	tiempo del job
#versionHelp=`empty_variable $versionHelp`       #	cuando introduce la option -v
#renice=`empty_variable $renice`                 #	variable para la prioridad del JOB
#GPU=`empty_variable $GPU`                       #	Usar GPU o no parametro numerico
#time_experiment=`empty_variable $time_experiment`       #	900 segundos = 15 minutos ## tiempo maximo que se hara dde una ejecucion de docking en el job (Por si se queda pillado que no fallen las demas ejecuciones)
#specialCommand=`empty_variable $specialCommand` #   comandos especiales de cada plantilla
#ext_query=`empty_variable $ext_query`           #   extension del query
#ext_target=`empty_variable $ext_target`               #   extension del receptor
#write_excel=`empty_variable $write_excel n`     #	escribir en el excel (n,, proyecto o vacio)
#resName=`empty_variable $resName`               #	nombre del residuo
#num_amino_acid=`empty_variable $num_amino_acid 0`   #   num Aminoacido
#gridSizeX=`empty_variable $gridSizeX`           #	Tamaño de la grid X
#gridSizeY=`empty_variable $gridSizeY`           #	Tamaño de la grid Y
#gridSizeZ=`empty_variable $gridSizeZ`           #	Tamaño de la grid Z
#queue=`empty_variable $queue`
#protocolP=`empty_variable $protocolP na`           #	protocolo de conversion target
#protocolL=`empty_variable $protocolL na`           #	protocolo de conversion liagndo
#kill_sm=`empty_variable $kill_sm`                   #	No utilizar, mata los procesos de lanzador
#chain=`empty_variable $chain`                   #	cadena de la target
#lanzTimeOut=`empty_variable $lanzTimeOut`       #	time out del comanto
#email=`empty_variable $email`                   #    se le envia un email cuando termian      #	email (Nop funcionara en todos los culter)
#number_execution=`empty_variable $number_execution -1` #numero de la ejecucion par aBD
#bd_exhaustiveness=`empty_variable $bd_exhaustiveness  1` #exaustiovidad en BD cada 1 CA, cada 2 CA, cada X CA ...
#bd_atom_default=`empty_variable $bd_atom_default CA`                       #residuo por defecto que se hace docking en BD
#target=`empty_variable $target no_target`
#mode_test=`empty_variable $mode_test`           # Ejecuta en secuencial un pequeño BD

#
#	parametros del usuario e internos
#
while (( $# )) #recorro todos los parametros y los asigno
 do

    #if [[ "$1" == \-* ]];then #esto quiere decir que sera un parametro con un delante -
    if [[ "$1" == \-[a-z]* ]] || [[ "$1" == \-[A-Z]* ]] || [[ "$1" == \-\-[a-a]* ]] || [[ "$1" == \--[A-Z]* ]];then #esto quiere decir que sera un parametro con un delante -
	   case `printf "%s" "$1" | tr '[:lower:]' '[:upper:]'`  in
		   	#entrada de datos desde shuttlemol
			#-F)   folder_analysis=${2%/}/;; # Carpeta donde se encuentran los resultados de MD 
#      -trj) TRAJ=$2
#      -tpr) TOP=$2
#      -gro) GRO=$2
      -F)   folder_analysis=${2%/};; # Carpeta donde se encuentran los resultados de MD
      -P)   profile=$2;; # Perfil (TARGET, TARGET_PROTEIN, TEST...)
      -D )  folder_experiment=${2%/}/;;	#directorio donde se guardaran los datos (no necesario)
       
      
      
      
#      -X )  x=$2;;			#posicion centro x
#			-Y )  y=$2;;			#posicion centro y
#			-Z )  z=$2;;			#posicion centro z
#			-D )  folder_experiment=${2%/}/;;	#directorio donde se guardaran los datos (no necesario)
#			-P|-NP ) echo -e "\nERROR: $1"
#			         echo "(-p) Opcion no valida ahora se usa -t"
#			        echo -e " ERROR: (-nt) Opcion no valida ahora se usa -nt\n"
#                 exit;;
#			-T )  target=$2;;		#fichero target o query
#			-L|-NL ) echo -e "\nERROR: $1"
#			         echo  "(-l) Opcion no valida ahora se usa -q"
#			         echo -e " ERROR: (-nl) Opcion no valida ahora se usa -nq\n"
#			         exit;;
#    		-Q )  query=$2;;		#directorio o fichero query
#			-SE) secuencial=$2;;	#para ejecutar shuttlemol en secuencial
#			-FX) flex=$2;;			#para indicarle que se queire flexibilidad (Vina) NO SE USA
#			-FL) flexFile=$2;; 		#Indicarle el fichero de flexibilidad previamente generado(prepare_flexreceptor4.py) cuando se divide en 2 -p 3s5z_rig -fl 3s5z_flex
#			-DY)  dinamyc="-dy";;	#para indicarle flexibilidad dinamica, tiene que existir un BD previamente (Vina)
#			-DD)  dirDinamyc=$2;;	#para indicarle el directorio de la flex (Se usa con dy)
#			-FI) ficheroBD=$2 ;;    #fichero de flexibilidad
#			-NC )  numPoses=$2 ;;	#numero de poses para la prueba defecto (1)
#			-AN)  t=$2;;  			##tamaño para el lado cuando se ahce BDC
#			-HI) histograms=$2;; 		##generar histograma despues de shuttlemol
#			-O )option=`echo $2 | awk '{print toupper($0)}'`   ;;			#bd | BDVS | BDC |o vs
#			-S )  software=`echo $2 | awk '{print toupper($0)}'`;;		 #software a utilizar
#			#-pr|-PR) nombreProject=$2;;		#nombre proyecto para indicar al cluster
#			-PN) proteinName=$2;;	#nombre de la target para subir al excell
#			-J )  num_per_job=$2;;	#numero de ejecuciones por job o en docking mutacion numero de distancia
#			-DE ) debug=$2;; 		#modo debug (del 1 hasta el 10 mirar shuttlemol/node_login/debug.sh)
#			-V  ) versionHelp="ve";;
#			-EL) write_excel=$2;;
#			-PRP) protocolP=$2;; 	#protocolo de conversion target
#			-PRL) protocolL=$2;; 	#protocolo de conversion liagndo
#            -NN )  nodos=$2;;
#			-EM ) email=$2;;
#			##uso interno
#			-RE) resName=$2;;		 #nombre del residuo
#			-NT) name_target=$2;;	 #nombre de la target a secas
#			-O ) option=$2;;		 #BD o VS
#			-I ) confAdicionalI=$2;; ##campos adicionales para VS o BD
#			-F ) confAdicionalF=$2;; ##campos adicionales para VS o BD
#			-C ) CWD=$2;;			 ##ruta donde se ejcuta
#			-NJ) name_job=$2;;		 #nombre del job identificativo
#			-BE) bd_exhaustiveness=$2;; # cada cuantos CA se hace docking
#			-UL) ultimo=$2;;		 #si se le pasa el -ul quiere decir que es el ultimo job
#			-NA) num_amino_acid=$2;;  #num de aminoacido de la prot
#			-CH) chain=$2;;			 #cadena de la target
#			-CK) check="ck";;
#			-TR) torsion=$2;;
#			-TD) time_experiment=$2;;    #tiempo que puede tardar como maximo cada docking, si tarda mas se cortara con la orden timeout "comando".
#			-NI) renice=$2;;         #sbatch a prioridad puede ir de 0 a 10000 siendo 10000 la más baja.    qsub es un número entre -1024 y +1023. Por defecto se asigna 0, y la prioridad es mayor cuanto más alto el número.
#			-CO) cores=$2;;			#nimero de cores a utlizar
#			-GP) GPU=$2;;			# GPU (SOLO EN MALAGA Y CON Gromacs)
#			-SC ) scoreCorte=$2;;		#corte para cuando de menos de x score eliminar la prueba
#			-MM ) mem=$2;;   		#memoria ram del docking
#			-EJ ) number_execution=$2 ;;  # #numero de la ejecucion, util para BD (el numeor de atomo no puede ser poe que algun pdb esta mal y pueden tener 2 carbonos alfas con el mismo num atomo)
#			-NQ ) name_query=$2;;		#nombre del liagndo
#			-IN ) ini=$2;;			    #iniicio para VS de los ficheros
#			-FN ) fin=$2;;      		#final para VS de los ficheros
#			-SA ) outJob=$2;;   		#nombre con el que se guardaran las salidas del job
#			-TJ ) time_job=$2;;			#tiempo en HH:MM:SS que puede durar el job
#			-LTO ) lanzTimeOut=$2;;
#			-PROJECT) project=$2;;  		#nombre del proyecto para el gesto de colas
#			-QUEUE  ) queue=$2;;
#			-GRID )grid=$2;;
#			-BDA ) bd_atom_default=$2;;
#			-EXP ) ext_target=$2;;
#			-EXL ) ext_query=$2;;
#			-GX ) gridSizeX=$2;;	#Tamaño de la grid X
# 	        -GY ) gridSizeY=$2;;	#Tamaño de la grid Y
#            -GZ ) gridSizeZ=$2;;	#Tamaño de la grid Z
#            -K )   kill_sm="Y";; #mata los procesos de shuttlemol, util cuando esta en secuendcial, ojo si esta ejecutando dos shuttlemol a la vez matara los dos
#            -PROFILE)  profile=$2;;
#			-TEST) mode_test="Y"
#			  option=${2};;
#			-h|-H ) f_help $2;;
#
#			*)
#                path_login_node=${pathSL}
#				read_file_template $1 $2
#				paraMError="para option $1"

	    esac
	fi
  shift #Modificamos $2 a $1, $3 a $2 ....
done

if [[ "${profile}" == "" ]];then
	profile="STANDAR_${option}"		#	Perfil para usar en get_histogram
fi

if [[ ${software} == "GR" ]] || [[ ${software} == "GRN" ]] || [[ ${software} == "GRC" ]] ;then
	
    if [[ "$query" == "no_query" ]] || [[ "$target" == "no_target" ]] || [[ "$query" == "" ]];then
    	
        if [[ $target != "no_target" ]] || [[ $target != "" ]];then
        	
        	aux="${target%.*}"
            mode_gr=`cat ${aux}*.conf |head -1 |grep  Profile | awk -F\: '{print $2}'`

        else
            mode_gr=`cat ${query}*.conf |head -1 |grep  Profile | awk -F\: '{print $2}'`
        fi
    fi
    
 
    if [[ "${mode_gr}" == "QUERIES" ]] || [[ "${mode_gr}" == "BIPHSIC_SYSTEMS" ]];then
        target="targets/test/queries.pdb"
        if [[ ! -f $target ]];then
            echo "ERROR: Para el modo QUERIES | BIPHSIC_SYSTEMS"
            echo "ERROR: debe existir un fichero pdb placebo en "$target
            exit
        fi
    elif [[ "${mode_gr}" == "TARGET" ]];then
	
        query="ShuttleMol/external_sw/gromacs/config_files/target/"
        aux=`ls ShuttleMol/external_sw/gromacs/config_files/target/*.mol2 |wc -l`
        if [[ ${aux}  != 1 ]];then
            echo "ERROR: Para el modo target"
            echo "ERROR: debe existir solo un fichero .mol2  "+query
            exit
        fi
    fi
fi




