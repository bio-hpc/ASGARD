#! /bin/bash																												

#======================================================================#
#======================================================================#
#
# Calculate average electron density map and structure factors from an 
# MD trajectory.
#
# Returns:
#			md_avg.map - average ASU electron density map in ccp4 map format
#			md_avg.mtz - structure factors of md_avg.mtz
#
# Author v1:James Holton
#        v2: Pawel Janowski
#
#Follow these four steps after editing setup variables in each script.
#1. Create a trajectory that is mass centered on the original crystal 
#   supercell. For this you can use fit.nc produced by XtalAnalyze.sh or
#   make your own in a similar way (usually XtalAnalyze.sh is run on a 
#   stripped trajectory, but for the density calculation one may wish to
#   include bulk water).
#
#2. MakePdb4Map.py - convert the trajectory into pdbs that will be used 
#   for the average. Add the corrrect cell parameters, remove virtual 
#   atoms, make sure element type column is correct.
#
#3. (optional) Prepare exp_pdb - this is the experimental map that 
#   you will be aligning against if align=1 below. This is not strictly 
#   necessary because the center of mass alignment results in optimal shifts 
#   of (0 0 0), ie no need to align if mass centered. If you decide to
#   use it the experimental map pdb needs to be similar to the pdb's 
#   from step two above: occupancy set to 1, column 78 filled, remove 
#   virtual atoms, add CRYST1 record.
#
#4. Set the following variable and execute md2map.sh.
#
#               *** USER SET THE FOLLOWING VARIABLES ***
#
# Directory with simulation frames (pdb format, must have element column,
#   CRYST1 of supercell, usually made by MakePdb4Map.py ). This directory 
#   should not contain any other files.
sim_pdb_dir=/home/pjanowsk/c/Case/4lzt/4lzt_ff12SB/analysis_i/average_density/test/PDBData/  
# supercell multiplicity (x y z)
MD_mult=( 3 2 2 )			         
# ccp4 space group symbol, no spaces                   
SG="P1"					                                  
# grid spacing for SFALL map 
GRID="GRID 90 90 90"		                    
# flat B-factor added to frames to avoid singularity in SFALL map calculation
B=15                                              
# structure factor calculation resolution limit
reso=0.9                                          
# set to 1 if request correlation alignment
align=0					                                  
# experimental supercell structure (CRYST1 of supercell); only used if align=1
exp_pdb=fav8_Ex2_4md2map.pdb
# if the simulation was run at constant volume, this can optionally be set to 1
# and the experimental CRYST1 pdb card provided. This obviates the need for 
# CRYST1 cards in the simulation frames in sim_pdb_dir.
constV=0
CRYST="CRYST1   27.240   31.870   34.230  88.52 108.53 111.89 P 1           1 "			                
#
#
#                *** DO NOT SET ANYTHING BELOW THIS LINE ***
#======================================================================#
#======================================================================#


function MakeMap {
	####################################################################
	#                                                                  #
	# CALCULATE ED MAP                                                 #
	#                                                                  #
	####################################################################
	
	#INPUT VARIABLES
	local pdb2map=$1 				#input pdb file name
	local mapname=$2				#output map name


  if [ $constV == 1 ]; then
cat <<EOF > awk.in
BEGIN{B=$B; printf("$CRYST\n");}
  /^ATOM/{
     RES= substr(\$0, 18, 9);
     X = substr(\$0, 31, 8);
     Y = substr(\$0, 39, 8);
     Z = substr(\$0, 47, 8);
     Ee = \$NF;
     printf("ATOM %6d %2s   %9s    %8.3f%8.3f%8.3f  1.00%6.2f%12s\n", \
            ++n,Ee,RES,X,Y,Z,B,Ee);}
  END{print "END"}
EOF
	else
cat <<EOF > awk.in
BEGIN{B=$B;na=${MD_mult[0]}; nb=${MD_mult[1]}; nc=${MD_mult[2]};} 
  /^CRYST/{a=\$2/na;b=\$3/nb;c=\$4/nc;al=\$5;be=\$6;ga=\$7;sg=(substr(\$0,56,15));
	printf "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %15s\n",a,b,c,al,be,ga,sg} 
  /^ATOM/{
     RES= substr(\$0, 18, 9);
     X = substr(\$0, 31, 8);
     Y = substr(\$0, 39, 8);
     Z = substr(\$0, 47, 8);
     Ee = \$NF;
     printf("ATOM %6d %2s   %9s    %8.3f%8.3f%8.3f  1.00%6.2f%12s\n", \
            ++n,Ee,RES,X,Y,Z,B,Ee);}
  END{print "END"}
EOF
  fi
awk -f awk.in $pdb2map > sfallme.pdb


	# GET UNIT CELL PARAMETERS
	pcell=( `head -1 sfallme.pdb | awk '{print $2,$3,$4,$5,$6,$7}'` )
	echo cell is ${pcell[*]}
	
	# ADD PDB SCALE MATRIX INFORMATION
	pdbset xyzin sfallme.pdb xyzout new.pdb << EOF > /dev/null
SPACE "$SG"
CELL ${pcell[*]}
EOF
	mv new.pdb sfallme.pdb


	# CALCULATE MAP
	sfall xyzin sfallme.pdb mapout ${mapname} << EOF > /dev/null
mode atmmap
CELL ${pcell[*]}
SYMM $SG
FORMFAC NGAUSS 5
$GRID
EOF
	#CLEAN
	rm -f params.txt 
	rm -f new.xyz 
	rm -f sfallme.xyz 
	#~ rm -f sfallme.pdb 
}


function AlignMaps {
	####################################################################
	#																   #
	# ALIGN TWO ED MAPS BY CONVOLUTION IN RECIPROCAL SPACE			   #
	#																   #
	####################################################################

	# The optimal x y z shift will be stored in $shift environmental variable
	
	# input variables
	right_map=$1 						# reference map (usually experimental)
	wrong_map=$2						# map to be shifted
	tempfile="$(mktemp -u tempXXXX_)"   # tempfile prefix
	
	# calculate SF from experimental map
	sfall mapin $right_map hklout ${tempfile}right.mtz << EOF > /dev/null
mode sfcalc mapin
resolution $reso
EOF

	# calculate SF from simulation map
	sfall mapin $wrong_map hklout ${tempfile}wrong.mtz << EOF > /dev/null
mode sfcalc mapin
resolution $reso
EOF

	# use deconvolution to find optimal shift
	# calculate Fexp/Fsim which are the SF's of the correlation function
	rm -f ${tempfile}del.mtz
	sftools << EOF > /dev/null
read ${tempfile}right.mtz
read ${tempfile}wrong.mtz
set labels
Fright
PHIright
Fwrong
PHIwrong
calc ( COL Fq PHIdel ) = ( COL Fright PHIright ) ( COL Fwrong PHIwrong ) /
calc COL W = COL Fq
select COL Fq > 1
calc COL W = 1 COL Fq /
select all
calc F COL Fdel = COL W 0.5 **
write ${tempfile}del.mtz col Fdel PHIdel
y
stop
EOF

	#calculate correlation map from correlation SF's
	fft hklin ${tempfile}del.mtz mapout ${tempfile}del.map << EOF > ${tempfile}.log
labin F1=Fdel PHI=PHIdel
reso $reso
EOF

	# make sure that we define "sigma" for the unmasked map
	echo "scale sigma 1 0" |\
	mapmask mapin ${tempfile}del.map mapout ${tempfile}zscore.map  > /dev/null

	#find highest peak on correlation map
	peakmax mapin ${tempfile}zscore.map xyzout ${tempfile}peak.pdb << EOF > ${tempfile}.log
numpeaks 10
EOF

	#print results
	cat ${tempfile}peak.pdb | awk '/^ATOM/{\
	 print substr($0,31,8),substr($0,39,8),substr($0,47,8),"   ",substr($0,61)+0}' |\
	 awk 'NR==1{max=$4} $4>max/3' | cat > ${tempfile}best_shift.txt
	zscore=`awk '{print $4;exit}' ${tempfile}best_shift.txt`
	echo "z-score: $zscore"
	shift=`awk '{print $1,$2,$3;exit}' ${tempfile}best_shift.txt`

	#clean
	rm -f ${tempfile}best_shift.txt
	rm -f ${tempfile}zscore.map
	rm -f ${tempfile}peak.pdb
	rm -f ${tempfile}del.map
	rm -f ${tempfile}.log
	rm -f ${tempfile}right.mtz
	rm -f ${tempfile}wrong.mtz
	rm -f ${tempfile}del.mtz
}


function AverageDensity {
	####################################################################
	#																   #
	# CALCULATE AVERAGE DENSITY MAP AND ITS STRUCTURE FACTORS   	   #
	#																   #
	####################################################################
		
	# SUM ALL MAPS
	echo "Summing maps"
	rm -f sum.map
	for ((frame=1; frame<=$frames; frame++)); do
		frame_map=`ls -1 maps | sed $frame'q;d'`
		if [ ! -e sum.map ]; then
			cp maps/${frame_map} sum.map

		else
			mapmask mapin1 sum.map mapin2 maps/${frame_map} mapout new.map <<EOF | grep -e "Mean density" |head -1
maps add
EOF
			mv new.map sum.map
		fi
	done
	
	#NORMALIZE SUM MAP
	scale=`echo $frames ${MD_mult[*]} $symops | awk '{print $1*$2*$3*$4*$5}'`	
	mapmask mapin sum.map mapout md_avg.map <<EOF | grep -e "Mean density"
scale factor `echo "1/$scale" | bc -l` 0
EOF
	echo -e "=======================\n"
	echo "Scaling map by 1/$scale"
	echo | mapdump mapin md_avg.map | egrep dens
	
	# CALCULATE SF FROM AVERAGE DENISTY MAP
	sfall mapin md_avg.map hklout md_avg_1.mtz << EOF > /dev/null
MODE SFCALC MAPIN
RESOLUTION $reso
EOF
	# EDIT MTZ FILE TO INCLUDE SIGFP AND FOM COLUMNS
	sftools << EOF > /dev/null
read md_avg_1.mtz
set labels
FP
PHI
calc Q COL SIGFP = 0.1
calc W COL FOM = 0.9
write md_avg_2.mtz col FP SIGFP PHI FOM
y
stop
EOF
	# REMOVE THE BFACTOR THAT WAS ADDED FOR SF CALCULATION TO AVOID SINGULARITY
	cad hklin1 md_avg_2.mtz hklout md_avg.mtz << EOF >/dev/null
scale file 1 1 -$B
labin file 1 all
EOF
	# CLEAN
	rm -f md_avg_1.mtz 
	rm -f md_avg_2.mtz
	rm -f sfall.mtz
	rm -f sum.map
}


########################################################################
#																	   #
# 	MAIN															   #
#																	   #
########################################################################

# SET-UP
frames=`ls -1 $sim_pdb_dir | wc -l`								  #no. of frames
symops=`awk -v SG=$SG '$4==SG{print $2;exit}' ${CLIBD}/symop.lib` #no. of sym operations
shift=" 0 0 0 "															
rm -rf shifts.txt maps
mkdir -p maps


# CALCULATE EXPERIMENTAL MAP
if [ $align == 1 ]; then
	echo "Calculating experimental map."
	MakeMap $exp_pdb right.map
	echo -e "=======================\n"
fi

# FOR EACH TRAJECTORY FRAME
start_time=`echo "puts [clock clicks -milliseconds]" | tclsh`
for ((frame=1; frame<=$frames; frame++)); do
	sim_pdb=`ls -1 $sim_pdb_dir | sed $frame'q;d'`
	
	# CALCULATE OPTIMAL SHIFT BY FOURIER CONVOLUTION
	if [ $align == 1 ]; then
		echo "Aligning ${sim_pdb_dir}${sim_pdb}"
		MakeMap ${sim_pdb_dir}${sim_pdb} wrong.map
		AlignMaps right.map wrong.map
		echo "$sim_pdb $shift" | awk '{printf("%-15s  %10s  %10s  %10s\n",$1,$2,$3,$4) }' | tee -a shifts.txt
	fi
	
	# CALCULATE MAP WITH SHIFT ($shift is now set to the optimal value)
	echo "Calculating map: ${sim_pdb_dir}${sim_pdb}"
	MakeMap ${sim_pdb_dir}${sim_pdb} maps/${sim_pdb%.*}.map
	
	# REPORT RUNNING TIME
	now_time=`echo "puts [clock clicks -milliseconds]" | tclsh`
	sofar=`echo $now_time $start_time | awk '{print ($1-$2)/1000}'`
	togo=`echo $sofar $frame $frames | awk '{print int($1/$2*$3-$1)}'`
	finish=`echo "puts [clock format [expr [clock seconds] + $togo]]" | tclsh`
	echo "Expect to finish at: $finish"
	echo -e "=======================\n"
done 

# CALCULATE SIMULATION AVERAGE ELECTRON DENSITY
AverageDensity





