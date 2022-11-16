#! /bin/bash
set -x
# 1. Run prepscript_nowat.sh to make nowat files first.
# 2. Also prepare unitcell "UC" files (cut out all atoms except 1 asym unit,
#    use the tleap script from prepscript_nowat.sh, modify it to make the 
#    prmtop and rst7. Change box dimensions.
# 3. Center UC.rst7 on pdb (translate with ptraj)
# 4. Modify arguments for the scripts: 
#		SplitTrajectory.py
#		RevSymm_netcdf.py
#		analyse_revsym.py
#		IndivUCrmsd.py
#       distance_matrix.py

# 5. List of scripts in FullAnalysisScripts: 
#       check_mergetraj.py
#       difference_vector.py
#       plotbfac.py
#       plotrmsd_permonomer.py
#       plotrmsd.py
#       plotvolume.py
#       process_mdout.perl
#       IndivUCrmsd.py
#       RevSymm_netcdf.py
#       SplitTrajectory.py
#       distance_matrix.py
#       analyse_revsym.py


#SET VARIABLES
WORKING_DIR=/home/pjanowsk/Case/4lzt/RunSh_1300ns/
SC_TOPO=4lztSh.prmtop 
SC_RST7=4lztSc_centonpdb_nowat.rst7  #supercell with pdb (experimental) coordinates for each asym unit and pdb box size
NOWAT_TOPO=4lztSc_nowat.prmtop
UC_RST7=UC_centonpdb.rst7 #centered on pdb!
UC_TOPO=UC.prmtop
#BFAC_CRYST=bfac_crystal.dat
OUT_ROOT=/home/case/xtal/4lzt_xtal/mdh
TRAJ_ROOT=/home/case/xtal/4lzt_xtal/mdh
TRAJ_EXT=nc
TRAJ_START=1
TRAJ_END=68
STRIP_MASK=":1669-9999"
CENT_ATOM=27
PLOT_NAME='4lzt_Sh'
ASYMUNITS=1
UNITCELLS=12


########################################################################

cd ${WORKING_DIR}

# strip waters and mergetrajectory
echo -e '\n############################\nmerging trajectory\n####################'
rm -rf ptraj_mergetrajectory mergetraj_cent_nowat.nc netcdf
#~ echo "trajin ${SC_RST7}" >> ptraj_mergetrajectory
for i in `seq ${TRAJ_START} ${TRAJ_END}`; do
	echo "trajin ${TRAJ_ROOT}${i}.${TRAJ_EXT}" >> ptraj_mergetrajectory
done
echo "strip ${STRIP_MASK}" >> ptraj_mergetrajectory
#~ echo "center mass origin" >> ptraj_mergetrajectory
echo "trajout mergetraj_nowat.nc netcdf" >> ptraj_mergetrajectory
cpptraj ${SC_TOPO} <ptraj_mergetrajectory

############################
# This is a 4lzt thing: because I couldn't add all the necessary waters at once,
# I added a part first, creating 4lztSc, equilibrated, than added more, to get
# 4lztSh. Thus, the pdb supercell is 4lztSc but I can't add this to the .nc files
# without first stripping waters because number of waters does not agree. So
# first I strip waters, than I add 4lztSc, than I center on pdb.
##########################
rm -rf ptraj_center
echo "trajin 4lztSc_centonpdb_nowat.rst7" >> ptraj_center
echo "trajin mergetraj_nowat.nc" >> ptraj_center
echo "center mass origin" >> ptraj_center
echo "trajout mergetraj_cent_nowat.nc netcdf" >> ptraj_center
cpptraj ${NOWAT_TOPO} <ptraj_center

# translate trajectory so centers of mass aligned with pdb crystal
echo -e '\n############################\ntranslating trajectory\n####################'
rm -rf ptraj_translate mergetraj_centonpdb_nowat.nc
./FullAnalysisScripts/difference_vector.py $UC_RST7 mergetraj_cent_nowat.nc $CENT_ATOM
tx=`cat tmp | awk '{print $1}'`
ty=`cat tmp | awk '{print $2}'`
tz=`cat tmp | awk '{print $3}'`
rm -rf tmp ptraj_translate
echo "trajin mergetraj_cent_nowat.nc" >>ptraj_translate
echo "translate x ${tx} y ${ty} z ${tz}" >>ptraj_translate
echo "trajout mergetraj_centonpdb_nowat.nc netcdf" >>ptraj_translate
ptraj ${NOWAT_TOPO} <ptraj_translate


# run some checks to make sure that the traj is centered on pdb
echo -e '\n############################\nchecking trajectory\n####################'
./FullAnalysisScripts/check_mergetraj.py ${NOWAT_TOPO} mergetraj_centonpdb_nowat.nc ${UC_RST7}
if [ $? -eq 1 ] ; then
	echo "ERROR in check_mergetraj.py."
	exit
	fi
if [ $? -eq 0 ] ; then
	echo "check_mergetraj.py passed!"
	fi


# split trajectories
echo -e '\n############################\nsplitting trajectory\n####################'
rm -rf splittrajectories
./FullAnalysisScripts/SplitTrajectory.py


# reverse symmetry trajectories and analyse
echo -e '\n############################\nreverse symmetry on split trajectories\n####################'
rm -rf revsym; mkdir revsym; cd revsym
../FullAnalysisScripts/RevSymm_netcdf.py


echo -e '\n############################\nanalyzing reversed split trajectories\n####################'
../FullAnalysisScripts/analyse_revsym.py
cd ..

# plot
echo -e '\n############################\nplotting\n####################'
./FullAnalysisScripts/plotrmsd.py $PLOT_NAME
./FullAnalysisScripts/plotrmsd_v4.py $PLOT_NAME
./FullAnalysisScripts/plotbfac.py $PLOT_NAME
./FullAnalysisScripts/plotrmsd_permonomer.py $PLOT_NAME $ASYMUNITS $UNITCELLS

# Avg Monomer RMSD

cd splittrajectories
../FullAnalysisScripts/IndivUCrmsd.py
cd ..

# volume
rm -rf volume
mkdir volume
cd volume
input=$(for i in `seq ${TRAJ_START} ${TRAJ_END}`; do echo -n "${OUT_ROOT}$i.out ";done)
../FullAnalysisScripts/process_mdout.perl $input
cd ..
./FullAnalysisScripts/plotvolume.py $SC_RST7 $PLOT_NAME

# distance matrix
echo -e '\n#############################\ncalculating distance matrices\n##############'
rm -rf distance_matrix
./FullAnalysisScripts/distance_matrix.py

# move plots to plots directory
rm -rf plots
mkdir plots
mv *png plots


