#! /bin/bash
set -x

#----------------------------------------------------------------------#
#                                                                      #
# A preparation script to merge and center of mass fit all frames from #
# a trajectory. Stripping is not performed.                            #
#----------------------------------------------------------------------#

#SET VARIABLES
WORKING_DIR=/home/pjanowsk/Case/4lzt/RunSi/average_density
TRAJ_ROOT=/home/case/xtal/4lzt_xtal/mdi
TRAJ_EXT=nc
TRAJ_START=1
TRAJ_END=75
OFFSET=10
SC_TOPO=/home/pjanowsk/Case/4lzt/RunSi/4lztSi.prmtop
SC_RST7=/home/pjanowsk/Case/4lzt/RunSi/4lztSi.rst7


########################################################################

cd ${WORKING_DIR}
rm -rf ctraj.merge.in merge_wat.nc

# mergetrajectory
echo
echo '##########################'
echo '# merging trajectory     #'
echo '##########################'
echo "parm ${SC_TOPO}" >> ctraj.merge.in
echo "trajin ${SC_RST7}" >> ctraj.merge.in
for i in `seq ${TRAJ_START} ${TRAJ_END}`; do
	echo "trajin ${TRAJ_ROOT}${i}.${TRAJ_EXT} 1 -1 ${OFFSET}" >> ctraj.merge.in
done
echo "center mass origin" >> ctraj.merge.in
echo "trajout merge_wat.nc netcdf" >> ctraj.merge.in
cpptraj <ctraj.merge.in
if [ $? -ne 0 ]; then
	exit
fi

# translate trajectory so centers of mass aligned with pdb crystal
echo
echo '##########################'
echo '# translating trajectory #'
echo '##########################'
cat <<EOF > ctraj.translate.in
parm ${SC_TOPO}
trajin merge_wat.nc
reference ../4lztSi_centonpdb.rst7
rmsd reference '!:WAT&!@H=' norotate
trajout fit_wat.nc netcdf
go
EOF
cpptraj < ctraj.translate.in
if [ $? -ne 0 ]; then
	exit
fi
rm ctraj.translate.in ctraj.merge.in

./check_mergetraj.py ${SC_TOPO} fit_wat.nc ${SC_RST7} 

# Now run MakePdb4Map.py


