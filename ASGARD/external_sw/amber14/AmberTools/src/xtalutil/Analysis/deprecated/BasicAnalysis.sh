#! /bin/sh
#set -x

# Super-brief instructions (DAC, 2/13): Copy this file to where your crystal
#      trajectory is.  In the copy, edit the variables below (per PAJ's 
#      instructions, also below), then run this script.  Outputs will contain 
#      a basic analysis of structures and B-factors.
#
#
# PAJ instructions: Feb.2013:
#
#  1. Prepare supercell trajectory for analysis. This does not have to be
#     centered in in any way. Such preparation could involve, for example,
#     stripping unnecessary atoms or filtering for a sparser selection of
#     frames to speed up analysis. The script StripMerge.sh in the auxillary
#     folder does this.
#
#  2. Prepare supercell topology and restart coordinate file. The script
#     prepscript_nowat.sh in the auxillary directory does this.
#
#  3. Prepare the topology, restart coordiantes and pdb file for the
#     single asymmetric unit.  This coordinates must be the same as the
#     original crystallographer's pdb. The pdb file should also contain
#     SMTRY/CRYST1 records. The script prepscript_asu.sh in the auxillary
#     folder does this.
#
#  4. Prepare files with the experimental B-factors. The analysis scripts
#     will calculate two sets of B-factors (by default named "calpha" and
#     "sdch" but you can specify the atom mask to whatever). These files
#     should contain the matching experimental B-factors used for the plots.
#     [dac note: need to specify the format for these files!]
#
#  5. All variables are now set in this master script. Copy this script to
#     the working directory (where the files in steps 1-4 reside), edit the
#     variables below, then run this script.
#
#  6. THINGS TO DO:
#     - add verbose flag to supress cpptraj output, erase intermediate
#     - scripts, etc. separate plotting from the analysis/data output.
#     - some variables that are set by hand now could be automated
#     - (unitcell, asus, residues)


# List of scripts used by this (master) script:
#       plotbfac.py
#       plotrmsd_permonomer.py
#       plotrmsd.py
#       plotvolume.py
#       IndivUCrmsd.py
#       RevSymm_netcdf.py
#       SplitTrajectory.py
#       distance_matrix.py
#       analyse_revsym.py


#======================================================================#
#                                                                      #
#  			SET ALL VARIABLES HERE                                     #
#                                                                      #

# location of the analysis scripts:
XTAL_ANALYSIS_PATH=$AMBERHOME/AmberTools/src/xtalutil/Analysis

# working directory variable, for use below:
WD=/Volumes/godel1/collagen/1cgd

#  Amber files for the supercell (might have waters, etc stripped):
SC_PRMTOP=$WD/supercell.parm7        # supercell topology
SC_RST7=$WD/supercell.rst7           # supercell with pdb (experimental) coords.
SC_TRAJECTORY=$WD/raw/md20-66.prot.nc       # trajectory to be analyzed

#  Amber files for the asymmetric unit:
ASU_RST7=$WD/asu.rst7                  # same coordinates as pdb!
ASU_PRMTOP=$WD/asu.parm7 
ASU_PDB=$WD/asu.pdb                    # must contain SMTRY and CRYST1 records 
                                       # from original pdb

#   Information on how the supercell was constructed
PROP="1 3 2"          # Propagation operations used to build supercell (a, b, c)
ASUS=4                # Number of ASU's in unit cell
UNITCELLS=24          # Number of unit cells in the supercell
RESIDUES=90           # Number of residues in the (stripped) ASU
PLOT_PFX="1cgd"       # filename prefix (for plots, etc.)

#                                                                      #
# 		SETTINGS FOR "ANALYSE_REVSYM.PY"                               #
#                                                                      # 
# This script will calculate two sets of B-factors, one which will be called
# by default "calpha" and one called "sidechain" calculated byres. Two sets
# of rmsd values will also be calculated, one by default called "heavy" and 
# one called "backbone".

BM1=":1-90@CA"                 # B factor Calpha mask
BM2=":1-90&!(@H=,CA,C,O,N)"    # B factor sidechain mask
RM1=":1-90&!(@H=)"             # RMSD heavy atom mask (also used for 
                               # RotFit B-factor calculation
RM2=":1-90@CA,C,N"             # RMSD backbone atom mask
Timestep=0.1   	               # Timestep for rmsd output
BFAC_CRYST_CALPHA=             # experimental Bfacs
BFAC_CRYST_SDCH=

#                                                                      #
#     DON'T SET ANYTHING BELOW THIS LINE (HOPEFULLY NOT NECESSARY)     #
#======================================================================#



########################################################################
cd ${WD}

# translate trajectory so centers of mass aligned with pdb crystal
echo
echo '##########################'
echo '# translating trajectory #'
echo '##########################'
cat <<EOF > ctraj.translate.in
parm ${SC_PRMTOP}
trajin ${SC_TRAJECTORY}
trajout fit.nc netcdf
reference ${SC_RST7}
rmsd reference '!@H=' norotate nofit out drift_nofit.dat
rmsd reference '!@H=' norotate out drift.dat
go
EOF
cpptraj < ctraj.translate.in
if [ $? -ne 0 ]; then
	exit
fi
rm ctraj.translate.in

#~ ./FullAnalysisScripts/check_mergetraj.py ${SC_PRMTOP} fit.nc ${SC_RST7} 
 

echo
echo '########################'
echo '# splitting trajectory #'
echo '########################'
rm -rf splittrajectories; mkdir splittrajectories; cd splittrajectories
${XTAL_ANALYSIS_PATH}/SplitTrajectory.py \
	-p ${SC_PRMTOP} -t ../fit.nc \
	-u ${UNITCELLS} -a ${ASUS} -r ${RESIDUES} 
cd ..

echo
echo '##########################################'
echo '# reverse symmetry on split trajectories #'
echo '##########################################'
rm -rf revsym; mkdir revsym; cd revsym
${XTAL_ANALYSIS_PATH}/RevSym_netcdf.py \
	-p ${ASU_PDB} -r ${SC_RST7} \
	-prop "$PROP"

echo
echo '#########################################'
echo '# analyzing reversed split trajectories #'
echo '#########################################'
${XTAL_ANALYSIS_PATH}/analyse_revsym.py  \
	-p ${ASU_PRMTOP} -pdb ${ASU_PDB}  \
	-u ${UNITCELLS} -a ${ASUS}  \
	-Bm1 ${BM1}  \
	-Bm2 ${BM2}  \
	-Rm1 ${RM1}  \
	-Rm2 ${RM2}  \
	-tm ${Timestep}
cd ..

# volume
${XTAL_ANALYSIS_PATH}/plotvolume.py \
	-t fit.nc \
	-r ${SC_RST7} \
	-Title $PLOT_PFX

	
# plot
echo
echo '############'
echo '# plotting #'
echo '############'
${XTAL_ANALYSIS_PATH}/plotrmsd.py $PLOT_PFX
${XTAL_ANALYSIS_PATH}/plotrmsd_v4.py $PLOT_PFX
${XTAL_ANALYSIS_PATH}/plotbfac.py $PLOT_PFX $BFAC_CRYST_CALPHA $BFAC_CRYST_SDCH
${XTAL_ANALYSIS_PATH}/plotrmsd_permonomer.py $PLOT_PFX $ASUS $UNITCELLS


# Avg Monomer RMSD
cd splittrajectories
${XTAL_ANALYSIS_PATH}/IndivUCrmsd.py \
	-p ${ASU_PRMTOP} -pdb ${ASU_PDB} -b $BFAC_CRYST_CALPHA \
	-u ${UNITCELLS} -a ${ASUS}  \
	-Bm1 ${BM1}  \
	-Rm1 ${RM1}  \
	-Title $PLOT_PFX
cd ..

 
echo
echo '#################################'
echo '# calculating distance matrices #'
echo '#################################'
rm -rf distance_matrix
${XTAL_ANALYSIS_PATH}/distance_matrix.py \
	-p ${SC_PRMTOP} -t ../fit.nc \
	-u ${UNITCELLS} -a ${ASUS}  \
	-tm ${Timestep}

# move plots to plots directory
rm -rf plots; mkdir plots; mv *.png volume.dat volume.txt plots
#~ rm -rf splittrajectories

