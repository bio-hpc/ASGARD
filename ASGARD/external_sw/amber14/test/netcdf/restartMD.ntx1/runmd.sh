#!/bin/bash
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

. ../MasterTest.sh

TOP=../AlaDipeptide.TIP3P.parm7
CleanFiles mdinfo md.out md.nc md?.rst7 md.in ene? cpptraj.out md.crd

# txtinpcrd.rst7 and ncinpcrd.rst7 are a text and netcdf file respectively
# generated from identical MD runs. Run 1 MD using the text restart, and
# another MD using the netcdf restart, reading only coordinates (not 
# velocities). The resulting restarts should match within a certain precision.

echo "  Netcdf MD restart read test, ntx=1"

CRD=../txtinpcrd.rst7
cat > md.in <<EOF
Ala2 explicit solvent MD using text restart coords only
&cntrl
   imin = 0, nstlim = 100, dt=0.001, 
   ntx = 1, irest = 0, ig = 1, ntxo = 1, 
   ntwx = 50, ioutfm = 1, ntpr = 10, ntwr = 50000, ntwe = 10,
   iwrap = 0, nscm = 0,
   ntc = 2, ntf = 2, ntb = 1, cut = 8.0, 
   ntt = 1, tautp = 5.0, temp0 = 300.0, tempi = 300.0,
   ntp = 0, taup = 5.0,
&end
&ewald
   nfft1 = 30, nfft2 = 30, nfft3 = 30
&end
EOF
$DO_PARALLEL $TESTsander -O -i md.in -c $CRD -p $TOP -o md.out -x md.nc -r md1.rst7 -e ene1
CheckError $?

CRD=../ncinpcrd.rst7
cat > md.in <<EOF
Ala2 explicit solvent MD using netcdf restart coords only
&cntrl
   imin = 0, nstlim = 100, dt=0.001, 
   ntx = 1, irest = 0, ig = 1, ntxo = 1, 
   ntwx = 50, ioutfm = 1, ntpr = 10, ntwr = 50000, ntwe = 10,
   iwrap = 0, nscm = 0,
   ntc = 2, ntf = 2, ntb = 1, cut = 8.0, 
   ntt = 1, tautp = 5.0, temp0 = 300.0, tempi = 300.0,
   ntp = 0, taup = 5.0,
&end
&ewald
   nfft1 = 30, nfft2 = 30, nfft3 = 30
&end
EOF
$DO_PARALLEL $TESTsander -O -i md.in -c $CRD -p $TOP -o md.out -x md.nc -r md2.rst7 -e ene2
CheckError $?

../../dacdif -a 0.00001 md1.rst7 md2.rst7 
../../dacdif -r 0.00001 ene1 ene2
# Traj check if cpptraj present
SetCpptraj
$CPPTRAJ -p $TOP -y md.nc -x md.crd > cpptraj.out
CheckError $?
../../dacdif md.crd.save md.crd
