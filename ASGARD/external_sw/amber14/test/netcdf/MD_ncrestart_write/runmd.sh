#!/bin/bash
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

. ../MasterTest.sh

TOP=../AlaDipeptide.TIP3P.parm7
CRD=../inpcrd
CleanFiles md.in md.out md.nc ncmd.rst7 mdinfo textmd.rst7

# Test writing of netcdf restart format

cat > md.in <<EOF
Ala2 explicit solvent MD with netcdf restart
&cntrl
   imin = 0, nstlim = 100, dt=0.001, 
   ntx = 5, irest = 1, ig = 1, ntxo = 2, 
   ntwx = 50, ioutfm = 1, ntpr = 10, ntwr = 50000,
   iwrap = 0, nscm = 0,
   ntc = 2, ntf = 2, ntb = 1, cut = 8.0, 
   ntt = 1, tautp = 5.0, temp0 = 300.0, tempi = 300.0,
   ntp = 0, taup = 5.0,
&end
&ewald
   nfft1 = 30, nfft2 = 30, nfft3 = 30
&end
EOF

echo "  Netcdf MD Restart Write Test"
$DO_PARALLEL $TESTsander -O -i md.in -c $CRD -p $TOP -o md.out -x md.nc -r ncmd.rst7
CheckError $?

# Second run to convert to text restart - comparing binary files is unreliable
cat > md.in <<EOF
Ala2 explicit solvent MD with netcdf restart
&cntrl
   imin = 0, nstlim = 1, dt=0.001, 
   ntx = 5, irest = 1, ig = 1, ntxo = 1, 
   ntwx = 50, ioutfm = 1, ntpr = 10, ntwr = 50000,
   iwrap = 0, nscm = 0,
   ntc = 2, ntf = 2, ntb = 1, cut = 8.0, 
   ntt = 1, tautp = 5.0, temp0 = 300.0, tempi = 300.0,
   ntp = 0, taup = 5.0,
&end
&ewald
   nfft1 = 30, nfft2 = 30, nfft3 = 30
&end
EOF
$DO_PARALLEL $TESTsander -O -i md.in -c ncmd.rst7 -p $TOP -o md.out -x md.nc -r textmd.rst7
CheckError $?

../../dacdif -a 0.0001 textmd.rst7.save textmd.rst7
