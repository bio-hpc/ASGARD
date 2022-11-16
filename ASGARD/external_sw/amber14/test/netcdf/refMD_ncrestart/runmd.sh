#!/bin/bash
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

. ../MasterTest.sh

TOP=../AlaDipeptide.TIP3P.parm7
CRD=../inpcrd
REF=ncrefc.rst7
CleanFiles md.out md.nc md.rst7 mdinfo md.in ene.dat

# Test use of netcdf restart file as reference coordinates

cat > md.in <<EOF
Ala2 explicit solvent MD with netcdf restart ref coords
&cntrl
   imin = 0, nstlim = 100, dt=0.001, 
   ntx = 5, irest = 1, ig = 1, ntxo = 1, 
   ntwx = 50, ioutfm = 1, ntpr = 10, ntwr = 50000, ntwe = 10, 
   iwrap = 0, nscm = 0,
   ntc = 2, ntf = 2, ntb = 1, cut = 8.0, 
   ntt = 1, tautp = 5.0, temp0 = 300.0, tempi = 300.0,
   ntp = 0, taup = 5.0,
   ntr = 1, restraintmask = ':1-3@C,CA,N', restraint_wt = 5.0, 
&end
&ewald
   nfft1 = 30, nfft2 = 30, nfft3 = 30
&end
EOF

echo "  Restrained MD with netcdf Restart Reference Coords Test"
$DO_PARALLEL $TESTsander -O -i md.in -c $CRD -ref $REF -p $TOP -o md.out -x md.nc -r md.rst7 -e ene.dat
CheckError $?

../../dacdif -a 0.00001 md.rst7.save md.rst7
if [ `basename $TESTsander` = "pmemd.mic_offload.MPI" ]; then
../../dacdif -r 1.0e-06 ene.dat.save ene.dat
else
../../dacdif ene.dat.save ene.dat
fi
