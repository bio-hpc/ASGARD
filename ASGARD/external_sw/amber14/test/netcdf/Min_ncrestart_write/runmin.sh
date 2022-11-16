#!/bin/bash
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

. ../MasterTest.sh

TOP=../AlaDipeptide.TIP3P.parm7
CRD=../inpcrd
CleanFiles min.out min.nc ncmin.rst7 mdinfo min.in textmin.rst7

# Test writing of netcdf restart during minimzation

cat > min.in <<EOF
Ala2 explicit solvent Min with netcdf restart
 &cntrl
   imin = 1, ntmin = 2, maxcyc = 100,
   ntx = 1, irest = 0, ig = 1, ntxo = 2,
   ntwx = 50, ioutfm = 1, ntpr = 10, ntwr = 50000, 
   ntc = 2, ntf = 2, ntb = 1, cut = 8.0, 
 &end
 &ewald
   nfft1=30, nfft2=30, nfft3=30,
 &end
EOF
echo "  Netcdf Minimization Restart Write Test"
$DO_PARALLEL $TESTsander -O -i min.in -c $CRD -p $TOP -o min.out -x min.nc -r ncmin.rst7
CheckError $?

# Second run to convert to text restart - comparing binary files is unreliable
cat > min.in <<EOF
Ala2 explicit solvent Min with netcdf restart
 &cntrl
   imin = 1, ntmin = 2, maxcyc = 0,
   ntx = 1, irest = 0, ig = 1, ntxo = 1,
   ntwx = 50, ioutfm = 1, ntpr = 10, ntwr = 50000, 
   ntc = 2, ntf = 2, ntb = 1, cut = 8.0, 
 &end
 &ewald
   nfft1=30, nfft2=30, nfft3=30,
 &end
EOF
$DO_PARALLEL $TESTsander -O -i min.in -c ncmin.rst7 -p $TOP -o min.out -x min.nc -r textmin.rst7
CheckError $?

../../dacdif -a 0.00001 textmin.rst7.save textmin.rst7 

