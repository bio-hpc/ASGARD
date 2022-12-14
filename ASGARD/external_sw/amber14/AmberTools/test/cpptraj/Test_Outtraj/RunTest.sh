#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles rms.in rmsd.dat test.crd T1.crd T2.crd maxmin.in maxmin.crd filter.dat

CheckNetcdf
TOP="../tz2.truncoct.parm7"

# Test 1
cat > rms.in <<EOF
noprogress
trajin ../tz2.truncoct.nc
outtraj test.crd title "trajectory generated by ptraj"
rms first :2-11 out rmsd.dat
outtraj T1.crd
trajout T2.crd
EOF
INPUT="rms.in"
RunCpptraj "RMSD Test with outtraj."
DoTest rmsd.dat.save rmsd.dat
DoTest T1.crd T2.crd
DoTest ../tz2.truncoct.crd test.crd

# Test 2
cat > maxmin.in <<EOF
noprogress
trajin ../tz2.truncoct.nc
rms R1 first :2-11
filter R1 min 0.7 max 0.8 out filter.dat
outtraj maxmin.crd 
EOF
INPUT="maxmin.in"
RunCpptraj "Outtraj Test with filtering."
DoTest maxmin.crd.save maxmin.crd

# Test 3
cat > maxmin.in <<EOF
trajin ../tz2.truncoct.nc
rms R1 first :2-11
outtraj maxmin.crd maxmin R1 min 0.7 max 0.8
EOF
INPUT="maxmin.in"
RunCpptraj "Outtraj Test with maxmin."
DoTest maxmin.crd.save maxmin.crd

CheckTest

EndTest

exit 0
