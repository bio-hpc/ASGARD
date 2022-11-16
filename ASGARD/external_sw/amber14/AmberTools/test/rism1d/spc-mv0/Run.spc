#!/bin/csh -f


cat > spc.inp <<EOF
&PARAMETERS
outlist='uxghcbtensq', THEORY='DRISM', closure='MV0',
!grid
NR=16384, DR=0.025, rout=100, kout=30,
!MDIIS
mdiis_nvec=20, mdiis_del=0.3, tolerance=1.e-12,
!iter
ksave=0, progress=1, maxstep=10000,
!ElStat
SMEAR=1, ADBCOR=0.5,
!bulk solvent properties
temperature=298, DIEps=78.497,
NSP=1
/
&SPECIES
!DENSITY=55.5,
!corresponds very closely to 0.0333 1/A3
DENSITY=55.296d0,
MODEL="../../../../dat/rism1d/model/SPC.mdl"
/
EOF

../../../bin/rism1d spc > spc.out || goto error

../check1d spc
/bin/rm -f spc.inp spc.out spc.sav

exit(0)

error:
echo "  ${0}:  Program error"
exit(1)
