#!/bin/csh -f
#TEST-PROGRAM sander
#TEST-DESCRIP test H, C and O parameters of PM3/MMX2
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

set sander = "../../../bin/sander"
if( $?TESTsander ) then
    set sander = $TESTsander
endif

if( ! $?DO_PARALLEL ) then
    setenv DO_PARALLEL " "
else
    echo "This test not set up for parallel"
    echo "need #nres>#nproc"
    exit 0
endif

cat > mdin <<EOF
formic acid min PM3-MMX2
&cntrl
 imin=1, maxcyc=20, ntmin=2, ncyc=10,
 ntwr=20, ntpr=1, ntb=0, igb=0,
 cut=999.0, ifqnt=1,
/
&qmmm
 iqmatoms=1,2,3,4,5,
 qm_theory='PM3',
 qmmm_int=4,
 verbosity=0,
/

EOF

set output = formic_acid_dimer.out

touch dummy
$DO_PARALLEL $sander -O -p formic_acid_dimer.ff99.top -c formic_acid_dimer.crd -o $output < dummy || goto error
../../dacdif -t 3 $output.save $output

/bin/rm -f restrt mdin mdinfo mdcrd dummy
exit(0)

error:
echo "  ${0}:  Program error"
exit(1)
