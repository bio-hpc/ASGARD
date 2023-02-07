#!/bin/csh -f
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
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

touch dummy
set MDIN = mdin.znb.md
set MDOUT = znme2.znb.md.out
$DO_PARALLEL $sander -O -i $MDIN -o $MDOUT < dummy || goto error
../../dacdif -t 1 $MDOUT.save $MDOUT

/bin/rm -f mdin mdinfo mdcrd dummy restrt
exit(0)

error:
echo "  ${0}:  Program error"
exit(1)
