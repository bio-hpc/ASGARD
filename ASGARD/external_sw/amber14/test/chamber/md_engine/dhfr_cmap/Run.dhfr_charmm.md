#!/bin/csh -f
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

if( ! $?TESTsander ) set TESTsander = "../../../../bin/sander"

if( ! $?DO_PARALLEL ) then
        setenv DO_PARALLEL " "
endif

cat > mdin <<EOF
 short md
 &cntrl
   ntx=1, irest=0,
   imin=0,nstlim=50,
   dt=0.002,ntc=2,ntf=2,
   ntt=1,tempi=300.0,temp0=300.0, 
   ntpr=1,igb=1,cut=9999.0,ntwx=0,
   ntwr=0,ntwe=0,ntb=0,
 /
EOF

set output = mdout.dhfr_charmm

touch dummy
$DO_PARALLEL $TESTsander -O -i mdin -c inpcrd -o $output < dummy || goto error

if ( `basename $TESTsander` == "pmemd.mic_offload.MPI" ) then
../../../dacdif -r 1.0e-07 $output.save $output
else
../../../dacdif $output.save $output
endif
/bin/rm -f mdin restrt mdinfo dummy
exit(0)

error:
echo "  ${0}:  Program error"
exit(1)
