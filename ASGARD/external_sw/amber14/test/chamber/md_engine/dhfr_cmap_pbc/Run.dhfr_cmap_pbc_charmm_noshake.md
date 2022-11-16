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
   dt=0.001,ntc=1,ntf=1,
   ntt=1,tempi=300.0,temp0=300.0, 
   ntpr=1,ntb=1,ntp=0,cut=9.0,ntwx=0,
   ntwr=0,ntwe=0,
 /
 &ewald
  ew_coeff=0.340,nfft1=96,nfft2=80,nfft3=64,order=4,vdwmeth=0,
 /
EOF

set output = mdout.dhfr_charmm_noshake.md

touch dummy
$DO_PARALLEL $TESTsander -O -i mdin -c inpcrd -o $output < dummy || goto error

if ( `basename $TESTsander` == "pmemd.mic_offload.MPI" ) then
../../../dacdif -r 1.0e-07 -t 1 $output.save $output
else
../../../dacdif -t 1 $output.save $output
endif
/bin/rm -f mdin restrt mdinfo dummy
exit(0)

error:
echo "  ${0}:  Program error"
exit(1)
