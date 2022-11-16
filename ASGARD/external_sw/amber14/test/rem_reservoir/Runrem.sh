#!/bin/bash
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, advanced
#TEST-STATE   not finished

# -- Files --
# cluster.info: RREMD clustering info
# rem.log.save, rst.1.02.save, mdout.1.03.save: Test save files
# rst.heat.0?: Input coordinates
# saveene: RREMD reservoir structure energies
# full.topo: Topology
# reserv/frame.XXXXXX: RREMD structure reservoir
# Runrem.sh: This file

. ../program_error.sh
sander=../../bin/sander.MPI
if [[ ! -z $TESTsander ]] ; then
  sander=$TESTsander
fi
if [[ -z $DO_PARALLEL ]] ; then
  echo "REM requires a parallel environment!"
  exit 0
fi
numprocs=`$DO_PARALLEL ../numprocs`
if [[ `basename $sander` = 'sander.MPI' ]] ; then
   if [ $numprocs -ne 4 -a $numprocs -ne 8 -a $numprocs -ne 12 -a $numprocs -ne 16 ] ; then
      echo "This test case requires 4, 8, 12, or 16 MPI threads!"
      echo "Only using $numprocs"
      exit 0
   fi
else
   echo "Reservoir Replica Exchange not yet supported."
   exit 0
   #if [ $numprocs -ne 8 -a $numprocs -ne 12 -a $numprocs -ne 16 ]; then
   #   echo "This test case requires 8, 12, or 16 MPI threads!"
   #   exit 0
   #fi
fi

# Generate groupfile and input
if [[ -e groupfile ]] ; then
  rm groupfile
fi
N=1
for TEMP in "300.00" "309.83" "319.91" "330.22" ; do
  EXT=`printf "%02i" $N`
  cat > mdin.$EXT <<EOF
Production 10 ps, NTV, Lang.,
 &cntrl
    ntpr=5000, ntwr=5000, ntwx=500, ntxo=2,
    ntf=2, ntc=2, ntp=0,
    ntt=3, gamma_ln=10, ig=785558, cut = 8.0,
    nstlim=100, numexchg=50, dt=0.002,
    imin=0, ntx=5, irest=1, ioutfm=1,
    temp0=$TEMP,
 /
EOF
  echo "-O -O -rem 1 -remlog rem.log -remtype rem.type -rremd 3 -saveene saveene -clusterinfo cluster.info -reservoir reserv/frame -i mdin.$EXT -p full.topo -c rst.heat.$EXT -r ncrst.$EXT -o mdout.$EXT -x traj.$EXT -inf mdinfo.$EXT" >> groupfile
  ((N++))
done

$DO_PARALLEL $sander -O -ng 4 -groupfile groupfile || error

../dacdif rem.log.save rem.log
../dacdif rem.type.save rem.type
../dacdif mdout.04.save mdout.04
# Remove files
/bin/rm mdout.01 mdout.02 mdout.03 mdin.?? mdinfo.?? ncrst.?? traj.??

exit 0
