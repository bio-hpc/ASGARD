#!/bin/sh
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

. ../common.sh

#
# don't run if there is no netCDF
#

set_NC_UTILS

if test $? -ne 0; then
   exit 0;
fi

set_SANDER sander.MPI
set_SANDER_CMD 4

JUNK="bbmd.log bbmd.?.txt umbrella.?.nc mdout.? inpcrd.?"
JUNK="${JUNK} mt19937.nc restrt.? mdcrd.? groups prmtop sander.out"

#
# remove the junk
#

/bin/rm -rf ${JUNK} junk.*

#
# prepare files
#

cp -p ../prmtop ../inpcrd.? ../groups ../mt19937.nc .

#
# run SANDER
#

${SANDER_CMD} -ng 4 -groupfile groups > sander.out 2>&1

../../dacdif -t 1 save/bbmd.log bbmd.log

../../dacdif save/mt19937.nc mt19937.nc

for i in 1 2 3 4; do
    ../../dacdif -t 2 save/bbmd.$i.txt bbmd.$i.txt
    ../../dacdif -t 1 save/mdout.$i mdout.$i
done

for i in 1 2 3; do
   do_ncdump '%14.10f' save/umbrella.$i.nc umbrella.$i.save.ncdump
   do_ncdump '%14.10f' umbrella.$i.nc umbrella.$i.ncdump
   ../../dacdif umbrella.$i.save.ncdump umbrella.$i.ncdump
done

#
# preserve the junk on failure
#

save_junk_on_failure ${JUNK}

#
# remove the junk
#

/bin/rm -f ${JUNK} *.ncdump profile_mpi
