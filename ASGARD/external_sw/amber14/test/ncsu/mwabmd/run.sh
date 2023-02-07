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

JUNK="umbrella.?.nc mdout.? restrt.? abmd.?.txt inpcrd.?"
JUNK="${JUNK} mdcrd.? groups prmtop sander.out *.log"

#
# remove the junk
#

/bin/rm -rf ${JUNK} junk.*

#
# prepare files
#

cp -p ../prmtop ../groups ../inpcrd.? .

#
# run SANDER
#

${SANDER_CMD} -ng 4 -rem 0 -groupfile groups > sander.out 2>&1

for i in 1 2 3 4; do
    ../../dacdif -t 1 save/abmd.$i.txt abmd.$i.txt
    ../../dacdif -t 1 save/mdout.$i mdout.$i
done

for i in 1 2 3 4; do
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
