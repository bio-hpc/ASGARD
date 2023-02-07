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
set_SANDER_CMD 1

JUNK="mdout mdinfo restrt monitor.txt umbrella.nc inpcrd prmtop sander.out"

#
# remove the junk
#

/bin/rm -rf ${JUNK} junk.*

#
# prepare files
#

cp -p ../prmtop .
cp -p ../inpcrd.4 ./inpcrd

#
# run SANDER
#

${SANDER_CMD} > sander.out 2>&1

../../dacdif -t 1 save/mdout mdout
../../dacdif -t 1 save/monitor.txt monitor.txt

do_ncdump '%14.10f' save/umbrella.nc umbrella.save.ncdump
do_ncdump '%14.10f' umbrella.nc umbrella.ncdump
../../dacdif umbrella.save.ncdump umbrella.ncdump

#
# preserve the junk on failure
#

save_junk_on_failure ${JUNK}

#
# remove the junk
#

/bin/rm -f ${JUNK} *.ncdump profile_mpi
