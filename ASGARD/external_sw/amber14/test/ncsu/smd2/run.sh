#!/bin/sh
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

. ../common.sh

#
# Test to make sure there are at most six processors.
#
numprocs=`${DO_PARALLEL} echo "Testing number of processors" | wc -l`
if [ "$numprocs" -gt 6 ] ; then
	echo "This test requires six or fewer processors. Skipping."
	exit
fi


set_SANDER sander.MPI
set_SANDER_CMD 1

JUNK="mdout mdinfo restrt work.txt sander.out"

#
# remove the junk
#

/bin/rm -rf ${JUNK} junk.*

#
# run SANDER
#

${SANDER_CMD} > sander.out 2>&1

../../dacdif -t 1 save/mdout mdout
../../dacdif -t 1 save/work.txt work.txt

#
# preserve the junk on failure
#

save_junk_on_failure ${JUNK}

#
# remove the junk
#

/bin/rm -f ${JUNK} profile_mpi
