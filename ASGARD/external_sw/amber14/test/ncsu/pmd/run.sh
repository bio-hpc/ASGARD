#!/bin/sh
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

. ../common.sh

set_SANDER sander.MPI
set_SANDER_CMD 1

JUNK="mdout mdinfo restrt pmd.txt inpcrd prmtop sander.out"

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
../../dacdif -t 1 save/pmd.txt pmd.txt

#
# preserve the junk on failure
#

save_junk_on_failure ${JUNK}

#
# remove the junk
#

/bin/rm -f ${JUNK} profile_mpi
