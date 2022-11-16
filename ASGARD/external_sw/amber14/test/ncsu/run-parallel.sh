#!/bin/sh
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

if test -z "${AMBERHOME}"; then
    export AMBERHOME="`pwd`/../../"
fi

if test -z "${TESTsander}"; then
    export TESTsander="${AMBERHOME}/bin/sander.MPI"
fi

TESTS_2="abmd_ANALYSIS abmd_FLOODING abmd_UMBRELLA smd pmd smd2"

# 4*n mpi threads
TESTS_4="bbmd abremd mwabmd premd"

TESTS="${TESTS_2} ${TESTS_4}"

for t in ${TESTS} ; do
    echo ">>>>>>> doing '$t'"
    (cd $t ; ./run.sh ; cd ..)
done
