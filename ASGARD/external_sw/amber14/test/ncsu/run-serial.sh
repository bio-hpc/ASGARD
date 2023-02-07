#!/bin/sh
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

if test -z "${AMBERHOME}"; then
    export AMBERHOME="`pwd`/../../"
fi

if test -z "${TESTsander}"; then
    export TESTsander="${AMBERHOME}/bin/sander"
fi

TESTS="abmd_ANALYSIS abmd_FLOODING abmd_UMBRELLA smd pmd smd2"

for t in ${TESTS} ; do
    echo ">>>>>>> doing '$t'"
    (cd $t ; ./run.sh ; cd ..)
done
