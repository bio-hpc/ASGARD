#!/bin/bash
 if [ "$#" -ne 1 ]; then
     echo "USAGE: $0 <sander>"
     echo
     echo "Checks if the <sander> executable does not support 3D-RISM"
     echo "calculations. If it does not, exits with a return value of 1 and"
     echo "an error message. <sander> maybe sander or sander.MPI. Other "
     echo "executables may give unexpected results, such as passing this test."
fi
touch foo
TESTsander=`which $1`
if [[ ! -x "$TESTsander" ]]; then
    echo "$1 is not an executable or does not exist."
    exit 1
fi
HAS_RISM=`$TESTsander -O -xvv foo 2> /dev/null | grep flag`
/bin/rm -f foo mdout mdin
if [ -n "$HAS_RISM" ] ; then
    echo "$TESTsander compiled without RISM support."
    exit 1
fi
exit 0
