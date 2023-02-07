#!/bin/sh

#This script is for testing the cuda implmentation via the test
#Makefiles. It supports the use of 2 command line flags to control
#testing:
#
# $1 = PREC_MODEL  - default = SPFP (Mixed Precision)
#
#Defaults
# PREC_MODEL = SPFP
#
#Example use:
#
# ./test_amber_cuda.x         (This will run with defaults)
# ./test_amber_cuda.x DPFP (This selects all double precision as the 
#                              precision model)
#
# ./test_amber_cuda.x SPFP  (This will run the hybrid SPFP precision model)

. ../AmberTools/test/test_check.sh
check_environment `dirname \`pwd\`` serial

if [ $# -lt 1 ]; then
    echo "Using default PREC_MODEL = SPFP"
    PREC_MODEL="SPFP"
else
    PREC_MODEL=$1
fi

if [ ! -d $AMBERHOME/test/cuda ]; then
    echo "CUDA tests not found; You probably only have AmberTools"
    exit 0
fi

date_string=`date +%Y-%m-%d_%H-%M-%S`
logdir="${AMBERHOME}/logs/test_amber_cuda"
logprefix="${logdir}/${date_string}"
logfile="${logprefix}.log"
difffile="${logprefix}.diff"

if [ ! -z "$DO_PARALLEL" ]; then
   echo "Warning. DO_PARALLEL is set to $DO_PARALLEL. This environment variable is being unset."
   unset DO_PARALLEL
fi

mkdir -p ${logdir}
(make -k -f Makefile test.serial.cuda PREC_MODEL=$PREC_MODEL 2>&1) | tee ${logfile}

passed_count=`grep PASS ${logfile} | wc -l`
questionable_count=`grep "FAILURE:" ${logfile} | wc -l`
error_count=`grep "[Ee]rror" ${logfile} | wc -l`

if [ "${passed_count}" -eq 1 ] ; then
	echo "       1 file comparison passed" | tee -a ${logfile}
else
	echo "${passed_count} file comparisons passed" | tee -a ${logfile}
fi

if [ "${questionable_count}" -eq 1 ] ; then
	echo "       1 file comparison failed" | tee -a ${logfile}
else
	echo "${questionable_count} file comparisons failed" | tee -a ${logfile}
fi

if [ "${error_count}" -eq 1 ] ; then
	echo "       1 test experienced an error" | tee -a ${logfile}
else
	echo "${error_count} tests experienced errors" | tee -a ${logfile}
fi

echo "Test log file saved as ${logfile}"

if [ -f TEST_FAILURES.diff ]; then
   mv TEST_FAILURES.diff ${difffile}
   echo "Test diffs file saved as ${difffile}"
else
   echo "No test diffs to save!"
fi

