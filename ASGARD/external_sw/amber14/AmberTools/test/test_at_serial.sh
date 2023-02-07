#!/bin/sh

. `dirname $0`/test_check.sh

date_string=`date +%Y-%m-%d_%H-%M-%S`
logdir="${AMBERHOME}/logs/test_at_serial"
logprefix="${logdir}/${date_string}"
logfile="${logprefix}.log"
difffile="${logprefix}.diff"

dir=`dirname \`pwd\``
check_environment `dirname $dir` serial

mkdir -p ${logdir}
(make -k -f Makefile test.serial 2>&1) | tee ${logfile}
(make -k -C ${AMBERHOME}/test test.serial 2>&1) | tee -a ${logfile}
(make -k -f Makefile finished 2>&1) | tee -a ${logfile}

passed_count=`grep PASS ${logfile} | wc -l`
questionable_count=`grep "FAILURE:" ${logfile} | wc -l`
error_count=`grep "Program error" ${logfile} | grep -v "echo" | wc -l`

echo "${passed_count} file comparisons passed" | tee -a ${logfile}
echo "${questionable_count} file comparisons failed" | tee -a ${logfile}
echo "${error_count} tests experienced errors" | tee -a ${logfile}

echo "Test log file saved as ${logfile}" | tee -a ${logfile}

if [ -f TEST_FAILURES.diff ]; then
   mv TEST_FAILURES.diff ${difffile}
   if [ -f ${AMBERHOME}/test/TEST_FAILURES.diff ]; then
      cat ${AMBERHOME}/test/TEST_FAILURES.diff >> ${difffile}
      /bin/rm ${AMBERHOME}/test/TEST_FAILURES.diff
   fi
   echo "Test diffs file saved as ${difffile}" | tee -a ${logfile}
else
   if [ -f ${AMBERHOME}/test/TEST_FAILURES.diff ]; then
      mv ${AMBERHOME}/test/TEST_FAILURES.diff ${difffile}
      echo "Test diffs file saved as ${difffile}" | tee -a ${logfile}
   else
      echo "No test diffs to save!" | tee -a ${logfile}
   fi
fi

# save summary for later reporting:
tail -5 ${logfile} > ${logdir}/at_summary
