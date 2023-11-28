#!/usr/bin/env bash
#
#

#-nt          int    0       Total number of threads to start (0 is guess)
#-ntmpi       int    0       Number of thread-MPI threads to start (0 is guess)
#-ntomp       int    0       Number of OpenMP threads per MPI process/thread

if [ "${nodos}" != "1" ];then
	cores=`expr ${cores} \* ${nodos}` 
	thrads="-ntmpi $cores "
else
	thrads="-ntomp $cores "
fi

if [ "${GPU}" != "N/A" ];then
	gpu="-nb gpu"
	#gpu="-nb gpu_cpu"
else
	gpu=""
fi