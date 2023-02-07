#!/bin/sh
# Daniel R. Roe
# 2010-04-02
# Certain tests require a stack that is at least 20480 kbytes.
# Check the current stack limit, then resize if necessary.
STACK=`ulimit -s`
if [ $STACK -lt 20480 ] ; then
  echo "This test requires a larger stack size; resizing." > /dev/stderr
  ulimit -s 20480
  STACK=`ulimit -s`
  if [ $STACK -lt 20480 ] ; then
    echo "Could not resize the stack." > /dev/stderr
    exit 1
  fi
fi
exit 0    
