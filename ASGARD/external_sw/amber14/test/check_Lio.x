#!/bin/sh
#checks if the Lio QM library is installed

if [ -z $LIOHOME ]; then
   # We can't find lio libs.
   echo 'Lio is not installed - Skipping Test...'
   echo 'Check your Lio installation'
   echo 'and the LIOHOME environment setting'
   echo ''
   exit 1 
fi

# Otherwise LIO exists.
exit 0
