#!/bin/csh -f
#checks if the file adf.exe exists
#if it does, it is assumed that we have a valid ADF installation

which adf.exe >& /dev/null
if( $status ) then
    # we didn't exit, so ADF seems not installed...
    echo 'ADF not installed - Skipping Test...'
    echo 'Check your ADF installation and make sure'
    echo 'that the ADF executable is called adf.exe'
    echo ''
    exit(1)
endif

exit(0)
