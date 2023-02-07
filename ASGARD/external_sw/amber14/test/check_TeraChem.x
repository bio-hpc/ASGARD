#!/bin/csh -f
#checks if the file terachem exists
#if it does, it is assumed that we have a valid TeraChem installation

which terachem >& /dev/null
if( $status ) then
    # we didn't exit, so TeraChem seems not installed...
    echo 'TeraChem not installed - Skipping Test...'
    echo 'Check your TeraChem installation and make sure'
    echo 'that the TeraChem executable is called terachem'
    echo ''
    exit(1)
endif

exit(0)

