#!/bin/csh -f
#checks if the file qchem exists
#if it does, it is assumed that we have a valid Q-Chem installation

which qchem >& /dev/null
if( $status ) then
    # we didn't exit, so Q-Chem seems not installed...
    echo 'Q-Chem not installed - Skipping Test...'
    echo 'Check your Q-Chem installation and make sure'
    echo 'that the Q-Chem executable is called qchem'
    echo ''
    exit(1)
endif

exit(0)
