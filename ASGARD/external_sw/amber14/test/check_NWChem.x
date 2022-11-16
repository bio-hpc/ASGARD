#!/bin/csh -f
#checks if the file nwchem exists and is executable
#if it does, it is assumed that we have a valid NWChem installation

which nwchem >& /dev/null
if( $status ) then
    # we didn't exit, so NWChem seems not installed...
    echo 'NWChem not installed - Skipping Test...'
    echo 'Check your NWChem installation and make sure'
    echo 'that the NWChem executable is called nwchem'
    echo ''
    exit(1)
endif

exit(0)

