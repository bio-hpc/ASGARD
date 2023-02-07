#!/bin/csh -f
#checks if the file rungms exists
#if it does, it is assumed that we have a valid GAMESS installation

which rungms >& /dev/null
if( $status ) then
    # we didn't exit, so ADF seems not installed...
    echo 'GAMESS not installed - Skipping Test...'
    echo 'Check your GAMESS installation and make sure'
    echo 'that the GAMESS executable is called rungms'
    echo ''
    exit(1)
endif

exit(0)

