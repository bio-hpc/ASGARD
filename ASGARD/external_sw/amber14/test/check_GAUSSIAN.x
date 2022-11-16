#!/bin/csh -f
#checks if the file g09 or g03 exists
#if it does, it is assumed that we have a valid GAUSSIAN installation

which g09 >& /dev/null
if( $status ) then
    which g03 >& /dev/null
    if( $status ) then
        # we didn't exit, so GAUSSIAN seems not installed...
        echo 'GAUSSIAN not installed - Skipping Test...'
        echo 'Check your GAUSSIAN installation and make sure'
        echo 'that the GAUSSIAN executable is called g09 or g03'
        echo ''
        exit(1)
    endif 
endif

exit(0)

