#!/bin/bash -f
#checks if the orca QM program is installed

if [ `which orca 2>/dev/null` ] ; then
	ORCA=`orca --help 2>/dev/null 2>/dev/null | grep -c Neese`
	if [ $ORCA -eq 1 ] ; then
		# orca exists and works
		exit 0
	fi
fi
# we didn't exit yet, so Orca seems not installed...
echo 'Orca not installed - Skipping Test...'
echo 'Check your Orca installation and make sure'
echo 'that the Orca executable is called orca'
echo ''
exit 1 

