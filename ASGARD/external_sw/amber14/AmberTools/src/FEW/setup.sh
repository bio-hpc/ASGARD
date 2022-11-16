#!/bin/sh

if perl < /dev/null > /dev/null 2>&1  ; then

	if perl -MChemistry::Mol -le 'print $INC{"Chemistry/Mol.pm"}' < /dev/null > /dev/null 2>&1; then
		exit 0
	else
		cd $AMBERHOME/AmberTools/src/FEW/additional_libs
		mkdir PerlMol
		cd $AMBERHOME/AmberTools/src/FEW/additional_libs/PerlMol-0.3500
		perl Makefile.PL PREFIX=$AMBERHOME/AmberTools/src/FEW/additional_libs/PerlMol
		make
		make install
		cd $AMBERHOME/AmberTools/src/FEW
	fi
	
else
	echo No Perl available on system. Please install Perl.
	exit 1
fi
