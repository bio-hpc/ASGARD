#!/bin/sh

# This script is invoked via [py_setup.sh BINDIR PYTHON], so $1 is the bin
# directory and $2 is the python interpreter that we're using.

if [ $# -ne 2 ]; then
	# We must not have a PYTHON interpreter, so we can't do any python programs
	exit 0
fi

# If our current Python does not have argparse installed (it is only part of the
# Python stdlib as of Python 2.7 and later), then put it in AMBERHOME/bin
$2 -c "import argparse" > /dev/null 2>&1

if [ $? -ne 0 ]; then
   /bin/cp argparse_amber.py $1/argparse.py
fi

# Launch with the appropriate python

# cpinutil.py with modules
sed -e "s@PYTHONEXE@$2@g" < cpinutil.py > $1/cpinutil.py
/bin/chmod +x $1/cpinutil.py

$2 -c "from cpinutils import residues, utilities" > /dev/null 2>&1

if [ $? -ne 0 ]; then
   echo "Error importing cpinutil.py python modules! cpinutil will not work."
   exit 1
fi

/bin/cp -LR cpinutils $1/

# softcore.py
sed -e "s@PYTHONEXE@$2@g" < softcore_setup.py > $1/softcore_setup.py
/bin/chmod +x $1/softcore_setup.py

# mdout_analyzer.py
sed -e "s@PYTHONEXE@$2@g" < mdout_analyzer.py > $1/mdout_analyzer.py
/bin/chmod +x $1/mdout_analyzer.py

# Do not try importing mdoutanalyzer since it relies on Tkinter, matplotlib, and
# numpy; none of which are pre-requisites for Amber.
/bin/cp -LR mdoutanalyzer $1/

# charmmlipid2amber.py
sed -e "s@PYTHONEXE@$2@g" < charmmlipid2amber.py > $1/charmmlipid2amber.py
/bin/chmod +x $1/charmmlipid2amber.py
