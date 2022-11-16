#!/bin/sh

# This script is invoked via [setup.sh $BINDIR $PYTHON] which passes us
# where the programs should be installed and which python to use. If either
# is missing, we just stop

# This means we don't have a Python
if [ $# -lt 2 ]; then
   exit 0
fi

sed -e "s@PYTHONEXE@$2@g" < pymdpbsa > $1/pymdpbsa
sed -e "s@PYTHONEXE@$2@g" < pytleap > $1/pytleap
sed -e "s@PYTHONEXE@$2@g" < pdb4amber05/pdb4amber > $1/pdb4amber

/bin/chmod +x $1/pymdpbsa $1/pytleap $1/pdb4amber
