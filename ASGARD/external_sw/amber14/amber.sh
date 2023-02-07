export AMBERHOME="/mnt/home/users/ac_001_um/jorgedlpg/gromacs/amber14"
export PATH="${AMBERHOME}/bin:${PATH}"

# Add location of Amber Python modules to default Python search path
if [ -z "$PYTHONPATH" ]; then
    export PYTHONPATH="${AMBERHOME}/lib/python2.7/site-packages"
else
    export PYTHONPATH="${AMBERHOME}/lib/python2.7/site-packages:${PYTHONPATH}"
fi
if [ -z "${LD_LIBRARY_PATH}" ]; then
   export LD_LIBRARY_PATH="${AMBERHOME}/lib"
else
   export LD_LIBRARY_PATH="${AMBERHOME}/lib:${LD_LIBRARY_PATH}"
fi
