setenv AMBERHOME "/mnt/home/users/ac_001_um/jorgedlpg/gromacs/amber14"
setenv PATH "${AMBERHOME}/bin:${PATH}"

# Add location of Amber Python modules to default Python search path
if( ! ($?PYTHONPATH) ) then
    setenv PYTHONPATH "${AMBERHOME}/lib/python2.7/site-packages"
else
    setenv PYTHONPATH "${AMBERHOME}/lib/python2.7/site-packages:${PYTHONPATH}"
endif
if( ! ($?LD_LIBRARY_PATH) ) then
   setenv LD_LIBRARY_PATH "${AMBERHOME}/lib"
else
   setenv LD_LIBRARY_PATH "${AMBERHOME}/lib:${LD_LIBRARY_PATH}"
endif
