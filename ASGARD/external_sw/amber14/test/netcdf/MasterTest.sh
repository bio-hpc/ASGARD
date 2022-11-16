# All test scripts in subdirectories should source this

# CleanFiles(): For every arg passed to the function, check for the file and rm it
CleanFiles() {
  while [[ ! -z $1 ]] ; do
    #for RMFILE in `find . -name "$1"` ; do
    if [[ -e $1 ]] ; then
      #echo "  Cleaning $1"
      rm $1
    fi
    #done
    shift
  done
  # If only cleaning requested no run needed, exit now
  if [[ $CLEAN -eq 1 ]] ; then
    exit 0
  fi
}

SetCpptraj() {
  CPPTRAJ="../../../bin/cpptraj"
  if [ ! -e "$CPPTRAJ" ] ; then
    CPPTRAJ=`which cpptraj`
    if [ -z "$CPPTRAJ" ] ; then
      echo "Warning: This check requires cpptraj to perform NetCDF to ASCII conversion."
      exit 0
    fi
  fi
}

CheckError() {
  if [[ $1 -ne 0 ]] ; then
    echo "${0}: Program error"
    exit 1
  fi
}

#==============================================================================
# If first argument is empty then NETCDF in config.h is empty, no BINTRAJ.
if [[ -z $1 ]] ; then
  echo "Warning: Amber compiled without NetCDF support. Skipping NetCDF test."
  exit 0
fi
# If the first argument is "clean" then no set-up is required. Script will
# exit when CleanFiles is called from sourcing script.
CLEAN=0
if [[ $1 = "clean" ]] ; then
  CLEAN=1
fi

# Get options, check environment
CPPTRAJ=""
if [[ $CLEAN -eq 0 ]] ; then
  # Check for TESTsander
  if [[ -z $TESTsander ]] ; then
    TESTsander="../../../bin/sander"
  fi
  # Check for DO_PARALLEL
  if [[ -z $DO_PARALLEL ]] ; then
    DO_PARALLEL=" "
  fi
fi

