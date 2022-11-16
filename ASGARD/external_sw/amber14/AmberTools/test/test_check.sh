# Provides a function to check that the Amber environment (specifically
# pertaining to LD_LIBRARY_PATH) is sufficient before running tests that will
# ultimately fail. This prevents spewing garbage and replaces it with a succinct
# error message

check_environment() {
   # check_environment should be called with 2 arguments. The first should be
   # the value of AMBERHOME. The second should be either "serial" or "parallel"

   # First check that AMBERHOME is set
   if [ -z "$AMBERHOME" ]; then
      echo ""
      echo "Error: AMBERHOME is not set!"
      echo ""
      echo "Set AMBERHOME to $1 and re-run the tests;"
      echo "alternatively, this and other commands for updating your environment"
      echo "are in \$AMBERHOME/amber.sh and \$AMBERHOME/amber.csh which can be sourced"
      echo "now or from your shell resource file (e.g., ~/.bashrc or ~/.cshrc)."
      echo ""
      exit 1
   fi

   # For parallel simulations, check that DO_PARALLEL is not empty. For serial
   # simulations, empty DO_PARALLEL and warn if it was set
   if [ "$2" = "serial" ]; then
      if [ ! -z "$DO_PARALLEL" ]; then
         echo ""
         echo "Warning: DO_PARALLEL is set to \"$DO_PARALLEL\"."
         echo "This environment variable is being unset for serial testing."
         echo ""
         unset DO_PARALLEL
      fi
   fi

   if [ "$2" = "parallel" ]; then
      if [ -z "$DO_PARALLEL" ]; then
         echo ""
         echo "Error: DO_PARALLEL must be set for parallel tests!"
         echo ""
         exit 1
      fi
      echo ""
      echo "Tests being run with DO_PARALLEL=\"$DO_PARALLEL\"."
      echo ""
   fi

   # Unset TESTsander, since that's used by most tests
   if [ ! -z "$TESTsander" ]; then
      echo ""
      echo "Warning: TESTsander set to \"$TESTsander\"."
      echo "Unsetting TESTsander for tests."
      echo ""
      unset TESTsander
   fi

   # Now make sure that AMBERHOME/lib is inside LD_LIBRARY_PATH if either the
   # NetCDF or FFTW3 shared libraries were built by Amber. This is not necessary
   # on macs.
   is_mac='no'
   uname -a | grep "Darwin" 2>&1 > /dev/null
   test $? -eq 0 && is_mac='yes'

   # We are done here for Macs
   test $is_mac = "yes" && return

   python << EOF
import os
import sys
ambhome = os.getenv('AMBERHOME')
real_ambhome = os.path.realpath(os.path.abspath(ambhome))
if ambhome is None:
    sys.stderr.write('AMBERHOME is not set. This must be set to run the tests!\\n')
    sys.exit(1)

def error():
    sys.stderr.write('We recommend adding the line:\\n\\n')
    sys.stderr.write('   test -f %s/amber.sh && source %s/amber.sh (sh/bash/zsh)\\n' % (ambhome, ambhome))
    sys.stderr.write('or\\n')
    sys.stderr.write('   test -f %s/amber.csh && source %s/amber.csh (csh/tcsh)\\n' % (ambhome, ambhome))
    sys.stderr.write('\\nto your login shell resource file (e.g., ~/.bashrc or ~/.cshrc)\\n')
    sys.exit(1)

try:
    import chemistry
except ImportError:
    sys.stderr.write('Could not import Amber Python modules. This likely means\\n'
        'that your Amber Python environment was not set up correctly\\n\\n')
    error()

if 'darwin' in sys.platform:
    sys.exit(0) # Nothing to check here

ldlib = os.getenv('LD_LIBRARY_PATH')
if ldlib is None:
    sys.stderr.write('LD_LIBRARY_PATH is not set. Make sure you\\n'
            'source %s/amber.sh for sh/bash/zsh\\n'
            '(or %s/amber.csh for csh/tcsh)\\n'
            'in order to properly set your environment up for using Amber\\n' %
            (ambhome, ambhome)
    )
for folder in ldlib.split(os.pathsep):
    folder = folder.strip()
    if not folder:
        continue
    realpath = os.path.realpath(os.path.abspath(folder))
    if realpath.startswith(os.path.join(real_ambhome, 'lib')):
        break
else:
    sys.stderr.write('Error: LD_LIBRARY_PATH does not include \$AMBERHOME/lib!\\n')
    sys.stderr.write('\\n')
    sys.stderr.write('Amber now requires \$AMBERHOME/lib to be added to your LD_LIBRARY_PATH\\n')
    sys.stderr.write('environment variable in order for all components to work.\\n')
    error()
sys.exit(0)
EOF
    test $? -eq 0 || exit 1
}
