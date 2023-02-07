
from distutils.core import setup
import os
import sys

# First the ParmedTools packages:
packages = ['MMPBSA_mods']

# Scripts
scripts = ['MMPBSA.py', 'ante-MMPBSA.py']

if __name__ == '__main__':

    parallel = False
    for i, arg in enumerate(sys.argv):
        # Parallel builds only build chemistry and compat24
        if arg == '--par':
            scripts = ['MMPBSA.py.MPI']
            import shutil
            shutil.copyfile('MMPBSA.py', 'MMPBSA.py.MPI')
            sys.argv.pop(i)
            parallel = True
            break

    # See if we need to install mpi4py or not
    try:
        from mpi4py import MPI
        if MPI.COMM_WORLD.Get_rank() != 0 or MPI.COMM_WORLD.Get_size() != 1:
            raise RuntimeError('mpi4py is not working')
        mpi4py_works = True
    except (ImportError, RuntimeError):
        mpi4py_works = False

    try:
        from distutils.command.build_py import build_py_2to3 as build_py
        from distutils.command.build_scripts import build_scripts_2to3 as build_scripts
        PY3 = True
    except ImportError:
        from distutils.command.build_py import build_py
        from distutils.command.build_scripts import build_scripts
        PY3 = False

    if parallel and not mpi4py_works:
        os.system('tar zxf mpi4py-1.2.2.tar.gz')
        os.chdir('mpi4py-1.2.2')
        try:
            print('Building mpi4py... This may take a few minutes.')
            ret = os.system('%s setup.py build > ../../mpi4py_install.log 2>&1'
                            % sys.executable)
            if ret != 0:
                sys.exit(' Error in mpi4py install. Check mpi4py_install.log')
            ret = os.system('%s setup.py install --prefix=%s' % (sys.executable,
                            os.getenv('AMBERHOME')))
            if ret != 0:
                sys.exit(' Error in mpi4py install. Check mpi4py_install.log')
        finally:
            os.chdir('..')

    setup(name='MMPBSA.py',
          version='15.0',
          description='Program for carrying out MM/PBSA-like end-state '
                      'free energy calculations',
          author='Jason Swails, T. Dwight McGee Jr., and Bill R. Miller III',
          author_email='jason.swails -at- gmail.com',
          url='http://ambermd.org',
          license='GPL v2 or later',
          packages=packages,
          cmdclass={'build_py': build_py, 'build_scripts': build_scripts},
          scripts=scripts)

    # Delete our temporary file
    if parallel:
        os.unlink('MMPBSA.py.MPI')
