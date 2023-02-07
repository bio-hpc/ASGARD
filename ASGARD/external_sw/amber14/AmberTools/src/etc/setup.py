
from distutils.core import setup
import os
import sys

# First the ParmedTools packages:
packages = ['cpinutils', 'mdoutanalyzer']

modules = []

# Scripts
scripts = ['cpinutil.py', 'charmmlipid2amber.py', 'mdout_analyzer.py',
           'softcore_setup.py']

if __name__ == '__main__':

    try:
        import argparse
    except ImportError:
        import shutil
        shutil.copyfile('argparse_amber.py', 'argparse.py')
        modules += ['argparse']

    try:
        from distutils.command.build_py import build_py_2to3 as build_py
        from distutils.command.build_scripts import build_scripts_2to3 as build_scripts
        PY3 = True
    except ImportError:
        from distutils.command.build_py import build_py
        from distutils.command.build_scripts import build_scripts
        PY3 = False

    setup(name='AmberTools',
          version='15.0',
          description='Various modules needed for AmberTools Python programs',
          author='Jason M. Swails, Ben Madej, and Thomas T. Joseph',
          author_email='jason.swails -at- gmail.com',
          url='http://ambermd.org',
          license='GPL v2 or later',
          packages=packages,
          py_modules=modules,
          cmdclass={'build_py': build_py, 'build_scripts': build_scripts},
          scripts=scripts)

    if os.path.exists('argparse.py'):
        os.unlink('argparse.py')
