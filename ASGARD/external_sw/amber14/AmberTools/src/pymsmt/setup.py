
from distutils.core import setup
import os
import sys

# Packages in MSMT toolkit
packages = ['interface', 'mcpb', 'pymsmtlib', 'pymsmtmol']

# Modules
modules = ['pymsmtexp']

# Scripts
scripts = ['tools/MCPB.py', 'tools/OptC4.py', 'tools/PdbSearcher.py']

# See if our Python version will support OpenMM. Of the ParmEd-supported
# Pythons, only 2.4 and 2.5 do not work with OpenMM
major, minor = sys.version_info[:2]

if __name__ == '__main__':

    try:
        from distutils.command.build_py import build_py_2to3 as build_py
        from distutils.command.build_scripts import build_scripts_2to3 as build_scripts
        PY3 = True
    except ImportError:
        from distutils.command.build_py import build_py
        from distutils.command.build_scripts import build_scripts
        PY3 = False

    setup(name='pyMSMT',
          version='15.0', # For AmberTools 15
          description='Amber parameter file editor',
          author='Pengfei Li',
          author_email='ldsoar1990 -at- gmail.com',
          license='GPL v2 or later',
          packages=packages,
          py_modules=modules,
          cmdclass={'build_py':build_py, 'build_scripts':build_scripts},
          scripts=scripts)
