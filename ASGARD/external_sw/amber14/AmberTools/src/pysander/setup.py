
from distutils.core import setup, Extension
import os
from os.path import join
import shutil
import sys

amberhome = os.getenv('AMBERHOME')
assert amberhome is not None

if amberhome is None:
    raise RuntimeError('AMBERHOME is not set! Cannot compile pysander')

packages = ['sander', 'sanderles']

os.system('/bin/rm -fr sanderles')
shutil.copytree('sander', 'sanderles')
incdir = [join(amberhome, 'include'),
          join(amberhome, 'AmberTools', 'src', 'include')]
libdir = [join(amberhome, 'lib')]

try:
    pysander = Extension('sander.pysander',
                         sources=['sander/src/pysandermodule.c'],
                         include_dirs=incdir, library_dirs=libdir,
                         libraries=['sander'],
                         depends=['sander/src/pysandermoduletypes.c',
                                  join(incdir[1], 'CompatibilityMacros.h')],
    )
    pysanderles = Extension('sanderles.pysander',
                            sources=['sanderles/src/pysandermodule.c'],
                            include_dirs=incdir, library_dirs=libdir,
                            libraries=['sanderles'],
                            depends=['sander/src/pysandermoduletypes.c',
                                     join(incdir[1], 'CompatibilityMacros.h')],
                            define_macros=[('LES', None)])
    setup(name='sander',
          version="15.0",
          license='GPL v2 or later',
          author='Jason Swails',
          author_email='jason.swails -at- gmail.com',
          description='SANDER energy/force evaluation',
          ext_modules=[pysander, pysanderles],
          packages=packages,
    )
finally:
    os.system('/bin/rm -fr sanderles')
