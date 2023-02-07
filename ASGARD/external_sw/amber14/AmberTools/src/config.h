#  Amber configuration file, created with: ./configure gnu

###############################################################################

# (1)  Location of the installation

BASEDIR=/mnt/home/users/ac_001_um/jorgedlpg/gromacs/amber14
BINDIR=/mnt/home/users/ac_001_um/jorgedlpg/gromacs/amber14/bin
LIBDIR=/mnt/home/users/ac_001_um/jorgedlpg/gromacs/amber14/lib
INCDIR=/mnt/home/users/ac_001_um/jorgedlpg/gromacs/amber14/include
DATDIR=/mnt/home/users/ac_001_um/jorgedlpg/gromacs/amber14/dat
LOGDIR=/mnt/home/users/ac_001_um/jorgedlpg/gromacs/amber14/logs

###############################################################################


#  (2) If you want NAB to search additional libraries by default, add them
#      to the FLIBS variable here.  (External libraries can also be linked into
#      NAB programs simply by including them on the command line; libraries
#      included in FLIBS are always searched.)

FLIBS=  -lsff -lpbsa -lrism -lfftw3 -larpack -llapack -lblas  -lnetcdf  -lgfortran -w 
FLIBS_PTRAJ= -larpack -llapack -lblas   -lgfortran -w
FLIBSF= -larpack -llapack -lblas -lxblas-amb  
FLIBS_FFTW3= -lfftw3
###############################################################################

#  (3)  Modify any of the following if you need to change, e.g. to use gcc
#        rather than cc, etc.

SHELL=/bin/sh
INSTALLTYPE=serial
BUILDAMBER=

#  Set the C compiler, etc. 

#  The configure script should be fine, but if you need to hand-edit,
#  here is some info:

#   Example:  CC-->gcc; LEX-->flex; YACC-->yacc (built in byacc)
#     Note: If your lexer is "really" flex, you need to set
#     LEX=flex below.  For example, on some distributions,
#     /usr/bin/lex is really just a pointer to /usr/bin/flex,
#     so LEX=flex is necessary.  In general, gcc seems to need flex.

#   The compiler flags CFLAGS and CXXFLAGS should always be used.
#   By contrast, *OPTFLAGS and *NOOPTFLAGS will only be used with
#   certain files, and usually at compile-time but not link-time.
#   Where *OPTFLAGS and *NOOPTFLAGS are requested (in Makefiles,
#   makedepend and depend), they should come before CFLAGS or
#   CXXFLAGS; this allows the user to override *OPTFLAGS and
#   *NOOPTFLAGS using the BUILDFLAGS variable.

#   AMBERBUILDFLAGS provides a hook into all stages of the build process.
#   It can be used to build debug versions, invoke special features, etc.
#   Example:  make AMBERBUILDFLAGS='-O0 -g' sander
#
CC=gcc
CFLAGS=-fPIC -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DBINTRAJ -DHASGZ -DHASBZ2 -D__PLUMED_HAS_DLOPEN  $(CUSTOMBUILDFLAGS)  $(AMBERBUILDFLAGS)
CNOOPTFLAGS=
COPTFLAGS=-O3 -mtune=native 
AMBERCFLAGS= $(AMBERBUILDFLAGS)
WARNFLAGS=-Wall -Wno-unused-function

CXX=g++
CPLUSPLUS=g++
CXXFLAGS=-fPIC  $(CUSTOMBUILDFLAGS) $(AMBERBUILDFLAGS)
CXXNOOPTFLAGS=
CXXOPTFLAGS=-fPIC -O3 
AMBERCXXFLAGS= $(AMBERBUILDFLAGS)

NABFLAGS= $(AMBERBUILDFLAGS)
PBSAFLAG=-DFFTW $(AMBERBUILDFLAGS)

LDFLAGS= $(CUSTOMBUILDFLAGS) $(AMBERBUILDFLAGS)
AMBERLDFLAGS=$(AMBERBUILDFLAGS)

LEX=   flex
YACC=  $(BINDIR)/yacc
AR=    ar rv
M4=    m4
RANLIB=ranlib

#  Set the C-preprocessor.  Code for a small preprocessor is in
#    ucpp-1.3;  it gets installed as $(BINDIR)/ucpp;

CPP=ucpp -l

#  These variables control whether we will use compiled versions of BLAS
#  and LAPACK (which are generally slower), or whether those libraries are
#  already available (presumably in an optimized form).

LAPACK=install
BLAS=install
F2C=skip

#  These variables determine whether builtin versions of certain components
#  can be used, or whether we need to compile our own versions.

UCPP=install
C9XCOMPLEX=skip

#  For Windows/cygwin, set SFX to ".exe"; for Unix/Linux leave it empty:
#  Set OBJSFX to ".obj" instead of ".o" on Windows:

SFX=
OSFX=.o
MV=mv
RM=rm
CP=cp

#  Information about Fortran compilation:

FC=gfortran
FFLAGS= -fPIC $(LOCALFLAGS) $(CUSTOMBUILDFLAGS) -I$(INCDIR) $(NETCDFINC)  $(AMBERBUILDFLAGS)
FNOOPTFLAGS= -O0
FOPTFLAGS= -O3 -mtune=native
AMBERFFLAGS=$(AMBERBUILDFLAGS)
FREEFORMAT_FLAG= -ffree-form
LM=-lm
FPP=cpp -traditional -P
FPPFLAGS= -DBINTRAJ -DEMIL  $(CUSTOMBUILDFLAGS) $(AMBERBUILDFLAGS)
AMBERFPPFLAGS=$(AMBERBUILDFLAGS)
FCREAL8=-fdefault-real-8
NOFORTRANMAIN=-lgfortran -w
FWARNFLAGS=-Wall -Wno-unused-function

XHOME= /usr/X11R6
XLIBS= -L/usr/X11R6/lib64 -L/usr/X11R6/lib
MAKE_XLEAP=install_xleap

NETCDF=$(INCDIR)/netcdf.mod
NETCDFLIB=$(LIBDIR)/libnetcdf.a
NETCDFLIBF=$(LIBDIR)/libnetcdff.a $(LIBDIR)/libnetcdf.a
NETCDFINC=-I$(INCDIR)
FFTWLIB=-lfftw3

EMIL=EMIL
EMILLIB=$(LIBDIR)/libemil.a -lstdc++

ZLIB=-lz
BZLIB=-lbz2

HASFC=yes
MTKPP=install_mtkpp
XBLAS=$(LIBDIR)/libxblas-amb.a
FFTW3=$(LIBDIR)/libfftw3.a
MDGX=yes

COMPILER=gnu
MKL=
MKL_PROCESSOR=

#CUDA Specific build flags
NVCC=
PMEMD_CU_INCLUDES=
PMEMD_CU_LIBS=
PMEMD_CU_DEFINES=

#PMEMD Specific build flags
PMEMD_F90=gfortran   -DBINTRAJ -DEMIL -DPUBFFT
PMEMD_FOPTFLAGS=-O3 -mtune=native $(AMBERBUILDFLAGS)
PMEMD_CC=gcc
PMEMD_COPTFLAGS=-O3 -mtune=native -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DBINTRAJ  $(AMBERBUILDFLAGS)
PMEMD_FLIBSF=  $(LIBDIR)/libemil.a -lstdc++
PMEMD_LD= gfortran  $(AMBERBUILDFLAGS)
LDOUT= -o 

#for NAB:
MPI=

#1D-RISM
RISM=yes

#3D-RISM NAB
RISMSFF=-DRISMSFF
SFF_RISM_INTERFACE=../rism/amber_rism_interface.NAB.o
TESTRISMSFF=testrism

#3D-RISM SANDER
RISMSANDER=-DRISMSANDER
SANDER_RISM_INTERFACE=../rism/amber_rism_interface.SANDER.o
FLIBS_RISMSANDER=-lrism
TESTRISMSANDER=testrism

#for EMIL:
EMIL_MPIFLAGS=

#PUPIL
PUPILLIBS=-lrt -lm -lc -L${PUPIL_PATH}/lib -lPUPIL -lPUPILBlind

#Python interpreter we are using
PYTHON=/mnt/home/soft/python/programs/x86_64/bin/python2.7

#Type of Amber_Phenix build and python include file
AMBERPHENIX=no 
INCLUDE_PY=/export/home_users/home/soft/python/programs/x86_64/include/python2.7
PYTHON_VER=2.7
PHENIX_INC=
PHENIX_REPO=

PYSANDER=install
PYTRAJ=no_pytraj

#For LIO QM GPU Library
LIOLIBS=

# OS-specific rules for making shared objects
SHARED_SUFFIX=.so
MAKE_SHARED=-shared

# PLUMED related variables:
PLUMED_INCLUDE_FILE=
PLUMED_LOAD=Plumed.o -ldl -Wl,-export-dynamic
PLUMED_DEPENDENCIES=Plumed.o
