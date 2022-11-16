#
# Basic include options.
# CC is the C compiler, normally invoked with options CFLAGS.
# LINKER and LDFLAGS function as CC and CFLAGS, but for linking stages.
#
CC = gcc
# For gfortran
#F2C_CONFIG=-DCONFIG_FC_UNDERSCORE
# For g77 or gfortran -ff2c
F2C_CONFIG=-DCONFIG_FC_DBL_UNDERSCORE -DCONFIG_FC_RETURNS_DBL_REAL
CFLAGS = -O3 -Dx86 -Wall $(F2C_CONFIG)
LINKER = $(CC)
LDFLAGS = 
EXTRA_LIBS = -lm

# M4 macro preprocessor
M4 	= m4
M4_OPTS = 
INDENT	= indent
INDENT_OPTS = -ce -i2 -nfc1 -br -brs -cs -npcs -nprs -npsl

#
#  The name of the libraries to be created/linked to
#
XBLASLIB = libxblas.a

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib
