# configured with [[ ./configure intel ]]

# Compilers
CXX = icpc
F90 = ifort
LD = ifort

# Flags
F90FLAGS = -warn all -O3 -xHost -ipo
CXXFLAGS = -Wall -O3 -xHost -ipo
LDFLAGS = -nofor_main -lstdc++ -lz -O3 -xHost -ipo

# Name of the program
PROGNAME = cphstats
PREFIX = /usr/local
