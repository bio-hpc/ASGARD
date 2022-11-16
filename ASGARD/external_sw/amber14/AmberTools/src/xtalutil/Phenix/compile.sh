#! /bin/bash

echo $LD_LIBRARY_PATH | tr ':' '\n' | grep amber


g++ -g -O0 tst_cdriver.cpp  -o tst_cdriver       \
		-I/home/pjanowsk/amberSD/include \
		-I/home/pawelrc/bin/phenix_svn/build/include \
        -I/home/pawelrc/bin/phenix_svn/source/cctbx_project \
        -I/usr/include/python2.7 \
        -I/home/pawelrc/bin/phenix_svn/source/amber/include \
        -L/home/pjanowsk/amberSD/lib \
        -lmdgx -lboost_python -lboost_system -lnetcdf -lfftw3 -lpython2.7 \
        ./phenix_amber_interface.so
