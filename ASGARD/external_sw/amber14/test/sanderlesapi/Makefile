include ../../config.h

test: Fortran C CPP Python

Fortran: testf
	./testf && /bin/rm -f testf

testf::
	@echo "Testing the Fortran LES API"
	$(FC) -I$(INCDIR) -L$(LIBDIR) test.F90 -o testf -lsanderles

C: testc
	(./testc && /bin/rm -f testc) || echo "Program error"

testc::
	@echo "Testing the C LES API"
	$(CC) -I$(INCDIR) -L$(LIBDIR) test.c -o testc -lsanderles

CPP: testcpp
	(./testcpp && /bin/rm -f testcpp) || echo "Program error"

testcpp::
	@echo "Testing the C++ LES API"
	$(CXX) -I$(INCDIR) -L$(LIBDIR) test.cpp -o testcpp -lsanderles

Python:
	@echo "Testing the Python LES API"
	($(PYTHON) test.py) || echo "Program error"

.SUFFIXES:
