include ../../config.h

test: Fortran Fortran2 C CPP Python clean

test.parm7: prepare.pmd
	$(BINDIR)/parmed.py -Oen -i prepare.pmd > /dev/null

Fortran: testf test.parm7
	@echo "Running testf"
	@./testf && /bin/rm -f testf

Fortran2: testf2
	@echo "Running testf2"
	@(./testf2 && /bin/rm -f testf2) || echo "Program error"

testf::
	@echo "Testing the Fortran API"
	$(FC) -I$(INCDIR) -L$(LIBDIR) test.F90 -o testf -lsander

testf2::
	@echo "Testing the Fortran file-less API"
	$(FC) -I$(INCDIR) -L$(LIBDIR) test2.F90 -o testf2 -lsander

C: testc test.parm7
	@echo "Running testc"
	@(./testc && /bin/rm -f testc) || echo "Program error"

testc::
	@echo "Testing the C API"
	$(CC) -I$(INCDIR) -L$(LIBDIR) test.c -o testc -lsander

CPP: testcpp test.parm7
	@echo "Running testcpp"
	@(./testcpp && /bin/rm -f testcpp) || echo "Program error"

testcpp::
	@echo "Testing the C++ API"
	$(CXX) -I$(INCDIR) -L$(LIBDIR) test.cpp -o testcpp -lsander

Python: test.parm7
	@echo "Testing the Python API"
	@echo "$(PYTHON) test.py"
	@$(PYTHON) test.py || echo "Program error"

clean:
	/bin/rm -f test.parm7 testf testf2 testc testcpp

.SUFFIXES:

.NOTPARALLEL: Fortran Fortran2 C CPP Python clean
