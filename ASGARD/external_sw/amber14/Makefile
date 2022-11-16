#  Simple top-level Makefile to point users to those hidden below:

CONFIG_FILE=config.h
include $(CONFIG_FILE)

install:
	cd AmberTools/src && $(MAKE) install
	@(if [ -n "$(BUILDAMBER)" ] ; then \
	    cd src && $(MAKE) install; \
	  else \
		echo "==============================================================" ;\
		echo "$(AMBERHOME)/src/Makefile not found, or -noamber was set." ;\
		echo "This is expected if you do not have Amber14." ;\
		echo "==============================================================" ;\
	fi ;\
	)

clean:
	-(cd AmberTools/src && $(MAKE) clean)
	-(cd src && $(MAKE) clean)

uninstall:
	-(cd AmberTools/src && $(MAKE) uninstall)
	-(cd src && $(MAKE) uninstall)

clean.test:
	-(cd AmberTools/test && $(MAKE) clean)
	-(cd test && $(MAKE) clean)

test.clean:  clean.test

# This target cleans all generated files including all executables;
# the intention is that distclean returns the distribution to a
# pristine state as though the user just untarred a fresh copy of the code and
# applied their updates
distclean: clean clean.test uninstall
	/bin/rm -rf include/* lib/* ./AmberTools/src/cpptraj/src/readline/libreadline.a
	/bin/rm -f ./$(CONFIG_FILE) ./AmberTools/src/$(CONFIG_FILE)

test::  test.$(INSTALLTYPE)

test.serial:
	(cd AmberTools/test && $(MAKE) test)
	@(if [ -n "$(BUILDAMBER)" ] ; then \
		cd test && $(MAKE) test; \
		echo ""; \
		echo "Summary of AmberTools serial tests:"; \
		echo ""; cat ../logs/test_at_serial/at_summary; \
	else \
		echo "==============================================================" ;\
		echo "$(AMBERHOME)/src/Makefile not found." ;\
		echo "This is expected if you have not installed Amber14." ;\
		echo "==============================================================" ;\
	fi ;\
	)

test.parallel:
	(cd AmberTools/test && $(MAKE) test.parallel)
	@(if [ -n "$(BUILDAMBER)" ] ; then \
		cd test && $(MAKE) test.parallel; \
		echo ""; \
		echo "Summary of AmberTools parallel tests:"; \
		echo ""; cat ../logs/test_at_parallel/at_summary; \
	else \
		echo "==============================================================" ;\
		echo "$(AMBERHOME)/src/Makefile not found." ;\
		echo "This is expected if you have not installed Amber14." ;\
		echo "==============================================================" ;\
	fi ;\
	)

test.cuda:
	-(cd test && $(MAKE) test.cuda)

test.cuda_parallel:
	-(cd test && $(MAKE) test.cuda.parallel)

test.mic_offload:
	-(cd test && $(MAKE) test.mic_offload)

nabonly:
	-(cd AmberTools/src && $(MAKE) nabonly)

openmp:
	-(cd AmberTools/src && $(MAKE) openmp)

help:
	@echo "configure is the first step in the installation.  Execute this: "
	@echo "./configure --help"
