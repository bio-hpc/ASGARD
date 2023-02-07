#   top-level $(MAKE)file for the AmberTools suite of programs
AMBERTOOLS=AmberTools14
CONFIG_FILE=config.h
CONFIG_COMMAND=./configure
PARALLEL_TOKEN=DMPI
MIC_TOKEN=mmic
MIC_OFFLOAD_TOKEN=DMIC_offload

include $(CONFIG_FILE)

.NOTPARALLEL:configured_serial THIRDPARTY $(MTKPP)
install: $(INSTALLTYPE)
.PHONY: INSTALL_MTKPP

serial: configured_serial THIRDPARTY $(MTKPP)
	@echo "Starting installation of ${AMBERTOOLS} serial at `date`".
#	utility routines and libraries:
	(cd ucpp-1.3 && $(MAKE) $(UCPP) )
	(cd byacc && $(MAKE) install )
	(cd arpack && $(MAKE) install );
	(cd lapack && $(MAKE) $(LAPACK) )
	(cd blas && $(MAKE) $(BLAS) )
	(cd c9x-complex && $(MAKE) $(C9XCOMPLEX) )
	(cd pbsa && $(MAKE) install )
	(cd gbnsr6 && $(MAKE) install )
	(cd cifparse && $(MAKE) install )

#	antechamber:
	(cd antechamber && $(MAKE) install )
	(cd sqm && $(MAKE) install )

#   cpptraj:
	$(MAKE) cpptraj

#   miscellaneous:
	(cd reduce && $(MAKE) install )

#   leap
	(cd leap && $(MAKE) install )

#	old Amber programs
	(cd lib && $(MAKE) install)
	(cd addles && $(MAKE) install)
	(cd sander && $(MAKE) serial)
	(cd nmr_aux && $(MAKE) install)
	(cd nmode && $(MAKE) install)

#   nab:
	(cd sff && $(MAKE) install )
	(cd pbsa && $(MAKE) libinstall )
	(cd rism && $(MAKE) $(RISM) )
	(cd nab && $(MAKE) install )
	(cd nss && $(MAKE) install )
	(cd etc && $(MAKE) install )

#   rism
	(if [ $(RISM) = "yes" ]; then \
		cd rism && $(MAKE) rism3d.snglpnt$(SFX); \
	fi;\
	)

#   mdgx:
	(cd mdgx && $(MAKE) $(MDGX) )

#   xtalutil:
	(cd xtalutil/Clibs/ && $(MAKE) install )
	(cd xtalutil/CPrograms/ && $(MAKE) install )
	(cd xtalutil/Phenix/ && $(MAKE) $(AMBERPHENIX) )

#   saxs and saxs_md:
	(cd saxs && $(MAKE) install )

#	Parmed
	(cd parmed && $(PYTHON) setup.py install --prefix=$(AMBERHOME) --install-scripts=$(BINDIR))

#	mmpbsa_py
	(cd mmpbsa_py && $(MAKE) install)

#	mm_pbsa
	(cd mm_pbsa && $(MAKE) install)

#	FEW
	(cd FEW && $(MAKE) install)

#	amberlite
	(cd amberlite && $(MAKE) install)

#	paramfit
	(cd paramfit && $(MAKE) install)

#	cphstats
	(cd cphstats && $(MAKE) install)

#	pysander
	(cd pysander && $(MAKE) $(PYSANDER))

#	pymsmt
	(cd pymsmt && $(PYTHON) setup.py install --prefix=$(AMBERHOME) --install-scripts=$(BINDIR))

	@echo ""
	@echo "Installation of ${AMBERTOOLS} serial is complete at `date`."
	@echo ""

nabonly: configured_serial $(XBLAS)  
#	utility routines and libraries:
	(cd ucpp-1.3 && $(MAKE) $(UCPP) )
	(cd byacc && $(MAKE) install )
	(cd arpack && $(MAKE) install );
	(cd lapack && $(MAKE) $(LAPACK) )
	(cd blas && $(MAKE) $(BLAS) )
	(cd c9x-complex && $(MAKE) $(C9XCOMPLEX) )
	(cd pbsa && $(MAKE) libinstall )
	(cd rism && $(MAKE) $(RISM) )
	(cd cifparse && $(MAKE) install )
	(cd sff && $(MAKE) install )
	(cd nab && $(MAKE) install )
	(cd nss && $(MAKE) install )

	@echo
	@echo "Completed installation of NAB"
	@echo

openmp: THIRDPARTY nabonly 
	$(MAKE) cpptraj
	(cd saxs && $(MAKE) install )
	(cd paramfit && $(MAKE) install )

	@echo ""
	@echo "Installation of ${AMBERTOOLS} openmp is complete at `date`."
	@echo ""

mtkpponly: configured_serial install_mtkpp

install_mtkpp:
	-(cd mtkpp && $(MAKE) && $(MAKE) -j 1 install )
	-(cd mtkpp/tests && $(MAKE) create )

.NOTPARALLEL:configured_parallel THIRDPARTY

parallel: configured_parallel THIRDPARTY
	@echo "Starting installation of ${AMBERTOOLS} MPI at `date`".
	(cd etc && $(MAKE) parallel)
	(cd sff && $(MAKE) libsff_mpi )
#	(cd pbsa && $(MAKE) clean && $(MAKE) install.parallel )
	(cd pbsa && $(MAKE) libinstall )
	$(MAKE) cpptraj
	(cd byacc && $(MAKE) install )
	(cd ucpp-1.3 && $(MAKE) $(UCPP) )
	(cd cifparse && $(MAKE) install )
	(cd arpack && $(MAKE) install );
	(cd lapack && $(MAKE) $(LAPACK) )
	(cd blas && $(MAKE) $(BLAS) )
	(cd mdgx && $(MAKE) $(MDGX) )
	(if [ "$(RISM)" = "yes" ]; then \
		cd rism && $(MAKE) install_mpi ; \
	fi; \
	)
	(cd nab && $(MAKE) install )
	(cd parmed && $(PYTHON) setup.py install --prefix=$(AMBERHOME))
	(cd mmpbsa_py && $(MAKE) parallel )
	(if [ "$(RISM)" = "yes" ]; then \
		cd rism && $(MAKE) rism3d.snglpnt.MPI$(SFX) ; \
	fi; \
	)
# Old Amber programs
	(cd addles && $(MAKE) install)
	(cd sander && $(MAKE) parallel)

	@echo ""
	@echo "Installation of ${AMBERTOOLS} MPI is complete at `date`."
	@echo ""

configured:
	@(if [ ! -f $(CONFIG_FILE) ] ; then \
		echo "Error: $(CONFIG_COMMAND) must be executed before $(MAKE) !" ;\
		exit 2 ;\ # $(CONFIG_COMMAND) ;\
	fi ;\
	)

configured_serial: configured
	@(if grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
		echo "Error: $(CONFIG_FILE) is of type parallel, not serial !" ;\
		echo "  Rerun $(CONFIG_COMMAND) without the -mpi flag." ;\
		exit 2 ;\
	fi ;\
	)

configured_parallel: configured
	@(grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null || \
		{ echo "Error: $(CONFIG_FILE) is of type serial, not parallel !" ;\
			echo "  Rerun $(CONFIG_COMMAND) with the -mpi flag." ;\
			exit 2 ;\
		} ;\
	)

configured_mic: configured
	@(grep $(MIC_TOKEN) $(CONFIG_FILE) > /dev/null || \
		{ echo "Error: $(CONFIG_FILE) is not of type mic_native !" ;\
			echo "Rerun $(CONFIG_FILE) WITH THE -mic_native flag." ;\
			exit 2 ;\
		} ;\
	)

configured_mic_offload: configured
	@(grep $(MIC_OFFLOAD_TOKEN) $(CONFIG_FILE) > /dev/null || \
                { echo "Error: $(CONFIG_FILE) is not of type mic_offload !" ;\
                        echo "Rerun $(CONFIG_FILE) WITH THE -mic_offload flag." ;\
                        exit 2 ;\
                } ;\
        )

netcdf_uninstall:
	-(cd netcdf-4.3.0 && $(MAKE) uninstall)
	-(cd netcdf-fortran-4.2 && $(MAKE) uninstall)
	-/bin/rm -f netcdf_config.log ../lib/libnetcdf* ../include/netcdf.*

clean::
	-(cd ucpp-1.3 && $(MAKE) clean )
	-(cd byacc && $(MAKE) clean )
	-(cd cifparse && $(MAKE) clean )
	-(cd nab && $(MAKE) clean )
	-(cd sff && $(MAKE) clean )
	-(cd nss && $(MAKE) clean )
	-(cd leap && $(MAKE) clean )
	-(cd arpack && $(MAKE) clean )
	-(cd blas && $(MAKE) clean )
	-(cd xblas && $(MAKE) clean )
	-(cd lapack && $(MAKE) clean )
	-(cd c9x-complex && $(MAKE) clean )
	-(cd etc && $(MAKE) clean )
	-(cd chamber && $(MAKE) clean )
	-(cd pbsa && $(MAKE) clean )
	-(cd gbnsr6 && $(MAKE) clean )
	-(cd antechamber && $(MAKE) clean )
	-(cd sqm && $(MAKE) clean )
	-(cd lib && $(MAKE) clean )
	-(cd ptraj && $(MAKE) clean )
	-(cd cpptraj && $(MAKE) clean)
	-(cd reduce && $(MAKE) clean )
	-(cd mtkpp && $(MAKE) clean )
	-(cd mdgx && $(MAKE) clean )
	-(cd xtalutil/CPrograms && $(MAKE) clean )
	-(cd saxs && $(MAKE) clean )
	-(cd paramfit && $(MAKE) clean )
	-(cd rism && $(MAKE) clean )
	-(cd fftw-3.3 && $(MAKE) clean )
	-($(RM) -f fftw-3.3/mpi/fftw3-mpi.f03 )
	-(cd mmpbsa_py && /bin/rm -fr MMPBSA.py.MPI build/ mpi4py-1.2.2/)
	-(cd mm_pbsa && $(MAKE) clean )
	-(cd FEW && $(MAKE) clean )
	-(cd parmed && /bin/rm -fr build/)
	-(cd cphstats && $(MAKE) clean )
	-(cd lib && $(MAKE) clean)
	-(cd addles && $(MAKE) clean)
	-(cd sander && $(MAKE) clean)
	-(cd nmr_aux && $(MAKE) clean)
	-(cd nmode && $(MAKE) clean)
	-(cd emil && $(MAKE) clean)
	-(cd sebomd && $(MAKE) clean)
	-(cd pymsmt && $(MAKE) clean)

uninstall: netcdf_uninstall
	-(cd ucpp-1.3 && $(MAKE) uninstall )
	-(cd byacc && $(MAKE) uninstall )
	-(cd cifparse && $(MAKE) uninstall )
	-(cd nab && $(MAKE) uninstall )
	-(cd sff && $(MAKE) uninstall )
	-(cd nss && $(MAKE) uninstall )
	-(cd leap && $(MAKE) uninstall )
	-(cd arpack && $(MAKE) uninstall )
	-(cd blas && $(MAKE) uninstall )
	-($(RM) -f $(LIBDIR)/libxblas.a )
	-(cd lapack && $(MAKE) uninstall )
	-(cd c9x-complex && $(MAKE) uninstall )
	-(cd etc && $(MAKE) uninstall )
	-(cd chamber && $(MAKE) uninstall )
	-(cd pbsa && $(MAKE) uninstall )
	-(cd gbnsr6 && $(MAKE) uninstall )
	-(cd antechamber && $(MAKE) uninstall )
	-(cd sqm && $(MAKE) uninstall )
	-(cd lib && $(MAKE) uninstall )
	-(cd ptraj && $(MAKE) uninstall )
	-(cd reduce && $(MAKE) uninstall )
	-(cd mdgx && $(MAKE) uninstall )
	-(cd xtalutil/CPrograms && $(MAKE) uninstall )
	-(cd saxs && $(MAKE) uninstall )
	-(cd mtkpp && $(MAKE) uninstall )
	-(cd paramfit && $(MAKE) uninstall )
	-(cd rism && $(MAKE) uninstall )
	-(cd fftw-3.3 && $(MAKE) uninstall)
	-(cd cpptraj && $(MAKE) uninstall)
	-(cd mmpbsa_py && $(MAKE) uninstall)
	-(cd mm_pbsa && $(MAKE) uninstall)
	-(cd FEW && $(MAKE) uninstall)
	-(cd parmed && /bin/rm -fr $(BINDIR)/*parmed.py)
	-(cd amberlite && $(MAKE) uninstall )
	-/bin/rm -f ../test/numprocs
	-/bin/rm -fr $(BINDIR)/ncmpidump $(BINDIR)/ncmpigen $(BINDIR)/ncvalid $(LIBDIR)/libnetcdf* $(INCDIR)/mtkpp
	-(cd $(LIBDIR) && /bin/rm -f libxblas-amb.a libfftw3_mpi.a libfftw3_mpi.la)
	-(cd $(LIBDIR) && /bin/rm -f libmdgx.so)
	-(cd $(LIBDIR) && /bin/rm -f libemil.a)
	-(cd $(INCDIR) && /bin/rm -f fftw3-mpi.f03 fftw3-mpi.h)
	-(cd lib && $(MAKE) uninstall)
	-(cd addles && $(MAKE) uninstall)
	-(cd sander && $(MAKE) uninstall.serial uninstall.parallel)
	-(cd nmr_aux && $(MAKE) uninstall)
	-(cd nmode && $(MAKE) uninstall)
	-(cd cphstats && $(MAKE) uninstall)
	-(cd pymsmt && $(MAKE) uninstall)
	-(cd $(LIBDIR) && /bin/rm -f python*/)

$(LIBDIR)/libxblas-amb.a:
	(cd xblas && $(MAKE) -j 1 lib-amb && mv libxblas-amb.a $(LIBDIR) )

$(LIBDIR)/libfftw3.a $(LIBDIR)/libfftw3_mpi.a: 
	cd fftw-3.3 && $(MAKE) && $(MAKE) -j 1 install;

THIRDPARTY: $(XBLAS) $(FFTW3)

cuda:
	@echo "$(AMBERTOOLS) has no CUDA-enabled components"

	-(cd ../../src && $(MAKE) cuda )

cuda_parallel:
	@echo "$(AMBERTOOLS) has no CUDA-enabled components"

	-(cd ../../src && $(MAKE) cuda_parallel )

mic: configured_mic 
	@echo "$(AMBERTOOLS) has no mic-enabled components"

	-(cd ../../src && $(MAKE) mic )

mic_parallel: configured_mic configured_parallel 
	@echo "$(AMBERTOOLS) has no mic-enabled components"

	-(cd ../../src && $(MAKE) mic_parallel )

mic_offload: configured_mic_offload 
	@echo "$(AMBERTOOLS) has no mic-enabled components"

	-(cd ../../src && $(MAKE) mic_offload )

# Rule for building cpptraj
cpptraj::
	(cd cpptraj && $(MAKE) $(INSTALLTYPE))

# Add a rule to build ptraj for anyone that wants to compile it easily.
ptraj::
	(cd ptraj && $(MAKE) install )

# Add a rule to build chamber for anyone that wants to compile it easily
chamber::
	(cd chamber && $(MAKE) install )
