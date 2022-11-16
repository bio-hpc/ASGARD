include ../src/config.h

######    Top-level targets, for general use:    ############################

test::
	./test_at_serial.sh

test.parallel::
	./test_at_parallel.sh

######    Intermediate-level targets:    ############################

test.serial: clean is_amberhome_defined \
	test.nab test.cpptraj test.antechamber test.mdgx \
	test.leap test.resp test.pbsa test.gbnsr6 test.mmpbsa test.parmed \
	test.elsize test.sqm test.rism1d test.amberlite \
	test.paramfit test.mm_pbsa test.FEW test.cphstats \
	test.charmmlipid2amber test.cpinutil

test.parallel2: clean is_amberhome_defined \
        test.nab test.mmpbsa

test.nab:
	cd nab && $(MAKE) -k test $(TESTRISMSFF)

test.ptraj:
	-cd ptraj_rmsa && ./Run.rms
	-cd ptraj_rms && ./Run.rms
#	-cd ptraj_les && ./Run.les
	-cd ptraj_matrix && ./Run.matrix
	-cd ptraj_order && ./Run.order
	-cd ptraj_trajectory && ./Run.ptraj
	-cd ptraj_remd_netcdf && ./Run.ptraj $(NETCDF)
	-cd ptraj_comprehensive && ./Run.comprehensive $(NETCDF)

test.antechamber:
	-cd antechamber/bondtype && ./Run.bondtype
	-cd antechamber/c60 && ./Run.c60
	-cd antechamber/charmm && ./Run.charmm
	-cd antechamber/chemokine && ./Run.1b2t
	-cd antechamber/tp && ./Run.tp
	-cd antechamber/ash && ./Run.ash
	-cd antechamber/sustiva && ./Run.sustiva
	-cd antechamber/fluorescein && ./Run.fluorescein
	-cd antechamber/guanine_amber && ./Run.guanine_amber
	-cd antechamber/fpph && ./Run.fpph_resp
	-cd antechamber/residuegen && ./Run.ala
	-cd antechamber/parmchk2 && ./Run.parmchk2

test.leap:
	-cd leap && ./Run.tleap
	-cd leap && ./Run.parm7
	-cd leap && ./Run.cmds
	-cd leap && ./Run.mbondi3
	-cd leap && ./Run.addions
	-cd leap/glycam/06j && ./Run.glycam
	-cd leap/glycam/06EPb && ./Run.glycam
	-cd leap/glycam/06j_10 && ./Run.glycam
	-cd leap/glycam/06j_12SB && ./Run.glycam

test.cpptraj:
	cd cpptraj && $(MAKE) -k test

test.paramfit:
	cd paramfit && $(MAKE) -k test

test.resp:
	cd resp && ./Run.ade

test.reduce:
	cd reduce && ./Run.reduce

test.mdgx:
	-cd mdgx/1p7e && ./Test.1p7e
	-cd mdgx/NVTwat && ./Test.nvtwat
	-cd mdgx/MLE && ./Test.mle
	-cd mdgx/Settle && ./Test.settle
	-cd mdgx/MCBarostat && ./Test.mcb
	-cd mdgx/TIP4P && ./Test.tip4p
	-cd mdgx/ForceCalc && ./Test.frcc
	-cd mdgx/AddEP && ./Test.addep
	-cd mdgx/QFit && ./Test.qfit
	-cd mdgx/Triclinic && ./Test.triclinic
#	-cd mdgx/Belly && ./Test.belly


test.mtkpp:
	cd mtkpp && ./Run

test.pbsa:
	-cd pbsa_dmp && ./Run.dmp.min
	-cd pbsa_p3m && ./Run.dmp.min
	-cd pbsa_bcopt && ./Run.dmp.min
	-cd pbsa_saopt && ./Run.dmp.min
	-cd pbsa_tri && ./Run.dmp.min
	-cd pbsa_rot && ./Run.charge.min
	-cd pbsa_lpb && ./Run.lsolver.min
	-cd pbsa_tsr && ./Run.tsrb.min
	-cd pbsa_radi && ./Run.ion.min
	-cd pbsa_npb && ./test
	-cd pbsa_frc && ./test
	-cd pbsa_ipb2 && ./Run.110D.min
	-cd pbsa_delphi && ./test
	-cd pbsa_IIM && ./Run.IIM.min
	-cd pbsa_ligand && ./Run.t4bnz.min
	-cd pbsa_mbfocus && ./test
	-cd pbsa_aug_c2 && ./Run.ion.min
	-cd pbsa_memb_test && ./Run.memb.test
	-cd pbsa_memb_test && ./Run.solv.test
	-cd pbsa_phiout && ./Run.dmp.min

test.pbsa.MPI:
	-cd pbsa_mbfocus && \
	   env TESTpbsa='../../bin/pbsa.MPI' ./Run.mb

test.gbnsr6:
	-cd gbnsr6 && ./Run.gbnsr6.min

test.ambpdb:
	-cd ambpdb && ./Run.ambpdb
	-cd ambpdb_first_protein && ./Run.ambpdb
	-cd ambpdb_first_rna && ./Run.ambpdb

test.elsize:
	cd elsize && ./Run.elsize

test.chamber:
	-cd ./chamber/trans-but-2-ene && ./Run
	-cd ./chamber/ala_ala_ala && ./Run
	-cd ./chamber/ala_ala_ala && ./Runx
	-cd ./chamber/poly_pro_gas && ./Run.poly_pro_gas_chamber
	-cd ./chamber/glucose && ./Run
	-cd ./chamber/dhfr_gas && ./Run.dhfr_gas_chamber
	-cd ./chamber/dhfr_cmap_pbc && ./Run.dhfr_cmap_pbc_chamber
	-cd ./chamber/2b1q-trlose-wat-c36 && ./Run
	-cd ./chamber/test-c27-vmd-psf && ./Run

test.cphstats:
	-cd cphstats && ./Run.cphstats

test.sqm:
	-cd ./sqm/c2h2 && ./Run
	-cd ./sqm/ch3f && ./Run
	-cd ./sqm/nma && ./Run
	-cd ./sqm/znh+ && ./Run
	-cd ./sqm/h2s && ./Run
	-cd ./sqm/dhplus && ./Run
	-cd ./sqm/fecp2 && ./Run
	-cd ./sqm/external_charges && ./Run

test.rism1d:
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/spc-kh && ./Run.spc; fi
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/spc-psen && ./Run.spc-nacl-3; fi
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/spc-hnc && ./Run.spc-lj; fi
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/spc-mv0 && ./Run.spc; fi
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/lj-py && ./Run.lj; fi
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/tip3p-kh && ./Run.tip3p; fi

test.rism1d.dev: test.rism1d
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/dev/spc-kh && ./Run.spc-nacl; fi
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/dev/spc-kh && ./Run.spc-lj; fi
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/dev/spc-psen && ./Run.spc-1; fi
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/dev/spc-psen && ./Run.spc-3; fi
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/dev/spc-psen && ./Run.spc-nacl-1; fi
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/dev/spc-psen && ./Run.spc-nacl-2; fi
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/dev/spc-hnc && ./Run.spc; fi
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/dev/tip3p-kh && ./Run.tip3p-nacl; fi
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/dev/tip3p-kk-2008 && ./Run.tip3p-nacl-1; fi
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/dev/tip3p-kk-2008 && ./Run.tip3p-nacl-2; fi
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/dev/tip3p-kk-2008 && ./Run.tip3p-nacl-3; fi
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/dev/tip3p-kk-2008 && ./Run.tip3p-nacl-4; fi
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/dev/dt && ./Run.all; fi
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism1d/dev/ && $(MAKE); fi

test.rism3d.thermo:
	-if [ -n "$(TESTRISMSFF)" ]; then cd ./rism3d.thermo/ala && ./Run.thermoala; fi

test.mmpbsa:
	cd mmpbsa_py && $(MAKE) test

test.mm_pbsa:
	cd ../src/mm_pbsa/Examples && ./Run.mmpbsa.test

test.FEW:
	cd ../src/FEW/examples/test && ./Run.few.test

test.parmed:
	cd parmed/normal_prmtop && ./Run.parmed
	cd parmed/chamber_prmtop && ./Run.parmed
	cd parmed/solvated_prmtop && ./Run.parmed
	cd parmed/change_parms && ./Run.parmed
	cd parmed/multifiles && ./Run.parmed
	cd parmed/add12_6_4 && ./Run.parmed
	cd parmed/ti_merge && ./Run.parmed
	cd parmed/chamber && ./Run.parmed

test.cpinutil:
	cd cpinutil && ./Run.cpin

test.charmmlipid2amber:
	cd charmmlipid2amber && ./Run.charmmlipid2amber

test.amberlite:
	cd amberlite && ./Run.pdb4amber

finished:
	@echo ""
	@echo "Finished test suite for AmberTools at `date`."
	@echo ""

is_amberhome_defined:
	@(if [ -z "$(AMBERHOME)" ] ; then \
	    echo "Error: AMBERHOME should be defined or else some tests will fail !" ;\
	    exit 2 ;\
	fi ;\
	)

clean:
	-./test_at_clean.sh
	-cd mmpbsa_py && $(MAKE) -k clean
	-cd ../src/mm_pbsa/Examples && ./Run.mmpbsa.clean
	-cd ../src/FEW/examples/test && ./Run.few.clean
	-cd cpptraj && $(MAKE) -k clean
	-cd leap/glycam/06j && ./Run.glycam clean
	-cd leap/glycam/06EPb && ./Run.glycam clean
	-cd leap/glycam/06j_10 && ./Run.glycam clean
	-cd leap/glycam/06j_12SB && ./Run.glycam clean

check:
	(find . -name "*.dif" ! -size 0c -print | \
	while read dif ;\
	do \
	    echo $$dif ;\
	    cat $$dif ;\
	    echo ;\
	done ;\
	)

