These tests check that CHAMBER and SANDER are correctly dealing with periodic 
boundary conditions within CHARMM restart files. 

The test system is a box of 3000 TIP3P waters that has been heated up to
~300K within charmm and a restart file saved of this state. See the ./build
directory for how this was done.

The comp_ene awk script in this directory extracts the various energy terms
out of the charmm and SANDER outputs for each test and parses these in a comparable
manner.

*Test 1 has no box information, hence is just a cube of water in vacuum.

*Test 2 has box information, however, PME is not used to evalulate its electrostatic 
energy. No switch function is used (eedmeth=4) is used for the direct sum Coulomb
interaction.

*Test 3 is the same as test 2, but PME is used. Identical PME settings are
used in both charmm and SANDER; the FFT grid and spline order needed to be
increased achieve compariable results.

No correction for vdw after cut off (vdwmeth=0) is set for all three tests,
but only has an effect when PME is on.


Known issues
============
This testcase is still work in progress since there are still minor
differences in the kinetic energy of the same restart file between charmm and
SANDER.

The results need to be parsed via DACDIF to give a discreet pass/fail result,
currently the result is dumped to stdout via the awk script.


Typical Output
==============

==============
Running test 1
==============
awk -f comp_ene.awk test1.out mdout.test1 mdinfo.test1
            amb         chm    difference
                               (amb - chm)
ELEC  -11682.0388 -11682.0385     -0.0003       amber:
                                         elec   -11682.0388
                                         ee14        0.0000
VDW     1288.0148   1288.0148      0.0000       amber:
                                         vdw      1288.0148
                                         vdw14       0.0000
EWALD
 self       0.0000      0.0000      0.0000
 rec        0.0000      0.0000      0.0000
 dir   -11682.0388 -11682.0385     -0.0003
 adj        0.0000      0.0000      0.0000

awk -f comp_ektot.awk test1.out mdinfo.test1
            amb         chm    difference
                               (amb - chm)
 Ektot     1768.4612   1765.5368      2.9244

==============
Running test 2
==============
awk -f comp_ene.awk test2.out mdout.test2 mdinfo.test2
            amb         chm    difference
                               (amb - chm)
ELEC  -18675.0800 -18675.0796     -0.0004       amber:
                                         elec   -18675.0800
                                         ee14        0.0000
VDW     1311.5689   1311.5689      0.0000       amber:
                                         vdw      1311.5689
                                         vdw14       0.0000
EWALD
 self       0.0000      0.0000      0.0000
 rec        0.0000      0.0000      0.0000
 dir   -18675.0800 -18675.0796     -0.0004
 adj        0.0000      0.0000      0.0000

awk -f comp_ektot.awk test2.out  mdinfo.test2
            amb         chm    difference
                               (amb - chm)
 Ektot     1773.5671   1773.2612      0.3059

==============
Running test 3
==============
awk -f comp_ene.awk test3.out mdout.test3 mdinfo.test3
            amb         chm    difference
                               (amb - chm)
ELEC  -11359.0028 -11359.0025     -0.0003       amber:
                                         elec   -11359.0028
                                         ee14        0.0000
VDW     1286.0903   1286.0903      0.0000       amber:
                                         vdw      1286.0903
                                         vdw14       0.0000
EWALD
 self  -60204.7678 -60204.7663     -0.0015
 rec       51.3174     51.3173      0.0000
 dir   -10508.2689 -10508.2686     -0.0003
 adj    59302.7165  59302.7151      0.0014

awk -f comp_ektot.awk test3.out mdinfo.test3
            amb         chm    difference
                               (amb - chm)
 Ektot     1773.6234   1773.6218      0.0016

