#!/usr/bin/env python
from __future__ import division

import os
import sys
try:
    import sander
    from chemistry.amber.readparm import Rst7
except ImportError:
    print('Could not import sander. Skipping test')
    sys.exit(0)

# For Py3 compatibility
try:
    xrange
except NameError:
    xrange = range

# For duck-typing test with set_positions
class Vec3(tuple):

    def __new__(cls, a, b, c):
        return tuple.__new__(cls, (a, b, c))

def compare(computed, regression, desc):
    """
    Compares expected and computed values. Returns True if the comparison
    failed, False otherwise
    """
    if abs(computed - regression) > 2e-4:
        print("%s failed: Expected %15.4f got %15.4f" %
              (desc, regression, computed))
        return True
    return False

def compare_forces(computed, regression):
    """ Compares two forces for the same components """
    if abs(computed - regression) > 3e-4:
        print('Force comparison failed: %s vs. %s (%s)' % (computed, regression,
            abs(computed-regression)))
        return True
    return False

def print_result(failed):
    if failed:
        print("Possible FAILURE")
    else:
        print("PASSED")
    print("="*62)

# Run the first test
print('Testing GB sander interface')
options = sander.gas_input(7)
options.cut = 9999.0
options.saltcon = 0.2
options.gbsa = 1
rst = Rst7.open("../gb7_trx/trxox.2.4ns.x")
sander.setup("../gb7_trx/prmtop_an", rst.coordinates, rst.box, options)
e, f = sander.energy_forces()

failed = compare(e.bond, 631.8993, 'Bond')
failed = compare(e.angle, 898.2543, 'Angle') or failed
failed = compare(e.dihedral, 566.4453, 'Dihedral') or failed
failed = compare(e.vdw_14, 348.8246, '1-4 vdW') or failed
failed = compare(e.elec_14, 5980.5047, '1-4 Elec') or failed
failed = compare(e.vdw, -768.3629, 'van der Waals') or failed
failed = compare(e.elec, -7874.4913, 'Electrostatic') or failed
failed = compare(e.gb, -1943.0838, 'EGB') or failed
failed = compare(e.surf, 33.8338, 'SASA (GBSA)') or failed

print_result(failed)
assert len(f) == 3 * sander.natom()
sander.cleanup()

failed2 = failed
failed = False

# Run the second test
print('Testing PME sander interface')
options = sander.pme_input()
options.cut = 8.0
rst = Rst7.open('../4096wat/eq1.x')
sander.setup("../4096wat/prmtop", rst.coordinates, rst.box, options)
e, f = sander.energy_forces()

failed = compare(e.bond, 0.0, "Bond") or failed
failed = compare(e.angle, 0.0, "Angle") or failed
failed = compare(e.dihedral, 0.0, "Dihedral") or failed
failed = compare(e.vdw_14, 0.0, "1-4 vdW") or failed
failed = compare(e.elec_14, 0.0, "1-4 Elec") or failed
failed = compare(e.vdw, 6028.9517, "van der Waals") or failed
failed = compare(e.elec, -45371.5995, "Electrostatic") or failed

print_result(failed)

assert len(f) == 3 * sander.natom()

failed2 = failed2 or failed
failed = False

# Run the third test
print('Checking for consistent forces')
fi = open("../4096wat/mdfrc_cmp.save", 'r')
try:
    line = fi.readline()
    while line and not line.strip().startswith('forces ='):
        line = fi.readline()

    for i in xrange(sander.natom()):
        i3 = i * 3
        frcs = [float(x) for x in fi.readline().split()]
        failed = failed or compare_forces(f[i3  ], frcs[0])
        failed = failed or compare_forces(f[i3+1], frcs[1])
        failed = failed or compare_forces(f[i3+2], frcs[2])
finally:
    fi.close()

print_result(failed)
sander.cleanup()

failed2 = failed2 or failed
failed = False

# Run the 4th test
print('Testing the QM/MM non-periodic interface')
options = sander.gas_input(1)
options.cut = 99.0
options.ifqnt = 1
qmmm_options = sander.qm_input()
qmmm_options.iqmatoms[:3] = [8, 9, 10]
qmmm_options.qm_theory = "PM3"
qmmm_options.qmcharge = 0
qmmm_options.qmgb = 2
qmmm_options.adjust_q = 0

rst = Rst7.open("../qmmm2/lysine_PM3_qmgb2/lysine.crd")
sander.setup("../qmmm2/lysine_PM3_qmgb2/prmtop",
             rst.coordinates, rst.box, options, qmmm_options)
e, f = sander.energy_forces()

failed = compare(e.bond, 0.0016, "Bond") or failed
failed = compare(e.angle, 0.3736, "Angle") or failed
failed = compare(e.dihedral, 0.0026, "Dihedral") or failed
failed = compare(e.vdw_14, 3.7051, "1-4 vdW") or failed
failed = compare(e.elec_14, 65.9137, "1-4 Elec") or failed
failed = compare(e.vdw, 0.1908, "van der Waals") or failed
failed = compare(e.elec, -4.1241, "Electrostatic") or failed
failed = compare(e.gb, -80.1406, "EGB") or failed
failed = compare(e.scf, -11.9100, "QM Escf") or failed


print_result(failed)
sander.cleanup()

failed2 = failed2 or failed
failed = False

# Run the 5th test
print('Testing the QM/MM periodic interface (PM3-PDDG)')
options = sander.pme_input()
options.cut = 8.0
options.ifqnt = 1
options.jfastw = 4

qmmm_options = sander.QmInputOptions()
qmmm_options.qm_theory = "PDDG-PM3"
qmmm_options.qmmask = ":1-2"
qmmm_options.qmcharge = 0
qmmm_options.scfconv = 1e-10
qmmm_tight_p_conv = 1
qmmm_options.qmmm_int = 5

rst = Rst7.open("../qmmm2/MechEm_nma-spcfwbox/inpcrd")
sander.setup("../qmmm2/MechEm_nma-spcfwbox/prmtop",
             rst.coordinates, rst.box, options, qmmm_options)
e, f = sander.energy_forces()

failed = compare(e.bond, 605.7349, "Bond") or failed
failed = compare(e.angle, 331.7679, "Angle") or failed
failed = compare(e.dihedral, 0.0000, "Dihedral") or failed
failed = compare(e.vdw_14, 0.0000, "1-4 vdW") or failed
failed = compare(e.elec_14, 0.0000, "1-4 Elec") or failed
failed = compare(e.vdw, 1281.8450, "van der Waals") or failed
failed = compare(e.elec, -7409.7167, "Electrostatic") or failed
failed = compare(e.scf, -37.1277, "QM Escf") or failed

print_result(failed)
sander.cleanup()

failed2 = failed2 or failed
failed = False

# Run 6th test
print('Testing the QM/MM periodic interface (DFTB)')
if not os.path.exists('../../dat/slko/C-C.skf'):
    print('Could not find the SLKO files. Skipping this test.')
else:
    qmmm_options.qm_theory = "DFTB"
    qmmm_options.dftb_3rd_order = "NONE"
    sander.setup("../qmmm2/MechEm_nma-spcfwbox/prmtop",
                 rst.coordinates, rst.box, options, qmmm_options)
    e, f = sander.energy_forces()
    failed = compare(e.bond, 605.7349, "Bond") or failed
    failed = compare(e.angle, 331.7679, "Angle") or failed
    failed = compare(e.dihedral, 0.0000, "Dihedral") or failed
    failed = compare(e.vdw_14, 0.0000, "1-4 vdW") or failed
    failed = compare(e.elec_14, 0.0000, "1-4 Elec") or failed
    failed = compare(e.vdw, 1281.8450, "van der Waals") or failed
    failed = compare(e.elec, -7409.7167, "Electrostatic") or failed
    failed = compare(e.scf, -1209.0254, "QM Escf") or failed
    
    print_result(failed)
    sander.cleanup()

# Run the 6th test
print('Testing the broader API functionality')

failed2 = failed or failed2
failed = False

options = sander.gas_input(7)
options.cut = 9999.0
options.saltcon = 0.2
options.gbsa = 1
rst = Rst7.open("../gb7_trx/trxox.2.4ns.x")
sander.setup("../gb7_trx/prmtop_an", rst.coordinates, rst.box, options)
e, f = sander.energy_forces()

failed = compare(e.bond, 631.8993, 'Bond')
failed = compare(e.angle, 898.2543, 'Angle') or failed
failed = compare(e.dihedral, 566.4453, 'Dihedral') or failed
failed = compare(e.vdw_14, 348.8246, '1-4 vdW') or failed
failed = compare(e.elec_14, 5980.5047, '1-4 Elec') or failed
failed = compare(e.vdw, -768.3629, 'van der Waals') or failed
failed = compare(e.elec, -7874.4913, 'Electrostatic') or failed
failed = compare(e.gb, -1943.0838, 'EGB') or failed
failed = compare(e.surf, 33.8338, 'SASA (GBSA)') or failed

inpcrd = Rst7.open('../gb7_trx/trxox.2.4pns.x')
sander.set_positions(inpcrd.coordinates)
e, f = sander.energy_forces()

failed = compare(e.bond, 342.0589, "Bond") or failed
failed = compare(e.angle, 877.9924, "Angle") or failed
failed = compare(e.dihedral, 580.6551, "Dihedral") or failed
failed = compare(e.vdw_14, 357.5908, "1-4 vdW") or failed
failed = compare(e.elec_14, 5973.9713, "1-4 Elec") or failed
failed = compare(e.vdw, -770.7973, "van der Waals") or failed
failed = compare(e.elec, -7883.3799, "Electrostatic") or failed
failed = compare(e.gb, -1938.5736, "EGB") or failed
failed = compare(e.surf, 33.8944, "SASA (GBSA)") or failed

sander.cleanup()

options = sander.pme_input()

rst = Rst7.open('../4096wat/eq1.x')
sander.setup('../4096wat/prmtop', rst.coordinates, rst.box, options)
e, f = sander.energy_forces()

failed = compare(e.bond, 0.0, "Bond") or failed
failed = compare(e.angle, 0.0, "Angle") or failed
failed = compare(e.dihedral, 0.0, "Dihedral") or failed
failed = compare(e.vdw_14, 0.0, "1-4 vdW") or failed
failed = compare(e.elec_14, 0.0, "1-4 Elec") or failed
failed = compare(e.vdw, 6028.9517, "van der Waals") or failed
failed = compare(e.elec, -45371.5995, "Electrostatic") or failed

sander.set_box(51, 51, 51, 90, 90, 90)

e, f = sander.energy_forces()

failed = compare(e.bond, 0., "Bond") or failed
failed = compare(e.angle, 0., "Angle") or failed
failed = compare(e.dihedral, 0., "Dihedral") or failed
failed = compare(e.vdw_14, 0., "1-4 vdW") or failed
failed = compare(e.elec_14, 0., "1-4 Elec") or failed
failed = compare(e.vdw, 12567.4552, "van der Waals") or failed
failed = compare(e.elec, -42333.7294, "Electrostatic") or failed

print_result(failed)

failed = False
print('Checking get_positions call')
coords = []
for i in range(sander.natom()):
    i3 = i * 3
    coords.append(Vec3(rst.coordinates[i3],
                       rst.coordinates[i3+1],
                       rst.coordinates[i3+2])
    )
sander.set_positions(coords)

for x, y in zip(rst.coordinates, sander.get_positions()):
    failed = compare(x, y, "Coordinates") or failed

print_result(failed)
