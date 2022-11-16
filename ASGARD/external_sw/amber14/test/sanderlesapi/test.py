#!/usr/bin/env python
from __future__ import division

import os
import sys
try:
    import sanderles as sander
    from chemistry.amber.readparm import Rst7
except ImportError:
    print('Could not import sanderles. Skipping test')
    sys.exit(0)

# For Py3 compatibility
try:
    xrange
except NameError:
    xrange = range

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

def print_result(failed):
    if failed:
        print("Possible FAILURE")
    else:
        print("PASSED")
    print("="*62)

# Run the first test
print('Testing GB sander interface (diffcoords w/ RDT)')
options = sander.gas_input(7)
options.cut = 9999.0
options.rgbmax = 100.0
options.rdt = 0.01
rst = Rst7.open("../LES_GB/les.diffcoords.r")
sander.setup("../LES_GB/les.prm", rst.coordinates, rst.box, options)
e, f = sander.energy_forces()

failed = compare(e.bond, 16.5749, "Bond")
failed = compare(e.angle, 21.5250, "Angle") or failed
failed = compare(e.dihedral, 35.5749, "Dihedral") or failed
failed = compare(e.vdw_14, 6.4411, "1-4 vdW") or failed
failed = compare(e.elec_14, 140.5502, "1-4 Elec") or failed
failed = compare(e.vdw, -4.6590, "van der Waals") or failed
failed = compare(e.elec, -198.7892, "Electrostatic") or failed
failed = compare(e.gb, -30.4884, "EGB") or failed
failed = compare(e.surf, 0.0, "SASA (GBSA)") or failed

print_result(failed)
assert len(f) == 3 * sander.natom()
sander.cleanup()

failed2 = failed
failed = False

# Run the second test
print('Testing GB sander interface (samecoords w/out RDT)')
options.rdt = 0.0
rst = Rst7.open("../LES_GB/les.samecoords.r")
sander.setup("../LES_GB/les.alt.prm", rst.coordinates, rst.box, options)
e, f = sander.energy_forces()

failed = compare(e.bond, 5.8375, "Bond")
failed = compare(e.angle, 19.0846, "Angle") or failed
failed = compare(e.dihedral, 32.7197, "Dihedral") or failed
failed = compare(e.vdw_14, 7.1039, "1-4 vdW") or failed
failed = compare(e.elec_14, 141.3377, "1-4 Elec") or failed
failed = compare(e.vdw, -3.0346, "van der Waals") or failed
failed = compare(e.elec, -202.2822, "Electrostatic") or failed
failed = compare(e.gb, -28.2130, "EGB") or failed

print_result(failed)

failed2 = failed
failed = False
sander.cleanup()

# Run the 3rd test
print('Testing PME sander interface')
options = sander.pme_input()
options.cut = 8.0
rst = Rst7.open("../LES/md.LES.x")
sander.setup("../LES/LES.prmtop.save", rst.coordinates, rst.box, options);
e, f = sander.energy_forces()

failed = compare(e.bond, 14.7095, "Bond")
failed = compare(e.angle, 34.6208, "Angle") or failed
failed = compare(e.dihedral, 35.3483, "Dihedral") or failed
failed = compare(e.vdw_14, 13.0097, "1-4 vdW") or failed
failed = compare(e.elec_14, 274.1453, "1-4 Elec") or failed
failed = compare(e.vdw, 545.6397, "van der Waals") or failed
failed = compare(e.elec, -4995.3084, "Electrostatic") or failed

print_result(failed)

assert len(f) == 3 * sander.natom()

failed2 = failed2 or failed
failed = False
sander.cleanup()
