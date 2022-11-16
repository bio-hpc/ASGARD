/* This program tests the sander C API by running the same tests as the Fortran
 * interface.
 */

#include "sander.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;

bool compare(const double, const double, const char*);
bool compare_forces(const double, const double);
bool slko_files_exist(void);

int main() {

    double *forces, *coordinates;
    double box[6];
    int natom;
    bool failed, failed2;

    sander_input options;
    qmmm_input_options qmmm_options;
    pot_ene energies;

    failed = false;
    
    cout << "Testing GB sander interface (diffcoords w/ RDT)" << endl;
    gas_sander_input(&options, 7);
    options.cut = 9999.0;
    options.rgbmax = 100.0;
    options.rdt = 0.01;

    natom = get_incprd_natom("../LES_GB/les.diffcoords.r");
    coordinates = new double[natom*3];
    read_inpcrd_file("../LES_GB/les.diffcoords.r", coordinates, box);
    sander_setup_mm("../LES_GB/les.prm", coordinates, box, &options);
    forces = new double[3*sander_natom()];

    energy_forces(&energies, forces);

/* Compare the energies to the output from the relevant test:
 NSTEP =        0   TIME(PS) = 1000010.000  TEMP(K) =     0.00  PRESS =     0.0
 Etot   =       -13.2705  EKtot   =         0.0000  EPtot      =       -13.2705
 BOND   =        16.5749  ANGLE   =        21.5250  DIHED      =        35.5749
 1-4 NB =         6.4411  1-4 EEL =       140.5502  VDWAALS    =        -4.6590
 EELEC  =      -198.7892  EGB     =       -30.4884  RESTRAINT  =         0.0000
*/

    failed = compare(energies.bond, 16.5749, "Bond");
    failed = compare(energies.angle, 21.5250, "Angle") || failed;
    failed = compare(energies.dihedral, 35.5749, "Dihedral") || failed;
    failed = compare(energies.vdw_14, 6.4411, "1-4 vdW") || failed;
    failed = compare(energies.elec_14, 140.5502, "1-4 Elec") || failed;
    failed = compare(energies.vdw, -4.6590, "van der Waals") || failed;
    failed = compare(energies.elec, -198.7892, "Electrostatic") || failed;
    failed = compare(energies.gb, -30.4884, "EGB") || failed;
    failed = compare(energies.surf, 0.0, "SASA (GBSA)") || failed;

    sander_cleanup();
    delete[] forces;
    delete[] coordinates;

    if (failed) {
        cout << "Possible FAILURE" << endl;
    } else {
        cout << "PASSED" << endl;
    }
    for (int i = 0; i < 62; i++) cout << "="; cout << endl;

    cout << "Testing GB sander interface (samecoords w/out RDT)" << endl;
    failed2 = failed;
    failed = false;
    options.rdt = 0.0;

    natom = get_inpcrd_natom("../LES_GB/les.samecoords.r");
    coordinates = new double[3*natom];
    read_inpcrd_file("../LES_GB/les.samecoords.r", coordinates, box);
    sander_setup_mm("../LES_GB/les.alt.prm", coordinates, box, &options);
    forces  = new double[3*natom];

    energy_forces(&energies, forces);

/* Compare the energies to the output from the relevant test:
 NSTEP =        0   TIME(PS) = 1000010.000  TEMP(K) =     0.00  PRESS =     0.0
 Etot   =       -27.4464  EKtot   =         0.0000  EPtot      =       -27.4464
 BOND   =         5.8375  ANGLE   =        19.0846  DIHED      =        32.7197
 1-4 NB =         7.1039  1-4 EEL =       141.3377  VDWAALS    =        -3.0346
 EELEC  =      -202.2822  EGB     =       -28.2130  RESTRAINT  =         0.0000
*/

    failed = compare(energies.bond, 5.8375, "Bond") || failed;
    failed = compare(energies.angle, 19.0846, "Angle") || failed;
    failed = compare(energies.dihedral, 32.7197, "Dihedral") || failed;
    failed = compare(energies.vdw_14, 7.1039, "1-4 vdW") || failed;
    failed = compare(energies.elec_14, 141.3377, "1-4 Elec") || failed;
    failed = compare(energies.vdw, -3.0346, "van der Waals") || failed;
    failed = compare(energies.elec, -202.2822, "Electrostatic") || failed;
    failed = compare(energies.gb, -28.2130, "EGB") || failed;

    sander_cleanup();
    delete[] forces;
    delete[] coordinates;

    if (failed > 0) {
        cout << "Possible FAILURE" << endl;
    } else {
        cout << "PASSED" << endl;
    }
    for (int i = 0; i < 62; i++) cout << "="; cout << endl;

    cout << "Testing PME sander interface" << endl;
    failed2 = failed || failed2;
    failed = false;
    pme_sander_input(&options);
    options.cut = 8.0;

    natom = get_inpcrd_natom("../LES/md.LES.x");
    coordinates = new double[3*natom];
    read_inpcrd_file("../LES/md.LES.x", coordinates, box);
    sander_setup_mm("../LES/LES.prmtop.save", coordinates, box, &options);
    forces = new double[3*sander_natom()];
    energy_forces(&energies, forces);

/* Compare the energies to the output from the relevant test:
 NSTEP =        1   TIME(PS) =      10.002  TEMP(K) =   289.54  PRESS =     0.0
 Etot   =     -3114.6572  EKtot   =       963.1779  EPtot      =     -4077.8351
 BOND   =        14.7095  ANGLE   =        34.6208  DIHED      =        35.3483
 1-4 NB =        13.0097  1-4 EEL =       274.1453  VDWAALS    =       545.6397
 EELEC  =     -4995.3084  EHBOND  =         0.0000  RESTRAINT  =         0.0000
*/
    failed = compare(energies.bond, 14.7095, "Bond");
    failed = compare(energies.angle, 34.6208, "Angle") || failed;
    failed = compare(energies.dihedral, 35.3483, "Dihedral") || failed;
    failed = compare(energies.vdw_14, 13.0097, "1-4 vdW") || failed;
    failed = compare(energies.elec_14, 274.1453, "1-4 Elec") || failed;
    failed = compare(energies.vdw, 545.6397, "van der Waals") || failed;
    failed = compare(energies.elec, -4995.3084, "Electrostatic") || failed;

    if (failed > 0) {
        printf("Possible FAILURE\n");
    } else {
        printf("PASSED\n");
    }
    for (int i = 0; i < 62; i++) cout << "="; cout << endl;
    delete[] forces;
    sander_cleanup();

    // Exit with the appropriate exit code
    if (failed || failed2) return 1;
    return 0;
}

bool compare(const double computed, const double regression, const char* desc) {

    // Compare to 4 decimal places

    if (fabs(computed - regression) > 2.0e-4) {
        printf("%s failed: Expected %15.4f got %15.4f\n",
               desc, regression, computed);
        return true;
    }
    return false;
}

bool compare_forces(const double computed, const double regression) {

    // Compare to 4 decimal places

    if (fabs(computed - regression) > 3.0e-4) {
        cout << "Force comparison failed: " << regression << " vs. "
             << computed << " (" << fabs(computed-regression) << ")" << endl;
        return true;
    }

    return false;
}

bool slko_files_exist(void) {

    if (ifstream("../../dat/slko/C-C.skf", ifstream::in))
        return true;
    return false;
}
