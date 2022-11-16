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

    double *forces;
    bool failed;
    double box[6];
    double *coordinates;
    int natom;

    sander_input options;
    qmmm_input_options qmmm_options;
    pot_ene energies;

    failed = false;
    
    cout << "Testing GB sander interface" << endl;
    gas_sander_input(&options, 7);
    options.cut = 9999.0;
    options.saltcon = 0.2;
    options.gbsa = 1;

    natom = get_inpcrd_natom("../gb7_trx/trxox.2.4ns.x");
    coordinates = new double[3*natom];
    read_inpcrd_file("../gb7_trx/trxox.2.4ns.x", coordinates, box);
    sander_setup_mm("../gb7_trx/prmtop_an", coordinates, box, &options);
    forces = new double[3*sander_natom()];

    energy_forces(&energies, forces);

/* Compare the energies to the output from the relevant test:
 BOND    =      631.8993  ANGLE   =      898.2543  DIHED      =      566.4453
 VDWAALS =     -768.3629  EEL     =    -7874.4913  EGB        =    -1943.0838
 1-4 VDW =      348.8246  1-4 EEL =     5980.5047  RESTRAINT  =        0.0000
 ESURF   =       33.8338
*/
    failed = compare(energies.bond, 631.8993, "Bond");
    failed = compare(energies.angle, 898.2543, "Angle") || failed;
    failed = compare(energies.dihedral, 566.4453, "Dihedral") || failed;
    failed = compare(energies.vdw_14, 348.8246, "1-4 vdW") || failed;
    failed = compare(energies.elec_14, 5980.5047, "1-4 Elec") || failed;
    failed = compare(energies.vdw, -768.3629, "van der Waals") || failed;
    failed = compare(energies.elec, -7874.4913, "Electrostatic") || failed;
    failed = compare(energies.gb, -1943.0838, "EGB") || failed;
    failed = compare(energies.surf, 33.8338, "SASA (GBSA)") || failed;

    sander_cleanup();

    if (failed) {
        cout << "Possible FAILURE" << endl;
    } else {
        cout << "PASSED" << endl;
    }
    for (int i = 0; i < 62; i++) cout << "="; cout << endl;

    cout << "Testing PME sander interface" << endl;
    failed = false;
    pme_sander_input(&options);
    options.cut = 8.0;

    delete[] forces;
    delete[] coordinates;

    natom = get_inpcrd_natom("../4096wat/eq1.x");
    coordinates = new double[3*natom];
    read_inpcrd_file("../4096wat/eq1.x", coordinates, box);
    sander_setup_mm("../4096wat/prmtop", coordinates, box, &options);
    forces = new double[3*sander_natom()];

    energy_forces(&energies, forces);

/* Compare the energies to the output from the relevant test:
 NSTEP =        1   TIME(PS) =       1.001  TEMP(K) =   298.28  PRESS =     0.0
 Etot   =    -32059.8471  EKtot   =      7282.8008  EPtot      =    -39342.6479
 BOND   =         0.0000  ANGLE   =         0.0000  DIHED      =         0.0000
 1-4 NB =         0.0000  1-4 EEL =         0.0000  VDWAALS    =      6028.9517
 EELEC  =    -45371.5995  EHBOND  =         0.0000  RESTRAINT  =         0.0000
*/
    failed = compare(energies.bond, 0.0, "Bond") || failed;
    failed = compare(energies.angle, 0.0, "Angle") || failed;
    failed = compare(energies.dihedral, 0.0, "Dihedral") || failed;
    failed = compare(energies.vdw_14, 0.0, "1-4 vdW") || failed;
    failed = compare(energies.elec_14, 0.0, "1-4 Elec") || failed;
    failed = compare(energies.vdw, 6028.9517, "van der Waals") || failed;
    failed = compare(energies.elec, -45371.5995, "Electrostatic") || failed;

    if (failed) {
        cout << "Possible FAILURE" << endl;
    } else {
        cout << "PASSED" << endl;
    }
    for (int i = 0; i < 62; i++) cout << "="; cout << endl;

    cout << "Checking for consistent forces" << endl;
    failed = false;

    // Now open up the saved file with all of the forces and check those
    ifstream frcfile("../4096wat/mdfrc_cmp.save", ifstream::in);

    string line;
    while (getline(frcfile, line)) {
        if (line != " forces =")
            continue;
        // Now we start parsing our forces
        for (int i = 0; i < 12288; i++) {
            getline(frcfile, line);
            double x, y, z;
            sscanf(line.c_str(), "%lf %lf %lf", &x, &y, &z);
            int i3 = i * 3;
            if (i == 43452) continue;
            failed = failed || compare_forces(forces[i3  ], x);
            failed = failed || compare_forces(forces[i3+1], y);
            failed = failed || compare_forces(forces[i3+2], z);
        }
        // If we got here, we're done parsing
        break;
    }
    delete[] forces;
    delete[] coordinates;
    // Close the frcfile
    frcfile.close();

    if (failed) {
        cout << "Possible FAILURE" << endl;
    } else {
        cout << "PASSED" << endl;
    }
    for (int i = 0; i < 62; i++) cout << "="; cout << endl;
    sander_cleanup();

    cout << "Testing the QM/MM non-periodic interface" << endl;
    gas_sander_input(&options, 1);
    options.cut = 99.0;
    options.ifqnt = 1;
    qm_sander_input(&qmmm_options);
    qmmm_options.iqmatoms[0] = 8;  // Ugh, C and C++ do not accept slicing...
    qmmm_options.iqmatoms[1] = 9;
    qmmm_options.iqmatoms[2] = 10;
    strncpy(qmmm_options.qm_theory, "PM3", 3);
    qmmm_options.qmcharge = 0;
    qmmm_options.qmgb = 2;
    qmmm_options.adjust_q = 0;

    natom = get_inpcrd_natom("../qmmm2/lysine_PM3_qmgb2/lysine.crd");
    coordinates = new double[3*natom];
    read_inpcrd_file("../qmmm2/lysine_PM3_qmgb2/lysine.crd", coordinates, box);
    sander_setup("../qmmm2/lysine_PM3_qmgb2/prmtop",
                 coordinates, box, &options, &qmmm_options);
    forces = new double[3*sander_natom()];

    energy_forces(&energies, forces);

    failed = compare(energies.bond, 0.0016, "Bond") || failed;
    failed = compare(energies.angle, 0.3736, "Angle") || failed;
    failed = compare(energies.dihedral, 0.0026, "Dihedral") || failed;
    failed = compare(energies.vdw_14, 3.7051, "1-4 vdW") || failed;
    failed = compare(energies.elec_14, 65.9137, "1-4 Elec") || failed;
    failed = compare(energies.vdw, 0.1908, "van der Waals") || failed;
    failed = compare(energies.elec, -4.1241, "Electrostatic") || failed;
    failed = compare(energies.gb, -80.1406, "EGB") || failed;
    failed = compare(energies.scf, -11.9100, "QM Escf") || failed;
    
    if (failed) {
        cout << "Possible FAILURE" << endl;
    } else {
        cout << "PASSED" << endl;
    }
    for (int i = 0; i < 62; i++) cout << "="; cout << endl;
    delete[] forces;
    delete[] coordinates;
    sander_cleanup();

    cout << "Testing the QM/MM periodic interface (PM3-PDDG)" << endl;

    pme_sander_input(&options);
    options.cut = 8.0;
    options.ifqnt = 1;
    options.jfastw = 4;

    qm_sander_input(&qmmm_options);
    strncpy(qmmm_options.qmmask, ":1-2", 4);
    strncpy(qmmm_options.qm_theory, "PDDG-PM3", 8);
    qmmm_options.qmcharge = 0;
    qmmm_options.scfconv = 1e-10;
    qmmm_options.tight_p_conv = 1;
    qmmm_options.qmmm_int = 5;

    natom = get_inpcrd_natom("../qmmm2/MechEm_nma-spcfwbox/inpcrd");
    coordinates = new double[3*natom];
    read_inpcrd_file("../qmmm2/MechEm_nma-spcfwbox/inpcrd", coordinates, box);
    sander_setup("../qmmm2/MechEm_nma-spcfwbox/prmtop",
                 coordinates, box, &options, &qmmm_options);
    forces = new double[3*sander_natom()];

    energy_forces(&energies, forces);

    failed = compare(energies.bond, 605.7349, "Bond") || failed;
    failed = compare(energies.angle, 331.7679, "Angle") || failed;
    failed = compare(energies.dihedral, 0.0000, "Dihedral") || failed;
    failed = compare(energies.vdw_14, 0.0000, "1-4 vdW") || failed;
    failed = compare(energies.elec_14, 0.0000, "1-4 Elec") || failed;
    failed = compare(energies.vdw, 1281.8450, "van der Waals") || failed;
    failed = compare(energies.elec, -7409.7167, "Electrostatic") || failed;
    failed = compare(energies.scf, -37.1277, "QM Escf") || failed;

    if (failed) {
        cout << "Possible FAILURE" << endl;
    } else {
        cout << "PASSED" << endl;
    }
    for (int i = 0; i < 62; i++) cout << "="; cout << endl;
    sander_cleanup();

    failed = 0;

    cout << "Testing the QM/MM periodic interface (DFTB)" << endl;

    if (slko_files_exist()) {
        qm_sander_input(&qmmm_options);
        strncpy(qmmm_options.qmmask, ":1-2", 4);
        strncpy(qmmm_options.qm_theory, "DFTB", 4);
        qmmm_options.qmcharge = 0;
        qmmm_options.scfconv = 1e-10;
        qmmm_options.tight_p_conv = 1;
        qmmm_options.qmmm_int = 5;
        sander_setup("../qmmm2/MechEm_nma-spcfwbox/prmtop",
                     coordinates, box, &options, &qmmm_options);
        energy_forces(&energies, forces);

        failed = compare(energies.bond, 605.7349, "Bond") || failed;
        failed = compare(energies.angle, 331.7679, "Angle") || failed;
        failed = compare(energies.dihedral, 0.0000, "Dihedral") || failed;
        failed = compare(energies.vdw_14, 0.0000, "1-4 vdW") || failed;
        failed = compare(energies.elec_14, 0.0000, "1-4 Elec") || failed;
        failed = compare(energies.vdw, 1281.8450, "van der Waals") || failed;
        failed = compare(energies.elec, -7409.7167, "Electrostatic") || failed;
        failed = compare(energies.scf, -1209.0254, "QM Escf") || failed;

        if (failed) {
            cout << "Possible FAILURE" << endl;
        } else {
            cout << "PASSED" << endl;
        }
    } else {
        cout << "Could not find the SLKO files. Skipping this test." << endl;
    }
    sander_cleanup();
    for (int i = 0; i < 62; i++) cout << "="; cout << endl;
    delete[] forces;
    delete[] coordinates;

    cout << "Testing the broader API functionality" << endl;

    failed = false;

    gas_sander_input(&options, 7);
    options.cut = 9999.0;
    options.saltcon = 0.2;
    options.gbsa = 1;

    natom = get_inpcrd_natom("../gb7_trx/trxox.2.4ns.x");
    coordinates = new double[3*natom];
    read_inpcrd_file("../gb7_trx/trxox.2.4ns.x", coordinates, box);
    sander_setup_mm("../gb7_trx/prmtop_an", coordinates, box, &options);
    forces = new double[3*sander_natom()];

    energy_forces(&energies, forces);

    failed = compare(energies.bond, 631.8993, "Bond") || failed;
    failed = compare(energies.angle, 898.2543, "Angle") || failed;
    failed = compare(energies.dihedral, 566.4453, "Dihedral") || failed;
    failed = compare(energies.vdw_14, 348.8246, "1-4 vdW") || failed;
    failed = compare(energies.elec_14, 5980.5047, "1-4 Elec") || failed;
    failed = compare(energies.vdw, -768.3629, "van der Waals") || failed;
    failed = compare(energies.elec, -7874.4913, "Electrostatic") || failed;
    failed = compare(energies.gb, -1943.0838, "EGB") || failed;
    failed = compare(energies.surf, 33.8338, "SASA (GBSA)") || failed;

    natom = get_inpcrd_natom("../gb7_trx/trxox.2.4pns.x");

    if (natom != 1654) {
        failed = true;
    } else {
        read_inpcrd_file("../gb7_trx/trxox.2.4pns.x", forces, box);
        failed = compare(forces[   0], -12.6153818, "Coordinate") || failed;
        failed = compare(forces[   1],   7.5431308, "Coordinate") || failed;
        failed = compare(forces[   2],  -3.4102499, "Coordinate") || failed;
        failed = compare(forces[   3], -12.7604696, "Coordinate") || failed;
        failed = compare(forces[   4],   7.4844904, "Coordinate") || failed;
        failed = compare(forces[   5],  -4.4080529, "Coordinate") || failed;
        failed = compare(forces[4959],  -0.3058399, "Coordinate") || failed;
        failed = compare(forces[4960], -15.0993313, "Coordinate") || failed;
        failed = compare(forces[4961],  13.7839947, "Coordinate") || failed;

        set_positions(forces);
        energy_forces(&energies, forces);

        failed = compare(energies.bond, 342.0589, "Bond") || failed;
        failed = compare(energies.angle, 877.9924, "Angle") || failed;
        failed = compare(energies.dihedral, 580.6551, "Dihedral") || failed;
        failed = compare(energies.vdw_14, 357.5908, "1-4 vdW") || failed;
        failed = compare(energies.elec_14, 5973.9713, "1-4 Elec") || failed;
        failed = compare(energies.vdw, -770.7973, "van der Waals") || failed;
        failed = compare(energies.elec, -7883.3799, "Electrostatic") || failed;
        failed = compare(energies.gb, -1938.5736, "EGB") || failed;
        failed = compare(energies.surf, 33.8944, "SASA (GBSA)") || failed;
    }
    delete[] forces;
    delete[] coordinates;
    sander_cleanup();
    
    pme_sander_input(&options);
    options.cut = 8.0;

    natom = get_inpcrd_natom("../4096wat/eq1.x");
    coordinates = new double[natom*3];
    read_inpcrd_file("../4096wat/eq1.x", coordinates, box);
    sander_setup_mm("../4096wat/prmtop", coordinates, box, &options);
    forces = new double[sander_natom()*3];

    energy_forces(&energies, forces);

    failed = compare(energies.bond, 0.0, "Bond") || failed;
    failed = compare(energies.angle, 0.0, "Angle") || failed;
    failed = compare(energies.dihedral, 0.0, "Dihedral") || failed;
    failed = compare(energies.vdw_14, 0.0, "1-4 vdW") || failed;
    failed = compare(energies.elec_14, 0.0, "1-4 Elec") || failed;
    failed = compare(energies.vdw, 6028.9517, "van der Waals") || failed;
    failed = compare(energies.elec, -45371.5995, "Electrostatic") || failed;

    set_box(51.0, 51.0, 51.0, 90.0, 90.0, 90.0);

    energy_forces(&energies, forces);

    failed = compare(energies.bond, 0., "Bond") || failed;
    failed = compare(energies.angle, 0., "Angle") || failed;
    failed = compare(energies.dihedral, 0., "Dihedral") || failed;
    failed = compare(energies.vdw_14, 0., "1-4 vdW") || failed;
    failed = compare(energies.elec_14, 0., "1-4 Elec") || failed;
    failed = compare(energies.vdw, 12567.4552, "van der Waals") || failed;
    failed = compare(energies.elec, -42333.7294, "Electrostatic") || failed;

    if (failed) {
        cout << "Possible FAILURE" << endl;
    } else {
        cout << "PASSED" << endl;
    }
    for (int i = 0; i < 62; i++) cout << "="; cout << endl;

    // Check get_positions
    failed = false;
    cout << "Checking get_positions call" << endl;
    get_positions(forces);

    for (int i = 0; i < sander_natom() * 3; i++) {
        failed = compare(forces[i], coordinates[i], "Coordinates") || failed;
    }

    if (failed) {
        cout << "Possible FAILURE" << endl;
    } else {
        cout << "PASSED" << endl;
    }
    for (int i = 0; i < 62; i++) cout << "="; cout << endl;

    sander_cleanup();
    delete[] forces;
    delete[] coordinates;

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
