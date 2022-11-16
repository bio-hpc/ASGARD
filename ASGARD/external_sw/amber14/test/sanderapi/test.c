/* This program tests the sander C API by running the same tests as the Fortran
 * interface.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sander.h"

int compare(const double, const double, const char*);
int slko_files_exist(void);

int main() {

    double *forces;
    double *coordinates;
    double box[6];
    int failed, i;
    int natom;

    sander_input options;
    qmmm_input_options qmmm_options;
    pot_ene energies;

    failed = 0;
    gas_sander_input(&options, 7);
    options.cut = 9999.0;
    options.saltcon = 0.2;
    options.gbsa = 1;

    natom = get_inpcrd_natom("../gb7_trx/trxox.2.4ns.x");
    coordinates = (double *)malloc(3*natom*sizeof(double));
    read_inpcrd_file("../gb7_trx/trxox.2.4ns.x", coordinates, box);
    printf("Testing proper treatment of a bad prmtop\n");
    if (sander_setup_mm("test.parm7", coordinates, box, &options) == 0)
        failed = 1;
    else
        failed = 0;
    failed += is_setup();
    if (failed > 0) {
        printf("Possible FAILURE\n");
        for (i = 0; i < 62; i++) printf("="); printf("\n");
        return 1;
    } else {
        printf("PASSED\n");
    }
    for (i = 0; i < 62; i++) printf("="); printf("\n");

    printf("Testing proper treatment of bad input\n");
    options.cut = -5.0;
    if (sander_setup_mm("../gb7_trx/prmtop_an", coordinates, box, &options) == 0)
        failed = 1;
    else
        failed = 0;
    failed += is_setup();

    if (failed > 0) {
        printf("Possible FAILURE\n");
        for (i = 0; i < 62; i++) printf("="); printf("\n");
        return 1;
    } else {
        printf("PASSED\n");
    }
    for (i = 0; i < 62; i++) printf("="); printf("\n");

    printf("Testing GB sander interface\n");
    options.cut = 9999.0;
    sander_setup_mm("../gb7_trx/prmtop_an", coordinates, box, &options);
    forces = (double*) malloc(3*sander_natom()*sizeof(double));
    energy_forces(&energies, forces);

/* Compare the energies to the output from the relevant test:
 BOND    =      631.8993  ANGLE   =      898.2543  DIHED      =      566.4453
 VDWAALS =     -768.3629  EEL     =    -7874.4913  EGB        =    -1943.0838
 1-4 VDW =      348.8246  1-4 EEL =     5980.5047  RESTRAINT  =        0.0000
 ESURF   =       33.8338
*/
    failed += compare(energies.bond, 631.8993, "Bond");
    failed += compare(energies.angle, 898.2543, "Angle");
    failed += compare(energies.dihedral, 566.4453, "Dihedral");
    failed += compare(energies.vdw_14, 348.8246, "1-4 vdW");
    failed += compare(energies.elec_14, 5980.5047, "1-4 Elec");
    failed += compare(energies.vdw, -768.3629, "van der Waals");
    failed += compare(energies.elec, -7874.4913, "Electrostatic");
    failed += compare(energies.gb, -1943.0838, "EGB");
    failed += compare(energies.surf, 33.8338, "SASA (GBSA)");

    sander_cleanup();

    if (failed > 0) {
        printf("Possible FAILURE\n");
    } else {
        printf("PASSED\n");
    }
    for (i = 0; i < 62; i++) printf("="); printf("\n");

    printf("Testing PME sander interface\n");
    failed = 0;
    pme_sander_input(&options);
    options.cut = 8.0;

    free(forces);
    free(coordinates);

    natom = get_inpcrd_natom("../4096wat/eq1.x");
    coordinates = (double *)malloc(natom*3*sizeof(double));
    read_inpcrd_file("../4096wat/eq1.x", coordinates, box);
    sander_setup_mm("../4096wat/prmtop", coordinates, box, &options);
    forces = (double*) malloc(3*sander_natom()*sizeof(double));
    energy_forces(&energies, forces);

/* Compare the energies to the output from the relevant test:
 NSTEP =        1   TIME(PS) =       1.001  TEMP(K) =   298.28  PRESS =     0.0
 Etot   =    -32059.8471  EKtot   =      7282.8008  EPtot      =    -39342.6479
 BOND   =         0.0000  ANGLE   =         0.0000  DIHED      =         0.0000
 1-4 NB =         0.0000  1-4 EEL =         0.0000  VDWAALS    =      6028.9517
 EELEC  =    -45371.5995  EHBOND  =         0.0000  RESTRAINT  =         0.0000
*/
    failed += compare(energies.bond, 0.0, "Bond");
    failed += compare(energies.angle, 0.0, "Angle");
    failed += compare(energies.dihedral, 0.0, "Dihedral");
    failed += compare(energies.vdw_14, 0.0, "1-4 vdW");
    failed += compare(energies.elec_14, 0.0, "1-4 Elec");
    failed += compare(energies.vdw, 6028.9517, "van der Waals");
    failed += compare(energies.elec, -45371.5995, "Electrostatic");

    if (failed > 0) {
        printf("Possible FAILURE\n");
    } else {
        printf("PASSED\n");
    }
    for (i = 0; i < 62; i++) printf("="); printf("\n");
    free(forces);
    free(coordinates);
    sander_cleanup();

    failed = 0;

    // Now test the various QM/MM capabilities

    printf("Testing the QM/MM non-periodic interface\n");
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
    coordinates = (double *)malloc(3*natom*sizeof(double));
    read_inpcrd_file("../qmmm2/lysine_PM3_qmgb2/lysine.crd", coordinates, box);
    sander_setup("../qmmm2/lysine_PM3_qmgb2/prmtop",
                 coordinates, box, &options, &qmmm_options);
    forces = (double *) malloc(3*sander_natom()*sizeof(double));

    energy_forces(&energies, forces);

    failed += compare(energies.bond, 0.0016, "Bond");
    failed += compare(energies.angle, 0.3736, "Angle");
    failed += compare(energies.dihedral, 0.0026, "Dihedral");
    failed += compare(energies.vdw_14, 3.7051, "1-4 vdW");
    failed += compare(energies.elec_14, 65.9137, "1-4 Elec");
    failed += compare(energies.vdw, 0.1908, "van der Waals");
    failed += compare(energies.elec, -4.1241, "Electrostatic");
    failed += compare(energies.gb, -80.1406, "EGB");
    failed += compare(energies.scf, -11.9100, "QM Escf");
    
    if (failed > 0) {
        printf("Possible FAILURE\n");
    } else {
        printf("PASSED\n");
    }
    for (i = 0; i < 62; i++) printf("="); printf("\n");
    free(forces);
    free(coordinates);
    sander_cleanup();

    printf("Testing the QM/MM periodic interface (PM3-PDDG)\n");

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
    coordinates = (double *)malloc(3*natom*sizeof(double));
    read_inpcrd_file("../qmmm2/MechEm_nma-spcfwbox/inpcrd", coordinates, box);
    sander_setup("../qmmm2/MechEm_nma-spcfwbox/prmtop",
                 coordinates, box, &options, &qmmm_options);
    forces = (double *) malloc(3*sander_natom()*sizeof(double));

    energy_forces(&energies, forces);

    failed += compare(energies.bond, 605.7349, "Bond");
    failed += compare(energies.angle, 331.7679, "Angle");
    failed += compare(energies.dihedral, 0.0000, "Dihedral");
    failed += compare(energies.vdw_14, 0.0000, "1-4 vdW");
    failed += compare(energies.elec_14, 0.0000, "1-4 Elec");
    failed += compare(energies.vdw, 1281.8450, "van der Waals");
    failed += compare(energies.elec, -7409.7167, "Electrostatic");
    failed += compare(energies.scf, -37.1277, "QM Escf");

    if (failed > 0) {
        printf("Possible FAILURE\n");
    } else {
        printf("PASSED\n");
    }
    for (i = 0; i < 62; i++) printf("="); printf("\n");
    sander_cleanup();

    failed = 0;

    printf("Testing the QM/MM periodic interface (DFTB)\n");

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

        failed += compare(energies.bond, 605.7349, "Bond");
        failed += compare(energies.angle, 331.7679, "Angle");
        failed += compare(energies.dihedral, 0.0000, "Dihedral");
        failed += compare(energies.vdw_14, 0.0000, "1-4 vdW");
        failed += compare(energies.elec_14, 0.0000, "1-4 Elec");
        failed += compare(energies.vdw, 1281.8450, "van der Waals");
        failed += compare(energies.elec, -7409.7167, "Electrostatic");
        failed += compare(energies.scf, -1209.0254, "QM Escf");

        if (failed > 0) {
            printf("Possible FAILURE\n");
        } else {
            printf("PASSED\n");
        }
        sander_cleanup();
    } else {
        printf("Could not find the SLKO files. Skipping this test.\n");
    }
    for (i = 0; i < 62; i++) printf("="); printf("\n");
    free(forces);

    printf("Testing the broader API functionality\n");

    failed = 0;
    gas_sander_input(&options, 7);
    options.cut = 9999.0;
    options.saltcon = 0.2;
    options.gbsa = 1;

    natom = get_inpcrd_natom("../gb7_trx/trxox.2.4ns.x");
    coordinates = (double *)malloc(3*natom*sizeof(double));
    read_inpcrd_file("../gb7_trx/trxox.2.4ns.x", coordinates, box);
    sander_setup_mm("../gb7_trx/prmtop_an", coordinates, box, &options);

    forces = (double*) malloc(3*sander_natom()*sizeof(double));
    energy_forces(&energies, forces);

/* Compare the energies to the output from the relevant test:
 BOND    =      631.8993  ANGLE   =      898.2543  DIHED      =      566.4453
 VDWAALS =     -768.3629  EEL     =    -7874.4913  EGB        =    -1943.0838
 1-4 VDW =      348.8246  1-4 EEL =     5980.5047  RESTRAINT  =        0.0000
 ESURF   =       33.8338
*/
    failed += compare(energies.bond, 631.8993, "Bond");
    failed += compare(energies.angle, 898.2543, "Angle");
    failed += compare(energies.dihedral, 566.4453, "Dihedral");
    failed += compare(energies.vdw_14, 348.8246, "1-4 vdW");
    failed += compare(energies.elec_14, 5980.5047, "1-4 Elec");
    failed += compare(energies.vdw, -768.3629, "van der Waals");
    failed += compare(energies.elec, -7874.4913, "Electrostatic");
    failed += compare(energies.gb, -1943.0838, "EGB");
    failed += compare(energies.surf, 33.8338, "SASA (GBSA)");

    natom = get_inpcrd_natom("../gb7_trx/trxox.2.4pns.x");

    if (natom != 1654) {
        failed += 1;
    } else {
        read_inpcrd_file("../gb7_trx/trxox.2.4pns.x", coordinates, box);
        failed += compare(coordinates[   0], -12.6153818, "Coordinate");
        failed += compare(coordinates[   1],   7.5431308, "Coordinate");
        failed += compare(coordinates[   2],  -3.4102499, "Coordinate");
        failed += compare(coordinates[   3], -12.7604696, "Coordinate");
        failed += compare(coordinates[   4],   7.4844904, "Coordinate");
        failed += compare(coordinates[   5],  -4.4080529, "Coordinate");
        failed += compare(coordinates[4959],  -0.3058399, "Coordinate");
        failed += compare(coordinates[4960], -15.0993313, "Coordinate");
        failed += compare(coordinates[4961],  13.7839947, "Coordinate");

        set_positions(coordinates);
        energy_forces(&energies, forces);

        failed += compare(energies.bond, 342.0589, "Bond");
        failed += compare(energies.angle, 877.9924, "Angle");
        failed += compare(energies.dihedral, 580.6551, "Dihedral");
        failed += compare(energies.vdw_14, 357.5908, "1-4 vdW");
        failed += compare(energies.elec_14, 5973.9713, "1-4 Elec");
        failed += compare(energies.vdw, -770.7973, "van der Waals");
        failed += compare(energies.elec, -7883.3799, "Electrostatic");
        failed += compare(energies.gb, -1938.5736, "EGB");
        failed += compare(energies.surf, 33.8944, "SASA (GBSA)");
    }
    free(forces);
    free(coordinates);
    sander_cleanup();

    pme_sander_input(&options);
    options.cut = 8.0;

    natom = get_inpcrd_natom("../4096wat/eq1.x");
    coordinates = (double *)malloc(natom*3*sizeof(double));
    read_inpcrd_file("../4096wat/eq1.x", coordinates, box);
    sander_setup_mm("../4096wat/prmtop", coordinates, box, &options);
    forces = (double *)malloc(sander_natom()*3*sizeof(double));

    energy_forces(&energies, forces);

    failed += compare(energies.bond, 0.0, "Bond");
    failed += compare(energies.angle, 0.0, "Angle");
    failed += compare(energies.dihedral, 0.0, "Dihedral");
    failed += compare(energies.vdw_14, 0.0, "1-4 vdW");
    failed += compare(energies.elec_14, 0.0, "1-4 Elec");
    failed += compare(energies.vdw, 6028.9517, "van der Waals");
    failed += compare(energies.elec, -45371.5995, "Electrostatic");

    set_box(51.0, 51.0, 51.0, 90.0, 90.0, 90.0);

    energy_forces(&energies, forces);

    failed += compare(energies.bond, 0., "Bond");
    failed += compare(energies.angle, 0., "Angle");
    failed += compare(energies.dihedral, 0., "Dihedral");
    failed += compare(energies.vdw_14, 0., "1-4 vdW");
    failed += compare(energies.elec_14, 0., "1-4 Elec");
    failed += compare(energies.vdw, 12567.4552, "van der Waals");
    failed += compare(energies.elec, -42333.7294, "Electrostatic");

    if (failed > 0) {
        printf("Possible FAILURE\n");
    } else {
        printf("PASSED\n");
    }
    for (i = 0; i < 62; i++) printf("="); printf("\n");

    // Now check that get_positions works as expected
    printf("Checking get_positions call\n");
    failed = 0;
    get_positions(forces);
    for (i = 0; i < sander_natom() * 3; i++) {
        failed += compare(forces[i], coordinates[i], "Coordinates");
    }
    sander_cleanup();
    free(forces);
    free(coordinates);

    if (failed > 0) {
        printf("Possible FAILURE\n");
    } else {
        printf("PASSED\n");
    }
    for (i = 0; i < 62; i++) printf("="); printf("\n");

    return 0;
}

int compare(const double computed, const double regression, const char* desc) {

    // Compare to 4 decimal places

    if (fabs(computed - regression) > 2.0e-4) {
        printf("%s failed: Expected %15.4f got %15.4f\n",
               desc, regression, computed);
        return 1;
    }
    return 0;
}

int slko_files_exist(void) {

    FILE *test;

    test = fopen("../../dat/slko/C-C.skf", "r");
    if (test == NULL)
        return 0;
    else
        fclose(test);
    return 1;
}
