/* This program tests the sander C API by running the same tests as the Fortran
 * interface.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sander.h"

int compare(const double, const double, const char*);

int main() {

    double *forces;
    double *coordinates;
    double box[6];
    int failed, failed2, i, natom;

    sander_input options;
    qmmm_input_options qmmm_options;
    pot_ene energies;

    failed = 0;
    printf("Testing GB sander interface (diffcoords w/ RDT)\n");
    gas_sander_input(&options, 7);
    options.cut = 9999.0;
    options.rgbmax = 100.0;
    options.rdt = 0.01;

    natom = get_inpcrd_natom("../LES_GB/les.diffcoords.r");
    coordinates = (double *)malloc(natom*3*sizeof(double));
    read_inpcrd_file("../LES_GB/les.diffcoords.r", coordinates, box);
    sander_setup_mm("../LES_GB/les.prm", coordinates, box, &options);

    forces = (double*) malloc(3*sander_natom()*sizeof(double));
    energy_forces(&energies, forces);

/* Compare the energies to the output from the relevant test:
 NSTEP =        0   TIME(PS) = 1000010.000  TEMP(K) =     0.00  PRESS =     0.0
 Etot   =       -13.2705  EKtot   =         0.0000  EPtot      =       -13.2705
 BOND   =        16.5749  ANGLE   =        21.5250  DIHED      =        35.5749
 1-4 NB =         6.4411  1-4 EEL =       140.5502  VDWAALS    =        -4.6590
 EELEC  =      -198.7892  EGB     =       -30.4884  RESTRAINT  =         0.0000
*/
    failed += compare(energies.bond, 16.5749, "Bond");
    failed += compare(energies.angle, 21.5250, "Angle");
    failed += compare(energies.dihedral, 35.5749, "Dihedral");
    failed += compare(energies.vdw_14, 6.4411, "1-4 vdW");
    failed += compare(energies.elec_14, 140.5502, "1-4 Elec");
    failed += compare(energies.vdw, -4.6590, "van der Waals");
    failed += compare(energies.elec, -198.7892, "Electrostatic");
    failed += compare(energies.gb, -30.4884, "EGB");
    failed += compare(energies.surf, 0.0, "SASA (GBSA)");

    sander_cleanup();
    free(coordinates);
    free(forces);

    if (failed > 0) {
        printf("Possible FAILURE\n");
    } else {
        printf("PASSED\n");
    }
    for (i = 0; i < 62; i++) printf("="); printf("\n");

    printf("Testing GB sander interface (samecoords w/out RDT)\n");
    failed2 = failed;
    failed = 0;
    options.rdt = 0.0;

    natom = get_inpcrd_natom("../LES_GB/les.samecoords.r");
    coordinates = (double *)malloc(natom*3*sizeof(double));
    read_inpcrd_file("../LES_GB/les.samecoords.r", coordinates, box);
    sander_setup_mm("../LES_GB/les.alt.prm", coordinates, box, &options);
    forces = (double *)malloc(sander_natom()*3*sizeof(double));

    energy_forces(&energies, forces);

/* Compare the energies to the output from the relevant test:
 NSTEP =        0   TIME(PS) = 1000010.000  TEMP(K) =     0.00  PRESS =     0.0
 Etot   =       -27.4464  EKtot   =         0.0000  EPtot      =       -27.4464
 BOND   =         5.8375  ANGLE   =        19.0846  DIHED      =        32.7197
 1-4 NB =         7.1039  1-4 EEL =       141.3377  VDWAALS    =        -3.0346
 EELEC  =      -202.2822  EGB     =       -28.2130  RESTRAINT  =         0.0000
*/

    failed += compare(energies.bond, 5.8375, "Bond");
    failed += compare(energies.angle, 19.0846, "Angle");
    failed += compare(energies.dihedral, 32.7197, "Dihedral");
    failed += compare(energies.vdw_14, 7.1039, "1-4 vdW");
    failed += compare(energies.elec_14, 141.3377, "1-4 Elec");
    failed += compare(energies.vdw, -3.0346, "van der Waals");
    failed += compare(energies.elec, -202.2822, "Electrostatic");
    failed += compare(energies.gb, -28.2130, "EGB");

    sander_cleanup();
    free(coordinates);
    free(forces);

    if (failed > 0) {
        printf("Possible FAILURE\n");
    } else {
        printf("PASSED\n");
    }
    for (i = 0; i < 62; i++) printf("="); printf("\n");

    printf("Testing PME sander interface\n");
    failed2 += failed;
    failed = 0;
    pme_sander_input(&options);
    options.cut = 8.0;

    natom = get_inpcrd_natom("../LES/md.LES.x");
    coordinates = (double *)malloc(3*natom*sizeof(double));
    read_inpcrd_file("../LES/md.LES.x", coordinates, box);
    sander_setup_mm("../LES/LES.prmtop.save", coordinates, box, &options);
    forces = (double*) malloc(3*sander_natom()*sizeof(double));
    energy_forces(&energies, forces);

/* Compare the energies to the output from the relevant test:
 NSTEP =        1   TIME(PS) =      10.002  TEMP(K) =   289.54  PRESS =     0.0
 Etot   =     -3114.6572  EKtot   =       963.1779  EPtot      =     -4077.8351
 BOND   =        14.7095  ANGLE   =        34.6208  DIHED      =        35.3483
 1-4 NB =        13.0097  1-4 EEL =       274.1453  VDWAALS    =       545.6397
 EELEC  =     -4995.3084  EHBOND  =         0.0000  RESTRAINT  =         0.0000
*/
    failed += compare(energies.bond, 14.7095, "Bond");
    failed += compare(energies.angle, 34.6208, "Angle");
    failed += compare(energies.dihedral, 35.3483, "Dihedral");
    failed += compare(energies.vdw_14, 13.0097, "1-4 vdW");
    failed += compare(energies.elec_14, 274.1453, "1-4 Elec");
    failed += compare(energies.vdw, 545.6397, "van der Waals");
    failed += compare(energies.elec, -4995.3084, "Electrostatic");

    if (failed > 0) {
        printf("Possible FAILURE\n");
    } else {
        printf("PASSED\n");
    }
    for (i = 0; i < 62; i++) printf("="); printf("\n");
    free(forces);
    free(coordinates);
    sander_cleanup();

    // Exit with the appropriate exit code
    if (failed > 0 || failed2 > 0) return 1;
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
