/* Python interface to sander functionality */

// Python includes
#include <Python.h>
#include "structmember.h"

// Support versions of Python older than 2.5 that didn't define Py_ssize_t
#if PY_VERSION_HEX < 0x02050000 && !defined(PY_SSIZE_T_MIN)
typedef int Py_ssize_t;
#   define PY_SSIZE_T_MAX INT_MAX
#   define PY_SSIZE_T_MIN INT_MIN
#endif

// A set of macros for use with Py3
#include "CompatibilityMacros.h"

// Standard C includes
#include <stdio.h>
#include <string.h>

// Amber-specific includes
#include "sander.h"

// Cordion off the type definitions, since they are large
#include "pysandermoduletypes.c"

/* sander can only be set up once, and must be cleaned up before being set up
 * again. Use this module-level variable to track whether an active system is
 * set up (0 means no system is currently set up, 1 means a system is currently
 * set up). You can only get energies and forces for a system that is set up,
 * you can only clean up a system that is set up, and you can only run setup
 * when a system is not already set up.
 */
static int IS_SETUP = 0;

/* Sander setup routine -- sets up a calculation to run with the given prmtop
 * file, inpcrd file, and input options. */
static PyObject*
pysander_setup(PyObject *self, PyObject *args) {

    char *prmtop;
    double *coordinates;
    double box[6];
    PyObject *arg2, *arg3, *arg4, *arg5;
    arg2 = NULL; arg3 = NULL; arg4 = NULL; arg5 = NULL;

    sander_input input;
    qmmm_input_options qm_input;

    // Needed to blank-out the strings
    qm_sander_input(&qm_input);

    // The passed arguments
    if (!PyArg_ParseTuple(args, "sOOO|O", &prmtop, &arg2, &arg3, &arg4, &arg5))
        return NULL;

    if (IS_SETUP) {
        // Raise a RuntimeError
        PyErr_SetString(PyExc_RuntimeError,
                        "A sander system is already set up!");
        return NULL;
    }

    pysander_InputOptions *mm_inp;
    pysander_QmInputOptions *qm_inp;

    if (!PyList_Check(arg2)) {
        PyErr_SetString(PyExc_TypeError, "2nd argument must be a list");
        return NULL;
    }

    if (!PyList_Check(arg3)) {
        PyErr_SetString(PyExc_TypeError, "3rd argument must be a list");
        return NULL;
    } else if (PyList_Size(arg3) != 6) {
        PyErr_SetString(PyExc_ValueError, "3rd argument must have 6 elements");
        return NULL;
    }

    if (!PyObject_TypeCheck(arg4, &pysander_InputOptionsType)) {
        PyErr_SetString(PyExc_TypeError,
                        "4th argument must be of type InputOptions");
        return NULL;
    }

    if (arg5 && !PyObject_TypeCheck(arg5, &pysander_QmInputOptionsType)) {
        PyErr_SetString(PyExc_TypeError,
                        "5th argument must be of type QmInputOptions");
        return NULL;
    }

    mm_inp = (pysander_InputOptions *) arg4;

    // Copy over values from mm_inp to input
    input.igb = (int) PyInt_AsLong(mm_inp->igb);
    input.alpb = (int) PyInt_AsLong(mm_inp->alpb);
    input.gbsa = (int) PyInt_AsLong(mm_inp->gbsa);
    input.lj1264 = (int) PyInt_AsLong(mm_inp->lj1264);
    input.ipb = (int) PyInt_AsLong(mm_inp->ipb);
    input.inp = (int) PyInt_AsLong(mm_inp->inp);
    input.vdwmeth = (int) PyInt_AsLong(mm_inp->vdwmeth);
    input.ew_type = (int) PyInt_AsLong(mm_inp->ew_type);
    input.ntb = (int) PyInt_AsLong(mm_inp->ntb);
    input.ifqnt = (int) PyInt_AsLong(mm_inp->ifqnt);
    input.jfastw = (int) PyInt_AsLong(mm_inp->jfastw);
    input.ntf = (int) PyInt_AsLong(mm_inp->ntf);
    input.ntc = (int) PyInt_AsLong(mm_inp->ntc);

    input.extdiel = PyFloat_AsDouble(mm_inp->extdiel);
    input.intdiel = PyFloat_AsDouble(mm_inp->intdiel);
    input.rgbmax = PyFloat_AsDouble(mm_inp->rgbmax);
    input.saltcon = PyFloat_AsDouble(mm_inp->saltcon);
    input.cut = PyFloat_AsDouble(mm_inp->cut);
    input.dielc = PyFloat_AsDouble(mm_inp->dielc);
    input.rdt = PyFloat_AsDouble(mm_inp->rdt);

    if (arg5) {
        qm_inp = (pysander_QmInputOptions *) arg5;
        // Copy over values from qm_inp to qm_input
        qm_input.qmgb = (int) PyInt_AsLong(qm_inp->qmgb);
        qm_input.lnk_atomic_no = (int) PyInt_AsLong(qm_inp->lnk_atomic_no);
        qm_input.ndiis_matrices = (int) PyInt_AsLong(qm_inp->ndiis_matrices);
        qm_input.ndiis_attempts = (int) PyInt_AsLong(qm_inp->ndiis_attempts);
        qm_input.lnk_method = (int) PyInt_AsLong(qm_inp->lnk_method);
        qm_input.qmcharge = (int) PyInt_AsLong(qm_inp->qmcharge);
        qm_input.corecharge = (int) PyInt_AsLong(qm_inp->corecharge);
        qm_input.buffercharge = (int) PyInt_AsLong(qm_inp->buffercharge);
        qm_input.spin = (int) PyInt_AsLong(qm_inp->spin);
        qm_input.qmqmdx = (int) PyInt_AsLong(qm_inp->qmqmdx);
        qm_input.verbosity = (int) PyInt_AsLong(qm_inp->verbosity);
        qm_input.printcharges = (int) PyInt_AsLong(qm_inp->printcharges);
        qm_input.printdipole = (int) PyInt_AsLong(qm_inp->printdipole);
        qm_input.print_eigenvalues = (int) PyInt_AsLong(qm_inp->print_eigenvalues);
        qm_input.peptide_corr = (int) PyInt_AsLong(qm_inp->peptide_corr);
        qm_input.itrmax = (int) PyInt_AsLong(qm_inp->itrmax);
        qm_input.printbondorders = (int) PyInt_AsLong(qm_inp->printbondorders);
        qm_input.qmshake = (int) PyInt_AsLong(qm_inp->qmshake);
        qm_input.qmmmrij_incore = (int) PyInt_AsLong(qm_inp->qmmmrij_incore);
        qm_input.qmqm_erep_incore = (int) PyInt_AsLong(qm_inp->qmqm_erep_incore);
        qm_input.pseudo_diag = (int) PyInt_AsLong(qm_inp->pseudo_diag);
        qm_input.qm_ewald = (int) PyInt_AsLong(qm_inp->qm_ewald);
        qm_input.qm_pme = (int) PyInt_AsLong(qm_inp->qm_pme);
        qm_input.kmaxqx = (int) PyInt_AsLong(qm_inp->kmaxqx);
        qm_input.kmaxqy = (int) PyInt_AsLong(qm_inp->kmaxqy);
        qm_input.kmaxqz = (int) PyInt_AsLong(qm_inp->kmaxqz);
        qm_input.ksqmaxq = (int) PyInt_AsLong(qm_inp->ksqmaxq);
        qm_input.qmmm_int = (int) PyInt_AsLong(qm_inp->qmmm_int);
        qm_input.adjust_q = (int) PyInt_AsLong(qm_inp->adjust_q);
        qm_input.tight_p_conv = (int) PyInt_AsLong(qm_inp->tight_p_conv);
        qm_input.diag_routine = (int) PyInt_AsLong(qm_inp->diag_routine);
        qm_input.density_predict = (int) PyInt_AsLong(qm_inp->density_predict);
        qm_input.fock_predict = (int) PyInt_AsLong(qm_inp->fock_predict);
        qm_input.vsolv = (int) PyInt_AsLong(qm_inp->vsolv);
        qm_input.dftb_maxiter = (int) PyInt_AsLong(qm_inp->dftb_maxiter);
        qm_input.dftb_disper = (int) PyInt_AsLong(qm_inp->dftb_disper);
        qm_input.dftb_chg = (int) PyInt_AsLong(qm_inp->dftb_chg);
        qm_input.abfqmmm = (int) PyInt_AsLong(qm_inp->abfqmmm);
        qm_input.hot_spot = (int) PyInt_AsLong(qm_inp->hot_spot);
        qm_input.qmmm_switch = (int) PyInt_AsLong(qm_inp->qmmm_switch);

        qm_input.qmcut = PyFloat_AsDouble(qm_inp->qmcut);
        qm_input.lnk_dis = PyFloat_AsDouble(qm_inp->lnk_dis);
        qm_input.scfconv = PyFloat_AsDouble(qm_inp->scfconv);
        qm_input.errconv = PyFloat_AsDouble(qm_inp->errconv);
        qm_input.dftb_telec = PyFloat_AsDouble(qm_inp->dftb_telec);
        qm_input.dftb_telec_step = PyFloat_AsDouble(qm_inp->dftb_telec_step);
        qm_input.fockp_d1 = PyFloat_AsDouble(qm_inp->fockp_d1);
        qm_input.fockp_d2 = PyFloat_AsDouble(qm_inp->fockp_d2);
        qm_input.fockp_d3 = PyFloat_AsDouble(qm_inp->fockp_d3);
        qm_input.fockp_d4 = PyFloat_AsDouble(qm_inp->fockp_d4);
        qm_input.damp = PyFloat_AsDouble(qm_inp->damp);
        qm_input.vshift = PyFloat_AsDouble(qm_inp->vshift);
        qm_input.kappa = PyFloat_AsDouble(qm_inp->kappa);
        qm_input.pseudo_diag_criteria = PyFloat_AsDouble(qm_inp->pseudo_diag_criteria);
        qm_input.min_heavy_mass = PyFloat_AsDouble(qm_inp->min_heavy_mass);
        qm_input.r_switch_hi = PyFloat_AsDouble(qm_inp->r_switch_hi);
        qm_input.r_switch_lo = PyFloat_AsDouble(qm_inp->r_switch_lo);


        // Error checking on the string input options
        size_t i;
        if (!PyObject_IS_STRING(qm_inp->qmmask)) {
            PyErr_SetString(PyExc_ValueError,
                            "qmmask must be a string");
            return NULL;
        }
        if (PyString_Size(qm_inp->qmmask) >= 8192) {
            PyErr_SetString(PyExc_ValueError,
                            "qmmask must be smaller than 8192 characters");
            return NULL;
        } else {
            strncpy(qm_input.qmmask, PyString_AsString(qm_inp->qmmask),
                    PyString_Size(qm_inp->qmmask));
            for (i = PyString_Size(qm_inp->qmmask); i < 8192; i++)
                qm_input.qmmask[i] = ' ';
        }

        if (!PyObject_IS_STRING(qm_inp->coremask)) {
            PyErr_SetString(PyExc_ValueError,
                            "coremask must be a string");
            return NULL;
        } else if (PyString_Size(qm_inp->coremask) >= 8192) {
            PyErr_SetString(PyExc_ValueError,
                            "coremask must be smaller than 8192 characters");
            return NULL;
        } else {
            strncpy(qm_input.coremask, PyString_AsString(qm_inp->coremask),
                    PyString_Size(qm_inp->coremask));
            for (i = PyString_Size(qm_inp->coremask); i < 8192; i++)
                qm_input.coremask[i] = ' ';
        }

        if (!PyObject_IS_STRING(qm_inp->buffermask)) {
            PyErr_SetString(PyExc_ValueError,
                            "buffermask must be a string");
            return NULL;
        } else if (PyString_Size(qm_inp->buffermask) >= 8192) {
            PyErr_SetString(PyExc_ValueError,
                            "buffermask must be smaller than 8192 characters");
            return NULL;
        } else {
            strncpy(qm_input.buffermask, PyString_AsString(qm_inp->buffermask),
                    PyString_Size(qm_inp->buffermask));
            for (i = PyString_Size(qm_inp->buffermask); i < 8192; i++)
                qm_input.buffermask[i] = ' ';
        }

        if (!PyObject_IS_STRING(qm_inp->centermask)) {
            PyErr_SetString(PyExc_ValueError,
                            "centermask must be a string");
            return NULL;
        } else if (PyString_Size(qm_inp->centermask) >= 8192) {
            PyErr_SetString(PyExc_ValueError,
                            "centermask must be smaller than 8192 characters");
            return NULL;
        } else {
            strncpy(qm_input.centermask, PyString_AsString(qm_inp->centermask),
                    PyString_Size(qm_inp->centermask));
            for (i = PyString_Size(qm_inp->centermask); i < 8192; i++)
                qm_input.centermask[i] = ' ';
        }

        if (!PyObject_IS_STRING(qm_inp->dftb_3rd_order)) {
            PyErr_SetString(PyExc_ValueError,
                            "dftb_3rd_order must be a string");
            return NULL;
        } else if (PyString_Size(qm_inp->dftb_3rd_order) >= 256) {
            PyErr_SetString(PyExc_ValueError,
                            "dftb_3rd_order must be smaller than 256 characters");
            return NULL;
        } else {
            strncpy(qm_input.dftb_3rd_order, PyString_AsString(qm_inp->dftb_3rd_order),
                    PyString_Size(qm_inp->dftb_3rd_order));
            for (i = PyString_Size(qm_inp->dftb_3rd_order); i < 256; i++)
                qm_input.dftb_3rd_order[i] = ' ';
        }

        if (!PyObject_IS_STRING(qm_inp->qm_theory)) {
            PyErr_SetString(PyExc_ValueError,
                            "qm_theory must be a string");
            return NULL;
        } else if (PyString_Size(qm_inp->qm_theory) >= 12) {
            PyErr_SetString(PyExc_ValueError,
                            "qm_theory must be smaller than 12 characters");
            return NULL;
        } else {
            strncpy(qm_input.qm_theory, PyString_AsString(qm_inp->qm_theory),
                    (int)PyString_Size(qm_inp->qm_theory));
            for (i = PyString_Size(qm_inp->qm_theory); i < 12; i++)
                qm_input.qm_theory[i] = ' ';
        }

        // Now copy over the arrays. Check that none of them are too large
        if (!PyList_Check(qm_inp->iqmatoms)) {
            PyErr_SetString(PyExc_ValueError,
                            "iqmatoms must be a list of integers");
            return NULL;
        } else if (PyList_Size(qm_inp->iqmatoms) > MAX_QUANTUM_ATOMS) {
            PyErr_SetString(PyExc_ValueError, "iqmatoms is too large");
            return NULL;
        } else {
            Py_ssize_t i;
            for (i = 0; i < PyList_Size(qm_inp->iqmatoms); i++) {
                qm_input.iqmatoms[i] = (int) 
                        PyInt_AsLong(PyList_GetItem(qm_inp->iqmatoms, i));
            }
            for (i = PyList_Size(qm_inp->iqmatoms); i < MAX_QUANTUM_ATOMS; i++)
                qm_input.iqmatoms[i] = 0;
        }

        if (!PyList_Check(qm_inp->core_iqmatoms)) {
            PyErr_SetString(PyExc_ValueError,
                            "core_iqmatoms must be a list");
            return NULL;
        } else if (PyList_Size(qm_inp->core_iqmatoms) > MAX_QUANTUM_ATOMS) {
            PyErr_SetString(PyExc_ValueError,
                            "core_iqmatoms is too large");
            return NULL;
        } else {
            Py_ssize_t i;
            for (i = 0; i < PyList_Size(qm_inp->core_iqmatoms); i++) {
                qm_input.core_iqmatoms[i] = (int)
                        PyInt_AsLong(PyList_GetItem(qm_inp->core_iqmatoms, i));
            }
            for (i = PyList_Size(qm_inp->core_iqmatoms); i < MAX_QUANTUM_ATOMS; i++)
                qm_input.core_iqmatoms[i] = 0;
        }

        if (!PyList_Check(qm_inp->buffer_iqmatoms)) {
            PyErr_SetString(PyExc_ValueError,
                            "buffer_iqmatoms must be a list");
            return NULL;
        } else if (PyList_Size(qm_inp->buffer_iqmatoms) > MAX_QUANTUM_ATOMS) {
            PyErr_SetString(PyExc_ValueError,
                            "buffer_iqmatoms is too large");
            return NULL;
        } else {
            Py_ssize_t i;
            for (i = 0; i < PyList_Size(qm_inp->buffer_iqmatoms); i++)
                qm_input.buffer_iqmatoms[i] = (int)
                        PyInt_AsLong(PyList_GetItem(qm_inp->buffer_iqmatoms, i));
            for (i = PyList_Size(qm_inp->buffer_iqmatoms); i < MAX_QUANTUM_ATOMS; i++)
                qm_input.buffer_iqmatoms[i] = 0;
        }
    }

    Py_ssize_t i;
    coordinates = (double *)malloc(PyList_Size(arg2)*sizeof(double));
    // Fill up the positions and box
    for (i = 0; i < PyList_Size(arg2); i++)
        coordinates[i] = PyFloat_AsDouble(PyList_GetItem(arg2, i));
    for (i = 0; i < 6; i++)
        box[i] = PyFloat_AsDouble(PyList_GetItem(arg3, i));

    if (sander_setup(prmtop, coordinates, box, &input, &qm_input)) {
        free(coordinates);
        PyErr_SetString(PyExc_RuntimeError, "Problem setting up sander");
        return NULL;
    }

    free(coordinates);
    IS_SETUP = 1;

    Py_RETURN_NONE;
}

static PyObject*
pysander_set_positions(PyObject *self, PyObject *args) {
    PyObject *pypositions;
    double *positions;

    if (!PyArg_ParseTuple(args, "O", &pypositions))
        return NULL;

    if (!IS_SETUP) {
        PyErr_SetString(PyExc_RuntimeError,
                        "No sander system is currently set up!");
        return NULL;
    }

    // Check that the passed positions is legitimate

    if (!PyList_Check(pypositions)) {
        PyErr_SetString(PyExc_TypeError,
                        "set_positions expects a list of coordinates");
        return NULL;
    }

    if (PyList_Size(pypositions) != 3 * sander_natom()) {
        PyErr_SetString(PyExc_ValueError,
                        "coordinate list must have length 3*natom");
        return NULL;
    }

    // Now allocate the positions and assign them all
    positions = (double *)malloc(sander_natom()*3*sizeof(double));

    Py_ssize_t i;
    for (i = 0; i < sander_natom()*3; i++) {
        positions[(size_t)i] = PyFloat_AsDouble(PyList_GetItem(pypositions, i));
    }

    set_positions(positions);
    free(positions);
    Py_RETURN_NONE;
}

static PyObject*
pysander_set_box(PyObject *self, PyObject *args) {

    double a, b, c, alpha, beta, gamma;

    if (!PyArg_ParseTuple(args, "dddddd", &a, &b, &c, &alpha, &beta, &gamma))
        return NULL;

    if (!IS_SETUP) {
        PyErr_SetString(PyExc_RuntimeError,
                        "No sander system is currently set up!");
        return NULL;
    }

    set_box(a, b, c, alpha, beta, gamma);

    Py_RETURN_NONE;
}

/* Deallocates the memory used by sander so sander can be set up and used again
 */
static PyObject*
pysander_cleanup(PyObject *self) {
    if (!IS_SETUP) {
        // Raise a RuntimeError
        PyErr_SetString(PyExc_RuntimeError,
                        "No sander system is currently set up!");
        return NULL;
    }
    sander_cleanup();
    IS_SETUP = 0;
    Py_RETURN_NONE;
}

#define ASSIGN_INT(var) Py_DECREF(ret->var); ret->var = PyInt_FromLong((long int)inp.var)
#define ASSIGN_FLOAT(var) Py_DECREF(ret->var); ret->var = PyFloat_FromDouble(inp.var)
/* Creates an input option struct with all of the options optimized for gas
 * phase or implicit solvent (i.e., aperiodic) calculations
 */
static PyObject*
pysander_gas_input(PyObject *self, PyObject *args) {

    long tmp = 6;
    if (!PyArg_ParseTuple(args, "|i", &tmp)) {
        return NULL;
    }
    int igb = (int) tmp;
    if (igb < 0 || igb > 10 || igb == 3 || igb == 9 || igb == 4) {
        PyErr_SetString(PyExc_ValueError,
                        "igb must be 0, 1, 2, 5, 6, 7, 8, or 10");
        return NULL;
    }
    sander_input inp;
    gas_sander_input(&inp, igb);
    pysander_InputOptions *ret = (pysander_InputOptions *)
            PyObject_CallObject((PyObject *) &pysander_InputOptionsType, NULL);
    if (ret == NULL)
        return NULL;
    // Integers
    ASSIGN_INT(igb);
    ASSIGN_INT(alpb);
    ASSIGN_INT(gbsa);
    ASSIGN_INT(lj1264);
    ASSIGN_INT(ipb);
    ASSIGN_INT(inp);
    ASSIGN_INT(vdwmeth);
    ASSIGN_INT(ew_type);
    ASSIGN_INT(ntb);
    ASSIGN_INT(ifqnt);
    ASSIGN_INT(jfastw);
    ASSIGN_INT(ntf);
    ASSIGN_INT(ntc);
    // Floats
    ASSIGN_FLOAT(extdiel);
    ASSIGN_FLOAT(intdiel);
    ASSIGN_FLOAT(rgbmax);
    ASSIGN_FLOAT(saltcon);
    ASSIGN_FLOAT(cut);
    ASSIGN_FLOAT(dielc);
    ASSIGN_FLOAT(rdt);

    return (PyObject *) ret;
}

/* Creates an input option struct with all of the options optimized for PME
 * calculations
 */
static PyObject *
pysander_pme_input(PyObject *self) {
    sander_input inp;
    pme_sander_input(&inp);
    pysander_InputOptions *ret = (pysander_InputOptions *)
            PyObject_CallObject((PyObject *) &pysander_InputOptionsType, NULL);
    if (ret == NULL)
        return NULL;
    // Integers
    ASSIGN_INT(igb);
    ASSIGN_INT(alpb);
    ASSIGN_INT(gbsa);
    ASSIGN_INT(lj1264);
    ASSIGN_INT(ipb);
    ASSIGN_INT(inp);
    ASSIGN_INT(vdwmeth);
    ASSIGN_INT(ew_type);
    ASSIGN_INT(ntb);
    ASSIGN_INT(ifqnt);
    ASSIGN_INT(jfastw);
    ASSIGN_INT(ntf);
    ASSIGN_INT(ntc);
    // Floats
    ASSIGN_FLOAT(extdiel);
    ASSIGN_FLOAT(intdiel);
    ASSIGN_FLOAT(rgbmax);
    ASSIGN_FLOAT(saltcon);
    ASSIGN_FLOAT(cut);
    ASSIGN_FLOAT(dielc);
    ASSIGN_FLOAT(rdt);

    return (PyObject *) ret;
}
#undef ASSIGN_INT
#undef ASSIGN_FLOAT

/* Returns the number of atoms in the currently set-up system. If the system is
 * not set up, raise RuntimeError
 */
static PyObject *
pysander_natom(PyObject *self) {
    if (IS_SETUP == 0) {
        PyErr_SetString(PyExc_RuntimeError,
                        "Cannot query number of atoms -- no system set up");
        return NULL;
    }

    return PyInt_FromLong((long int)sander_natom());
}

static PyObject *
pysander_energy_forces(PyObject *self) {

    if (IS_SETUP == 0) {
        PyErr_SetString(PyExc_RuntimeError,
                        "Cannot compute energies and forces -- no system set up");
        return NULL;
    }

    pot_ene energies;
    int natom3 = 3 * sander_natom();
    double *forces = (double *) malloc(natom3*sizeof(double));

    energy_forces(&energies, forces);

    // Now construct the return values

    pysander_EnergyTerms *py_energies = (pysander_EnergyTerms *)
            PyObject_CallObject((PyObject *) &pysander_EnergyTermsType, NULL);
    PyObject *py_forces = PyList_New(natom3);

    py_energies->tot = PyFloat_FromDouble(energies.tot);
    py_energies->vdw = PyFloat_FromDouble(energies.vdw);
    py_energies->elec = PyFloat_FromDouble(energies.elec);
    py_energies->gb = PyFloat_FromDouble(energies.gb);
    py_energies->bond = PyFloat_FromDouble(energies.bond);
    py_energies->angle = PyFloat_FromDouble(energies.angle);
    py_energies->dihedral = PyFloat_FromDouble(energies.dihedral);
    py_energies->vdw_14 = PyFloat_FromDouble(energies.vdw_14);
    py_energies->elec_14 = PyFloat_FromDouble(energies.elec_14);
    py_energies->constraint = PyFloat_FromDouble(energies.constraint);
    py_energies->polar = PyFloat_FromDouble(energies.polar);
    py_energies->hbond = PyFloat_FromDouble(energies.hbond);
    py_energies->surf = PyFloat_FromDouble(energies.surf);
    py_energies->scf = PyFloat_FromDouble(energies.scf);
    py_energies->disp = PyFloat_FromDouble(energies.disp);
    py_energies->dvdl = PyFloat_FromDouble(energies.dvdl);
    py_energies->angle_ub = PyFloat_FromDouble(energies.angle_ub);
    py_energies->imp = PyFloat_FromDouble(energies.imp);
    py_energies->cmap = PyFloat_FromDouble(energies.cmap);
    py_energies->emap = PyFloat_FromDouble(energies.emap);
    py_energies->les = PyFloat_FromDouble(energies.les);
    py_energies->noe = PyFloat_FromDouble(energies.noe);
    py_energies->pb = PyFloat_FromDouble(energies.pb);
    py_energies->rism = PyFloat_FromDouble(energies.rism);
    py_energies->ct = PyFloat_FromDouble(energies.ct);
    py_energies->amd_boost = PyFloat_FromDouble(energies.amd_boost);

    Py_ssize_t i;
    for (i = 0; i < (Py_ssize_t) natom3; i++)
        PyList_SET_ITEM(py_forces, i, PyFloat_FromDouble(forces[i]));
    free(forces);

    PyObject *ret = PyTuple_New(2);
    PyTuple_SET_ITEM(ret, 0, (PyObject *)py_energies);
    PyTuple_SET_ITEM(ret, 1, py_forces);

    return ret;
}

static PyObject *
pysander_get_positions(PyObject *self) {

    if (IS_SETUP == 0) {
        PyErr_SetString(PyExc_RuntimeError,
                        "Cannot get positions when no system is set up.");
        return NULL;
    }

    int natom3 = 3 * sander_natom();
    double *positions = (double *) malloc(natom3*sizeof(double));

    PyObject *py_positions = PyList_New(natom3);

    get_positions(positions);

    Py_ssize_t i;
    for (i = 0; i < (Py_ssize_t) natom3; i++)
        PyList_SET_ITEM(py_positions, i, PyFloat_FromDouble(positions[i]));
    free(positions);

    return py_positions;
}

static PyObject *
pysander_is_setup(PyObject *self) {
    if (IS_SETUP == 0)
        Py_RETURN_FALSE;
    Py_RETURN_TRUE;
}

/* Python module initialization */

static PyMethodDef
pysanderMethods[] = {
    { "setup", (PyCFunction) pysander_setup, METH_VARARGS,
            "Sets up sander calc (private)"},
    { "cleanup", (PyCFunction) pysander_cleanup, METH_NOARGS,
            "Cleans up sander calc"},
    { "gas_input", (PyCFunction) pysander_gas_input, METH_VARARGS,
            "Returns a populated InputOptions instance optimized for gas-phase\n"
            "or implicit solvent calculations. Optional argument is the GB model\n"
            "you wish to use (1, 2, 5, 7, or 8)"},
    { "pme_input", (PyCFunction) pysander_pme_input, METH_NOARGS,
            "Returns a populated InputOptions instance optimized with Amber\n"
            "defaults for PME calculations.\n"},
    { "natom", (PyCFunction) pysander_natom, METH_NOARGS,
            "Returns the number of atoms in the currently set-up system"},
    { "energy_forces", (PyCFunction) pysander_energy_forces, METH_NOARGS,
            "Computes energies and forces from the given set of coordinates.\n"
            "\n"
            "Returns\n"
            "-------\n"
            "   energy : type EnergyTerms\n"
            "       An EnergyTerms instance populated with the energy components\n"
            "       in kilocalories per mole\n"
            "\n"
            "   forces : list\n"
            "       A list of all forces in kilocalories/mole/Angstroms"},
    { "set_positions", (PyCFunction) pysander_set_positions, METH_VARARGS,
            "Sets the active positions to the passed list of positions (private)"},
    { "get_positions", (PyCFunction) pysander_get_positions, METH_NOARGS,
            "Returns the currently active positions as a list"},
    { "set_box", (PyCFunction) pysander_set_box, METH_VARARGS,
            "Sets the box dimensions of the active system.\n"
            "\n"
            "Parameters\n"
            "----------\n"
            "a : float\n"
            "    Length of the first side of the unit cell\n"
            "b : float\n"
            "    Length of the second side of the unit cell\n"
            "c : float\n"
            "    Length of the third side of the unit cell\n"
            "alpha : float\n"
            "    Angle between sides b and c of the unit cell\n"
            "beta : float\n"
            "    Angle between sides a and c of the unit cell\n"
            "gamma : float\n"
            "    Angle between sides a and b of the unit cell\n"},
    { "is_setup", (PyCFunction) pysander_is_setup, METH_NOARGS,
            "Returns True if sander is set up and False otherwise"},
    {NULL}, // sentinel
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "pysander",                                                 // m_name
    "Python interface into sander energy and force evaluation", // m_doc
    -1,                                                         // m_size
    pysanderMethods,                                            // m_methods
    NULL,
    NULL,
    NULL,
    NULL,
};
#endif

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC
PyInit_pysander(void) {
    // Type declarations
    if (PyType_Ready(&pysander_InputOptionsType) < 0)
        return NULL;
    if (PyType_Ready(&pysander_EnergyTermsType) < 0)
        return NULL;
    if (PyType_Ready(&pysander_QmInputOptionsType) < 0)
        return NULL;
    PyObject* m = PyModule_Create(&moduledef);
#else
PyMODINIT_FUNC
initpysander(void) {
    if (PyType_Ready(&pysander_InputOptionsType))
        return;
    if (PyType_Ready(&pysander_EnergyTermsType))
        return;
    if (PyType_Ready(&pysander_QmInputOptionsType))
        return;
    PyObject* m = Py_InitModule3("pysander", pysanderMethods,
                "Python interface into sander energy and force evaluation");
#endif

    // Now add the types
    Py_INCREF(&pysander_InputOptionsType);
    PyModule_AddObject(m, "InputOptions", (PyObject*) &pysander_InputOptionsType);
    Py_INCREF(&pysander_EnergyTermsType);
    PyModule_AddObject(m, "EnergyTerms", (PyObject *) &pysander_EnergyTermsType);
    Py_INCREF(&pysander_QmInputOptionsType);
    PyModule_AddObject(m, "QmInputOptions", (PyObject *) &pysander_QmInputOptionsType);

#if PY_MAJOR_VERSION >= 3
    return m;
#endif
}
