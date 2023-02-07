/* Define the types that we want to make available. These include the input
 * structs and the energy struct. Make sure these stay up-to-date with the
 * structs defined in the sander API (sander.h), which in turn have to remain
 * up-to-date with the types defined in the Fortran module
 */

// Input options
typedef struct {
    PyObject_HEAD
    PyObject *igb;      // int
    PyObject *alpb;     // int
    PyObject *gbsa;     // int
    PyObject *lj1264;   // int
    PyObject *ipb;      // int
    PyObject *inp;      // int
    PyObject *vdwmeth;  // int
    PyObject *ew_type;  // int
    PyObject *ntb;      // int
    PyObject *ifqnt;    // int
    PyObject *jfastw;   // int
    PyObject *ntf;      // int
    PyObject *ntc;      // int

    PyObject *extdiel;  // double
    PyObject *intdiel;  // double
    PyObject *rgbmax;   // double
    PyObject *saltcon;  // double
    PyObject *cut;      // double
    PyObject *dielc;    // double
    PyObject *rdt;      // double
} pysander_InputOptions;

static void
pysander_InputOptions_dealloc(pysander_InputOptions* self) {
    Py_DECREF(self->igb);
    Py_DECREF(self->alpb);
    Py_DECREF(self->gbsa);
    Py_DECREF(self->lj1264);
    Py_DECREF(self->ipb);
    Py_DECREF(self->inp);
    Py_DECREF(self->vdwmeth);
    Py_DECREF(self->ew_type);
    Py_DECREF(self->ntb);
    Py_DECREF(self->ifqnt);
    Py_DECREF(self->jfastw);
    Py_DECREF(self->ntf);
    Py_DECREF(self->ntc);

    Py_DECREF(self->extdiel);
    Py_DECREF(self->intdiel);
    Py_DECREF(self->rgbmax);
    Py_DECREF(self->saltcon);
    Py_DECREF(self->cut);
    Py_DECREF(self->dielc);
    Py_DECREF(self->rdt);
    PY_DESTROY_TYPE;
}

static PyObject *
pysander_InputOptions_new(PyTypeObject *type) {
    pysander_InputOptions *self;
    self = (pysander_InputOptions *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->igb = PyInt_FromLong(0);
        self->alpb = PyInt_FromLong(0);
        self->gbsa = PyInt_FromLong(0);
        self->lj1264 = PyInt_FromLong(0);
        self->ipb = PyInt_FromLong(0);
        self->inp = PyInt_FromLong(0);
        self->vdwmeth = PyInt_FromLong(0);
        self->ew_type = PyInt_FromLong(0);
        self->ntb = PyInt_FromLong(0);
        self->ifqnt = PyInt_FromLong(0);
        self->jfastw = PyInt_FromLong(0);
        self->ntf = PyInt_FromLong(0);
        self->ntc = PyInt_FromLong(0);

        self->extdiel = PyFloat_FromDouble(0.0);
        self->intdiel = PyFloat_FromDouble(0.0);
        self->rgbmax = PyFloat_FromDouble(0.0);
        self->saltcon = PyFloat_FromDouble(0.0);
        self->cut = PyFloat_FromDouble(0.0);
        self->dielc = PyFloat_FromDouble(0.0);
        self->rdt = PyFloat_FromDouble(0.0);
    }

    return (PyObject *) self;
}

static PyMemberDef pysander_InputOptionMembers[] = {
    {"igb", T_OBJECT_EX, offsetof(pysander_InputOptions, igb), 0,
                "GB model to use"},
    {"alpb", T_OBJECT_EX, offsetof(pysander_InputOptions, alpb), 0,
                "Whether to use ALPB"},
    {"gbsa", T_OBJECT_EX, offsetof(pysander_InputOptions, gbsa), 0,
                "Whether to use a SASA term with GB"},
    {"lj1264", T_OBJECT_EX, offsetof(pysander_InputOptions, lj1264), 0,
                "Use the 12-6-4 potential"},
    {"ipb", T_OBJECT_EX, offsetof(pysander_InputOptions, ipb), 0,
                "Use PB"},
    {"inp", T_OBJECT_EX, offsetof(pysander_InputOptions, inp), 0,
                "PB SASA model to use"},
    {"vdwmeth", T_OBJECT_EX, offsetof(pysander_InputOptions, vdwmeth), 0,
                "Whether to use long-range dispersion correction"},
    {"ew_type", T_OBJECT_EX, offsetof(pysander_InputOptions, ew_type), 0,
                "Determines whether to use PME or Ewald for long-range electrostatics\n"
                "0 - Use PME\n"
                "1 - Use Ewald\n"},
    {"ntb", T_OBJECT_EX, offsetof(pysander_InputOptions, ntb), 0,
                "Whether PBC are present"},
    {"ifqnt", T_OBJECT_EX, offsetof(pysander_InputOptions, ifqnt), 0,
                "Whether to use QM/MM"},
    {"jfastw", T_OBJECT_EX, offsetof(pysander_InputOptions, jfastw), 0,
                "Whether to use analytical constraint algo. for 3-pt. waters"},
    {"ntf", T_OBJECT_EX, offsetof(pysander_InputOptions, ntf), 0,
                "Which (if any) potential energy terms are omitted"},
    {"ntc", T_OBJECT_EX, offsetof(pysander_InputOptions, ntc), 0,
                "Flag to set whether or not SHAKE is used to constrain bonds"},

    {"extdiel", T_OBJECT_EX, offsetof(pysander_InputOptions, extdiel), 0,
                "External dielectric constant for GB"},
    {"intdiel", T_OBJECT_EX, offsetof(pysander_InputOptions, intdiel), 0,
                "Internal dielectric constant for GB"},
    {"rgbmax", T_OBJECT_EX, offsetof(pysander_InputOptions, rgbmax), 0,
                "Effective radii cutoff"},
    {"saltcon", T_OBJECT_EX, offsetof(pysander_InputOptions, saltcon), 0,
                "GB salt concentration (M)"},
    {"cut", T_OBJECT_EX, offsetof(pysander_InputOptions, cut), 0,
                "Nonbonded cutoff"},
    {"dielc", T_OBJECT_EX, offsetof(pysander_InputOptions, dielc), 0,
                "dielectric constant"},
    {"rdt", T_OBJECT_EX, offsetof(pysander_InputOptions, rdt), 0,
                "Cutoff determining when only a single effective GB radius will\n"
                "be used when computing energies with LES"},
    {NULL} /* sentinel */
};

static PyTypeObject pysander_InputOptionsType = {
#if PY_MAJOR_VERSION >= 3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,                              // ob_size
#endif
    "sander.pysander.InputOptions", // tp_name
    sizeof(pysander_InputOptions),  // tp_basicsize
    0,                              // tp_itemsize
    (destructor)pysander_InputOptions_dealloc, // tp_dealloc
    0,                              // tp_print
    0,                              // tp_getattr
    0,                              // tp_setattr
    0,                              // tp_compare
    0,                              // tp_repr
    0,                              // tp_as_number
    0,                              // tp_as_sequence
    0,                              // tp_as_mapping
    0,                              // tp_hash
    0,                              // tp_call
    0,                              // tp_str
    0,                              // tp_getattro
    0,                              // tp_setattro
    0,                              // tp_as_buffer
    Py_TPFLAGS_DEFAULT,             // tp_flags
    "List of sander input options", // tp_doc
    0,		                        // tp_traverse
    0,		                        // tp_clear
    0,		                        // tp_richcompare
    0,		                        // tp_weaklistoffset
    0,		                        // tp_iter
    0,		                        // tp_iternext
    0,                              // tp_methods
    pysander_InputOptionMembers,    // tp_members
    0,                              // tp_getset
    0,                              // tp_base
    0,                              // tp_dict
    0,                              // tp_descr_get
    0,                              // tp_descr_set
    0,                              // tp_dictoffset
    0,                              // tp_init
    0,                              // tp_alloc
    (newfunc)pysander_InputOptions_new,// tp_new

};

// Energy struct
typedef struct {
    PyObject_HEAD
    PyObject *tot;
    PyObject *vdw;
    PyObject *elec;
    PyObject *gb;
    PyObject *bond;
    PyObject *angle;
    PyObject *dihedral;
    PyObject *vdw_14;
    PyObject *elec_14;
    PyObject *constraint;
    PyObject *polar;
    PyObject *hbond;
    PyObject *surf;
    PyObject *scf;
    PyObject *disp;
    PyObject *dvdl;
    PyObject *angle_ub;
    PyObject *imp;
    PyObject *cmap;
    PyObject *emap;
    PyObject *les;
    PyObject *noe;
    PyObject *pb;
    PyObject *rism;
    PyObject *ct;
    PyObject *amd_boost;
} pysander_EnergyTerms;

#define ASSIGN(var) self->var = PyFloat_FromDouble(0.0)

static PyObject *
pysander_EnergyTerms_new(PyTypeObject *type) {
    pysander_EnergyTerms *self;
    self = (pysander_EnergyTerms *)type->tp_alloc(type, 0);
    if (self != NULL) {
        ASSIGN(tot);
        ASSIGN(vdw);
        ASSIGN(elec);
        ASSIGN(gb);
        ASSIGN(bond);
        ASSIGN(angle);
        ASSIGN(dihedral);
        ASSIGN(vdw_14);
        ASSIGN(elec_14);
        ASSIGN(constraint);
        ASSIGN(polar);
        ASSIGN(hbond);
        ASSIGN(surf);
        ASSIGN(scf);
        ASSIGN(disp);
        ASSIGN(dvdl);
        ASSIGN(angle_ub);
        ASSIGN(imp);
        ASSIGN(cmap);
        ASSIGN(emap);
        ASSIGN(les);
        ASSIGN(noe);
        ASSIGN(pb);
        ASSIGN(rism);
        ASSIGN(ct);
        ASSIGN(amd_boost);
    }

    return (PyObject *) self;
}

#undef ASSIGN

static void
pysander_EnergyTerms_dealloc(pysander_EnergyTerms* self) {
    Py_DECREF(self->tot);
    Py_DECREF(self->vdw);
    Py_DECREF(self->elec);
    Py_DECREF(self->gb);
    Py_DECREF(self->bond);
    Py_DECREF(self->angle);
    Py_DECREF(self->dihedral);
    Py_DECREF(self->vdw_14);
    Py_DECREF(self->elec_14);
    Py_DECREF(self->constraint);
    Py_DECREF(self->polar);
    Py_DECREF(self->hbond);
    Py_DECREF(self->surf);
    Py_DECREF(self->scf);
    Py_DECREF(self->disp);
    Py_DECREF(self->dvdl);
    Py_DECREF(self->angle_ub);
    Py_DECREF(self->imp);
    Py_DECREF(self->cmap);
    Py_DECREF(self->emap);
    Py_DECREF(self->les);
    Py_DECREF(self->noe);
    Py_DECREF(self->pb);
    Py_DECREF(self->rism);
    Py_DECREF(self->ct);
    Py_DECREF(self->amd_boost);
    PY_DESTROY_TYPE;
}

static PyMemberDef pysander_EnergyTermsMembers[] = {
    {"tot", T_OBJECT_EX, offsetof(pysander_EnergyTerms, tot), 0,
                "Total potential energy"},
    {"vdw", T_OBJECT_EX, offsetof(pysander_EnergyTerms, vdw), 0,
                "van der Waals energy (excluding 1-4)"},
    {"elec", T_OBJECT_EX, offsetof(pysander_EnergyTerms, elec), 0,
                "Electrostatic energy (excluding 1-4)"},
    {"gb", T_OBJECT_EX, offsetof(pysander_EnergyTerms, gb), 0,
                "Generalized Born polar solvation energy"},
    {"bond", T_OBJECT_EX, offsetof(pysander_EnergyTerms, bond), 0,
                "Bond energy"},
    {"angle", T_OBJECT_EX, offsetof(pysander_EnergyTerms, angle), 0,
                "Angle energy"},
    {"dihedral", T_OBJECT_EX, offsetof(pysander_EnergyTerms, dihedral), 0,
                "Dihedral energy (including impropers)"},
    {"vdw_14", T_OBJECT_EX, offsetof(pysander_EnergyTerms, vdw_14), 0,
                "1-4 van der Waals energy"},
    {"elec_14", T_OBJECT_EX, offsetof(pysander_EnergyTerms, elec_14), 0,
                "1-4 electrostatic energy"},
    {"constraint", T_OBJECT_EX, offsetof(pysander_EnergyTerms, constraint), 0,
                "Restraint energy"},
    {"polar", T_OBJECT_EX, offsetof(pysander_EnergyTerms, polar), 0,
                "Polarization energy (for polarized force fields)"},
    {"hbond", T_OBJECT_EX, offsetof(pysander_EnergyTerms, hbond), 0,
                "Hydrogen bond (10-12 potential) energy"},
    {"surf", T_OBJECT_EX, offsetof(pysander_EnergyTerms, surf), 0,
                "Nonpolar solvation energy for implicit solvent"},
    {"scf", T_OBJECT_EX, offsetof(pysander_EnergyTerms, scf), 0,
                "QM energy"},
    {"disp", T_OBJECT_EX, offsetof(pysander_EnergyTerms, disp), 0,
                "Dispersion nonpolar solvation energy from PB"},
    {"dvdl", T_OBJECT_EX, offsetof(pysander_EnergyTerms, dvdl), 0,
                "DV/DL from TI"},
    {"angle_ub", T_OBJECT_EX, offsetof(pysander_EnergyTerms, angle_ub), 0,
                "Urey-Bradley energy (CHARMM FF only)"},
    {"imp", T_OBJECT_EX, offsetof(pysander_EnergyTerms, imp), 0,
                "Improper torsion energy (CHARMM FF only)"},
    {"cmap", T_OBJECT_EX, offsetof(pysander_EnergyTerms, cmap), 0,
                "Coupled torsion correction map energy (CHARMM only)"},
    {"emap", T_OBJECT_EX, offsetof(pysander_EnergyTerms, emap), 0,
                "Energy map restraint energy"},
    {"les", T_OBJECT_EX, offsetof(pysander_EnergyTerms, les), 0,
                "LES energy"},
    {"noe", T_OBJECT_EX, offsetof(pysander_EnergyTerms, noe), 0,
                "NOE restraint energy"},
    {"pb", T_OBJECT_EX, offsetof(pysander_EnergyTerms, pb), 0,
                "PB polar solvation energy"},
    {"rism", T_OBJECT_EX, offsetof(pysander_EnergyTerms, rism), 0,
                "3D-RISM energy"},
    {"ct", T_OBJECT_EX, offsetof(pysander_EnergyTerms, ct), 0,
                "Charge-transfer energy (from charge-relocation module)"},
    {"amd_boost", T_OBJECT_EX, offsetof(pysander_EnergyTerms, amd_boost), 0,
                "accelerated MD boost energy"},
    {NULL} /* sentinel */
};

static PyTypeObject pysander_EnergyTermsType = {
#if PY_MAJOR_VERSION >= 3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,                              // ob_size
#endif
    "sander.pysander.EnergyTerms",  // tp_name
    sizeof(pysander_EnergyTerms),   // tp_basicsize
    0,                              // tp_itemsize
    (destructor)pysander_EnergyTerms_dealloc, // tp_dealloc
    0,                              // tp_print
    0,                              // tp_getattr
    0,                              // tp_setattr
    0,                              // tp_compare
    0,                              // tp_repr
    0,                              // tp_as_number
    0,                              // tp_as_sequence
    0,                              // tp_as_mapping
    0,                              // tp_hash 
    0,                              // tp_call
    0,                              // tp_str
    0,                              // tp_getattro
    0,                              // tp_setattro
    0,                              // tp_as_buffer
    Py_TPFLAGS_DEFAULT,             // tp_flags
    "List of sander energy terms",  // tp_doc 
    0,		                        // tp_traverse
    0,		                        // tp_clear
    0,		                        // tp_richcompare
    0,		                        // tp_weaklistoffset
    0,		                        // tp_iter
    0,		                        // tp_iternext
    0,                              // tp_methods
    pysander_EnergyTermsMembers,    // tp_members
    0,                              // tp_getset
    0,                              // tp_base
    0,                              // tp_dict
    0,                              // tp_descr_get
    0,                              // tp_descr_set
    0,                              // tp_dictoffset
    0,                              // tp_init
    0,                              // tp_alloc
    (newfunc)pysander_EnergyTerms_new,// tp_new
    
};

// QM/MM options
typedef struct {
    PyObject_HEAD
    PyObject *iqmatoms; // List, length MAX_QUANTUM_ATOMS
    PyObject *qmgb;
    PyObject *lnk_atomic_no;
    PyObject *ndiis_matrices;
    PyObject *ndiis_attempts;
    PyObject *lnk_method;
    PyObject *qmcharge;
    PyObject *corecharge;
    PyObject *buffercharge;
    PyObject *spin;
    PyObject *qmqmdx;
    PyObject *verbosity;
    PyObject *printcharges;
    PyObject *printdipole;
    PyObject *print_eigenvalues;
    PyObject *peptide_corr;
    PyObject *itrmax;
    PyObject *printbondorders;
    PyObject *qmshake;
    PyObject *qmmmrij_incore;
    PyObject *qmqm_erep_incore;
    PyObject *pseudo_diag;
    PyObject *qm_ewald;
    PyObject *qm_pme;
    PyObject *kmaxqx;
    PyObject *kmaxqy;
    PyObject *kmaxqz;
    PyObject *ksqmaxq;
    PyObject *qmmm_int;
    PyObject *adjust_q;
    PyObject *tight_p_conv;
    PyObject *diag_routine;
    PyObject *density_predict;
    PyObject *fock_predict;
    PyObject *vsolv;
    PyObject *dftb_maxiter;
    PyObject *dftb_disper;
    PyObject *dftb_chg;
    PyObject *abfqmmm;
    PyObject *hot_spot;
    PyObject *qmmm_switch;
    PyObject *core_iqmatoms; // List, length MAX_QUANTUM_ATOMS
    PyObject *buffer_iqmatoms; // List, length MAX_QUANTUM_ATOMS
    PyObject *qmcut;
    PyObject *lnk_dis;
    PyObject *scfconv;
    PyObject *errconv;
    PyObject *dftb_telec;
    PyObject *dftb_telec_step;
    PyObject *fockp_d1;
    PyObject *fockp_d2;
    PyObject *fockp_d3;
    PyObject *fockp_d4;
    PyObject *damp;
    PyObject *vshift;
    PyObject *kappa;
    PyObject *pseudo_diag_criteria;
    PyObject *min_heavy_mass;
    PyObject *r_switch_hi;
    PyObject *r_switch_lo;
    PyObject *qmmask;           // String, length 8192
    PyObject *coremask;         // String, length 8192
    PyObject *buffermask;       // String, length 8192
    PyObject *centermask;       // String, length 8192
    PyObject *dftb_3rd_order;   // String, length 256
    PyObject *qm_theory;        // String, length 12
} pysander_QmInputOptions;

#define ASSIGN_FLOAT(var) self->var = PyFloat_FromDouble(inp.var)
#define ASSIGN_LIST(var, len) self->var = PyList_New(len)
#if PY_MAJOR_VERSION >= 3
#   define ASSIGN_INT(var) self->var = PyInt_FromLong(inp.var)
#   define ASSIGN_STRING(var, val) self->var = PyUnicode_FromString(val)
#else
#   define ASSIGN_INT(var) self->var = PyLong_FromLong(inp.var)
#   define ASSIGN_STRING(var, val) self->var = PyString_FromString(val)
#endif

static PyObject *
pysander_QmInputOptions_new(PyTypeObject *type) {
    Py_ssize_t i;
    qmmm_input_options inp;
    qm_sander_input(&inp);
    pysander_QmInputOptions *self;
    self = (pysander_QmInputOptions *)type->tp_alloc(type, 0);
    if (self != NULL) {
        ASSIGN_LIST(iqmatoms, MAX_QUANTUM_ATOMS);
        ASSIGN_INT(qmgb);
        ASSIGN_INT(lnk_atomic_no);
        ASSIGN_INT(ndiis_matrices);
        ASSIGN_INT(ndiis_attempts);
        ASSIGN_INT(lnk_method);
        ASSIGN_INT(qmcharge);
        ASSIGN_INT(corecharge);
        ASSIGN_INT(buffercharge);
        ASSIGN_INT(spin);
        ASSIGN_INT(qmqmdx);
        ASSIGN_INT(verbosity);
        ASSIGN_INT(printcharges);
        ASSIGN_INT(printdipole);
        ASSIGN_INT(print_eigenvalues);
        ASSIGN_INT(peptide_corr);
        ASSIGN_INT(itrmax);
        ASSIGN_INT(printbondorders);
        ASSIGN_INT(qmshake);
        ASSIGN_INT(qmmmrij_incore);
        ASSIGN_INT(qmqm_erep_incore);
        ASSIGN_INT(pseudo_diag);
        ASSIGN_INT(qm_ewald);
        ASSIGN_INT(qm_pme);
        ASSIGN_INT(kmaxqx);
        ASSIGN_INT(kmaxqy);
        ASSIGN_INT(kmaxqz);
        ASSIGN_INT(ksqmaxq);
        ASSIGN_INT(qmmm_int);
        ASSIGN_INT(adjust_q);
        ASSIGN_INT(tight_p_conv);
        ASSIGN_INT(diag_routine);
        ASSIGN_INT(density_predict);
        ASSIGN_INT(fock_predict);
        ASSIGN_INT(vsolv);
        ASSIGN_INT(dftb_maxiter);
        ASSIGN_INT(dftb_disper);
        ASSIGN_INT(dftb_chg);
        ASSIGN_INT(abfqmmm);
        ASSIGN_INT(hot_spot);
        ASSIGN_INT(qmmm_switch);
        ASSIGN_LIST(core_iqmatoms, MAX_QUANTUM_ATOMS);
        ASSIGN_LIST(buffer_iqmatoms, MAX_QUANTUM_ATOMS);

        ASSIGN_FLOAT(qmcut);
        ASSIGN_FLOAT(lnk_dis);
        ASSIGN_FLOAT(scfconv);
        ASSIGN_FLOAT(errconv);
        ASSIGN_FLOAT(dftb_telec);
        ASSIGN_FLOAT(dftb_telec_step);
        ASSIGN_FLOAT(fockp_d1);
        ASSIGN_FLOAT(fockp_d2);
        ASSIGN_FLOAT(fockp_d3);
        ASSIGN_FLOAT(fockp_d4);
        ASSIGN_FLOAT(damp);
        ASSIGN_FLOAT(vshift);
        ASSIGN_FLOAT(kappa);
        ASSIGN_FLOAT(pseudo_diag_criteria);
        ASSIGN_FLOAT(min_heavy_mass);
        ASSIGN_FLOAT(r_switch_hi);
        ASSIGN_FLOAT(r_switch_lo);

        ASSIGN_STRING(qmmask, "");
        ASSIGN_STRING(coremask, "");
        ASSIGN_STRING(buffermask, "");
        ASSIGN_STRING(centermask, "");
        ASSIGN_STRING(dftb_3rd_order, "NONE");
        ASSIGN_STRING(qm_theory, "");

        // Now assign all of the lists to zeros
        for (i = 0; i < MAX_QUANTUM_ATOMS; i++) {
            PyList_SetItem(self->iqmatoms, i, PyInt_FromLong(0));
            PyList_SetItem(self->core_iqmatoms, i, PyInt_FromLong(0));
            PyList_SetItem(self->buffer_iqmatoms, i, PyInt_FromLong(0));
        }
    }

    return (PyObject *) self;
}

#undef ASSIGN_INT
#undef ASSIGN_FLOAT
#undef ASSIGN_LIST
#undef ASSIGN_STRING

static void pysander_QmInputOptions_dealloc(pysander_QmInputOptions *self) {
    Py_DECREF(self->iqmatoms);
    Py_DECREF(self->qmgb);
    Py_DECREF(self->lnk_atomic_no);
    Py_DECREF(self->ndiis_matrices);
    Py_DECREF(self->ndiis_attempts);
    Py_DECREF(self->lnk_method);
    Py_DECREF(self->qmcharge);
    Py_DECREF(self->corecharge);
    Py_DECREF(self->buffercharge);
    Py_DECREF(self->spin);
    Py_DECREF(self->qmqmdx);
    Py_DECREF(self->verbosity);
    Py_DECREF(self->printcharges);
    Py_DECREF(self->printdipole);
    Py_DECREF(self->print_eigenvalues);
    Py_DECREF(self->peptide_corr);
    Py_DECREF(self->itrmax);
    Py_DECREF(self->printbondorders);
    Py_DECREF(self->qmshake);
    Py_DECREF(self->qmmmrij_incore);
    Py_DECREF(self->qmqm_erep_incore);
    Py_DECREF(self->pseudo_diag);
    Py_DECREF(self->qm_ewald);
    Py_DECREF(self->qm_pme);
    Py_DECREF(self->kmaxqx);
    Py_DECREF(self->kmaxqy);
    Py_DECREF(self->kmaxqz);
    Py_DECREF(self->ksqmaxq);
    Py_DECREF(self->qmmm_int);
    Py_DECREF(self->adjust_q);
    Py_DECREF(self->tight_p_conv);
    Py_DECREF(self->diag_routine);
    Py_DECREF(self->density_predict);
    Py_DECREF(self->fock_predict);
    Py_DECREF(self->vsolv);
    Py_DECREF(self->dftb_maxiter);
    Py_DECREF(self->dftb_disper);
    Py_DECREF(self->dftb_chg);
    Py_DECREF(self->abfqmmm);
    Py_DECREF(self->hot_spot);
    Py_DECREF(self->qmmm_switch);
    Py_DECREF(self->core_iqmatoms);
    Py_DECREF(self->buffer_iqmatoms);
    Py_DECREF(self->qmcut);
    Py_DECREF(self->lnk_dis);
    Py_DECREF(self->scfconv);
    Py_DECREF(self->errconv);
    Py_DECREF(self->dftb_telec);
    Py_DECREF(self->dftb_telec_step);
    Py_DECREF(self->fockp_d1);
    Py_DECREF(self->fockp_d2);
    Py_DECREF(self->fockp_d3);
    Py_DECREF(self->fockp_d4);
    Py_DECREF(self->damp);
    Py_DECREF(self->vshift);
    Py_DECREF(self->kappa);
    Py_DECREF(self->pseudo_diag_criteria);
    Py_DECREF(self->min_heavy_mass);
    Py_DECREF(self->r_switch_hi);
    Py_DECREF(self->r_switch_lo);
    Py_DECREF(self->qmmask);
    Py_DECREF(self->coremask);
    Py_DECREF(self->buffermask);
    Py_DECREF(self->centermask);
    Py_DECREF(self->dftb_3rd_order);
    Py_DECREF(self->qm_theory);
    PY_DESTROY_TYPE;
}
static PyMemberDef pysander_QmInputOptionsMembers[] = {
    {"iqmatoms", T_OBJECT_EX, offsetof(pysander_QmInputOptions, iqmatoms), 0,
        "List of atom indexes (starting from 1) to be treated using QM"},
    {"qmgb", T_OBJECT_EX, offsetof(pysander_QmInputOptions, qmgb), 0,
        "GB model to use for QM region (use InputOptions.igb instead)"},
    {"lnk_atomic_no", T_OBJECT_EX, offsetof(pysander_QmInputOptions, lnk_atomic_no), 0,
        "Atomic number of element to use as link atoms"},
    {"ndiis_matrices", T_OBJECT_EX, offsetof(pysander_QmInputOptions, ndiis_matrices), 0,
        "Number of previous error matrices to use in DIIS convergence"},
    {"ndiis_attempts", T_OBJECT_EX, offsetof(pysander_QmInputOptions, ndiis_attempts), 0,
        "Number of DIIS attempts"},
    {"lnk_method", T_OBJECT_EX, offsetof(pysander_QmInputOptions, lnk_method), 0,
        "Link atom method"},
    {"qmcharge", T_OBJECT_EX, offsetof(pysander_QmInputOptions, qmcharge), 0,
        "Charge of QM region (integer)"},
    {"corecharge", T_OBJECT_EX, offsetof(pysander_QmInputOptions, corecharge), 0,
        "Charge of QM core region (integer)"},
    {"buffercharge", T_OBJECT_EX, offsetof(pysander_QmInputOptions, buffercharge), 0,
        "Charge of QM buffer region (integer)"},
    {"spin", T_OBJECT_EX, offsetof(pysander_QmInputOptions, spin), 0,
        "Spin multiplicity"},
    {"qmqmdx", T_OBJECT_EX, offsetof(pysander_QmInputOptions, qmqmdx), 0,
        "Analytical (1) or numerical (2) QM-QM derivatives"},
    {"verbosity", T_OBJECT_EX, offsetof(pysander_QmInputOptions, verbosity), 0,
        "QM/MM verbosity (should always be 0)"},
    {"printcharges", T_OBJECT_EX, offsetof(pysander_QmInputOptions, printcharges), 0,
        "Whether to print QM charges to stdout"},
    {"printdipole", T_OBJECT_EX, offsetof(pysander_QmInputOptions, printdipole), 0,
        "Whether to print QM dipoles to stdout"},
    {"print_eigenvalues", T_OBJECT_EX, offsetof(pysander_QmInputOptions, print_eigenvalues), 0,
        "Whether to print QM eigenvalues"},
    {"peptide_corr", T_OBJECT_EX, offsetof(pysander_QmInputOptions, peptide_corr), 0,
        "Don't (0) or Do (1) apply a correction to peptide linkages"},
    {"itrmax", T_OBJECT_EX, offsetof(pysander_QmInputOptions, itrmax), 0,
        "Maximum number of SCF iterations"},
    {"printbondorders", T_OBJECT_EX, offsetof(pysander_QmInputOptions, printbondorders), 0,
        "Whether to print bond orders to stdout"},
    {"qmshake", T_OBJECT_EX, offsetof(pysander_QmInputOptions, qmshake), 0,
        "Whether to constrain H-heavy bonds in QM region using SHAKE"},
    {"qmmmrij_incore", T_OBJECT_EX, offsetof(pysander_QmInputOptions, qmmmrij_incore), 0,
        "????"},
    {"qmqm_erep_incore", T_OBJECT_EX, offsetof(pysander_QmInputOptions, qmqm_erep_incore), 0,
        "????"},
    {"pseudo_diag", T_OBJECT_EX, offsetof(pysander_QmInputOptions, pseudo_diag), 0,
        "Whether to use pseudo-diagonalizer for Fock matrix"},
    {"qm_ewald", T_OBJECT_EX, offsetof(pysander_QmInputOptions, qm_ewald), 0,
        "Whether to use Ewald in QM/MM (overridden by other settings)"},
    {"qm_pme", T_OBJECT_EX, offsetof(pysander_QmInputOptions, qm_pme), 0,
        "Whether to use Particle mesh Ewald in QM/MM (overridden by other settings)"},
    {"kmaxqx", T_OBJECT_EX, offsetof(pysander_QmInputOptions, kmaxqx), 0,
        "Max number of k-space vectors to use in X dimension for QM PME/Ewald"},
    {"kmaxqy", T_OBJECT_EX, offsetof(pysander_QmInputOptions, kmaxqy), 0,
        "Max number of k-space vectors to use in Y dimension for QM PME/Ewald"},
    {"kmaxqz", T_OBJECT_EX, offsetof(pysander_QmInputOptions, kmaxqz), 0,
        "Max number of k-space vectors to use in Z dimension for QM PME/Ewald"},
    {"ksqmaxq", T_OBJECT_EX, offsetof(pysander_QmInputOptions, ksqmaxq), 0,
        "Max number of k^2 values for spherical cutoff in recip. space (QM PME/Ewald)"},
    {"qmmm_int", T_OBJECT_EX, offsetof(pysander_QmInputOptions, qmmm_int), 0,
        "Controls QM-MM interactions (see AmberTools manual)"},
    {"adjust_q", T_OBJECT_EX, offsetof(pysander_QmInputOptions, adjust_q), 0,
        "Method to control how charges are adjusted to conserve charge"},
    {"tight_p_conv", T_OBJECT_EX, offsetof(pysander_QmInputOptions, tight_p_conv), 0,
        "Controls tightness of convergence criteria on density matrix in SCF"},
    {"diag_routine", T_OBJECT_EX, offsetof(pysander_QmInputOptions, diag_routine), 0,
        "Controls which diagonalization routine is used for the Fock matrix"},
    {"density_predict", T_OBJECT_EX, offsetof(pysander_QmInputOptions, density_predict), 0,
        "Initial guess method (??? -- just use default)"},
    {"fock_predict", T_OBJECT_EX, offsetof(pysander_QmInputOptions, fock_predict), 0,
        "Initial guess for the Fock matrix"},
    {"vsolv", T_OBJECT_EX, offsetof(pysander_QmInputOptions, vsolv), 0,
        "Controls inclusion of solvent in QM region in adaptive QM/MM"},
    {"dftb_maxiter", T_OBJECT_EX, offsetof(pysander_QmInputOptions, dftb_maxiter), 0,
        "Max number of SCF iterations for DFTB calculations"},
    {"dftb_disper", T_OBJECT_EX, offsetof(pysander_QmInputOptions, dftb_disper), 0,
        "Whether to use a dispersion correction for SCC-DFTB"},
    {"dftb_chg", T_OBJECT_EX, offsetof(pysander_QmInputOptions, dftb_chg), 0,
        "Type of charges to report (0 -- Mulliken or 2 -- CM3 charges)"},
    {"abfqmmm", T_OBJECT_EX, offsetof(pysander_QmInputOptions, abfqmmm), 0,
        "Whether to do adaptive biased force QM/MM"},
    {"hot_spot", T_OBJECT_EX, offsetof(pysander_QmInputOptions, hot_spot), 0,
        "Whether to use hot spot-like adaptive QM/MM"},
    {"qmmm_switch", T_OBJECT_EX, offsetof(pysander_QmInputOptions, qmmm_switch), 0,
        "Whether to switch QM/MM non-bonded interactions"},
    {"core_iqmatoms", T_OBJECT_EX, offsetof(pysander_QmInputOptions, core_iqmatoms), 0,
        "List of QM atom indexes (starting from 1) in the core QM region"},
    {"buffer_iqmatoms", T_OBJECT_EX, offsetof(pysander_QmInputOptions, buffer_iqmatoms), 0,
        "List of QM atom indexes (starting from 1) in the buffer QM region"},
    {"qmcut", T_OBJECT_EX, offsetof(pysander_QmInputOptions, qmcut), 0,
        "Cutoff for QM-MM nonbonded interactions"},
    {"lnk_dis", T_OBJECT_EX, offsetof(pysander_QmInputOptions, lnk_dis), 0,
        "Bond distance for bonds with link atoms"},
    {"scfconv", T_OBJECT_EX, offsetof(pysander_QmInputOptions, scfconv), 0,
        "SCF convergence criteria"},
    {"errconv", T_OBJECT_EX, offsetof(pysander_QmInputOptions, errconv), 0,
        "SCF tolerance on the error matrix"},
    {"dftb_telec", T_OBJECT_EX, offsetof(pysander_QmInputOptions, dftb_telec), 0,
        "Electronic temperature (K) used to accelerate SCC-DFTB convergence"},
    {"dftb_telec_step", T_OBJECT_EX, offsetof(pysander_QmInputOptions, dftb_telec_step), 0,
        "Step size for dftb_telec changes"},
    {"fockp_d1", T_OBJECT_EX, offsetof(pysander_QmInputOptions, fockp_d1), 0,
        "????"},
    {"fockp_d2", T_OBJECT_EX, offsetof(pysander_QmInputOptions, fockp_d2), 0,
        "????"},
    {"fockp_d3", T_OBJECT_EX, offsetof(pysander_QmInputOptions, fockp_d3), 0,
        "????"},
    {"fockp_d4", T_OBJECT_EX, offsetof(pysander_QmInputOptions, fockp_d4), 0,
        "????"},
    {"damp", T_OBJECT_EX, offsetof(pysander_QmInputOptions, damp), 0,
        "Damping factor"},
    {"vshift", T_OBJECT_EX, offsetof(pysander_QmInputOptions, vshift), 0,
        "Level shifting control"},
    {"kappa", T_OBJECT_EX, offsetof(pysander_QmInputOptions, kappa), 0,
        "Set automatically to control salt concentration in GB calcs (do not use)"},
    {"pseudo_diag_criteria", T_OBJECT_EX, offsetof(pysander_QmInputOptions, pseudo_diag_criteria), 0,
        "Pseudo diagonalization 'convergence' criteria"},
    {"min_heavy_mass", T_OBJECT_EX, offsetof(pysander_QmInputOptions, min_heavy_mass), 0,
        "Smallest atomic mass to consider as a \"heavy\" atom"},
    {"r_switch_hi", T_OBJECT_EX, offsetof(pysander_QmInputOptions, r_switch_hi), 0,
        "Distance at which switched interactions are 0 (should be qmcut)"},
    {"r_switch_lo", T_OBJECT_EX, offsetof(pysander_QmInputOptions, r_switch_lo), 0,
        "Distance at which switching function turns on (default qmcut-2)"},
    {"qmmask", T_OBJECT_EX, offsetof(pysander_QmInputOptions, qmmask), 0,
        "Amber-style atom mask to specify QM region"},
    {"coremask", T_OBJECT_EX, offsetof(pysander_QmInputOptions, coremask), 0,
        "Amber-style atom mask to specify core QM region"},
    {"buffermask", T_OBJECT_EX, offsetof(pysander_QmInputOptions, buffermask), 0,
        "Amber-style atom mask to specify buffer QM region"},
    {"centermask", T_OBJECT_EX, offsetof(pysander_QmInputOptions, centermask), 0,
        "Amber-style atom mask to specify center of QM region"},
    {"dftb_3rd_order", T_OBJECT_EX, offsetof(pysander_QmInputOptions, dftb_3rd_order), 0,
        "Whether to use DFTB 3rd-order correction"},
    {"qm_theory", T_OBJECT_EX, offsetof(pysander_QmInputOptions, qm_theory), 0,
        "Level of QM theory to use (see AmberTools manual for options)"},
    {NULL} /* sentinel */
};

static PyTypeObject pysander_QmInputOptionsType = {
#if PY_MAJOR_VERSION >= 3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(NULL)
    0,                              // ob_size
#endif
    "sander.pysander.QmInputOptions",// tp_name
    sizeof(pysander_QmInputOptions),  // tp_basicsize
    0,                              // tp_itemsize
    (destructor)pysander_QmInputOptions_dealloc, // tp_dealloc
    0,                              // tp_print
    0,                              // tp_getattr
    0,                              // tp_setattr
    0,                              // tp_compare
    0,                              // tp_repr
    0,                              // tp_as_number
    0,                              // tp_as_sequence
    0,                              // tp_as_mapping
    0,                              // tp_hash
    0,                              // tp_call
    0,                              // tp_str
    0,                              // tp_getattro
    0,                              // tp_setattro
    0,                              // tp_as_buffer
    Py_TPFLAGS_DEFAULT,             // tp_flags
    "List of QM/MM input options",  // tp_doc
    0,		                        // tp_traverse
    0,		                        // tp_clear
    0,		                        // tp_richcompare
    0,		                        // tp_weaklistoffset
    0,		                        // tp_iter
    0,		                        // tp_iternext
    0,                              // tp_methods
    pysander_QmInputOptionsMembers, // tp_members
    0,                              // tp_getset
    0,                              // tp_base
    0,                              // tp_dict
    0,                              // tp_descr_get
    0,                              // tp_descr_set
    0,                              // tp_dictoffset
    0,                              // tp_init
    0,                              // tp_alloc
    (newfunc)pysander_QmInputOptions_new,// tp_new

};
