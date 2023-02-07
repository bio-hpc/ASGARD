GFORTRAN module created from qmmm_module.F90 on Thu Feb 11 10:37:22 2016
MD5:70c7e668cb8f530baed5e65632bbdd9d -- If you edit this, you'll get what you deserve.

(() () () () () ()
() () () () () () () () () () () () () () () () () () () () ())

()

()

()

()

(2 'allocate_qmmm' 'qmmm_module' 'allocate_qmmm' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 3 0 (4 5 6) () 0 () () 0 0)
7 'alph_mm' 'qmmm_module' 'alph_mm' 1 ((PARAMETER UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE) (REAL 8 0 0 REAL ()) 0 0 () (
CONSTANT (REAL 8 0 0 REAL ()) 0 '0.50000000000000@1') () 0 () () 0 0)
8 'axis_tol' 'qmmm_module' 'axis_tol' 1 ((PARAMETER UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE) (REAL 8 0 0 REAL ()) 0 0 () (
CONSTANT (REAL 8 0 0 REAL ()) 0 '0.2af31dc4611874@-6') () 0 () () 0 0)
9 'deallocate_qmmm' 'qmmm_module' 'deallocate_qmmm' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 10 0 (11 12 13 14) () 0 () () 0 0)
15 'default_qmmm_input_options' 'qmmm_module' 'default_qmmm_input_options'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN SUBROUTINE) (
UNKNOWN 0 0 0 UNKNOWN ()) 16 0 (17) () 0 () () 0 0)
18 'exponential_cutoff' 'qmmm_module' 'exponential_cutoff' 1 ((
PARAMETER UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE) (REAL 8 0 0
REAL ()) 0 0 () (CONSTANT (REAL 8 0 0 REAL ()) 0 '0.1e000000000000@2') ()
0 () () 0 0)
19 'get_atomic_number' 'qmmm_module' 'get_atomic_number' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN SUBROUTINE ALWAYS_EXPLICIT) (
UNKNOWN 0 0 0 UNKNOWN ()) 20 0 (21 22 23 24) () 0 () () 0 0)
25 'overlap_cutoff' 'qmmm_module' 'overlap_cutoff' 1 ((PARAMETER
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN IMPLICIT-SAVE) (REAL 8 0 0 REAL ())
0 0 () (CONSTANT (REAL 8 0 0 REAL ()) 0 '0.1651b3f1314df4@3') () 0 () ()
0 0)
26 'qm2_params' 'qmmm_module' 'qm2_params' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN) (DERIVED 27 0 0 DERIVED ()) 0 0 () () 0 ()
() 0 0)
28 'qm2_rij_eqns' 'qmmm_module' 'qm2_rij_eqns' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN) (DERIVED 29 0 0 DERIVED ())
0 0 () () 0 () () 0 0)
29 'qm2_rij_eqns_structure' 'qmmm_module' 'qm2_rij_eqns_structure' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((30 'qmmmrijdata' (REAL 8 0 0
REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (31
'qmmmrij_allocated' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()))
PUBLIC () 0 0)
32 'qm2_struct' 'qmmm_module' 'qm2_struct' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN EXPLICIT-SAVE) (DERIVED 33 0 0 DERIVED ()) 0 0 () ()
0 () () 0 0)
33 'qm2_structure' 'qmmm_module' 'qm2_structure' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (UNKNOWN 0 0 0
UNKNOWN ()) 0 0 () () 0 ((34 'den_matrix' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(35 'old_den_matrix' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (36 'old2_density' (
REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0)) (37 'md_den_mat_guess1' (REAL 8 0 0 REAL ())
(1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ())
0)) (38 'md_den_mat_guess2' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1
0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (39
'fock_mat_final4' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (40 'fock_mat_final3'
(REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0)) (41 'fock_mat_final2' (REAL 8 0 0 REAL ())
(1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ())
0)) (42 'fock_mat_final1' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (43 'fock_matrix' (
REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0)) (44 'qm_mm_e_repul' (REAL 8 0 0 REAL ()) (
2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (45 'qm_qm_2e_repul' (REAL 8 0 0 REAL ()) (1 DEFERRED ()
()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (46
'hmatrix' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (47 'qm_qm_e_repul' (REAL 8 0 0 REAL
()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (48 'fock2_ptot2' (REAL 8 0 0 REAL ()) (2 DEFERRED () ()
() ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (49
'eigen_vectors' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (50 'eigen_values' (
REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0)) (51 'scf_mchg' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(52 'diis_fock' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (53 'diis_errmat' (
REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS (NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0)) (54 'diis_mat' (REAL 8 0 0 REAL ()) (2
DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (55 'matsize' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS
()) (56 'n2el' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
57 'norbs' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (58
'nclosed' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (59
'nopenclosed' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
60 'qm_mm_e_repul_allocated' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (61 'n_peptide_links' (INTEGER 4 0 0 INTEGER ()) () 0
0 0 UNKNOWN-ACCESS ()) (62 'peptide_links' (INTEGER 4 0 0 INTEGER ()) (
2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (63 'calc_mchg_scf' (
LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ())) PUBLIC () 0 0)
64 'qm_ewald_structure' 'qmmm_module' 'qm_ewald_structure' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (UNKNOWN 0 0 0
UNKNOWN ()) 0 0 () () 0 ((65 'kvec' (REAL 8 0 0 REAL ()) (1 DEFERRED ()
()) 1 1 0 UNKNOWN-ACCESS ()) (66 'dkvec' (REAL 8 0 0 REAL ()) (2
DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (67 'dmkv' (REAL 8 0 0
REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (68 'ktable'
(REAL 8 0 0 REAL ()) (3 DEFERRED () () () () () ()) 1 1 0 UNKNOWN-ACCESS
()) (69 'qmktable' (REAL 8 0 0 REAL ()) (3 DEFERRED () () () () () ()) 1
1 0 UNKNOWN-ACCESS ()) (70 'mmpot' (REAL 8 0 0 REAL ()) (1 DEFERRED () ())
1 1 0 UNKNOWN-ACCESS ()) (71 'qmpot' (REAL 8 0 0 REAL ()) (1 DEFERRED ()
()) 1 1 0 UNKNOWN-ACCESS ()) (72 'coulpot' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (73 'd_ewald_mm' (REAL 8 0 0
REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (74
'ewald_core' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (75
'mm_recip_e' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (76 'kappa'
(REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (77 'totkq' (INTEGER 4
0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (78 'natom' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (79 'ewald_startup' (LOGICAL 4 0
0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ())) PUBLIC () 0 0)
80 'qm_gb' 'qmmm_module' 'qm_gb' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN) (DERIVED 81 0 0 DERIVED ()) 0 0 () () 0 ()
() 0 0)
81 'qm_gb_structure' 'qmmm_module' 'qm_gb_structure' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (UNKNOWN 0 0 0
UNKNOWN ()) 0 0 () () 0 ((82 'qmqm_onefij' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (83 'qmqm_kappafij' (REAL 8 0 0
REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (84 'gb_mmpot' (
REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (85
'gb_qmpot' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS
()) (86 'intdieli' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (87
'extdieli' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (88 'kappa'
(REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (89 'mmcut2' (REAL 8 0
0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (90 'one_arad_beta' (REAL 8 0 0
REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (91 'qmqm_gb_list' (INTEGER 4 0 0
INTEGER ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (92
'saltcon_on' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (93
'alpb_on' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ())) PUBLIC
() 0 0)
94 'qmewald' 'qmmm_module' 'qmewald' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN) (DERIVED 64 0 0 DERIVED ()) 0 0 () () 0 ()
() 0 0)
95 'qmmm_div' 'qmmm_module' 'qmmm_div' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN) (DERIVED 96 0 0 DERIVED ()) 0 0 () () 0 ()
() 0 0)
96 'qmmm_div_structure' 'qmmm_module' 'qmmm_div_structure' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (UNKNOWN 0 0 0
UNKNOWN ()) 0 0 () () 0 ((97 'ntotatm' (INTEGER 4 0 0 INTEGER ()) () 0 0
0 UNKNOWN-ACCESS ()) (98 'all_atom_numbers' (INTEGER 4 0 0 INTEGER ()) (
1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ())) PUBLIC () 0 0)
99 'qmmm_input_options' 'qmmm_module' 'qmmm_input_options' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN SEQUENCE) (UNKNOWN 0 0 0
UNKNOWN ()) 0 0 () () 0 ((100 'qmcut' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (101 'lnk_dis' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (102 'scfconv' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (103 'errconv' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (104 'dftb_telec' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (105 'dftb_telec_step' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (106 'fockp_d1' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (107 'fockp_d2' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (108 'fockp_d3' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (109 'fockp_d4' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (110 'damp' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (111 'vshift' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (112 'kappa' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (113 'pseudo_diag_criteria' (REAL 8 0 0 REAL ()) () 0
0 0 UNKNOWN-ACCESS ()) (114 'min_heavy_mass' (REAL 8 0 0 REAL ()) () 0 0
0 UNKNOWN-ACCESS ()) (115 'r_switch_hi' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (116 'r_switch_lo' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (117 'iqmatoms' (INTEGER 4 0 0 INTEGER ()) (1
EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '10000')) 1 0 0 UNKNOWN-ACCESS ()) (118 'qmgb' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (119 'lnk_atomic_no'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (120
'ndiis_matrices' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(121 'ndiis_attempts' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS
()) (122 'lnk_method' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS
()) (123 'qmcharge' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(124 'corecharge' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(125 'buffercharge' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(126 'spin' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (127
'qmqmdx' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (128
'verbosity' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (129
'printcharges' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
130 'printdipole' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(131 'print_eigenvalues' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (132 'peptide_corr' (INTEGER 4 0 0 INTEGER ()) () 0 0
0 UNKNOWN-ACCESS ()) (133 'itrmax' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (134 'printbondorders' (INTEGER 4 0 0 INTEGER ()) ()
0 0 0 UNKNOWN-ACCESS ()) (135 'qmshake' (INTEGER 4 0 0 INTEGER ()) () 0
0 0 UNKNOWN-ACCESS ()) (136 'qmmmrij_incore' (INTEGER 4 0 0 INTEGER ())
() 0 0 0 UNKNOWN-ACCESS ()) (137 'qmqm_erep_incore' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (138 'pseudo_diag' (INTEGER 4 0
0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (139 'qm_ewald' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (140 'qm_pme' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (141 'kmaxqx' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (142 'kmaxqy' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (143 'kmaxqz' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (144 'ksqmaxq' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (145 'qmmm_int' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (146 'adjust_q' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (147 'tight_p_conv' (INTEGER 4 0
0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (148 'diag_routine' (INTEGER 4
0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (149 'density_predict' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (150 'fock_predict'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (151 'vsolv' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (152 'dftb_maxiter'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (153 'dftb_disper'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (154 'dftb_chg' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (155 'abfqmmm' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (156 'hot_spot' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (157 'qmmm_switch'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (158
'core_iqmatoms' (INTEGER 4 0 0 INTEGER ()) (1 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'10000')) 1 0 0 UNKNOWN-ACCESS ()) (159 'buffer_iqmatoms' (INTEGER 4 0 0
INTEGER ()) (1 EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '10000')) 1 0 0 UNKNOWN-ACCESS ())
(160 'qmmask' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '8192'))) () 0 0 0 UNKNOWN-ACCESS ()) (161 'coremask' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192')))
() 0 0 0 UNKNOWN-ACCESS ()) (162 'buffermask' (CHARACTER 1 0 0 CHARACTER
((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () 0 0 0
UNKNOWN-ACCESS ()) (163 'centermask' (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () 0 0 0 UNKNOWN-ACCESS
()) (164 'dftb_3rd_order' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '256'))) () 0 0 0 UNKNOWN-ACCESS ()) (165
'qm_theory' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER
()) 0 '12'))) () 0 0 0 UNKNOWN-ACCESS ())) PUBLIC () 0 0)
166 'qmmm_mpi' 'qmmm_module' 'qmmm_mpi' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN UNKNOWN) (DERIVED 167 0 0 DERIVED ()) 0 0 () () 0 ()
() 0 0)
167 'qmmm_mpi_structure' 'qmmm_module' 'qmmm_mpi_structure' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (UNKNOWN 0 0 0
UNKNOWN ()) 0 0 () () 0 ((168 'commqmmm' (INTEGER 4 0 0 INTEGER ()) () 0
0 0 UNKNOWN-ACCESS ()) (169 'numthreads' (INTEGER 4 0 0 INTEGER ()) () 0
0 0 UNKNOWN-ACCESS ()) (170 'mytaskid' (INTEGER 4 0 0 INTEGER ()) () 0 0
0 UNKNOWN-ACCESS ()) (171 'openmp_numthreads' (INTEGER 4 0 0 INTEGER ())
() 0 0 0 UNKNOWN-ACCESS ()) (172 'natom_start' (INTEGER 4 0 0 INTEGER ())
() 0 0 0 UNKNOWN-ACCESS ()) (173 'natom_end' (INTEGER 4 0 0 INTEGER ())
() 0 0 0 UNKNOWN-ACCESS ()) (174 'nquant_nlink_start' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (175 'nquant_nlink_end' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (176 'totkq_count'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (177 'kvec_start'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (178 'kvec_end' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (179 'two_e_offset'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (180
'nquant_nlink_istart' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS
()) (181 'nquant_nlink_iend' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (182 'nquant_nlink_loop_extent_begin' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (183
'nquant_nlink_loop_extent_end' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (184 'nquant_nlink_jrange' (INTEGER 4 0 0 INTEGER ())
(2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (185 'commqmmm_master'
(LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ())) PUBLIC () 0 0)
186 'qmmm_nml' 'qmmm_module' 'qmmm_nml' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN EXPLICIT-SAVE) (DERIVED 187 0 0 DERIVED ()) 0 0 ()
() 0 () () 0 0)
188 'qmmm_opnq' 'qmmm_module' 'qmmm_opnq' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN EXPLICIT-SAVE) (DERIVED 189 0 0 DERIVED ()) 0 0 ()
() 0 () () 0 0)
189 'qmmm_opnq_structure' 'qmmm_module' 'qmmm_opnq_structure' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((190 'useopnq' (LOGICAL 4 0 0
LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS (CONSTANT (LOGICAL 4 0 0 LOGICAL ())
0 0)) (191 'opnqcorrection' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS
()) (192 'vdwcorrection' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ())
(193 'switching' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS (
CONSTANT (LOGICAL 4 0 0 LOGICAL ()) 0 1)) (194 'nb_cutoff' (REAL 8 0 0
REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (195 'switch_cutoff1' (REAL 8 0 0
REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (196 'switch_cutoff2' (REAL 8 0 0
REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (197 'mm_atomtype' (INTEGER 4 0 0
INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (198 'supported'
(LOGICAL 4 0 0 LOGICAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (
199 'atomic_number' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (200 'lj_r' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1
1 0 UNKNOWN-ACCESS ()) (201 'lj_epsilon' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ())) PUBLIC () 0 0)
202 'qmmm_scratch' 'qmmm_module' 'qmmm_scratch' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN) (DERIVED 203 0 0 DERIVED ())
0 0 () () 0 () () 0 0)
203 'qmmm_scratch_structure' 'qmmm_module' 'qmmm_scratch_structure' 1 (
(DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((204 'matsize_red_scratch' (REAL
8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (205
'qm_pme_scratch' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (206 'mat_diag_workspace' (REAL 8 0 0 REAL ()) (2
DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (207
'pdiag_scr_norbs_norbs' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1
1 0 UNKNOWN-ACCESS ()) (208 'pdiag_scr_noccupied_norbs' (REAL 8 0 0 REAL
()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (209 'pdiag_vectmp1'
(REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (210
'pdiag_vectmp2' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (211 'pdiag_vectmp3' (REAL 8 0 0 REAL ()) (1 DEFERRED
() ()) 1 1 0 UNKNOWN-ACCESS ()) (212 'pdiag_vecjs' (INTEGER 4 0 0
INTEGER ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (213
'lapack_dc_real_scr' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (214 'lapack_dc_int_scr' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (215 'qm_real_scratch' (REAL 8
0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (216
'qm_int_scratch' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (217 'lapack_dc_real_scr_aloc' (INTEGER 4 0 0 INTEGER
()) () 0 0 0 UNKNOWN-ACCESS ()) (218 'lapack_dc_int_scr_aloc' (INTEGER 4
0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (219 'qm_mm_pairs_allocated'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())) PUBLIC () 0 0)
220 'qmmm_struct' 'qmmm_module' 'qmmm_struct' 1 ((VARIABLE
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN EXPLICIT-SAVE) (DERIVED 221 0 0
DERIVED ()) 0 0 () () 0 () () 0 0)
222 'qmmm_vsolv' 'qmmm_module' 'qmmm_vsolv' 1 ((VARIABLE UNKNOWN-INTENT
UNKNOWN-PROC UNKNOWN EXPLICIT-SAVE) (DERIVED 223 0 0 DERIVED ()) 0 0 ()
() 0 () () 0 0)
224 'qmsort' 'qmmm_module' 'qmsort' 1 ((PROCEDURE UNKNOWN-INTENT
MODULE-PROC DECL UNKNOWN SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 225 0 (
226) () 0 () () 0 0)
227 'validate_qm_atoms' 'qmmm_module' 'validate_qm_atoms' 1 ((PROCEDURE
UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN SUBROUTINE) (UNKNOWN 0 0 0
UNKNOWN ()) 228 0 (229 230 231) () 0 () () 0 0)
221 'qmmm_struct_type' 'qmmm_struct_module' 'qmmm_struct_type' 1 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((232 'enuclr_qmqm' (REAL 8 0 0
REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (233 'enuclr_qmmm' (REAL 8 0 0 REAL
()) () 0 0 0 UNKNOWN-ACCESS ()) (234 'elec_eng' (REAL 8 0 0 REAL ()) ()
0 0 0 UNKNOWN-ACCESS ()) (235 'coulombic_eng' (REAL 8 0 0 REAL ()) () 0
0 0 UNKNOWN-ACCESS ()) (236 'dcorrection' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (237 'hcorrection' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (238 'qm_resp_charges' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(239 'qm_resp_charge_sum' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ())
(240 'mm_link_pair_resp_charges' (REAL 8 0 0 REAL ()) (1 DEFERRED () ())
1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (241
'mm_link_pair_saved_coords' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ())
1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (242 'qm_coords'
(REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (243 'scaled_mm_charges' (REAL 8 0 0
REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (244 'dxyzqm' (REAL 8 0 0 REAL ()) (2 DEFERRED () () ()
()) 1 1 0 UNKNOWN-ACCESS ()) (245 'dxyzcl' (REAL 8 0 0 REAL ()) (2
DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (246 'qm_xcrd' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1
0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (247
'switched_mmpot' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (248 'natom' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (249 'nquant' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (250 'core_nquant'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (251
'buffer_nquant' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
252 'nlink' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (253
'nquant_nlink' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
254 'qm_ntypes' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
255 'noshake_overlap' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS
()) (256 'qm_type_id' (INTEGER 4 0 0 INTEGER ()) (1 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'86')) 1 0 0 UNKNOWN-ACCESS ()) (257 'qm_atom_type' (INTEGER 4 0 0
INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (258 'link_pairs' (INTEGER 4 0 0 INTEGER ()) (2 DEFERRED
() () () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
259 'iqmatoms' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (260 'core_iqmatoms'
(INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (261 'buffer_iqmatoms' (INTEGER 4 0
0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0)) (262 'qm_nsubset' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (263 'core_nsubset' (INTEGER 4 0 0 INTEGER ()) () 0 0
0 UNKNOWN-ACCESS ()) (264 'buffer_nsubset' (INTEGER 4 0 0 INTEGER ()) ()
0 0 0 UNKNOWN-ACCESS ()) (265 'qm_subsetatoms' (INTEGER 4 0 0 INTEGER ())
(1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ())
0)) (266 'core_subsetatoms' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ())
1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (267
'buffer_subsetatoms' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (268 'center_nsubset'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (269
'center_subsetatoms' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (270
'iqm_atomic_numbers' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (271 'qm_mm_pairs' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (272
'qm_mm_pair_list' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (273
'qm_mm_pair_atom_numbers' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ())
1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (274
'num_qmmm_calls' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(275 'atom_mask' (LOGICAL 4 0 0 LOGICAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (276 'mm_link_mask'
(LOGICAL 4 0 0 LOGICAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (277 'qm_mm_first_call' (LOGICAL 4 0
0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (278 'fock_first_call' (
LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (279
'fock2_2atm_first_call' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (280 'qm2_allocate_e_repul_first_call' (LOGICAL 4 0 0
LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (281
'qm2_calc_rij_eqns_first_call' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (282 'qm2_scf_first_call' (LOGICAL 4 0 0 LOGICAL ())
() 0 0 0 UNKNOWN-ACCESS ()) (283 'zero_link_charges_first_call' (
LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (284
'adj_mm_link_pair_crd_first_call' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (285 'am1_or_pm3' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (286 'pddg_in_use' (LOGICAL 4 0 0 LOGICAL ()) () 0 0
0 UNKNOWN-ACCESS ()) (287 'mmcoords_contains_lnk_coords' (LOGICAL 4 0 0
LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (288 'pm3mmx_interface' (
LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (289 'abfqmmm' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 INTEGER ()) 0 '0')) (290 'hot_spot' (INTEGER 4 0 0 INTEGER ()) () 0 0
0 UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (291
'r_core_in' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (292
'r_core_out' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (293
'r_qm_in' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (294 'r_qm_out'
(REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (295 'r_buffer_in' (
REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (296 'r_buffer_out' (
REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (297 'cut_bond_list_file'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '256')))
() 0 0 0 UNKNOWN-ACCESS ()) (298 'oxidation_number_list_file' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '256')))
() 0 0 0 UNKNOWN-ACCESS ()) (299 'mom_cons_type' (INTEGER 4 0 0 INTEGER
()) () 0 0 0 UNKNOWN-ACCESS ()) (300 'mom_cons_region' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (301 'fix_atom_list' (INTEGER 4
0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (302 'solvent_atom_number' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (303
'selection_type' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(304 'center_type' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(305 'initial_selection_type' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (306 'max_bonds_per_atom' (INTEGER 4 0 0 INTEGER ())
() 0 0 0 UNKNOWN-ACCESS ()) (307 'n_max_recursive' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (308 'min_heavy_mass' (REAL 8 0
0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (309 'gamma_ln_qm' (REAL 8 0 0
REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (310 'read_idrst_file' (CHARACTER 1
0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '256'))) () 0 0 0
UNKNOWN-ACCESS ()) (311 'write_idrst_file' (CHARACTER 1 0 0 CHARACTER (
(CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '256'))) () 0 0 0 UNKNOWN-ACCESS
()) (312 'ntwidrst' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(313 'pdb_file' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '256'))) () 0 0 0 UNKNOWN-ACCESS ()) (314 'ntwpdb' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())) PUBLIC () 0 0)
17 'options' '' 'options' 16 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN UNKNOWN
DUMMY) (DERIVED 99 0 0 DERIVED ()) 0 0 () () 0 () () 0 0)
4 'qmmm_nml' '' 'qmmm_nml' 3 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN DUMMY) (DERIVED 315 0 0 DERIVED ()) 0 0 () () 0 () () 0 0)
5 'qmmm_struct' '' 'qmmm_struct' 3 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN DUMMY) (DERIVED 316 0 0 DERIVED ()) 0 0 () () 0 () () 0 0)
6 'natom' '' 'natom' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN DUMMY)
(INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () 0 0)
316 'qmmm_struct_type' 'qmmm_struct_module' 'qmmm_struct_type' 3 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((317 'enuclr_qmqm' (REAL 8 0 0
REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (318 'enuclr_qmmm' (REAL 8 0 0 REAL
()) () 0 0 0 UNKNOWN-ACCESS ()) (319 'elec_eng' (REAL 8 0 0 REAL ()) ()
0 0 0 UNKNOWN-ACCESS ()) (320 'coulombic_eng' (REAL 8 0 0 REAL ()) () 0
0 0 UNKNOWN-ACCESS ()) (321 'dcorrection' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (322 'hcorrection' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (323 'qm_resp_charges' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(324 'qm_resp_charge_sum' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ())
(325 'mm_link_pair_resp_charges' (REAL 8 0 0 REAL ()) (1 DEFERRED () ())
1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (326
'mm_link_pair_saved_coords' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ())
1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (327 'qm_coords'
(REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (328 'scaled_mm_charges' (REAL 8 0 0
REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (329 'dxyzqm' (REAL 8 0 0 REAL ()) (2 DEFERRED () () ()
()) 1 1 0 UNKNOWN-ACCESS ()) (330 'dxyzcl' (REAL 8 0 0 REAL ()) (2
DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (331 'qm_xcrd' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1
0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (332
'switched_mmpot' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (333 'natom' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (334 'nquant' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (335 'core_nquant'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (336
'buffer_nquant' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
337 'nlink' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (338
'nquant_nlink' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
339 'qm_ntypes' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
340 'noshake_overlap' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS
()) (341 'qm_type_id' (INTEGER 4 0 0 INTEGER ()) (1 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'86')) 1 0 0 UNKNOWN-ACCESS ()) (342 'qm_atom_type' (INTEGER 4 0 0
INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (343 'link_pairs' (INTEGER 4 0 0 INTEGER ()) (2 DEFERRED
() () () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
344 'iqmatoms' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (345 'core_iqmatoms'
(INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (346 'buffer_iqmatoms' (INTEGER 4 0
0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0)) (347 'qm_nsubset' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (348 'core_nsubset' (INTEGER 4 0 0 INTEGER ()) () 0 0
0 UNKNOWN-ACCESS ()) (349 'buffer_nsubset' (INTEGER 4 0 0 INTEGER ()) ()
0 0 0 UNKNOWN-ACCESS ()) (350 'qm_subsetatoms' (INTEGER 4 0 0 INTEGER ())
(1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ())
0)) (351 'core_subsetatoms' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ())
1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (352
'buffer_subsetatoms' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (353 'center_nsubset'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (354
'center_subsetatoms' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (355
'iqm_atomic_numbers' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (356 'qm_mm_pairs' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (357
'qm_mm_pair_list' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (358
'qm_mm_pair_atom_numbers' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ())
1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (359
'num_qmmm_calls' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(360 'atom_mask' (LOGICAL 4 0 0 LOGICAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (361 'mm_link_mask'
(LOGICAL 4 0 0 LOGICAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (362 'qm_mm_first_call' (LOGICAL 4 0
0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (363 'fock_first_call' (
LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (364
'fock2_2atm_first_call' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (365 'qm2_allocate_e_repul_first_call' (LOGICAL 4 0 0
LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (366
'qm2_calc_rij_eqns_first_call' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (367 'qm2_scf_first_call' (LOGICAL 4 0 0 LOGICAL ())
() 0 0 0 UNKNOWN-ACCESS ()) (368 'zero_link_charges_first_call' (
LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (369
'adj_mm_link_pair_crd_first_call' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (370 'am1_or_pm3' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (371 'pddg_in_use' (LOGICAL 4 0 0 LOGICAL ()) () 0 0
0 UNKNOWN-ACCESS ()) (372 'mmcoords_contains_lnk_coords' (LOGICAL 4 0 0
LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (373 'pm3mmx_interface' (
LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (374 'abfqmmm' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 INTEGER ()) 0 '0')) (375 'hot_spot' (INTEGER 4 0 0 INTEGER ()) () 0 0
0 UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (376
'r_core_in' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (377
'r_core_out' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (378
'r_qm_in' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (379 'r_qm_out'
(REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (380 'r_buffer_in' (
REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (381 'r_buffer_out' (
REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (382 'cut_bond_list_file'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '256')))
() 0 0 0 UNKNOWN-ACCESS ()) (383 'oxidation_number_list_file' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '256')))
() 0 0 0 UNKNOWN-ACCESS ()) (384 'mom_cons_type' (INTEGER 4 0 0 INTEGER
()) () 0 0 0 UNKNOWN-ACCESS ()) (385 'mom_cons_region' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (386 'fix_atom_list' (INTEGER 4
0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (387 'solvent_atom_number' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (388
'selection_type' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(389 'center_type' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(390 'initial_selection_type' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (391 'max_bonds_per_atom' (INTEGER 4 0 0 INTEGER ())
() 0 0 0 UNKNOWN-ACCESS ()) (392 'n_max_recursive' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (393 'min_heavy_mass' (REAL 8 0
0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (394 'gamma_ln_qm' (REAL 8 0 0
REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (395 'read_idrst_file' (CHARACTER 1
0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '256'))) () 0 0 0
UNKNOWN-ACCESS ()) (396 'write_idrst_file' (CHARACTER 1 0 0 CHARACTER (
(CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '256'))) () 0 0 0 UNKNOWN-ACCESS
()) (397 'ntwidrst' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(398 'pdb_file' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '256'))) () 0 0 0 UNKNOWN-ACCESS ()) (399 'ntwpdb' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())) PUBLIC () 0 0)
315 'qmmm_nml_type' 'qmmm_nml_module' 'qmmm_nml_type' 3 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (UNKNOWN 0 0 0
UNKNOWN ()) 0 0 () () 0 ((400 'qmcut' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (401 'qmcut2' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (402 'lnk_dis' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (403 'scfconv' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (404 'density_conv' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (405 'errconv' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (406 'ndiis_matrices' (INTEGER 4 0 0 INTEGER ()) () 0
0 0 UNKNOWN-ACCESS ()) (407 'ndiis_attempts' (INTEGER 4 0 0 INTEGER ())
() 0 0 0 UNKNOWN-ACCESS ()) (408 'pseudo_diag_criteria' (REAL 8 0 0 REAL
()) () 0 0 0 UNKNOWN-ACCESS ()) (409 'lnk_atomic_no' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (410 'lnk_method' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (411 'qmgb' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (412 'qmtheory' (DERIVED 413 0 0
DERIVED ()) () 0 0 0 UNKNOWN-ACCESS ()) (414 'qmcharge' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (415 'corecharge' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (416 'buffercharge' (INTEGER 4 0
0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (417 'spin' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (418 'verbosity' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (419 'itrmax' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (420 'qmshake' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (421 'kmaxqx' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (422 'kmaxqy' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (423 'kmaxqz' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (424 'ksqmaxq' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (425 'kappa' (REAL 8 0 0 REAL ())
() 0 0 0 UNKNOWN-ACCESS ()) (426 'qm_ewald' (INTEGER 4 0 0 INTEGER ()) ()
0 0 0 UNKNOWN-ACCESS ()) (427 'qmmm_int' (INTEGER 4 0 0 INTEGER ()) () 0
0 0 UNKNOWN-ACCESS ()) (428 'qmmm_switch' (LOGICAL 4 0 0 LOGICAL ()) ()
0 0 0 UNKNOWN-ACCESS ()) (429 'r_switch_lo' (REAL 8 0 0 REAL ()) () 0 0
0 UNKNOWN-ACCESS ()) (430 'r_switch_hi' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (431 'adjust_q' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (432 'diag_routine' (INTEGER 4 0 0 INTEGER ()) () 0 0
0 UNKNOWN-ACCESS ()) (433 'density_predict' (INTEGER 4 0 0 INTEGER ()) ()
0 0 0 UNKNOWN-ACCESS ()) (434 'vsolv' (INTEGER 4 0 0 INTEGER ()) () 0 0
0 UNKNOWN-ACCESS ()) (435 'fock_predict' (INTEGER 4 0 0 INTEGER ()) () 0
0 0 UNKNOWN-ACCESS ()) (436 'fockp_d1' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (437 'fockp_d2' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (438 'fockp_d3' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (439 'fockp_d4' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (440 'idc' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (441 'divpb' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (442 'nquant' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (443 'iqmatoms' (INTEGER 4 0 0 INTEGER ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(444 'chg_lambda' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (445
'dftb_maxiter' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
446 'dftb_disper' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(447 'dftb_chg' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
448 'printdipole' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(449 'print_eigenvalues' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (450 'dftb_telec' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (451 'dftb_telec_step' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (452 'dftb_3rd_order' (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '256'))) () 0 0 0 UNKNOWN-ACCESS ())
(453 'ifqnt' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (
454 'qmqm_analyt' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ())
(455 'tight_p_conv' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ())
(456 'printcharges' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ())
(457 'printbondorders' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (458 'peptide_corr' (LOGICAL 4 0 0 LOGICAL ()) () 0 0
0 UNKNOWN-ACCESS ()) (459 'qmqm_erep_incore' (LOGICAL 4 0 0 LOGICAL ())
() 0 0 0 UNKNOWN-ACCESS ()) (460 'allow_pseudo_diag' (LOGICAL 4 0 0
LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (461 'qmmmrij_incore' (LOGICAL 4
0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (462 'writepdb' (LOGICAL 4 0
0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (463 'qm_pme' (LOGICAL 4 0 0
LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (464 'damp' (REAL 8 0 0 REAL ())
() 0 0 0 UNKNOWN-ACCESS ()) (465 'vshift' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ())) PUBLIC () 0 0)
13 'qmmm_vsolv' '' 'qmmm_vsolv' 10 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN DUMMY) (DERIVED 466 0 0 DERIVED ()) 0 0 () () 0 () () 0 0)
14 'qm2_params' '' 'qm2_params' 10 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN DUMMY) (DERIVED 467 0 0 DERIVED ()) 0 0 () () 0 () () 0 0)
226 'iqmatoms' '' 'iqmatoms' 225 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1
ASSUMED_SIZE (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () 0 0)
22 'atom_mass' '' 'atom_mass' 20 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN DUMMY) (REAL 8 0 0 REAL ()) 0 0 () () 0 () () 0 0)
21 'atom_name' '' 'atom_name' 20 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN DUMMY) (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '4'))) 0 0 () () 0 () () 0 0)
23 'atomic_number' '' 'atomic_number' 20 ((VARIABLE OUT UNKNOWN-PROC
UNKNOWN UNKNOWN DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () 0 0)
24 'errorflag' '' 'errorflag' 20 ((VARIABLE OUT UNKNOWN-PROC UNKNOWN
UNKNOWN OPTIONAL DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () () 0 0)
229 'iqmatoms' '' 'iqmatoms' 228 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 EXPLICIT (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (VARIABLE (INTEGER 4 0 0
INTEGER ()) 0 230 ())) 0 () () 0 0)
230 'nquant' '' 'nquant' 228 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () 0 0)
231 'natom' '' 'natom' 228 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () 0 0)
187 'qmmm_nml_type' 'qmmm_nml_module' 'qmmm_nml_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (UNKNOWN 0 0 0
UNKNOWN ()) 0 0 () () 0 ((468 'qmcut' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (469 'qmcut2' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (470 'lnk_dis' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (471 'scfconv' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (472 'density_conv' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (473 'errconv' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (474 'ndiis_matrices' (INTEGER 4 0 0 INTEGER ()) () 0
0 0 UNKNOWN-ACCESS ()) (475 'ndiis_attempts' (INTEGER 4 0 0 INTEGER ())
() 0 0 0 UNKNOWN-ACCESS ()) (476 'pseudo_diag_criteria' (REAL 8 0 0 REAL
()) () 0 0 0 UNKNOWN-ACCESS ()) (477 'lnk_atomic_no' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (478 'lnk_method' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (479 'qmgb' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (480 'qmtheory' (DERIVED 481 0 0
DERIVED ()) () 0 0 0 UNKNOWN-ACCESS ()) (482 'qmcharge' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (483 'corecharge' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (484 'buffercharge' (INTEGER 4 0
0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (485 'spin' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (486 'verbosity' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (487 'itrmax' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (488 'qmshake' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (489 'kmaxqx' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (490 'kmaxqy' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (491 'kmaxqz' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (492 'ksqmaxq' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (493 'kappa' (REAL 8 0 0 REAL ())
() 0 0 0 UNKNOWN-ACCESS ()) (494 'qm_ewald' (INTEGER 4 0 0 INTEGER ()) ()
0 0 0 UNKNOWN-ACCESS ()) (495 'qmmm_int' (INTEGER 4 0 0 INTEGER ()) () 0
0 0 UNKNOWN-ACCESS ()) (496 'qmmm_switch' (LOGICAL 4 0 0 LOGICAL ()) ()
0 0 0 UNKNOWN-ACCESS ()) (497 'r_switch_lo' (REAL 8 0 0 REAL ()) () 0 0
0 UNKNOWN-ACCESS ()) (498 'r_switch_hi' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (499 'adjust_q' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (500 'diag_routine' (INTEGER 4 0 0 INTEGER ()) () 0 0
0 UNKNOWN-ACCESS ()) (501 'density_predict' (INTEGER 4 0 0 INTEGER ()) ()
0 0 0 UNKNOWN-ACCESS ()) (502 'vsolv' (INTEGER 4 0 0 INTEGER ()) () 0 0
0 UNKNOWN-ACCESS ()) (503 'fock_predict' (INTEGER 4 0 0 INTEGER ()) () 0
0 0 UNKNOWN-ACCESS ()) (504 'fockp_d1' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (505 'fockp_d2' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (506 'fockp_d3' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (507 'fockp_d4' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (508 'idc' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (509 'divpb' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (510 'nquant' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (511 'iqmatoms' (INTEGER 4 0 0 INTEGER ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(512 'chg_lambda' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (513
'dftb_maxiter' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
514 'dftb_disper' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(515 'dftb_chg' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
516 'printdipole' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(517 'print_eigenvalues' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (518 'dftb_telec' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (519 'dftb_telec_step' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (520 'dftb_3rd_order' (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '256'))) () 0 0 0 UNKNOWN-ACCESS ())
(521 'ifqnt' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (
522 'qmqm_analyt' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ())
(523 'tight_p_conv' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ())
(524 'printcharges' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ())
(525 'printbondorders' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (526 'peptide_corr' (LOGICAL 4 0 0 LOGICAL ()) () 0 0
0 UNKNOWN-ACCESS ()) (527 'qmqm_erep_incore' (LOGICAL 4 0 0 LOGICAL ())
() 0 0 0 UNKNOWN-ACCESS ()) (528 'allow_pseudo_diag' (LOGICAL 4 0 0
LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (529 'qmmmrij_incore' (LOGICAL 4
0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (530 'writepdb' (LOGICAL 4 0
0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (531 'qm_pme' (LOGICAL 4 0 0
LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (532 'damp' (REAL 8 0 0 REAL ())
() 0 0 0 UNKNOWN-ACCESS ()) (533 'vshift' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ())) PUBLIC () 0 0)
481 'qmtheorytype' 'qmmm_qmtheorymodule' 'qmtheorytype' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((534 'pm3' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (535 'am1' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (536 'am1d' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (537 'mndo' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (538 'mndod' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (539 'pddgpm3' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (540 'pddgmndo' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (541 'pm3carb1' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (542 'pm3znb' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (543 'dftb' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (544 'rm1' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (545 'pddgpm3_08' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (546 'pm6' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (547 'dispersion' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (548 'dispersion_hydrogenplus' (LOGICAL 4 0 0 LOGICAL
()) () 0 0 0 UNKNOWN-ACCESS ()) (549 'pm3mais' (LOGICAL 4 0 0 LOGICAL ())
() 0 0 0 UNKNOWN-ACCESS ()) (550 'extern' (LOGICAL 4 0 0 LOGICAL ()) ()
0 0 0 UNKNOWN-ACCESS ()) (551 'sebomd' (LOGICAL 4 0 0 LOGICAL ()) () 0 0
0 UNKNOWN-ACCESS ())) PUBLIC () 0 0)
223 'qmmm_vsolv_type' 'qmmm_vsolv_module' 'qmmm_vsolv_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (UNKNOWN 0 0 0
UNKNOWN ()) 0 0 () () 0 ((552 'debug' (LOGICAL 4 0 0 LOGICAL ()) () 0 0
0 UNKNOWN-ACCESS ()) (553 'verbosity' (INTEGER 4 0 0 INTEGER ()) () 0 0
0 UNKNOWN-ACCESS ()) (554 'recalculate' (LOGICAL 4 0 0 LOGICAL ()) () 0
0 0 UNKNOWN-ACCESS ()) (555 'nearest_qm_solvent' (INTEGER 4 0 0 INTEGER
()) () 0 0 0 UNKNOWN-ACCESS ()) (556 'nearest_qm_solvent_fq' (INTEGER 4
0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (557
'nearest_qm_solvent_resname' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '4'))) () 0 0 0 UNKNOWN-ACCESS ()) (558
'nearest_qm_solvent_center_id' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (559 'qm_center_atom_id' (INTEGER 4 0 0 INTEGER ()) ()
0 0 0 UNKNOWN-ACCESS ()) (560 'fixed_nquant' (INTEGER 4 0 0 INTEGER ())
() 0 0 0 UNKNOWN-ACCESS ()) (561 'nsolv_res' (INTEGER 4 0 0 INTEGER ())
() 0 0 0 UNKNOWN-ACCESS ()) (562 'natom_solv_res' (INTEGER 4 0 0 INTEGER
()) () 0 0 0 UNKNOWN-ACCESS ()) (563 'fixed_iqmatoms' (INTEGER 4 0 0
INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (564 'solvent_pointers' (INTEGER 4 0 0 INTEGER ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(565 'nearest_solvent_pointers' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED ()
()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (566
'nearest_solvent_pointers_prev' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED ()
()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (567
'nearest_solvent_distances' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1
0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (568 'prmtop_numbnd'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (569 'nbonh' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (570 'iibh' (
INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0)) (571 'ijbh' (INTEGER 4 0 0 INTEGER ()) (
1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ())
0)) (572 'icbh' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (573 'nbona' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (574 'iiba' (
INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0)) (575 'ijba' (INTEGER 4 0 0 INTEGER ()) (
1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ())
0)) (576 'icba' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (577 'ntheth' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (578 'iith' (
INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0)) (579 'ijth' (INTEGER 4 0 0 INTEGER ()) (
1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ())
0)) (580 'ikth' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (581 'icth' (
INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0)) (582 'ntheta' (INTEGER 4 0 0 INTEGER ())
() 0 0 0 UNKNOWN-ACCESS ()) (583 'iita' (INTEGER 4 0 0 INTEGER ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(584 'ijta' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (585 'ikta' (
INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0)) (586 'icta' (INTEGER 4 0 0 INTEGER ()) (
1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ())
0)) (587 'nphih' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(588 'iiph' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (589 'ijph' (
INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0)) (590 'ikph' (INTEGER 4 0 0 INTEGER ()) (
1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ())
0)) (591 'ilph' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (592 'icph' (
INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0)) (593 'nphia' (INTEGER 4 0 0 INTEGER ()) ()
0 0 0 UNKNOWN-ACCESS ()) (594 'iipa' (INTEGER 4 0 0 INTEGER ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(595 'ijpa' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (596 'ikpa' (
INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0)) (597 'ilpa' (INTEGER 4 0 0 INTEGER ()) (
1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ())
0)) (598 'icpa' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC () 0 0)
413 'qmtheorytype' 'qmmm_qmtheorymodule' 'qmtheorytype' 3 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((599 'pm3' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (600 'am1' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (601 'am1d' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (602 'mndo' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (603 'mndod' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (604 'pddgpm3' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (605 'pddgmndo' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (606 'pm3carb1' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (607 'pm3znb' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (608 'dftb' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (609 'rm1' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (610 'pddgpm3_08' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (611 'pm6' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (612 'dispersion' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (613 'dispersion_hydrogenplus' (LOGICAL 4 0 0 LOGICAL
()) () 0 0 0 UNKNOWN-ACCESS ()) (614 'pm3mais' (LOGICAL 4 0 0 LOGICAL ())
() 0 0 0 UNKNOWN-ACCESS ()) (615 'extern' (LOGICAL 4 0 0 LOGICAL ()) ()
0 0 0 UNKNOWN-ACCESS ()) (616 'sebomd' (LOGICAL 4 0 0 LOGICAL ()) () 0 0
0 UNKNOWN-ACCESS ())) PUBLIC () 0 0)
11 'qmmm_nml' '' 'qmmm_nml' 10 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN DUMMY) (DERIVED 617 0 0 DERIVED ()) 0 0 () () 0 () () 0 0)
12 'qmmm_struct' '' 'qmmm_struct' 10 ((VARIABLE INOUT UNKNOWN-PROC
UNKNOWN UNKNOWN DUMMY) (DERIVED 618 0 0 DERIVED ()) 0 0 () () 0 () () 0
0)
618 'qmmm_struct_type' 'qmmm_struct_module' 'qmmm_struct_type' 10 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((619 'enuclr_qmqm' (REAL 8 0 0
REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (620 'enuclr_qmmm' (REAL 8 0 0 REAL
()) () 0 0 0 UNKNOWN-ACCESS ()) (621 'elec_eng' (REAL 8 0 0 REAL ()) ()
0 0 0 UNKNOWN-ACCESS ()) (622 'coulombic_eng' (REAL 8 0 0 REAL ()) () 0
0 0 UNKNOWN-ACCESS ()) (623 'dcorrection' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (624 'hcorrection' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (625 'qm_resp_charges' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(626 'qm_resp_charge_sum' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ())
(627 'mm_link_pair_resp_charges' (REAL 8 0 0 REAL ()) (1 DEFERRED () ())
1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (628
'mm_link_pair_saved_coords' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ())
1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (629 'qm_coords'
(REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (630 'scaled_mm_charges' (REAL 8 0 0
REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (631 'dxyzqm' (REAL 8 0 0 REAL ()) (2 DEFERRED () () ()
()) 1 1 0 UNKNOWN-ACCESS ()) (632 'dxyzcl' (REAL 8 0 0 REAL ()) (2
DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (633 'qm_xcrd' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1
0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (634
'switched_mmpot' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (635 'natom' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (636 'nquant' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (637 'core_nquant'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (638
'buffer_nquant' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
639 'nlink' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (640
'nquant_nlink' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
641 'qm_ntypes' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
642 'noshake_overlap' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS
()) (643 'qm_type_id' (INTEGER 4 0 0 INTEGER ()) (1 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'86')) 1 0 0 UNKNOWN-ACCESS ()) (644 'qm_atom_type' (INTEGER 4 0 0
INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (645 'link_pairs' (INTEGER 4 0 0 INTEGER ()) (2 DEFERRED
() () () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (
646 'iqmatoms' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (647 'core_iqmatoms'
(INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (648 'buffer_iqmatoms' (INTEGER 4 0
0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0
0 UNKNOWN ()) 0)) (649 'qm_nsubset' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (650 'core_nsubset' (INTEGER 4 0 0 INTEGER ()) () 0 0
0 UNKNOWN-ACCESS ()) (651 'buffer_nsubset' (INTEGER 4 0 0 INTEGER ()) ()
0 0 0 UNKNOWN-ACCESS ()) (652 'qm_subsetatoms' (INTEGER 4 0 0 INTEGER ())
(1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ())
0)) (653 'core_subsetatoms' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ())
1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (654
'buffer_subsetatoms' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (655 'center_nsubset'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (656
'center_subsetatoms' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (657
'iqm_atomic_numbers' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (658 'qm_mm_pairs' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (659
'qm_mm_pair_list' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (660
'qm_mm_pair_atom_numbers' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ())
1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (661
'num_qmmm_calls' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(662 'atom_mask' (LOGICAL 4 0 0 LOGICAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (663 'mm_link_mask'
(LOGICAL 4 0 0 LOGICAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (664 'qm_mm_first_call' (LOGICAL 4 0
0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (665 'fock_first_call' (
LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (666
'fock2_2atm_first_call' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (667 'qm2_allocate_e_repul_first_call' (LOGICAL 4 0 0
LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (668
'qm2_calc_rij_eqns_first_call' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (669 'qm2_scf_first_call' (LOGICAL 4 0 0 LOGICAL ())
() 0 0 0 UNKNOWN-ACCESS ()) (670 'zero_link_charges_first_call' (
LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (671
'adj_mm_link_pair_crd_first_call' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (672 'am1_or_pm3' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (673 'pddg_in_use' (LOGICAL 4 0 0 LOGICAL ()) () 0 0
0 UNKNOWN-ACCESS ()) (674 'mmcoords_contains_lnk_coords' (LOGICAL 4 0 0
LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (675 'pm3mmx_interface' (
LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (676 'abfqmmm' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0
0 INTEGER ()) 0 '0')) (677 'hot_spot' (INTEGER 4 0 0 INTEGER ()) () 0 0
0 UNKNOWN-ACCESS (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '0')) (678
'r_core_in' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (679
'r_core_out' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (680
'r_qm_in' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (681 'r_qm_out'
(REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (682 'r_buffer_in' (
REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (683 'r_buffer_out' (
REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (684 'cut_bond_list_file'
(CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '256')))
() 0 0 0 UNKNOWN-ACCESS ()) (685 'oxidation_number_list_file' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '256')))
() 0 0 0 UNKNOWN-ACCESS ()) (686 'mom_cons_type' (INTEGER 4 0 0 INTEGER
()) () 0 0 0 UNKNOWN-ACCESS ()) (687 'mom_cons_region' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (688 'fix_atom_list' (INTEGER 4
0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (689 'solvent_atom_number' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (690
'selection_type' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(691 'center_type' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(692 'initial_selection_type' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (693 'max_bonds_per_atom' (INTEGER 4 0 0 INTEGER ())
() 0 0 0 UNKNOWN-ACCESS ()) (694 'n_max_recursive' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (695 'min_heavy_mass' (REAL 8 0
0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (696 'gamma_ln_qm' (REAL 8 0 0
REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (697 'read_idrst_file' (CHARACTER 1
0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '256'))) () 0 0 0
UNKNOWN-ACCESS ()) (698 'write_idrst_file' (CHARACTER 1 0 0 CHARACTER (
(CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '256'))) () 0 0 0 UNKNOWN-ACCESS
()) (699 'ntwidrst' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(700 'pdb_file' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '256'))) () 0 0 0 UNKNOWN-ACCESS ()) (701 'ntwpdb' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())) PUBLIC () 0 0)
466 'qmmm_vsolv_type' 'qmmm_vsolv_module' 'qmmm_vsolv_type' 10 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((702 'debug' (LOGICAL 4 0 0
LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (703 'verbosity' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (704 'recalculate' (LOGICAL 4 0
0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (705 'nearest_qm_solvent' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (706
'nearest_qm_solvent_fq' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (707 'nearest_qm_solvent_resname' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '4'))) () 0 0 0
UNKNOWN-ACCESS ()) (708 'nearest_qm_solvent_center_id' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (709 'qm_center_atom_id' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (710 'fixed_nquant'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (711 'nsolv_res'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (712
'natom_solv_res' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(713 'fixed_iqmatoms' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1
0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (714
'solvent_pointers' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (715
'nearest_solvent_pointers' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ())
1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (716
'nearest_solvent_pointers_prev' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED ()
()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (717
'nearest_solvent_distances' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1
0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (718 'prmtop_numbnd'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (719 'nbonh' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (720 'iibh' (
INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0)) (721 'ijbh' (INTEGER 4 0 0 INTEGER ()) (
1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ())
0)) (722 'icbh' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (723 'nbona' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (724 'iiba' (
INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0)) (725 'ijba' (INTEGER 4 0 0 INTEGER ()) (
1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ())
0)) (726 'icba' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (727 'ntheth' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (728 'iith' (
INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0)) (729 'ijth' (INTEGER 4 0 0 INTEGER ()) (
1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ())
0)) (730 'ikth' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (731 'icth' (
INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0)) (732 'ntheta' (INTEGER 4 0 0 INTEGER ())
() 0 0 0 UNKNOWN-ACCESS ()) (733 'iita' (INTEGER 4 0 0 INTEGER ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(734 'ijta' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (735 'ikta' (
INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0)) (736 'icta' (INTEGER 4 0 0 INTEGER ()) (
1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ())
0)) (737 'nphih' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(738 'iiph' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (739 'ijph' (
INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0)) (740 'ikph' (INTEGER 4 0 0 INTEGER ()) (
1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ())
0)) (741 'ilph' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (742 'icph' (
INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0)) (743 'nphia' (INTEGER 4 0 0 INTEGER ()) ()
0 0 0 UNKNOWN-ACCESS ()) (744 'iipa' (INTEGER 4 0 0 INTEGER ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(745 'ijpa' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (746 'ikpa' (
INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0)) (747 'ilpa' (INTEGER 4 0 0 INTEGER ()) (
1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ())
0)) (748 'icpa' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))) PUBLIC () 0 0)
27 'qm2_params_type' 'qm2_params_module' 'qm2_params_type' 1 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (UNKNOWN 0 0 0
UNKNOWN ()) 0 0 () () 0 ((749 'tot_heat_form' (REAL 8 0 0 REAL ()) () 0
0 0 UNKNOWN-ACCESS ()) (750 'sp_quantum_number' (INTEGER 4 0 0 INTEGER ())
(1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (751 'd_quantum_number' (
INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (
752 'gss' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ())
(753 'hsp' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS
()) (754 'hpp' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (755 'dd' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ())
1 1 0 UNKNOWN-ACCESS ()) (756 'po' (REAL 8 0 0 REAL ()) (2 DEFERRED () ()
() ()) 1 1 0 UNKNOWN-ACCESS ()) (757 'core_chg' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (758 'orb_elec_ke' (REAL 8 0 0
REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (759 'betasas'
(REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (
760 'betasap' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0
UNKNOWN-ACCESS ()) (761 'betasad' (REAL 8 0 0 REAL ()) (2 DEFERRED () ()
() ()) 1 1 0 UNKNOWN-ACCESS ()) (762 'betapap' (REAL 8 0 0 REAL ()) (2
DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (763 'betapad' (REAL 8 0
0 REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (764
'betadad' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0
UNKNOWN-ACCESS ()) (765 'gnn' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1
1 0 UNKNOWN-ACCESS ()) (766 'rho_core' (REAL 8 0 0 REAL ()) (1 DEFERRED
() ()) 1 1 0 UNKNOWN-ACCESS ()) (767 'f0sd' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (768 'g2sd' (REAL 8 0 0 REAL ())
(1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (769 'fn1' (REAL 8 0 0 REAL
()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (770 'fn2' (REAL 8
0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (771 'fn3'
(REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (
772 'onec2elec_params' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1
0 UNKNOWN-ACCESS ()) (773 'multip_2c_elec_params' (REAL 8 0 0 REAL ()) (
2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (774 'cc_exp_params' (
REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (775
'pm6_alpab' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0
UNKNOWN-ACCESS ()) (776 'pm6_xab' (REAL 8 0 0 REAL ()) (2 DEFERRED () ()
() ()) 1 1 0 UNKNOWN-ACCESS ()) (777 'pm3mais_alpab' (REAL 8 0 0 REAL ())
(3 DEFERRED () () () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (778
'pm3mais_betab' (REAL 8 0 0 REAL ()) (3 DEFERRED () () () () () ()) 1 1
0 UNKNOWN-ACCESS ()) (779 'pm3mais_gamab' (REAL 8 0 0 REAL ()) (3
DEFERRED () () () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (780
's_orb_exp_by_type' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (781 'p_orb_exp_by_type' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (782 'd_orb_exp_by_type' (REAL
8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (783
's_orb_exp_tail_by_type' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (784 'p_orb_exp_tail_by_type' (REAL 8 0 0 REAL ()) (
1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (785 'd_orb_exp_tail_by_type'
(REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (786
'pddge1' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ())
(787 'pddge2' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (788 'scale_factor1_pm3mmx' (REAL 8 0 0 REAL ()) (2
DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (789 'scale_factor2_pm3mmx'
(REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (
790 'rho_pm3mmx' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (791 'atom_orb_zz_sxs_over_sas' (REAL 8 0 0 REAL ())
(4 DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (792
'atom_orb_zz_sxp_over_sap' (REAL 8 0 0 REAL ()) (4 DEFERRED () () () ()
() () () ()) 1 1 0 UNKNOWN-ACCESS ()) (793 'atom_orb_zz_sxd_over_sad' (
REAL 8 0 0 REAL ()) (4 DEFERRED () () () () () () () ()) 1 1 0
UNKNOWN-ACCESS ()) (794 'atom_orb_zz_pxp_over_pap' (REAL 8 0 0 REAL ())
(4 DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (795
'atom_orb_zz_pxd_over_pad' (REAL 8 0 0 REAL ()) (4 DEFERRED () () () ()
() () () ()) 1 1 0 UNKNOWN-ACCESS ()) (796 'atom_orb_zz_dxd_over_dad' (
REAL 8 0 0 REAL ()) (4 DEFERRED () () () () () () () ()) 1 1 0
UNKNOWN-ACCESS ()) (797 'atom_orb_ss_eqn' (REAL 8 0 0 REAL ()) (4
DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (798
'atom_orb_sp_ovlp' (REAL 8 0 0 REAL ()) (4 DEFERRED () () () () () () ()
()) 1 1 0 UNKNOWN-ACCESS ()) (799 'atom_orb_sd_ovlp' (REAL 8 0 0 REAL ())
(4 DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (800
'atom_orb_pd_ovlp' (REAL 8 0 0 REAL ()) (4 DEFERRED () () () () () () ()
()) 1 1 0 UNKNOWN-ACCESS ()) (801 'atom_orb_pp_ovlp_inj' (REAL 8 0 0
REAL ()) (4 DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (
802 'atom_orb_pp_ovlp_ieqj1' (REAL 8 0 0 REAL ()) (4 DEFERRED () () () ()
() () () ()) 1 1 0 UNKNOWN-ACCESS ()) (803 'atom_orb_pp_ovlp_ieqj2' (
REAL 8 0 0 REAL ()) (4 DEFERRED () () () () () () () ()) 1 1 0
UNKNOWN-ACCESS ()) (804 'atom_orb_dd_ovlp_inj' (REAL 8 0 0 REAL ()) (4
DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (805
'atom_orb_dd_ovlp_ieqj1' (REAL 8 0 0 REAL ()) (4 DEFERRED () () () () ()
() () ()) 1 1 0 UNKNOWN-ACCESS ()) (806 'atom_orb_dd_ovlp_ieqj2' (REAL 8
0 0 REAL ()) (4 DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS ())
(807 'atom_orb_ss_eqn_adb' (REAL 8 0 0 REAL ()) (4 DEFERRED () () () ()
() () () ()) 1 1 0 UNKNOWN-ACCESS ()) (808 'atom_orb_sp_eqn_xy' (REAL 8
0 0 REAL ()) (4 DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS ())
(809 'atom_orb_sp_eqn_xx1' (REAL 8 0 0 REAL ()) (4 DEFERRED () () () ()
() () () ()) 1 1 0 UNKNOWN-ACCESS ()) (810 'atom_orb_sp_eqn_xx2' (REAL 8
0 0 REAL ()) (4 DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS ())
(811 'atom_orb_pp_eqn_xxy1' (REAL 8 0 0 REAL ()) (4 DEFERRED () () () ()
() () () ()) 1 1 0 UNKNOWN-ACCESS ()) (812 'atom_orb_pp_eqn_xxy2' (REAL
8 0 0 REAL ()) (4 DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS
()) (813 'pddg_term1' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1
0 UNKNOWN-ACCESS ()) (814 'pddg_term2' (REAL 8 0 0 REAL ()) (2 DEFERRED
() () () ()) 1 1 0 UNKNOWN-ACCESS ()) (815 'pddg_term3' (REAL 8 0 0 REAL
()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (816 'pddg_term4'
(REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (
817 'natomic_orbs' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (818 'orb_loc' (INTEGER 4 0 0 INTEGER ()) (2 DEFERRED
() () () ()) 1 1 0 UNKNOWN-ACCESS ()) (819 'pascal_tri1' (INTEGER 4 0 0
INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (820 'pascal_tri2'
(INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (
821 'num_fn' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (822 'qxd_supported' (LOGICAL 4 0 0 LOGICAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (823 'qxd_s' (REAL 8 0 0 REAL ())
(1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (824 'qxd_z0' (REAL 8 0 0
REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (825 'qxd_zq' (
REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (826
'qxd_d0' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ())
(827 'qxd_dq' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (828 'qxd_q0' (REAL 8 0 0 REAL ()) (1 DEFERRED () ())
1 1 0 UNKNOWN-ACCESS ()) (829 'qxd_qq' (REAL 8 0 0 REAL ()) (1 DEFERRED
() ()) 1 1 0 UNKNOWN-ACCESS ()) (830 'qxd_neff' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ())) PUBLIC () 0 0)
617 'qmmm_nml_type' 'qmmm_nml_module' 'qmmm_nml_type' 10 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (UNKNOWN 0 0 0
UNKNOWN ()) 0 0 () () 0 ((831 'qmcut' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (832 'qmcut2' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (833 'lnk_dis' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (834 'scfconv' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (835 'density_conv' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (836 'errconv' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (837 'ndiis_matrices' (INTEGER 4 0 0 INTEGER ()) () 0
0 0 UNKNOWN-ACCESS ()) (838 'ndiis_attempts' (INTEGER 4 0 0 INTEGER ())
() 0 0 0 UNKNOWN-ACCESS ()) (839 'pseudo_diag_criteria' (REAL 8 0 0 REAL
()) () 0 0 0 UNKNOWN-ACCESS ()) (840 'lnk_atomic_no' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (841 'lnk_method' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (842 'qmgb' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (843 'qmtheory' (DERIVED 844 0 0
DERIVED ()) () 0 0 0 UNKNOWN-ACCESS ()) (845 'qmcharge' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (846 'corecharge' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (847 'buffercharge' (INTEGER 4 0
0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (848 'spin' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (849 'verbosity' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (850 'itrmax' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (851 'qmshake' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (852 'kmaxqx' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (853 'kmaxqy' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (854 'kmaxqz' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (855 'ksqmaxq' (INTEGER 4 0 0
INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (856 'kappa' (REAL 8 0 0 REAL ())
() 0 0 0 UNKNOWN-ACCESS ()) (857 'qm_ewald' (INTEGER 4 0 0 INTEGER ()) ()
0 0 0 UNKNOWN-ACCESS ()) (858 'qmmm_int' (INTEGER 4 0 0 INTEGER ()) () 0
0 0 UNKNOWN-ACCESS ()) (859 'qmmm_switch' (LOGICAL 4 0 0 LOGICAL ()) ()
0 0 0 UNKNOWN-ACCESS ()) (860 'r_switch_lo' (REAL 8 0 0 REAL ()) () 0 0
0 UNKNOWN-ACCESS ()) (861 'r_switch_hi' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (862 'adjust_q' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (863 'diag_routine' (INTEGER 4 0 0 INTEGER ()) () 0 0
0 UNKNOWN-ACCESS ()) (864 'density_predict' (INTEGER 4 0 0 INTEGER ()) ()
0 0 0 UNKNOWN-ACCESS ()) (865 'vsolv' (INTEGER 4 0 0 INTEGER ()) () 0 0
0 UNKNOWN-ACCESS ()) (866 'fock_predict' (INTEGER 4 0 0 INTEGER ()) () 0
0 0 UNKNOWN-ACCESS ()) (867 'fockp_d1' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (868 'fockp_d2' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (869 'fockp_d3' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (870 'fockp_d4' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (871 'idc' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (872 'divpb' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (873 'nquant' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (874 'iqmatoms' (INTEGER 4 0 0 INTEGER ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(875 'chg_lambda' (REAL 8 0 0 REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (876
'dftb_maxiter' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
877 'dftb_disper' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(878 'dftb_chg' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
879 'printdipole' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(880 'print_eigenvalues' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (881 'dftb_telec' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (882 'dftb_telec_step' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (883 'dftb_3rd_order' (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '256'))) () 0 0 0 UNKNOWN-ACCESS ())
(884 'ifqnt' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (
885 'qmqm_analyt' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ())
(886 'tight_p_conv' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ())
(887 'printcharges' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ())
(888 'printbondorders' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (889 'peptide_corr' (LOGICAL 4 0 0 LOGICAL ()) () 0 0
0 UNKNOWN-ACCESS ()) (890 'qmqm_erep_incore' (LOGICAL 4 0 0 LOGICAL ())
() 0 0 0 UNKNOWN-ACCESS ()) (891 'allow_pseudo_diag' (LOGICAL 4 0 0
LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (892 'qmmmrij_incore' (LOGICAL 4
0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (893 'writepdb' (LOGICAL 4 0
0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (894 'qm_pme' (LOGICAL 4 0 0
LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (895 'damp' (REAL 8 0 0 REAL ())
() 0 0 0 UNKNOWN-ACCESS ()) (896 'vshift' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ())) PUBLIC () 0 0)
844 'qmtheorytype' 'qmmm_qmtheorymodule' 'qmtheorytype' 10 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN) (UNKNOWN 0 0 0 UNKNOWN ())
0 0 () () 0 ((897 'pm3' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (898 'am1' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (899 'am1d' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (900 'mndo' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (901 'mndod' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (902 'pddgpm3' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (903 'pddgmndo' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (904 'pm3carb1' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (905 'pm3znb' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (906 'dftb' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (907 'rm1' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (908 'pddgpm3_08' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (909 'pm6' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (910 'dispersion' (LOGICAL 4 0 0 LOGICAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (911 'dispersion_hydrogenplus' (LOGICAL 4 0 0 LOGICAL
()) () 0 0 0 UNKNOWN-ACCESS ()) (912 'pm3mais' (LOGICAL 4 0 0 LOGICAL ())
() 0 0 0 UNKNOWN-ACCESS ()) (913 'extern' (LOGICAL 4 0 0 LOGICAL ()) ()
0 0 0 UNKNOWN-ACCESS ()) (914 'sebomd' (LOGICAL 4 0 0 LOGICAL ()) () 0 0
0 UNKNOWN-ACCESS ())) PUBLIC () 0 0)
467 'qm2_params_type' 'qm2_params_module' 'qm2_params_type' 10 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((915 'tot_heat_form' (REAL 8 0 0
REAL ()) () 0 0 0 UNKNOWN-ACCESS ()) (916 'sp_quantum_number' (INTEGER 4
0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (917
'd_quantum_number' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (918 'gss' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1
1 0 UNKNOWN-ACCESS ()) (919 'hsp' (REAL 8 0 0 REAL ()) (1 DEFERRED () ())
1 1 0 UNKNOWN-ACCESS ()) (920 'hpp' (REAL 8 0 0 REAL ()) (1 DEFERRED ()
()) 1 1 0 UNKNOWN-ACCESS ()) (921 'dd' (REAL 8 0 0 REAL ()) (2 DEFERRED
() () () ()) 1 1 0 UNKNOWN-ACCESS ()) (922 'po' (REAL 8 0 0 REAL ()) (2
DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (923 'core_chg' (REAL 8 0
0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (924 'orb_elec_ke'
(REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (
925 'betasas' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0
UNKNOWN-ACCESS ()) (926 'betasap' (REAL 8 0 0 REAL ()) (2 DEFERRED () ()
() ()) 1 1 0 UNKNOWN-ACCESS ()) (927 'betasad' (REAL 8 0 0 REAL ()) (2
DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (928 'betapap' (REAL 8 0
0 REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (929
'betapad' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0
UNKNOWN-ACCESS ()) (930 'betadad' (REAL 8 0 0 REAL ()) (2 DEFERRED () ()
() ()) 1 1 0 UNKNOWN-ACCESS ()) (931 'gnn' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (932 'rho_core' (REAL 8 0 0
REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (933 'f0sd' (REAL 8
0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (934 'g2sd' (
REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (935 'fn1'
(REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (
936 'fn2' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0
UNKNOWN-ACCESS ()) (937 'fn3' (REAL 8 0 0 REAL ()) (2 DEFERRED () () ()
()) 1 1 0 UNKNOWN-ACCESS ()) (938 'onec2elec_params' (REAL 8 0 0 REAL ())
(2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (939
'multip_2c_elec_params' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1
1 0 UNKNOWN-ACCESS ()) (940 'cc_exp_params' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (941 'pm6_alpab' (REAL 8 0 0
REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (942 'pm6_xab'
(REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (
943 'pm3mais_alpab' (REAL 8 0 0 REAL ()) (3 DEFERRED () () () () () ())
1 1 0 UNKNOWN-ACCESS ()) (944 'pm3mais_betab' (REAL 8 0 0 REAL ()) (3
DEFERRED () () () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (945 'pm3mais_gamab'
(REAL 8 0 0 REAL ()) (3 DEFERRED () () () () () ()) 1 1 0 UNKNOWN-ACCESS
()) (946 's_orb_exp_by_type' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1
0 UNKNOWN-ACCESS ()) (947 'p_orb_exp_by_type' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (948 'd_orb_exp_by_type' (REAL
8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (949
's_orb_exp_tail_by_type' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (950 'p_orb_exp_tail_by_type' (REAL 8 0 0 REAL ()) (
1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (951 'd_orb_exp_tail_by_type'
(REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (952
'pddge1' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ())
(953 'pddge2' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (954 'scale_factor1_pm3mmx' (REAL 8 0 0 REAL ()) (2
DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (955 'scale_factor2_pm3mmx'
(REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (
956 'rho_pm3mmx' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (957 'atom_orb_zz_sxs_over_sas' (REAL 8 0 0 REAL ())
(4 DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (958
'atom_orb_zz_sxp_over_sap' (REAL 8 0 0 REAL ()) (4 DEFERRED () () () ()
() () () ()) 1 1 0 UNKNOWN-ACCESS ()) (959 'atom_orb_zz_sxd_over_sad' (
REAL 8 0 0 REAL ()) (4 DEFERRED () () () () () () () ()) 1 1 0
UNKNOWN-ACCESS ()) (960 'atom_orb_zz_pxp_over_pap' (REAL 8 0 0 REAL ())
(4 DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (961
'atom_orb_zz_pxd_over_pad' (REAL 8 0 0 REAL ()) (4 DEFERRED () () () ()
() () () ()) 1 1 0 UNKNOWN-ACCESS ()) (962 'atom_orb_zz_dxd_over_dad' (
REAL 8 0 0 REAL ()) (4 DEFERRED () () () () () () () ()) 1 1 0
UNKNOWN-ACCESS ()) (963 'atom_orb_ss_eqn' (REAL 8 0 0 REAL ()) (4
DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (964
'atom_orb_sp_ovlp' (REAL 8 0 0 REAL ()) (4 DEFERRED () () () () () () ()
()) 1 1 0 UNKNOWN-ACCESS ()) (965 'atom_orb_sd_ovlp' (REAL 8 0 0 REAL ())
(4 DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (966
'atom_orb_pd_ovlp' (REAL 8 0 0 REAL ()) (4 DEFERRED () () () () () () ()
()) 1 1 0 UNKNOWN-ACCESS ()) (967 'atom_orb_pp_ovlp_inj' (REAL 8 0 0
REAL ()) (4 DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (
968 'atom_orb_pp_ovlp_ieqj1' (REAL 8 0 0 REAL ()) (4 DEFERRED () () () ()
() () () ()) 1 1 0 UNKNOWN-ACCESS ()) (969 'atom_orb_pp_ovlp_ieqj2' (
REAL 8 0 0 REAL ()) (4 DEFERRED () () () () () () () ()) 1 1 0
UNKNOWN-ACCESS ()) (970 'atom_orb_dd_ovlp_inj' (REAL 8 0 0 REAL ()) (4
DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (971
'atom_orb_dd_ovlp_ieqj1' (REAL 8 0 0 REAL ()) (4 DEFERRED () () () () ()
() () ()) 1 1 0 UNKNOWN-ACCESS ()) (972 'atom_orb_dd_ovlp_ieqj2' (REAL 8
0 0 REAL ()) (4 DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS ())
(973 'atom_orb_ss_eqn_adb' (REAL 8 0 0 REAL ()) (4 DEFERRED () () () ()
() () () ()) 1 1 0 UNKNOWN-ACCESS ()) (974 'atom_orb_sp_eqn_xy' (REAL 8
0 0 REAL ()) (4 DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS ())
(975 'atom_orb_sp_eqn_xx1' (REAL 8 0 0 REAL ()) (4 DEFERRED () () () ()
() () () ()) 1 1 0 UNKNOWN-ACCESS ()) (976 'atom_orb_sp_eqn_xx2' (REAL 8
0 0 REAL ()) (4 DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS ())
(977 'atom_orb_pp_eqn_xxy1' (REAL 8 0 0 REAL ()) (4 DEFERRED () () () ()
() () () ()) 1 1 0 UNKNOWN-ACCESS ()) (978 'atom_orb_pp_eqn_xxy2' (REAL
8 0 0 REAL ()) (4 DEFERRED () () () () () () () ()) 1 1 0 UNKNOWN-ACCESS
()) (979 'pddg_term1' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1
0 UNKNOWN-ACCESS ()) (980 'pddg_term2' (REAL 8 0 0 REAL ()) (2 DEFERRED
() () () ()) 1 1 0 UNKNOWN-ACCESS ()) (981 'pddg_term3' (REAL 8 0 0 REAL
()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (982 'pddg_term4'
(REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (
983 'natomic_orbs' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (984 'orb_loc' (INTEGER 4 0 0 INTEGER ()) (2 DEFERRED
() () () ()) 1 1 0 UNKNOWN-ACCESS ()) (985 'pascal_tri1' (INTEGER 4 0 0
INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (986 'pascal_tri2'
(INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (
987 'num_fn' (INTEGER 4 0 0 INTEGER ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (988 'qxd_supported' (LOGICAL 4 0 0 LOGICAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (989 'qxd_s' (REAL 8 0 0 REAL ())
(1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (990 'qxd_z0' (REAL 8 0 0
REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (991 'qxd_zq' (
REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (992
'qxd_d0' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ())
(993 'qxd_dq' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (994 'qxd_q0' (REAL 8 0 0 REAL ()) (1 DEFERRED () ())
1 1 0 UNKNOWN-ACCESS ()) (995 'qxd_qq' (REAL 8 0 0 REAL ()) (1 DEFERRED
() ()) 1 1 0 UNKNOWN-ACCESS ()) (996 'qxd_neff' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ())) PUBLIC () 0 0)
)

('allocate_qmmm' 0 2 'alph_mm' 0 7 'axis_tol' 0 8 'deallocate_qmmm' 0 9
'default_qmmm_input_options' 0 15 'exponential_cutoff' 0 18
'get_atomic_number' 0 19 'overlap_cutoff' 0 25 'qm2_params' 0 26
'qm2_rij_eqns' 0 28 'qm2_rij_eqns_structure' 0 29 'qm2_struct' 0 32
'qm2_structure' 0 33 'qm_ewald_structure' 0 64 'qm_gb' 0 80
'qm_gb_structure' 0 81 'qmewald' 0 94 'qmmm_div' 0 95 'qmmm_div_structure'
0 96 'qmmm_input_options' 0 99 'qmmm_mpi' 0 166 'qmmm_mpi_structure' 0
167 'qmmm_nml' 0 186 'qmmm_opnq' 0 188 'qmmm_opnq_structure' 0 189
'qmmm_scratch' 0 202 'qmmm_scratch_structure' 0 203 'qmmm_struct' 0 220
'qmmm_vsolv' 0 222 'qmsort' 0 224 'validate_qm_atoms' 0 227)
