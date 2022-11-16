GFORTRAN module created from qm2_read_nm_and_alloc.F90 on Thu Feb 11 10:40:20 2016
MD5:b0afb95a0c04faa8f21784729788c620 -- If you edit this, you'll get what you deserve.

(
() () () () () () () () () () () () () () () () () () () () () () () ()
() () ())

()

()

()

()

(2 'read_qmmm_nm_and_alloc' 'qmmm_read_and_alloc' 'read_qmmm_nm_and_alloc'
1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL UNKNOWN SUBROUTINE
ALWAYS_EXPLICIT) (UNKNOWN 0 0 0 UNKNOWN ()) 3 0 (4 5 6 7 8 9 10 11 12 13
14 15) () 0 () () 0 0)
4 'igb' '' 'igb' 3 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN
UNKNOWN DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () 0 0)
5 'ih' '' 'ih' 3 ((VARIABLE UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN
DIMENSION DUMMY) (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '4'))) 0 0 () (1 ASSUMED_SIZE (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '1') ()) 0 () () 0 0)
6 'ix' '' 'ix' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN DIMENSION
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 ASSUMED_SIZE (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () 0 0)
7 'x' '' 'x' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN DIMENSION
DUMMY) (REAL 8 0 0 REAL ()) 0 0 () (1 ASSUMED_SIZE (CONSTANT (INTEGER 4
0 0 INTEGER ()) 0 '1') ()) 0 () () 0 0)
8 'cut' '' 'cut' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN DUMMY) (
REAL 8 0 0 REAL ()) 0 0 () () 0 () () 0 0)
9 'use_pme' '' 'use_pme' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () 0 0)
10 'ntb' '' 'ntb' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN DUMMY) (
INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () 0 0)
11 'qmstep' '' 'qmstep' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () 0 0)
12 'isabfqm' '' 'isabfqm' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
DIMENSION DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () (1 ASSUMED_SIZE (
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') ()) 0 () () 0 0)
13 'abfqmcharge' '' 'abfqmcharge' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () 0 0)
14 'read_file' '' 'read_file' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () () 0 0)
15 'options' '' 'options' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
OPTIONAL DUMMY) (DERIVED 16 0 0 DERIVED ()) 0 0 () () 0 () () 0 0)
16 'qmmm_input_options' 'qmmm_module' 'qmmm_input_options' 3 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN SEQUENCE) (UNKNOWN 0 0 0
UNKNOWN ()) 0 0 () () 0 ((17 'qmcut' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (18 'lnk_dis' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (19 'scfconv' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (20 'errconv' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (21 'dftb_telec' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (22 'dftb_telec_step' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (23 'fockp_d1' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (24 'fockp_d2' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (25 'fockp_d3' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (26 'fockp_d4' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (27 'damp' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (28 'vshift' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (29 'kappa' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (30 'pseudo_diag_criteria' (REAL 8 0 0 REAL ()) () 0
0 0 UNKNOWN-ACCESS ()) (31 'min_heavy_mass' (REAL 8 0 0 REAL ()) () 0 0
0 UNKNOWN-ACCESS ()) (32 'r_switch_hi' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (33 'r_switch_lo' (REAL 8 0 0 REAL ()) () 0 0 0
UNKNOWN-ACCESS ()) (34 'iqmatoms' (INTEGER 4 0 0 INTEGER ()) (1 EXPLICIT
(CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '10000')) 1 0 0 UNKNOWN-ACCESS ()) (35 'qmgb' (INTEGER 4 0
0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (36 'lnk_atomic_no' (INTEGER 4
0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (37 'ndiis_matrices' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (38 'ndiis_attempts'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (39 'lnk_method'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (40 'qmcharge' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (41 'corecharge' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (42 'buffercharge'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (43 'spin' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (44 'qmqmdx' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (45 'verbosity' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (46 'printcharges'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (47 'printdipole'
(INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (48
'print_eigenvalues' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(49 'peptide_corr' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(50 'itrmax' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (51
'printbondorders' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(52 'qmshake' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
53 'qmmmrij_incore' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ())
(54 'qmqm_erep_incore' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (55 'pseudo_diag' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (56 'qm_ewald' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (57 'qm_pme' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (58 'kmaxqx' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (59 'kmaxqy' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (60 'kmaxqz' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (61 'ksqmaxq' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (62 'qmmm_int' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (63 'adjust_q' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (64 'tight_p_conv' (INTEGER 4 0 0 INTEGER ()) () 0 0
0 UNKNOWN-ACCESS ()) (65 'diag_routine' (INTEGER 4 0 0 INTEGER ()) () 0
0 0 UNKNOWN-ACCESS ()) (66 'density_predict' (INTEGER 4 0 0 INTEGER ())
() 0 0 0 UNKNOWN-ACCESS ()) (67 'fock_predict' (INTEGER 4 0 0 INTEGER ())
() 0 0 0 UNKNOWN-ACCESS ()) (68 'vsolv' (INTEGER 4 0 0 INTEGER ()) () 0
0 0 UNKNOWN-ACCESS ()) (69 'dftb_maxiter' (INTEGER 4 0 0 INTEGER ()) ()
0 0 0 UNKNOWN-ACCESS ()) (70 'dftb_disper' (INTEGER 4 0 0 INTEGER ()) ()
0 0 0 UNKNOWN-ACCESS ()) (71 'dftb_chg' (INTEGER 4 0 0 INTEGER ()) () 0
0 0 UNKNOWN-ACCESS ()) (72 'abfqmmm' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (73 'hot_spot' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (74 'qmmm_switch' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (75 'core_iqmatoms' (INTEGER 4 0 0 INTEGER ()) (1
EXPLICIT (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER
4 0 0 INTEGER ()) 0 '10000')) 1 0 0 UNKNOWN-ACCESS ()) (76
'buffer_iqmatoms' (INTEGER 4 0 0 INTEGER ()) (1 EXPLICIT (CONSTANT (
INTEGER 4 0 0 INTEGER ()) 0 '1') (CONSTANT (INTEGER 4 0 0 INTEGER ()) 0
'10000')) 1 0 0 UNKNOWN-ACCESS ()) (77 'qmmask' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () 0 0 0
UNKNOWN-ACCESS ()) (78 'coremask' (CHARACTER 1 0 0 CHARACTER ((CONSTANT
(INTEGER 4 0 0 INTEGER ()) 0 '8192'))) () 0 0 0 UNKNOWN-ACCESS ()) (79
'buffermask' (CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0
INTEGER ()) 0 '8192'))) () 0 0 0 UNKNOWN-ACCESS ()) (80 'centermask' (
CHARACTER 1 0 0 CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '8192')))
() 0 0 0 UNKNOWN-ACCESS ()) (81 'dftb_3rd_order' (CHARACTER 1 0 0
CHARACTER ((CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '256'))) () 0 0 0
UNKNOWN-ACCESS ()) (82 'qm_theory' (CHARACTER 1 0 0 CHARACTER ((
CONSTANT (INTEGER 4 0 0 INTEGER ()) 0 '12'))) () 0 0 0 UNKNOWN-ACCESS ()))
PUBLIC () 0 0)
)

('read_qmmm_nm_and_alloc' 0 2)
