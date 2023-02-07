GFORTRAN module created from qm2_diagonalizer_module.F90 on Thu Feb 11 10:37:33 2016
MD5:f1b53d167e1de69e5327353acf95e421 -- If you edit this, you'll get what you deserve.

(() () () () () () () ()
() () () () () () () () () () () () () () () () () () ())

()

()

()

()

(2 'qm2_diagonalizer_setup' 'qm2_diagonalizer_module'
'qm2_diagonalizer_setup' 1 ((PROCEDURE UNKNOWN-INTENT MODULE-PROC DECL
UNKNOWN SUBROUTINE) (UNKNOWN 0 0 0 UNKNOWN ()) 3 0 (4 5 6 7 8 9 10) () 0
() () 0 0)
4 'diag_routine' '' 'diag_routine' 3 ((VARIABLE INOUT UNKNOWN-PROC
UNKNOWN UNKNOWN DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () 0 0)
5 'allow_pseudo_diag' '' 'allow_pseudo_diag' 3 ((VARIABLE IN
UNKNOWN-PROC UNKNOWN UNKNOWN DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () ()
0 () () 0 0)
6 'verbosity' '' 'verbosity' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN
UNKNOWN DUMMY) (INTEGER 4 0 0 INTEGER ()) 0 0 () () 0 () () 0 0)
7 'master' '' 'master' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () () 0 0)
8 'qm2_struct' '' 'qm2_struct' 3 ((VARIABLE INOUT UNKNOWN-PROC UNKNOWN
UNKNOWN DUMMY) (DERIVED 11 0 0 DERIVED ()) 0 0 () () 0 () () 0 0)
9 'qmmm_scratch' '' 'qmmm_scratch' 3 ((VARIABLE INOUT UNKNOWN-PROC
UNKNOWN UNKNOWN DUMMY) (DERIVED 12 0 0 DERIVED ()) 0 0 () () 0 () () 0 0)
10 'silence' '' 'silence' 3 ((VARIABLE IN UNKNOWN-PROC UNKNOWN UNKNOWN
DUMMY) (LOGICAL 4 0 0 LOGICAL ()) 0 0 () () 0 () () 0 0)
11 'qm2_structure' 'qmmm_module' 'qm2_structure' 3 ((DERIVED
UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (UNKNOWN 0 0 0
UNKNOWN ()) 0 0 () () 0 ((13 'den_matrix' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(14 'old_den_matrix' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (15 'old2_density' (
REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0)) (16 'md_den_mat_guess1' (REAL 8 0 0 REAL ())
(1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ())
0)) (17 'md_den_mat_guess2' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1
0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (18
'fock_mat_final4' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (19 'fock_mat_final3'
(REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0)) (20 'fock_mat_final2' (REAL 8 0 0 REAL ())
(1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ())
0)) (21 'fock_mat_final1' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (22 'fock_matrix' (
REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0)) (23 'qm_mm_e_repul' (REAL 8 0 0 REAL ()) (
2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (24 'qm_qm_2e_repul' (REAL 8 0 0 REAL ()) (1 DEFERRED ()
()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (25
'hmatrix' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (
NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (26 'qm_qm_e_repul' (REAL 8 0 0 REAL
()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0
UNKNOWN ()) 0)) (27 'fock2_ptot2' (REAL 8 0 0 REAL ()) (2 DEFERRED () ()
() ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (28
'eigen_vectors' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (29 'eigen_values' (
REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (
UNKNOWN 0 0 0 UNKNOWN ()) 0)) (30 'scf_mchg' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0))
(31 'diis_fock' (REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0
UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN ()) 0)) (32 'diis_errmat' (
REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS (NULL
(UNKNOWN 0 0 0 UNKNOWN ()) 0)) (33 'diis_mat' (REAL 8 0 0 REAL ()) (2
DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS (NULL (UNKNOWN 0 0 0 UNKNOWN
()) 0)) (34 'matsize' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS
()) (35 'n2el' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
36 'norbs' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (37
'nclosed' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (38
'nopenclosed' (INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (
39 'qm_mm_e_repul_allocated' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (40 'n_peptide_links' (INTEGER 4 0 0 INTEGER ()) () 0
0 0 UNKNOWN-ACCESS ()) (41 'peptide_links' (INTEGER 4 0 0 INTEGER ()) (
2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (42 'calc_mchg_scf' (
LOGICAL 4 0 0 LOGICAL ()) () 0 0 0 UNKNOWN-ACCESS ())) PUBLIC () 0 0)
12 'qmmm_scratch_structure' 'qmmm_module' 'qmmm_scratch_structure' 3 ((
DERIVED UNKNOWN-INTENT UNKNOWN-PROC UNKNOWN UNKNOWN POINTER_COMP) (
UNKNOWN 0 0 0 UNKNOWN ()) 0 0 () () 0 ((43 'matsize_red_scratch' (REAL 8
0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (44
'qm_pme_scratch' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (45 'mat_diag_workspace' (REAL 8 0 0 REAL ()) (2
DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (46 'pdiag_scr_norbs_norbs'
(REAL 8 0 0 REAL ()) (2 DEFERRED () () () ()) 1 1 0 UNKNOWN-ACCESS ()) (
47 'pdiag_scr_noccupied_norbs' (REAL 8 0 0 REAL ()) (2 DEFERRED () () ()
()) 1 1 0 UNKNOWN-ACCESS ()) (48 'pdiag_vectmp1' (REAL 8 0 0 REAL ()) (
1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (49 'pdiag_vectmp2' (REAL 8 0
0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (50 'pdiag_vectmp3'
(REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (51
'pdiag_vecjs' (INTEGER 4 0 0 INTEGER ()) (2 DEFERRED () () () ()) 1 1 0
UNKNOWN-ACCESS ()) (52 'lapack_dc_real_scr' (REAL 8 0 0 REAL ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (53 'lapack_dc_int_scr' (REAL 8
0 0 REAL ()) (1 DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (54
'qm_real_scratch' (REAL 8 0 0 REAL ()) (1 DEFERRED () ()) 1 1 0
UNKNOWN-ACCESS ()) (55 'qm_int_scratch' (INTEGER 4 0 0 INTEGER ()) (1
DEFERRED () ()) 1 1 0 UNKNOWN-ACCESS ()) (56 'lapack_dc_real_scr_aloc' (
INTEGER 4 0 0 INTEGER ()) () 0 0 0 UNKNOWN-ACCESS ()) (57
'lapack_dc_int_scr_aloc' (INTEGER 4 0 0 INTEGER ()) () 0 0 0
UNKNOWN-ACCESS ()) (58 'qm_mm_pairs_allocated' (INTEGER 4 0 0 INTEGER ())
() 0 0 0 UNKNOWN-ACCESS ())) PUBLIC () 0 0)
)

('qm2_diagonalizer_setup' 0 2)
