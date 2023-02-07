dnl Generates test code for hemv2
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"
#include "blas_extended_test.h"

/* 0 -- 1 */
#define UPLO_START 0
#define UPLO_END   1

/* 0 -- 1 */
#define ORDER_START  0
#define ORDER_END    1

/* 0 -- 2 */
#define ALPHA_START  0
#define ALPHA_END    2

/* 0 -- 2 */
#define BETA_START   0
#define BETA_END     2

/* -1 -- 1 */
#define NORM_START   -1
#define NORM_END     1

/* 0 -- 2 */
#define LDA_START    0
#define LDA_END      2

/* 0 -- 2 */
#define PREC_START   0
#define PREC_END     2

/* 0 -- 1 */
#define RANDOMIZE_START 0
#define RANDOMIZE_END   1

/* -2 -- 2 (Stride) */
#define INCX_START -2
#define INCX_END 2

/* -2 -- 2 (Stride) */
#define INCY_START -2
#define INCY_END 2

#define NUM_DATA 7

include(cblas.m4)dnl
include(test-common.m4)dnl
include(hemv2/hemv2-common.m4)dnl
dnl
dnl
define(`DO_TEST_HEMV2_NAME', 
  `ifelse(`$1&&$1', `$2&&$3', `do_test_$1hemv2$4', `do_test_$1hemv2_$2_$3$4')')
dnl
dnl
define(`DO_TEST_HEMV2_PARAMS', 
  `int n,  
   int ntests, int *seed, double thresh, int debug, float test_prob, 
   double *min_ratio, double *max_ratio, int *num_bad_ratio, int *num_tests')
dnl
dnl
define(`DO_TEST_HEMV2_HEAD', 
  `void DO_TEST_HEMV2_NAME($1, $2, $3, $4) 
           (DO_TEST_HEMV2_PARAMS($1, $2, $3, $4))')
dnl
dnl
define(`DO_TEST_HEMV2', 
  `DO_TEST_HEMV2_HEAD($1, $2, $3, $4) {
   DO_TEST_HEMV2_BODY($1, $2, $3, $4)
  }')
dnl
dnl
dnl
define(`TESTGEN_HEMV2_NAME', `ifelse(
       `$2&&$3', `$1&&$1', `BLAS_$1hemv2_testgen', `BLAS_$1hemv2_$2_$3_testgen')')
dnl
dnl
define(`DO_TEST_HEMV2_BODY', 
  `
  /* Function name */
  const char fname[] = "DO_TEST_HEMV2_NAME($1, $2, $3, $4)";
  int i;
  int yi;
  int incyi, yi0;
  int test_count;
  int bad_ratio_count;
  int ri;
  int incri = 1;
  int incx, incy;
  double ratio;
  double ratio_min, ratio_max;
  double eps_int;  /* internal machine epsilon     */
  double un_int;   /* internal underflow threshold */

  DECLARE(rin, $1_type)
  DECLARE(rout, $1_type)
  DECLARE(r_true_elem, EXTRA_TYPE($1_type))

  enum blas_order_type order_type;
  enum blas_uplo_type uplo_type;
  enum blas_prec_type prec;
 
  int order_val, uplo_val;
  int lda_val, incx_val, incy_val;
  int alpha_val, beta_val;
  
  ifelse(`$4', `_x', `int prec_val;', `')
 
  int lda;
  int alpha_flag, beta_flag;
  int saved_seed;
  int norm;
  int test_no;

  DECLARE(alpha, $1_type)
  DECLARE(beta, $1_type)
  DECLARE_VECTOR(a, $2_type)
  DECLARE_VECTOR(head_x, $3_type)
  DECLARE_VECTOR(tail_x, $3_type)
  DECLARE_VECTOR(y, $1_type)
  DECLARE_VECTOR(a_vec, $2_type)
  DECLARE_VECTOR(y_gen, $1_type)
  DECLARE_VECTOR(head_x_gen, $3_type)
  DECLARE_VECTOR(tail_x_gen, $3_type)
  DECLARE_VECTOR(ratios, real_D)

  /* true result calculated by testgen, in double-double */
  DECLARE_VECTOR(r_true, EXTRA_TYPE($1_type))

  FPU_FIX_DECL;
  
  if (n < 0) BLAS_error(fname, -1, n, NULL);
  if (ntests < 0) BLAS_error(fname, -2, ntests, NULL);

  /* initialization */
  saved_seed = *seed;
  ratio = 0.0;
  ratio_min = 1e308;
  ratio_max = 0.0;

  *num_tests = 0;
  *num_bad_ratio = 0;
  *min_ratio = 0.0;
  *max_ratio = 0.0;

  if (n == 0) return;
  INC_ADJUST(incri, EXTRA_TYPE($1_type))

  FPU_FIX_START;
 
  MALLOC_VECTOR(y, $1_type, 2 * n)
  MALLOC_VECTOR(y_gen, $1_type, n)
  MALLOC_VECTOR(head_x_gen, $3_type, n)
  MALLOC_VECTOR(tail_x_gen, $3_type, n)
  MALLOC_VECTOR(a, $2_type, 2 * n * n)
  MALLOC_VECTOR(head_x, $3_type, 2 * n)
  MALLOC_VECTOR(tail_x, $3_type, 2 * n)
  MALLOC_VECTOR(a_vec, $2_type, n)
  MALLOC_VECTOR(r_true, EXTRA_TYPE($1_type), n)
  MALLOC_VECTOR(ratios, real_D, n)

  test_count = 0;
  bad_ratio_count = 0;

  /* vary alpha */
  for (alpha_val = ALPHA_START; alpha_val <= ALPHA_END; alpha_val++) {
    
    SET_ALPHA($1_type)

    /* vary beta */
    for (beta_val = BETA_START; beta_val <= BETA_END; beta_val++) {
      SET_BETA($1_type)

      ifelse($4, _x, `
      /* varying extra precs */
      for (prec_val = PREC_START; prec_val <= PREC_END; prec_val++) {')
      SET_INTERNAL_PARAMS($1_type, $4)

      /* vary norm -- underflow, approx 1, overflow */
      for (norm = NORM_START; norm <= NORM_END; norm++) {
       
        /* number of tests */
        for (test_no = 0; test_no < ntests; test_no++) {

          /* vary storage format */
          for (order_val = ORDER_START; order_val <= ORDER_END; order_val++) {

            order_type = (order_val == 0) ? blas_colmajor : blas_rowmajor;

            /* vary upper / lower variation */
            for (uplo_val = UPLO_START; uplo_val <= UPLO_END; uplo_val++) {

              uplo_type = (uplo_val == 0) ? blas_upper : blas_lower;
      
                /* vary lda = n, n+1, 2*n */
                for (lda_val = LDA_START; lda_val <= LDA_END; lda_val++) {
                
                  lda = (lda_val == 0) ? n : (lda_val == 1) ? n+1 : 2*n;
                    
                  saved_seed = *seed;
                  /* For the sake of speed, we throw out this case at random */
                  if ( xrand(seed) >= test_prob ) continue;

                  SET_ALPHA($1_type)
                  SET_BETA($1_type)

                  /* finally we are here to generate the test case */
                  TESTGEN_HEMV2_NAME($1, $2, $3) (norm, order_type, 
                    uplo_type, n, &alpha, alpha_flag, &beta, beta_flag,
                    a, lda, head_x_gen, tail_x_gen, y_gen, seed, HEAD(r_true),
                    TAIL(r_true));
                  test_count++;

                  /* vary incx = -2, -1, 1, 2 */
                  for (incx_val = INCX_START; incx_val <= INCX_END; incx_val++) {
                    
                    incx = incx_val;
                    if (0 == incx) continue;
                    
                    /* vary incy = -2, -1, 1, 2 */
                    for (incy_val = INCY_START; incy_val <= INCY_END; incy_val++) {
                      
                      incy = incy_val;
                      if (0 == incy) continue;

                      /* copy generated vector with appropriate incs. */
                      $1copy_vector(y_gen, n, 1, y, incy);
                      $3copy_vector(head_x_gen, n, 1, head_x, incx);
                      $3copy_vector(tail_x_gen, n, 1, tail_x, incx);

                      /* call hemv2 routines to be tested */ 
                      FPU_FIX_STOP;
                      HEMV2_NAME($1, $2, $3, $4)(order_type,
                        uplo_type, n, alpha, a, lda, head_x, tail_x, incx, beta, 
                        y, incy ifelse(`$4', `_x', `, prec'));
                      FPU_FIX_START;

                      /* now compute the ratio using test_BLAS_xdot */
                      /* copy a row from A, use x, run dot test */

                      incyi = incy;
                      INC_ADJUST(incyi, $1_type)
                      yi0 = (incy > 0) ? 0 : (-n + 1) * incyi;

                      for (i = 0, yi = yi0, ri = 0; 
                           i < n; i++, yi += incyi, ri += incri) {
                        $2he_copy_row(order_type, uplo_type, blas_left_side,
                                      n, a, lda, a_vec, i);   
                          
                        /* just use the x vector - it was unchanged (in theory) */
                        GET_VECTOR_ELEMENT(rin, y_gen, i, $1_type)
                        GET_VECTOR_ELEMENT(rout, y, yi, $1_type)
                        GET_VECTOR_ELEMENT(r_true_elem, r_true, ri, EXTRA_TYPE($1_type))

                        TEST_DOT2_NAME($1, $2, $3, $4)(n, blas_no_conj, alpha, beta,
                          rin, rout, HEAD(r_true)_elem, TAIL(r_true)_elem, a_vec, 1,
                          head_x, tail_x, incx, eps_int, un_int, &ratios[i]);

                        /* take the max ratio */
                        if (i == 0) {
                          ratio = ratios[0];

                          /* The !<= below causes NaN errors to be included.
                           * Note that (NaN > 0) is false */
                        } else if (!(ratios[i] <= ratio)) {
                          ratio = ratios[i];
                        }
                        
                      }  /* end of dot-test loop */


                      /* The !<= below causes NaN errors to be included.
                       * Note that (NaN > 0) is false */
                      if (!(ratio <= thresh)) {

                        if (debug == 3) {
                          printf("\n\t\tTest # %d\n", test_count);
                          printf("y type : $1, a type : $2, x type : $3\n");
                          printf("Seed = %d\t", saved_seed);
                          printf("n %d\n", n);
                          printf("LDA %d  INCX %d  INCY %d\n", lda, incx, incx);

                          if (order_type == blas_rowmajor)
                            printf("row ");
                          else
                            printf("col ");

                          if (uplo_type == blas_upper)
                            printf("upper ");
                          else 
                            printf("lower ");

                          printf("NORM %d, ALPHA %d, BETA %d\n", 
                                 norm, alpha_val, beta_val);

                          /* print out info */
                          PRINT_VAR(alpha, $1_type);
                          printf("   ");
                          PRINT_VAR(beta, $1_type);
                          printf("\n");

                          printf("a\n");
                          $2he_print_matrix(a, n, lda, order_type, uplo_type);
                          $3print_vector(head_x, n, incx, "head_x");
                          $3print_vector(tail_x, n, incx, "tail_x");
                          $1print_vector(y_gen, n, incy, "y_gen");
                          $1print_vector(y, n, incy, "y");
                          DOUBLE_ABBREV($1)print_vector(HEAD(r_true), n, 1, "HEAD(r_true)");
                          dprint_vector(ratios, n, 1, "ratios");
                          printf("ratio = %g\n", ratio);
                        }
                        bad_ratio_count++;
                        if (bad_ratio_count >= MAX_BAD_TESTS) {
                          printf("\ntoo many failures, exiting....");
                          printf("\nTesting and compilation");
                          printf(" are incomplete\n\n");
                          goto end;
                        }
                        if ( !(ratio <= TOTAL_FAILURE_THRESHOLD) ) {
                          printf("\nFlagrant ratio error, exiting...");
                          printf("\nTesting and compilation");
                          printf(" are incomplete\n\n");
                          goto end;
                        }
                      }

                      if (!(ratio <= ratio_max)) ratio_max = ratio;
                      if (ratio != 0.0 && !(ratio >= ratio_min))
                        ratio_min = ratio;

                    } /* end of incy loop */

                  } /* end of incx loop */

                } /* end of lda loop */     
              
            } /* end of uplo loop */

          } /* end of order loop */

        } /* end of nr test loop */
        
      } /* end of norm loop */

dnl ---------------------------------------------------------
ifelse(`$4', `_x', `
        } /* end of prec loop */  ')
dnl ---------------------------------------------------------

    } /* end of beta loop */

  } /* end of alpha loop */

end:
  FPU_FIX_STOP;

  FREE_VECTOR(y, $1_type)
  FREE_VECTOR(a, $2_type)
  FREE_VECTOR(y_gen, $1_type)
  FREE_VECTOR(head_x, $3_type)
  FREE_VECTOR(tail_x, $3_type)
  FREE_VECTOR(head_x_gen, $3_type)
  FREE_VECTOR(tail_x_gen, $3_type)
  FREE_VECTOR(r_true, EXTRA_TYPE($1_type))
  FREE_VECTOR(ratios, real_D)
  FREE_VECTOR(a_vec, $2_type)

  *max_ratio = ratio_max;
  *min_ratio = ratio_min;
  *num_tests = test_count;
  *num_bad_ratio = bad_ratio_count;
') 
dnl
dnl
dnl
define(`CALL_DO_TEST_HEMV2', 
 `fname = "HEMV2_NAME($1, $2, $3, $4)";
  printf("Testing %s...\n", fname);
  total_tests = 0;
  total_bad_ratios = 0;
  total_min_ratio = 1e308;
  total_max_ratio = 0.0;
  for (i = 0; i < nsizes; i++) {
    n = n_data[i][0];

    DO_TEST_HEMV2_NAME($1, $2, $3, $4)(n, ntests, &seed, thresh, debug,
             test_prob,   
             &min_ratio, &max_ratio, &num_bad_ratio, &num_tests);

    if (debug == 2 || (debug == 1 && num_bad_ratio > 0)) {
      printf("   [%d %d]: ", n, n);
      printf("bad/total = %d/%d, min_ratio = %g, max_ratio = %g\n", 
             num_bad_ratio, num_tests, min_ratio, max_ratio);
    }

    total_tests += num_tests;
    total_bad_ratios += num_bad_ratio;
    if (total_min_ratio > min_ratio)
      total_min_ratio = min_ratio;
    if (total_max_ratio < max_ratio)
      total_max_ratio = max_ratio;
  }

  nr_routines++;
  if (total_bad_ratios == 0)
    printf("PASS> ");
  else {
    printf("FAIL> ");
    nr_failed_routines++;
  }
  printf("%-24s: bad/total = %d/%d, max_ratio = %.2e\n\n",
      fname, total_bad_ratios, total_tests, max_ratio);
 ') 

FOREACH(`HEMV2_ARGS', `
DO_TEST_HEMV2(arg)')dnl

MAIN(`
  int i;
  int n_data[NUM_DATA][1] = {{4}, {2}, {3}, {8}, {10}, {1}, {7}};', `
  if (nsizes < 0 || nsizes > NUM_DATA)
    BLAS_error("do_test_hemv2", -1, nsizes, NULL);
FOREACH(`HEMV2_ARGS', `
CALL_DO_TEST_HEMV2(arg)')')dnl

