#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_dsbmv_s_d(enum blas_order_type order, enum blas_uplo_type uplo,
		    int n, int k, double alpha, const float *a, int lda,
		    const double *x, int incx, double beta,
		    double *y, int incy)

/* 
 * Purpose
 * =======
 *
 * This routines computes the matrix product:
 *
 *     y  <-  alpha * A * x  +  beta * y
 * 
 * where A is a symmetric band matrix.
 *
 * Arguments
 * =========
 *
 * order   (input) enum blas_order_type
 *         Storage format of input symmetric matrix A.
 * 
 * uplo    (input) enum blas_uplo_type
 *         Determines which half of matrix A (upper or lower triangle)
 *         is accessed.
 *
 * n       (input) int
 *         Dimension of A and size of vectors x, y.
 *
 * k       (input) int
 *         Number of subdiagonals ( = number of superdiagonals)
 *      
 * alpha   (input) double
 * 
 * a       (input) float*
 *         Matrix A.
 *
 * lda     (input) int
 *         Leading dimension of matrix A.
 *
 * x       (input) double*
 *         Vector x.
 *   
 * incx    (input) int
 *         Stride for vector x.
 *
 * beta    (input) double
 * 
 * y       (input/output) double*
 *         Vector y.
 *
 * incy    (input) int
 *         Stride for vector y. 
 *
 *
 *  Notes on storing a symmetric band matrix:
 * 
 *      Integers in the below arrays represent values of
 *              type float.
 *
 *    if we have a symettric matrix:
 *
 *      1  2  3  0  0
 *      2  4  5  6  0
 *      3  5  7  8  9
 *      0  6  8  10 11
 *      0  0  9  11 12
 *
 *     This matrix has n == 5, and k == 2. It can be stored in the
 *      following ways:
 *
 *      Notes for the examples:
 *      Each column below represents a contigous vector.
 *      Columns are strided by lda.
 *      An asterisk (*) represents a position in the 
 *       matrix that is not used.
 *      Note that the minimum lda (size of column) is 3 (k+1).
 *       lda may be arbitrarily large; an lda > 3 would mean
 *       there would be unused data at the bottom of the below
 *       columns.        
 *
 *    blas_colmajor and blas_upper:
 *      *  *  3  6  9  
 *      *  2  5  8  11 
 *      1  4  7  10 12
 *
 *
 *    blas_colmajor and blas_lower
 *      1  4  7  10  12
 *      2  5  8  11  *
 *      3  6  9  *   *
 *
 *
 *    blas_rowmajor and blas_upper 
 *      Columns here also represent contigous arrays.
 *      1  4  7  10  12
 *      2  5  8  11  *
 *      3  6  9  *   *
 *
 *
 *    blas_rowmajor and blas_lower
 *      Columns here also represent contigous arrays.
 *      *  *  3  6   9
 *      *  2  5  8   11
 *      1  4  7  10  12
 *
 */
{
  static const char routine_name[] = "BLAS_dsbmv_s_d";

  /* Integer Index Variables */
  int i, j;
  int xi, yi;
  int aij, astarti, x_starti, y_starti;
  int incaij, incaij2;
  int n_i;
  int maxj_first, maxj_second;

  /* Input Matrices */
  const float *a_i = a;
  const double *x_i = x;

  /* Output Vector */
  double *y_i = y;

  /* Input Scalars */
  double alpha_i = alpha;
  double beta_i = beta;

  /* Temporary Floating-Point Variables */
  float a_elem;
  double x_elem;
  double y_elem;
  double prod;
  double sum;
  double tmp1;
  double tmp2;


  /* Test for no-op */
  if (n <= 0) {
    return;
  }
  if (alpha_i == 0.0 && beta_i == 1.0) {
    return;
  }

  /* Check for error conditions. */
  if (order != blas_colmajor && order != blas_rowmajor) {
    BLAS_error(routine_name, -1, order, 0);
  }
  if (uplo != blas_upper && uplo != blas_lower) {
    BLAS_error(routine_name, -2, uplo, 0);
  }
  if (n < 0) {
    BLAS_error(routine_name, -3, n, 0);
  }
  if (k < 0 || k > n) {
    BLAS_error(routine_name, -4, k, 0);
  }
  if ((lda < k + 1) || (lda < 1)) {
    BLAS_error(routine_name, -7, lda, 0);
  }
  if (incx == 0) {
    BLAS_error(routine_name, -9, incx, 0);
  }
  if (incy == 0) {
    BLAS_error(routine_name, -12, incy, 0);
  }

  /* Set Index Parameters */
  n_i = n;

  if (((uplo == blas_upper) && (order == blas_colmajor)) ||
      ((uplo == blas_lower) && (order == blas_rowmajor))) {
    incaij = 1;			/* increment in first loop */
    incaij2 = lda - 1;		/* increment in second loop */
    astarti = k;		/* does not start on zero element */
  } else {
    incaij = lda - 1;
    incaij2 = 1;
    astarti = 0;		/* start on first element of array */
  }
  /* Adjustment to increments (if any) */
  if (incx < 0) {
    x_starti = (-n_i + 1) * incx;
  } else {
    x_starti = 0;
  }
  if (incy < 0) {
    y_starti = (-n_i + 1) * incy;
  } else {
    y_starti = 0;
  }



  /* alpha = 0.  In this case, just return beta * y */
  if (alpha_i == 0.0) {
    for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
      y_elem = y_i[yi];
      tmp1 = y_elem * beta_i;
      y_i[yi] = tmp1;
    }
  } else {
    /*  determine the loop interation counts */
    /* number of elements done in first loop 
       (this will increase by one over each column up to a limit) */
    maxj_first = 0;
    /* number of elements done in second loop the first time */
    maxj_second = MIN(k + 1, n_i);

    /* Case alpha == 1. */
    if (alpha_i == 1.0) {

      if (beta_i == 0.0) {
	/* Case alpha = 1, beta = 0.  We compute  y <--- A * x */
	for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	  sum = 0.0;
	  for (j = 0, aij = astarti, xi = x_starti;
	       j < maxj_first; j++, aij += incaij, xi += incx) {
	    a_elem = a_i[aij];
	    x_elem = x_i[xi];
	    prod = a_elem * x_elem;
	    sum = sum + prod;
	  }
	  for (j = 0; j < maxj_second; j++, aij += incaij2, xi += incx) {
	    a_elem = a_i[aij];
	    x_elem = x_i[xi];
	    prod = a_elem * x_elem;
	    sum = sum + prod;
	  }
	  y_i[yi] = sum;
	  if (i + 1 >= (n_i - k))
	    maxj_second--;
	  if (i >= k) {
	    astarti += (incaij + incaij2);
	    x_starti += incx;
	  } else {
	    maxj_first++;
	    astarti += incaij2;
	  }
	}
      } else {
	/* Case alpha = 1, but beta != 0. 
	   We compute  y  <--- A * x + beta * y */
	for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	  sum = 0.0;

	  for (j = 0, aij = astarti, xi = x_starti;
	       j < maxj_first; j++, aij += incaij, xi += incx) {
	    a_elem = a_i[aij];
	    x_elem = x_i[xi];
	    prod = a_elem * x_elem;
	    sum = sum + prod;
	  }
	  for (j = 0; j < maxj_second; j++, aij += incaij2, xi += incx) {
	    a_elem = a_i[aij];
	    x_elem = x_i[xi];
	    prod = a_elem * x_elem;
	    sum = sum + prod;
	  }
	  y_elem = y_i[yi];
	  tmp2 = y_elem * beta_i;
	  tmp1 = sum;
	  tmp1 = tmp2 + tmp1;
	  y_i[yi] = tmp1;
	  if (i + 1 >= (n_i - k))
	    maxj_second--;
	  if (i >= k) {
	    astarti += (incaij + incaij2);
	    x_starti += incx;
	  } else {
	    maxj_first++;
	    astarti += incaij2;
	  }
	}
      }
    } else {
      if (beta_i == 0.0) {
	/* Case alpha != 1, but beta == 0. 
	   We compute  y  <--- A * x * a */
	for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	  sum = 0.0;

	  for (j = 0, aij = astarti, xi = x_starti;
	       j < maxj_first; j++, aij += incaij, xi += incx) {
	    a_elem = a_i[aij];
	    x_elem = x_i[xi];
	    prod = a_elem * x_elem;
	    sum = sum + prod;
	  }
	  for (j = 0; j < maxj_second; j++, aij += incaij2, xi += incx) {
	    a_elem = a_i[aij];
	    x_elem = x_i[xi];
	    prod = a_elem * x_elem;
	    sum = sum + prod;
	  }
	  y_elem = y_i[yi];
	  tmp1 = sum * alpha_i;
	  y_i[yi] = tmp1;
	  if (i + 1 >= (n_i - k))
	    maxj_second--;
	  if (i >= k) {
	    astarti += (incaij + incaij2);
	    x_starti += incx;
	  } else {
	    maxj_first++;
	    astarti += incaij2;
	  }
	}
      } else {
	/* The most general form,   y <--- alpha * A * x + beta * y */
	for (i = 0, yi = y_starti; i < n_i; i++, yi += incy) {
	  sum = 0.0;

	  for (j = 0, aij = astarti, xi = x_starti;
	       j < maxj_first; j++, aij += incaij, xi += incx) {
	    a_elem = a_i[aij];
	    x_elem = x_i[xi];
	    prod = a_elem * x_elem;
	    sum = sum + prod;
	  }
	  for (j = 0; j < maxj_second; j++, aij += incaij2, xi += incx) {
	    a_elem = a_i[aij];
	    x_elem = x_i[xi];
	    prod = a_elem * x_elem;
	    sum = sum + prod;
	  }
	  y_elem = y_i[yi];
	  tmp2 = y_elem * beta_i;
	  tmp1 = sum * alpha_i;
	  tmp1 = tmp2 + tmp1;
	  y_i[yi] = tmp1;
	  if (i + 1 >= (n_i - k))
	    maxj_second--;
	  if (i >= k) {
	    astarti += (incaij + incaij2);
	    x_starti += incx;
	  } else {
	    maxj_first++;
	    astarti += incaij2;
	  }
	}
      }
    }
  }


}				/* end BLAS_dsbmv_s_d */
