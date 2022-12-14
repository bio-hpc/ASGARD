#include <math.h>
#include "blas_extended.h"
#include "blas_extended_private.h"

void BLAS_dtrsv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, double alpha, const double *T, int ldt,
		  double *x, int incx, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 * 
 * This routine solve :
 * 
 *     x <- alpha * inverse(T) * x
 * 
 * Arguments
 * =========
 * 
 * order  (input) enum blas_order_type
 *        column major, row major
 *
 * uplo   (input) enum blas_uplo_type
 *        upper, lower
 *
 * trans  (input) enum blas_trans_type
 *        no trans, trans, conj trans
 * 
 * diag   (input) enum blas_diag_type
 *        unit, non unit
 *
 * n      (input) int
 *        the dimension of T
 * 
 * alpha  (input) double
 * 
 * T      (input) double*
 *        Triangular matrix
 *
 * x      (input) const double*
 *           Array of length n.
 * 
 * incx   (input) int
 *           The stride used to access components x[i].
 *
 * prec   (input) enum blas_prec_type
 *        Specifies the internal precision to be used.
 *        = blas_prec_single: single precision.
 *        = blas_prec_double: double precision.
 *        = blas_prec_extra : anything at least 1.5 times as accurate
 *                            than double, and wider than 80-bits.
 *                            We use double-double in our implementation.
 *
 */
{
  char *routine_name = "BLAS_dtrsv";

  int i, j;			/* used to idx matrix */
  int ix, jx;			/* used to idx vector x */
  int start_x;			/* used as the starting idx to vector x */
  const double *T_i = T;	/* internal matrix T */
  double *x_i = x;		/* internal x */
  double alpha_i = alpha;	/* internal alpha */
  double T_element;		/* temporary variable for an element of matrix A */
  int incT = 1;			/* internal ldt */

  if ((order != blas_rowmajor && order != blas_colmajor) ||
      (uplo != blas_upper && uplo != blas_lower) ||
      (trans != blas_trans && trans !=
       blas_no_trans && trans != blas_conj_trans) ||
      (diag != blas_non_unit_diag && diag != blas_unit_diag) ||
      (ldt < n) || (incx == 0)) {
    BLAS_error(routine_name, 0, 0, NULL);
  }

  if (n <= 0)
    return;



  /* configuring the vector starting idx */
  if (incx <= 0) {
    start_x = -(n - 1) * incx;
  } else {
    start_x = 0;
  }

  /* if alpha is zero, then return x as a zero vector */
  if (alpha_i == 0.0) {
    ix = start_x;
    for (i = 0; i < n; i++) {
      x_i[ix] = 0.0;
      ix += incx;
    }
    return;
  }
  switch (prec) {
  case blas_prec_single:
  case blas_prec_double:
  case blas_prec_indigenous:

    {
      double temp1;		/* temporary variable for calculations */
      double temp2;		/* temporary variable for calculations */
      double temp3;		/* temporary variable for calculations */

      if ((order == blas_rowmajor &&
	   trans == blas_no_trans && uplo == blas_upper) ||
	  (order == blas_colmajor &&
	   trans != blas_no_trans && uplo == blas_lower)) {

	jx = start_x + (n - 1) * incx;
	for (j = n - 1; j >= 0; j--) {

	  /* compute Xj = alpha*Xj - SUM Tij(or Tji) * Xi
	     i=j+1 to n-1           */
	  temp3 = x_i[jx];
	  temp1 = temp3 * alpha_i;

	  ix = start_x + (n - 1) * incx;
	  for (i = n - 1; i >= j + 1; i--) {
	    T_element = T_i[i * incT + j * ldt * incT];

	    temp3 = x_i[ix];
	    temp2 = temp3 * T_element;
	    temp1 = temp1 + (-temp2);
	    ix -= incx;
	  }			/* for j<n */

	  /* if the diagonal entry is not equal to one, then divide Xj by 
	     the entry */
	  if (diag == blas_non_unit_diag) {
	    T_element = T_i[j * incT + j * ldt * incT];


	    temp1 = temp1 / T_element;

	  }
	  /* if (diag == blas_non_unit_diag) */
	  x_i[jx] = temp1;

	  jx -= incx;
	}			/* for j>=0 */
      } else if ((order == blas_rowmajor &&
		  trans == blas_no_trans && uplo == blas_lower) ||
		 (order == blas_colmajor &&
		  trans != blas_no_trans && uplo == blas_upper)) {

	jx = start_x;
	for (j = 0; j < n; j++) {

	  /* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
	     i=j+1 to n-1           */
	  temp3 = x_i[jx];
	  /* multiply by alpha */
	  temp1 = temp3 * alpha_i;

	  ix = start_x;
	  for (i = 0; i < j; i++) {
	    T_element = T_i[i * incT + j * ldt * incT];

	    temp3 = x_i[ix];
	    temp2 = temp3 * T_element;
	    temp1 = temp1 + (-temp2);
	    ix += incx;
	  }			/* for i<j */

	  /* if the diagonal entry is not equal to one, then divide Xj by 
	     the entry */
	  if (diag == blas_non_unit_diag) {
	    T_element = T_i[j * incT + j * ldt * incT];


	    temp1 = temp1 / T_element;

	  }
	  /* if (diag == blas_non_unit_diag) */
	  x_i[jx] = temp1;
	  jx += incx;
	}			/* for j<n */
      } else if ((order == blas_rowmajor &&
		  trans != blas_no_trans && uplo == blas_lower) ||
		 (order == blas_colmajor &&
		  trans == blas_no_trans && uplo == blas_upper)) {

	jx = start_x + (n - 1) * incx;
	for (j = n - 1; j >= 0; j--) {

	  /* compute Xj = alpha*Xj - SUM Tij(or Tji) * Xi
	     i=j+1 to n-1           */
	  temp3 = x_i[jx];
	  temp1 = temp3 * alpha_i;

	  ix = start_x + (n - 1) * incx;
	  for (i = n - 1; i >= j + 1; i--) {
	    T_element = T_i[j * incT + i * ldt * incT];

	    temp3 = x_i[ix];
	    temp2 = temp3 * T_element;
	    temp1 = temp1 + (-temp2);
	    ix -= incx;
	  }			/* for j<n */

	  /* if the diagonal entry is not equal to one, then divide Xj by 
	     the entry */
	  if (diag == blas_non_unit_diag) {
	    T_element = T_i[j * incT + j * ldt * incT];


	    temp1 = temp1 / T_element;

	  }
	  /* if (diag == blas_non_unit_diag) */
	  x_i[jx] = temp1;

	  jx -= incx;
	}			/* for j>=0 */
      } else if ((order == blas_rowmajor &&
		  trans != blas_no_trans && uplo == blas_upper) ||
		 (order == blas_colmajor &&
		  trans == blas_no_trans && uplo == blas_lower)) {

	jx = start_x;
	for (j = 0; j < n; j++) {

	  /* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
	     i=j+1 to n-1           */
	  temp3 = x_i[jx];
	  /* multiply by alpha */
	  temp1 = temp3 * alpha_i;

	  ix = start_x;
	  for (i = 0; i < j; i++) {
	    T_element = T_i[j * incT + i * ldt * incT];

	    temp3 = x_i[ix];
	    temp2 = temp3 * T_element;
	    temp1 = temp1 + (-temp2);
	    ix += incx;
	  }			/* for i<j */

	  /* if the diagonal entry is not equal to one, then divide Xj by 
	     the entry */
	  if (diag == blas_non_unit_diag) {
	    T_element = T_i[j * incT + j * ldt * incT];


	    temp1 = temp1 / T_element;

	  }
	  /* if (diag == blas_non_unit_diag) */
	  x_i[jx] = temp1;
	  jx += incx;
	}			/* for j<n */
      }
    }
    break;
  case blas_prec_extra:
    {
      FPU_FIX_DECL;
      FPU_FIX_START;
      {
	{
	  int inc_intx;		/* inc for intx */
	  double head_temp1, tail_temp1;	/* temporary variable for calculations */
	  double head_temp2, tail_temp2;	/* temporary variable for calculations */
	  double head_temp3, tail_temp3;	/* temporary variable for calculations */
	  double *head_intx, *tail_intx;	/* copy of x used for calculations */

	  /* allocate space for intx */
	  head_intx = (double *) blas_malloc(n * sizeof(double));
	  tail_intx = (double *) blas_malloc(n * sizeof(double));
	  if (n > 0 && (head_intx == NULL || tail_intx == NULL)) {
	    BLAS_error("blas_malloc", 0, 0, "malloc failed.\n");
	  }

	  /* since intx is for internal usage, set it to 1 and then adjust
	     it if necessary */
	  inc_intx = 1;


	  /* copy x to intx */
	  ix = start_x;
	  jx = 0;
	  for (i = 0; i < n; i++) {
	    head_temp1 = x_i[ix];
	    tail_temp1 = 0.0;
	    head_intx[jx] = head_temp1;
	    tail_intx[jx] = tail_temp1;
	    ix += incx;
	    jx += inc_intx;
	  }

	  if ((order == blas_rowmajor &&
	       trans == blas_no_trans && uplo == blas_upper) ||
	      (order == blas_colmajor &&
	       trans != blas_no_trans && uplo == blas_lower)) {

	    jx = (n - 1) * inc_intx;
	    for (j = n - 1; j >= 0; j--) {

	      /* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
	         i=j+1 to n-1           */
	      head_temp3 = head_intx[jx];
	      tail_temp3 = tail_intx[jx];
	      /* multiply by alpha */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_temp3 * split;
		a11 = con - head_temp3;
		a11 = con - a11;
		a21 = head_temp3 - a11;
		con = alpha_i * split;
		b1 = con - alpha_i;
		b1 = con - b1;
		b2 = alpha_i - b1;

		c11 = head_temp3 * alpha_i;
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_temp3 * alpha_i;
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_temp1 = t1 + t2;
		tail_temp1 = t2 - (head_temp1 - t1);
	      }

	      ix = (n - 1) * inc_intx;
	      for (i = n - 1; i >= j + 1; i--) {
		T_element = T_i[i * incT + j * ldt * incT];

		head_temp3 = head_intx[ix];
		tail_temp3 = tail_intx[ix];
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_temp3 * split;
		  a11 = con - head_temp3;
		  a11 = con - a11;
		  a21 = head_temp3 - a11;
		  con = T_element * split;
		  b1 = con - T_element;
		  b1 = con - b1;
		  b2 = T_element - b1;

		  c11 = head_temp3 * T_element;
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_temp3 * T_element;
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_temp2 = t1 + t2;
		  tail_temp2 = t2 - (head_temp2 - t1);
		}
		{
		  double head_bt, tail_bt;
		  head_bt = -head_temp2;
		  tail_bt = -tail_temp2;
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_temp1 + head_bt;
		    bv = s1 - head_temp1;
		    s2 = ((head_bt - bv) + (head_temp1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_temp1 + tail_bt;
		    bv = t1 - tail_temp1;
		    t2 = ((tail_bt - bv) + (tail_temp1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_temp1 = t1 + t2;
		    tail_temp1 = t2 - (head_temp1 - t1);
		  }
		}
		ix -= inc_intx;
	      }			/* for j<n */

	      /* if the diagonal entry is not equal to one, then divide Xj by 
	         the entry */
	      if (diag == blas_non_unit_diag) {
		T_element = T_i[j * incT + j * ldt * incT];


		{
		  /* Compute double-double = double-double / double,
		     using a Newton iteration scheme. */
		  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

		  /* Compute a DP approximation to the quotient. */
		  t1 = head_temp1 / T_element;

		  /* Split t1 and b into two parts with at most 26 bits each,
		     using the Dekker-Veltkamp method. */
		  con = t1 * split;
		  t11 = con - (con - t1);
		  t21 = t1 - t11;
		  con = T_element * split;
		  b1 = con - (con - T_element);
		  b2 = T_element - b1;

		  /* Compute t1 * b using Dekker method. */
		  t12 = t1 * T_element;
		  t22 = (((t11 * b1 - t12) + t11 * b2) + t21 * b1) + t21 * b2;

		  /* Compute dda - (t12, t22) using Knuth trick. */
		  t11 = head_temp1 - t12;
		  e = t11 - head_temp1;
		  t21 =
		    ((-t12 - e) + (head_temp1 - (t11 - e))) + tail_temp1 -
		    t22;

		  /* Compute high-order word of (t11, t21) and divide by b. */
		  t2 = (t11 + t21) / T_element;

		  /* The result is t1 + t2, after normalization. */
		  head_temp1 = t1 + t2;
		  tail_temp1 = t2 - (head_temp1 - t1);
		}

	      }
	      /* if (diag == blas_non_unit_diag) */
	      head_intx[jx] = head_temp1;
	      tail_intx[jx] = tail_temp1;

	      jx -= inc_intx;
	    }			/* for j>=0 */
	  } else if ((order == blas_rowmajor &&
		      trans == blas_no_trans && uplo == blas_lower) ||
		     (order == blas_colmajor &&
		      trans != blas_no_trans && uplo == blas_upper)) {

	    jx = 0;
	    for (j = 0; j < n; j++) {

	      /* compute Xj = Xj - SUM Aij(or Aji) * Xi
	         i=j+1 to n-1           */
	      head_temp3 = head_intx[jx];
	      tail_temp3 = tail_intx[jx];
	      /* multiply by alpha */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_temp3 * split;
		a11 = con - head_temp3;
		a11 = con - a11;
		a21 = head_temp3 - a11;
		con = alpha_i * split;
		b1 = con - alpha_i;
		b1 = con - b1;
		b2 = alpha_i - b1;

		c11 = head_temp3 * alpha_i;
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_temp3 * alpha_i;
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_temp1 = t1 + t2;
		tail_temp1 = t2 - (head_temp1 - t1);
	      }

	      ix = 0;
	      for (i = 0; i < j; i++) {
		T_element = T_i[i * incT + j * ldt * incT];

		head_temp3 = head_intx[ix];
		tail_temp3 = tail_intx[ix];
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_temp3 * split;
		  a11 = con - head_temp3;
		  a11 = con - a11;
		  a21 = head_temp3 - a11;
		  con = T_element * split;
		  b1 = con - T_element;
		  b1 = con - b1;
		  b2 = T_element - b1;

		  c11 = head_temp3 * T_element;
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_temp3 * T_element;
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_temp2 = t1 + t2;
		  tail_temp2 = t2 - (head_temp2 - t1);
		}
		{
		  double head_bt, tail_bt;
		  head_bt = -head_temp2;
		  tail_bt = -tail_temp2;
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_temp1 + head_bt;
		    bv = s1 - head_temp1;
		    s2 = ((head_bt - bv) + (head_temp1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_temp1 + tail_bt;
		    bv = t1 - tail_temp1;
		    t2 = ((tail_bt - bv) + (tail_temp1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_temp1 = t1 + t2;
		    tail_temp1 = t2 - (head_temp1 - t1);
		  }
		}
		ix += inc_intx;
	      }			/* for i<j */

	      /* if the diagonal entry is not equal to one, then divide Xj by 
	         the entry */
	      if (diag == blas_non_unit_diag) {
		T_element = T_i[j * incT + j * ldt * incT];


		{
		  /* Compute double-double = double-double / double,
		     using a Newton iteration scheme. */
		  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

		  /* Compute a DP approximation to the quotient. */
		  t1 = head_temp1 / T_element;

		  /* Split t1 and b into two parts with at most 26 bits each,
		     using the Dekker-Veltkamp method. */
		  con = t1 * split;
		  t11 = con - (con - t1);
		  t21 = t1 - t11;
		  con = T_element * split;
		  b1 = con - (con - T_element);
		  b2 = T_element - b1;

		  /* Compute t1 * b using Dekker method. */
		  t12 = t1 * T_element;
		  t22 = (((t11 * b1 - t12) + t11 * b2) + t21 * b1) + t21 * b2;

		  /* Compute dda - (t12, t22) using Knuth trick. */
		  t11 = head_temp1 - t12;
		  e = t11 - head_temp1;
		  t21 =
		    ((-t12 - e) + (head_temp1 - (t11 - e))) + tail_temp1 -
		    t22;

		  /* Compute high-order word of (t11, t21) and divide by b. */
		  t2 = (t11 + t21) / T_element;

		  /* The result is t1 + t2, after normalization. */
		  head_temp1 = t1 + t2;
		  tail_temp1 = t2 - (head_temp1 - t1);
		}

	      }
	      /* if (diag == blas_non_unit_diag) */
	      head_intx[jx] = head_temp1;
	      tail_intx[jx] = tail_temp1;
	      jx += inc_intx;
	    }			/* for j<n */
	  } else if ((order == blas_rowmajor &&
		      trans != blas_no_trans && uplo == blas_lower) ||
		     (order == blas_colmajor &&
		      trans == blas_no_trans && uplo == blas_upper)) {

	    jx = (n - 1) * inc_intx;
	    for (j = n - 1; j >= 0; j--) {

	      /* compute Xj = alpha*Xj - SUM Aij(or Aji) * Xi
	         i=j+1 to n-1           */
	      head_temp3 = head_intx[jx];
	      tail_temp3 = tail_intx[jx];
	      /* multiply by alpha */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_temp3 * split;
		a11 = con - head_temp3;
		a11 = con - a11;
		a21 = head_temp3 - a11;
		con = alpha_i * split;
		b1 = con - alpha_i;
		b1 = con - b1;
		b2 = alpha_i - b1;

		c11 = head_temp3 * alpha_i;
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_temp3 * alpha_i;
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_temp1 = t1 + t2;
		tail_temp1 = t2 - (head_temp1 - t1);
	      }

	      ix = (n - 1) * inc_intx;
	      for (i = n - 1; i >= j + 1; i--) {
		T_element = T_i[j * incT + i * ldt * incT];

		head_temp3 = head_intx[ix];
		tail_temp3 = tail_intx[ix];
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_temp3 * split;
		  a11 = con - head_temp3;
		  a11 = con - a11;
		  a21 = head_temp3 - a11;
		  con = T_element * split;
		  b1 = con - T_element;
		  b1 = con - b1;
		  b2 = T_element - b1;

		  c11 = head_temp3 * T_element;
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_temp3 * T_element;
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_temp2 = t1 + t2;
		  tail_temp2 = t2 - (head_temp2 - t1);
		}
		{
		  double head_bt, tail_bt;
		  head_bt = -head_temp2;
		  tail_bt = -tail_temp2;
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_temp1 + head_bt;
		    bv = s1 - head_temp1;
		    s2 = ((head_bt - bv) + (head_temp1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_temp1 + tail_bt;
		    bv = t1 - tail_temp1;
		    t2 = ((tail_bt - bv) + (tail_temp1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_temp1 = t1 + t2;
		    tail_temp1 = t2 - (head_temp1 - t1);
		  }
		}
		ix -= inc_intx;
	      }			/* for j<n */

	      /* if the diagonal entry is not equal to one, then divide Xj by 
	         the entry */
	      if (diag == blas_non_unit_diag) {
		T_element = T_i[j * incT + j * ldt * incT];


		{
		  /* Compute double-double = double-double / double,
		     using a Newton iteration scheme. */
		  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

		  /* Compute a DP approximation to the quotient. */
		  t1 = head_temp1 / T_element;

		  /* Split t1 and b into two parts with at most 26 bits each,
		     using the Dekker-Veltkamp method. */
		  con = t1 * split;
		  t11 = con - (con - t1);
		  t21 = t1 - t11;
		  con = T_element * split;
		  b1 = con - (con - T_element);
		  b2 = T_element - b1;

		  /* Compute t1 * b using Dekker method. */
		  t12 = t1 * T_element;
		  t22 = (((t11 * b1 - t12) + t11 * b2) + t21 * b1) + t21 * b2;

		  /* Compute dda - (t12, t22) using Knuth trick. */
		  t11 = head_temp1 - t12;
		  e = t11 - head_temp1;
		  t21 =
		    ((-t12 - e) + (head_temp1 - (t11 - e))) + tail_temp1 -
		    t22;

		  /* Compute high-order word of (t11, t21) and divide by b. */
		  t2 = (t11 + t21) / T_element;

		  /* The result is t1 + t2, after normalization. */
		  head_temp1 = t1 + t2;
		  tail_temp1 = t2 - (head_temp1 - t1);
		}

	      }
	      /* if (diag == blas_non_unit_diag) */
	      head_intx[jx] = head_temp1;
	      tail_intx[jx] = tail_temp1;

	      jx -= inc_intx;
	    }			/* for j>=0 */
	  } else if ((order == blas_rowmajor &&
		      trans != blas_no_trans && uplo == blas_upper) ||
		     (order == blas_colmajor &&
		      trans == blas_no_trans && uplo == blas_lower)) {

	    jx = 0;
	    for (j = 0; j < n; j++) {

	      /* compute Xj = Xj - SUM Aij(or Aji) * Xi
	         i=j+1 to n-1           */
	      head_temp3 = head_intx[jx];
	      tail_temp3 = tail_intx[jx];
	      /* multiply by alpha */
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_temp3 * split;
		a11 = con - head_temp3;
		a11 = con - a11;
		a21 = head_temp3 - a11;
		con = alpha_i * split;
		b1 = con - alpha_i;
		b1 = con - b1;
		b2 = alpha_i - b1;

		c11 = head_temp3 * alpha_i;
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_temp3 * alpha_i;
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_temp1 = t1 + t2;
		tail_temp1 = t2 - (head_temp1 - t1);
	      }

	      ix = 0;
	      for (i = 0; i < j; i++) {
		T_element = T_i[j * incT + i * ldt * incT];

		head_temp3 = head_intx[ix];
		tail_temp3 = tail_intx[ix];
		{
		  /* Compute double-double = double-double * double. */
		  double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		  con = head_temp3 * split;
		  a11 = con - head_temp3;
		  a11 = con - a11;
		  a21 = head_temp3 - a11;
		  con = T_element * split;
		  b1 = con - T_element;
		  b1 = con - b1;
		  b2 = T_element - b1;

		  c11 = head_temp3 * T_element;
		  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		  c2 = tail_temp3 * T_element;
		  t1 = c11 + c2;
		  t2 = (c2 - (t1 - c11)) + c21;

		  head_temp2 = t1 + t2;
		  tail_temp2 = t2 - (head_temp2 - t1);
		}
		{
		  double head_bt, tail_bt;
		  head_bt = -head_temp2;
		  tail_bt = -tail_temp2;
		  {
		    /* Compute double-double = double-double + double-double. */
		    double bv;
		    double s1, s2, t1, t2;

		    /* Add two hi words. */
		    s1 = head_temp1 + head_bt;
		    bv = s1 - head_temp1;
		    s2 = ((head_bt - bv) + (head_temp1 - (s1 - bv)));

		    /* Add two lo words. */
		    t1 = tail_temp1 + tail_bt;
		    bv = t1 - tail_temp1;
		    t2 = ((tail_bt - bv) + (tail_temp1 - (t1 - bv)));

		    s2 += t1;

		    /* Renormalize (s1, s2)  to  (t1, s2) */
		    t1 = s1 + s2;
		    s2 = s2 - (t1 - s1);

		    t2 += s2;

		    /* Renormalize (t1, t2)  */
		    head_temp1 = t1 + t2;
		    tail_temp1 = t2 - (head_temp1 - t1);
		  }
		}
		ix += inc_intx;
	      }			/* for i<j */

	      /* if the diagonal entry is not equal to one, then divide Xj by 
	         the entry */
	      if (diag == blas_non_unit_diag) {
		T_element = T_i[j * incT + j * ldt * incT];


		{
		  /* Compute double-double = double-double / double,
		     using a Newton iteration scheme. */
		  double b1, b2, con, e, t1, t2, t11, t21, t12, t22;

		  /* Compute a DP approximation to the quotient. */
		  t1 = head_temp1 / T_element;

		  /* Split t1 and b into two parts with at most 26 bits each,
		     using the Dekker-Veltkamp method. */
		  con = t1 * split;
		  t11 = con - (con - t1);
		  t21 = t1 - t11;
		  con = T_element * split;
		  b1 = con - (con - T_element);
		  b2 = T_element - b1;

		  /* Compute t1 * b using Dekker method. */
		  t12 = t1 * T_element;
		  t22 = (((t11 * b1 - t12) + t11 * b2) + t21 * b1) + t21 * b2;

		  /* Compute dda - (t12, t22) using Knuth trick. */
		  t11 = head_temp1 - t12;
		  e = t11 - head_temp1;
		  t21 =
		    ((-t12 - e) + (head_temp1 - (t11 - e))) + tail_temp1 -
		    t22;

		  /* Compute high-order word of (t11, t21) and divide by b. */
		  t2 = (t11 + t21) / T_element;

		  /* The result is t1 + t2, after normalization. */
		  head_temp1 = t1 + t2;
		  tail_temp1 = t2 - (head_temp1 - t1);
		}

	      }
	      /* if (diag == blas_non_unit_diag) */
	      head_intx[jx] = head_temp1;
	      tail_intx[jx] = tail_temp1;
	      jx += inc_intx;
	    }			/* for j<n */
	  }

	  /* copy the final results from intx to x */
	  ix = start_x;
	  jx = 0;
	  for (i = 0; i < n; i++) {
	    head_temp1 = head_intx[jx];
	    tail_temp1 = tail_intx[jx];
	    x_i[ix] = head_temp1;
	    ix += incx;
	    jx += inc_intx;
	  }

	  blas_free(head_intx);
	  blas_free(tail_intx);
	}
      }
      FPU_FIX_STOP;
    }
    break;
  }
}
