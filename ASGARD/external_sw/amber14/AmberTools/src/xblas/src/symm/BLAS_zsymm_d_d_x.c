#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_zsymm_d_d_x(enum blas_order_type order, enum blas_side_type side,
		      enum blas_uplo_type uplo, int m, int n,
		      const void *alpha, const double *a, int lda,
		      const double *b, int ldb, const void *beta,
		      void *c, int ldc, enum blas_prec_type prec)

/* 
 * Purpose
 * =======
 *
 * This routines computes one of the matrix product:
 *
 *     C  <-  alpha * A * B  +  beta * C
 *     C  <-  alpha * B * A  +  beta * C
 * 
 * where A is a symmetric matrix.
 *
 * Arguments
 * =========
 *
 * order   (input) enum blas_order_type
 *         Storage format of input matrices A, B, and C.
 * 
 * side    (input) enum blas_side_type
 *         Determines which side of matrix B is matrix A is multiplied.
 *
 * uplo    (input) enum blas_uplo_type
 *         Determines which half of matrix A (upper or lower triangle)
 *         is accessed.
 *
 * m n     (input) int
 *         Size of matrices A, B, and C.
 *         Matrix A is m-by-m if it is multiplied on the left, 
 *                     n-by-n otherwise.
 *         Matrices B and C are m-by-n.
 *
 * alpha   (input) const void*
 * 
 * a       (input) const double*
 *         Matrix A.
 *
 * lda     (input) int
 *         Leading dimension of matrix A.
 *
 * b       (input) const double*
 *         Matrix B.
 *   
 * ldb     (input) int
 *         Leading dimension of matrix B.
 *
 * beta    (input) const void*
 * 
 * c       (input/output) void*
 *         Matrix C.
 * 
 * ldc     (input) int
 *         Leading dimension of matrix C.
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
  switch (prec) {

  case blas_prec_single:
  case blas_prec_double:
  case blas_prec_indigenous:{

      /* Integer Index Variables */
      int i, j, k;

      int ai, bj, ci;
      int aik, bkj, cij;

      int incai, incbj, incci;
      int incaik1, incaik2, incbkj, inccij;

      int m_i, n_i;

      /* Input Matrices */
      const double *a_i = a;
      const double *b_i = b;

      /* Output Matrix */
      double *c_i = (double *) c;

      /* Input Scalars */
      double *alpha_i = (double *) alpha;
      double *beta_i = (double *) beta;

      /* Temporary Floating-Point Variables */
      double a_elem;
      double b_elem;
      double c_elem[2];
      double prod;
      double sum;
      double tmp1[2];
      double tmp2[2];



      /* Check for error conditions. */
      if (m <= 0 || n <= 0) {
	return;
      }

      if (order == blas_colmajor && (ldb < m || ldc < m)) {
	return;
      }
      if (order == blas_rowmajor && (ldb < n || ldc < n)) {
	return;
      }
      if (side == blas_left_side && lda < m) {
	return;
      }
      if (side == blas_right_side && lda < n) {
	return;
      }

      /* Test for no-op */
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0
	  && (beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
	return;
      }

      /* Set Index Parameters */
      if (side == blas_left_side) {
	m_i = m;
	n_i = n;
      } else {
	m_i = n;
	n_i = m;
      }

      if ((order == blas_colmajor && side == blas_left_side) ||
	  (order == blas_rowmajor && side == blas_right_side)) {
	incbj = ldb;
	incbkj = 1;
	incci = 1;
	inccij = ldc;
      } else {
	incbj = 1;
	incbkj = ldb;
	incci = ldc;
	inccij = 1;
      }

      if ((order == blas_colmajor && uplo == blas_upper) ||
	  (order == blas_rowmajor && uplo == blas_lower)) {
	incai = lda;
	incaik1 = 1;
	incaik2 = lda;
      } else {
	incai = 1;
	incaik1 = lda;
	incaik2 = 1;
      }



      /* Adjustment to increments (if any) */
      incci *= 2;
      inccij *= 2;






      /* alpha = 0.  In this case, just return beta * C */
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
	for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
	  for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {
	    c_elem[0] = c_i[cij];
	    c_elem[1] = c_i[cij + 1];
	    {
	      tmp1[0] =
		(double) c_elem[0] * beta_i[0] -
		(double) c_elem[1] * beta_i[1];
	      tmp1[1] =
		(double) c_elem[0] * beta_i[1] +
		(double) c_elem[1] * beta_i[0];
	    }
	    c_i[cij] = tmp1[0];
	    c_i[cij + 1] = tmp1[1];
	  }
	}
      } else if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {

	/* Case alpha == 1. */

	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* Case alpha = 1, beta = 0.  We compute  C <--- A * B   or  B * A */
	  for (i = 0, ci = 0, ai = 0; i < m_i; i++, ci += incci, ai += incai) {
	    for (j = 0, cij = ci, bj = 0; j < n_i;
		 j++, cij += inccij, bj += incbj) {

	      sum = 0.0;

	      for (k = 0, aik = ai, bkj = bj; k < i;
		   k++, aik += incaik1, bkj += incbkj) {
		a_elem = a_i[aik];
		b_elem = b_i[bkj];
		prod = a_elem * b_elem;
		sum = sum + prod;
	      }
	      for (; k < m_i; k++, aik += incaik2, bkj += incbkj) {
		a_elem = a_i[aik];
		b_elem = b_i[bkj];
		prod = a_elem * b_elem;
		sum = sum + prod;
	      }
	      tmp1[0] = sum;
	      tmp1[1] = 0.0;
	      c_i[cij] = tmp1[0];
	      c_i[cij + 1] = tmp1[1];
	    }
	  }
	} else {
	  /* Case alpha = 1, but beta != 0. 
	     We compute  C  <--- A * B + beta * C 
	     or  C  <--- B * A + beta * C  */

	  for (i = 0, ci = 0, ai = 0; i < m_i; i++, ci += incci, ai += incai) {
	    for (j = 0, cij = ci, bj = 0; j < n_i;
		 j++, cij += inccij, bj += incbj) {

	      sum = 0.0;

	      for (k = 0, aik = ai, bkj = bj; k < i;
		   k++, aik += incaik1, bkj += incbkj) {
		a_elem = a_i[aik];
		b_elem = b_i[bkj];
		prod = a_elem * b_elem;
		sum = sum + prod;
	      }
	      for (; k < m_i; k++, aik += incaik2, bkj += incbkj) {
		a_elem = a_i[aik];
		b_elem = b_i[bkj];
		prod = a_elem * b_elem;
		sum = sum + prod;
	      }
	      c_elem[0] = c_i[cij];
	      c_elem[1] = c_i[cij + 1];
	      {
		tmp2[0] =
		  (double) c_elem[0] * beta_i[0] -
		  (double) c_elem[1] * beta_i[1];
		tmp2[1] =
		  (double) c_elem[0] * beta_i[1] +
		  (double) c_elem[1] * beta_i[0];
	      }
	      tmp1[0] = sum;
	      tmp1[1] = 0.0;
	      tmp1[0] = tmp2[0] + tmp1[0];
	      tmp1[1] = tmp2[1] + tmp1[1];
	      c_i[cij] = tmp1[0];
	      c_i[cij + 1] = tmp1[1];
	    }
	  }
	}

      } else {
	/* The most general form,   C <--- alpha * A * B + beta * C  
	   or   C <--- alpha * B * A + beta * C  */

	for (i = 0, ci = 0, ai = 0; i < m_i; i++, ci += incci, ai += incai) {
	  for (j = 0, cij = ci, bj = 0; j < n_i;
	       j++, cij += inccij, bj += incbj) {

	    sum = 0.0;

	    for (k = 0, aik = ai, bkj = bj; k < i;
		 k++, aik += incaik1, bkj += incbkj) {
	      a_elem = a_i[aik];
	      b_elem = b_i[bkj];
	      prod = a_elem * b_elem;
	      sum = sum + prod;
	    }
	    for (; k < m_i; k++, aik += incaik2, bkj += incbkj) {
	      a_elem = a_i[aik];
	      b_elem = b_i[bkj];
	      prod = a_elem * b_elem;
	      sum = sum + prod;
	    }
	    {
	      tmp1[0] = alpha_i[0] * sum;
	      tmp1[1] = alpha_i[1] * sum;
	    }
	    c_elem[0] = c_i[cij];
	    c_elem[1] = c_i[cij + 1];
	    {
	      tmp2[0] =
		(double) c_elem[0] * beta_i[0] -
		(double) c_elem[1] * beta_i[1];
	      tmp2[1] =
		(double) c_elem[0] * beta_i[1] +
		(double) c_elem[1] * beta_i[0];
	    }
	    tmp1[0] = tmp1[0] + tmp2[0];
	    tmp1[1] = tmp1[1] + tmp2[1];
	    c_i[cij] = tmp1[0];
	    c_i[cij + 1] = tmp1[1];
	  }
	}
      }



      break;
    }

  case blas_prec_extra:{

      /* Integer Index Variables */
      int i, j, k;

      int ai, bj, ci;
      int aik, bkj, cij;

      int incai, incbj, incci;
      int incaik1, incaik2, incbkj, inccij;

      int m_i, n_i;

      /* Input Matrices */
      const double *a_i = a;
      const double *b_i = b;

      /* Output Matrix */
      double *c_i = (double *) c;

      /* Input Scalars */
      double *alpha_i = (double *) alpha;
      double *beta_i = (double *) beta;

      /* Temporary Floating-Point Variables */
      double a_elem;
      double b_elem;
      double c_elem[2];
      double head_prod, tail_prod;
      double head_sum, tail_sum;
      double head_tmp1[2], tail_tmp1[2];
      double head_tmp2[2], tail_tmp2[2];

      FPU_FIX_DECL;

      /* Check for error conditions. */
      if (m <= 0 || n <= 0) {
	return;
      }

      if (order == blas_colmajor && (ldb < m || ldc < m)) {
	return;
      }
      if (order == blas_rowmajor && (ldb < n || ldc < n)) {
	return;
      }
      if (side == blas_left_side && lda < m) {
	return;
      }
      if (side == blas_right_side && lda < n) {
	return;
      }

      /* Test for no-op */
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0
	  && (beta_i[0] == 1.0 && beta_i[1] == 0.0)) {
	return;
      }

      /* Set Index Parameters */
      if (side == blas_left_side) {
	m_i = m;
	n_i = n;
      } else {
	m_i = n;
	n_i = m;
      }

      if ((order == blas_colmajor && side == blas_left_side) ||
	  (order == blas_rowmajor && side == blas_right_side)) {
	incbj = ldb;
	incbkj = 1;
	incci = 1;
	inccij = ldc;
      } else {
	incbj = 1;
	incbkj = ldb;
	incci = ldc;
	inccij = 1;
      }

      if ((order == blas_colmajor && uplo == blas_upper) ||
	  (order == blas_rowmajor && uplo == blas_lower)) {
	incai = lda;
	incaik1 = 1;
	incaik2 = lda;
      } else {
	incai = 1;
	incaik1 = lda;
	incaik2 = 1;
      }

      FPU_FIX_START;

      /* Adjustment to increments (if any) */
      incci *= 2;
      inccij *= 2;






      /* alpha = 0.  In this case, just return beta * C */
      if (alpha_i[0] == 0.0 && alpha_i[1] == 0.0) {
	for (i = 0, ci = 0; i < m_i; i++, ci += incci) {
	  for (j = 0, cij = ci; j < n_i; j++, cij += inccij) {
	    c_elem[0] = c_i[cij];
	    c_elem[1] = c_i[cij + 1];
	    {
	      /* Compute complex-extra = complex-double * complex-double. */
	      double head_t1, tail_t1;
	      double head_t2, tail_t2;
	      /* Real part */
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = c_elem[0] * split;
		a1 = con - c_elem[0];
		a1 = con - a1;
		a2 = c_elem[0] - a1;
		con = beta_i[0] * split;
		b1 = con - beta_i[0];
		b1 = con - b1;
		b2 = beta_i[0] - b1;

		head_t1 = c_elem[0] * beta_i[0];
		tail_t1 =
		  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = c_elem[1] * split;
		a1 = con - c_elem[1];
		a1 = con - a1;
		a2 = c_elem[1] - a1;
		con = beta_i[1] * split;
		b1 = con - beta_i[1];
		b1 = con - b1;
		b2 = beta_i[1] - b1;

		head_t2 = c_elem[1] * beta_i[1];
		tail_t2 =
		  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_t2 = -head_t2;
	      tail_t2 = -tail_t2;
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_tmp1[0] = head_t1;
	      tail_tmp1[0] = tail_t1;
	      /* Imaginary part */
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = c_elem[1] * split;
		a1 = con - c_elem[1];
		a1 = con - a1;
		a2 = c_elem[1] - a1;
		con = beta_i[0] * split;
		b1 = con - beta_i[0];
		b1 = con - b1;
		b2 = beta_i[0] - b1;

		head_t1 = c_elem[1] * beta_i[0];
		tail_t1 =
		  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = c_elem[0] * split;
		a1 = con - c_elem[0];
		a1 = con - a1;
		a2 = c_elem[0] - a1;
		con = beta_i[1] * split;
		b1 = con - beta_i[1];
		b1 = con - b1;
		b2 = beta_i[1] - b1;

		head_t2 = c_elem[0] * beta_i[1];
		tail_t2 =
		  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_tmp1[1] = head_t1;
	      tail_tmp1[1] = tail_t1;
	    }
	    c_i[cij] = head_tmp1[0];
	    c_i[cij + 1] = head_tmp1[1];
	  }
	}
      } else if ((alpha_i[0] == 1.0 && alpha_i[1] == 0.0)) {

	/* Case alpha == 1. */

	if (beta_i[0] == 0.0 && beta_i[1] == 0.0) {
	  /* Case alpha = 1, beta = 0.  We compute  C <--- A * B   or  B * A */
	  for (i = 0, ci = 0, ai = 0; i < m_i; i++, ci += incci, ai += incai) {
	    for (j = 0, cij = ci, bj = 0; j < n_i;
		 j++, cij += inccij, bj += incbj) {

	      head_sum = tail_sum = 0.0;

	      for (k = 0, aik = ai, bkj = bj; k < i;
		   k++, aik += incaik1, bkj += incbkj) {
		a_elem = a_i[aik];
		b_elem = b_i[bkj];
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = a_elem * split;
		  a1 = con - a_elem;
		  a1 = con - a1;
		  a2 = a_elem - a1;
		  con = b_elem * split;
		  b1 = con - b_elem;
		  b1 = con - b1;
		  b2 = b_elem - b1;

		  head_prod = a_elem * b_elem;
		  tail_prod =
		    (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_sum + head_prod;
		  bv = s1 - head_sum;
		  s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_sum + tail_prod;
		  bv = t1 - tail_sum;
		  t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_sum = t1 + t2;
		  tail_sum = t2 - (head_sum - t1);
		}
	      }
	      for (; k < m_i; k++, aik += incaik2, bkj += incbkj) {
		a_elem = a_i[aik];
		b_elem = b_i[bkj];
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = a_elem * split;
		  a1 = con - a_elem;
		  a1 = con - a1;
		  a2 = a_elem - a1;
		  con = b_elem * split;
		  b1 = con - b_elem;
		  b1 = con - b1;
		  b2 = b_elem - b1;

		  head_prod = a_elem * b_elem;
		  tail_prod =
		    (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_sum + head_prod;
		  bv = s1 - head_sum;
		  s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_sum + tail_prod;
		  bv = t1 - tail_sum;
		  t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_sum = t1 + t2;
		  tail_sum = t2 - (head_sum - t1);
		}
	      }
	      head_tmp1[0] = head_sum;
	      tail_tmp1[0] = tail_sum;
	      head_tmp1[1] = tail_tmp1[1] = 0.0;
	      c_i[cij] = head_tmp1[0];
	      c_i[cij + 1] = head_tmp1[1];
	    }
	  }
	} else {
	  /* Case alpha = 1, but beta != 0. 
	     We compute  C  <--- A * B + beta * C 
	     or  C  <--- B * A + beta * C  */

	  for (i = 0, ci = 0, ai = 0; i < m_i; i++, ci += incci, ai += incai) {
	    for (j = 0, cij = ci, bj = 0; j < n_i;
		 j++, cij += inccij, bj += incbj) {

	      head_sum = tail_sum = 0.0;

	      for (k = 0, aik = ai, bkj = bj; k < i;
		   k++, aik += incaik1, bkj += incbkj) {
		a_elem = a_i[aik];
		b_elem = b_i[bkj];
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = a_elem * split;
		  a1 = con - a_elem;
		  a1 = con - a1;
		  a2 = a_elem - a1;
		  con = b_elem * split;
		  b1 = con - b_elem;
		  b1 = con - b1;
		  b2 = b_elem - b1;

		  head_prod = a_elem * b_elem;
		  tail_prod =
		    (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_sum + head_prod;
		  bv = s1 - head_sum;
		  s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_sum + tail_prod;
		  bv = t1 - tail_sum;
		  t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_sum = t1 + t2;
		  tail_sum = t2 - (head_sum - t1);
		}
	      }
	      for (; k < m_i; k++, aik += incaik2, bkj += incbkj) {
		a_elem = a_i[aik];
		b_elem = b_i[bkj];
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = a_elem * split;
		  a1 = con - a_elem;
		  a1 = con - a1;
		  a2 = a_elem - a1;
		  con = b_elem * split;
		  b1 = con - b_elem;
		  b1 = con - b1;
		  b2 = b_elem - b1;

		  head_prod = a_elem * b_elem;
		  tail_prod =
		    (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_sum + head_prod;
		  bv = s1 - head_sum;
		  s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_sum + tail_prod;
		  bv = t1 - tail_sum;
		  t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_sum = t1 + t2;
		  tail_sum = t2 - (head_sum - t1);
		}
	      }
	      c_elem[0] = c_i[cij];
	      c_elem[1] = c_i[cij + 1];
	      {
		/* Compute complex-extra = complex-double * complex-double. */
		double head_t1, tail_t1;
		double head_t2, tail_t2;
		/* Real part */
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = c_elem[0] * split;
		  a1 = con - c_elem[0];
		  a1 = con - a1;
		  a2 = c_elem[0] - a1;
		  con = beta_i[0] * split;
		  b1 = con - beta_i[0];
		  b1 = con - b1;
		  b2 = beta_i[0] - b1;

		  head_t1 = c_elem[0] * beta_i[0];
		  tail_t1 =
		    (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = c_elem[1] * split;
		  a1 = con - c_elem[1];
		  a1 = con - a1;
		  a2 = c_elem[1] - a1;
		  con = beta_i[1] * split;
		  b1 = con - beta_i[1];
		  b1 = con - b1;
		  b2 = beta_i[1] - b1;

		  head_t2 = c_elem[1] * beta_i[1];
		  tail_t2 =
		    (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		head_t2 = -head_t2;
		tail_t2 = -tail_t2;
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_t1 + head_t2;
		  bv = s1 - head_t1;
		  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_t1 + tail_t2;
		  bv = t1 - tail_t1;
		  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		head_tmp2[0] = head_t1;
		tail_tmp2[0] = tail_t1;
		/* Imaginary part */
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = c_elem[1] * split;
		  a1 = con - c_elem[1];
		  a1 = con - a1;
		  a2 = c_elem[1] - a1;
		  con = beta_i[0] * split;
		  b1 = con - beta_i[0];
		  b1 = con - b1;
		  b2 = beta_i[0] - b1;

		  head_t1 = c_elem[1] * beta_i[0];
		  tail_t1 =
		    (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		{
		  /* Compute double_double = double * double. */
		  double a1, a2, b1, b2, con;

		  con = c_elem[0] * split;
		  a1 = con - c_elem[0];
		  a1 = con - a1;
		  a2 = c_elem[0] - a1;
		  con = beta_i[1] * split;
		  b1 = con - beta_i[1];
		  b1 = con - b1;
		  b2 = beta_i[1] - b1;

		  head_t2 = c_elem[0] * beta_i[1];
		  tail_t2 =
		    (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
		}
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_t1 + head_t2;
		  bv = s1 - head_t1;
		  s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_t1 + tail_t2;
		  bv = t1 - tail_t1;
		  t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t1 = t1 + t2;
		  tail_t1 = t2 - (head_t1 - t1);
		}
		head_tmp2[1] = head_t1;
		tail_tmp2[1] = tail_t1;
	      }
	      head_tmp1[0] = head_sum;
	      tail_tmp1[0] = tail_sum;
	      head_tmp1[1] = tail_tmp1[1] = 0.0;
	      {
		double head_t, tail_t;
		double head_a, tail_a;
		double head_b, tail_b;
		/* Real part */
		head_a = head_tmp2[0];
		tail_a = tail_tmp2[0];
		head_b = head_tmp1[0];
		tail_b = tail_tmp1[0];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_tmp1[0] = head_t;
		tail_tmp1[0] = tail_t;
		/* Imaginary part */
		head_a = head_tmp2[1];
		tail_a = tail_tmp2[1];
		head_b = head_tmp1[1];
		tail_b = tail_tmp1[1];
		{
		  /* Compute double-double = double-double + double-double. */
		  double bv;
		  double s1, s2, t1, t2;

		  /* Add two hi words. */
		  s1 = head_a + head_b;
		  bv = s1 - head_a;
		  s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		  /* Add two lo words. */
		  t1 = tail_a + tail_b;
		  bv = t1 - tail_a;
		  t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		  s2 += t1;

		  /* Renormalize (s1, s2)  to  (t1, s2) */
		  t1 = s1 + s2;
		  s2 = s2 - (t1 - s1);

		  t2 += s2;

		  /* Renormalize (t1, t2)  */
		  head_t = t1 + t2;
		  tail_t = t2 - (head_t - t1);
		}
		head_tmp1[1] = head_t;
		tail_tmp1[1] = tail_t;
	      }
	      c_i[cij] = head_tmp1[0];
	      c_i[cij + 1] = head_tmp1[1];
	    }
	  }
	}

      } else {
	/* The most general form,   C <--- alpha * A * B + beta * C  
	   or   C <--- alpha * B * A + beta * C  */

	for (i = 0, ci = 0, ai = 0; i < m_i; i++, ci += incci, ai += incai) {
	  for (j = 0, cij = ci, bj = 0; j < n_i;
	       j++, cij += inccij, bj += incbj) {

	    head_sum = tail_sum = 0.0;

	    for (k = 0, aik = ai, bkj = bj; k < i;
		 k++, aik += incaik1, bkj += incbkj) {
	      a_elem = a_i[aik];
	      b_elem = b_i[bkj];
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = a_elem * split;
		a1 = con - a_elem;
		a1 = con - a1;
		a2 = a_elem - a1;
		con = b_elem * split;
		b1 = con - b_elem;
		b1 = con - b1;
		b2 = b_elem - b1;

		head_prod = a_elem * b_elem;
		tail_prod =
		  (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_sum + head_prod;
		bv = s1 - head_sum;
		s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sum + tail_prod;
		bv = t1 - tail_sum;
		t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sum = t1 + t2;
		tail_sum = t2 - (head_sum - t1);
	      }
	    }
	    for (; k < m_i; k++, aik += incaik2, bkj += incbkj) {
	      a_elem = a_i[aik];
	      b_elem = b_i[bkj];
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = a_elem * split;
		a1 = con - a_elem;
		a1 = con - a1;
		a2 = a_elem - a1;
		con = b_elem * split;
		b1 = con - b_elem;
		b1 = con - b1;
		b2 = b_elem - b1;

		head_prod = a_elem * b_elem;
		tail_prod =
		  (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_sum + head_prod;
		bv = s1 - head_sum;
		s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_sum + tail_prod;
		bv = t1 - tail_sum;
		t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_sum = t1 + t2;
		tail_sum = t2 - (head_sum - t1);
	      }
	    }
	    {
	      /* Compute complex-extra = complex-double * real. */
	      double head_t, tail_t;
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_sum * split;
		a11 = con - head_sum;
		a11 = con - a11;
		a21 = head_sum - a11;
		con = alpha_i[0] * split;
		b1 = con - alpha_i[0];
		b1 = con - b1;
		b2 = alpha_i[0] - b1;

		c11 = head_sum * alpha_i[0];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_sum * alpha_i[0];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_tmp1[0] = head_t;
	      tail_tmp1[0] = tail_t;
	      {
		/* Compute double-double = double-double * double. */
		double a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		con = head_sum * split;
		a11 = con - head_sum;
		a11 = con - a11;
		a21 = head_sum - a11;
		con = alpha_i[1] * split;
		b1 = con - alpha_i[1];
		b1 = con - b1;
		b2 = alpha_i[1] - b1;

		c11 = head_sum * alpha_i[1];
		c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		c2 = tail_sum * alpha_i[1];
		t1 = c11 + c2;
		t2 = (c2 - (t1 - c11)) + c21;

		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_tmp1[1] = head_t;
	      tail_tmp1[1] = tail_t;
	    }
	    c_elem[0] = c_i[cij];
	    c_elem[1] = c_i[cij + 1];
	    {
	      /* Compute complex-extra = complex-double * complex-double. */
	      double head_t1, tail_t1;
	      double head_t2, tail_t2;
	      /* Real part */
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = c_elem[0] * split;
		a1 = con - c_elem[0];
		a1 = con - a1;
		a2 = c_elem[0] - a1;
		con = beta_i[0] * split;
		b1 = con - beta_i[0];
		b1 = con - b1;
		b2 = beta_i[0] - b1;

		head_t1 = c_elem[0] * beta_i[0];
		tail_t1 =
		  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = c_elem[1] * split;
		a1 = con - c_elem[1];
		a1 = con - a1;
		a2 = c_elem[1] - a1;
		con = beta_i[1] * split;
		b1 = con - beta_i[1];
		b1 = con - b1;
		b2 = beta_i[1] - b1;

		head_t2 = c_elem[1] * beta_i[1];
		tail_t2 =
		  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      head_t2 = -head_t2;
	      tail_t2 = -tail_t2;
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_tmp2[0] = head_t1;
	      tail_tmp2[0] = tail_t1;
	      /* Imaginary part */
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = c_elem[1] * split;
		a1 = con - c_elem[1];
		a1 = con - a1;
		a2 = c_elem[1] - a1;
		con = beta_i[0] * split;
		b1 = con - beta_i[0];
		b1 = con - b1;
		b2 = beta_i[0] - b1;

		head_t1 = c_elem[1] * beta_i[0];
		tail_t1 =
		  (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double_double = double * double. */
		double a1, a2, b1, b2, con;

		con = c_elem[0] * split;
		a1 = con - c_elem[0];
		a1 = con - a1;
		a2 = c_elem[0] - a1;
		con = beta_i[1] * split;
		b1 = con - beta_i[1];
		b1 = con - b1;
		b2 = beta_i[1] - b1;

		head_t2 = c_elem[0] * beta_i[1];
		tail_t2 =
		  (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_t1 + head_t2;
		bv = s1 - head_t1;
		s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_t1 + tail_t2;
		bv = t1 - tail_t1;
		t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t1 = t1 + t2;
		tail_t1 = t2 - (head_t1 - t1);
	      }
	      head_tmp2[1] = head_t1;
	      tail_tmp2[1] = tail_t1;
	    }
	    {
	      double head_t, tail_t;
	      double head_a, tail_a;
	      double head_b, tail_b;
	      /* Real part */
	      head_a = head_tmp1[0];
	      tail_a = tail_tmp1[0];
	      head_b = head_tmp2[0];
	      tail_b = tail_tmp2[0];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_tmp1[0] = head_t;
	      tail_tmp1[0] = tail_t;
	      /* Imaginary part */
	      head_a = head_tmp1[1];
	      tail_a = tail_tmp1[1];
	      head_b = head_tmp2[1];
	      tail_b = tail_tmp2[1];
	      {
		/* Compute double-double = double-double + double-double. */
		double bv;
		double s1, s2, t1, t2;

		/* Add two hi words. */
		s1 = head_a + head_b;
		bv = s1 - head_a;
		s2 = ((head_b - bv) + (head_a - (s1 - bv)));

		/* Add two lo words. */
		t1 = tail_a + tail_b;
		bv = t1 - tail_a;
		t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

		s2 += t1;

		/* Renormalize (s1, s2)  to  (t1, s2) */
		t1 = s1 + s2;
		s2 = s2 - (t1 - s1);

		t2 += s2;

		/* Renormalize (t1, t2)  */
		head_t = t1 + t2;
		tail_t = t2 - (head_t - t1);
	      }
	      head_tmp1[1] = head_t;
	      tail_tmp1[1] = tail_t;
	    }
	    c_i[cij] = head_tmp1[0];
	    c_i[cij + 1] = head_tmp1[1];
	  }
	}
      }

      FPU_FIX_STOP;

      break;
    }
  }
}
