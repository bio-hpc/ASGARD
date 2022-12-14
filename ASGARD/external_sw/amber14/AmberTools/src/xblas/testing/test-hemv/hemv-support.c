#include <stdio.h>
#include "blas_extended.h"

void che_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_side_type side, int n, void *a, int lda,
		  void *a_vec, int row)

/* Copies the given row of matrix a into the supplied vector. */
{
  int conj_flag;
  int i;
  int ai, incai1, incai2;
  int vi, incvi = 1;

  float a_elem[2];
  const float *a_i = (float *) a;
  float *a_vec_i = (float *) a_vec;

  if ((side == blas_left_side && uplo == blas_upper) ||
      (side == blas_right_side && uplo == blas_lower))
    conj_flag = 1;
  else
    conj_flag = 0;

  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai1 = 1;
    incai2 = lda;
    ai = lda * row;
  } else {
    incai1 = lda;
    incai2 = 1;
    ai = row;
  }

  incai1 *= 2;
  incai2 *= 2;
  ai *= 2;
  incvi *= 2;

  for (i = 0, vi = 0; i < row; i++, ai += incai1, vi += incvi) {
    a_elem[0] = a_i[ai];
    a_elem[1] = a_i[ai + 1];
    if (conj_flag == 1) {
      a_elem[1] = -a_elem[1];
    }
    a_vec_i[vi] = a_elem[0];
    a_vec_i[vi + 1] = a_elem[1];
  }
  for (; i < n; i++, ai += incai2, vi += incvi) {
    a_elem[0] = a_i[ai];
    a_elem[1] = a_i[ai + 1];
    if (conj_flag == 0) {
      a_elem[1] = -a_elem[1];
    }
    a_vec_i[vi] = a_elem[0];
    a_vec_i[vi + 1] = a_elem[1];
  }
}
void zhe_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_side_type side, int n, void *a, int lda,
		  void *a_vec, int row)

/* Copies the given row of matrix a into the supplied vector. */
{
  int conj_flag;
  int i;
  int ai, incai1, incai2;
  int vi, incvi = 1;

  double a_elem[2];
  const double *a_i = (double *) a;
  double *a_vec_i = (double *) a_vec;

  if ((side == blas_left_side && uplo == blas_upper) ||
      (side == blas_right_side && uplo == blas_lower))
    conj_flag = 1;
  else
    conj_flag = 0;

  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai1 = 1;
    incai2 = lda;
    ai = lda * row;
  } else {
    incai1 = lda;
    incai2 = 1;
    ai = row;
  }

  incai1 *= 2;
  incai2 *= 2;
  ai *= 2;
  incvi *= 2;

  for (i = 0, vi = 0; i < row; i++, ai += incai1, vi += incvi) {
    a_elem[0] = a_i[ai];
    a_elem[1] = a_i[ai + 1];
    if (conj_flag == 1) {
      a_elem[1] = -a_elem[1];
    }
    a_vec_i[vi] = a_elem[0];
    a_vec_i[vi + 1] = a_elem[1];
  }
  for (; i < n; i++, ai += incai2, vi += incvi) {
    a_elem[0] = a_i[ai];
    a_elem[1] = a_i[ai + 1];
    if (conj_flag == 0) {
      a_elem[1] = -a_elem[1];
    }
    a_vec_i[vi] = a_elem[0];
    a_vec_i[vi + 1] = a_elem[1];
  }
}

void che_print_matrix(void *a, int n, int lda,
		      enum blas_order_type order, enum blas_uplo_type uplo)
{
  int ai, aij;
  int incai, incaij1, incaij2;
  int i, j;

  float a_elem[2];
  const float *a_i = (float *) a;

  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai = lda;
    incaij1 = 1;
    incaij2 = lda;
  } else {
    incai = 1;
    incaij1 = lda;
    incaij2 = 1;
  }

  incai *= 2;
  incaij1 *= 2;
  incaij2 *= 2;

  for (i = 0, ai = 0; i < n; i++, ai += incai) {

    for (j = 0, aij = ai; j < i; j++, aij += incaij1) {
      a_elem[0] = a_i[aij];
      a_elem[1] = a_i[aij + 1];
      if (uplo == blas_upper) {
	a_elem[1] = -a_elem[1];
      }
      printf("(%16.8e, %16.8e)", a_elem[0], a_elem[1]);
    }
    a_elem[0] = a_i[ai];
    a_elem[1] = 0.0;
    printf("(%16.8e, %16.8e)", a_elem[0], a_elem[1]);
    j++;
    aij += incaij2;
    for (; j < n; j++, aij += incaij2) {
      a_elem[0] = a_i[aij];
      a_elem[1] = a_i[aij + 1];
      if (uplo == blas_lower) {
	a_elem[1] = -a_elem[1];
      }
      printf("(%16.8e, %16.8e)", a_elem[0], a_elem[1]);
    }
    printf("\n");
  }

}
void zhe_print_matrix(void *a, int n, int lda,
		      enum blas_order_type order, enum blas_uplo_type uplo)
{
  int ai, aij;
  int incai, incaij1, incaij2;
  int i, j;

  double a_elem[2];
  const double *a_i = (double *) a;

  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai = lda;
    incaij1 = 1;
    incaij2 = lda;
  } else {
    incai = 1;
    incaij1 = lda;
    incaij2 = 1;
  }

  incai *= 2;
  incaij1 *= 2;
  incaij2 *= 2;

  for (i = 0, ai = 0; i < n; i++, ai += incai) {

    for (j = 0, aij = ai; j < i; j++, aij += incaij1) {
      a_elem[0] = a_i[aij];
      a_elem[1] = a_i[aij + 1];
      if (uplo == blas_upper) {
	a_elem[1] = -a_elem[1];
      }
      printf("(%24.16e, %24.16e)", a_elem[0], a_elem[1]);
    }
    a_elem[0] = a_i[ai];
    a_elem[1] = 0.0;
    printf("(%24.16e, %24.16e)", a_elem[0], a_elem[1]);
    j++;
    aij += incaij2;
    for (; j < n; j++, aij += incaij2) {
      a_elem[0] = a_i[aij];
      a_elem[1] = a_i[aij + 1];
      if (uplo == blas_lower) {
	a_elem[1] = -a_elem[1];
      }
      printf("(%24.16e, %24.16e)", a_elem[0], a_elem[1]);
    }
    printf("\n");
  }

}

void sskew_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      enum blas_side_type side, int n, float *a, int lda,
		      float *a_vec, int row)

/*
 *  Copies the given vector into the given row of symmetric matrix a.
 */
{
  int conj_flag;
  int i;
  int ai, incai1, incai2;
  int vi, incvi = 1;

  float a_elem;
  float *a_i = a;
  const float *a_vec_i = a_vec;

  if ((side == blas_left_side && uplo == blas_upper) ||
      (side == blas_right_side && uplo == blas_lower))
    conj_flag = 1;
  else
    conj_flag = 0;

  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai1 = 1;
    incai2 = lda;
    ai = lda * row;
  } else {
    incai1 = lda;
    incai2 = 1;
    ai = row;
  }






  for (i = 0, vi = 0; i < row; i++, ai += incai1, vi += incvi) {
    a_elem = a_vec_i[vi];
    if (conj_flag == 1) {
      a_elem = -a_elem;
    }
    a_i[ai] = a_elem;
  }

  for (; i < n; i++, ai += incai2, vi += incvi) {
    a_elem = a_vec_i[vi];
    if (conj_flag == 0) {
      a_elem = -a_elem;
    }
    a_i[ai] = a_elem;
  }
}
void dskew_commit_row(enum blas_order_type order, enum blas_uplo_type uplo,
		      enum blas_side_type side, int n, double *a, int lda,
		      double *a_vec, int row)

/*
 *  Copies the given vector into the given row of symmetric matrix a.
 */
{
  int conj_flag;
  int i;
  int ai, incai1, incai2;
  int vi, incvi = 1;

  double a_elem;
  double *a_i = a;
  const double *a_vec_i = a_vec;

  if ((side == blas_left_side && uplo == blas_upper) ||
      (side == blas_right_side && uplo == blas_lower))
    conj_flag = 1;
  else
    conj_flag = 0;

  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai1 = 1;
    incai2 = lda;
    ai = lda * row;
  } else {
    incai1 = lda;
    incai2 = 1;
    ai = row;
  }






  for (i = 0, vi = 0; i < row; i++, ai += incai1, vi += incvi) {
    a_elem = a_vec_i[vi];
    if (conj_flag == 1) {
      a_elem = -a_elem;
    }
    a_i[ai] = a_elem;
  }

  for (; i < n; i++, ai += incai2, vi += incvi) {
    a_elem = a_vec_i[vi];
    if (conj_flag == 0) {
      a_elem = -a_elem;
    }
    a_i[ai] = a_elem;
  }
}

void sskew_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_side_type side, int n, float *a, int lda,
		    float *a_vec, int row)

/*
 *  Copies the given row of skew matrix a into the supplied vector.
 */
{
  int conj_flag;
  int i;
  int ai, incai1, incai2;
  int vi, incvi = 1;

  float a_elem;
  const float *a_i = a;
  float *a_vec_i = a_vec;

  if ((side == blas_left_side && uplo == blas_upper) ||
      (side == blas_right_side && uplo == blas_lower))
    conj_flag = 1;
  else
    conj_flag = 0;

  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai1 = 1;
    incai2 = lda;
    ai = lda * row;
  } else {
    incai1 = lda;
    incai2 = 1;
    ai = row;
  }






  for (i = 0, vi = 0; i < row; i++, ai += incai1, vi += incvi) {
    a_elem = a_i[ai];
    if (conj_flag == 1) {
      a_elem = -a_elem;
    }
    a_vec_i[vi] = a_elem;
  }

  for (; i < n; i++, ai += incai2, vi += incvi) {
    a_elem = a_i[ai];
    if (conj_flag == 0) {
      a_elem = -a_elem;
    }
    a_vec_i[vi] = a_elem;
  }
}
void dskew_copy_row(enum blas_order_type order, enum blas_uplo_type uplo,
		    enum blas_side_type side, int n, double *a, int lda,
		    double *a_vec, int row)

/*
 *  Copies the given row of skew matrix a into the supplied vector.
 */
{
  int conj_flag;
  int i;
  int ai, incai1, incai2;
  int vi, incvi = 1;

  double a_elem;
  const double *a_i = a;
  double *a_vec_i = a_vec;

  if ((side == blas_left_side && uplo == blas_upper) ||
      (side == blas_right_side && uplo == blas_lower))
    conj_flag = 1;
  else
    conj_flag = 0;

  if ((order == blas_colmajor && uplo == blas_upper) ||
      (order == blas_rowmajor && uplo == blas_lower)) {
    incai1 = 1;
    incai2 = lda;
    ai = lda * row;
  } else {
    incai1 = lda;
    incai2 = 1;
    ai = row;
  }






  for (i = 0, vi = 0; i < row; i++, ai += incai1, vi += incvi) {
    a_elem = a_i[ai];
    if (conj_flag == 1) {
      a_elem = -a_elem;
    }
    a_vec_i[vi] = a_elem;
  }

  for (; i < n; i++, ai += incai2, vi += incvi) {
    a_elem = a_i[ai];
    if (conj_flag == 0) {
      a_elem = -a_elem;
    }
    a_vec_i[vi] = a_elem;
  }
}
