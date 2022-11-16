//
//  Copyright (C) 2002, 2003 Si-Lab b.v.b.a and Toon Knapen
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_BLAS_BLAS_H
#define BOOST_NUMERIC_BINDINGS_BLAS_BLAS_H

/*
 * const-correct prototypes for BLAS functions
 *
 */

#include <boost/numeric/bindings/blas/blas_names.h>
#include <boost/numeric/bindings/traits/type.h>

extern "C"
{
  //
  // Level 1
  //
  void   BLAS_SSCAL(const integer_t *n, const float*      alpha, float*      x, const integer_t* incx);
  void   BLAS_DSCAL(const integer_t *n, const double*     alpha, double*     x, const integer_t* incx);
  void   BLAS_CSCAL(const integer_t *n, const fcomplex_t* alpha, fcomplex_t* x, const integer_t* incx);
  void   BLAS_ZSCAL(const integer_t *n, const dcomplex_t* alpha, dcomplex_t* x, const integer_t* incx);

  void   BLAS_SAXPY(const integer_t *n, const float*      alpha, const float*      x, const integer_t* incx,  float*      y, const integer_t* incy);
  void   BLAS_DAXPY(const integer_t *n, const double*     alpha, const double*     x, const integer_t* incx,  double*     y, const integer_t* incy);
  void   BLAS_CAXPY(const integer_t *n, const fcomplex_t* alpha, const fcomplex_t* x, const integer_t* incx,  fcomplex_t* y, const integer_t* incy);
  void   BLAS_ZAXPY(const integer_t *n, const dcomplex_t* alpha, const dcomplex_t* x, const integer_t* incx,  dcomplex_t* y, const integer_t* incy);

#ifndef BIND_FORTRAN_F2C_RETURN_CONVENTIONS
  float  BLAS_SDOT (const integer_t *n, const float  *x, const integer_t *incx, const float  *y, const integer_t *incy);
#else
  double BLAS_SDOT (const integer_t *n, const float  *x, const integer_t *incx, const float  *y, const integer_t *incy);
#endif
  double BLAS_DDOT (const integer_t *n, const double *x, const integer_t *incx, const double *y, const integer_t *incy);

#ifndef BIND_FORTRAN_F2C_RETURN_CONVENTIONS
  fcomplex_t BLAS_CDOTU(const integer_t *n, const fcomplex_t *x, const integer_t *incx, const fcomplex_t *y, const integer_t *incy);
#else
  void   BLAS_CDOTU(fcomplex_t* ret, const integer_t *n, const fcomplex_t *x, const integer_t *incx, const fcomplex_t *y, const integer_t *incy);
#endif
#ifndef BIND_FORTRAN_F2C_RETURN_CONVENTIONS
  dcomplex_t BLAS_ZDOTU(const integer_t *n, const dcomplex_t *x, const integer_t *incx, const dcomplex_t *y, const integer_t *incy);
#else
  void   BLAS_ZDOTU(dcomplex_t* ret, const integer_t *n, const dcomplex_t *x, const integer_t *incx, const dcomplex_t *y, const integer_t *incy);
#endif

#ifndef BIND_FORTRAN_F2C_RETURN_CONVENTIONS
  fcomplex_t BLAS_CDOTC(const integer_t *n, const fcomplex_t *x, const integer_t *incx, const fcomplex_t *y, const integer_t *incy);
#else
  void   BLAS_CDOTC(fcomplex_t* ret, const integer_t *n, const fcomplex_t *x, const integer_t *incx, const fcomplex_t *y, const integer_t *incy);
#endif
#ifndef BIND_FORTRAN_F2C_RETURN_CONVENTIONS
  dcomplex_t BLAS_ZDOTC(const integer_t *n, const dcomplex_t *x, const integer_t *incx, const dcomplex_t *y, const integer_t *incy);
#else
  void   BLAS_ZDOTC(dcomplex_t* ret, const integer_t *n, const dcomplex_t *x, const integer_t *incx, const dcomplex_t *y, const integer_t *incy);
#endif

#ifndef BIND_FORTRAN_F2C_RETURN_CONVENTIONS
  float   BLAS_SNRM2(const integer_t *n, const float  *x, const integer_t *incx);
#else
  double  BLAS_SNRM2(const integer_t *n, const float  *x, const integer_t *incx);
#endif
  double  BLAS_DNRM2(const integer_t *n, const double *x, const integer_t *incx);
#ifndef BIND_FORTRAN_F2C_RETURN_CONVENTIONS
  float   BLAS_SCNRM2(const integer_t *n, const fcomplex_t  *x, const integer_t *incx);
#else
  double  BLAS_SCNRM2(const integer_t *n, const fcomplex_t  *x, const integer_t *incx);
#endif
  double  BLAS_DZNRM2(const integer_t *n, const dcomplex_t *x, const integer_t *incx);

#ifndef BIND_FORTRAN_F2C_RETURN_CONVENTIONS
  float   BLAS_SASUM(const integer_t *n, const float  *x, const integer_t *incx);
#else
  double  BLAS_SASUM(const integer_t *n, const float  *x, const integer_t *incx);
#endif
  double  BLAS_DASUM(const integer_t *n, const double *x, const integer_t *incx);
#ifndef BIND_FORTRAN_F2C_RETURN_CONVENTIONS
  float   BLAS_SCASUM(const integer_t *n, const fcomplex_t  *x, const integer_t *incx);
#else
  double  BLAS_SCASUM(const integer_t *n, const fcomplex_t  *x, const integer_t *incx);
#endif
  double  BLAS_DZASUM(const integer_t *n, const dcomplex_t *x, const integer_t *incx);

  void BLAS_SCOPY( const integer_t *n, const float  *x, const integer_t *incx, float  *y, const integer_t *incy);
  void BLAS_DCOPY( const integer_t *n, const double  *x, const integer_t *incx, double  *y, const integer_t *incy);
  void BLAS_CCOPY( const integer_t *n, const fcomplex_t  *x, const integer_t *incx, fcomplex_t  *y, const integer_t *incy);
  void BLAS_ZCOPY( const integer_t *n, const dcomplex_t  *x, const integer_t *incx, dcomplex_t  *y, const integer_t *incy);


  //
  // Level 2
  //
  void   BLAS_SGEMV(const char *trans, const integer_t *m, const integer_t *n, const float    *alpha, const float    *a, const integer_t *lda, const float    *x, const integer_t *incx, const float    *beta, float    *y, const integer_t *incy) ;
  void   BLAS_DGEMV(const char *trans, const integer_t *m, const integer_t *n, const double   *alpha, const double   *a, const integer_t *lda, const double   *x, const integer_t *incx, const double   *beta, double   *y, const integer_t *incy) ;
  void   BLAS_CGEMV(const char *trans, const integer_t *m, const integer_t *n, const fcomplex_t *alpha, const fcomplex_t *a, const integer_t *lda, const fcomplex_t *x, const integer_t *incx, const fcomplex_t *beta, fcomplex_t *y, const integer_t *incy) ;
  void   BLAS_ZGEMV(const char *trans, const integer_t *m, const integer_t *n, const dcomplex_t *alpha, const dcomplex_t *a, const integer_t *lda, const dcomplex_t *x, const integer_t *incx, const dcomplex_t *beta, dcomplex_t *y, const integer_t *incy) ;

  void   BLAS_SGER (const integer_t *m, const integer_t *n, const float * alpha, const float * x, const integer_t *incx, const float * y, const integer_t *incy, float *a,  const integer_t *lda);
  void   BLAS_DGER (const integer_t *m, const integer_t *n, const double *alpha, const double *x, const integer_t *incx, const double *y, const integer_t *incy, double *a, const integer_t *lda);

  void   BLAS_CGERU(const integer_t *m, const integer_t *n, const fcomplex_t *alpha, const fcomplex_t *x, const integer_t *incx, const fcomplex_t *y, const integer_t *incy, fcomplex_t *a, const integer_t *lda);
  void   BLAS_ZGERU(const integer_t *m, const integer_t *n, const dcomplex_t *alpha, const dcomplex_t *x, const integer_t *incx, const dcomplex_t *y, const integer_t *incy, dcomplex_t *a, const integer_t *lda);

  void   BLAS_CGERC(const integer_t *m, const integer_t *n, const fcomplex_t *alpha, const fcomplex_t *x, const integer_t *incx, const fcomplex_t *y, const integer_t *incy, fcomplex_t *a, const integer_t *lda);
  void   BLAS_ZGERC(const integer_t *m, const integer_t *n, const dcomplex_t *alpha, const dcomplex_t *x, const integer_t *incx, const dcomplex_t *y, const integer_t *incy, dcomplex_t *a, const integer_t *lda);


  //
  // Level 3
  //
  void   BLAS_SGEMM(const char *transa, const char *transb, const integer_t *m, const integer_t *n, const integer_t *k, const float      *alpha, const float      *a, const integer_t *lda, const float      *b, const integer_t *ldb, const float      *beta, float      *c, const integer_t *ldc);
  void   BLAS_DGEMM(const char *transa, const char *transb, const integer_t *m, const integer_t *n, const integer_t *k, const double     *alpha, const double     *a, const integer_t *lda, const double     *b, const integer_t *ldb, const double     *beta, double     *c, const integer_t *ldc);
  void   BLAS_CGEMM(const char *transa, const char *transb, const integer_t *m, const integer_t *n, const integer_t *k, const fcomplex_t *alpha, const fcomplex_t *a, const integer_t *lda, const fcomplex_t *b, const integer_t *ldb, const fcomplex_t *beta, fcomplex_t *c, const integer_t *ldc);
  void   BLAS_ZGEMM(const char *transa, const char *transb, const integer_t *m, const integer_t *n, const integer_t *k, const dcomplex_t *alpha, const dcomplex_t *a, const integer_t *lda, const dcomplex_t *b, const integer_t *ldb, const dcomplex_t *beta, dcomplex_t *c, const integer_t *ldc);

  void   BLAS_SSYRK ( const char* uplo, const char* trans, const integer_t* n, const integer_t* k, const float* alpha,
                      const float* a, const integer_t* lda, const float* beta, float* c, const integer_t* ldc );
  void   BLAS_DSYRK ( const char* uplo, const char* trans, const integer_t* n, const integer_t* k, const double* alpha,
                      const double* a, const integer_t* lda, const double* beta, double* c, const integer_t* ldc );
  void   BLAS_CSYRK ( const char* uplo, const char* trans, const integer_t* n, const integer_t* k, const fcomplex_t* alpha,
                      const fcomplex_t* a, const integer_t* lda, const fcomplex_t* beta, fcomplex_t* c, const integer_t* ldc );
  void   BLAS_ZSYRK ( const char* uplo, const char* trans, const integer_t* n, const integer_t* k, const dcomplex_t* alpha,
                      const dcomplex_t* a, const integer_t* lda, const dcomplex_t* beta, dcomplex_t* c, const integer_t* ldc );
  void   BLAS_CHERK ( const char* uplo, const char* trans, const integer_t* n, const integer_t* k, const float* alpha,
                      const fcomplex_t* a, const integer_t* lda, const float* beta, fcomplex_t* c, const integer_t* ldc );
  void   BLAS_ZHERK ( const char* uplo, const char* trans, const integer_t* n, const integer_t* k, const double* alpha,
                      const dcomplex_t* a, const integer_t* lda, const double* beta, dcomplex_t* c, const integer_t* ldc );

  void BLAS_STRSM( const char* side, const char* uplo, const char* transa, const char* diag, const integer_t* m,
                   const integer_t* n, float const* alpha, float const* a, integer_t const* lda, float* b, integer_t const* ldb );
  void BLAS_DTRSM( const char* side, const char* uplo, const char* transa, const char* diag, const integer_t* m,
                   const integer_t* n, double const* alpha, double const* a, integer_t const* lda, double* b, integer_t const* ldb );
  void BLAS_CTRSM( const char* side, const char* uplo, const char* transa, const char* diag, const integer_t* m,
                   const integer_t* n, fcomplex_t const* alpha, fcomplex_t const* a, integer_t const* lda, fcomplex_t* b, integer_t const* ldb );
  void BLAS_ZTRSM( const char* side, const char* uplo, const char* transa, const char* diag, const integer_t* m,
                   const integer_t* n, dcomplex_t const* alpha, dcomplex_t const* a, integer_t const* lda, dcomplex_t* b, integer_t const* ldb );

}

#endif // BOOST_NUMERIC_BINDINGS_BLAS_BLAS_H
