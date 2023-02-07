/*
 *
 * Copyright (c) Toon Knapen & Kresimir Fresl 2003
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * KF acknowledges the support of the Faculty of Civil Engineering,
 * University of Zagreb, Croatia.
 *
 */

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_LAPACK_H
#define BOOST_NUMERIC_BINDINGS_LAPACK_LAPACK_H

#include <boost/numeric/bindings/traits/type.h>
#include <boost/numeric/bindings/lapack/lapack_names.h>

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
#  define BOOST_NUMERIC_BINDINGS_FORTRAN
#endif

extern "C" {

  /********************************************************************/
  /*                        linear systems                            */
  /********************************************************************/

  /* general */

  void LAPACK_SGESV (integer_t const* n, integer_t const* nrhs,
                     float* a, integer_t const* lda, integer_t* ipiv,
                     float* b, integer_t const* ldb, integer_t* info);
  void LAPACK_DGESV (integer_t const* n, integer_t const* nrhs,
                     double* a, integer_t const* lda, integer_t* ipiv,
                     double* b, integer_t const* ldb, integer_t* info);
  void LAPACK_CGESV (integer_t const* n, integer_t const* nrhs,
                     fcomplex_t* a, integer_t const* lda, integer_t* ipiv,
                     fcomplex_t* b, integer_t const* ldb, integer_t* info);
  void LAPACK_ZGESV (integer_t const* n, integer_t const* nrhs,
                     dcomplex_t* a, integer_t const* lda, integer_t* ipiv,
                     dcomplex_t* b, integer_t const* ldb, integer_t* info);

  void LAPACK_SGETRF (integer_t const* n, integer_t const* nrhs,
                      float* a, integer_t const* lda, integer_t* ipiv, integer_t* info);
  void LAPACK_DGETRF (integer_t const* n, integer_t const* nrhs,
                      double* a, integer_t const* lda, integer_t* ipiv, integer_t* info);
  void LAPACK_CGETRF (integer_t const* n, integer_t const* nrhs,
                      fcomplex_t* a, integer_t const* lda,
                      integer_t* ipiv, integer_t* info);
  void LAPACK_ZGETRF (integer_t const* n, integer_t const* nrhs,
                      dcomplex_t* a, integer_t const* lda,
                      integer_t* ipiv, integer_t* info);

  void LAPACK_SGETRS (char const* trans, integer_t const* n, integer_t const* nrhs,
                      float const* a, integer_t const* lda, integer_t const* ipiv,
                      float* b, integer_t const* ldb, integer_t* info);
  void LAPACK_DGETRS (char const* trans, integer_t const* n, integer_t const* nrhs,
                      double const* a, integer_t const* lda, integer_t const* ipiv,
                      double* b, integer_t const* ldb, integer_t* info);
  void LAPACK_CGETRS (char const* trans, integer_t const* n, integer_t const* nrhs,
                      fcomplex_t const* a, integer_t const* lda, integer_t const* ipiv,
                      fcomplex_t* b, integer_t const* ldb, integer_t* info);
  void LAPACK_ZGETRS (char const* trans, integer_t const* n, integer_t const* nrhs,
                      dcomplex_t const* a, integer_t const* lda, integer_t const* ipiv,
                      dcomplex_t* b, integer_t const* ldb, integer_t* info);

  void LAPACK_SGETRI (integer_t const* n, float* a, integer_t const* lda, integer_t const* ipiv,
                      float* work, integer_t const* lwork, integer_t* info);
  void LAPACK_DGETRI (integer_t const* n, double* a, integer_t const* lda, integer_t const* ipiv,
                      double* work, integer_t const* lwork, integer_t* info);
  void LAPACK_CGETRI (integer_t const* n, fcomplex_t* a, integer_t const* lda, integer_t const* ipiv,
                      fcomplex_t* work, integer_t const* lwork, integer_t* info);
  void LAPACK_ZGETRI (integer_t const* n, dcomplex_t* a, integer_t const* lda, integer_t const* ipiv,
                      dcomplex_t* work, integer_t const* lwork, integer_t* info);

  /* symmetric/Hermitian positive definite */

  void LAPACK_SPOSV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     float* a, integer_t const* lda,
                     float* b, integer_t const* ldb, integer_t* info);
  void LAPACK_DPOSV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     double* a, integer_t const* lda,
                     double* b, integer_t const* ldb, integer_t* info);
  void LAPACK_CPOSV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     fcomplex_t* a, integer_t const* lda,
                     fcomplex_t* b, integer_t const* ldb, integer_t* info);
  void LAPACK_ZPOSV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     dcomplex_t* a, integer_t const* lda,
                     dcomplex_t* b, integer_t const* ldb, integer_t* info);

  void LAPACK_SPOTRF (char const* uplo, integer_t const* n,
                      float* a, integer_t const* lda, integer_t* info);
  void LAPACK_DPOTRF (char const* uplo, integer_t const* n,
                      double* a, integer_t const* lda, integer_t* info);
  void LAPACK_CPOTRF (char const* uplo, integer_t const* n,
                      fcomplex_t* a, integer_t const* lda, integer_t* info);
  void LAPACK_ZPOTRF (char const* uplo, integer_t const* n,
                      dcomplex_t* a, integer_t const* lda, integer_t* info);

  void LAPACK_SPOTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      float const* a, integer_t const* lda,
                      float* b, integer_t const* ldb, integer_t* info);
  void LAPACK_DPOTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      double const* a, integer_t const* lda,
                      double* b, integer_t const* ldb, integer_t* info);
  void LAPACK_CPOTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      fcomplex_t const* a, integer_t const* lda,
                      fcomplex_t* b, integer_t const* ldb, integer_t* info);
  void LAPACK_ZPOTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      dcomplex_t const* a, integer_t const* lda,
                      dcomplex_t* b, integer_t const* ldb, integer_t* info);


  /* symmetric/Hermitian positive definite in packed storage */

  void LAPACK_SPPSV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     float* ap, float* b, integer_t const* ldb, integer_t* info);
  void LAPACK_DPPSV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     double* ap, double* b, integer_t const* ldb, integer_t* info);
  void LAPACK_CPPSV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     fcomplex_t* ap, fcomplex_t* b, integer_t const* ldb, integer_t* info);
  void LAPACK_ZPPSV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     dcomplex_t* ap, dcomplex_t* b, integer_t const* ldb, integer_t* info);

  void LAPACK_SPPTRF (char const* uplo, integer_t const* n, float* ap, integer_t* info);
  void LAPACK_DPPTRF (char const* uplo, integer_t const* n, double* ap, integer_t* info);
  void LAPACK_CPPTRF (char const* uplo, integer_t const* n,
                      fcomplex_t* ap, integer_t* info);
  void LAPACK_ZPPTRF (char const* uplo, integer_t const* n,
                      dcomplex_t* ap, integer_t* info);

  void LAPACK_SPPTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      float const* ap, float* b, integer_t const* ldb, integer_t* info);
  void LAPACK_DPPTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      double const* ap, double* b, integer_t const* ldb, integer_t* info);
  void LAPACK_CPPTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      fcomplex_t const* ap,
                      fcomplex_t* b, integer_t const* ldb, integer_t* info);
  void LAPACK_ZPPTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      dcomplex_t const* ap,
                      dcomplex_t* b, integer_t const* ldb, integer_t* info);

  void LAPACK_SPPTRI (char const* uplo, integer_t const* n, float* ap, integer_t* info);
  void LAPACK_DPPTRI (char const* uplo, integer_t const* n, double* ap, integer_t* info);
  void LAPACK_CPPTRI (char const* uplo, integer_t const* n,
                      fcomplex_t* ap, integer_t* info);
  void LAPACK_ZPPTRI (char const* uplo, integer_t const* n,
                      dcomplex_t* ap, integer_t* info);

  /* symmetric/Hermitian positive definite tridiagonal */

  void LAPACK_SPTSV (integer_t const* N, integer_t const* NRHS, float* D, float* E
                    , float* B, integer_t const* LDB, integer_t* INFO
                    );
  void LAPACK_DPTSV (integer_t const* N, integer_t const* NRHS, double* D, double* E
                    , double* B, integer_t const* LDB, integer_t* INFO
                    );
  void LAPACK_CPTSV (integer_t const* N, integer_t const* NRHS, float* D, fcomplex_t* E
                    , fcomplex_t* B, integer_t const* LDB, integer_t* INFO
                    );
  void LAPACK_ZPTSV (integer_t const* N, integer_t const* NRHS, double* D, dcomplex_t* E
                    , dcomplex_t* B, integer_t const* LDB, integer_t* INFO
                    );

  void LAPACK_SPTTRF ( integer_t const* n, float* d, float* e, integer_t* info);
  void LAPACK_DPTTRF ( integer_t const* n, double* d, double* e, integer_t* info);
  void LAPACK_CPTTRF ( integer_t const* n, float* d, fcomplex_t* e, integer_t* info);
  void LAPACK_ZPTTRF ( integer_t const* n, double* d, dcomplex_t* e, integer_t* info);

  void LAPACK_SPTTRS ( integer_t const* n, integer_t const* nrhs,
                      float const* d, float const* e,
                      float* b, integer_t const* ldb, integer_t* info);
  void LAPACK_DPTTRS ( integer_t const* n, integer_t const* nrhs,
                      double const* d, double const* e,
                      double* b, integer_t const* ldb, integer_t* info);
  void LAPACK_CPTTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      float const* d, fcomplex_t const* e,
                      fcomplex_t* b, integer_t const* ldb, integer_t* info);
  void LAPACK_ZPTTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      double const* d, dcomplex_t const* e,
                      dcomplex_t* b, integer_t const* ldb, integer_t* info);



  /* symmetric/Hermitian indefinite and complex symmetric */

  void LAPACK_SSYSV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     float* a, integer_t const* lda, integer_t* ipiv,
                     float* b, integer_t const* ldb,
                     float* w, integer_t const* lw, integer_t* info);
  void LAPACK_DSYSV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     double* a, integer_t const* lda, integer_t* ipiv,
                     double* b, integer_t const* ldb,
                     double* w, integer_t const* lw, integer_t* info);
  void LAPACK_CSYSV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     fcomplex_t* a, integer_t const* lda, integer_t* ipiv,
                     fcomplex_t* b, integer_t const* ldb,
                     fcomplex_t* w, integer_t const* lw, integer_t* info);
  void LAPACK_ZSYSV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     dcomplex_t* a, integer_t const* lda, integer_t* ipiv,
                     dcomplex_t* b, integer_t const* ldb,
                     dcomplex_t* w, integer_t const* lw, integer_t* info);

  void LAPACK_CHESV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     fcomplex_t* a, integer_t const* lda, integer_t* ipiv,
                     fcomplex_t* b, integer_t const* ldb,
                     fcomplex_t* w, integer_t const* lw, integer_t* info);
  void LAPACK_ZHESV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     dcomplex_t* a, integer_t const* lda, integer_t* ipiv,
                     dcomplex_t* b, integer_t const* ldb,
                     dcomplex_t* w, integer_t const* lw, integer_t* info);

  void LAPACK_SSYTRF (char const* uplo, integer_t const* n,
                      float* a, integer_t const* lda, integer_t* ipiv,
                      float* w, integer_t const* lw, integer_t* info);
  void LAPACK_DSYTRF (char const* uplo, integer_t const* n,
                      double* a, integer_t const* lda, integer_t* ipiv,
                      double* w, integer_t const* lw, integer_t* info);
  void LAPACK_CSYTRF (char const* uplo, integer_t const* n,
                      fcomplex_t* a, integer_t const* lda, integer_t* ipiv,
                      fcomplex_t* w, integer_t const* lw, integer_t* info);
  void LAPACK_ZSYTRF (char const* uplo, integer_t const* n,
                      dcomplex_t* a, integer_t const* lda, integer_t* ipiv,
                      dcomplex_t* w, integer_t const* lw, integer_t* info);

  void LAPACK_CHETRF (char const* uplo, integer_t const* n,
                      fcomplex_t* a, integer_t const* lda, integer_t* ipiv,
                      fcomplex_t* w, integer_t const* lw, integer_t* info);
  void LAPACK_ZHETRF (char const* uplo, integer_t const* n,
                      dcomplex_t* a, integer_t const* lda, integer_t* ipiv,
                      dcomplex_t* w, integer_t const* lw, integer_t* info);

  void LAPACK_SSYTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      float const* a, integer_t const* lda, integer_t const* ipiv,
                      float* b, integer_t const* ldb, integer_t* info);
  void LAPACK_DSYTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      double const* a, integer_t const* lda, integer_t const* ipiv,
                      double* b, integer_t const* ldb, integer_t* info);
  void LAPACK_CSYTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      fcomplex_t const* a, integer_t const* lda, integer_t const* ipiv,
                      fcomplex_t* b, integer_t const* ldb, integer_t* info);
  void LAPACK_ZSYTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      dcomplex_t const* a, integer_t const* lda, integer_t const* ipiv,
                      dcomplex_t* b, integer_t const* ldb, integer_t* info);

  void LAPACK_SSYTRI (char const* uplo, integer_t const* n, float* a,
                      integer_t const* lda, integer_t const* ipiv, float* work,
                      integer_t* info);
  void LAPACK_DSYTRI (char const* uplo, integer_t const* n, double* a,
                      integer_t const* lda, integer_t const* ipiv, double* work,
                      integer_t* info);
  void LAPACK_CSYTRI (char const* uplo, integer_t const* n, fcomplex_t* a,
                      integer_t const* lda, integer_t const* ipiv, fcomplex_t* work,
                      integer_t* info);
  void LAPACK_ZSYTRI (char const* uplo, integer_t const* n, dcomplex_t* a,
                      integer_t const* lda, integer_t const* ipiv, dcomplex_t* work,
                      integer_t* info);
 
  void LAPACK_CHETRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      fcomplex_t const* a, integer_t const* lda, integer_t const* ipiv,
                      fcomplex_t* b, integer_t const* ldb, integer_t* info);
  void LAPACK_ZHETRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      dcomplex_t const* a, integer_t const* lda, integer_t const* ipiv, 
                      dcomplex_t* b, integer_t const* ldb, integer_t* info);


  /* symmetric/Hermitian indefinite and complex symmetric in packed storage */

  void LAPACK_SSPSV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     float* ap, integer_t* ipiv,
                     float* b, integer_t const* ldb, integer_t* info);
  void LAPACK_DSPSV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     double* ap, integer_t* ipiv,
                     double* b, integer_t const* ldb, integer_t* info);
  void LAPACK_CSPSV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     fcomplex_t* ap, integer_t* ipiv,
                     fcomplex_t* b, integer_t const* ldb, integer_t* info);
  void LAPACK_ZSPSV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     dcomplex_t* ap, integer_t* ipiv,
                     dcomplex_t* b, integer_t const* ldb, integer_t* info);

  void LAPACK_CHPSV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     fcomplex_t* ap, integer_t* ipiv,
                     fcomplex_t* b, integer_t const* ldb, integer_t* info);
  void LAPACK_ZHPSV (char const* uplo, integer_t const* n, integer_t const* nrhs,
                     dcomplex_t* ap, integer_t* ipiv,
                     dcomplex_t* b, integer_t const* ldb, integer_t* info);

  void LAPACK_SSPTRF (char const* uplo, integer_t const* n,
                      float* ap, integer_t* ipiv, integer_t* info);
  void LAPACK_DSPTRF (char const* uplo, integer_t const* n,
                      double* ap, integer_t* ipiv, integer_t* info);
  void LAPACK_CSPTRF (char const* uplo, integer_t const* n,
                      fcomplex_t* ap, integer_t* ipiv, integer_t* info);
  void LAPACK_ZSPTRF (char const* uplo, integer_t const* n,
                      dcomplex_t* ap, integer_t* ipiv, integer_t* info);

  void LAPACK_CHPTRF (char const* uplo, integer_t const* n,
                      fcomplex_t* ap, integer_t* ipiv, integer_t* info);
  void LAPACK_ZHPTRF (char const* uplo, integer_t const* n,
                      dcomplex_t* ap, integer_t* ipiv, integer_t* info);

  void LAPACK_SSPTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      float const* ap, integer_t const* ipiv,
                      float* b, integer_t const* ldb, integer_t* info);
  void LAPACK_DSPTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      double const* ap, integer_t const* ipiv,
                      double* b, integer_t const* ldb, integer_t* info);
  void LAPACK_CSPTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      fcomplex_t const* ap, integer_t const* ipiv,
                      fcomplex_t* b, integer_t const* ldb, integer_t* info);
  void LAPACK_ZSPTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      dcomplex_t const* ap, integer_t const* ipiv,
                      dcomplex_t* b, integer_t const* ldb, integer_t* info);

  void LAPACK_SSPTRI (char const* uplo, integer_t const* n, float const* ap,
                      integer_t const* ipiv, float* work, integer_t* info);
  void LAPACK_DSPTRI (char const* uplo, integer_t const* n, double const* ap,
                      integer_t const* ipiv, double* work, integer_t* info);
  void LAPACK_CSPTRI (char const* uplo, integer_t const* n, fcomplex_t const* ap,
                      integer_t const* ipiv, fcomplex_t* work, integer_t* info);
  void LAPACK_ZSPTRI (char const* uplo, integer_t const* n, dcomplex_t const* ap,
                      integer_t const* ipiv, dcomplex_t* work, integer_t* info);

  void LAPACK_CHPTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      fcomplex_t const* ap, integer_t const* ipiv,
                      fcomplex_t* b, integer_t const* ldb, integer_t* info);
  void LAPACK_ZHPTRS (char const* uplo, integer_t const* n, integer_t const* nrhs,
                      dcomplex_t const* ap, integer_t const* ipiv,
                      dcomplex_t* b, integer_t const* ldb, integer_t* info);

  /* banded */

  void LAPACK_SGBTRF (integer_t const* n, integer_t const* m, integer_t const* kl, integer_t const* ku,
                      float* ab, integer_t const* ldab, integer_t* ipiv, integer_t* info);
  void LAPACK_DGBTRF (integer_t const* n, integer_t const* m, integer_t const* kl, integer_t const* ku,
                      double* ab, integer_t const* ldab, integer_t* ipiv, integer_t* info);
  void LAPACK_CGBTRF (integer_t const* n, integer_t const* m, integer_t const* kl, integer_t const* ku,
                      fcomplex_t* ab, integer_t const* ldab, integer_t* ipiv, integer_t* info);
  void LAPACK_ZGBTRF (integer_t const* n, integer_t const* m, integer_t const* kl, integer_t const* ku,
                      dcomplex_t* ab, integer_t const* ldab, integer_t* ipiv, integer_t* info);

  void LAPACK_SGBTRS (char const* trans, integer_t const* n, integer_t const* kl, integer_t const* ku, integer_t const* nrhs,
                      float const* ab, integer_t const* ldab, integer_t const* ipiv,
                      float* b, integer_t const* ldb, integer_t* info);
  void LAPACK_DGBTRS (char const* trans, integer_t const* n, integer_t const* kl, integer_t const* ku, integer_t const* nrhs,
                      double const* ab, integer_t const* ldab, integer_t const* ipiv,
                      double* b, integer_t const* ldb, integer_t* info);
  void LAPACK_CGBTRS (char const* trans, integer_t const* n, integer_t const* kl, integer_t const* ku, integer_t const* nrhs,
                      fcomplex_t const* ab, integer_t const* ldab, integer_t const* ipiv,
                      fcomplex_t* b, integer_t const* ldb, integer_t* info);
  void LAPACK_ZGBTRS (char const* trans, integer_t const* n, integer_t const* kl, integer_t const* ku, integer_t const* nrhs,
                      dcomplex_t const* ab, integer_t const* ldab, integer_t const* ipiv,
                      dcomplex_t* b, integer_t const* ldb, integer_t* info);


  /**********************************************************************/
  /*                         eigenproblems                              */
  /**********************************************************************/

  void LAPACK_SGEES (const char* jobvs, const char* sort, logical_t* select, const integer_t* n,
                     float* a, const integer_t * lda, const integer_t* sdim, float* wr, float* wi,
                     float* vs, const integer_t * ldvs, float* work, const integer_t * lwork,
                     bool* bwork, integer_t* info);

  void LAPACK_DGEES (const char* jobvs, const char* sort, logical_t* select, const integer_t* n,
                     double* a, const integer_t * lda, const integer_t* sdim, double* wr, double* wi,
                     double* vs, const integer_t * ldvs, double* work, const integer_t * lwork,
                     bool* bwork, integer_t* info);

  void LAPACK_CGEES( const char* jobvs, const char* sort, logical_t* select, const integer_t *n,
                     fcomplex_t* a, const integer_t * lda, integer_t * sdim, fcomplex_t* w, fcomplex_t* vs,
                     const integer_t * ldvs, fcomplex_t* work, const integer_t * lwork, float* rwork,
                     bool* bwork, integer_t* info );

  void LAPACK_ZGEES( const char* jobvs, const char* sort, const logical_t* select, const integer_t *n,
                     dcomplex_t* a, const integer_t * lda, integer_t * sdim, dcomplex_t* w, dcomplex_t* vs,
                     const integer_t * ldvs, dcomplex_t* work, const integer_t * lwork, double* rwork,
                     bool* bwork, integer_t* info );


  void LAPACK_SGEEV( const char* jobvl, const char* jobvr, const integer_t* n, float* a,
                    const integer_t* lda, float* wr, float* wi, float* vl, const integer_t* ldvl,
                    float* vr, const integer_t* ldvr, float* work, const integer_t* lwork, integer_t* info );

  void LAPACK_DGEEV( const char* jobvl, const char* jobvr, const integer_t* n, double* a,
                    const integer_t* lda, double* wr, double* wi, double* vl, const integer_t* ldvl,
                    double* vr, const integer_t* ldvr, double* work, const integer_t* lwork, integer_t* info );

  void LAPACK_CGEEV( const char* jobvl, const char* jobvr, const integer_t* n, fcomplex_t* a,
                    const integer_t* lda, fcomplex_t* w, fcomplex_t* vl, const integer_t* ldvl,
                    fcomplex_t* vr, const integer_t* ldvr, fcomplex_t* work, const integer_t* lwork,
                    float* rwork, integer_t* info );

  void LAPACK_ZGEEV( const char* jobvl, const char* jobvr, const integer_t* n, dcomplex_t* a,
                    const integer_t* lda, dcomplex_t* w, dcomplex_t* vl, const integer_t* ldvl,
                    dcomplex_t* vr, const integer_t* ldvr, dcomplex_t* work, const integer_t* lwork,
                    double* rwork, integer_t* info );


  void LAPACK_SSYEV( const char* jobz, const char* uplo, const integer_t *n,
                     float* a, const integer_t * lda, float* w,
                     float* work, const integer_t * lwork, integer_t* info );

  void LAPACK_DSYEV( const char* jobz, const char* uplo, const integer_t *n,
                     double* a, const integer_t * lda, double* w,
                     double* work, const integer_t * lwork, integer_t* info );

  void LAPACK_CHEEV( const char* jobz, const char* uplo, const integer_t *n,
                     fcomplex_t* a, const integer_t * lda, float* w,
                     fcomplex_t* work, const integer_t * lwork, float* rwork,
                     integer_t* info );

  void LAPACK_ZHEEV( const char* jobz, const char* uplo, const integer_t *n,
                     dcomplex_t* a, const integer_t * lda, double* w,
                     dcomplex_t* work, const integer_t * lwork, double* rwork,
                     integer_t* info );


  void LAPACK_SSYEVD( const char* jobz, const char* uplo, const integer_t* n,
                      float* a, const integer_t* lda, float* w,
                      float* work, const integer_t* lwork,
                      integer_t* iwork, const integer_t* liwork, integer_t* info);

  void LAPACK_DSYEVD( const char* jobz, const char* uplo, const integer_t* n,
                      double* a, const integer_t* lda, double* w,
                      double* work, const integer_t* lwork,
                      integer_t* iwork, const integer_t* liwork, integer_t* info);

  void LAPACK_CHEEVD( const char* jobz, const char* uplo, const integer_t* n,
                      fcomplex_t* a, const integer_t* lda, float* w,
                      fcomplex_t* work, const integer_t* lwork, float* rwork, const integer_t* lrwork,
                      integer_t* iwork, const integer_t* liwork, integer_t* info);

  void LAPACK_ZHEEVD( const char* jobz, const char* uplo, const integer_t* n,
                      dcomplex_t* a, const integer_t* lda, double* w,
                      dcomplex_t* work, const integer_t* lwork, double* rwork, const integer_t* lrwork,
                      integer_t* iwork, const integer_t* liwork, integer_t* info);


  void LAPACK_SSYEVX( const char* jobz, const char* range, const char* uplo, const integer_t* n,
                      float* a, const integer_t* lda, const float* vl, const float* vu, const integer_t* il, const integer_t* iu,
                      const float* abstol, integer_t* m, float* w, float* z, const integer_t* ldz,
                      float* work, const integer_t* lwork,
                      integer_t* iwork, integer_t* ifail, integer_t* info);

  void LAPACK_DSYEVX( const char* jobz, const char* range, const char* uplo, const integer_t* n,
                      double* a, const integer_t* lda, const double* vl, const double* vu, const integer_t* il, const integer_t* iu,
                      const double* abstol, integer_t* m, double* w, double* z, const integer_t* ldz,
                      double* work, const integer_t* lwork,
                      integer_t* iwork, integer_t* ifail, integer_t* info);

  void LAPACK_CHEEVX( const char* jobz, const char* range, const char* uplo, const integer_t* n,
                      fcomplex_t* a, const integer_t* lda, const float* vl, const float* vu, const integer_t* il, const integer_t* iu,
                      const float* abstol, integer_t* m, float* w, fcomplex_t* z, const integer_t* ldz,
                      fcomplex_t* work, const integer_t* lwork, float* rwork,
                      integer_t* iwork, integer_t* ifail, integer_t* info);

  void LAPACK_ZHEEVX( const char* jobz, const char* range, const char* uplo, const integer_t* n,
                      dcomplex_t* a, const integer_t* lda, const double* vl, const double* vu, const integer_t* il, const integer_t* iu,
                      const double* abstol, integer_t* m, double* w, dcomplex_t* z, const integer_t* ldz,
                      dcomplex_t* work, const integer_t* lwork, double* rwork,
                      integer_t* iwork, integer_t* ifail, integer_t* info);


  void LAPACK_CTREVC( const char* side, const char* howmny, const logical_t* select, const integer_t *n,
                     fcomplex_t* t, const integer_t * ldt, fcomplex_t* vl, const integer_t* ldvl,
                     fcomplex_t* vr, const integer_t * ldvr, const integer_t * nm, integer_t* m, fcomplex_t* work,
                     float* rwork, integer_t* info );

  void LAPACK_ZTREVC( const char* side, const char* howmny, const logical_t* select, const integer_t *n,
                     dcomplex_t* t, const integer_t * ldt, dcomplex_t* vl, const integer_t* ldvl,
                     dcomplex_t* vr, const integer_t * ldvr, const integer_t * nm, integer_t* m, dcomplex_t* work,
                     double* rwork, integer_t* info );

  void LAPACK_STREVC( const char* side, const char* howmny, const logical_t* select, const integer_t *n,
                     float* t, const integer_t * ldt, float* vl, const integer_t* ldvl,
                     float* vr, const integer_t * ldvr, const integer_t * nm, integer_t* m, float* work,
                     integer_t* info );

  void LAPACK_DTREVC( const char* side, const char* howmny, const logical_t* select, const integer_t *n,
                     double* t, const integer_t * ldt, double* vl, const integer_t* ldvl,
                     double* vr, const integer_t * ldvr, const integer_t * nm, integer_t* m, double* work,
                     integer_t* info );


  void LAPACK_STREXC( const char* compq, const integer_t *n,
                     float* t, const integer_t * ldt, float* q, const integer_t* ldq,
                     integer_t* ifst, const integer_t * ilst, float* work, integer_t* info );

  void LAPACK_DTREXC( const char* compq, const integer_t *n,
                     double* t, const integer_t * ldt, double* q, const integer_t* ldq,
                     integer_t* ifst, const integer_t * ilst, double* work, integer_t* info );

  void LAPACK_CTREXC( const char* compq, const integer_t *n,
                     fcomplex_t* t, const integer_t * ldt, fcomplex_t* q, const integer_t* ldq,
                     integer_t* ifst, const integer_t * ilst, integer_t* info );

  void LAPACK_ZTREXC( const char* compq, const integer_t *n,
                     dcomplex_t* t, const integer_t * ldt, dcomplex_t* q, const integer_t* ldq,
                     integer_t* ifst, const integer_t * ilst, integer_t* info );


  /* Hessenberg matrices */

  void LAPACK_SHSEQR( const char* JOB, const char* COMPZ, const integer_t* N, const integer_t* ILO, const integer_t* IHI, float* H,
                      const integer_t* LDH, float* WR, float* WI, float* Z, integer_t const* LDZ,
                      float* WORK, const integer_t* LWORK, integer_t* INFO ) ;

  void LAPACK_CHSEQR( const char* JOB, const char* COMPZ, const integer_t* N, const integer_t* ILO, const integer_t* IHI, fcomplex_t* H,
                      const integer_t* LDH, fcomplex_t* W, fcomplex_t* Z, integer_t const* LDZ,
                      fcomplex_t* WORK, const integer_t* LWORK, integer_t* INFO ) ;

  void LAPACK_DHSEQR( const char* JOB, const char* COMPZ, const integer_t* N, const integer_t* ILO, const integer_t* IHI, double* H,
                      const integer_t* LDH, double* WR, double* WI, double* Z, integer_t const* LDZ,
                      double* WORK, const integer_t* LWORK, integer_t* INFO ) ;

  void LAPACK_ZHSEQR( const char* JOB, const char* COMPZ, const integer_t* N, const integer_t* ILO, const integer_t* IHI, dcomplex_t* H,
                      const integer_t* LDH, dcomplex_t* W, dcomplex_t* Z, integer_t const* LDZ,
                      dcomplex_t* WORK, const integer_t* LWORK, integer_t* INFO ) ;

  /* Hermitian tridiagonal matrices */
  
  void LAPACK_SSTEQR( char const* compz, integer_t const* n, float* d, float* E, float* z, integer_t const* ldz, float* work, integer_t* info ) ;
  void LAPACK_DSTEQR( char const* compz, integer_t const* n, double* d, double* E, double* z, integer_t const* ldz, double* work, integer_t* info ) ;

  /* Hermitian banded matrices */
  
  void LAPACK_SSBEV( char const* jobz, char const* uplo, integer_t const* n,
                     integer_t const* kd, float* ab, integer_t const* ldab, float* w,
                     float* z, integer_t const* ldz, float* work, integer_t* info );

  void LAPACK_DSBEV( char const* jobz, char const* uplo, integer_t const* n,
                     integer_t const* kd, double* ab, integer_t const* ldab, double* w,
                     double* z, integer_t const* ldz, double* work, integer_t* info );

  void LAPACK_CHBEV( char const* jobz, char const* uplo, integer_t const* n,
                     integer_t const* kd, fcomplex_t* ab, integer_t const* ldab, float* w,
                     fcomplex_t* z, integer_t const* ldz, fcomplex_t* work,
                     float* rwork, integer_t* info );

  void LAPACK_ZHBEV( char const* jobz, char const* uplo, integer_t const* n,
                     integer_t const* kd, dcomplex_t* ab, integer_t const* ldab, double* w,
                     dcomplex_t* z, integer_t const* ldz, dcomplex_t* work,
                     double* rwork, integer_t* info );


  void LAPACK_SSBEVX( char const* jobz, char const* range, char const* uplo, integer_t const* n,
                      integer_t const* kd, float* ab, integer_t const* ldab, float* q, integer_t const* ldq,
                      const float* vl, const float* vu, const integer_t* il, const integer_t* iu,
                      const float* abstol, integer_t* m,
                      float* w, float* z, integer_t const* ldz, float* work,
                      integer_t* iwork, integer_t* ifail, integer_t* info );

  void LAPACK_DSBEVX( char const* jobz, char const* range, char const* uplo, integer_t const* n,
                      integer_t const* kd, double* ab, integer_t const* ldab, double* q, integer_t const* ldq,
                      const double* vl, const double* vu, const integer_t* il, const integer_t* iu,
                      const double* abstol, integer_t* m,
                      double* w, double* z, integer_t const* ldz, double* work,
                      integer_t* iwork, integer_t* ifail, integer_t* info );

  void LAPACK_CHBEVX( char const* jobz, char const* range, char const* uplo, integer_t const* n,
                      integer_t const* kd, fcomplex_t* ab, integer_t const* ldab, fcomplex_t* q, integer_t const* ldq,
                      const float* vl, const float* vu, const integer_t* il, const integer_t* iu,
                      const float* abstol, integer_t* m,
                      float* w, fcomplex_t* z, integer_t const* ldz, fcomplex_t* work, float* rwork,
                      integer_t* iwork, integer_t* ifail, integer_t* info );

  void LAPACK_ZHBEVX( char const* jobz, char const* range, char const* uplo, integer_t const* n,
                      integer_t const* kd, dcomplex_t* ab, integer_t const* ldab, dcomplex_t* q, integer_t const* ldq,
                      const double* vl, const double* vu, const integer_t* il, const integer_t* iu,
                      const double* abstol, integer_t* m,
                      double* w, dcomplex_t* z, integer_t const* ldz, dcomplex_t* work, double* rwork,
                      integer_t* iwork, integer_t* ifail, integer_t* info );


  /*********************************************************************/
  /*       Auxiliary routines for eigenvalue problems                  */
  /*********************************************************************/

  void LAPACK_SSYTRD( char const* uplo, integer_t const* n, float* a, integer_t const* lda, float* d,
                      float* e, float* tau, float* work, integer_t const* lwork, integer_t* INFO ) ;

  void LAPACK_DSYTRD( char const* uplo, integer_t const* n, double* a, integer_t const* lda, double* d,
                      double* e, double* tau, double* work, integer_t const* lwork, integer_t* INFO ) ;


  /**********************************************************************/
  /*               generalized eigenvalue / eigenvector                 */
  /**********************************************************************/

   void LAPACK_SSYGV(integer_t const* itype, char const* jobz, char const* uplo, integer_t const* n,
                     float* a, integer_t const* lda, float* b, integer_t const* ldb,
                     float* w, float* work, integer_t const* lwork, integer_t* info);

   void LAPACK_DSYGV(integer_t const* itype, char const* jobz, char const* uplo, integer_t const* n,
                     double* a, integer_t const* lda, double* b, integer_t const* ldb,
                     double* w, double* work, integer_t const* lwork, integer_t* info);

   void LAPACK_CHEGV(integer_t const* itype, char const* jobz, char const* uplo, integer_t const* n,
                     fcomplex_t* a, integer_t const* lda, fcomplex_t* b, integer_t const* ldb,
                     float* w, fcomplex_t* work, integer_t const* lwork, float* rwork, integer_t* info);

   void LAPACK_ZHEGV(integer_t const *itype, char const* jobz, char const *uplo, integer_t const * n,
                     dcomplex_t *a, integer_t const *lda, dcomplex_t *b, integer_t const *ldb,
                     double *w, dcomplex_t *work, integer_t const *lwork, double* rwork, integer_t* info);


  /*********************************************************************/
  /*                             SVD                                   */
  /*********************************************************************/

  void LAPACK_SGESVD (char const* jobu, char const* jobvt,
                      integer_t const* m, integer_t const* n, float* a, integer_t const* lda,
                      float* s, float* u, integer_t const* ldu,
                      float* vt, integer_t const* ldvt,
                      float* work, integer_t const* lwork, integer_t* info);
  void LAPACK_DGESVD (char const* jobu, char const* jobvt,
                      integer_t const* m, integer_t const* n, double* a, integer_t const* lda,
                      double* s, double* u, integer_t const* ldu,
                      double* vt, integer_t const* ldvt,
                      double* work, integer_t const* lwork, integer_t* info);
  void LAPACK_CGESVD (char const* jobu, char const* jobvt,
                      integer_t const* m, integer_t const* n,
                      fcomplex_t* a, integer_t const* lda,
                      float* s, fcomplex_t* u, integer_t const* ldu,
                      fcomplex_t* vt, integer_t const* ldvt,
                      fcomplex_t* work, integer_t const* lwork,
                      float* rwork, integer_t* info);
  void LAPACK_ZGESVD (char const* jobu, char const* jobvt,
                      integer_t const* m, integer_t const* n,
                      dcomplex_t* a, integer_t const* lda,
                      double* s, dcomplex_t* u, integer_t const* ldu,
                      dcomplex_t* vt, integer_t const* ldvt,
                      dcomplex_t* work, integer_t const* lwork,
                      double* rwork, integer_t* info);

  void LAPACK_SGESDD (char const* jobz, integer_t const* m, integer_t const* n,
                      float* a, integer_t const* lda,
                      float* s, float* u, integer_t const* ldu,
                      float* vt, integer_t const* ldvt,
                      float* work, integer_t const* lwork, integer_t* iwork, integer_t* info);
  void LAPACK_DGESDD (char const* jobz, integer_t const* m, integer_t const* n,
                      double* a, integer_t const* lda,
                      double* s, double* u, integer_t const* ldu,
                      double* vt, integer_t const* ldvt,
                      double* work, integer_t const* lwork, integer_t* iwork, integer_t* info);
  void LAPACK_CGESDD (char const* jobz, integer_t const* m, integer_t const* n,
                      fcomplex_t* a, integer_t const* lda,
                      float* s, fcomplex_t* u, integer_t const* ldu,
                      fcomplex_t* vt, integer_t const* ldvt,
                      fcomplex_t* work, integer_t const* lwork,
                      float* rwork, integer_t* iwork, integer_t* info);
  void LAPACK_ZGESDD (char const* jobz, integer_t const* m, integer_t const* n,
                      dcomplex_t* a, integer_t const* lda,
                      double* s, dcomplex_t* u, integer_t const* ldu,
                      dcomplex_t* vt, integer_t const* ldvt,
                      dcomplex_t* work, integer_t const* lwork,
                      double* rwork, integer_t* iwork, integer_t* info);


  /*********************************************************************/
  /*                    QR factorization                               */
  /*********************************************************************/

  void LAPACK_SGEQRF( const integer_t* m, const integer_t* n, float* a, const integer_t* lda,
                      float* tau, float* work, const integer_t* lwork, integer_t* info );
  void LAPACK_DGEQRF( const integer_t* m, const integer_t* n, double* a, const integer_t* lda,
                      double* tau, double* work, const integer_t* lwork, integer_t* info );
  void LAPACK_CGEQRF( const integer_t* m, const integer_t* n, fcomplex_t* a, const integer_t* lda,
                      fcomplex_t* tau, fcomplex_t* work, const integer_t* lwork, integer_t* info );
  void LAPACK_ZGEQRF( const integer_t* m, const integer_t* n, dcomplex_t* a, const integer_t* lda,
                      dcomplex_t* tau, dcomplex_t* work, const integer_t* lwork, integer_t* info );



  void LAPACK_SORMQR( const char* side, const char* trans, const integer_t* m,
                      const integer_t* n, const integer_t* k, const float* a,
                      const integer_t* lda, const float* tau,
                      float* c, const integer_t* ldc, float* work,
                      const integer_t* lwork, integer_t* info );
  void LAPACK_DORMQR( const char* side, const char* trans, const integer_t* m,
                      const integer_t* n, const integer_t* k, const double* a,
                      const integer_t* lda, const double* tau,
                      double* c, const integer_t* ldc, double* work,
                      const integer_t* lwork, integer_t* info );
  void LAPACK_CUNMQR( const char* side, const char* trans, const integer_t* m,
                      const integer_t* n, const integer_t* k, const fcomplex_t* a,
                      const integer_t* lda, const fcomplex_t* tau,
                      fcomplex_t* c, const integer_t* ldc, fcomplex_t* work,
                      const integer_t* lwork, integer_t* info );
  void LAPACK_ZUNMQR( const char* side, const char* trans, const integer_t* m,
                      const integer_t* n, const integer_t* k, const dcomplex_t* a,
                      const integer_t* lda, const dcomplex_t* tau,
                      dcomplex_t* c, const integer_t* ldc, dcomplex_t* work,
                      const integer_t* lwork, integer_t* info );

  void LAPACK_SORGQR( const integer_t* m, const integer_t* n, const integer_t* k,
                      float* a, const integer_t* lda, float* tau,
                      float* work, const integer_t* lwork, const integer_t* info);
  void LAPACK_DORGQR( const integer_t* m, const integer_t* n, const integer_t* k,
                      double* a, const integer_t* lda, double* tau,
                      double* work, const integer_t* lwork, const integer_t* info);
  void LAPACK_CUNGQR( const integer_t* m, const integer_t* n, const integer_t* k,
                      fcomplex_t* a, const integer_t* lda, fcomplex_t* tau,
                      fcomplex_t* work, const integer_t* lwork, const integer_t* info);
  void LAPACK_ZUNGQR( const integer_t* m, const integer_t* n, const integer_t* k,
                      dcomplex_t* a, const integer_t* lda, dcomplex_t* tau,
                      dcomplex_t* work, const integer_t* lwork, const integer_t* info);


  /********************************************************************/
  /*                          Least Squares                           */
  /********************************************************************/

  void LAPACK_SGELS(const char* trans, const integer_t* m, const integer_t* n,
                    const integer_t *nrhs, float* a, const integer_t* lda,
                    float* b, const integer_t* ldb, float* work,
                    const integer_t* lwork, integer_t* info);
  void LAPACK_DGELS(const char* trans, const integer_t* m, const integer_t* n,
                    const integer_t *nrhs, double* a, const integer_t* lda,
                    double* b, const integer_t* ldb, double* work,
                    const integer_t* lwork, integer_t* info);
  void LAPACK_CGELS(const char* trans, const integer_t* m, const integer_t* n,
                    const integer_t *nrhs, fcomplex_t* a, const integer_t* lda,
                    fcomplex_t* b, const integer_t* ldb, fcomplex_t* work,
                    const integer_t* lwork, integer_t* info);
  void LAPACK_ZGELS(const char* trans, const integer_t* m, const integer_t* n,
                    const integer_t *nrhs, dcomplex_t* a, const integer_t* lda,
                    dcomplex_t* b, const integer_t* ldb, dcomplex_t* work,
                    const integer_t* lwork, integer_t* info);


  void LAPACK_SGELSS(const integer_t *m, const integer_t *n, const integer_t *nrhs,
                     float *a, const integer_t *lda, float *b, const integer_t *ldb,
                     float *s, const float *rcond, integer_t *rank, float *work,
                     const integer_t *lwork, integer_t *info);
  void LAPACK_DGELSS(const integer_t *m, const integer_t *n, const integer_t *nrhs,
                     double *a, const integer_t *lda, double *b, const integer_t *ldb,
                     double *s, const double *rcond, integer_t *rank, double *work,
                     const integer_t *lwork, integer_t *info);
  void LAPACK_CGELSS(const integer_t *m, const integer_t *n, const integer_t *nrhs,
                     fcomplex_t *a, const integer_t *lda, fcomplex_t *b, const integer_t *ldb,
                     float *s, const float *rcond, integer_t *rank, fcomplex_t *work,
                     const integer_t *lwork, float *rwork, integer_t *info);
  void LAPACK_ZGELSS(const integer_t *m, const integer_t *n, const integer_t *nrhs,
                     dcomplex_t *a, const integer_t *lda, dcomplex_t *b, const integer_t *ldb,
                     double *s, const double *rcond, integer_t *rank, dcomplex_t *work,
                     const integer_t *lwork, double *rwork, integer_t *info);


  void LAPACK_SGELSD(const integer_t *m, const integer_t *n, const integer_t *nrhs,
                     float *a, const integer_t *lda, float *b, const integer_t *ldb,
                     float *s, const float *rcond, integer_t *rank, float *work,
                     const integer_t *lwork, integer_t *iwork, integer_t *info);
  void LAPACK_DGELSD(const integer_t *m, const integer_t *n, const integer_t *nrhs,
                     double *a, const integer_t *lda, double *b, const integer_t *ldb,
                     double *s, const double *rcond, integer_t *rank, double *work,
                     const integer_t *lwork, integer_t *iwork, integer_t *info);
  void LAPACK_CGELSD(const integer_t *m, const integer_t *n, const integer_t *nrhs,
                     fcomplex_t *a, const integer_t *lda, fcomplex_t *b, const integer_t *ldb,
                     float *s, const float *rcond, integer_t *rank, fcomplex_t *work,
                     const integer_t *lwork, float *rwork, integer_t *iwork, integer_t *info);
  void LAPACK_ZGELSD(const integer_t *m, const integer_t *n, const integer_t *nrhs,
                     dcomplex_t *a, const integer_t *lda, dcomplex_t *b, const integer_t *ldb,
                     double *s, const double *rcond, integer_t *rank, dcomplex_t *work,
                     const integer_t *lwork, double *rwork, integer_t *iwork, integer_t *info);



  /********************************************************************/
  /*                          auxiliary                               */
  /********************************************************************/

  integer_t LAPACK_ILAENV (integer_t const* ispec, const char* name, const char* opt,
                     integer_t const* n1, integer_t const* n2, integer_t const* n3, 
                     integer_t const* n4, integer_t, integer_t);

}

#endif
