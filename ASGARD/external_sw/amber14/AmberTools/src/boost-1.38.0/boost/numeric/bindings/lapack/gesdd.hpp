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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_GESDD_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_GESDD_HPP

#include <boost/numeric/bindings/traits/type_traits.hpp>
#include <boost/numeric/bindings/traits/traits.hpp>
#include <boost/numeric/bindings/lapack/lapack.h>
#include <boost/numeric/bindings/traits/detail/array.hpp>
#include <boost/numeric/bindings/traits/detail/utils.hpp>

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
#  include <boost/static_assert.hpp>
#  include <boost/type_traits/is_same.hpp>
#endif

#include <cassert>

namespace boost { namespace numeric { namespace bindings {

  namespace lapack {

    ///////////////////////////////////////////////////////////////////
    //
    // singular value decomposition 
    //
    ///////////////////////////////////////////////////////////////////

    /*
     * (divide and conquer driver)
     * gesdd() computes the singular value decomposition (SVD) of
     * M-by-N matrix A, optionally computing the left and/or right
     * singular vectors, by using divide-and-conquer method.
     * The SVD is written
     *
     *     A = U * S * V^T    or    A = U * S * V^H
     *
     * where S is an M-by-N matrix which is zero except for its min(m,n)
     * diagonal elements, U is an M-by-M orthogonal/unitary matrix, and V
     * is an N-by-N orthogonal/unitary matrix. The diagonal elements of S
     * are the singular values of A; they are real and non-negative, and
     * are returnede in descending  order. The first min(m,n) columns of
     * U and V are the left and right singular vectors of A. (Note that
     * the routine returns V^T or V^H, not V.
     */

    namespace detail {

      inline
      void gesdd (char const jobz, integer_t const m, integer_t const n,
                  float* a, integer_t const lda,
                  float* s, float* u, integer_t const ldu,
                  float* vt, integer_t const ldvt,
                  float* work, integer_t const lwork, float* /* dummy */,
                  integer_t* iwork, integer_t* info)
      {
        LAPACK_SGESDD (&jobz, &m, &n, a, &lda, s,
                       u, &ldu, vt, &ldvt, work, &lwork, iwork, info);
      }

      inline
      void gesdd (char const jobz, integer_t const m, integer_t const n,
                  double* a, integer_t const lda,
                  double* s, double* u, integer_t const ldu,
                  double* vt, integer_t const ldvt,
                  double* work, integer_t const lwork, double* /* dummy */,
                  integer_t* iwork, integer_t* info)
      {
        LAPACK_DGESDD (&jobz, &m, &n, a, &lda, s,
                       u, &ldu, vt, &ldvt, work, &lwork, iwork, info);
      }

      inline
      void gesdd (char const jobz, integer_t const m, integer_t const n,
                  traits::complex_f* a, integer_t const lda,
                  float* s, traits::complex_f* u, integer_t const ldu,
                  traits::complex_f* vt, integer_t const ldvt,
                  traits::complex_f* work, integer_t const lwork,
                  float* rwork, integer_t* iwork, integer_t* info)
      {
        LAPACK_CGESDD (&jobz, &m, &n,
                       traits::complex_ptr (a), &lda, s,
                       traits::complex_ptr (u), &ldu,
                       traits::complex_ptr (vt), &ldvt,
                       traits::complex_ptr (work), &lwork,
                       rwork, iwork, info);
      }

      inline
      void gesdd (char const jobz, integer_t const m, integer_t const n,
                  traits::complex_d* a, integer_t const lda,
                  double* s, traits::complex_d* u, integer_t const ldu,
                  traits::complex_d* vt, integer_t const ldvt,
                  traits::complex_d* work, integer_t const lwork,
                  double* rwork, integer_t* iwork, integer_t* info)
      {
        LAPACK_ZGESDD (&jobz, &m, &n,
                       traits::complex_ptr (a), &lda, s,
                       traits::complex_ptr (u), &ldu,
                       traits::complex_ptr (vt), &ldvt,
                       traits::complex_ptr (work), &lwork,
                       rwork, iwork, info);
      }


      inline
      integer_t gesdd_min_work (float, char jobz, integer_t m, integer_t n) {
        integer_t minmn = m < n ? m : n;
        integer_t maxmn = m < n ? n : m;
        integer_t m3 = 3 * minmn;
        integer_t m4 = 4 * minmn;
        integer_t minw;
        if (jobz == 'N') {
          // leading comments:
          //   LWORK >= 3*min(M,N) + max(max(M,N), 6*min(M,N))
          // code:
          //   LWORK >= 3*min(M,N) + max(max(M,N), 7*min(M,N))
          integer_t m7 = 7 * minmn;
          minw = maxmn < m7 ? m7 : maxmn;
          minw += m3;
        }
        if (jobz == 'O') {
          // LWORK >= 3*min(M,N)*min(M,N)
          //          + max(max(M,N), 5*min(M,N)*min(M,N)+4*min(M,N))
          integer_t m5 = 5 * minmn * minmn + m4;
          minw = maxmn < m5 ? m5 : maxmn;
          minw += m3 * minmn;
        }
        if (jobz == 'S' || jobz == 'A') {
          // LWORK >= 3*min(M,N)*min(M,N)
          //          + max(max(M,N), 4*min(M,N)*min(M,N)+4*min(M,N)).
          integer_t m44 = m4 * minmn + m4;
          minw = maxmn < m44 ? m44 : maxmn;
          minw += m3 * minmn;
        }
        return minw; 
      }
      inline
      integer_t gesdd_min_work (double, char jobz, integer_t m, integer_t n) {
        integer_t minmn = m < n ? m : n;
        integer_t maxmn = m < n ? n : m;
        integer_t m3 = 3 * minmn;
        integer_t m4 = 4 * minmn;
        integer_t minw;
        assert(jobz == 'N' || jobz == 'O' || jobz == 'S' || jobz == 'A');
        if (jobz == 'N') {
          // leading comments:
          //   LWORK >= 3*min(M,N) + max(max(M,N), 6*min(M,N))
          // code:
          //   LWORK >= 3*min(M,N) + max(max(M,N), 7*min(M,N))
          integer_t m7 = 7 * minmn;
          minw = maxmn < m7 ? m7 : maxmn;
          minw += m3;
        }
        else if (jobz == 'O') {
          // LWORK >= 3*min(M,N)*min(M,N)
          //          + max(max(M,N), 5*min(M,N)*min(M,N)+4*min(M,N))
          integer_t m5 = 5 * minmn * minmn + m4;
          minw = maxmn < m5 ? m5 : maxmn;
          minw += m3 * minmn;
        }
        else if (jobz == 'S' || jobz == 'A') {
          // LWORK >= 3*min(M,N)*min(M,N)
          //          + max(max(M,N), 4*min(M,N)*min(M,N)+4*min(M,N)).
          integer_t m44 = m4 * minmn + m4;
          minw = maxmn < m44 ? m44 : maxmn;
          minw += m3 * minmn;
        }
        else {
          assert(false);
          return 0 ;
        }
        return minw;
      }
      inline
      integer_t gesdd_min_work (traits::complex_f, char jobz, integer_t m, integer_t n) {
        integer_t minmn = m < n ? m : n;
        integer_t maxmn = m < n ? n : m;
        integer_t m2 = 2 * minmn;
        integer_t minw = m2 + maxmn;
        if (jobz == 'N')
          // LWORK >= 2*min(M,N)+max(M,N)
          ;
        if (jobz == 'O')
          // LWORK >= 2*min(M,N)*min(M,N) + 2*min(M,N) + max(M,N)
          minw += m2 * minmn;
        if (jobz == 'S' || jobz == 'A')
          // LWORK >= min(M,N)*min(M,N) + 2*min(M,N) + max(M,N)
          minw += minmn * minmn;
        return minw;
      }
      inline
      integer_t gesdd_min_work (traits::complex_d, char jobz, integer_t m, integer_t n) {
        integer_t minmn = m < n ? m : n;
        integer_t maxmn = m < n ? n : m;
        integer_t m2 = 2 * minmn;
        integer_t minw = m2 + maxmn;
        if (jobz == 'N')
          // LWORK >= 2*min(M,N)+max(M,N)
          ;
        if (jobz == 'O')
          // LWORK >= 2*min(M,N)*min(M,N) + 2*min(M,N) + max(M,N)
          minw += m2 * minmn;
        if (jobz == 'S' || jobz == 'A')
          // LWORK >= min(M,N)*min(M,N) + 2*min(M,N) + max(M,N)
          minw += minmn * minmn;
        return minw;
      }

      inline
      integer_t gesdd_rwork (float, char, integer_t, integer_t) { return 1; }
      inline
      integer_t gesdd_rwork (double, char, integer_t, integer_t) { return 1; }
      inline
      integer_t gesdd_rwork (traits::complex_f, char jobz, integer_t m, integer_t n) {
        integer_t minmn = m < n ? m : n;
        integer_t minw;
        if (jobz == 'N')
          // LWORK >= 7*min(M,N)
          minw = 7 * minmn;
        else
          // LRWORK >= 5*min(M,N)*min(M,N) + 5*min(M,N)
          minw = 5 * (minmn * minmn + minmn);
        return minw;
      }
      inline
      integer_t gesdd_rwork (traits::complex_d, char jobz, integer_t m, integer_t n) {
        integer_t minmn = m < n ? m : n;
        integer_t minw;
        if (jobz == 'N')
          // LWORK >= 7*min(M,N)
          minw = 7 * minmn;
        else
          // LRWORK >= 5*min(M,N)*min(M,N) + 5*min(M,N)
          minw = 5 * (minmn * minmn + minmn);
        return minw;
      }

      inline
      integer_t gesdd_iwork (integer_t m, integer_t n) {
        integer_t minmn = m < n ? m : n;
        return 8 * minmn;
      }

    } // detail


    template <typename MatrA>
    integer_t gesdd_work (char const q, char const jobz, MatrA const& a)
    {

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrA>::matrix_structure,
        traits::general_t
      >::value));
#endif

#ifdef BOOST_NUMERIC_BINDINGS_LAPACK_2
      assert (q == 'M');
#else
      assert (q == 'M' || q == 'O');
#endif 
      assert (jobz == 'N' || jobz == 'O' || jobz == 'A' || jobz == 'S');

      integer_t const m = traits::matrix_size1 (a);
      integer_t const n = traits::matrix_size2 (a);
      integer_t lw = -13;

#ifndef BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS
      typedef typename traits::matrix_traits<MatrA>::value_type val_t;
#else
      typedef typename MatrA::value_type val_t;
#endif

      if (q == 'M')
        lw = detail::gesdd_min_work (val_t(), jobz, m, n);

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_2
      MatrA& a2 = const_cast<MatrA&> (a);
      if (q == 'O') {
        // traits::detail::array<val_t> w (1);
        val_t w;
        integer_t info;
        detail::gesdd (jobz, m, n,
                       traits::matrix_storage (a2),
                       traits::leading_dimension (a2),
                       0, // traits::vector_storage (s),
                       0, // traits::matrix_storage (u),
                       m, // traits::leading_dimension (u),
                       0, // traits::matrix_storage (vt),
                       n, // traits::leading_dimension (vt),
                       &w, // traits::vector_storage (w),
                       -1, // traits::vector_size (w),
                       0, // traits::vector_storage (rw),
                       0, // traits::vector_storage (iw),
                       &info);
        assert (info == 0);

        lw = traits::detail::to_int (w);
        // // lw = traits::detail::to_int (w[0]);
        /*
         * is there a bug in LAPACK? or in Mandrake's .rpm?
         * if m == 3, n == 4 and jobz == 'N' (real A),
         * gesdd() returns optimal size == 1 while minimum size == 27
         */
        // integer_t lwo = traits::detail::to_int (w);
        // integer_t lwmin = detail::gesdd_min_work (val_t(), jobz, m, n);
        // lw = lwo < lwmin ? lwmin : lwo;
      }
#endif

      return lw;
    }


    template <typename MatrA>
    integer_t gesdd_rwork (char jobz, MatrA const& a) {

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrA>::matrix_structure,
        traits::general_t
      >::value));
#endif

      assert (jobz == 'N' || jobz == 'O' || jobz == 'A' || jobz == 'S');

#ifndef BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS
      typedef typename traits::matrix_traits<MatrA>::value_type val_t;
#else
      typedef typename MatrA::value_type val_t;
#endif

      return detail::gesdd_rwork (val_t(), jobz,
                                  traits::matrix_size1 (a),
                                  traits::matrix_size2 (a));
    }


    template <typename MatrA>
    integer_t gesdd_iwork (MatrA const& a) {

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrA>::matrix_structure,
        traits::general_t
      >::value));
#endif

      return detail::gesdd_iwork (traits::matrix_size1 (a),
                                  traits::matrix_size2 (a));
    }


    template <typename MatrA, typename VecS,
              typename MatrU, typename MatrV, typename VecW, typename VecIW>
    int gesdd (char const jobz, MatrA& a,
               VecS& s, MatrU& u, MatrV& vt, VecW& w, VecIW& iw)
    {

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrA>::matrix_structure,
        traits::general_t
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrU>::matrix_structure,
        traits::general_t
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrV>::matrix_structure,
        traits::general_t
      >::value));

      BOOST_STATIC_ASSERT(
        (boost::is_same<
          typename traits::matrix_traits<MatrA>::value_type, float
        >::value
        ||
        boost::is_same<
          typename traits::matrix_traits<MatrA>::value_type, double
        >::value));
#endif

      integer_t const m = traits::matrix_size1 (a);
      integer_t const n = traits::matrix_size2 (a);
      integer_t const minmn = m < n ? m : n;

      assert (minmn == traits::vector_size (s));
      assert ((jobz == 'N')
              || ((jobz == 'O' || jobz == 'A') && m >= n)
              || ((jobz == 'O' || jobz == 'A')
                  && m < n
                  && m == traits::matrix_size2 (u))
              || (jobz == 'S' && minmn == traits::matrix_size2 (u)));
      assert ((jobz == 'N' && traits::leading_dimension (u) >= 1)
              || (jobz == 'O'
                  && m >= n
                  && traits::leading_dimension (u) >= 1)
              || (jobz == 'O'
                  && m < n
                  && traits::leading_dimension (u) >= m)
              || (jobz == 'A' && traits::leading_dimension (u) >= m)
              || (jobz == 'S' && traits::leading_dimension (u) >= m));
      assert (n == traits::matrix_size2 (vt));
      assert ((jobz == 'N' && traits::leading_dimension (vt) >= 1)
              || (jobz == 'O'
                  && m < n
                  && traits::leading_dimension (vt) >= 1)
              || (jobz == 'O'
                  && m >= n
                  && traits::leading_dimension (vt) >= n)
              || (jobz == 'A' && traits::leading_dimension (vt) >= n)
              || (jobz == 'S' && traits::leading_dimension (vt) >= minmn));
#ifndef BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS
      typedef typename traits::matrix_traits<MatrA>::value_type val_t;
#else
      typedef typename MatrA::value_type val_t;
#endif
      assert (traits::vector_size (w)
              >= detail::gesdd_min_work (val_t(), jobz, m, n));
      assert (traits::vector_size (iw) >= detail::gesdd_iwork (m, n));

      integer_t info;
      detail::gesdd (jobz, m, n,
                     traits::matrix_storage (a),
                     traits::leading_dimension (a),
                     traits::vector_storage (s),
                     traits::matrix_storage (u),
                     traits::leading_dimension (u),
                     traits::matrix_storage (vt),
                     traits::leading_dimension (vt),
                     traits::vector_storage (w),
                     traits::vector_size (w),
                     0, // dummy argument
                     traits::vector_storage (iw),
                     &info);
      return info;
    }


    template <typename MatrA, typename VecS, typename MatrU,
              typename MatrV, typename VecW, typename VecRW, typename VecIW>
    int gesdd (char const jobz, MatrA& a,
               VecS& s, MatrU& u, MatrV& vt, VecW& w, VecRW& rw, VecIW& iw)
    {

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrA>::matrix_structure,
        traits::general_t
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrU>::matrix_structure,
        traits::general_t
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrV>::matrix_structure,
        traits::general_t
      >::value));
#endif

      integer_t const m = traits::matrix_size1 (a);
      integer_t const n = traits::matrix_size2 (a);
      integer_t const minmn = m < n ? m : n;

      assert (minmn == traits::vector_size (s));
      assert ((jobz == 'N')
              || ((jobz == 'O' || jobz == 'A') && m >= n)
              || ((jobz == 'O' || jobz == 'A')
                  && m < n
                  && m == traits::matrix_size2 (u))
              || (jobz == 'S' && minmn == traits::matrix_size2 (u)));
      assert ((jobz == 'N' && traits::leading_dimension (u) >= 1)
              || (jobz == 'O'
                  && m >= n
                  && traits::leading_dimension (u) >= 1)
              || (jobz == 'O'
                  && m < n
                  && traits::leading_dimension (u) >= m)
              || (jobz == 'A' && traits::leading_dimension (u) >= m)
              || (jobz == 'S' && traits::leading_dimension (u) >= m));
      assert (n == traits::matrix_size2 (vt));
      assert ((jobz == 'N' && traits::leading_dimension (vt) >= 1)
              || (jobz == 'O'
                  && m < n
                  && traits::leading_dimension (vt) >= 1)
              || (jobz == 'O'
                  && m >= n
                  && traits::leading_dimension (vt) >= n)
              || (jobz == 'A' && traits::leading_dimension (vt) >= n)
              || (jobz == 'S' && traits::leading_dimension (vt) >= minmn));
#ifndef BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS
      typedef typename traits::matrix_traits<MatrA>::value_type val_t;
#else
      typedef typename MatrA::value_type val_t;
#endif
      assert (traits::vector_size (w)
              >= detail::gesdd_min_work (val_t(), jobz, m, n));
      assert (traits::vector_size (rw)
              >= detail::gesdd_rwork (val_t(), jobz, m, n));
      assert (traits::vector_size (iw) >= detail::gesdd_iwork (m, n));

      integer_t info;
      detail::gesdd (jobz, m, n,
                     traits::matrix_storage (a),
                     traits::leading_dimension (a),
                     traits::vector_storage (s),
                     traits::matrix_storage (u),
                     traits::leading_dimension (u),
                     traits::matrix_storage (vt),
                     traits::leading_dimension (vt),
                     traits::vector_storage (w),
                     traits::vector_size (w),
                     traits::vector_storage (rw),
                     traits::vector_storage (iw),
                     &info);
      return info;
    }


    template <typename MatrA, typename VecS, typename MatrU, typename MatrV>
    int gesdd (char const opt, char const jobz,
               MatrA& a, VecS& s, MatrU& u, MatrV& vt)
    {

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrA>::matrix_structure,
        traits::general_t
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrU>::matrix_structure,
        traits::general_t
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrV>::matrix_structure,
        traits::general_t
      >::value));
#endif

      integer_t const m = traits::matrix_size1 (a);
      integer_t const n = traits::matrix_size2 (a);
#ifndef NDEBUG
      integer_t const minmn = m < n ? m : n;
#endif // NDEBUG

      assert (minmn == traits::vector_size (s));
      assert ((jobz == 'N')
              || ((jobz == 'O' || jobz == 'A') && m >= n)
              || ((jobz == 'O' || jobz == 'A')
                  && m < n
                  && m == traits::matrix_size2 (u))
              || (jobz == 'S' && minmn == traits::matrix_size2 (u)));
      assert ((jobz == 'N' && traits::leading_dimension (u) >= 1)
              || (jobz == 'O'
                  && m >= n
                  && traits::leading_dimension (u) >= 1)
              || (jobz == 'O'
                  && m < n
                  && traits::leading_dimension (u) >= m)
              || (jobz == 'A' && traits::leading_dimension (u) >= m)
              || (jobz == 'S' && traits::leading_dimension (u) >= m));
      assert (n == traits::matrix_size2 (vt));
      assert ((jobz == 'N' && traits::leading_dimension (vt) >= 1)
              || (jobz == 'O'
                  && m < n
                  && traits::leading_dimension (vt) >= 1)
              || (jobz == 'O'
                  && m >= n
                  && traits::leading_dimension (vt) >= n)
              || (jobz == 'A' && traits::leading_dimension (vt) >= n)
              || (jobz == 'S' && traits::leading_dimension (vt) >= minmn));

#ifdef BOOST_NUMERIC_BINDINGS_LAPACK_2
      assert (opt == 'M');
#else
      assert (opt == 'M' || opt == 'O');
#endif

#ifndef BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS
      typedef typename traits::matrix_traits<MatrA>::value_type val_t;
#else 
      typedef typename MatrA::value_type val_t;
#endif 
      typedef typename traits::type_traits<val_t>::real_type real_t;

      integer_t const lw = gesdd_work (opt, jobz, a);
      traits::detail::array<val_t> w (lw);
      if (!w.valid()) return -101;

      integer_t const lrw = gesdd_rwork (jobz, a);
      traits::detail::array<real_t> rw (lrw);
      if (!rw.valid()) return -102;

      integer_t const liw = gesdd_iwork (a);
      traits::detail::array<integer_t> iw (liw);
      if (!iw.valid()) return -103;

      integer_t info;
      detail::gesdd (jobz, m, n,
                     traits::matrix_storage (a),
                     traits::leading_dimension (a),
                     traits::vector_storage (s),
                     traits::matrix_storage (u),
                     traits::leading_dimension (u),
                     traits::matrix_storage (vt),
                     traits::leading_dimension (vt),
                     traits::vector_storage (w),
                     lw, //traits::vector_size (w),
                     traits::vector_storage (rw),
                     traits::vector_storage (iw),
                     &info);
      return info;
    }


#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_2

    template <typename MatrA, typename VecS, typename MatrU, typename MatrV>
    inline
    int gesdd (char const jobz, MatrA& a, VecS& s, MatrU& u, MatrV& vt) {
      return gesdd ('O', jobz, a, s, u, vt);
    }

    template <typename MatrA, typename VecS, typename MatrU, typename MatrV>
    inline
    int gesdd (MatrA& a, VecS& s, MatrU& u, MatrV& vt) {
      return gesdd ('O', 'S', a, s, u, vt);
    }

    template <typename MatrA, typename VecS>
    int gesdd (MatrA& a, VecS& s) {

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrA>::matrix_structure,
        traits::general_t
      >::value));
#endif 

      integer_t const m = traits::matrix_size1 (a);
      integer_t const n = traits::matrix_size2 (a);
      integer_t const minmn = m < n ? m : n;

      assert (minmn == traits::vector_size (s));

#ifndef BOOST_NUMERIC_BINDINGS_POOR_MANS_TRAITS
      typedef typename traits::matrix_traits<MatrA>::value_type val_t;
#else
      typedef typename MatrA::value_type val_t;
#endif 
      typedef typename traits::type_traits<val_t>::real_type real_t;

      integer_t const lw = gesdd_work ('O', 'N', a);
      traits::detail::array<val_t> w (lw);
      if (!w.valid()) return -101;

      integer_t const lrw = gesdd_rwork ('N', a);
      traits::detail::array<real_t> rw (lrw);
      if (!rw.valid()) return -102;

      integer_t const liw = gesdd_iwork (a);
      traits::detail::array<integer_t> iw (liw);
      if (!iw.valid()) return -103;

      integer_t info;
      detail::gesdd ('N', m, n,
                     traits::matrix_storage (a),
                     traits::leading_dimension (a),
                     traits::vector_storage (s),
                     0, // traits::matrix_storage (u),
                     1, // traits::leading_dimension (u),
                     0, // traits::matrix_storage (vt),
                     1, // traits::leading_dimension (vt),
                     traits::vector_storage (w),
                     traits::vector_size (w),
                     traits::vector_storage (rw),
                     traits::vector_storage (iw),
                     &info);
      return info;
    }

#endif

  } // namespace lapack

}}}

#endif
