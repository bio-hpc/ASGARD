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

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_TREVC_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_TREVC_HPP

#include <complex>
#include <boost/numeric/bindings/traits/traits.hpp>
#include <boost/numeric/bindings/traits/type_traits.hpp>
#include <boost/numeric/bindings/lapack/lapack.h>
#include <boost/numeric/bindings/traits/detail/array.hpp>

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
#  include <boost/static_assert.hpp>
#  include <boost/type_traits.hpp>
#endif


namespace boost { namespace numeric { namespace bindings {

  namespace lapack {

    ///////////////////////////////////////////////////////////////////
    //
    // Compute eigenvectors of Schur matrix (computed by gees).
    //
    ///////////////////////////////////////////////////////////////////

    /*
     * trevc() computes the eigenvectors using the Schur factorization
     *
     * Let  A = U * S * herm(U), then trecv computes the eigenvectors of
     * S, and may optionally apply them to U.
     *
     * To compute the Schur factorization, see gees.
     */

    namespace detail {
      inline
      void trevc (char const side, char const howmny, const logical_t* select, integer_t const n,
                 float* t, integer_t const ldt, float* vl, integer_t const ldvl, float* vr, integer_t const ldvr,
                 integer_t const mm, integer_t& m, float* work, integer_t& info)
      {
        LAPACK_STREVC (&side, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, &mm, &m, work, &info);
      }

      inline
      void trevc (char const side, char const howmny, const logical_t* select, integer_t const n,
                 double* t, integer_t const ldt, double* vl, integer_t const ldvl, double* vr, integer_t const ldvr,
                 integer_t const mm, integer_t& m, double* work, integer_t& info)
      {
        LAPACK_DTREVC (&side, &howmny, select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, &mm, &m, work, &info);
      }

      inline
      void trevc (char const side, char const howmny, const logical_t* select, integer_t const n,
                 traits::complex_f* t, integer_t const ldt, traits::complex_f* vl, integer_t const ldvl, traits::complex_f* vr, integer_t const ldvr,
                 integer_t const mm, integer_t& m, traits::complex_f* work, integer_t& info)
      {
        LAPACK_CTREVC (&side, &howmny, select, &n, traits::complex_ptr(t), &ldt, traits::complex_ptr(vl), &ldvl,
                        traits::complex_ptr(vr), &ldvr, &mm, &m, traits::complex_ptr(work+n), reinterpret_cast<float*>(work), &info);
      }

      inline
      void trevc (char const side, char const howmny, const logical_t* select, integer_t const n,
                  traits::complex_d* t, integer_t const ldt, traits::complex_d* vl, integer_t const ldvl, traits::complex_d* vr, integer_t const ldvr,
                  integer_t const mm, integer_t& m, traits::complex_d* work, integer_t& info)
      {
        LAPACK_ZTREVC (&side, &howmny, select, &n, traits::complex_ptr(t), &ldt,
                       traits::complex_ptr(vl), &ldvl, traits::complex_ptr(vr), &ldvr,
                       &mm, &m, traits::complex_ptr(work+n), reinterpret_cast<double*>(work), &info);
      }

    }

    // Compute Schur factorization with Schur vectors
    template <typename MatrT, typename VL, typename VR, typename Work>
    int trevc (char const side, char const howmny, MatrT& t, VL& vl, VR& vr, Work& work) {

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrT>::matrix_structure,
        traits::general_t
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<VL>::matrix_structure,
        traits::general_t
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<VR>::matrix_structure,
        traits::general_t
      >::value));
#endif

      integer_t const n = traits::matrix_size1 (t);
      assert (n == traits::matrix_size2 (t));
      assert (n == traits::matrix_size1 (vl));
      assert (n == traits::matrix_size2 (vl));
      assert (n == traits::matrix_size1 (vr));
      assert (n == traits::matrix_size2 (vr));
      assert (3*n <= traits::vector_size (work));

      logical_t* select=0;

      integer_t mm=n;
      integer_t m;
      integer_t info;
      detail::trevc (side, howmny, select, n,
                    traits::matrix_storage (t),
                    traits::leading_dimension (t),
                    traits::matrix_storage (vl),
                    traits::leading_dimension (vl),
                    traits::matrix_storage (vr),
                    traits::leading_dimension (vr),
                    mm,
                    m,
                    traits::vector_storage (work),
                    info);
      return info;
    }

  }

}}}

#endif
