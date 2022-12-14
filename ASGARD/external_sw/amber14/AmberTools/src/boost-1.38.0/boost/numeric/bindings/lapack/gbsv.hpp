//
// Copyright Vardan Akopian 2007
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_GBSV_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_GBSV_HPP

#include <boost/numeric/bindings/traits/traits.hpp>
#include <boost/numeric/bindings/traits/type_traits.hpp>
#include <boost/numeric/bindings/lapack/lapack.h>


#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK
#  include <boost/static_assert.hpp>
#  include <boost/type_traits/is_same.hpp>
#endif


namespace boost { namespace numeric { namespace bindings {

  namespace lapack {


    namespace detail {
      inline
      void gbtrf (integer_t const n, integer_t const m, integer_t const kl, integer_t const ku,
                  float* ab, integer_t const ldab, integer_t* ipiv, integer_t* info)
      {
        LAPACK_SGBTRF (&n, &m, &kl, &ku, ab, &ldab, ipiv, info);
      }
      inline
      void gbtrf (integer_t const n, integer_t const m, integer_t const kl, integer_t const ku,
                  double* ab, integer_t const ldab, integer_t* ipiv, integer_t* info)
      {
        LAPACK_DGBTRF (&n, &m, &kl, &ku, ab, &ldab, ipiv, info);
      }
      inline
      void gbtrf (integer_t const n, integer_t const m, integer_t const kl, integer_t const ku,
                  traits::complex_f* ab, integer_t const ldab, integer_t* ipiv, integer_t* info)
      {
        LAPACK_CGBTRF (&n, &m, &kl, &ku, traits::complex_ptr(ab), &ldab, ipiv, info);
      }
      inline
      void gbtrf (integer_t const n, integer_t const m, integer_t const kl, integer_t const ku,
                  traits::complex_d* ab, integer_t const ldab, integer_t* ipiv, integer_t* info)
      {
        LAPACK_ZGBTRF (&n, &m, &kl, &ku, traits::complex_ptr(ab), &ldab, ipiv, info);
      }
    }

    template <typename MatrA, typename IVec>
    int gbtrf (MatrA& a, IVec& ipiv) {

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK 
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrA>::matrix_structure,
        traits::banded_t
      >::value));
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrA>::ordering_type,
        traits::column_major_t
      >::value));
#endif

      integer_t const n = traits::matrix_size1 (a);
      integer_t const m = traits::matrix_size2 (a);
      assert (traits::vector_size (ipiv) == (m < n ? m : n));

      // if the matrix has kl lower and ku upper diagonals, then we should have
      // allocated kl lower and kl+ku upper diagonals
      integer_t const kl = traits::matrix_lower_bandwidth (a);
      integer_t const ku = traits::matrix_upper_bandwidth (a) - kl;
      integer_t const ld = traits::leading_dimension (a);

      assert(ku >= 0);

      integer_t info;
      detail::gbtrf (n, m, kl, ku,
                     traits::matrix_storage (a),
                     ld,
                     traits::vector_storage (ipiv),
                     &info);
      return info;
    }


    namespace detail {
      inline
      void gbtrs (char const trans, integer_t const n, integer_t const kl, integer_t const ku, integer_t const m,
                  float const* ab, integer_t const ldab, integer_t const* ipiv,
                  float* b, integer_t const ldb, integer_t* info)
      {
        LAPACK_SGBTRS (&trans, &n, &kl, &ku, &m, ab, &ldab, ipiv, b, &ldb, info);
      }
      inline
      void gbtrs (char const trans, integer_t const n, integer_t const kl, integer_t const ku, integer_t const m,
                  double const* ab, integer_t const ldab, integer_t const* ipiv,
                  double* b, integer_t const ldb, integer_t* info)
      {
        LAPACK_DGBTRS (&trans, &n, &kl, &ku, &m, ab, &ldab, ipiv, b, &ldb, info);
      }
      inline
      void gbtrs (char const trans, integer_t const n, integer_t const kl, integer_t const ku, integer_t const m,
                  traits::complex_f const* ab, integer_t const ldab, integer_t const* ipiv,
                  traits::complex_f* b, integer_t const ldb, integer_t* info)
      {
        LAPACK_CGBTRS (&trans, &n, &kl, &ku, &m, traits::complex_ptr(ab), &ldab, ipiv, traits::complex_ptr(b), &ldb, info);
      }
      inline
      void gbtrs (char const trans, integer_t const n, integer_t const kl, integer_t const ku, integer_t const m,
                  traits::complex_d const* ab, integer_t const ldab, integer_t const* ipiv,
                  traits::complex_d* b, integer_t const ldb, integer_t* info)
      {
        LAPACK_ZGBTRS (&trans, &n, &kl, &ku, &m, traits::complex_ptr(ab), &ldab, ipiv, traits::complex_ptr(b), &ldb, info);
      }
    }


    template <typename MatrA, typename MatrB, typename IVec>
    int gbtrs (char const trans, MatrA const& a, IVec const& ipiv, MatrB& b)
    {
      assert (trans == 'N' || trans == 'T' || trans == 'C');

#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK 
      BOOST_STATIC_ASSERT((boost::is_same<
        typename traits::matrix_traits<MatrA>::matrix_structure,
        traits::banded_t
      >::value));
#endif

      integer_t const n = traits::matrix_size1 (a);
      assert (n == traits::matrix_size2 (a));
      assert (n == traits::matrix_size1 (b));
      assert (n == traits::vector_size (ipiv));

      // if the matrix has kl lower and ku upper diagonals, then we should have
      // allocated kl lower and kl+ku upper diagonals
      integer_t const kl = traits::matrix_lower_bandwidth (a);
      integer_t const ku = traits::matrix_upper_bandwidth (a) - kl;
      integer_t const ld = traits::leading_dimension (a);

      assert(ku >= 0);

      integer_t info;
      detail::gbtrs (trans, n, kl, ku, traits::matrix_size2 (b),
#ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
                     traits::matrix_storage (a),
#else
                     traits::matrix_storage_const (a),
#endif
                     ld,
#ifndef BOOST_NO_FUNCTION_TEMPLATE_ORDERING
                     traits::vector_storage (ipiv),
#else
                     traits::vector_storage_const (ipiv),
#endif
                     traits::matrix_storage (b),
                     traits::leading_dimension (b),
                     &info);
      return info;
    }


  }

}}}

#endif
