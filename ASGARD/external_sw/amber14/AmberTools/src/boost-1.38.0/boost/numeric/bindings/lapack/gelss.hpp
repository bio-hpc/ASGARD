//
// Copyright Jesse Manning 2007
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_GELSS_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_GELSS_HPP

#include <algorithm>

#include <boost/numeric/bindings/traits/type.hpp>
#include <boost/numeric/bindings/traits/traits.hpp>
#include <boost/numeric/bindings/traits/type_traits.hpp>
#include <boost/numeric/bindings/lapack/lapack.h>
#include <boost/numeric/bindings/lapack/workspace.hpp>
#include <boost/numeric/bindings/traits/detail/array.hpp>
#include <boost/numeric/bindings/traits/detail/utils.hpp>


#ifndef BOOST_NUMERIC_BINDINGS_NO_STRUCTURE_CHECK 
#  include <boost/static_assert.hpp>
#  include <boost/type_traits.hpp>
#endif

namespace boost { namespace numeric { namespace bindings {

    namespace lapack {

        namespace detail {

            inline void gelss(const integer_t m, const integer_t n, const integer_t nrhs,
                              float *a, const integer_t lda, float *b, const integer_t ldb,
                              float *s, const float rcond, integer_t *rank, float *work,
                              const integer_t lwork, integer_t *info)
            {
                LAPACK_SGELSS(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, info);
            }

            inline void gelss(const integer_t m, const integer_t n, const integer_t nrhs,
                              double *a, const integer_t lda, double *b, const integer_t ldb,
                              double *s, const double rcond, integer_t *rank, double *work,
                              const integer_t lwork, integer_t *info)
            {
                LAPACK_DGELSS(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, rank, work, &lwork, info);
            }

            inline void gelss(const integer_t m, const integer_t n, const integer_t nrhs,
                              traits::complex_f *a, const integer_t lda, traits::complex_f *b,
                              const integer_t ldb, float *s, const float rcond, integer_t *rank,
                              traits::complex_f *work, const integer_t lwork, float *rwork, integer_t *info)
            {
                LAPACK_CGELSS(&m, &n, &nrhs, traits::complex_ptr(a),
                              &lda, traits::complex_ptr(b), &ldb, s,
                              &rcond, rank, traits::complex_ptr(work),
                              &lwork, rwork, info);
            }

            inline void gelss(const integer_t m, const integer_t n, const integer_t nrhs,
                              traits::complex_d *a, const integer_t lda, traits::complex_d *b,
                              const integer_t ldb, double *s, const double rcond, integer_t *rank,
                              traits::complex_d *work, const integer_t lwork, double *rwork, integer_t *info)
            {
                LAPACK_ZGELSS(&m, &n, &nrhs, traits::complex_ptr(a),
                              &lda, traits::complex_ptr(b), &ldb, s,
                              &rcond, rank, traits::complex_ptr(work),
                              &lwork, rwork, info);
            }

            // gelss for real type
            template <typename MatrA, typename MatrB, typename VecS, typename Work>
            int gelss(MatrA& A, MatrB& B, VecS& s, Work& work)
            {
                typedef typename MatrA::value_type val_t;
                typedef typename traits::type_traits<val_t>::real_type real_t;

                const std::ptrdiff_t m = traits::matrix_size1(A);
                const std::ptrdiff_t n = traits::matrix_size2(A);
                const std::ptrdiff_t nrhs = traits::matrix_size2(B);
                const std::ptrdiff_t maxmn = std::max(m, n);
                const std::ptrdiff_t minmn = std::min(m, n);

                // sanity checks
                assert(m >= 0 && n >= 0);
                assert(nrhs >= 0);
                assert(traits::leading_dimension(A) >= std::max<std::ptrdiff_t>(1, m));
                assert(traits::leading_dimension(B) >= std::max<std::ptrdiff_t>(1, maxmn));
                assert(traits::vector_size(work) >= 1);
                assert(traits::vector_size(s) >= std::max<std::ptrdiff_t>(1, minmn));

                integer_t info;
                const real_t rcond = -1;    // use machine precision
                integer_t rank;

                detail::gelss(traits::matrix_size1(A),
                              traits::matrix_size2(A),
                              traits::matrix_size2(B),
                              traits::matrix_storage(A),
                              traits::leading_dimension(A),
                              traits::matrix_storage(B),
                              traits::leading_dimension(B),
                              traits::vector_storage(s),
                              rcond,
                              &rank,
                              traits::vector_storage(work),
                              traits::vector_size(work),
                              &info);

                return info;
            }

            // gelss for complex type
            template <typename MatrA, typename MatrB, typename VecS, typename Work, typename RWork>
            int gelss(MatrA& A, MatrB& B, VecS& s, Work& work, RWork& rwork)
            {
                typedef typename MatrA::value_type val_t;
                typedef typename traits::type_traits<val_t>::real_type real_t;

                const std::ptrdiff_t m = traits::matrix_size1(A);
                const std::ptrdiff_t n = traits::matrix_size2(A);
                const std::ptrdiff_t nrhs = traits::matrix_size2(B);
                const std::ptrdiff_t maxmn = std::max(m, n);
                const std::ptrdiff_t minmn = std::min(m, n);

                // sanity checks
                assert(m >= 0 && n >= 0);
                assert(nrhs >= 0);
                assert(traits::leading_dimension(A) >= std::max<std::ptrdiff_t>(1, m));
                assert(traits::leading_dimension(B) >= std::max<std::ptrdiff_t>(1, maxmn));
                assert(traits::vector_size(work) >= 1);
                assert(traits::vector_size(s) >= std::max<std::ptrdiff_t>(1, minmn));

                integer_t info;
                const real_t rcond = -1;    // use machine precision
                integer_t rank;

                detail::gelss(traits::matrix_size1(A),
                              traits::matrix_size2(A),
                              traits::matrix_size2(B),
                              traits::matrix_storage(A),
                              traits::leading_dimension(A),
                              traits::matrix_storage(B),
                              traits::leading_dimension(B),
                              traits::vector_storage(s),
                              rcond,
                              &rank,
                              traits::vector_storage(work),
                              traits::vector_size(work),
                              traits::vector_storage(rwork),
                              &info);

                return info;
            }

            // default minimal workspace functor
            template <int N>
            struct Gelss { };

            // specialization for gelss (sgelss, dgelss)
            template <>
            struct Gelss<1>
            {
                template <typename MatrA, typename MatrB, typename VecS>
                int operator() (MatrA& A, MatrB& B, VecS& s, minimal_workspace) const
                {
                    typedef typename traits::matrix_traits<MatrA>::value_type val_t;

                    const std::ptrdiff_t m = traits::matrix_size1(A);
                    const std::ptrdiff_t n = traits::matrix_size2(A);
                    const std::ptrdiff_t rhs = traits::matrix_size2(B);

                    const std::ptrdiff_t minmn = std::min(m, n);           // minmn = m < n ? m : n
                    const std::ptrdiff_t maxmn = std::max(m, n);           // maxmn = m > n ? m : n
                    const std::ptrdiff_t maxmnr = std::max(maxmn, rhs);    // maxmnr = maxmn > rhs ? maxmn : rhs

                    traits::detail::array<val_t> work(3*minmn + std::max(2*minmn, maxmnr));

                    return gelss(A, B, s, work);
                }

                template <typename MatrA, typename MatrB, typename VecS>
                int operator() (MatrA& A, MatrB& B, VecS& s, optimal_workspace) const
                {
                    typedef typename traits::matrix_traits<MatrA>::value_type val_t;
                    typedef typename traits::type_traits<val_t>::real_type real_t;

                    //const std::ptrdiff_t m = traits::matrix_size1(A);
                    //const std::ptrdiff_t n = traits::matrix_size2(A);
                    //const std::ptrdiff_t rhs = traits::matrix_size2(B);

                    //const std::ptrdiff_t minmn = std::min(m, n);           // minmn = m < n ? m : n
                    //const std::ptrdiff_t maxmn = std::max(m, n);           // maxmn = m > n ? m : n
                    //const std::ptrdiff_t maxmnr = std::max(maxmn, rhs);    // maxmnr = maxmn > rhs ? maxmn : rhs

                    val_t temp_work;

                    const real_t rcond = -1;
                    integer_t rank;
                    integer_t info;

                    // query for optimal workspace size
                    detail::gelss(traits::matrix_size1(A),
                        traits::matrix_size2(A),
                        traits::matrix_size2(B),
                        traits::matrix_storage(A),
                        traits::leading_dimension(A),
                        traits::matrix_storage(B),
                        traits::leading_dimension(B),
                        traits::vector_storage(s),
                        rcond,
                        &rank,
                        &temp_work, //traits::vector_storage(work),
                        -1,         //traits::vector_size(work),
                        &info);

                    assert(info == 0);

                    const integer_t lwork = traits::detail::to_int(temp_work);

                    traits::detail::array<val_t> work(lwork);

                    return gelss(A, B, s, work);
                }

                template <typename MatrA, typename MatrB, typename VecS, typename Work>
                int operator() (MatrA& A, MatrB& B, VecS& s, detail::workspace1<Work> workspace) const
                {
                    return gelss(A, B, s, workspace.select(typename traits::matrix_traits<MatrA>::value_type()));
                }
            };

            // specialization for gelss (cgelss, zgelss)
            template <>
            struct Gelss<2>
            {
                template <typename MatrA, typename MatrB, typename VecS>
                int operator() (MatrA& A, MatrB& B, VecS& s, minimal_workspace) const
                {
                    typedef typename traits::matrix_traits<MatrA>::value_type val_t;
                    typedef typename traits::type_traits<val_t>::real_type real_t;

                    const std::ptrdiff_t m = traits::matrix_size1(A);
                    const std::ptrdiff_t n = traits::matrix_size2(A);
                    const std::ptrdiff_t rhs = traits::matrix_size2(B);

                    const std::ptrdiff_t minmn = std::min(m, n);           // minmn = m < n ? m : n
                    const std::ptrdiff_t maxmn = std::max(m, n);           // maxmn = m > n ? m : n
                    const std::ptrdiff_t maxmnr = std::max(maxmn, rhs);    // maxmnr = maxmn > rhs ? maxmn : rhs

                    traits::detail::array<val_t> work(2*minmn + maxmnr);
                    traits::detail::array<real_t> rwork(std::max<std::ptrdiff_t>(1, (5*minmn)));

                    return gelss(A, B, s, work, rwork);
                }

                template <typename MatrA, typename MatrB, typename VecS>
                int operator() (MatrA& A, MatrB& B, VecS& s, optimal_workspace) const
                {
                    typedef typename MatrA::value_type val_t;
                    typedef typename traits::type_traits<val_t>::real_type real_t;

                    const std::ptrdiff_t m = traits::matrix_size1(A);
                    const std::ptrdiff_t n = traits::matrix_size2(A);
                    //const std::ptrdiff_t rhs = traits::matrix_size2(B);

                    const std::ptrdiff_t minmn = std::min(m, n);           // minmn = m < n ? m : n
                    //const std::ptrdiff_t maxmn = std::max(m, n);           // maxmn = m > n ? m : n
                    //const std::ptrdiff_t maxmnr = std::max(maxmn, rhs);    // maxmnr = maxmn > rhs ? maxmn : rhs

                    val_t temp_work;
                    real_t temp_rwork;

                    const real_t rcond = -1;
                    integer_t rank;
                    integer_t info;

                    // query for optimal workspace size
                    detail::gelss(traits::matrix_size1(A),
                                  traits::matrix_size2(A),
                                  traits::matrix_size2(B),
                                  traits::matrix_storage(A),
                                  traits::leading_dimension(A),
                                  traits::matrix_storage(B),
                                  traits::leading_dimension(B),
                                  traits::vector_storage(s),
                                  rcond,
                                  &rank,
                                  &temp_work,   //traits::vector_storage(work),
                                  -1,           //traits::vector_size(work),
                                  &temp_rwork,
                                  &info);

                    assert(info == 0);

                    const integer_t lwork = traits::detail::to_int(temp_work);

                    traits::detail::array<val_t> work(lwork);
                    traits::detail::array<real_t> rwork(std::max<std::ptrdiff_t>(1, (5*minmn)));

                    return gelss(A, B, s, work, rwork);
                }

                template <typename MatrA, typename MatrB, typename VecS, typename Work, typename RWork>
                int operator() (MatrA& A, MatrB& B, VecS& s, detail::workspace2<Work, RWork> workspace) const
                {
                  typedef typename traits::matrix_traits<MatrA>::value_type    value_type ;
                  typedef typename traits::type_traits<value_type>::real_type real_type ;
                    return gelss(A, B, s, workspace.select(value_type()), workspace.select(real_type()));
                }
            };

        } // detail

        // gelss
        // Parameters:
        //  A:          matrix of coefficients
        //  B:          matrix of solutions (stored column-wise)
        //  s:          vector to store singular values on output, length >= max(1, min(m,n))
        //  workspace:  either optimal, minimal, or user supplied
        //
        template <typename MatrA, typename MatrB, typename VecS, typename Work>
        int gelss(MatrA& A, MatrB& B, VecS& s, Work workspace)
        {
            typedef typename traits::matrix_traits<MatrA>::value_type val_t;

            return detail::Gelss<n_workspace_args<val_t>::value>() (A, B, s, workspace);
        }

        // gelss, no singular values are returned
        // Parameters:
        //  A:          matrix of coefficients
        //  B:          matrix of solutions (stored column-wise)
        //  workspace:  either optimal, minimal, or user supplied
        //
        template <typename MatrA, typename MatrB, typename Work>
        int gelss(MatrA& A, MatrB& B, Work workspace)
        {
            typedef typename traits::matrix_traits<MatrA>::value_type val_t;
            typedef typename traits::type_traits<val_t>::real_type real_t;

            const std::ptrdiff_t m = traits::matrix_size1(A);
            const std::ptrdiff_t n = traits::matrix_size2(A);

            const std::ptrdiff_t s_size = std::max<std::ptrdiff_t>(1, std::min(m,n));
            traits::detail::array<real_t> s(s_size);

            return detail::Gelss<n_workspace_args<val_t>::value>() (A, B, s, workspace);
        }

    } // namespace lapack

}}}

#endif
