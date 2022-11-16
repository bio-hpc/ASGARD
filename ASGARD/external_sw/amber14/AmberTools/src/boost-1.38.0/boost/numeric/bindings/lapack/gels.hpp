//
// Copyright Jesse Manning 2007
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_GELS_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_GELS_HPP

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

        ////////////////////////////////////////////////////////////////////////
        //
        // Linear Least Squares of an underdetermined or overdetermined matrix
        //
        ///////////////////////////////////////////////////////////////////////

        /*  gels - uses the LQ or QR factorization to solve an overdetermined
         *         or underdetermined linear system.  A full rank matrix is
         *         assumed.
         *
         *  The linear least squares system is defined by A*x=b.  A is the m-by-n
         *  coefficients matrix and b is the m-by-nrhs matrix.  Several
         *  right hand side vectors b and solution vectors x can be handled in
         *  a single call; they are stored as the columns of the m-by-nrhs right
         *  hand side matrix B and the n-by-nrh solution matrix x.
         *
         *  If trans = 'N' and m >= n:  find least squares solution of overdetermined system
         *      minimizes || b - A x ||2
         *
         *  if trans = 'N' and m < n: find minimum norm solution of underdetermined system
         *      A*X = B
         *
         *  if trans = 'T' or 'C' and m >= n: find minimum norm soltion of underdetermined system
         *      A^H*X = B
         *
         *  if trans = 'T' or 'C' and m < n: find least squares solution of overdetermined system
         *      minimize || b - A^H x ||2
         *
         *  Workspace is organized following the arguments in the calling sequence.
         *  optimal_workspace() : for optimizing use of blas 3 kernels
         *  minimal_workspace() : minimum size of work arrays, but does not allow for optimization
         *                        of blas 3 kernels
         *  workspace( work )   : where work size must be at least min(m, n) + max(1, m, n, nrhs)
         *
         *
         *  Each solution vector x (length n) is returned column-wise in the rhs matrix B.
         */

        namespace detail {

            inline void gels(char const trans, const integer_t m, const integer_t n,
                             const integer_t nrhs, float *a, const integer_t lda,
                             float *b, const integer_t ldb, float *work,
                             const integer_t lwork, integer_t *info)
            {
                LAPACK_SGELS(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, info);
            }

            inline void gels(char const trans, const integer_t m, const integer_t n,
                             const integer_t nrhs, double *a, const integer_t lda,
                             double *b, const integer_t ldb, double *work,
                             const integer_t lwork, integer_t *info)
            {
                LAPACK_DGELS(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, info);
            }

            inline void gels(char const trans, const integer_t m, const integer_t n,
                             const integer_t nrhs, traits::complex_f *a, const integer_t lda,
                             traits::complex_f *b, const integer_t ldb, traits::complex_f *work,
                             const integer_t lwork, integer_t *info)
            {
                LAPACK_CGELS(&trans, &m, &n, &nrhs,
                             traits::complex_ptr(a), &lda,
                             traits::complex_ptr(b), &ldb,
                             traits::complex_ptr(work), &lwork, info);
            }

            inline void gels(char const trans, const integer_t m, const integer_t n,
                             const integer_t nrhs, traits::complex_d *a, const integer_t lda,
                             traits::complex_d *b, const integer_t ldb, traits::complex_d *work,
                             const integer_t lwork, integer_t *info)
            {
                LAPACK_ZGELS(&trans, &m, &n, &nrhs,
                             traits::complex_ptr(a), &lda,
                             traits::complex_ptr(b), &ldb,
                             traits::complex_ptr(work), &lwork, info);
            }

            // generic function that calls more detailed lapack function
            template <typename MatrA, typename VecB, typename Work>
            int gels(const char trans, MatrA& A, VecB& b, Work& work)
            {
#ifndef NDEBUG
                const std::ptrdiff_t m = traits::matrix_size1(A);
                const std::ptrdiff_t n = traits::matrix_size2(A);
                const std::ptrdiff_t mrhs = traits::matrix_size1(b);
                const std::ptrdiff_t nrhs = traits::matrix_size2(b);
#endif

                // sanity checks
                assert(trans == 'N' || trans == 'T' || trans == 'C');
                assert(m >= 0 && n >= 0);
                assert(mrhs >= 0 && nrhs >= 0);
                assert(traits::leading_dimension(A) >= 1 && traits::leading_dimension(b) >= 1);
                assert(mrhs == std::max(m, n));

                integer_t info;
                detail::gels(trans,
                             traits::matrix_size1(A),
                             traits::matrix_size2(A),
                             traits::matrix_size2(b),
                             traits::matrix_storage(A),
                             traits::leading_dimension(A),
                             traits::matrix_storage(b),
                             traits::leading_dimension(b),
                             traits::vector_storage(work),
                             traits::vector_size(work),
                             &info);

                return info;
            }

            // query for recommended workspace
            template <typename MatrA, typename VecB>
            integer_t gels_optimal_work(const char trans, MatrA& A, VecB& b)
            {
                typename MatrA::value_type work;
                integer_t info;
                detail::gels(trans,
                    traits::matrix_size1(A),
                    traits::matrix_size2(A),
                    traits::matrix_size2(b),
                    traits::matrix_storage(A),
                    traits::leading_dimension(A),
                    traits::matrix_storage(b),
                    traits::leading_dimension(b),
                    &work, //traits::vector_storage(work),
                    -1, // traits::vector_size(work),
                    &info);

                assert(info == 0);

                integer_t lwork = traits::detail::to_int(work);

                return lwork;
            }
        } // namespace detail


        template <typename MatrA, typename VecB>
        int gels(const char trans, MatrA& A, VecB& b, optimal_workspace)
        {
            // query optimal workspace size
            integer_t work_size = detail::gels_optimal_work(trans, A, b);
            traits::detail::array<typename MatrA::value_type> work(work_size);

            return detail::gels(trans, A, b, work);
        }

        template <typename MatrA, typename VecB>
        int gels(const char trans, MatrA& A, VecB& b, minimal_workspace)
        {
            const std::ptrdiff_t m = traits::matrix_size1(A);
            const std::ptrdiff_t n = traits::matrix_size2(A);
            const std::ptrdiff_t r = traits::matrix_size2(b);

            const std::ptrdiff_t minmn = std::min(m, n);       //m < n ? m : n;
            const std::ptrdiff_t maxmn = std::max(m, n);       // m > n ? m : n;
            const std::ptrdiff_t maxdim = std::max(maxmn, r);  // maxmn > r ? maxmn : r;

            traits::detail::array<typename MatrA::value_type> work(minmn + std::max<std::ptrdiff_t>(1, maxdim));

            return detail::gels(trans, A, b, work);
        }

        template <typename MatrA, typename VecB, typename Work>
        int gels(const char trans, MatrA& A, VecB& b, detail::workspace1<Work> workspace)
        {
            return detail::gels(trans, A, b, workspace.select(typename traits::vector_traits<VecB>::value_type()));
        }

    } // namespace lapack

}}}

#endif
