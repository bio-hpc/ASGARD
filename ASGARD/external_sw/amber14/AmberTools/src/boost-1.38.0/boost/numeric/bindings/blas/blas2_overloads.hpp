//
//  Copyright (C) Toon Knapen 2003
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_BLAS_BLAS2_OVERLOADS_HPP
#define BOOST_NUMERIC_BINDINGS_BLAS_BLAS2_OVERLOADS_HPP

#include <boost/numeric/bindings/blas/blas.h>
#include <boost/numeric/bindings/traits/type_traits.hpp>

namespace boost { namespace numeric { namespace bindings { namespace blas { namespace detail {

  using namespace boost::numeric::bindings::traits ;

  inline
  void gemv( char TRANS, const integer_t& m, const integer_t& n, const float    & alpha, const float    * a_ptr, const integer_t& lda, const float    * x_ptr, const integer_t& incx, const float    & beta, float    * y_ptr, const integer_t& incy ) { BLAS_SGEMV( &TRANS, &m, &n,          ( &alpha ),          ( a_ptr ), &lda,          ( x_ptr ), &incx,          ( &beta ),          ( y_ptr ), &incy ) ; }
  inline
  void gemv( char TRANS, const integer_t& m, const integer_t& n, const double   & alpha, const double   * a_ptr, const integer_t& lda, const double   * x_ptr, const integer_t& incx, const double   & beta, double   * y_ptr, const integer_t& incy ) { BLAS_DGEMV( &TRANS, &m, &n,          ( &alpha ),          ( a_ptr ), &lda,          ( x_ptr ), &incx,          ( &beta ),          ( y_ptr ), &incy ) ; }
  inline
  void gemv( char TRANS, const integer_t& m, const integer_t& n, const complex_f& alpha, const complex_f* a_ptr, const integer_t& lda, const complex_f* x_ptr, const integer_t& incx, const complex_f& beta, complex_f* y_ptr, const integer_t& incy ) { BLAS_CGEMV( &TRANS, &m, &n, complex_ptr( &alpha ), complex_ptr( a_ptr ), &lda, complex_ptr( x_ptr ), &incx, complex_ptr( &beta ), complex_ptr( y_ptr ), &incy ) ; }
  inline
  void gemv( char TRANS, const integer_t& m, const integer_t& n, const complex_d& alpha, const complex_d* a_ptr, const integer_t& lda, const complex_d* x_ptr, const integer_t& incx, const complex_d& beta, complex_d* y_ptr, const integer_t& incy ) { BLAS_ZGEMV( &TRANS, &m, &n, complex_ptr( &alpha ), complex_ptr( a_ptr ), &lda, complex_ptr( x_ptr ), &incx, complex_ptr( &beta ), complex_ptr( y_ptr ), &incy ) ; }

  inline
  void ger( const integer_t& m, const integer_t& n, const float  & alpha, const float  * x_ptr, const integer_t& incx, const float  * y_ptr, const integer_t& incy, float  * a_ptr, const integer_t& lda ) { BLAS_SGER( &m, &n, &alpha, x_ptr, &incx, y_ptr, &incy, a_ptr, &lda ) ; }
  inline
  void ger( const integer_t& m, const integer_t& n, const double & alpha, const double * x_ptr, const integer_t& incx, const double * y_ptr, const integer_t& incy, double * a_ptr, const integer_t& lda ) { BLAS_DGER( &m, &n, &alpha, x_ptr, &incx, y_ptr, &incy, a_ptr, &lda ) ; }
/*  
  inline
  void geru( const integer_t& m, const integer_t& n, const complex_f & alpha, const complex_f * x_ptr, const integer_t& incx, complex_f * y_ptr, const integer_t& incy, complex_f * a_ptr, const integer_t& lda ) { BLAS_CGERU( &m, &n, complex_ptr( &alpha ), complex_ptr( x_ptr ), &incx, complex_ptr( y_ptr ), &incy, complex_ptr( a_ptr ), &lda ) ; }
  inline
  void geru( const integer_t& m, const integer_t& n, const complex_d & alpha, const complex_d * x_ptr, const integer_t& incx, complex_d * y_ptr, const integer_t& incy, complex_d * a_ptr, const integer_t& lda ) { BLAS_ZGERU( &m, &n, complex_ptr( &alpha ), complex_ptr( x_ptr ), &incx, complex_ptr( y_ptr ), &incy, complex_ptr( a_ptr ), &lda ) ; }
*/
}}}}}

#endif // BOOST_NUMERIC_BINDINGS_BLAS_BLAS2_OVERLOADS_HPP

