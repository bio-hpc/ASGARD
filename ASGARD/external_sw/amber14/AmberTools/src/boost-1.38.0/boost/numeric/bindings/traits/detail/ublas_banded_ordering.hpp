/*
 * 
 * Copyright (c) Kresimir Fresl 2002 
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * Author acknowledges the support of the Faculty of Civil Engineering, 
 * University of Zagreb, Croatia.
 *
 */

#ifndef BOOST_NUMERIC_BINDINGS_TRAITS_DETAIL_UBLAS_BANDED_ORDERING_H
#define BOOST_NUMERIC_BINDINGS_TRAITS_DETAIL_UBLAS_BANDED_ORDERING_H

#include <boost/numeric/ublas/fwd.hpp> 

namespace boost { namespace numeric { namespace bindings { namespace traits {

  namespace detail {

    template <typename StOrdTag>
    struct ublas_banded_ordering {};
    
    template<> 
    struct ublas_banded_ordering<boost::numeric::ublas::row_major_tag> {
      // When orientation_category==row_major_tag then the ublas banded format corresponds to
      // the LAPACK band format, which really is a column_major format.
      typedef column_major_t                        type; 

      template <typename M>
      static typename M::size_type leading_dimension( M const& m ) {
        return m.lower() + m.upper() + 1 ;
      }

      template <typename M>
      static typename M::size_type stride1( M const& m ) {
        return 1 ;
      }

      template <typename M>
      static typename M::size_type stride2( M const& m ) {
        return leading_dimension(m)-1 ;
      }
    };
    
    template<> 
    struct ublas_banded_ordering<boost::numeric::ublas::column_major_tag> {
      // The type row_major_t is just used to indicate that this is not a column_major format.
      typedef row_major_t                        type; 

      template <typename M>
      static typename M::size_type leading_dimension( M const& m ) {
        return m.size2() ;
      }

      template <typename M>
      static typename M::size_type stride1( M const& m ) {
        return leading_dimension(m) ;
      }

      template <typename M>
      static typename M::size_type stride2( M const& m ) {
        return 1-leading_dimension(m) ;
      }
    };
  }

}}}}

#endif 
