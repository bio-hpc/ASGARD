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

#ifndef BOOST_NUMERIC_BINDINGS_TRAITS_UBLAS_UPLO_H
#define BOOST_NUMERIC_BINDINGS_TRAITS_UBLAS_UPLO_H

#include <boost/numeric/ublas/fwd.hpp> 

namespace boost { namespace numeric { namespace bindings { namespace traits {

  namespace detail {

    template <typename UpLoTag>
    struct ublas_uplo {};
    
    template<> 
    struct ublas_uplo<boost::numeric::ublas::lower_tag> {
      typedef lower_t type; 
    };

    template<> 
    struct ublas_uplo<boost::numeric::ublas::upper_tag> {
      typedef upper_t type; 
    };

    template<typename I> 
    struct ublas_uplo< boost::numeric::ublas::basic_upper<I> > {
      typedef upper_t type; 
    };

    template<typename I> 
    struct ublas_uplo< boost::numeric::ublas::basic_lower<I> > {
      typedef lower_t type; 
    };

  }

}}}}

#endif 
