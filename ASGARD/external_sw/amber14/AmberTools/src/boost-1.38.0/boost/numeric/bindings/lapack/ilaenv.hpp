/*
 *
 * Copyright (c) Toon Knapen, Kresimir Fresl and Matthias Troyer 2003
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * KF acknowledges the support of the Faculty of Civil Engineering,
 * University of Zagreb, Croatia.
 *
 */

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_ILAENV_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_ILAENV_HPP

#include <cstring>
#include <boost/numeric/bindings/lapack/lapack.h>

namespace boost { namespace numeric { namespace bindings { namespace lapack {

  /*
   * ilaenv() is called from the LAPACK routines to choose
   * problem-dependent parameters such as the block sizes
   * for the local environment.
   */
  
  inline
  integer_t ilaenv (integer_t const ispec, const char* name, const char* opts,
              integer_t const n1 = -1, integer_t const n2 = -1,
              integer_t const n3 = -1, integer_t const n4 = -1)
  {
    return ::LAPACK_ILAENV (&ispec, name, opts, &n1, &n2, &n3, &n4,
                          std::strlen (name), std::strlen (opts));
  }


}}}}

#endif
