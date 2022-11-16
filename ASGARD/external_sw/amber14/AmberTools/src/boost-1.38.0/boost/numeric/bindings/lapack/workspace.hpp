/*
 *
 * Copyright (c) Karl Meerbergen & Kresimir Fresl 2003
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * KF acknowledges the support of the Faculty of Civil Engineering,
 * University of Zagreb, Croatia.
 *
 */

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_WORKSPACE_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_WORKSPACE_HPP

#include <boost/numeric/bindings/traits/detail/array_impl.hpp>
#include <boost/numeric/bindings/traits/type.hpp>
#include <boost/numeric/bindings/traits/type_traits.hpp>
#include <boost/numeric/bindings/traits/vector_traits.hpp>
#include <memory>

namespace boost { namespace numeric { namespace bindings {

  /*
   * Organization of workspace in Lapack.
   * We allow one of the following arguments in a number of Lapack functions
   * - minimal_workspace() : the function allocates the minimum workspace required for the function
   * - optimal_workspace() : the function allocates the amount of workspace that allows optimal
   *                         execution.
   * - workspace( work )   : the function uses the workspace array in work.
   * - workspace( rwork, work ) : the function uses a real array rwork and a compolex array work as
   *                              workspace. (There are Lapack functions for complex matrices
   *                              that require two workarrays)
   * */

  namespace lapack {


     // Four classes are introduced to distinguish between the different type of memory allocations

     struct minimal_workspace {} ;

     struct optimal_workspace {} ;

     namespace detail {

        template <typename W>
        class workspace1 {
          public:
            workspace1(W& w)
            : w_( w )
            {}

          public:
            typedef typename traits::vector_traits<W>::value_type value_type ;
            W& select( value_type const& ) { return w_ ; }

          private:
            W& w_ ;
        }; // struct workspace1

        template <typename W, typename WRI>
        class workspace2 {
          public:
            workspace2(W& w, WRI& wri)
            : w_(w)
            , wri_(wri)
            {}

          public:
            typedef typename traits::vector_traits<W>::value_type w_value_type ;
            W& select( w_value_type const& ) { return w_ ; }

            typedef typename traits::vector_traits<WRI>::value_type wri_value_type ;
            WRI& select( wri_value_type const& ) { return wri_ ; }

          private:
            W& w_ ;
            WRI& wri_ ;
        }; // struct workspace2

        template <typename W, typename WR, typename WI>
        class workspace3 {
          public:
            workspace3(W& w, WR& wr, WI& wi)
            : w_(w)
            , wr_(wr)
            , wi_(wi)
            {}

          public:
            typedef typename traits::vector_traits<W>::value_type w_value_type ;
            W& select( w_value_type const& ) { return w_ ; }

            typedef typename traits::vector_traits<WR>::value_type wr_value_type ;
            WR& select( wr_value_type const& ) { return wr_ ; }

            typedef typename traits::vector_traits<WI>::value_type wi_value_type ;
            WI& select( wi_value_type const& ) { return wi_ ; }

          private:
            W& w_ ;
            WR& wr_ ;
            WI& wi_ ;
        }; // struct workspace3

        template <typename W, typename WR, typename WI, typename WB>
        class workspace4 {
          public:
            workspace4(W& w, WR& wr, WI& wi, WB& wb)
            : w_(w)
            , wr_(wr)
            , wi_(wi)
            , wb_(wb)
            {}

          public:
            typedef typename traits::vector_traits<W>::value_type w_value_type ;
            W& select( w_value_type const& ) { return w_ ; }

            typedef typename traits::vector_traits<WR>::value_type wr_value_type ;
            WR& select( wr_value_type const& ) { return wr_ ; }

            typedef typename traits::vector_traits<WI>::value_type wi_value_type ;
            WI& select( wi_value_type const& ) { return wi_ ; }

            typedef typename traits::vector_traits<WB>::value_type wb_value_type ;
            WB& select( wb_value_type const& ) { return wb_ ; }

          private:
            W& w_ ;
            WR& wr_ ;
            WI& wi_ ;
            WB& wb_ ;
        }; // struct workspace4

     }

     template <typename W>
     detail::workspace1<W> workspace(W& w) {
        return detail::workspace1<W>(w) ;
     } // workspace()

     //
     // Two situations:
     //   Real valued: workspace( real array, integer array )
     //   Complex valued: workspace( complex array, real array )
     //
     template <typename W, typename WRI>
     detail::workspace2<W,WRI> workspace(W& w, WRI& wri) {
        return detail::workspace2<W,WRI>(w, wri) ;
     } // workspace()

     //
     //   Complex valued: workspace( complex array, real array, integer array )
     //
     template <typename W, typename WR, typename WI>
     detail::workspace3<W,WR,WI> workspace(W& w, WR& wr, WI& wi) {
        return detail::workspace3<W,WR,WI>(w, wr, wi) ;
     } // workspace()

     //
     //   Complex valued: workspace( complex array, real array, integer array, bool array )
     //
     template <typename W, typename WR, typename WI, typename WB>
     detail::workspace4<W,WR,WI,WB> workspace(W& w, WR& wr, WI& wi, WB& wb) {
        return detail::workspace4<W,WR,WI,WB>(w, wr, wi, wb) ;
     } // workspace()

     /// Select the number of workspaces depending on the value_type
     template <typename T>
     struct n_workspace_args { };

     template <>
     struct n_workspace_args< float > {
        static const int value = 1 ;
     };

     template <>
     struct n_workspace_args< double > {
        static const int value = 1 ;
     };

     template <>
     struct n_workspace_args< traits::complex_f > {
        static const int value = 2 ;
     };

     template <>
     struct n_workspace_args< traits::complex_d > {
        static const int value = 2 ;
     };

  }

}}}

#endif
