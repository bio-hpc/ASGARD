
#ifndef BOOST_MPL_SET_AUX_SIZE_IMPL_HPP_INCLUDED
#define BOOST_MPL_SET_AUX_SIZE_IMPL_HPP_INCLUDED

// Copyright Aleksey Gurtovoy 2003-2004
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/mpl for documentation.

// $Id: size_impl.hpp,v 1.1 2010/03/29 19:23:44 peters Exp $
// $Date: 2010/03/29 19:23:44 $
// $Revision: 1.1 $

#include <boost/mpl/size_fwd.hpp>
#include <boost/mpl/set/aux_/tag.hpp>

namespace boost { namespace mpl {

template<>
struct size_impl< aux::set_tag >
{
    template< typename Set > struct apply
        : Set::size
    {
    };
};

}}

#endif // BOOST_MPL_SET_AUX_SIZE_IMPL_HPP_INCLUDED
