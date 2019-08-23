//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/special_functions/factorials.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/special_functions/prime.hpp>
//
// Note this header includes no other headers, this is
// important if this test is to be meaningful:
//
#include "test_compile_result.hpp"

void compile_and_link_test()
{
   check_result<boost::uint32_t>(boost::math::prime(u));
   //
   // Add constexpr tests here:
   //
#ifdef BOOST_MATH_HAVE_CONSTEXPR_TABLES
   constexpr boost::uint32_t ce_f = boost::math::prime(2);
   check_result<boost::uint32_t>(ce_f);
#endif
}

