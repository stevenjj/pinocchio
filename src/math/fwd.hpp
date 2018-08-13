//
// Copyright (c) 2016,2018 CNRS
//
// This file is part of Pinocchio
// Pinocchio is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// Pinocchio is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// Pinocchio If not, see
// <http://www.gnu.org/licenses/>.

#ifndef __se3_math_fwd_hpp__
#define __se3_math_fwd_hpp__

#include "pinocchio/fwd.hpp"
#include <math.h>
#include <boost/math/constants/constants.hpp>

#ifdef PINOCCHIO_WITH_CPPAD_SUPPORT
namespace boost
{
  namespace math
  {
    namespace constants
    {
      namespace detail
      {
        template<typename Scalar>
        struct constant_pi< CppAD::AD<Scalar> > : constant_pi<Scalar> {};
      }
    }
  }
}
#endif

namespace se3
{
  ///
  /// \brief Returns the value of PI according to the template parameters Scalar
  ///
  /// \tparam Scalar The scalar type of the return pi value
  ///
  template<typename Scalar>
  const Scalar PI()
  { return boost::math::constants::pi<Scalar>(); }
  
  /// The value of PI for double scalar type
  const double PId = PI<double>();
  
  namespace math
  {
    using std::fabs;
    using std::sqrt;
    
#ifdef PINOCCHIO_WITH_CPPAD_SUPPORT
    using CppAD::fabs;
    using CppAD::sqrt;
#endif
  }
}

#endif //#ifndef __se3_math_fwd_hpp__
