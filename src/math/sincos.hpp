//
// Copyright (c) 2015-2016,2018 CNRS
// Copyright (c) 2015 Wandercraft, 86 rue de Paris 91400 Orsay, France.
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

#ifndef __math_sincos_hpp__
#define __math_sincos_hpp__

#include <cmath>

namespace se3
{
  /// Forward declaration
  template<typename Scalar> struct SINCOSAlgo;
  
  
  ///
  /// \brief Computes sin/cos values of a given input scalar.
  ///
  /// \tparam Scalar Type of the input/output variables
  ///
  /// \param[in] a The input scalar from which we evalute the sin and cos.
  /// \param[inout] sa Variable containing the sin of a.
  /// \param[inout] ca Variable containing the cos of a.
  ///
  template<typename Scalar>
  void SINCOS(const Scalar & a, Scalar * sa, Scalar * ca)
  {
    SINCOSAlgo<Scalar>::run(a,sa,ca);
  }
  
  /// Generic evaluation of sin/cos functions.
  template<typename Scalar>
  struct SINCOSAlgo
  {
    static void run(const Scalar & a, Scalar * sa, Scalar * ca)
    {
      (*sa) = std::sin(a); (*ca) = std::cos(a);
    }
  };
  
#ifdef PINOCCHIO_WITH_CPPAD_SUPPORT
  
#include <cppad/cppad.hpp>
  
  /// Implementation for CppAD scalar types.
  template<typename Scalar>
  struct SINCOSAlgo< CppAD::AD<Scalar> >
  {
    static void run(const CppAD::AD<Scalar> & a, CppAD::AD<Scalar> * sa, CppAD::AD<Scalar> * ca)
    {
      (*sa) = CppAD::sin(a); (*ca) = CppAD::cos(a);
    }
  };
  
#endif
  
}

#endif //#ifndef __math_sincos_hpp__
