//
// Copyright (c) 2015-2018 CNRS
// Copyright (c) 2015 Wandercraft, 86 rue de Paris 91400 Orsay, France.
//

#ifndef __pinocchio_python_python_hpp__
#define __pinocchio_python_python_hpp__

#include "pinocchio/fwd.hpp"

namespace pinocchio
{
  namespace python
  {
    // Expose spatial classes
    void exposeSE3();
    void exposeForce();
    void exposeMotion();
    void exposeInertia();
    void exposeExplog();
    
    // Expose multibody classes
    void exposeJoints();
    void exposeModel();
    void exposeFrame();
    void exposeData();
    
    // Expose geometry module
    void exposeGeometry();
    
    // Expose parsers
    void exposeParsers();
    
    // Expose algorithms
    void exposeAlgorithms();
    
#ifdef PINOCCHIO_WITH_HPP_FCL
    void exposeFCL();
#endif // PINOCCHIO_WITH_HPP_FCL

  } // namespace python
} // namespace pinocchio

#endif // ifndef __pinocchio_python_python_hpp__

