//
// Copyright (c) 2016 CNRS
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

#ifndef __se3_python_frame_hpp__
#define __se3_python_frame_hpp__

#include <eigenpy/exception.hpp>
#include <eigenpy/eigenpy.hpp>

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "pinocchio/multibody/frame.hpp"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/container/aligned-vector.hpp"
#include "pinocchio/bindings/python/utils/copyable.hpp"

namespace se3
{
  namespace python
  {
    namespace bp = boost::python;

    struct FramePythonVisitor
      : public boost::python::def_visitor< FramePythonVisitor >
    {
      typedef Model::Index Index;
      typedef Model::JointIndex JointIndex;
      typedef Model::FrameIndex FrameIndex;

    public:

      template<class PyClass>
      void visit(PyClass& cl) const 
      {
        cl
          .def(bp::init< const std::string&,const JointIndex, const FrameIndex, const SE3&,FrameType> ((bp::arg("name (string)"),bp::arg("index of parent joint"), bp::args("index of parent frame"), bp::arg("SE3 placement"), bp::arg("type (FrameType)")),
                "Initialize from name, parent joint id, parent frame id and placement wrt parent joint."))

          .def_readwrite("name", &Frame::name, "name  of the frame")
          .def_readwrite("parent", &Frame::parent, "id of the parent joint")
          .def_readwrite("previousFrame", &Frame::previousFrame, "id of the previous frame") 
          .add_property("placement", 
                        &FramePythonVisitor::getPlacementWrtParentJoint, 
                        &FramePythonVisitor::setPlacementWrtParentJoint, 
                        "placement in the parent joint local frame")
          .def_readwrite("type", &Frame::type, "type of the frame")

          .def(bp::self_ns::str(bp::self_ns::self))
          .def(bp::self_ns::repr(bp::self_ns::self))
          ;
      }


      static SE3 getPlacementWrtParentJoint(const Frame & self) { return self.placement; }
      static void setPlacementWrtParentJoint(Frame & self, const SE3 & placement) { self.placement = placement; }

      static void expose()
      {
        bp::enum_<FrameType>("FrameType")
            .value("OP_FRAME",OP_FRAME)
            .value("JOINT",JOINT)
            .value("FIXED_JOINT",FIXED_JOINT)
            .value("BODY",BODY)
            .value("SENSOR",SENSOR)
            ;

        bp::class_<Frame>("Frame",
                           "A Plucker coordinate frame related to a parent joint inside a kinematic tree.\n\n",
	                         bp::no_init
                         )
        .def(FramePythonVisitor())
        .def(CopyableVisitor<Frame>())
        ;
    
        bp::class_< std::vector<Frame> >("StdVec_Frame")
        .def(bp::vector_indexing_suite< container::aligned_vector<Frame> >());
      }


    };
    

  } // namespace python
} // namespace se3

#endif // ifndef __se3_python_frame_hpp__
