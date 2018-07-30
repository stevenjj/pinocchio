//
// Copyright (c) 2015-2018 CNRS
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

#ifndef __se3_jacobian_hxx__
#define __se3_jacobian_hxx__

#include "pinocchio/multibody/visitor.hpp"
#include "pinocchio/algorithm/check.hpp"

/// @cond DEV

namespace se3
{
  
  template<typename JointCollection, typename ConfigVectorType>
  struct JointJacobiansForwardStep
  : public fusion::JointVisitorBase< JointJacobiansForwardStep<JointCollection,ConfigVectorType> >
  {
    typedef ModelTpl<JointCollection> Model;
    typedef DataTpl<JointCollection> Data;
    
    typedef boost::fusion::vector<const Model &,
                                  Data &,
                                  const ConfigVectorType &> ArgsType;
    
    template<typename JointModel>
    static void algo(const JointModelBase<JointModel> & jmodel,
                     JointDataBase<typename JointModel::JointDataDerived> & jdata,
                     const Model & model,
                     Data & data,
                     const Eigen::MatrixBase<ConfigVectorType> & q)
    {
      typedef typename Model::JointIndex JointIndex;
      
      const JointIndex & i = (JointIndex) jmodel.id();
      const JointIndex & parent = model.parents[i];
      
      jmodel.calc(jdata.derived(),q);
      
      data.liMi[i] = model.jointPlacements[i]*jdata.M();
      if(parent>0) data.oMi[i] = data.oMi[parent]*data.liMi[i];
      else         data.oMi[i] = data.liMi[i];
      
      jmodel.jointCols(data.J) = data.oMi[i].act(jdata.S());
    }
  
  };
  
  template<typename JointCollection, typename ConfigVectorType>
  inline const typename DataTpl<JointCollection>::Matrix6x &
  computeJointJacobians(const ModelTpl<JointCollection> & model,
                        DataTpl<JointCollection> & data,
                        const Eigen::MatrixBase<ConfigVectorType> & q)
  {
    assert(model.check(data) && "data is not consistent with model.");
    assert(q.size() == model.nq && "The configuration vector is not of right size");
    
    typedef ModelTpl<JointCollection> Model;
    typedef typename Model::JointIndex JointIndex;
    
    typedef JointJacobiansForwardStep<JointCollection,ConfigVectorType> Pass;
    for(JointIndex i=1; i<(JointIndex) model.njoints; ++i)
    {
      Pass::run(model.joints[i],data.joints[i],
                typename Pass::ArgsType(model,data,q));
    }
  
    return data.J;
  }
  
  template<typename JointCollection>
  struct JointJacobiansForwardStep2
  : public fusion::JointVisitorBase< JointJacobiansForwardStep2<JointCollection> >
  {
    typedef ModelTpl<JointCollection> Model;
    typedef DataTpl<JointCollection> Data;
    
    typedef boost::fusion::vector<Data &> ArgsType;
    
    template<typename JointModel>
    static void algo(const JointModelBase<JointModel> & jmodel,
                     JointDataBase<typename JointModel::JointDataDerived> & jdata,
                     Data & data)
    {
      typedef typename Model::JointIndex JointIndex;
      
      const JointIndex & i = jmodel.id();
      jmodel.jointCols(data.J) = data.oMi[i].act(jdata.S());
    }
    
  };
  
  template<typename JointCollection>
  inline const typename DataTpl<JointCollection>::Matrix6x &
  computeJointJacobians(const ModelTpl<JointCollection> & model,
                        DataTpl<JointCollection> & data)
  {
    assert(model.check(data) && "data is not consistent with model.");
    
    typedef ModelTpl<JointCollection> Model;
    typedef typename Model::JointIndex JointIndex;
    
    typedef JointJacobiansForwardStep2<JointCollection> Pass;
    for(JointIndex i=1; i< (JointIndex)model.njoints; ++i)
    {
      Pass::run(model.joints[i],data.joints[i],
                typename Pass::ArgsType(data));
    }
    
    return data.J;
  }
  
  /* Return the jacobian of the output frame attached to joint <jointId> in the
   world frame or in the local frame depending on the template argument. The
   function computeJacobians should have been called first. */
  template<typename JointCollection, typename Matrix6Like>
  inline void getJointJacobian(const ModelTpl<JointCollection> & model,
                               const DataTpl<JointCollection> & data,
                               const typename ModelTpl<JointCollection>::JointIndex jointId,
                               const ReferenceFrame rf,
                               const Eigen::MatrixBase<Matrix6Like> & J)
  {
    assert( J.rows() == 6 );
    assert( J.cols() == model.nv );
    assert(model.check(data) && "data is not consistent with model.");
    
    typedef DataTpl<JointCollection> Data;
    
    Matrix6Like & J_ = EIGEN_CONST_CAST(Matrix6Like,J);
    
    const typename Data::SE3 & oMjoint = data.oMi[jointId];
    int colRef = nv(model.joints[jointId])+idx_v(model.joints[jointId])-1;
    for(int j=colRef;j>=0;j=data.parents_fromRow[(Model::Index)j])
    {
      if(rf == WORLD)   J_.col(j) = data.J.col(j);
      else              J_.col(j) = oMjoint.actInv(Motion(data.J.col(j))).toVector(); // TODO: use MotionRef
    }
  }
  
  template<typename JointCollection, typename ConfigVectorType, typename Matrix6Like>
  struct JointJacobianForwardStep
  : public fusion::JointVisitorBase< JointJacobianForwardStep<JointCollection,ConfigVectorType,Matrix6Like> >
  {
    typedef ModelTpl<JointCollection> Model;
    typedef DataTpl<JointCollection> Data;
    
    typedef boost::fusion::vector<const Model &,
                                  Data &,
                                  const ConfigVectorType &,
                                  Matrix6Like &
                                  > ArgsType;
    
    template<typename JointModel>
    static void algo(const JointModelBase<JointModel> & jmodel,
                     JointDataBase<typename JointModel::JointDataDerived> & jdata,
                     const Model & model,
                     Data & data,
                     const Eigen::MatrixBase<ConfigVectorType> & q,
                     const Eigen::MatrixBase<Matrix6Like> & J)
    {
      typedef typename Model::JointIndex JointIndex;
      const JointIndex & i = jmodel.id();
      const JointIndex & parent = model.parents[i];
      
      jmodel.calc(jdata.derived(),q);
      
      data.liMi[i] = model.jointPlacements[i]*jdata.M();
      data.iMf[parent] = data.liMi[i]*data.iMf[i];
      
      Matrix6Like & J_ = EIGEN_CONST_CAST(Matrix6Like,J);
      jmodel.jointCols(J_) = data.iMf[i].inverse().act(jdata.S()); // TODO: use MotionRef
    }
  
  };
  
  template<typename JointCollection, typename ConfigVectorType, typename Matrix6Like>
  inline void jointJacobian(const ModelTpl<JointCollection> & model,
                            DataTpl<JointCollection> & data,
                            const Eigen::MatrixBase<ConfigVectorType> & q,
                            const JointIndex jointId,
                            const Eigen::MatrixBase<Matrix6Like> & J)
  {
    assert(model.check(data) && "data is not consistent with model.");
    assert(q.size() == model.nq && "The configuration vector is not of right size");
    
    typedef ModelTpl<JointCollection> Model;
    typedef typename Model::JointIndex JointIndex;
    
    data.iMf[jointId].setIdentity();
    typedef JointJacobianForwardStep<JointCollection,ConfigVectorType,Matrix6Like> Pass;
    for(JointIndex i=jointId; i>0; i=model.parents[i])
    {
      Pass::run(model.joints[i],data.joints[i],
                typename Pass::ArgsType(model,data,q,EIGEN_CONST_CAST(Matrix6Like,J)));
    }
  }
  
  template<typename JointCollection, typename ConfigVectorType, typename TangentVectorType>
  struct JointJacobiansTimeVariationForwardStep
  : public fusion::JointVisitorBase< JointJacobiansTimeVariationForwardStep<JointCollection,ConfigVectorType,TangentVectorType> >
  {
    typedef ModelTpl<JointCollection> Model;
    typedef DataTpl<JointCollection> Data;
    
    typedef boost::fusion::vector<const Model &,
                                  Data &,
                                  const ConfigVectorType &,
                                  const TangentVectorType &> ArgsType;
    
    template<typename JointModel>
    static void algo(const JointModelBase<JointModel> & jmodel,
                     JointDataBase<typename JointModel::JointDataDerived> & jdata,
                     const Model & model,
                     Data & data,
                     const Eigen::MatrixBase<ConfigVectorType> & q,
                     const Eigen::MatrixBase<TangentVectorType> & v)
    {
      typedef typename Model::JointIndex JointIndex;
      typedef typename Data::SE3 SE3;
      typedef typename Data::Motion Motion;
      
      const JointIndex & i = (JointIndex) jmodel.id();
      const JointIndex & parent = model.parents[i];
      
      SE3 & oMi = data.oMi[i];
      Motion & vJ = data.v[i];
      
      jmodel.calc(jdata.derived(),q,v);
      
      vJ = jdata.v();
      
      data.liMi[i] = model.jointPlacements[i]*jdata.M();
      if(parent>0)
      {
        oMi = data.oMi[parent]*data.liMi[i];
        vJ += data.liMi[i].actInv(data.v[parent]);
      }
      else
      {
        oMi = data.liMi[i];
      }
      
      jmodel.jointCols(data.J) = oMi.act(jdata.S());
      
      // Spatial velocity of joint i expressed in the global frame o
      data.ov[i] = oMi.act(vJ);
      
      typedef typename SizeDepType<JointModel::NV>::template ColsReturn<typename Data::Matrix6x>::Type ColsBlock;
      ColsBlock dJcols = jmodel.jointCols(data.dJ);
      ColsBlock Jcols = jmodel.jointCols(data.J);
      
      motionSet::motionAction(data.ov[i],Jcols,dJcols);
    }
    
  };
  
  template<typename JointCollection, typename ConfigVectorType, typename TangentVectorType>
  inline const typename DataTpl<JointCollection>::Matrix6x &
  computeJointJacobiansTimeVariation(const ModelTpl<JointCollection> & model,
                                     DataTpl<JointCollection> & data,
                                     const Eigen::MatrixBase<ConfigVectorType> & q,
                                     const Eigen::MatrixBase<TangentVectorType> & v)
  {
    assert(model.check(data) && "data is not consistent with model.");
    assert(q.size() == model.nq && "The configuration vector is not of right size");
    assert(v.size() == model.nv && "The velocity vector is not of right size");
    
    typedef ModelTpl<JointCollection> Model;
    typedef typename Model::JointIndex JointIndex;
    
    typedef JointJacobiansTimeVariationForwardStep<JointCollection,ConfigVectorType,TangentVectorType> Pass;
    for(JointIndex i=1; i<(JointIndex)model.njoints; ++i)
    {
      Pass::run(model.joints[i],data.joints[i],
                typename Pass::ArgsType(model,data,q,v));
    }
    
    return data.dJ;
  }
  
  template<typename JointCollection, typename Matrix6Like>
  inline void getJointJacobianTimeVariation(const ModelTpl<JointCollection> & model,
                                            const DataTpl<JointCollection> & data,
                                            const JointIndex jointId,
                                            const ReferenceFrame rf,
                                            const Eigen::MatrixBase<Matrix6Like> & dJ)
  {
    assert( dJ.rows() == 6 );
    assert( dJ.cols() == model.nv );
    assert(model.check(data) && "data is not consistent with model.");
    
    typedef DataTpl<JointCollection> Data;
    
    Matrix6Like & dJ_ = EIGEN_CONST_CAST(Matrix6Like,dJ);
    
    const typename Data::SE3 & oMjoint = data.oMi[jointId];
    int colRef = nv(model.joints[jointId])+idx_v(model.joints[jointId])-1;
    for(int j=colRef;j>=0;j=data.parents_fromRow[(size_t)j])
    {
      if(rf == WORLD)   dJ_.col(j) = data.dJ.col(j);
      else              dJ_.col(j) = oMjoint.actInv(Motion(data.dJ.col(j))).toVector(); // TODO: use MotionRef
    }
  }
  
  
} // namespace se3

/// @endcond

#endif // ifndef __se3_jacobian_hxx__
