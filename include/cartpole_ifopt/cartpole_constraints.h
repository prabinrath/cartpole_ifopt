#pragma once
#include <ifopt/constraint_set.h>
#include "cartpole_ifopt/cartpole.h"
#include "cartpole_ifopt/CartpoleSystemLinearized.h"

using namespace ifopt;

namespace cartpole {

class CartpoleConstraint : public ConstraintSet {
public:
  const size_t N_Nodes;
  static const size_t STATE_DIM = 4;
  static const size_t CONTROL_DIM = 1;
  typedef ct::core::ControlVector<CONTROL_DIM, double> control_vector_t;
  typedef ct::core::StateVector<STATE_DIM, double> state_vector_t;
  typedef ct::core::generated::CartpoleSystemLinearized CartpoleSystemLinearized;

  // This constraint set contains multiple related constraints.
  CartpoleConstraint(size_t n_nodes, double dt, CartpoleSystemPtr system, const std::string& name = "dynamics-constraint") : 
  ConstraintSet((n_nodes-1)*STATE_DIM, name), system_(system), dt_(dt), N_Nodes(n_nodes),
  system_jac_(new CartpoleSystemLinearized()) {}

  void GetStateNodes(VectorXd &state_set, std::vector<state_vector_t> &state_nodes) const
  {
    size_t idx = 0;
    for(state_vector_t &node : state_nodes)
    {
      for(size_t i = 0; i < STATE_DIM; i++)
      {
        node(i) = state_set(idx++);
      }
    }
  }

  void GetControlNodes(VectorXd &control_set, std::vector<control_vector_t> &control_nodes) const
  {
    size_t idx = 0;
    for(control_vector_t &node : control_nodes)
    {
      for(size_t i = 0; i < CONTROL_DIM; i++)
      {
        node(i) = control_set(idx++);
      }
    }
  }

  // The constraint value minus the constant value "c", moved to bounds.
  VectorXd GetValues() const override
  {
    VectorXd g(GetRows());
    
    VectorXd state_set = GetVariables()->GetComponent("state-set")->GetValues();
    std::vector<state_vector_t> state_nodes(N_Nodes);
    GetStateNodes(state_set, state_nodes);

    VectorXd control_set = GetVariables()->GetComponent("control-set")->GetValues();
    std::vector<control_vector_t> control_nodes(N_Nodes);
    GetControlNodes(control_set, control_nodes);

    size_t row_count = 0;
    for(size_t idx = 0; idx < N_Nodes - 1; idx++)
    {
      state_vector_t derivative0, derivative1;
      system_->computeControlledDynamics(state_nodes[idx], 0, control_nodes[idx], derivative0);
      system_->computeControlledDynamics(state_nodes[idx+1], 0, control_nodes[idx+1], derivative1);
      state_vector_t constraint = state_nodes[idx+1] - (state_nodes[idx] + 0.5 * dt_ * (derivative0 + derivative1));
      for(int i = 0; i < STATE_DIM; i++)
      {
        g(row_count++) = constraint(i);
      }
    }
    return g;
  };

  // The only constraint in this set is an equality constraint to "c".
  // Constant values should always be put into GetBounds(), not GetValues().
  // For inequality constraints (<,>), use Bounds(x, inf) or Bounds(-inf, x).
  VecBound GetBounds() const override
  {
    VecBound bounds(GetRows());
    for(Bounds &b : bounds)
    {
      b = Bounds(0, 0);
    }
    return bounds;
  }

  // This function provides the first derivative of the constraints.
  // In case this is too difficult to write, you can also tell the solvers to
  // approximate the derivatives by finite differences and not overwrite this
  // function, e.g. in ipopt.cc::use_jacobian_approximation_ = true
  // Attention: see the parent class function for important information on sparsity pattern.
  void FillJacobianBlock (std::string var_set, Jacobian& jac_block) const override
  {
    VectorXd state_set = GetVariables()->GetComponent("state-set")->GetValues();
    std::vector<state_vector_t> state_nodes(N_Nodes);
    GetStateNodes(state_set, state_nodes);

    VectorXd control_set = GetVariables()->GetComponent("control-set")->GetValues();
    std::vector<control_vector_t> control_nodes(N_Nodes);
    GetControlNodes(control_set, control_nodes);

    // must fill only that submatrix of the overall Jacobian that relates
    // to this constraint. 
    if (var_set == "state-set") {
      for(size_t idx = 0; idx < N_Nodes - 1; idx++)
      {
        CartpoleSystemLinearized::state_matrix_t A0 = system_jac_->getDerivativeState(state_nodes[idx], control_nodes[idx]);
        CartpoleSystemLinearized::state_matrix_t A1 = system_jac_->getDerivativeState(state_nodes[idx+1], control_nodes[idx+1]);
        
        CartpoleSystemLinearized::state_matrix_t jac_constraint0 = 
        -1 * (Eigen::MatrixXd::Identity(STATE_DIM, STATE_DIM) + 0.5 * dt_ * A0);
        CartpoleSystemLinearized::state_matrix_t jac_constraint1 = 
        Eigen::MatrixXd::Identity(STATE_DIM, STATE_DIM) - 0.5 * dt_ * A1;
        
        for(size_t jdx = idx*STATE_DIM, i = 0; jdx < (idx+1)*STATE_DIM; jdx++, i++)
        {
          for(size_t kdx = idx*STATE_DIM, j = 0; kdx < (idx+1)*STATE_DIM; kdx++, j++)
          {
            jac_block.coeffRef(jdx, kdx) = jac_constraint0(i, j);
          }          
        }

        for(size_t jdx = idx*STATE_DIM, i = 0; jdx < (idx+1)*STATE_DIM; jdx++, i++)
        {
          for(size_t kdx = (idx+1)*STATE_DIM, j = 0; kdx < (idx+2)*STATE_DIM; kdx++, j++)
          {
            jac_block.coeffRef(jdx, kdx) = jac_constraint1(i, j);
          }          
        }
      }
    }

    if (var_set == "control-set") {
      for(size_t idx = 0; idx < N_Nodes - 1; idx++)
      {
        CartpoleSystemLinearized::state_control_matrix_t B0 = system_jac_->getDerivativeControl(state_nodes[idx], control_nodes[idx]);
        CartpoleSystemLinearized::state_control_matrix_t B1 = system_jac_->getDerivativeControl(state_nodes[idx+1], control_nodes[idx+1]);
        CartpoleSystemLinearized::state_control_matrix_t jac_constraint0 = -1 * (0.5 * dt_ * B0);        
        CartpoleSystemLinearized::state_control_matrix_t jac_constraint1 = -1 * (0.5 * dt_ * B1);

        for(size_t jdx = idx*STATE_DIM, i = 0; jdx < (idx+1)*STATE_DIM; jdx++, i++)
        {
          for(size_t kdx = idx*CONTROL_DIM, j = 0; kdx < (idx+1)*CONTROL_DIM; kdx++, j++)
          {
            jac_block.coeffRef(jdx, kdx) = jac_constraint0(i, j);
          }          
        }

        for(size_t jdx = idx*STATE_DIM, i = 0; jdx < (idx+1)*STATE_DIM; jdx++, i++)
        {
          for(size_t kdx = (idx+1)*CONTROL_DIM, j = 0; kdx < (idx+2)*CONTROL_DIM; kdx++, j++)
          {
            jac_block.coeffRef(jdx, kdx) = jac_constraint1(i, j);
          }          
        }
      }
    }
  }

private:
  CartpoleSystemPtr system_;
  std::shared_ptr<CartpoleSystemLinearized> system_jac_;
  double dt_;
};

} //namespace cartpole