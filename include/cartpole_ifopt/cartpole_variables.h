#pragma once
#include <ifopt/variable_set.h>
#include <ct/core/core.h>

using namespace ifopt;

namespace cartpole {

class StateVariables : public VariableSet 
{
public:
  const size_t N_Nodes;
  static const size_t STATE_DIM = 4;
  typedef ct::core::StateVector<STATE_DIM, double> state_vector_t;

  // Every variable set has a name, here "state-set". this allows the constraints
  // and costs to define values and Jacobians specifically w.r.t this variable set.
  StateVariables(size_t n_nodes, double d_max, double a_max, 
  state_vector_t target, const std::string& name = "state-set") : 
  VariableSet(n_nodes*STATE_DIM, name), N_Nodes(n_nodes), d_max_(d_max), a_max_(a_max), target_(target)
  {
    nodes_.resize(N_Nodes);
    // the initial values where the NLP starts iterating from
    for(size_t n = 0; n < N_Nodes; n++)
    {
      //subject to change for better convergence
      nodes_[n](0) = (n/N_Nodes) * target_(0);
      nodes_[n](1) = (n/N_Nodes) * target_(1);
      nodes_[n](2) = 0;
      nodes_[n](3) = 0;
    }
  }

  // Here is where you can transform the Eigen::Vector into whatever
  // internal representation of your variables
  void SetVariables(const VectorXd& x) override
  {
    size_t idx = 0;

    for(state_vector_t &node : nodes_)
    {
      for(size_t i = 0; i < STATE_DIM; i++)
      {
        node(i) = x(idx++);
      }
    }
  };

  // Here is the reverse transformation from the internal representation to
  // to the Eigen::Vector
  VectorXd GetValues() const override
  {
    VectorXd x(N_Nodes*STATE_DIM);
    size_t idx = 0;

    for(state_vector_t node : nodes_)
    {
      for(size_t i = 0; i < STATE_DIM; i++)
      {
        x(idx++) = node(i);
      }
    }
    
    return x;
  };

  // Each variable has an upper and lower bound set here
  VecBound GetBounds() const override
  {
    VecBound bounds(GetRows());
    size_t idx = 0;

    for(int i = 0; i < STATE_DIM; i++)
    {
      bounds[idx++] = Bounds(0, 0);
    }

    while(idx < STATE_DIM*(N_Nodes -1))
    {
      bounds[idx++] = Bounds(-d_max_, d_max_);
      bounds[idx++] = Bounds(-a_max_, a_max_);
      bounds[idx++] = NoBound;
      bounds[idx++] = NoBound;
    }

    for(int i = 0; i < STATE_DIM; i++)
    {
      bounds[idx++] = Bounds(target_(i), target_(i));
    }

    return bounds;
  }

private:
  std::vector<state_vector_t> nodes_;
  state_vector_t target_;
  double d_max_, a_max_;
};

class ControlVariables : public VariableSet 
{
public:
  const size_t N_Nodes;
  static const size_t CONTROL_DIM = 1;
  typedef ct::core::ControlVector<CONTROL_DIM, double> control_vector_t;

  // Every variable set has a name, here "control-set". this allows the constraints
  // and costs to define values and Jacobians specifically w.r.t this variable set.
  ControlVariables(size_t n_nodes, double f_max, const std::string& name = "control-set") : 
  VariableSet(n_nodes*CONTROL_DIM, name), N_Nodes(n_nodes), f_max_(f_max)
  {
    nodes_.resize(N_Nodes);
    // the initial values where the NLP starts iterating from
    for(size_t n = 0; n < N_Nodes; n++)
    {
      //subject to change for better convergence
      nodes_[n](0) = 0;
    }
  }

  // Here is where you can transform the Eigen::Vector into whatever
  // internal representation of your variables
  void SetVariables(const VectorXd& x) override
  {
    size_t idx = 0;

    for(control_vector_t &node : nodes_)
    {
      for(size_t i = 0; i < CONTROL_DIM; i++)
      {
        node(i) = x(idx++);
      }
    }
  };

  // Here is the reverse transformation from the internal representation to
  // to the Eigen::Vector
  VectorXd GetValues() const override
  {
    VectorXd x(N_Nodes*CONTROL_DIM);
    size_t idx = 0;

    for(control_vector_t node : nodes_)
    {
      for(size_t i = 0; i < CONTROL_DIM; i++)
      {
        x(idx++) = node(i);
      }
    }

    return x;
  };

  // Each variable has an upper and lower bound set here
  VecBound GetBounds() const override
  {
    VecBound bounds(GetRows());
    
    for(size_t idx = 0; idx < N_Nodes; idx++)
    {
      bounds[idx] = Bounds(-f_max_, f_max_);
    }

    return bounds;
  }

private:
  std::vector<control_vector_t> nodes_;
  double f_max_;
};

} //namespace cartpole