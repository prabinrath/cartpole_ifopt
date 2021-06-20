#pragma once
#include <ifopt/cost_term.h>
#include <ct/core/core.h>

using namespace ifopt;

namespace cartpole {

class CartpoleCost: public CostTerm {
public:
  const size_t N_Nodes;
  static const size_t CONTROL_DIM = 1;
  typedef ct::core::ControlVector<CONTROL_DIM, double> control_vector_t;

  CartpoleCost(size_t n_nodes, double dt, const std::string& name = "force-squared-cost") : CostTerm(name), N_Nodes(n_nodes), dt_(dt) {}

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

  double GetCost() const override
  {
    VectorXd control_set = GetVariables()->GetComponent("control-set")->GetValues();
    std::vector<control_vector_t> control_nodes(N_Nodes);
    GetControlNodes(control_set, control_nodes);

    double cost = 0;
    size_t idx = 0;
    cost += std::pow(control_nodes[idx++](0), 2);
    while(idx < CONTROL_DIM*(N_Nodes - 1))
    {
      cost += 2*std::pow(control_nodes[idx++](0), 2);
    }
    cost += std::pow(control_nodes[idx++](0), 2);
    return 0.5 * dt_ * cost;
  };

  void FillJacobianBlock (std::string var_set, Jacobian& jac) const override
  {
    if (var_set == "control-set") {
      VectorXd control_set = GetVariables()->GetComponent("control-set")->GetValues();
      std::vector<control_vector_t> control_nodes(N_Nodes);
      GetControlNodes(control_set, control_nodes);

      size_t idx = 0;
      
      jac.coeffRef(0, idx) = dt_ * control_nodes[idx](0);
      idx++;
      while(idx < CONTROL_DIM*(N_Nodes - 1))
      {
        jac.coeffRef(0, idx) = 2 * dt_ * control_nodes[idx](0);
        idx++;
      }
      jac.coeffRef(0, idx) = dt_ * control_nodes[idx](0);
      idx++;
    }
  }

private:
  double dt_;
};

} //namespace cartpole