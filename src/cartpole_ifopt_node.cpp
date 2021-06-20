#include <cmath>
#include <ros/package.h>
#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>
#include "cartpole_ifopt/cartpole.h"
#include "cartpole_ifopt/cartpole_variables.h"
#include "cartpole_ifopt/cartpole_constraints.h"
#include "cartpole_ifopt/cartpole_cost.h"

int main(int argc, char** argv)
{
  // create cartpole system
  double m1 = 5;  // (kg) Cart mass
  double m2 = 1;  // (kg) pole mass
  double l = 2;   // (m) pendulum (pole) length
  double dist = 0.8;  // How far must the cart translate during its swing-up
  double maxForce = 200;  // Maximum actuator forces
  double duration = 2;
  double dt = 0.1;

  size_t n_nodes = duration/dt;
  double d_max = 2*dist;
  double a_max = 2*M_PI;
  cartpole::CartpoleSystem::state_vector_t target;
  target << dist, M_PI, 0, 0;

  cartpole::CartpoleSystemPtr system(new cartpole::CartpoleSystem(m1, m2, l));

  // 1. define the problem
  ifopt::Problem nlp;
  nlp.AddVariableSet  (std::make_shared<cartpole::StateVariables>(n_nodes, d_max, a_max, target));
  nlp.AddVariableSet  (std::make_shared<cartpole::ControlVariables>(n_nodes, maxForce));
  nlp.AddConstraintSet(std::make_shared<cartpole::CartpoleConstraint>(n_nodes, dt, system));
  nlp.AddCostSet      (std::make_shared<cartpole::CartpoleCost>(n_nodes, dt));
  nlp.PrintCurrent();

  // 2. choose solver and options
  ifopt::IpoptSolver ipopt;
  ipopt.SetOption("linear_solver", "mumps");
  ipopt.SetOption("jacobian_approximation", "exact");
  // ipopt.SetOption("jacobian_approximation", "finite-difference-values");

  // 3 . solve
  ipopt.Solve(nlp);
  Eigen::VectorXd x = nlp.GetOptVariables()->GetValues();

  // print results for quick check on matlab
  std::cout << "X = [";
  size_t idx = 0;
  while(idx < 4 * n_nodes)
  {
    std::cout << x[idx++] << " ";
    std::cout << x[idx++] << " ";
    std::cout << x[idx++] << " ";
    std::cout << x[idx++] << ";";
  }
  std::cout << "]" << std::endl;
  std::cout << "U = [";
  while(idx < 5 * n_nodes)
  {
    std::cout << x[idx++] << " ";
  }
  std::cout << "]" << std::endl;
}