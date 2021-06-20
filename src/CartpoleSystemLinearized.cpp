/**********************************************************************************************************************
This file is part of the Control Toolbox (https://github.com/ethz-adrl/control-toolbox), copyright by ETH Zurich.
Licensed under the BSD-2 license (see LICENSE file in main directory)
**********************************************************************************************************************/

// clang-format off

#include "cartpole_ifopt/CartpoleSystemLinearized.h"

namespace ct {
namespace core {
namespace generated {


const CartpoleSystemLinearized::state_matrix_t& CartpoleSystemLinearized::getDerivativeState(
    const state_vector_t& x,
    const control_vector_t& u,
    const double t)
{
    double* jac = dFdx_.data();
    Eigen::Matrix<double, 4 + 1, 1> x_in;
    x_in << x, u;

        vX_[0] = cos(x_in[1]);
    vX_[1] = cos(x_in[1]);
    vX_[2] = 9.81 * vX_[1];
    vX_[3] = -1 * sin(x_in[1]);
    vX_[4] = sin(x_in[1]);
    vX_[5] = 2. * vX_[4];
    vX_[6] = vX_[5] * x_in[3];
    vX_[7] = vX_[1] * vX_[1];
    vX_[8] = 6. - vX_[7];
    vX_[9] = vX_[1] * vX_[3] + vX_[3] * vX_[1];
    jac[6] = (2. * vX_[0] * x_in[3] * x_in[3] + vX_[2] * vX_[0] + 9.81 * vX_[3] * vX_[4] - (x_in[4] + vX_[6] * x_in[3] + vX_[2] * vX_[4]) / vX_[8] * (0 - vX_[9])) / vX_[8];
    vX_[2] = 2. * vX_[1];
    vX_[10] = vX_[2] * vX_[4];
    vX_[11] = vX_[10] * x_in[3];
    vX_[7] = 10. + 2. * (1 - vX_[7]);
    jac[7] = (-1 * (x_in[4] * vX_[3] + (vX_[2] * vX_[0] + 2. * vX_[3] * vX_[4]) * x_in[3] * x_in[3] + 58.86 * vX_[0]) - (-1 * (x_in[4] * vX_[1] + vX_[11] * x_in[3] + 58.86 * vX_[4])) / vX_[7] * 2. * (- vX_[9])) / vX_[7];
    jac[14] = (vX_[6] + vX_[5] * x_in[3]) / vX_[8];
    jac[15] = (-1 * (vX_[11] + vX_[10] * x_in[3])) / vX_[7];
    // dependent variables without operations
    jac[8] = 1;
    jac[13] = 1;


    return dFdx_;
}

const CartpoleSystemLinearized::state_control_matrix_t& CartpoleSystemLinearized::getDerivativeControl(
    const state_vector_t& x,
    const control_vector_t& u,
    const double t)
{
    double* jac = dFdu_.data();
    Eigen::Matrix<double, 4 + 1, 1> x_in;
    x_in << x, u;

        vU_[0] = cos(x_in[1]);
    vU_[1] = vU_[0] * vU_[0];
    jac[2] = 1 / (6. - vU_[1]);
    jac[3] = (-1 * vU_[0]) / (10. + 2. * (1 - vU_[1]));


    return dFdu_;
}

} // namespace generated
} // namespace core
} // namespace ct

// clang-format on
