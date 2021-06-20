#pragma once
#include <cmath>
#include <memory>
#include <iostream>
#include <ct/core/core.h>  // as usual, include CT

namespace cartpole {
namespace tpl {

template <typename SCALAR>
class CartpoleSystem : public ct::core::ControlledSystem<4, 1, SCALAR>
{
public:
    static const size_t STATE_DIM = 4;
    static const size_t CONTROL_DIM = 1;

    typedef ct::core::ControlVector<CONTROL_DIM, SCALAR> control_vector_t;
    typedef ct::core::StateVector<STATE_DIM, SCALAR> state_vector_t;
    typedef ct::core::tpl::TraitSelector<SCALAR> trait;

    CartpoleSystem() = delete;

    CartpoleSystem(SCALAR m1, SCALAR m2, SCALAR l, 
    std::shared_ptr<ct::core::Controller<STATE_DIM, CONTROL_DIM, SCALAR>> controller = nullptr)
    :ct::core::ControlledSystem<STATE_DIM, CONTROL_DIM, SCALAR>(controller, ct::core::SYSTEM_TYPE::GENERAL),
    m1_(m1), m2_(m2), l_(l)
    {
    }
    CartpoleSystem(const CartpoleSystem& arg) 
    : ct::core::ControlledSystem<STATE_DIM, CONTROL_DIM, SCALAR>(arg),
    m1_(arg.m1_), m2_(arg.m2_), l_(arg.l_)
    {        
    }
    virtual ~CartpoleSystem() {}
    CartpoleSystem* clone() const override { return new CartpoleSystem(*this); }

    virtual void computeControlledDynamics(const state_vector_t& state,
        const SCALAR& t,
        const control_vector_t& control,
        state_vector_t& derivative) override
    {
        derivative(0) = state(2);
        derivative(1) = state(3);
        derivative(2) = (l_*m2_*trait::Trait::sin(state(1))*state(3)*state(3) + control(0) + m2_*g_*trait::Trait::cos(state(1))*trait::Trait::sin(state(1))) / 
                        (m1_ + m2_*(1-trait::Trait::cos(state(1))*trait::Trait::cos(state(1))));
        derivative(3) = -1 * (l_*m2_*trait::Trait::cos(state(1))*trait::Trait::sin(state(1))*state(3)*state(3) + control(0)*trait::Trait::cos(state(1)) + (m1_+m2_)*g_*trait::Trait::sin(state(1))) /
                        (l_*m1_ + l_*m2_*(1-trait::Trait::cos(state(1))*trait::Trait::cos(state(1))));
        // std::cout << "State" << state << std::endl;
        // std::cout << "Derivative" << derivative << std::endl;
    }

private:
    SCALAR m1_;
    SCALAR m2_;
    SCALAR l_;
    const SCALAR g_ = SCALAR(9.81);
};
}

typedef tpl::CartpoleSystem<double> CartpoleSystem;
typedef std::shared_ptr<CartpoleSystem> CartpoleSystemPtr;

}  // namespace cartpole