#include "element.h"

// -------- Base element class --------
template<typename ODE>
Element1D<ODE>::Element1D(double l, const ODE &ode): _l(l), _eq(ode) {}

template<typename ODE>
double Element1D<ODE>::length() const { return _l; }

template<>
double Element1D<LinearODE<2>>::equation_coeff(unsigned int index) const {
    assert(index >= 0 && index <= nodes());
    return _eq.coeffs[index];
}

// -------- Linear element --------
template<unsigned Order>
LinearElement1DL<Order>::LinearElement1DL(double l, const LinearODE<Order> &ode):
    Element1D<LinearODE<Order>>(l, ode) {}

template<unsigned Order>
unsigned int LinearElement1DL<Order>::nodes() const { return 2; }

// -------- Cubic element --------
template<unsigned Order>
CubicElement1DL<Order>::CubicElement1DL(double l, const LinearODE<Order> &ode):
    Element1D<LinearODE<Order>>(l, ode) {}

template<unsigned Order>
unsigned int CubicElement1DL<Order>::nodes() const { return 4; }

// -------- Exact instantation for second order --------
template<>
MathMatrix LinearElement1DL<2>::dirichlet() const {
    Eigen::Matrix2d _dirichletA, _dirichletB, _dirichletC;
    _dirichletA << (1 / _l), (-1 / _l), (-1 / _l), (1 / _l);
    _dirichletB << (1. / 2), (-1. / 2), (1. / 2), (-1. / 2);
    _dirichletC << (_l / 3), (_l / 6), (_l / 6), (_l / 3);
    _dirichletA *= -_eq.coeffs[0];
    _dirichletB *= -_eq.coeffs[1];
    _dirichletC *= _eq.coeffs[2];

    return _dirichletA + _dirichletB + _dirichletC;
}

template<>
MathVector LinearElement1DL<2>::loads() const {
    MathVector _loads(2);
    double d = _eq.coeffs[3];
    _loads << -d * _l / 2, -d * _l / 2;

    return _loads;
}

template<>
MathMatrix CubicElement1DL<2>::dirichlet() const {
    const double a1 = (-37 / (10 * _l));
    const double a2 = ((189 / (40 * _l)));
    const double a3 = (-27 / (20 * _l));
    const double a4 = (13 / (40 * _l));
    const double a5 = (-54 / (5 * _l));
    const double a6 = (297 / (40 * _l));
    Eigen::Matrix4d _dirichletA, _dirichletB, _dirichletC;
    _dirichletA << a1, a2, a3, a4,
        a2, a5, a6, a3,
        a3, a6, a5, a2,
        a4, a3, a2, a1;

    constexpr double b1 = -0.5;
    constexpr double b2 = 57. / 80;
    constexpr double b3 = -0.3;
    constexpr double b4 = 7. / 80;
    constexpr double b5 = 0.;
    constexpr double b6 = 81. / 80;
    _dirichletB << b1, b2, b3, b4,
        -b2, b5, b6, b3,
        -b3, -b6, -b5, -b2,
        -b4, -b3, -b2, -b1;

    _dirichletA *= _eq.coeffs[0];
    _dirichletB *= _eq.coeffs[1];

    return _dirichletA + _dirichletB;
}

template<>
MathVector CubicElement1DL<2>::loads() const {
    Eigen::Vector4d _loads;
    const double d1 = -_l / 8;
    const double d2 = -3. * _l / 8;
    _loads << d1, d2, d2, d1;
    _loads *= _eq.coeffs[3];

    return _loads;
}
