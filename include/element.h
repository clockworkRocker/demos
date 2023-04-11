#ifndef LAB2_V2_INCLUDE_ELEMENT_H_
#define LAB2_V2_INCLUDE_ELEMENT_H_

#include <Eigen/Dense>
using MathMatrix = Eigen::MatrixXd;
using MathVector = Eigen::VectorXd;

#include "equation.h"

template <typename ODE>
class Element1D {
 protected:
  const double _l;
  const ODE& _eq;
 public:
  Element1D(double l, const ODE& ode);
  virtual ~Element1D() = default;
  double length() const;
  virtual unsigned nodes() const = 0;
  virtual double equation_coeff(unsigned index) const;
  virtual MathMatrix dirichlet() const = 0;
  virtual MathVector loads() const = 0;
};

enum ElemType {
  Linear1D, Cubic1D
};

template<unsigned Order>
class LinearElement1DL: public Element1D<LinearODE<Order>> {
 public:
  LinearElement1DL(double l, const LinearODE<Order> &ode);
  unsigned nodes() const override;
  MathMatrix dirichlet() const override;
  MathVector loads() const override;
};

template<unsigned Order>
class CubicElement1DL: public Element1D<LinearODE<Order>> {
 public:
  CubicElement1DL(double l, const LinearODE<Order> &ode);
  unsigned nodes() const override;
  MathMatrix dirichlet() const override;
  MathVector loads() const override;
};

#include "../src/element.cpp"

#endif //LAB2_V2_INCLUDE_ELEMENT_H_
/* ==================================================================================================== */

