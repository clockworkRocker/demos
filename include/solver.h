#ifndef LAB2_V2_INCLUDE_SOLVER_H_
#define LAB2_V2_INCLUDE_SOLVER_H_

#include "element.h"

#include <Eigen/Sparse>
#include <iostream>
#include <utility>
#include <vector>

using SparseMatrix = Eigen::SparseMatrix<double>;
template <typename T>
using Interval = std::pair<T, T>;

enum boundary_t {
  FIRST, SECOND
};

struct BoundaryCondition {
  boundary_t type;
  double x;
  double value;

  BoundaryCondition(boundary_t new_type, double new_x, double new_value);
};

template<typename ODE>
class Solver1D {
 protected:
  Interval<double> _interval = {0, 0};
  unsigned int _total_elements = 0;
  unsigned int _total_nodes = 0;

  MathMatrix _dirichlet;
  MathVector _loads;
  MathVector _unknowns;

  const ElemType _element_type;
  Element1D<ODE> *_element;
  std::ostream & _output;

  Element1D<ODE>* _make_element(const ODE& equation, double length);
  void _assemble();
  void _apply_boundaries(const std::vector<BoundaryCondition>& conditions);
  void _apply_first(const BoundaryCondition& condition);
  void _apply_second(const BoundaryCondition& condition);

 public:
  explicit Solver1D(ElemType elem_type, std::ostream& os = std::cout);
  int solve(const ODE &equation, const std::vector<BoundaryCondition> &conditions, const Interval<double> &interval, unsigned total_elems);
  void print_unknowns(std::ostream &stream = std::cout);
};

#include "../src/solver.cpp"
#endif //LAB2_V2_INCLUDE_SOLVER_H_
