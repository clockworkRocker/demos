#include "solver.h"

#include <Eigen/SparseCholesky>

BoundaryCondition::BoundaryCondition(boundary_t new_type, double new_x, double new_value):
    type(new_type),
    x(new_x),
    value(new_value)
    {}

template<typename ODE>
Solver1D<ODE>::Solver1D(ElemType elem_type, std::ostream &os): _element_type(elem_type), _output(os) {}

template<typename ODE>
void Solver1D<ODE>::_assemble() {
    const unsigned int dim = (_element->nodes() - 1) * _total_elements + 1;
    _total_nodes = dim;
    _dirichlet = MathMatrix::Zero(dim, dim);
    _loads = MathVector::Zero(dim);
    _unknowns = MathVector::Zero(dim);

    for (unsigned i = 0; i < _total_elements; ++i) {
        const unsigned pos = (_element->nodes() - 1) * i;
        _dirichlet.block(pos, pos, _element->nodes(), _element->nodes()) += _element->dirichlet();
        _loads.segment(pos, _element->nodes()) += _element->loads();
    }
}

template<typename ODE>
void Solver1D<ODE>::_apply_first(const BoundaryCondition &condition) {
    unsigned index = (condition.x - _interval.first) / (_interval.second - _interval.first) * (_total_nodes - 1);
    assert(index >= 0 && index < _dirichlet.rows());

    _dirichlet.row(index).setZero();
    //_dirichlet.col(index).setZero();
    _dirichlet(index, index) = 1;
    _loads(index) = condition.value;
}

template<typename ODE>
void Solver1D<ODE>::_apply_second(const BoundaryCondition &condition) {
    unsigned index = (condition.x - _interval.first) / (_interval.second - _interval.first) * (_total_nodes - 1);
    assert(index >= 0 && index < _dirichlet.rows());

    _loads(index) -= _element->equation_coeff(0) * condition.value;
}

template<typename ODE>
void Solver1D<ODE>::_apply_boundaries(const std::vector<BoundaryCondition> &conditions) {
    for (auto &condition: conditions)
        switch (condition.type) {
            case FIRST: _apply_first(condition);
                break;
            case SECOND: _apply_second(condition);
                break;
            default: break;
        }
}

template<>
Element1D<LinearODE<2>> *Solver1D<LinearODE<2>>::_make_element(const LinearODE<2> &equation, double length) {
    switch (_element_type) {
        case Linear1D: return new LinearElement1DL<2>(length, equation);
        case Cubic1D: return new CubicElement1DL<2>(length, equation);
    }

    return nullptr;
}

template<typename ODE>
void Solver1D<ODE>::print_unknowns(std::ostream &stream) {
    stream << "# X\tU\n";
    const double start = _interval.first;
    const double step = (_interval.second - start) / _total_nodes;

    for (int i = 0; i < _unknowns.size(); ++i)
        stream << start + i * step
               << '\t' << _unknowns[i]
               << '\n';
}

template<typename ODE>
int Solver1D<ODE>::solve(const ODE &equation,
                         const std::vector<BoundaryCondition> &conditions,
                         const Interval<double> &interval,
                         unsigned int total_elems) {
    _total_elements = total_elems;
    _interval = interval;
    _element = _make_element(equation, (interval.second - interval.first) / total_elems);
    _assemble();
    _apply_boundaries(conditions);

    Eigen::ColPivHouseholderQR<MathMatrix> decomposition(_dirichlet);
    _unknowns = decomposition.template solve(_loads);

    print_unknowns(_output);

    delete _element;
    _element = nullptr;

    return 0;
}