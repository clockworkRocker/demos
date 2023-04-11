#include "thermal.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <cstring>

#include <omp.h>

/* ------------------------------------------------------ NODE
 * ------------------------------------------------------ */
ThermoNode2d::ThermoNode2d() : Point2d(), temperature(0. / 0.) {}

ThermoNode2d::ThermoNode2d(double x, double y, double val)
    : Point2d(x, y), temperature(val) {}

/* ------------------------------------------------------ MESH
 * ------------------------------------------------------ */
MFDMesh2d::MFDMesh2d(Shape2d& geometry, double dx, double dy)
    : _dx(dx),
      _dy(dy),
      _rows((size_t)((geometry.maxY() - geometry.minY()) / dy) + 1),
      _cols((size_t)((geometry.maxX() - geometry.minX()) / dx) + 1),
      _nodes(new ThermoNode2d[_rows * _cols]),
      _row_offsets(new size_t[_rows]),
      _row_ends(new size_t[_rows]) {
  double startX = geometry.minX();
  double startY = geometry.minY();

  for (int i = 0; i < _rows; ++i) {
    bool foundLeftEdge = false;
    bool foundRightEdge = false;
    _row_offsets[i] = 0;
    _row_ends[i] = _cols;
    for (int j = 0; j < _cols; ++j) {
      at(i, j).x = startX + j * dx;
      at(i, j).y = startY + i * dy;
      at(i, j).temperature = 0;
      if (!foundLeftEdge) {
        if (!geometry.hasPoint(at(i, j)))
          _row_offsets[i]++;
        else
          foundLeftEdge = true;
      } else {
        if (!foundRightEdge && !geometry.hasPoint(at(i, j))) {
          foundRightEdge = true;
          _row_ends[i] = j;
        }
      }
    }
  }
}

MFDMesh2d::MFDMesh2d(const MFDMesh2d& rhs)
    : _dx(rhs._dx),
      _dy(rhs._dy),
      _rows(rhs._rows),
      _cols(rhs._cols),
      _nodes(new ThermoNode2d[_rows * _cols]),
      _row_offsets(new size_t[_rows]),
      _row_ends(new size_t[_rows]) {
  memcpy(_nodes, rhs._nodes, _rows * _cols * sizeof(ThermoNode2d));
  memcpy(_row_offsets, rhs._row_offsets, _rows * sizeof(double));
  memcpy(_row_ends, rhs._row_ends, _rows * sizeof(double));
}

MFDMesh2d::~MFDMesh2d() {
  delete[] _nodes;
  delete[] _row_offsets;
  delete[] _row_ends;
}

double MFDMesh2d::dx() const {
  return _dx;
}

double MFDMesh2d::dy() const {
  return _dy;
}

size_t MFDMesh2d::rows() const {
  return _rows;
}

size_t MFDMesh2d::rowStart(size_t row) {
  return _row_offsets[row];
}

size_t MFDMesh2d::rowEnd(size_t row) {
  return _row_ends[row];
}

ThermoNode2d& MFDMesh2d::operator[](size_t index) {
  return _nodes[index];
}

ThermoNode2d& MFDMesh2d::at(size_t row, size_t col) {
  return _nodes[row * _cols + col];
}

bool MFDMesh2d::indexOfNode(double x, double y, size_t* index) {
  for (int i = 0; i < _rows; ++i)
    for (int j = 0; j < _cols; ++j)
      if (at(i, j).x == x && at(i, j).y == y) {
        if (index)
          *index = i * _cols + j;
        return true;
      }

  return false;
}

void MFDMesh2d::debugPrint() {
  for (int i = _rows - 1; i >= 0; --i) {
    for (int j = 0; j < _row_ends[i]; ++j) {
      if (j < _row_offsets[i])
        std::cout << "  ";
      else
        std::cout << "* ";
    }
    std::cout << std::endl;
  }
}

/* ------------------------------------------------- THERMAL SHAPE
 * -------------------------------------------------- */
ThermoShape2d::ThermoShape2d(std::vector<Line2d*> edges)
    : Shape2d(edges), _mesh(nullptr) {}

ThermoShape2d::~ThermoShape2d() {
  delete _mesh;
}

ThermoShape2d::BoundaryCondition::BoundaryCondition(size_t type,
                                                    size_t geom_type,
                                                    size_t geom_index,
                                                    double value)
    : type(type),
      geometry_type(geom_type),
      geometry_index(geom_index),
      value(value) {}

void ThermoShape2d::generateMesh(double dx, double dy) {
  _mesh = new MFDMesh2d(*this, dx, dy);
}

bool ThermoShape2d::lookupNode(double x, double y, size_t* index) {
  if (_mesh)
    return _mesh->indexOfNode(x, y, index);

  return false;
}

bool ThermoShape2d::addBoundary(size_t geometry_type,
                                size_t geometry_index,
                                size_t type,
                                double value) {
  _boundary_conditions.emplace_back(type, geometry_type, geometry_index, value);
  return true;
}

/* ---------------------------------------------------- SOLVER
 * ------------------------------------------------------ */
MFDSolver2d::~MFDSolver2d() {
  delete _next_step;
}

double MFDSolver2d::step() const {
  return _time_step;
}

bool MFDSolver2d::solve(ThermoShape2d& shape,
                        double total_time,
                        std::fstream& output_file) {
  MFDMesh2d* tmp = nullptr;

  _task = &shape;
  if (!_task->_mesh)
    return false;

  _next_step = new MFDMesh2d(*_task->_mesh);
  _time_step = _task->_mesh->_dx * _task->_mesh->_dy / 5;

#pragma omp parallel num_threads(2)
  {
#pragma omp single
    std::cout << "Hallo, I am using " << omp_get_num_threads() << " threads\n";
    for (int i = 0; i < (int)(total_time / _time_step) + 1; ++i) {
      _applyBoundaries();

#pragma omp for collapse(2)
      for (size_t j = 1; j < _next_step->_rows - 1; ++j)
        for (size_t k = 1; k < _next_step->_cols - 1; ++k)
          _calculateNode(j, k);

#pragma omp single
      {
        tmp = _task->_mesh;
        _task->_mesh = _next_step;
        _next_step = tmp;
      }
    }
  }
  _printTimeLayerForPlot(output_file);
  return true;
}

void MFDSolver2d::_applyBoundaries() {
  for (auto& boundary : _task->_boundary_conditions) {
    switch (boundary.geometry_type) {
      case GM_NODE:
#pragma omp single
      {
        _applyBoundaryToNode(boundary);
      } break;
      case GM_EDGE:
        _applyBoundaryToEdge(boundary);
        break;
    }
  }
}

void MFDSolver2d::_applyBoundaryToNode(
    ThermoShape2d::BoundaryCondition& condition) {
  size_t i = condition.geometry_index / _task->_mesh->_cols;
  size_t j = condition.geometry_index % _task->_mesh->_cols;

  switch (condition.type) {
    case TM_TEMPERATURE:
      _task->_mesh->at(i, j).temperature = condition.value;
      break;
    case TM_HEAT_FLUX:
      if (j < _task->_mesh->_row_ends[i] - 1)
        _task->_mesh->at(i, j).temperature =
            _task->_mesh->at(i, j + 1).temperature -
            condition.value * _task->_mesh->_dx;
      else
        _task->_mesh->at(i, j).temperature =
            _task->_mesh->at(i, j - 1).temperature +
            condition.value * _task->_mesh->_dx;
      break;
    case TM_CONVECTION:
      double n = sqrt(_next_step->_dx * _next_step->_dx +
                      _next_step->_dy * _next_step->_dy);
      double nx = _next_step->_dx / n;
      double ny = _next_step->_dy / n;
      double num =
          condition.value -
          _task->_mesh->at(i, j - 1).temperature * nx / _task->_mesh->_dx +
          _task->_mesh->at(i + 1, j).temperature * ny / _task->_mesh->_dy;
      double den = -nx / _task->_mesh->_dx + ny / _task->_mesh->_dy + 1;

      _task->_mesh->at(i, j).temperature =
          (_next_step->at(i - 1, j - 1).temperature + condition.value * n) /
          (1 + n);
      break;
  }
}

void MFDSolver2d::_applyBoundaryToEdge(
    ThermoShape2d::BoundaryCondition& condition) {
  Line2d* edge = _task->edge(condition.geometry_index);
  double precision = _task->_mesh->_dx * (edge->maxY() - edge->minY()) * 0.5;

#pragma omp for collapse(2)
  for (size_t i = 0; i < _task->_mesh->_rows; ++i) {
    for (size_t j = 0; j < _task->_mesh->_cols; ++j)
      if (edge->hasPoint(_task->_mesh->at(i, j), precision))
        switch (condition.type) {
          case TM_TEMPERATURE:
            _task->_mesh->at(i, j).temperature = condition.value;
            break;
          case TM_HEAT_FLUX:
            _task->_mesh->at(i, j).temperature =
                _task->_mesh->_dx * condition.value +
                _task->_mesh->at(i, j + 1).temperature;
            break;
        };
  }
}

double MFDSolver2d::_d2dx2(size_t i, size_t j) {
  double next = _task->_mesh->at(i, j + 1).temperature;
  double curr = _task->_mesh->at(i, j).temperature;
  double prev = _task->_mesh->at(i, j - 1).temperature;
  double dx = _task->_mesh->_dx;

  return (next - 2 * curr + prev) / (dx * dx);
}

double MFDSolver2d::_d2dy2(size_t i, size_t j) {
  double next = _task->_mesh->at(i + 1, j).temperature;
  double curr = _task->_mesh->at(i, j).temperature;
  double prev = _task->_mesh->at(i - 1, j).temperature;
  double dy = _task->_mesh->_dy;

  return (next - 2 * curr + prev) / (dy * dy);
}

void MFDSolver2d::_calculateNode(size_t i, size_t j) {
  _next_step->at(i, j).temperature = _task->_mesh->at(i, j).temperature +
                                     _time_step * (_d2dx2(i, j) + _d2dy2(i, j));
}

void MFDSolver2d::_calculateTimeLayer() {
#pragma omp for
  for (size_t i = 1; i < _next_step->_rows - 1; ++i)
    for (size_t j = _next_step->rowStart(i) + 1; j < _next_step->rowEnd(i) - 1;
         ++j)
      _calculateNode(i, j);
}

void MFDSolver2d::_printTimeLayerForPlot(std::fstream& file) {
  for (int i = 0; i < _next_step->rows(); ++i) {
    for (size_t j = _next_step->rowStart(i); j < _next_step->rowEnd(i); ++j)
      file << _next_step->at(i, j).x << ' ' << _next_step->at(i, j).y << ' '
           << _next_step->at(i, j).temperature << std::endl;

    file << std::endl;
  }
  file << std::endl;
}

void MFDSolver2d::printTimeLayer(std::fstream& file) {
  file.precision(4);
  for (int i = _next_step->_rows; i >= 0; --i) {
    for (int j = _next_step->rowStart(i); j < _next_step->rowEnd(i); ++j)
      file << std::setw(5) << std::setfill(' ') << ' '
           << _next_step->at(i, j).temperature;
    file << std::endl;
  }
}