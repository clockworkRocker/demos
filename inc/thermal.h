#ifndef LAB1_INC_THERMAL_H_
#define LAB1_INC_THERMAL_H_

#include <fstream>
#include "geometry.h"

enum geometry_types { GM_NODE = 0, GM_EDGE = 1 };

enum boundary_types {
  TM_HEAT_PROOF = 0,
  TM_TEMPERATURE = 1,
  TM_HEAT_FLUX = 2,
  TM_CONVECTION = 3
};

class MFDSolver2d;
struct ThermoNode2d : Point2d {
  double temperature;
  ThermoNode2d();
  ThermoNode2d(double x, double y, double val = 0.);
};

class MFDMesh2d {
 private:
  size_t _rows, _cols;
  double _dx, _dy;
  ThermoNode2d* _nodes;
  size_t* _row_offsets;
  size_t* _row_ends;

 public:
  MFDMesh2d(Shape2d& geometry, double dx, double dy);
  MFDMesh2d(const MFDMesh2d& rhs);
  // MFDMesh2d(Shape2d &geometry, size_t total_nodes);
  ~MFDMesh2d();

  double dx() const;
  double dy() const;
  ThermoNode2d& operator[](size_t index);
  ThermoNode2d& at(size_t row, size_t col);
  size_t rows() const;
  size_t rowStart(size_t row);
  size_t rowEnd(size_t row);
  bool indexOfNode(double x, double y, size_t* index);

  void debugPrint();

  friend MFDSolver2d;
};

class ThermoShape2d : public Shape2d {
 protected:
  struct BoundaryCondition {
    size_t type;
    double value;
    size_t geometry_type;
    size_t geometry_index;

    BoundaryCondition() = default;
    BoundaryCondition(size_t type,
                      size_t geom_type,
                      size_t geom_index,
                      double value);
  };

  std::vector<BoundaryCondition> _boundary_conditions;
  MFDMesh2d* _mesh;

 public:
  explicit ThermoShape2d(std::vector<Line2d*> edges);
  ~ThermoShape2d() override;

  void generateMesh(double dx, double dy);
  bool lookupNode(double x, double y, size_t* index = nullptr);
  bool addBoundary(size_t geometry_type,
                   size_t geometry_index,
                   size_t type,
                   double value);

  friend MFDSolver2d;
};

class MFDSolver2d {
 private:
  ThermoShape2d* _task = nullptr;
  MFDMesh2d* _next_step = nullptr;
  double _time_step;

  inline void _applyBoundaries();
  inline void _applyBoundaryToNode(ThermoShape2d::BoundaryCondition& condition);
  inline void _applyBoundaryToEdge(ThermoShape2d::BoundaryCondition& condition);
  inline void _calculateNode(size_t i, size_t j);
  inline double _ddx(size_t i, size_t j);
  inline double _ddy(size_t i, size_t j);
  inline double _d2dx2(size_t i, size_t j);
  inline double _d2dy2(size_t i, size_t j);
  void _calculateTimeLayer();
  void _printTimeLayerForPlot(std::fstream& file);

 public:
  bool solve(ThermoShape2d& shape,
             double total_time,
             std::fstream& output_file);
  double step() const;
  void printTimeLayer(std::fstream& file);
  ~MFDSolver2d();
};

#endif  // LAB1_INC_THERMAL_H_
