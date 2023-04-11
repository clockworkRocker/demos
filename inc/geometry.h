#ifndef LAB1_INC_GEOMETRY_H_
#define LAB1_INC_GEOMETRY_H_

#include <vector>

const double EPS = 1e-5;

struct Point2d {
  double x, y;

  Point2d(double _x = 0., double _y = 0.);
};

class Line2d {
 public:
  virtual ~Line2d() = default;

  virtual double length() = 0;

  virtual Point2d start() = 0;
  virtual Point2d end() = 0;

  virtual double maxX() = 0;
  virtual double maxY() = 0;
  virtual double minX() = 0;
  virtual double minY() = 0;

  virtual bool hasPoint(const Point2d &point, double precision = EPS) = 0;
  virtual bool intersectsLine(Line2d &line) = 0;
};

class StraightLine2d : public Line2d {
 protected:
  Point2d _start, _end;

 public:
  StraightLine2d(double x1, double y1, double x2, double y2);
  StraightLine2d(const Point2d &point1, const Point2d &point2);

  Point2d start() override;
  Point2d end() override;

  double length() override;

  double maxX() override;
  double maxY() override;
  double minX() override;
  double minY() override;

  bool hasPoint(const Point2d &point, double precision) override;
  bool intersectsLine(Line2d &line) override;
};

class Arc2d : public Line2d {
 protected:
  Point2d _center;
  double _radius;
  double _ang_start, _ang_end;
 public:
  Arc2d(double x, double y, double r, double a1, double a2);
  Arc2d(const Point2d &center, double r, double a1, double a2);

  Point2d center() { return _center; }

  double radius() const { return _radius; }

  double angle1() const { return _ang_start; }

  double angle2() const { return _ang_end; }

  Point2d start() override;
  Point2d end() override;

  double length() override;

  double maxX() override;
  double maxY() override;
  double minX() override;
  double minY() override;

  bool hasPoint(const Point2d &point, double precision) override;
  bool intersectsLine(Line2d &line) override;
};

class Shape2d {
 protected:
  std::vector<Line2d *> _edges;

 public:
  Shape2d() = default;
  explicit Shape2d(std::vector<Line2d *> &edges);
  virtual ~Shape2d() = default;

  Line2d *edge(size_t index) { return index < _edges.size() ? _edges[index] : nullptr; }

  virtual double maxX();
  virtual double maxY();
  virtual double minX();
  virtual double minY();

  bool hasPoint(Point2d &point);
};

#endif //LAB1_INC_GEOMETRY_H_
