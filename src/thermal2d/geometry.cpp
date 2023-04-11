#include <cmath>
#include "geometry.h"

#define MIN(x1, x2) ((x1) < (x2) ? (x1) : (x2))
#define MAX(x1, x2) ((x1) > (x2) ? (x1) : (x2))

/* ===================================================== POINTS ===================================================== */
Point2d::Point2d(double _x, double _y) :
    x(_x),
    y(_y) {}

/* ===================================================== LINES ====================================================== */
/* ------------------------------------------------- Straight line -------------------------------------------------- */
StraightLine2d::StraightLine2d(double x1, double y1, double x2, double y2) :
    _start(Point2d(x1, y1)),
    _end(Point2d(x2, y2)) {}

StraightLine2d::StraightLine2d(const Point2d &point1, const Point2d &point2) :
    _start(point1),
    _end(point2) {}

Point2d StraightLine2d::start() { return _start; }

Point2d StraightLine2d::end() { return _end; }

double StraightLine2d::length() {
    return sqrt(
        (_end.x - _start.x) * (_end.x - _start.x) +
            (_end.y - _start.y) * (_end.y - _start.y)
    );
}

double StraightLine2d::maxX() { return MAX(_start.x, _end.x); }

double StraightLine2d::minX() { return MIN(_start.x, _end.x); }

double StraightLine2d::maxY() { return MAX(_start.y, _end.y); }

double StraightLine2d::minY() { return MIN(_start.y, _end.y); }

bool StraightLine2d::hasPoint(const Point2d &point, double precision) {
    if (point.x > MAX(_start.x, _end.x) || point.x < MIN(_start.x, _end.x))
        return false;
    if (point.y > MAX(_start.y, _end.y) || point.y < MIN(_start.y, _end.y))
        return false;
    return fabs((point.y - _start.y) * (_end.x - _start.x) - (point.x - _start.x) * (_end.y - _start.y)) <= precision;
}

bool StraightLine2d::intersectsLine(Line2d &line) {
    double x1 = _start.x;
    double x2 = _end.x;
    double x3 = line.start().x;
    double x4 = line.end().x;
    double y1 = _start.y;
    double y2 = _end.y;
    double y3 = line.start().y;
    double y4 = line.end().y;

    double d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

    if (!d)
        return false;

    double p1 = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / d;
    double p2 = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / d;

    return (p1 >= 0 && p1 <= 1 && p2 >= 0 && p2 <= 1);
}

/* ------------------------------------------------------ Arc ------------------------------------------------------- */
Arc2d::Arc2d(double x, double y, double r, double a1, double a2) :
    _center(x, y),
    _radius(r),
    _ang_start(a1),
    _ang_end(a2) {}

Arc2d::Arc2d(const Point2d &center, double r, double a1, double a2) :
    _center(center),
    _radius(r),
    _ang_start(a1),
    _ang_end(a2) {}

Point2d Arc2d::start() {
    return {
        _center.x + _radius * cos(_ang_start),
        _center.y + _radius * sin(_ang_start)
    };
}

Point2d Arc2d::end() {
    return {
        _center.x + _radius * cos(_ang_end),
        _center.y + _radius * sin(_ang_end)
    };
}

double Arc2d::length() {
    return _radius * (_ang_end - _ang_start);
}

double Arc2d::maxX() {
    if (fabs(0 - _ang_start) / fabs(_ang_end - _ang_start) < 1)
        return _center.x + _radius;

    if (fabs(0 - _ang_start) < fabs(0 - _ang_end))
        return _center.x + _radius * cos(_ang_start);
    else
        return _center.x + _radius * cos(_ang_end);
}

double Arc2d::minX() {
    if (fabs(M_PI - _ang_start) / fabs(_ang_end - _ang_start) < 1)
        return _center.x - _radius;

    if (fabs(M_PI - _ang_start) < fabs(M_PI - _ang_end))
        return _center.x + _radius * cos(_ang_start);
    else
        return _center.x + _radius * cos(_ang_end);
}

double Arc2d::maxY() {
    if (fabs(M_PI_2 - _ang_start) / fabs(_ang_end - _ang_start) < 1)
        return _center.y + _radius;

    if (fabs(M_PI_2 - _ang_start) < fabs(M_PI_2 - _ang_end))
        return _center.y + _radius * sin(_ang_start);
    else
        return _center.y + _radius * sin(_ang_end);
}

double Arc2d::minY() {
    if (fabs(-M_PI_2 - _ang_start) / fabs(_ang_end - _ang_start) < 1)
        return _center.y - _radius;

    if (fabs(-M_PI_2 - _ang_start) < fabs(-M_PI_2 - _ang_end))
        return _center.y + _radius * sin(_ang_start);
    else
        return _center.y + _radius * sin(_ang_end);
}

bool Arc2d::hasPoint(const Point2d &point, double precision) {
    double p1 = (acos((point.x - _center.x) / _radius) - _ang_start) / (_ang_end - _ang_start);
    double p2 = (asin((point.y - _center.y) / _radius) - _ang_start) / (_ang_end - _ang_start);
    return fabs(p1 - p2) <= EPS;
}

bool Arc2d::intersectsLine(Line2d &line) {
    // FIXME
    return false;
}

/* ==================================================== SHAPES ====================================================== */
Shape2d::Shape2d(std::vector<Line2d *> &edges) : _edges(edges) {}

double Shape2d::maxX() {
    double maxX = _edges[0]->maxX();
    for (auto *edge : _edges)
        if (edge->maxX() > maxX)
            maxX = edge->maxX();

    return maxX;
}

double Shape2d::minX() {
    double minX = _edges[0]->minX();
    for (auto *edge : _edges)
        if (edge->minX() < minX)
            minX = edge->minX();

    return minX;
}

double Shape2d::maxY() {
    double maxY = _edges[0]->maxY();
    for (auto *edge : _edges)
        if (edge->maxY() > maxY)
            maxY = edge->maxY();

    return maxY;
}

double Shape2d::minY() {
    double minY = _edges[0]->minY();
    for (auto *edge : _edges)
        if (edge->minY() < minY)
            minY = edge->minY();

    return minY;
}

bool Shape2d::hasPoint(Point2d &point) {
    int intersections = 0;
    StraightLine2d ray(point, Point2d(10 * point.x + 1e6, 10 * point.y + 1e6));

    for (auto &edge : _edges){
        if (edge->hasPoint(point))
            return true;
        if (edge->intersectsLine(ray))
            intersections++;
        }

    return (intersections % 2);
}