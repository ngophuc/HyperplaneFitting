#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>
#include <numeric>
#include <algorithm>
#include <iterator>

using namespace std;

#include "DGtal/base/BasicTypes.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/PointListReader.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/Color.h"

using namespace DGtal;

namespace Z4i {
    typedef DGtal::int32_t Integer;
    typedef SpaceND< 4, Integer > Space;
    typedef Space Z4;
    typedef Space::Point Point;
    typedef Space::Vector Vector;
    typedef Space::RealPoint RealPoint;
    typedef Space::RealVector RealVector;
}
typedef DGtal::Viewer3D<> MyViewer;


//Common functions
vector<size_t> sort_indexes(const vector<size_t> &v);

//Line fitting functions
template <typename TPoint2D, typename TPoint3D>
TPoint3D lineEquation(TPoint2D p1, TPoint2D p2);
template <typename TPoint2D, typename TPoint3D>
TPoint3D lineEquation(int a, int b, TPoint2D p);
template <typename TPoint>
Z2i::RealPoint projectPoint(const double aa, const double bb, const double cc, const TPoint aM);
template <typename TPoint>
double distancePointSegment(TPoint p, double a, double b, double c);
double distancePoints(Z2i::RealPoint p1, Z2i::RealPoint p2);
template <typename TPoint2D, typename TPoint3D>
double widthTriangle(TPoint2D p1, TPoint2D p2, TPoint2D p3);
template <typename TPoint2D, typename TPoint3D>
double widthTriangle(TPoint2D p1, TPoint2D p2, TPoint2D p3, int index);
template <typename TPoint>
bool belongDL(int a, int b, int mu, int omega, TPoint p);
template <typename TPoint>
vector<TPoint> countFitting(int a, int b, int mu, int omega, const vector<TPoint>& vecP);
template <typename TPoint2D, typename TPoint3D>
vector<TPoint2D> fittingTriangle(TPoint2D t1, TPoint2D t2, TPoint2D t3, const vector<TPoint2D>& vecP, double width, int& index);
template <typename TPoint2D, typename TPoint3D>
void drawFittingTriangle(TPoint2D t1, TPoint2D t2, TPoint2D t3, TPoint2D bmin, TPoint2D bmax, int index, Z2i::RealPoint& p1, Z2i::RealPoint& p2, Z2i::RealPoint& p3, Z2i::RealPoint& p4);

//Plane fitting functions
template <typename TPoint3D, typename TPoint4D>
TPoint4D planeEquation(TPoint3D p1, TPoint3D p2, TPoint3D p3);
template <typename TPoint3D, typename TPoint4D>
TPoint4D planeEquation(int a, int b, int c, TPoint3D p);
template <typename TPoint>
Z3i::RealPoint projectPoint(const double aa, const double bb, const double cc, const double dd, const TPoint aM);
template <typename TPoint>
double distancePointPlane(TPoint p, double a, double b, double c, double d);
double distancePoints(Z3i::RealPoint p1, Z3i::RealPoint p2);
template <typename TPoint3D, typename TPoint4D>
double widthTetradron(TPoint3D p1, TPoint3D p2, TPoint3D p3, TPoint3D p4);
template <typename TPoint3D, typename TPoint4D>
double widthTetradron(TPoint3D p1, TPoint3D p2, TPoint3D p3, TPoint3D p4, int index);
template <typename TPoint>
bool belongDP(int a, int b, int c, int mu, int omega, TPoint p);
template <typename TPoint>
vector<TPoint> countFitting(int a, int b, int c, int mu, int omega, const vector<TPoint>& vecP);
template <typename TPoint3D, typename TPoint4D>
vector<TPoint3D> fittingTetradron(TPoint3D t1, TPoint3D t2, TPoint3D t3, TPoint3D t4, const vector<TPoint3D>& vecP, double width, int& index);
template <typename TPoint3D, typename TPoint4D>
void drawFittingTetradron(TPoint3D t1, TPoint3D t2, TPoint3D t3, TPoint3D t4, int index,
                          int minX, int minY, int minZ, int maxX, int maxY, int maxZ,
                          Z3i::RealPoint& p11, Z3i::RealPoint& p12, Z3i::RealPoint& p13, Z3i::RealPoint& p14,
                          Z3i::RealPoint& p21, Z3i::RealPoint& p22, Z3i::RealPoint& p23, Z3i::RealPoint& p24);
#endif // FUNCTIONS_H