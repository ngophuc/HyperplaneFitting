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

//Common functions
vector<size_t> sort_indexes(const vector<size_t> &v);

//Line fitting functions
template<typename TPoint>
struct FittingLineFct
{
public:
    static void findBoundingBox(const vector<TPoint>& vec, TPoint& min, TPoint& max);
    template <typename TPoint2D, typename TPoint3D>
    static TPoint3D lineEquation(TPoint2D p1, TPoint2D p2);
    template <typename TPoint2D, typename TPoint3D>
    static TPoint3D lineEquation(int a, int b, TPoint2D p);
    static Z2i::RealPoint projectPoint(const double aa, const double bb, const double cc, const TPoint aM);
    static double distancePointSegment(TPoint p, double a, double b, double c);
    static double distancePoints(Z2i::RealPoint p1, Z2i::RealPoint p2);
    template <typename TPoint2D, typename TPoint3D>
    static double widthTriangle(TPoint2D p1, TPoint2D p2, TPoint2D p3);
    template <typename TPoint2D, typename TPoint3D>
    static double widthTriangle(TPoint2D p1, TPoint2D p2, TPoint2D p3, int index);
    static bool belongDL(int a, int b, int mu, int omega, TPoint p);
    static vector<TPoint> countFitting(int a, int b, int mu, int omega, const vector<TPoint>& vecP);
    template <typename TPoint2D, typename TPoint3D>
    static vector<TPoint2D> fittingTriangle(TPoint2D t1, TPoint2D t2, TPoint2D t3, const vector<TPoint2D>& vecP, double width, int& index);
    template <typename TPoint2D, typename TPoint3D>
    static void drawFittingTriangle(TPoint2D t1, TPoint2D t2, TPoint2D t3, TPoint2D bmin, TPoint2D bmax, int index, Z2i::RealPoint& p1, Z2i::RealPoint& p2, Z2i::RealPoint& p3, Z2i::RealPoint& p4);
};

//Plane fitting functions
template<typename TPoint>
struct FittingPlaneFct
{
public:
    static void findBoundingBox(const vector<TPoint>& vec, TPoint& min, TPoint& max);
    template <typename TPoint3D, typename TPoint4D>
    static TPoint4D planeEquation(TPoint3D p1, TPoint3D p2, TPoint3D p3);
    template <typename TPoint3D, typename TPoint4D>
    static TPoint4D planeEquation(int a, int b, int c, TPoint3D p);
    static Z3i::RealPoint projectPoint(const double aa, const double bb, const double cc, const double dd, const TPoint aM);
    static double distancePointPlane(TPoint p, double a, double b, double c, double d);
    template <typename TPoint1, typename TPoint2>
    static double distancePoints(TPoint1 p1, TPoint2 p2);
    template <typename TPoint3D, typename TPoint4D>
    static double widthTetradron(TPoint3D p1, TPoint3D p2, TPoint3D p3, TPoint3D p4);
    template <typename TPoint3D, typename TPoint4D>
    static double widthTetradron(TPoint3D p1, TPoint3D p2, TPoint3D p3, TPoint3D p4, int index);
    static bool belongDP(int a, int b, int c, int mu, int omega, TPoint p);
    static vector<TPoint> countFitting(int a, int b, int c, int mu, int omega, const vector<TPoint>& vecP);
    template <typename TPoint3D, typename TPoint4D>
    static vector<TPoint3D> fittingTetradron(TPoint3D t1, TPoint3D t2, TPoint3D t3, TPoint3D t4, const vector<TPoint3D>& vecP, double width, int& index);
    template <typename TPoint3D, typename TPoint4D>
    static void drawFittingTetradron(TPoint3D t1, TPoint3D t2, TPoint3D t3, TPoint3D t4, int index,
                              int minX, int minY, int minZ, int maxX, int maxY, int maxZ,
                              Z3i::RealPoint& p11, Z3i::RealPoint& p12, Z3i::RealPoint& p13, Z3i::RealPoint& p14,
                              Z3i::RealPoint& p21, Z3i::RealPoint& p22, Z3i::RealPoint& p23, Z3i::RealPoint& p24);
};

template<typename TPoint>
struct FittingHyperplaneFct
{
public:
    static void findBoundingBox(const vector<TPoint>& vec, TPoint& min, TPoint& max);
    template <typename TPoint4D>
    static vector<long long int> hyperplaneEquation(TPoint4D p1, TPoint4D p2, TPoint4D p3, TPoint4D p4);
    template <typename TPoint4D>
    static vector<long long int> hyperplaneEquation(long long int a, long long int b, long long int c, long long int d, TPoint4D p);
    static Z4i::RealPoint projectPoint(const double aa, const double bb, const double cc, const double dd, const double ee, const TPoint aM);
    static double distancePointHyperplane(TPoint p, double a, double b, double c, double d, double e);
    template <typename TPoint4D>
    static double widthFullCell(TPoint4D p1, TPoint4D p2, TPoint4D p3, TPoint4D p4, TPoint4D p5);
    template <typename TPoint4D>
    static double widthFullCell(TPoint4D p1, TPoint4D p2, TPoint4D p3, TPoint4D p4, TPoint4D p5, int index);
    static bool belongHP(int a, int b, int c, int d, int mu, int omega, TPoint p);
    static vector<TPoint> countFitting(int a, int b, int c, int d, int mu, int omega, const vector<TPoint>& vecP);
    template <typename TPoint4D>
    static vector<TPoint4D> fittingFullCell(TPoint4D t1, TPoint4D t2, TPoint4D t3, TPoint4D t4, TPoint4D t5, const vector<TPoint4D>& vecP, double width, int& index);
};

#include "Functions.ih" // Includes inline functions.
#endif // FUNCTIONS_H
