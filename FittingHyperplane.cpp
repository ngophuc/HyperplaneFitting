#include <CGAL/config.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/algorithm.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>

const int D = 4;    // Dimension
typedef CGAL::Epick_d< CGAL::Dimension_tag<D> > K;
typedef CGAL::Delaunay_triangulation<K> Delaunay;
typedef Delaunay::Point Point;
typedef Delaunay::size_type size_type;
typedef Delaunay::Face Face;
typedef std::vector<Face> Faces;
typedef Delaunay::Full_cell_iterator Full_cell_iterator;
typedef Delaunay::Facet Facet;
typedef Delaunay::Full_cell_handle Full_cell_handle;
typedef std::vector<Full_cell_handle> Full_cells;
typedef Delaunay::Vertex_handle     Vertex_handle;
typedef Delaunay::Vertex_iterator   Vertex_iterator;
typedef Delaunay::Facet_iterator    Facet_iterator;
typedef Delaunay::Locate_type       Locate_type;

#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>
#include <time.h>

#include "DGtal/base/BasicTypes.h"
#include "DGtal/helpers/StdDefs.h"

using namespace std;
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

template <typename TPoint>
void findBoundingBox(const vector<TPoint>& vec, TPoint& min, TPoint& max)
{
    auto minX=vec.at(0)[0], minY=vec.at(0)[1], minZ=vec.at(0)[2], minT=vec.at(0)[3];
    auto maxX=vec.at(0)[0], maxY=vec.at(0)[1], maxZ=vec.at(0)[2], maxT=vec.at(0)[3];
    for(size_t it=1; it<vec.size(); it++) {
        if(vec.at(it)[0]<minX)
            minX=vec.at(it)[0];
        if(vec.at(it)[1]<minY)
            minY=vec.at(it)[1];
        if(vec.at(it)[2]<minZ)
            minZ=vec.at(it)[2];
        if(vec.at(it)[3]<minT)
            minT=vec.at(it)[3];
        
        if(vec.at(it)[0]>maxX)
            maxX=vec.at(it)[0];
        if(vec.at(it)[1]>maxY)
            maxY=vec.at(it)[1];
        if(vec.at(it)[2]>maxZ)
            maxZ=vec.at(it)[2];
        if(vec.at(it)[3]>maxT)
            maxT=vec.at(it)[3];
    }
    min[0]=minX;
    min[1]=minY;
    min[2]=minZ;
    min[3]=minT;
    max[0]=maxX;
    max[1]=maxY;
    max[2]=maxZ;
    max[3]=maxT;
}

template <typename TType>
TType findDet3x3(
        TType a11, TType a12, TType a13,
        TType a21, TType a22, TType a23,
        TType a31, TType a32, TType a33 )
{
    return( a11*a22*a33 + a12*a23*a31 + a13*a21*a32 -
            a13*a22*a31 - a12*a21*a33 - a11*a23*a32 );
}

template <typename TType>
TType findDet4x4(
        TType a11, TType a12, TType a13, TType a14,
        TType a21, TType a22, TType a23, TType a24,
        TType a31, TType a32, TType a33, TType a34,
        TType a41, TType a42, TType a43, TType a44 )
{
    return( a11*findDet3x3(a22, a23, a24, a32, a33, a34, a42, a43, a44) -
            a12*findDet3x3(a21, a23, a24, a31, a33, a34, a41, a43, a44) +
            a13*findDet3x3(a21, a22, a24, a31, a32, a34, a41, a42, a44) -
            a14*findDet3x3(a21, a22, a23, a31, a32, a33, a41, a42, a43));
}

vector<double> solve4x4LinearEquation(
        double a11, double a12, double a13, double a14,
        double a21, double a22, double a23, double a24,
        double a31, double a32, double a33, double a34,
        double a41, double a42, double a43, double a44,
        double k1,  double k2,  double k3,  double k4 )
{
    vector<double> answer; //double[] answer = new double[4];
    double det = findDet4x4( a11, a12, a13, a14, a21, a22, a23, a24,
                             a31, a32, a33, a34, a41, a42, a43, a44 );
    if( det == 0.0 )
    {
        return answer;    //No unique solution
    }

    double a = findDet4x4(  k1, a12, a13, a14, k2,  a22, a23, a24,
                             k3, a32, a33, a34, k4,  a42, a43, a44 )/det;
    double b = findDet4x4( a11,  k1, a13, a14, a21,  k2, a23, a24,
                            a31,  k3, a33, a34, a41,  k4, a43, a44 )/det;
    double c = findDet4x4( a11, a12,  k1, a14, a21, a22,  k2, a24,
                            a31, a32,  k3, a34, a41, a42,  k4, a44 )/det;
    double d = findDet4x4( a11, a12, a13,  k1, a21, a22, a23,  k2,
                            a31, a32, a33,  k3, a41, a42, a43,  k4 )/det;
    answer.push_back(a);
    answer.push_back(b);
    answer.push_back(c);
    answer.push_back(d);
    return( answer );
}

vector<int> solve4x4LinearEquation(
        int a11, int a12, int a13, int a14,
        int a21, int a22, int a23, int a24,
        int a31, int a32, int a33, int a34,
        int a41, int a42, int a43, int a44,
        int k1,  int k2,  int k3,  int k4 )
{
    vector<int> answer; //int[] answer = new int[5];
    int det = findDet4x4( a11, a12, a13, a14, a21, a22, a23, a24,
                             a31, a32, a33, a34, a41, a42, a43, a44 );
    if( det == 0.0 )
    {
        return answer;    //No unique solution
    }
    int e = det;
    int a = findDet4x4(  k1, a12, a13, a14, k2,  a22, a23, a24,
                          k3, a32, a33, a34, k4,  a42, a43, a44 );
    int b = findDet4x4( a11,  k1, a13, a14, a21,  k2, a23, a24,
                            a31,  k3, a33, a34, a41,  k4, a43, a44 );
    int c = findDet4x4( a11, a12,  k1, a14, a21, a22,  k2, a24,
                            a31, a32,  k3, a34, a41, a42,  k4, a44 );
    int d = findDet4x4( a11, a12, a13,  k1, a21, a22, a23,  k2,
                            a31, a32, a33,  k3, a41, a42, a43,  k4 );
    answer.push_back(a);
    answer.push_back(b);
    answer.push_back(c);
    answer.push_back(d);
    answer.push_back(e);
    return( answer );
}

//Equation of line passing through three points: ax+by+cz+dt+e=0
template <typename TPoint4D> //FIXME
vector<long long int> hyperplaneEquation(TPoint4D p1, TPoint4D p2, TPoint4D p3, TPoint4D p4)
{
    auto x1=p1[0], y1=p1[1], z1=p1[2], t1=p1[3];
    auto x2=p2[0], y2=p2[1], z2=p2[2], t2=p2[3];
    auto x3=p3[0], y3=p3[1], z3=p3[2], t3=p3[3];
    auto x4=p4[0], y4=p4[1], z4=p4[2], t4=p4[3];
    vector<int> answer = solve4x4LinearEquation(
                                x1,  y1,  z1,  t1,
                                x2,  y2,  z2,  t2,
                                x3,  y3,  z3,  t3,
                                x4,  y4,  z4,  t4,
                                -1,  -1,  -1,  -1);
    long long int a=0, b=0, c=0, d=0, e=0;
    if(answer.size()!=0) {
        a=answer[0];
        b=answer[1];
        c=answer[2];
        d=answer[3];
        e=answer[4];
    }
    
    vector<long long int> hp;
    hp.push_back(a);
    hp.push_back(b);
    hp.push_back(c);
    hp.push_back(d);
    hp.push_back(e);
    return hp;
}

//Equation of line of direction (a,b,c,d) passing through a point p: e=-ax-by-cz-dt
template <typename TPoint4D>
vector<long long int> hyperplaneEquation(long long int a, long long int b, long long int c, long long int d, TPoint4D p)
{
    vector<long long int> hp;
    hp.push_back(a);
    hp.push_back(b);
    hp.push_back(c);
    hp.push_back(d);
    hp.push_back(-(a*p[0]+b*p[1]+c*p[2]+d*p[3]));
    return hp;
}

template <typename TPoint>
Z4i::RealPoint projectPoint(const double aa, const double bb, const double cc, const double dd, const double ee, const TPoint aM)
{
    double t = -(aa*aM[0]+bb*aM[1]+cc*aM[2]+dd*aM[3]+ee)/(aa*aa+bb*bb+cc*cc+dd*dd);
    return Z4i::RealPoint( aa*t+aM[0], bb*t+aM[1], cc*t+aM[2], dd*t+aM[3] );
}

//Plane equation: a x + by + cz +d=0
template <typename TPoint>
double distancePointHyperplane(TPoint p, double a, double b, double c, double d, double e)
{
    if(a==0 && b==0 && c==0 && d==0)
        return 0;
    return fabs(a*p[0] + b*p[1] + c*p[2] + d*p[3] + e)/sqrt(a*a+b*b+c*c+d*d);
}

double distancePoints(Z4i::RealPoint p1, Z4i::RealPoint p2)
{
    return sqrt((p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]) + (p1[2]-p2[2])*(p1[2]-p2[2]) + (p1[3]-p2[3])*(p1[3]-p2[3]));
}

template <typename TPoint4D>
double widthFullCell(TPoint4D p1, TPoint4D p2, TPoint4D p3, TPoint4D p4, TPoint4D p5) {
    vector<long long int> pl1=hyperplaneEquation<TPoint4D>(p2,p3,p4,p5);//plane 1 : t1t2t3t5->t1
    vector<long long int> pl2=hyperplaneEquation<TPoint4D>(p1,p3,p4,p5);//plane 2 : t1t3t4t5->t2
    vector<long long int> pl3=hyperplaneEquation<TPoint4D>(p1,p2,p4,p5);//plane 3 : t1t2t4t5->t3
    vector<long long int> pl4=hyperplaneEquation<TPoint4D>(p1,p2,p3,p5);//plane 4 : t1t2t3t5->t4
    vector<long long int> pl5=hyperplaneEquation<TPoint4D>(p1,p2,p3,p4);//plane 5 : t1t2t3t4->t5
    
    Z4i::RealPoint pp1=projectPoint((double)pl1[0],(double)pl1[1],(double)pl1[2],(double)pl1[3],(double)pl1[4],p1);
    Z4i::RealPoint pp2=projectPoint((double)pl2[0],(double)pl2[1],(double)pl2[2],(double)pl2[3],(double)pl2[4],p2);
    Z4i::RealPoint pp3=projectPoint((double)pl3[0],(double)pl3[1],(double)pl3[2],(double)pl3[3],(double)pl3[4],p3);
    Z4i::RealPoint pp4=projectPoint((double)pl4[0],(double)pl4[1],(double)pl4[2],(double)pl4[3],(double)pl4[4],p4);
    Z4i::RealPoint pp5=projectPoint((double)pl5[0],(double)pl5[1],(double)pl5[2],(double)pl5[3],(double)pl5[4],p5);
    
    double d1=distancePoints(p1,pp1);
    double d2=distancePoints(p2,pp2);
    double d3=distancePoints(p3,pp3);
    double d4=distancePoints(p4,pp4);
    double d5=distancePoints(p5,pp5);
    
    return std::min(d1,std::min(d2,std::min(d3,std::min(d4,d5))));
}

template <typename TPoint4D>
double widthFullCell(TPoint4D p1, TPoint4D p2, TPoint4D p3, TPoint4D p4, TPoint4D p5, int index) {
    vector<long long int> pl1=hyperplaneEquation<TPoint4D>(p2,p3,p4,p5);//plane 1 : t1t2t3t5->t1
    vector<long long int> pl2=hyperplaneEquation<TPoint4D>(p1,p3,p4,p5);//plane 2 : t1t3t4t5->t2
    vector<long long int> pl3=hyperplaneEquation<TPoint4D>(p1,p2,p4,p5);//plane 3 : t1t2t4t5->t3
    vector<long long int> pl4=hyperplaneEquation<TPoint4D>(p1,p2,p3,p5);//plane 4 : t1t2t3t5->t4
    vector<long long int> pl5=hyperplaneEquation<TPoint4D>(p1,p2,p3,p4);//plane 5 : t1t2t3t4->t5
    
    Z4i::RealPoint pp1=projectPoint((double)pl1[0],(double)pl1[1],(double)pl1[2],(double)pl1[3],(double)pl1[4],p1);
    Z4i::RealPoint pp2=projectPoint((double)pl2[0],(double)pl2[1],(double)pl2[2],(double)pl2[3],(double)pl2[4],p2);
    Z4i::RealPoint pp3=projectPoint((double)pl3[0],(double)pl3[1],(double)pl3[2],(double)pl3[3],(double)pl3[4],p3);
    Z4i::RealPoint pp4=projectPoint((double)pl4[0],(double)pl4[1],(double)pl4[2],(double)pl4[3],(double)pl4[4],p4);
    Z4i::RealPoint pp5=projectPoint((double)pl5[0],(double)pl5[1],(double)pl5[2],(double)pl5[3],(double)pl5[4],p5);
    
    if(index==0)
        return distancePoints(p1,pp1);
    if(index==1)
        return distancePoints(p2,pp2);
    if(index==2)
        return distancePoints(p3,pp3);
    if(index==3)
        return distancePoints(p4,pp4);
    return distancePoints(p5,pp5);
}

//Digital plane : mu <= ax + by + cz + dt <= omega
template <typename TPoint>
bool belongHP(int a, int b, int c, int d, int mu, int omega, TPoint p)
{
    long int r1=a*p[0]+b*p[1]+c*p[2]+d*p[3]+mu;
    long int r2=a*p[0]+b*p[1]+c*p[2]+d*p[3]+omega;
    return (r1*r2)<=0;
}

template <typename TPoint>
vector<TPoint> countFitting(int a, int b, int c, int d, int mu, int omega, const vector<TPoint>& vecP)
{
    vector<TPoint> p;
    for(size_t it=0; it<vecP.size(); it++)
        if(belongHP(a,b,c,d,mu,omega,vecP.at(it)))
            p.push_back(vecP.at(it));
    return p;
}

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v)
{
    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);
    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});
    return idx;
}

//Given a tetradron (t1t2t3t4), compute the fitting plane
template <typename TPoint4D>
vector<TPoint4D> fittingFullCell(TPoint4D t1, TPoint4D t2, TPoint4D t3, TPoint4D t4, TPoint4D t5, const vector<TPoint4D>& vecP, double width, int& index)
{
    vector<TPoint4D> vP1,vP2,vP3,vP4,vP5;
    double w1=widthFullCell<TPoint4D>(t1,t2,t3,t4,t5,0);
    if(w1<=width) {//plane 1 : t2t3t4t5->t1
        vector<long long int> pl1=hyperplaneEquation<TPoint4D>(t2, t3, t4, t5);
        vector<long long int> pl2=hyperplaneEquation<TPoint4D>(pl1[0], pl1[1], pl1[2], pl1[3], t1);
        auto a=pl1[0];
        auto b=pl1[1];
        auto c=pl1[2];
        auto d=pl1[3];
        auto mu=std::min(pl1[4],pl2[4]),omega=std::max(pl1[4],pl2[4]);
        vP1=countFitting(a,b,c,d,mu,omega,vecP);
    }
    double w2=widthFullCell<TPoint4D>(t1,t2,t3,t4,t5,1);
    if(w2<=width) {//plane 2 : t1t3t4t5->t2
        vector<long long int> pl1=hyperplaneEquation<TPoint4D>(t1, t3, t4, t5);
        vector<long long int> pl2=hyperplaneEquation<TPoint4D>(pl1[0], pl1[1], pl1[2], pl1[3], t2);
        auto a=pl1[0];
        auto b=pl1[1];
        auto c=pl1[2];
        auto d=pl1[3];
        auto mu=std::min(pl1[4],pl2[4]),omega=std::max(pl1[4],pl2[4]);
        vP2=countFitting(a,b,c,d,mu,omega,vecP);
    }
    double w3=widthFullCell<TPoint4D>(t1,t2,t3,t4,t5,2);
    if(w3<=width) {//plane 3 : t1t2t4t5->t3
        vector<long long int> pl1=hyperplaneEquation<TPoint4D>(t1, t2, t4, t5);
        vector<long long int> pl2=hyperplaneEquation<TPoint4D>(pl1[0], pl1[1], pl1[2], pl1[3], t3);
        auto a=pl1[0];
        auto b=pl1[1];
        auto c=pl1[2];
        auto d=pl1[3];
        auto mu=std::min(pl1[4],pl2[4]),omega=std::max(pl1[4],pl2[4]);
        vP3=countFitting(a,b,c,d,mu,omega,vecP);
    }
    double w4=widthFullCell<TPoint4D>(t1,t2,t3,t4,t5,3);
    if(w4<=width) {//plane 4 : t1t2t3t5->t4
        vector<long long int> pl1=hyperplaneEquation<TPoint4D>(t1, t2, t3, t5);
        vector<long long int> pl2=hyperplaneEquation<TPoint4D>(pl1[0], pl1[1], pl1[2], pl1[3], t4);
        auto a=pl1[0];
        auto b=pl1[1];
        auto c=pl1[2];
        auto d=pl1[3];
        auto mu=std::min(pl1[4],pl2[4]),omega=std::max(pl1[4],pl2[4]);
        vP4=countFitting(a,b,c,d,mu,omega,vecP);
    }
    double w5=widthFullCell<TPoint4D>(t1,t2,t3,t4,t5,4);
    if(w5<=width) {//plane 5 : t1t2t3t4->t5
        vector<long long int> pl1=hyperplaneEquation<TPoint4D>(t1, t2, t3, t4);
        vector<long long int> pl2=hyperplaneEquation<TPoint4D>(pl1[0], pl1[1], pl1[2], pl1[3], t5);
        auto a=pl1[0];
        auto b=pl1[1];
        auto c=pl1[2];
        auto d=pl1[3];
        auto mu=std::min(pl1[4],pl2[4]),omega=std::max(pl1[4],pl2[4]);
        vP5=countFitting(a,b,c,d,mu,omega,vecP);
    }
    //Sort index of vP1...vP4 to find max fitting
    vector<size_t> s;
    s.push_back(vP1.size());
    s.push_back(vP2.size());
    s.push_back(vP3.size());
    s.push_back(vP4.size());
    s.push_back(vP5.size());
    vector<size_t> ss=sort_indexes(s);
    index=ss.at(0);
    if(index==0) return vP1;
    if(index==1) return vP2;
    if(index==2) return vP3;
    if(index==3) return vP4;
    return vP5;
}

int main(int argc, char** argv)
{
	if(argc!=2) {
        std::cout<<"Usage: FittingHyperplane <filename> \n";
        exit(EXIT_FAILURE);
    }
    string baseInputName=argv[1];
    string singleName = baseInputName.substr(baseInputName.find_last_of("/")+1);
    string outputExt = baseInputName.substr(baseInputName.find_last_of(".")+1);
    singleName = singleName.substr(0, singleName.find_last_of("."));
    /*****  Read input points *****/
    string filename=singleName;
    string input=baseInputName;
    ifstream inFile;
    inFile.open(input.c_str());
    int nbPoint=0, x=0, y=0, z=0, t=0;
    inFile >> nbPoint;
    std::vector<Z4i::Point> tL;
    for(int it=0; it<nbPoint; it++){
        inFile >> x >> y >> z >> t;
        tL.push_back(Z4i::Point(x,y,z,t));
    }
    std::vector<Point> L;
    for(vector<Z4i::Point>::const_iterator it=tL.begin(); it!=tL.end(); it++)
        L.push_back(Point((*it)[0],(*it)[1],(*it)[2],(*it)[3]));
    /*****  Read input points *****/
    
    /***** Delaunay Triangulation *****/
    Delaunay T(D);                      // create triangulation
    CGAL_assertion(T.empty());
    T.insert(L.begin(), L.end());  // compute triangulation
    CGAL_assertion( T.is_valid() );
    std::size_t n = T.number_of_vertices(); //Nb of vertices
    std::size_t m = T.number_of_full_cells() ; //Nb of full cells
    /***** Delaunay Triangulation *****/
    
    /***** Write/print the triangulation *****/
    // Write the vertices
    std::vector<Vertex_handle> TV(n+1);
    size_type i = 0;
    for (Vertex_iterator it = T.vertices_begin(), end = T.vertices_end(); it != end; ++it)
        TV[i++] = it;
    CGAL_triangulation_assertion( i == n+1 );
    CGAL_triangulation_assertion( T.is_infinite(TV[0]) );
    CGAL::Unique_hash_map<Vertex_handle, std::size_t > V;
    V[T.infinite_vertex()] = 0;
    for (i=1; i <= n; i++)
        V[TV[i]] = i;
    // write the cells
    i=0;
    CGAL::Unique_hash_map<Full_cell_handle, std::size_t > C;
    for(Full_cell_iterator cit = T.full_cells_begin(); cit != T.full_cells_end(); ++cit )
        C[cit] = i++;
    CGAL_triangulation_assertion( i == m );
    /***** Write/print the triangulation *****/
    
    /***** copy points cgal -> dgtal *****/
    vector<Z4i::Point> vecPoints;
    for (Vertex_iterator it = T.vertices_begin(), end = T.vertices_end(); it != end; ++it)
        if(it != T.vertices_begin())//do not consider infinite vertex
            vecPoints.push_back(Z4i::Point((*it).point()[0],(*it).point()[1],(*it).point()[2],(*it).point()[3]));
    vector<Z4i::Point> vecVertices=vecPoints;
    vector<vector<int> > vecFullCells;
    for(Full_cell_iterator it = T.full_cells_begin(); it != T.full_cells_end(); ++it)
        if(V[it->vertex(0)]!=0 && V[it->vertex(1)]!=0 && V[it->vertex(2)]!=0 && V[it->vertex(3)]!=0 && V[it->vertex(4)]!=0) { //do not consider infinite vertex
            vector<int> cell;
            cell.push_back(V[it->vertex(0)]-1);
            cell.push_back(V[it->vertex(1)]-1);
            cell.push_back(V[it->vertex(2)]-1);
            cell.push_back(V[it->vertex(3)]-1);
            cell.push_back(V[it->vertex(4)]-1);
            vecFullCells.push_back(cell);
        }
    /***** copy points cgal -> dgtal *****/
    
    /***** Filter over tetradron width *****/
    int omega=1;
    double width=omega;//1
    vector<vector<int> > vecFullCellFilter;
    for(vector<vector<int> >::const_iterator it=vecFullCells.begin(); it!=vecFullCells.end(); it++) {
        vector<int> cell = *it;//vecTetra.at(it)
        int t1=cell[0], t2=cell[1], t3=cell[2], t4=cell[3], t5=cell[4];
        Z4i::Point p1=vecVertices.at(t1), p2=vecVertices.at(t2), p3=vecVertices.at(t3), p4=vecVertices.at(t4), p5=vecVertices.at(t5);
        double w = widthFullCell<Z4i::Point>(p1, p2, p3, p4,p5);
        if(w<=width && w>1e-3) {
            vecFullCellFilter.push_back(cell);
        }
    }
    /***** Filter over triangle width *****/
    
    /***** Fitting plane from tetradra *****/
    //For each triangle admissible compute the number of fitting
    int max=0, idMax=-1, index;
    for(size_t idT=0; idT<vecFullCellFilter.size(); idT++) {
        Z4i::Point p1=vecVertices.at(vecFullCellFilter.at(idT)[0]);
        Z4i::Point p2=vecVertices.at(vecFullCellFilter.at(idT)[1]);
        Z4i::Point p3=vecVertices.at(vecFullCellFilter.at(idT)[2]);
        Z4i::Point p4=vecVertices.at(vecFullCellFilter.at(idT)[3]);
        Z4i::Point p5=vecVertices.at(vecFullCellFilter.at(idT)[4]);
        vector<Z4i::Point> res=fittingFullCell<Z4i::Point>(p1,p2,p3,p4,p5,vecVertices,width,index);
        if(res.size()>=max) {
            max=res.size();
            idMax=idT;
        }
    }
    Z4i::Point p1=vecVertices.at(vecFullCellFilter.at(idMax)[0]);
    Z4i::Point p2=vecVertices.at(vecFullCellFilter.at(idMax)[1]);
    Z4i::Point p3=vecVertices.at(vecFullCellFilter.at(idMax)[2]);
    Z4i::Point p4=vecVertices.at(vecFullCellFilter.at(idMax)[3]);
    Z4i::Point p5=vecVertices.at(vecFullCellFilter.at(idMax)[4]);
    vector<Z4i::Point> res=fittingFullCell<Z4i::Point>(p1,p2,p3,p4,p5,vecVertices,width,index);
    /***** Fitting plane from tetradra *****/

    /****Write the result to file ***/
    ofstream outfile;
    std::string output=filename+"_Inliers.txt";
    outfile.open (output.c_str());
    outfile << res.size() << std::endl;
    for (vector<Z4i::Point>::iterator it = res.begin(), end = res.end(); it != end; ++it)
        outfile << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << " " << (*it)[3] <<std::endl;
    outfile.close();
    /****Write the result to file ***/

    return EXIT_SUCCESS;
}
