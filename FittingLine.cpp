#include "Functions.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef K::Segment_2 Segment_2;
typedef K::Ray_2 Ray_2;
typedef K::Line_2 Line_2;
typedef CGAL::Delaunay_triangulation_2<K>  Delaunay_triangulation_2;
typedef Delaunay_triangulation_2::Vertex_handle    Vertex_handle;
typedef Delaunay_triangulation_2::Face_handle      Face_handle;
typedef Delaunay_triangulation_2::Vertex_iterator  Vertex_iterator;
typedef Delaunay_triangulation_2::Edge_iterator    Edge_iterator;
typedef Delaunay_triangulation_2::Face_iterator    Face_iterator;

//A class to recover Voronoi diagram from stream
struct Cropped_voronoi_from_delaunay{
    std::list<Segment_2> m_cropped_vd;
    Iso_rectangle_2 m_bbox;
    
    Cropped_voronoi_from_delaunay(const Iso_rectangle_2& bbox):m_bbox(bbox){}
    
    template <class RSL>
    void crop_and_extract_segment(const RSL& rsl){
        CGAL::Object obj = CGAL::intersection(rsl,m_bbox);
        const Segment_2* s=CGAL::object_cast<Segment_2>(&obj);
        if (s) m_cropped_vd.push_back(*s);
    }
    
    void operator<<(const Ray_2& ray)    { crop_and_extract_segment(ray); }
    void operator<<(const Line_2& line)  { crop_and_extract_segment(line); }
    void operator<<(const Segment_2& seg){ crop_and_extract_segment(seg); }
};

template <typename TPoint>
void findBoundingBox(const vector<TPoint>& vec, TPoint& min, TPoint& max)
{
    int minX=vec.at(0)[0], minY=vec.at(0)[1];
    int maxX=vec.at(0)[0], maxY=vec.at(0)[1];
    //for(vector<TPoint>::const_iterator it=vec.begin(); it!=vec.end(); it++)
    for(size_t it=1; it<vec.size(); it++) {
        if(vec.at(it)[0]<minX)
            minX=vec.at(it)[0];
        if(vec.at(it)[1]<minY)
            minY=vec.at(it)[1];
        
        if(vec.at(it)[0]>maxX)
            maxX=vec.at(it)[0];
        if(vec.at(it)[1]>maxY)
            maxY=vec.at(it)[1];
    }
    min[0]=minX;
    min[1]=minY;
    max[0]=maxX;
    max[1]=maxY;
}

//Equation of line passing through two points: ax+by+c=0
template <typename TPoint2D, typename TPoint3D>
TPoint3D lineEquation(TPoint2D p1, TPoint2D p2)
{
    return TPoint3D(p1[1]-p2[1],p2[0]-p1[0],p1[0]*p2[1]-p2[0]*p1[1]);
}
//Equation of line of direction (a,b) passing through a point p: c=-ax-by
template <typename TPoint2D, typename TPoint3D>
TPoint3D lineEquation(int a, int b, TPoint2D p)
{
    return TPoint3D(a,b,-(a*p[0]+b*p[1]));
}

template <typename TPoint>
Z2i::RealPoint projectPoint(const double aa, const double bb, const double cc, const TPoint aM)
{
    double xm=aM[0];
    double ym=aM[1];
    //computation line: ax+by+c=0
    double d2=-( aa * aa + bb * bb );
    double s=-bb * xm + aa * ym;
    double xp=( aa * cc + bb * s ) / d2;
    double yp=( bb * cc - aa * s ) / d2;
    return Z2i::RealPoint( xp, yp );
}

//Segment equation: a x + by + c=0
template <typename TPoint>
double distancePointSegment(TPoint p, double a, double b, double c) 
{
    if(a==0 && b==0)
        return 0;
    return fabs(a*p[0] + b*p[1] + c)/sqrt(a*a+b*b);
}

double distancePoints(Z2i::RealPoint p1, Z2i::RealPoint p2)
{
    return sqrt((p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]));
}

template <typename TPoint2D, typename TPoint3D>
double widthTriangle(TPoint2D p1, TPoint2D p2, TPoint2D p3) {
    TPoint3D l1=lineEquation<TPoint2D,TPoint3D>(p2,p3);
    TPoint3D l2=lineEquation<TPoint2D,TPoint3D>(p3,p1);
    TPoint3D l3=lineEquation<TPoint2D,TPoint3D>(p2,p1);
    
    Z2i::RealPoint pp3=projectPoint((double)l1[0],(double)l1[1],(double)l1[2],p1);
    Z2i::RealPoint pp2=projectPoint((double)l2[0],(double)l2[1],(double)l2[2],p2);
    Z2i::RealPoint pp1=projectPoint((double)l3[0],(double)l3[1],(double)l3[2],p3);
    
    double d1=distancePoints(p1,pp1);//distancePointSegment(pp1, (double)l1[0],(double)l1[1],(double)l1[2]);
    double d2=distancePoints(p2,pp2);//distancePointSegment(pp2, (double)l2[0],(double)l2[1],(double)l2[2]);
    double d3=distancePoints(p3,pp3);//distancePointSegment(pp3, (double)l3[0],(double)l3[1],(double)l3[2]);
    
    return std::min(d1,std::min(d2,d3));
}

template <typename TPoint2D, typename TPoint3D>
double widthTriangle(TPoint2D p1, TPoint2D p2, TPoint2D p3, int index) {
    TPoint3D l1=lineEquation<TPoint2D,TPoint3D>(p2,p3);
    TPoint3D l2=lineEquation<TPoint2D,TPoint3D>(p3,p1);
    TPoint3D l3=lineEquation<TPoint2D,TPoint3D>(p2,p1);

    Z2i::RealPoint pp3=projectPoint((double)l1[0],(double)l1[1],(double)l1[2],p1);
    Z2i::RealPoint pp2=projectPoint((double)l2[0],(double)l2[1],(double)l2[2],p2);
    Z2i::RealPoint pp1=projectPoint((double)l3[0],(double)l3[1],(double)l3[2],p3);

    double d1=distancePoints(p1,pp1);//distancePointSegment(pp1, (double)l1[0],(double)l1[1],(double)l1[2]);
    double d2=distancePoints(p2,pp2);//distancePointSegment(pp2, (double)l2[0],(double)l2[1],(double)l2[2]);
    double d3=distancePoints(p3,pp3);//distancePointSegment(pp3, (double)l3[0],(double)l3[1],(double)l3[2]);
    
    if(index==0) return d1;
    if(index==1) return d2;
    return d3;
}

//Digital Straight Line : mu <= ax + by <= omega
template <typename TPoint>
bool belongDL(int a, int b, int mu, int omega, TPoint p)
{
    int r1=a*p[0]+b*p[1]+mu;
    int r2=a*p[0]+b*p[1]+omega;
    return r1*r2<=0;
}
template <typename TPoint>
vector<TPoint> countFitting(int a, int b, int mu, int omega, const vector<TPoint>& vecP)
{
    vector<TPoint> p;
    for(size_t it=0; it<vecP.size(); it++)
        if(belongDL(a,b,mu,omega,vecP.at(it)))
            p.push_back(vecP.at(it));
    return p;
}

//Given a triangle (t1t2t3), compute the fitting line
template <typename TPoint2D, typename TPoint3D>
vector<TPoint2D> fittingTriangle(TPoint2D t1, TPoint2D t2, TPoint2D t3, const vector<TPoint2D>& vecP, double width, int& index)
{
    vector<TPoint2D> vP1,vP2,vP3;
    double w1=widthTriangle<TPoint2D,TPoint3D>(t1,t2,t3,0);
    if(w1<=width) { //Line 1 : t2t3->t1
        TPoint3D l1=lineEquation<TPoint2D,TPoint3D>(t3, t2);
        TPoint3D l2=lineEquation<TPoint2D,TPoint3D>(l1[0], l1[1], t1);
        int a=l1[0];
        int b=l1[1];
        int mu=l1[2],omega=l2[2];
        vP1=countFitting(a,b,mu,omega,vecP);
    }
    double w2=widthTriangle<TPoint2D,TPoint3D>(t1,t2,t3,1);
    if(w2<=width) {//Line 2 : t1t3->t2
        TPoint3D l1=lineEquation<TPoint2D,TPoint3D>(t1, t3);
        TPoint3D l2=lineEquation<TPoint2D,TPoint3D>(l1[0], l1[1], t2);
        int a=l1[0];
        int b=l1[1];
        int mu=l1[2],omega=l2[2];
        vP2=countFitting(a,b,mu,omega,vecP);
    }
    double w3=widthTriangle<TPoint2D,TPoint3D>(t1,t2,t3,2);
    if(w3<=width) {//Line 3 : t1t2->t3
        TPoint3D l1=lineEquation<TPoint2D,TPoint3D>(t2, t1);
        TPoint3D l2=lineEquation<TPoint2D,TPoint3D>(l1[0], l1[1], t3);
        int a=l1[0];
        int b=l1[1];
        int mu=l1[2],omega=l2[2];
        vP3=countFitting(a,b,mu,omega,vecP);
    }
    //Sort index of vP1...vP4 to find max fitting
    vector<size_t> s;
    s.push_back(vP1.size());
    s.push_back(vP2.size());
    s.push_back(vP3.size());
    vector<size_t> ss=sort_indexes(s);
    index=ss.at(0);
    if(index==0)
        return vP1;
    if(index==1)
        return vP2;
    return vP3;
}

template <typename TPoint2D, typename TPoint3D>
void drawFittingTriangle(TPoint2D t1, TPoint2D t2, TPoint2D t3, TPoint2D bmin, TPoint2D bmax, int index, Z2i::RealPoint& p1, Z2i::RealPoint& p2, Z2i::RealPoint& p3, Z2i::RealPoint& p4)
{
    int minX=bmin[0];
    int minY=bmin[1];
    int maxX=bmax[0];
    int maxY=bmax[1];
    
    TPoint3D l1,l2;
    if(index==0)//Line 1 : t2t3->t1
    {
        l1=lineEquation<TPoint2D,TPoint3D>(t3, t2);
        l2=lineEquation<TPoint2D,TPoint3D>(l1[0], l1[1], t1);
    }
    else if(index==1)//Line 2 : t1t3->t2
    {
        l1=lineEquation<TPoint2D,TPoint3D>(t1, t3);
        l2=lineEquation<TPoint2D,TPoint3D>(l1[0], l1[1], t2);
    }
    else//Line 3 : t1t2->t3
    {
        l1=lineEquation<TPoint2D,TPoint3D>(t2, t1);
        l2=lineEquation<TPoint2D,TPoint3D>(l1[0], l1[1], t3);
    }
    
    int a=l1[0];
    int b=l1[1];
    int mu=l1[2],omega=l2[2];
    if(a!=0) {
        //ax+by+mu=0
        p1[1]=minY;
        p1[0]=-(mu+b*p1[1])/a;
        p2[1]=maxY;
        p2[0]=-(mu+b*p2[1])/a;
        //ax+by+omega=0
        p3[1]=minY;
        p3[0]=-(omega+b*p3[1])/a;
        p4[1]=maxY;
        p4[0]=-(omega+b*p4[1])/a;
    }
    else if(b!=0) {
        //ax+by+mu=0
        p1[0]=minX;
        p1[1]=-(mu+a*p1[0])/b;
        p2[0]=maxX;
        p2[1]=-(mu+a*p2[0])/b;
        //ax+by+omega=0
        p3[0]=minX;
        p3[1]=-(omega+a*p3[0])/b;
        p4[0]=maxX;
        p4[1]=-(omega+a*p4[0])/b;
    }
}

int main(){
    /*****  Read input points *****/
    string filename="test1";//cube2_d8_boundary6
    string input=filename+".txt";
    ifstream inFile;
    inFile.open(input.c_str());
    int nbPoint=0, x=0, y=0;
    inFile >> nbPoint;
    std::vector<Z2i::Point> tL;
    for(int it=0; it<nbPoint; it++){
        inFile >> x >> y;
        tL.push_back(Z2i::Point(x,y));
    }
    std::list<Point_2> L;
    for(vector<Z2i::Point>::const_iterator it=tL.begin(); it!=tL.end(); it++)
        L.push_front(Point_2((*it)[0],(*it)[1]));
    /*****  Read input points *****/

    /***** Delaunay Triangulation *****/
    Delaunay_triangulation_2 dt2;
    dt2.insert(L.begin(),L.end());
    /***** Delaunay Triangulation *****/
    
    Board2D aBoard;
    Z2i::Point bmin,bmax;
    findBoundingBox(tL,bmin,bmax);
    bmin=bmin-Z2i::Point(2,2);
    bmax=bmax+Z2i::Point(2,2);
    //Draw the input points
    Delaunay_triangulation_2::size_type n = dt2.number_of_vertices();
    std::vector<Vertex_handle> TV(n);
    size_t i = 0;
    vector<Z2i::Point> vecVertices;
    for (Vertex_iterator it = dt2.vertices_begin(), end = dt2.vertices_end(); it != end; ++it) {
        TV[i++] = it;
        vecVertices.push_back(Z2i::Point((*it).point()[0],(*it).point()[1]));
        aBoard.drawCircle((*it).point()[0],(*it).point()[1],0.05);
    }
    string output=filename+"_Input.svg";
    aBoard.saveSVG(output.c_str());

    //Get triangles from Delaunay triangulation
    CGAL::Unique_hash_map<Vertex_handle, std::size_t > V;
    V[dt2.infinite_vertex()] = 0;
    for (i=0; i < n; i++)
        V[TV[i]] = i;
    vector<Z3i::Point> vecFaces;
    for(Face_iterator it = dt2.faces_begin(); it != dt2.faces_end(); ++it) 
        vecFaces.push_back(Z3i::Point(V[it->vertex(0)],V[it->vertex(1)],V[it->vertex(2)]));
    
    //Filter over triangle width
    double width=2.5;
    vector<Z3i::Point> vecTriangle;
    aBoard.setLineWidth(1.0);
    aBoard.setPenColor(DGtal::Color::Blue);
    for(vector<Z3i::Point>::const_iterator it=vecFaces.begin(); it!= vecFaces.end(); it++) {
        double w=widthTriangle<Z2i::Point,Z3i::Point>(vecVertices.at((*it)[0]),vecVertices.at((*it)[1]),vecVertices.at((*it)[2]));
        if(w<width)
            vecTriangle.push_back(*it);
    }
    
    //For each triangle admissible compute the number of fitting
    int max=0, idMax=-1, index;
    for(size_t idT=0; idT<vecTriangle.size(); idT++) {
        Z2i::Point p1=vecVertices.at(vecTriangle.at(idT)[0]);
        Z2i::Point p2=vecVertices.at(vecTriangle.at(idT)[1]);
        Z2i::Point p3=vecVertices.at(vecTriangle.at(idT)[2]);
        vector<Z2i::Point> res=fittingTriangle<Z2i::Point,Z3i::Point>(p1,p2,p3,vecVertices,width,index);
        //cout<<"fitting of traingle "<<idT<<" contains "<<res.size()<<" points"<<endl;
        if(res.size()>=max) {
            max=res.size();
            idMax=idT;
        }
    }
    //cout<<"Max fitting is "<<idMax<<" with "<<max<<" points"<<endl;
    Z2i::Point p1=vecVertices.at(vecTriangle.at(idMax)[0]);
    Z2i::Point p2=vecVertices.at(vecTriangle.at(idMax)[1]);
    Z2i::Point p3=vecVertices.at(vecTriangle.at(idMax)[2]);
    vector<Z2i::Point> res=fittingTriangle<Z2i::Point,Z3i::Point>(p1,p2,p3,vecVertices,width,index);

    //Draw the best fitting line
    Z2i::RealPoint dp1,dp2,dp3,dp4;
    drawFittingTriangle<Z2i::Point,Z3i::Point>(p1,p2,p3,bmin,bmax,index,dp1,dp2,dp3,dp4);
    aBoard.setLineWidth(3.0);
    aBoard.setPenColor(DGtal::Color::Red);
    aBoard.drawLine(dp1[0],dp1[1],dp2[0],dp2[1]);
    aBoard.drawLine(dp3[0],dp3[1],dp4[0],dp4[1]);
    output=filename+"_Output.svg";
    aBoard.saveSVG(output.c_str());
    
    return 0;
}
