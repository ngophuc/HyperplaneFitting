#include "Functions.h"

#include<CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<CGAL/point_generators_3.h>
#include<CGAL/algorithm.h>
#include<CGAL/random_selection.h>
#include<CGAL/Delaunay_triangulation_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Delaunay_triangulation_3<Kernel> Delaunay;
typedef Delaunay::size_type        size_type;
typedef Delaunay::Vertex_handle    Vertex_handle;
typedef Delaunay::Vertex_iterator  Vertex_iterator;
typedef Delaunay::Cell_handle      Cell_handle;
typedef Delaunay::Cell_iterator    Cell_iterator;
typedef Delaunay::Edge_iterator    Edge_iterator;
typedef Delaunay::Facet_iterator   Facet_iterator;
typedef Delaunay::Locate_type      Locate_type;
typedef Delaunay::Point            Point;

template <typename TPoint>
void findBoundingBox(const vector<TPoint>& vec, TPoint& min, TPoint& max)
{
    int minX=vec.at(0)[0], minY=vec.at(0)[1], minZ=vec.at(0)[2];
    int maxX=vec.at(0)[0], maxY=vec.at(0)[1], maxZ=vec.at(0)[2];
    //for(vector<TPoint>::const_iterator it=vec.begin(); it!=vec.end(); it++)
    for(size_t it=1; it<vec.size(); it++) {
        if(vec.at(it)[0]<minX)
            minX=vec.at(it)[0];
        if(vec.at(it)[1]<minY)
            minY=vec.at(it)[1];
        if(vec.at(it)[2]<minZ)
            minZ=vec.at(it)[2];
        
        if(vec.at(it)[0]>maxX)
            maxX=vec.at(it)[0];
        if(vec.at(it)[1]>maxY)
            maxY=vec.at(it)[1];
        if(vec.at(it)[2]>maxZ)
            maxZ=vec.at(it)[2];
    }
    min[0]=minX;
    min[1]=minY;
    min[2]=minZ;
    max[0]=maxX;
    max[1]=maxY;
    max[2]=maxZ;
}

//Equation of line passing through three points: ax+by+cz+d=0
template <typename TPoint3D, typename TPoint4D>
TPoint4D planeEquation(TPoint3D p1, TPoint3D p2, TPoint3D p3)
{
    TPoint3D p1p2 (p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2]);//AB
    TPoint3D p1p3 (p3[0]-p1[0],p3[1]-p1[1],p3[2]-p1[2]);//AC
    int a=p1p2[1]*p1p3[2]-p1p2[2]*p1p3[1];
    int b=p1p2[2]*p1p3[0]-p1p2[0]*p1p3[2];
    int c=p1p2[0]*p1p3[1]-p1p2[1]*p1p3[0];
    int d=-(a*p1[0]+b*p1[1]+c*p1[2]);
    return TPoint4D(a,b,c,d);
    //return TPoint4D(p1p2[1]*p1p3[2]-p1p2[2]*p1p3[1],p1p2[2]*p1p3[0]-p1p2[0]*p1p3[2],p1p2[0]*p1p3[1]-p1p2[1]*p1p3[0],-(a*p1[0]+b*p1[1]+c*p1[2]));
}
//Equation of line of direction (a,b,c) passing through a point p: d=-ax-by-cz
template <typename TPoint3D, typename TPoint4D>
TPoint4D planeEquation(int a, int b, int c, TPoint3D p)
{
    return TPoint4D(a,b,c,-(a*p[0]+b*p[1]+c*p[2]));
}

template <typename TPoint>
Z3i::RealPoint projectPoint(const double aa, const double bb, const double cc, const double dd, const TPoint aM)
{
    double t = -(aa*aM[0]+bb*aM[1]+cc*aM[2]+dd)/(aa*aa+bb*bb+cc*cc);
    return Z3i::RealPoint( aa*t+aM[0], bb*t+aM[1], cc*t+aM[2] );
}

//Plane equation: a x + by + cz +d=0
template <typename TPoint>
double distancePointPlane(TPoint p, double a, double b, double c, double d)
{
    if(a==0 && b==0 && c==0)
        return 0;
    return fabs(a*p[0] + b*p[1] + c*p[2])/sqrt(a*a+b*b+c*c);
}

double distancePoints(Z3i::RealPoint p1, Z3i::RealPoint p2)
{
    return sqrt((p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]) + (p1[2]-p2[2])*(p1[2]-p2[2]));
}

template <typename TPoint3D, typename TPoint4D>
double widthTetradron(TPoint3D p1, TPoint3D p2, TPoint3D p3, TPoint3D p4)
{
    TPoint4D pl1=planeEquation<TPoint3D,TPoint4D>(p1,p2,p3);
    TPoint4D pl2=planeEquation<TPoint3D,TPoint4D>(p1,p2,p4);
    TPoint4D pl3=planeEquation<TPoint3D,TPoint4D>(p2,p3,p4);
    TPoint4D pl4=planeEquation<TPoint3D,TPoint4D>(p1,p3,p4);
    
    Z3i::RealPoint pp1=projectPoint((double)pl3[0],(double)pl3[1],(double)pl3[2],(double)pl3[3],p1);
    Z3i::RealPoint pp2=projectPoint((double)pl4[0],(double)pl4[1],(double)pl4[2],(double)pl4[3],p2);
    Z3i::RealPoint pp3=projectPoint((double)pl2[0],(double)pl2[1],(double)pl2[2],(double)pl2[3],p3);
    Z3i::RealPoint pp4=projectPoint((double)pl1[0],(double)pl1[1],(double)pl1[2],(double)pl1[3],p4);
    
    double d1=distancePoints(p1,pp1);//distancePointPlane(pp1, (double)l1[0],(double)l1[1],(double)l1[2]);
    double d2=distancePoints(p2,pp2);//distancePointPlane(pp2, (double)l2[0],(double)l2[1],(double)l2[2]);
    double d3=distancePoints(p3,pp3);//distancePointPlane(pp3, (double)l3[0],(double)l3[1],(double)l3[2]);
    double d4=distancePoints(p4,pp4);//distancePointPlane(pp4, (double)l4[0],(double)l4[1],(double)l4[2]);
    
    return std::min(d1,std::min(d2,std::min(d3,d4)));
}

template <typename TPoint3D, typename TPoint4D>
double widthTetradron(TPoint3D p1, TPoint3D p2, TPoint3D p3, TPoint3D p4, int index)
{
    TPoint4D pl1=planeEquation<TPoint3D,TPoint4D>(p1,p2,p3);
    TPoint4D pl2=planeEquation<TPoint3D,TPoint4D>(p1,p2,p4);
    TPoint4D pl3=planeEquation<TPoint3D,TPoint4D>(p2,p3,p4);
    TPoint4D pl4=planeEquation<TPoint3D,TPoint4D>(p1,p3,p4);
    
    Z3i::RealPoint pp1=projectPoint((double)pl3[0],(double)pl3[1],(double)pl3[2],(double)pl3[3],p1);
    Z3i::RealPoint pp2=projectPoint((double)pl4[0],(double)pl4[1],(double)pl4[2],(double)pl4[3],p2);
    Z3i::RealPoint pp3=projectPoint((double)pl2[0],(double)pl2[1],(double)pl2[2],(double)pl2[3],p3);
    Z3i::RealPoint pp4=projectPoint((double)pl1[0],(double)pl1[1],(double)pl1[2],(double)pl1[3],p4);
    
    if(index==0)
        return distancePoints(p1,pp1);
    if(index==1)
        return distancePoints(p2,pp2);
    if(index==2)
        return distancePoints(p3,pp3);
    return distancePoints(p4,pp4);
}

//Digital plane : mu <= ax + by + cz <= omega
template <typename TPoint>
bool belongDP(int a, int b, int c, int mu, int omega, TPoint p)
{
    long int r1=a*p[0]+b*p[1]+c*p[2]+mu;
    long int r2=a*p[0]+b*p[1]+c*p[2]+omega;
    return (r1*r2)<=0;
}

template <typename TPoint>
vector<TPoint> countFitting(int a, int b, int c, int mu, int omega, const vector<TPoint>& vecP)
{
    vector<TPoint> p;
    for(size_t it=0; it<vecP.size(); it++)
        if(belongDP(a,b,c,mu,omega,vecP.at(it)))
            p.push_back(vecP.at(it));
    return p;
}

//Given a tetradron (t1t2t3t4), compute the fitting plane
template <typename TPoint3D, typename TPoint4D>
vector<TPoint3D> fittingTetradron(TPoint3D t1, TPoint3D t2, TPoint3D t3, TPoint3D t4, const vector<TPoint3D>& vecP, double width, int& index)
{
    vector<TPoint3D> vP1,vP2,vP3,vP4;
    double w1=widthTetradron<TPoint3D,TPoint4D>(t1,t2,t3,t4,0);
    if(w1<=width) {//plane 1 : t2t3t4->t1
        TPoint4D pl1=planeEquation<TPoint3D,TPoint4D>(t2, t3, t4);
        TPoint4D pl2=planeEquation<TPoint3D,TPoint4D>(pl1[0], pl1[1], pl1[2], t1);
        int a=pl1[0];
        int b=pl1[1];
        int c=pl1[2];
        int mu=pl1[3],omega=pl2[3];
        vP1=countFitting(a,b,c,mu,omega,vecP);
    }
    double w2=widthTetradron<TPoint3D,TPoint4D>(t1,t2,t3,t4,1);
    if(w2<=width) {//plane 2 : t1t3t4->t2
        TPoint4D pl1=planeEquation<TPoint3D,TPoint4D>(t1, t3, t4);
        TPoint4D pl2=planeEquation<TPoint3D,TPoint4D>(pl1[0], pl1[1], pl1[2], t2);
        int a=pl1[0];
        int b=pl1[1];
        int c=pl1[2];
        int mu=pl1[3],omega=pl2[3];
        vP2=countFitting(a,b,c,mu,omega,vecP);
    }
    double w3=widthTetradron<TPoint3D,TPoint4D>(t1,t2,t3,t4,2);
    if(w3<=width) {//plane 3 : t1t2t4->t3
        TPoint4D pl1=planeEquation<TPoint3D,TPoint4D>(t1, t2, t4);
        TPoint4D pl2=planeEquation<TPoint3D,TPoint4D>(pl1[0], pl1[1], pl1[2], t3);
        int a=pl1[0];
        int b=pl1[1];
        int c=pl1[2];
        int mu=pl1[3],omega=pl2[3];
        vP3=countFitting(a,b,c,mu,omega,vecP);
    }
    double w4=widthTetradron<TPoint3D,TPoint4D>(t1,t2,t3,t4,3);
    if(w4<=width) {//plane 4 : t1t2t3->t4
        TPoint4D pl1=planeEquation<TPoint3D,TPoint4D>(t1, t2, t3);
        TPoint4D pl2=planeEquation<TPoint3D,TPoint4D>(pl1[0], pl1[1], pl1[2], t4);
        int a=pl1[0];
        int b=pl1[1];
        int c=pl1[2];
        int mu=pl1[3],omega=pl2[3];
        vP4=countFitting(a,b,c,mu,omega,vecP);
    }
    
    //Sort index of vP1...vP4 to find max fitting
    vector<size_t> s;
    s.push_back(vP1.size());
    s.push_back(vP2.size());
    s.push_back(vP3.size());
    s.push_back(vP4.size());
    vector<size_t> ss=sort_indexes(s);
    index=ss.at(0);
    if(index==0)
        return vP1;
    if(index==1)
        return vP2;
    if(index==2)
        return vP3;
    return vP4;
}

template <typename TPoint3D, typename TPoint4D>
void drawFittingTetradron(TPoint3D t1, TPoint3D t2, TPoint3D t3, TPoint3D t4, int index,
                          int minX, int minY, int minZ, int maxX, int maxY, int maxZ,
                          Z3i::RealPoint& p11, Z3i::RealPoint& p12, Z3i::RealPoint& p13, Z3i::RealPoint& p14,
                          Z3i::RealPoint& p21, Z3i::RealPoint& p22, Z3i::RealPoint& p23, Z3i::RealPoint& p24)
{
    TPoint4D pl1,pl2;
    if(index==0) {//plane 1 : t2t3t4->t1
        pl1=planeEquation<TPoint3D,TPoint4D>(t2, t3, t4);
        pl2=planeEquation<TPoint3D,TPoint4D>(pl1[0], pl1[1], pl1[2], t1);
    }
    else if(index==1) {//plane 2 : t1t3t4->t2
        pl1=planeEquation<TPoint3D,TPoint4D>(t1, t3, t4);
        pl2=planeEquation<TPoint3D,TPoint4D>(pl1[0], pl1[1], pl1[2], t2);
    }
    else if(index==2) {//plane 3 : t1t2t4->t3
        pl1=planeEquation<TPoint3D,TPoint4D>(t1, t2, t4);
        pl2=planeEquation<TPoint3D,TPoint4D>(pl1[0], pl1[1], pl1[2], t3);
    }
    else {//plane 4 : t1t2t3->t4
        pl1=planeEquation<TPoint3D,TPoint4D>(t1, t2, t3);
        pl2=planeEquation<TPoint3D,TPoint4D>(pl1[0], pl1[1], pl1[2], t4);
    }
    
    int a=pl1[0];
    int b=pl1[1];
    int c=pl1[2];
    int mu=pl1[3],omega=pl2[3];
    if(c!=0) {
        //ax+by+cz+mu=0
        p11[0]=minX;
        p11[1]=minY;
        p11[2]=-(mu+a*p11[0]+b*p11[1])/c;
        p12[0]=minX;
        p12[1]=maxY;
        p12[2]=-(mu+a*p12[0]+b*p12[1])/c;
        p13[0]=maxX;
        p13[1]=maxY;
        p13[2]=-(mu+a*p13[0]+b*p13[1])/c;
        p14[0]=maxX;
        p14[1]=minY;
        p14[2]=-(mu+a*p14[0]+b*p14[1])/c;
        //ax+by+cz+omega=0
        p21[0]=minX;
        p21[1]=minY;
        p21[2]=-(omega+a*p21[0]+b*p21[1])/c;
        p22[0]=minX;
        p22[1]=maxY;
        p22[2]=-(omega+a*p22[0]+b*p22[1])/c;
        p23[0]=maxX;
        p23[1]=maxY;
        p23[2]=-(omega+a*p23[0]+b*p23[1])/c;
        p24[0]=maxX;
        p24[1]=minY;
        p24[2]=-(omega+a*p24[0]+b*p24[1])/c;
    }
    else if(b!=0) {
        //ax+by+cz+mu=0
        p11[0]=minX;
        p11[2]=minZ;
        p11[1]=-(mu+a*p11[0]+c*p11[2])/b;
        p12[0]=minX;
        p12[2]=maxZ;
        p12[1]=-(mu+a*p12[0]+c*p12[2])/b;
        p13[0]=maxX;
        p13[2]=maxZ;
        p13[1]=-(mu+a*p13[0]+c*p13[2])/b;
        p14[0]=maxX;
        p14[2]=minZ;
        p14[1]=-(mu+a*p14[0]+c*p14[2])/b;
        //ax+by+cz+omega=0
        p21[0]=minX;
        p21[2]=minZ;
        p21[1]=-(omega+a*p21[0]+c*p21[2])/b;
        p22[0]=minX;
        p22[2]=maxZ;
        p22[1]=-(omega+a*p22[0]+c*p22[2])/b;
        p23[0]=maxX;
        p23[2]=maxZ;
        p23[1]=-(omega+a*p23[0]+c*p23[2])/b;
        p24[0]=maxX;
        p24[2]=minZ;
        p24[1]=-(omega+a*p24[0]+c*p24[2])/b;
    }
    else if(a!=0) {
        //ax+by+cz+mu=0
        p11[2]=minZ;
        p11[1]=minY;
        p11[0]=-(mu+c*p11[2]+b*p11[1])/a;
        p12[2]=minZ;
        p12[1]=maxY;
        p12[0]=-(mu+c*p12[2]+b*p12[1])/a;
        p13[2]=maxZ;
        p13[1]=maxY;
        p13[0]=-(mu+c*p13[2]+b*p13[1])/a;
        p14[2]=maxZ;
        p14[1]=minY;
        p14[0]=-(mu+c*p14[2]+b*p14[1])/a;
        //ax+by+cz+omega=0
        p21[2]=minZ;
        p21[1]=minY;
        p21[0]=-(omega+c*p21[2]+b*p21[1])/a;
        p22[2]=minZ;
        p22[1]=maxY;
        p22[0]=-(omega+c*p22[2]+b*p22[1])/a;
        p23[2]=maxZ;
        p23[1]=maxY;
        p23[0]=-(omega+c*p23[2]+b*p23[1])/a;
        p24[2]=maxZ;
        p24[1]=minY;
        p24[0]=-(omega+c*p24[2]+b*p24[1])/a;
    }
}

int main(int argc, char** argv)
{
    /*****  Read input points *****/
    string filename="test1.txt";//cube2_d8_boundary6
    ifstream inFile;
    inFile.open(filename.c_str());
    int nbPoint=0, x=0, y=0, z=0;
    inFile >> nbPoint;
    std::vector<Z3i::Point> tL;
    for(int it=0; it<nbPoint; it++){
        inFile >> x >> y >> z;
        tL.push_back(Z3i::Point(x,y,z));
    }
    std::list<Point> L;
    for(vector<Z3i::Point>::const_iterator it=tL.begin(); it!=tL.end(); it++)
        L.push_front(Point((*it)[0],(*it)[1],(*it)[2]));
    /*****  Read input points *****/
    
    /***** Delaunay Triangulation *****/
    Delaunay T(L.begin(), L.end());
    std::size_t n = T.number_of_vertices(); //Nb of vertices
    std::size_t m = T.number_of_cells(); //Nb of facets
    if (n == 0) return -1;
    /***** Triangulation *****/
    
    /***** Write the triangulation *****/
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
    Cell_iterator it;
    i = 0;
    CGAL::Unique_hash_map<Cell_handle, std::size_t > C;
    for(it = T.cells_begin(); it != T.cells_end(); ++it)
        C[it] = i++;
    CGAL_triangulation_assertion( i == m );
    /***** Write the triangulation *****/
    
    /***** copy points cgal -> dgtal *****/
    vector<Z3i::Point> vecPoints;
    for (Vertex_iterator it = T.vertices_begin(), end = T.vertices_end(); it != end; ++it)
        if(it != T.vertices_begin())//do not consider infinite vertex
            vecPoints.push_back(Z3i::Point((*it).point()[0],(*it).point()[1],(*it).point()[2]));
    vector<Z3i::Point> vecVertices=vecPoints;
    vector<Z4i::Point> vecTetra;//cells
    for(Cell_iterator it = T.cells_begin(); it != T.cells_end(); ++it)
        if(V[it->vertex(0)]!=0 && V[it->vertex(1)]!=0 && V[it->vertex(2)]!=0 && V[it->vertex(3)]!=0) //do not consider infinite vertex
            vecTetra.push_back(Z4i::Point(V[it->vertex(0)]-1,V[it->vertex(1)]-1,V[it->vertex(2)]-1,V[it->vertex(3)]-1));
    vector<Z3i::Point> vecTri;//facets
    for(size_t it = 0; it != vecTetra.size(); it++) {
        vecTri.push_back(Z3i::Point(vecTetra.at(it)[0],vecTetra.at(it)[1],vecTetra.at(it)[2]));//012
        vecTri.push_back(Z3i::Point(vecTetra.at(it)[0],vecTetra.at(it)[2],vecTetra.at(it)[3]));//023
        vecTri.push_back(Z3i::Point(vecTetra.at(it)[1],vecTetra.at(it)[2],vecTetra.at(it)[3]));//123
        vecTri.push_back(Z3i::Point(vecTetra.at(it)[0],vecTetra.at(it)[1],vecTetra.at(it)[3]));//013
    }
    //cout<<"vecTetra="<<vecTetra.size()<<endl;
    /***** copy points cgal -> dgtal *****/
    
    /***** Filter over tetradron width *****/
    double width=0.5, w=0;
    vector<Z4i::Point> vecTetraFilter;
    vector<Z3i::Point> vecTriFilter;//facets
    for(vector<Z4i::Point>::const_iterator it=vecTetra.begin(); it!=vecTetra.end(); it++)
    {
        Z4i::Point tetra = *it;//vecTetra.at(it)
        int t1=tetra[0], t2=tetra[1], t3=tetra[2], t4=tetra[3];
        Z3i::Point p1=vecVertices.at(t1), p2=vecVertices.at(t2), p3=vecVertices.at(t3), p4=vecVertices.at(t4);
        w = widthTetradron<Z3i::Point,Z4i::Point>(p1, p2, p3, p4);
        if(w<=width)
        {
            vecTetraFilter.push_back(tetra);
            vecTriFilter.push_back(Z3i::Point(t1,t2,t3));//123
            vecTriFilter.push_back(Z3i::Point(t1,t3,t4));//134
            vecTriFilter.push_back(Z3i::Point(t2,t3,t4));//234
            vecTriFilter.push_back(Z3i::Point(t1,t2,t4));//124
        }
    }
    cout<<"vecTetraFilter.size="<<vecTetraFilter.size()<<endl;
    /***** Filter over triangle width *****/
    
    /***** Fitting plane from tetradra *****/
    //For each triangle admissible compute the number of fitting
    int max=0, idMax=-1, index;
    for(size_t idT=0; idT<vecTetraFilter.size(); idT++)
    {
        Z3i::Point p1=vecVertices.at(vecTetraFilter.at(idT)[0]);
        Z3i::Point p2=vecVertices.at(vecTetraFilter.at(idT)[1]);
        Z3i::Point p3=vecVertices.at(vecTetraFilter.at(idT)[2]);
        Z3i::Point p4=vecVertices.at(vecTetraFilter.at(idT)[3]);
        
        vector<Z3i::Point> res=fittingTetradron<Z3i::Point,Z4i::Point>(p1,p2,p3,p4,vecVertices,width,index);
        //cout<<"fitting of tetradron "<<idT<<" (index="<<index<<") contains "<<res.size()<<" points"<<endl;
        if(res.size()>=max) {
            max=res.size();
            idMax=idT;
        }
    }
    //cout<<"Max fitting is "<<idMax<<" with "<<max<<"/"<<vecVertices.size()<<" points"<<endl;
    /***** Fitting plane from tetradra *****/
    
    /***** Display the result of fitting *****/
    QApplication application(argc,argv);
    MyViewer viewer;
    viewer.show();
    viewer.setBackgroundColor(QColor(255, 255, 255, 0));
    viewer << SetMode3D("PointVector", "Both");//Grid
    viewer << CustomColors3D(Color(0,250,0,255),Color(0,250,0,255));
    //Display points
    viewer << SetMode3D("PointVector", "Both");//Grid
    viewer << CustomColors3D(Color(0,250,0,255),Color(0,250,0,255));
    for (size_t i=0; i<vecPoints.size(); i++) {
        Z3i::Point p(vecPoints[i][0],vecPoints[i][1],vecPoints[i][2]);
        viewer.addBall(p,0.1);//viewer << p;
    }
    
    //display the best fitting plane
    viewer << CustomColors3D(Color(250,0,0,100),Color(250,0,0,100));
    Z3i::Point p1=vecVertices.at(vecTetraFilter.at(idMax)[0]);
    Z3i::Point p2=vecVertices.at(vecTetraFilter.at(idMax)[1]);
    Z3i::Point p3=vecVertices.at(vecTetraFilter.at(idMax)[2]);
    Z3i::Point p4=vecVertices.at(vecTetraFilter.at(idMax)[3]);
    cout<<"Best fitting :"<<p1<<p2<<p3<<p4;
    vector<Z3i::Point> res=fittingTetradron<Z3i::Point,Z4i::Point>(p1,p2,p3,p4,vecVertices,width,index);
    
    Z3i::Point bmin,bmax;
    findBoundingBox(vecVertices,bmin,bmax);
    Z3i::RealPoint p11,p12,p13,p14,p21,p22,p23,p24;
    drawFittingTetradron<Z3i::Point,Z4i::Point>(p1,p2,p3,p4,index,bmin[0],bmin[1],bmin[2],bmax[0],bmax[1],bmax[2],p11,p12,p13,p14,p21,p22,p23,p24);
    viewer.addTriangle(p1,p2,p3);
    viewer.addTriangle(p1,p2,p4);
    viewer.addTriangle(p4,p2,p3);
    viewer.addTriangle(p1,p4,p3);
    viewer.addQuad(p11, p12, p13, p14);
    viewer.addQuad(p21, p22, p23, p24);
    /***** Display the result of fitting *****/
    
    viewer << CustomColors3D(DGtal::Color::Gray, DGtal::Color::Gray);
    viewer << MyViewer::updateDisplay;
    return application.exec();
}
