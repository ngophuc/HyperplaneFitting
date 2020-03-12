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

int main(int argc, char** argv)
{
    if(argc!=2) {
        std::cout<<"Usage: FittingLines <filename> \n";
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
        w = FittingPlaneFct<Z3i::Point>::widthTetradron<Z3i::Point,Z4i::Point>(p1, p2, p3, p4);
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
        
        vector<Z3i::Point> res=FittingPlaneFct<Z3i::Point>::fittingTetradron<Z3i::Point,Z4i::Point>(p1,p2,p3,p4,vecVertices,width,index);
        if(res.size()>=max) {
            max=res.size();
            idMax=idT;
        }
    }
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
        viewer.addBall(p,0.1);
    }
    
    //display the best fitting plane
    Z3i::Point p1=vecVertices.at(vecTetraFilter.at(idMax)[0]);
    Z3i::Point p2=vecVertices.at(vecTetraFilter.at(idMax)[1]);
    Z3i::Point p3=vecVertices.at(vecTetraFilter.at(idMax)[2]);
    Z3i::Point p4=vecVertices.at(vecTetraFilter.at(idMax)[3]);
    vector<Z3i::Point> res=FittingPlaneFct<Z3i::Point>::fittingTetradron<Z3i::Point,Z4i::Point>(p1,p2,p3,p4,vecVertices,width,index);
    viewer << CustomColors3D(Color(0,0,250),Color(0,0,250));
    for (vector<Z3i::Point>::iterator it = res.begin(), end = res.end(); it != end; ++it)
        viewer.addBall(*it,0.15);
    Z3i::Point bmin,bmax;
    FittingPlaneFct<Z3i::Point>::findBoundingBox(vecVertices,bmin,bmax);
    Z3i::RealPoint p11,p12,p13,p14,p21,p22,p23,p24;
    FittingPlaneFct<Z3i::Point>::drawFittingTetradron<Z3i::Point,Z4i::Point>(p1,p2,p3,p4,index,bmin[0],bmin[1],bmin[2],bmax[0],bmax[1],bmax[2],p11,p12,p13,p14,p21,p22,p23,p24);
    viewer << CustomColors3D(Color(250,0,0,100),Color(250,0,0,100));
    viewer.addQuad(p11, p12, p13, p14);
    viewer.addQuad(p21, p22, p23, p24);
    /***** Display the result of fitting *****/
    
    /****Write the result to file ***/
    ofstream outfile;
    std::string output=filename+"_Inliers.txt";
    outfile.open (output.c_str());
    outfile << res.size() << std::endl;
    for (vector<Z3i::Point>::iterator it = res.begin(), end = res.end(); it != end; ++it)
        outfile << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] <<std::endl;
    outfile.close();
    /****Write the result to file ***/
    
    viewer << CustomColors3D(DGtal::Color::Gray, DGtal::Color::Gray);
    viewer << MyViewer::updateDisplay;
    return application.exec();
}
