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


int main(int argc, char** argv){
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
    FittingLineFct<Z2i::Point>::findBoundingBox(tL,bmin,bmax);
    //Draw the input points
    Delaunay_triangulation_2::size_type n = dt2.number_of_vertices();
    std::vector<Vertex_handle> TV(n);
    size_t i = 0;
    vector<Z2i::Point> vecVertices;
    for (Vertex_iterator it = dt2.vertices_begin(), end = dt2.vertices_end(); it != end; ++it) {
        TV[i++] = it;
        vecVertices.push_back(Z2i::Point((*it).point()[0],(*it).point()[1]));
        aBoard.fillCircle((*it).point()[0],(*it).point()[1],0.3);
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
        double w=FittingLineFct<Z2i::Point>::widthTriangle<Z2i::Point,Z3i::Point>(vecVertices.at((*it)[0]),vecVertices.at((*it)[1]),vecVertices.at((*it)[2]));
        if(w<width)
            vecTriangle.push_back(*it);
    }
    
    //For each triangle admissible compute the number of fitting
    int max=0, idMax=-1, index;
    for(size_t idT=0; idT<vecTriangle.size(); idT++) {
        Z2i::Point p1=vecVertices.at(vecTriangle.at(idT)[0]);
        Z2i::Point p2=vecVertices.at(vecTriangle.at(idT)[1]);
        Z2i::Point p3=vecVertices.at(vecTriangle.at(idT)[2]);
        vector<Z2i::Point> res=FittingLineFct<Z2i::Point>::fittingTriangle<Z2i::Point,Z3i::Point>(p1,p2,p3,vecVertices,width,index);
        if(res.size()>=max) {
            max=res.size();
            idMax=idT;
        }
    }
    Z2i::Point p1=vecVertices.at(vecTriangle.at(idMax)[0]);
    Z2i::Point p2=vecVertices.at(vecTriangle.at(idMax)[1]);
    Z2i::Point p3=vecVertices.at(vecTriangle.at(idMax)[2]);
    vector<Z2i::Point> res=FittingLineFct<Z2i::Point>::fittingTriangle<Z2i::Point,Z3i::Point>(p1,p2,p3,vecVertices,width,index);
    
    //Draw the best fitting line
    Z2i::RealPoint dp1,dp2,dp3,dp4;
    FittingLineFct<Z2i::Point>::drawFittingTriangle<Z2i::Point,Z3i::Point>(p1,p2,p3,bmin,bmax,index,dp1,dp2,dp3,dp4);
    aBoard.setLineWidth(3.0);
    aBoard.setPenColor(DGtal::Color::Red);
    aBoard.drawLine(dp1[0],dp1[1],dp2[0],dp2[1]);
    aBoard.drawLine(dp3[0],dp3[1],dp4[0],dp4[1]);
    aBoard.setPenColor(DGtal::Color::Blue);
    for (vector<Z2i::Point>::iterator it = res.begin(), end = res.end(); it != end; ++it)
        aBoard.fillCircle((*it)[0],(*it)[1],0.4);
    output=filename+"_Output.svg";
    aBoard.saveSVG(output.c_str());
    
    return 0;
}
