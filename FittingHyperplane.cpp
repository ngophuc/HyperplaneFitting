#include "Functions.h"

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
        double w = FittingHyperplaneFct<Z4i::Point>::widthFullCell<Z4i::Point>(p1, p2, p3, p4, p5);
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
        vector<Z4i::Point> res=FittingHyperplaneFct<Z4i::Point>::fittingFullCell<Z4i::Point>(p1,p2,p3,p4,p5,vecVertices,width,index);
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
    vector<Z4i::Point> res=FittingHyperplaneFct<Z4i::Point>::fittingFullCell<Z4i::Point>(p1,p2,p3,p4,p5,vecVertices,width,index);
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
