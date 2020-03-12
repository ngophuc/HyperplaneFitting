## HyperplaneFitting

Code source associated to IWCIA'20 paper: <b> Digital hyperplane fitting </b>

## Dependencies
The program uses some C++ 11 feature, so we recommend the use of gcc 4.2 or later to compile. The program requires these libraries to be installed:

* Boost
* DGtal (To install DGtal see https://github.com/DGtal-team/DGtal/blob/master/README.md)
* Eigen3
* CGAL
* CMake 3.15.2

# Installation 
To install the program see <a href="https://github.com/ngophuc/HyperplaneFitting/blob/master/Install.txt">INSTALL.txt</a> file

# Examples
The program takes as input a file containing a list of points with the first line indicating the total number of point, and following by the coordinates of each point. 
* Example for a file containing 2D points 
```
10
1 10
1 26
3 3
3 42
4 18
4 34
8 6
8 39
9 25
10 11
```
* Example for a file containing 3D points 
```
10
8 1 1
9 1 1
10 2 1
11 3 1
12 4 1
11 12 1
12 12 1
7 1 2
8 2 2
9 2 2
```
You can find more examples in <b>Data/Data2D</b> and <b>Data/Data3D</b>

* Run the program <b>FittingLine</b> 
```
  ./FittingLine -i points2D.txt
```
* Run the program <b>FittingPlane</b> 
```
  ./FittingPlane -i points3D.txt
```
