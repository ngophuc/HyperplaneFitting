# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( fittingHyperplane )

cmake_minimum_required(VERSION 2.8.11)
aux_source_directory(. SRC_LIST)

#Required in DGtal
FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

#Required in Egen3
FIND_PACKAGE(Eigen3 REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
LINK_DIRECTORIES(${EIGEN3_LIBRARY_DIRS})

#Required CGAL and its components
find_package( CGAL QUIET COMPONENTS  )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()

endif()

# include for local directory

# include for local package

add_definitions(-g)
# Creating entries for target: convexHull3D
# ############################

add_executable(FittingLine FittingLine Functions.h Functions.ih)
target_link_libraries(FittingLine ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${EIGEN3_INCLUDE_DIRS} ${OpenCV_LIBS})

add_executable(FittingPlane FittingPlane Functions.h Functions.ih)
target_link_libraries(FittingPlane ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${EIGEN3_INCLUDE_DIRS} ${OpenCV_LIBS})

add_executable(FittingHyperplane FittingHyperplane Functions.h Functions.ih)
target_link_libraries(FittingHyperplane ${DGTAL_LIBRARIES} ${DGtalToolsLibDependencies} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} ${EIGEN3_INCLUDE_DIRS} ${OpenCV_LIBS})
