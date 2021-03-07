cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
# project( gplines++ VERSION 0.0.1 DESCRIPTION "Generalized Splines in c++")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(PYTHON_VERSION
    3.6
    CACHE STRING "Build bindings for Python version")
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
add_library(gsplines SHARED src/BasisLegendre.cpp)
target_link_libraries(gsplines ${GSL_LIBRARIES} ${Eigen3_LIBRARIES})
target_include_directories(gsplines PUBLIC include ${GSL_INCLUDE_DIRS}
                                           ${Eigen3_INCLUDE_DIRS})

find_package(pybind11 REQUIRED)

pybind11_add_module(pygsplines bindings/pybasis.cpp)
message(WARNING ${Eigen3_INCLUDE_DIRS})
target_include_directories(pygsplines PUBLIC include ${GSL_INCLUDE_DIRS}
                                             /usr/include/eigen3)
target_link_libraries(pygsplines ${GSL_LIBRARIES} ${Eigen3_LIBRARIES})
# enable_testing() add_subdirectory(tests)
