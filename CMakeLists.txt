project(gsplines)
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
set(CMAKE_BUILD_TYPE Release)
# project( gplines++ VERSION 0.0.1 DESCRIPTION "Generalized Splines in c++")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(PYTHON_VERSION
    3.6
    CACHE STRING "Build bindings for Python version")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
add_library(
  gsplines SHARED
  ${CMAKE_SOURCE_DIR}/src/BasisLegendre.cpp
  ${CMAKE_SOURCE_DIR}/src/Interpolator.cpp
  ${CMAKE_SOURCE_DIR}/src/PiecewiseFunction.cpp
  ${CMAKE_SOURCE_DIR}/src/Sobolev.cpp)

target_link_libraries(gsplines ${GSL_LIBRARIES} ${Eigen3_LIBRARIES})
target_include_directories(
  gsplines PUBLIC ${CMAKE_SOURCE_DIR}/include ${GSL_INCLUDE_DIRS}
                  ${Eigen3_INCLUDE_DIRS})
# set_target_properties(gsplines PROPERTIES LIBRARY_OUTPUT_DIRECTORY lib)

set(PYBIND11_PYTHON_VERSION "3.6")
link_directories(${CMAKE_SOURCE_DIR}/build)
add_subdirectory(${CMAKE_SOURCE_DIR}/modules/pybind11)
pybind11_add_module(pygsplines ${CMAKE_SOURCE_DIR}/bindings/pybasis.cpp)
target_include_directories(pygsplines PUBLIC include ${GSL_INCLUDE_DIRS}
                                             /usr/include/eigen3)
target_link_libraries(pygsplines PUBLIC ${GSL_LIBRARIES} ${Eigen3_LIBRARIES}
                                        gsplines)
enable_testing()
add_subdirectory(tests)
