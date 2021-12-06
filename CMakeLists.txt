project(gsplines)
include(CTest)
find_program(CTEST_MEMORYCHECK_COMMAND dsfsdf)
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE
    "-O3 -DNDEBUG -funroll-loops   -mfpmath=sse  -fopenmp")
# project( gplines++ VERSION 0.0.1 DESCRIPTION "Generalized Splines in c++")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(PYTHON_VERSION
    3.6
    CACHE STRING "Build bindings for Python version")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# include ifopt with custon ROS installation
set(ifopt_DIR /opt/ros/noetic/share/ifopt/cmake)
find_package(ifopt REQUIRED)
find_package(Eigen3 REQUIRED)
add_library(
  gsplines SHARED
  ${PROJECT_SOURCE_DIR}/src/Basis/BasisLegendre.cpp
  ${PROJECT_SOURCE_DIR}/src/Basis/BasisLagrange.cpp
  ${PROJECT_SOURCE_DIR}/src/Basis/Basis.cpp
  ${PROJECT_SOURCE_DIR}/src/Interpolator.cpp
  ${PROJECT_SOURCE_DIR}/src/GSpline.cpp
  ${PROJECT_SOURCE_DIR}/src/FunctionalAnalysis/Sobolev.cpp
  ${PROJECT_SOURCE_DIR}/src/FunctionalAnalysis/Integral.cpp
  ${PROJECT_SOURCE_DIR}/src/Optimization/ipopt_interface.cpp
  ${PROJECT_SOURCE_DIR}/src/Optimization/ipopt_solver.cpp
  ${PROJECT_SOURCE_DIR}/src/Functions/FunctionSum.cpp
  ${PROJECT_SOURCE_DIR}/src/Functions/FunctionMul.cpp
  ${PROJECT_SOURCE_DIR}/src/Functions/FunctionsComp.cpp
  ${PROJECT_SOURCE_DIR}/src/Functions/FunctionsConcat.cpp
  ${PROJECT_SOURCE_DIR}/src/Functions/Functions.cpp
  ${PROJECT_SOURCE_DIR}/src/Functions/FunctionExpression.cpp
  ${PROJECT_SOURCE_DIR}/src/Functions/FunctionGenericOperations.cpp
  ${PROJECT_SOURCE_DIR}/src/Functions/FunctionBase.cpp
  ${PROJECT_SOURCE_DIR}/src/Functions/ElementalFunctions.cpp
  ${PROJECT_SOURCE_DIR}/src/Collocation/GaussLobattoPointsWeights.cpp)

target_link_libraries(gsplines ${Eigen3_LIBRARIES} ${ifopt_LIBRARIES})
target_include_directories(
  gsplines PUBLIC ${PROJECT_SOURCE_DIR}/include ${Eigen3_INCLUDE_DIRS}
                  ${ifopt_INCLUDE_DIRS} /usr/include/eigen3)
# set_target_properties(gsplines PROPERTIES LIBRARY_OUTPUT_DIRECTORY lib)

link_directories(${PROJECT_SOURCE_DIR}/build)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  set(PYBIND11_PYTHON_VERSION "3.6")
  add_subdirectory(${PROJECT_SOURCE_DIR}/modules/pybind11)
endif()

pybind11_add_module(pygsplines ${PROJECT_SOURCE_DIR}/bindings/bindings.cpp)
target_include_directories(pygsplines PUBLIC include /usr/include/eigen3)
target_link_libraries(pygsplines PUBLIC ${Eigen3_LIBRARIES} ${ifopt_LIBRARIES}
                                        gsplines)

set_target_properties(pygsplines PROPERTIES OUTPUT_NAME gsplines)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  enable_testing()
  add_subdirectory(tests)
endif()

if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json")
  execute_process(
    COMMAND
      ${CMAKE_COMMAND} -E copy_if_different
      ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
      ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json)
endif()
