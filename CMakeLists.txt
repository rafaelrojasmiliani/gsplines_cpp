cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(gsplines VERSION "0.0.1")
include(CTest)
find_program(CTEST_MEMORYCHECK_COMMAND valgrind)
# -mfpmath avx2
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -funroll-loops -mfpmath=sse ")
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -pthread ")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_INSTALL_PREFIX /usr)
set(PYTHON_VERSION
    3.8
    CACHE STRING "Build bindings for Python version")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# include ifopt with custon ROS installation
find_package(ifopt REQUIRED)
find_package(Eigen3 REQUIRED)
add_library(
  gsplines SHARED
  ${PROJECT_SOURCE_DIR}/src/Basis/BasisLegendre.cpp
  ${PROJECT_SOURCE_DIR}/src/Basis/BasisLagrange.cpp
  ${PROJECT_SOURCE_DIR}/src/Basis/Basis0101.cpp
  ${PROJECT_SOURCE_DIR}/src/Basis/Basis.cpp
  ${PROJECT_SOURCE_DIR}/src/Interpolator.cpp
  ${PROJECT_SOURCE_DIR}/src/GSpline.cpp
  ${PROJECT_SOURCE_DIR}/src/Tools.cpp
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
  ${PROJECT_SOURCE_DIR}/src/Collocation/GaussLobattoPointsWeights.cpp
  ${PROJECT_SOURCE_DIR}/src/Collocation/GaussLobattoLagrange.cpp
  ${PROJECT_SOURCE_DIR}/src/Collocation/GaussLobattoLagrangeFunctionals.cpp
  # --
)

target_link_libraries(gsplines ${EIGEN3_LIBRARIES} ${ifopt_LIBRARIES})
target_include_directories(
  gsplines
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:include> ${EIGEN3_INCLUDE_DIRS}
         ${ifopt_INCLUDE_DIRS})

set_target_properties(
  gsplines
  PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON
             VERSION ${gsplines_VERSION}
             SOVERSION ${gsplines_VERSION_MAJOR}.${gsplines_VERSION_MINOR}
             EXPORT_NAME GSplines)

link_directories(${PROJECT_SOURCE_DIR}/build)

if(NOT TARGET pybind11)
  set(PYBIND11_PYTHON_VERSION "3.8")
  add_subdirectory(${PROJECT_SOURCE_DIR}/modules/pybind11)
endif()

pybind11_add_module(pygsplines ${PROJECT_SOURCE_DIR}/bindings/bindings.cpp)

target_include_directories(
  pygsplines PUBLIC ${EIGEN3_INCLUDE_DIRS}
                    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/modules/pybind11>)

target_link_libraries(pygsplines PUBLIC ${EIGEN3_LIBRARIES} ${ifopt_LIBRARIES}
                                        gsplines)

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

include(GNUInstallDirs)
set(INSTALL_CMAKE_CONFIG_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/GSplines)

install(
  TARGETS gsplines
  EXPORT GSplinesTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c
          "import site; print(site.getsitepackages()[0])"
  OUTPUT_VARIABLE PYTHON_INSTALL_PATH
  OUTPUT_STRIP_TRAILING_WHITESPACE)

install(DIRECTORY python/ DESTINATION ${PYTHON_INSTALL_PATH})
install(TARGETS pygsplines DESTINATION ${PYTHON_INSTALL_PATH}/gsplines)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

export(EXPORT GSplinesTargets
       FILE ${CMAKE_CURRENT_BINARY_DIR}/GSplinesTargets.cmake)

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/GSplinesConfigVersion.cmake
  COMPATIBILITY SameMajorVersion)
configure_package_config_file(
  cmake/GSplinesConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/GSplinesConfig.cmake
  INSTALL_DESTINATION ${INSTALL_CMAKE_CONFIG_DIR})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/GSplinesConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/GSplinesConfigVersion.cmake
        DESTINATION ${INSTALL_CMAKE_CONFIG_DIR})
install(EXPORT GSplinesTargets DESTINATION ${INSTALL_CMAKE_CONFIG_DIR})

# Packaging
set(CPACK_PACKAGE_VENDOR "Rafael A. Rojas")
set(CPACK_GENERATOR "DEB")
set(CPACK_PACKAGE_VERSION ${gsplines_VERSION})
set(CPACK_SYSTEM_NAME ${CMAKE_HOST_SYSTEM_PROCESSOR})

# Debian
find_program(DPKG_PROG dpkg DOC "'dpkg' executable")
if(DPKG_PROG)
  execute_process(
    COMMAND ${DPKG_PROG} --print-architecture
    OUTPUT_VARIABLE CPACK_DEBIAN_PACKAGE_ARCHITECTURE
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  # Change system name to use the correct architecture in file name
  set(CPACK_SYSTEM_NAME ${CPACK_DEBIAN_PACKAGE_ARCHITECTURE})
endif()
# Debian versions require a dash
set(CPACK_DEBIAN_PACKAGE_VERSION ${CPACK_PACKAGE_VERSION}-1)
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Rafael A. Rojas")
set(CPACK_DEBIAN_PACKAGE_DEPENDS
    "coinor-libipopt-dev, ros-noetic-ifopt, libeigen3-dev")

include(CPack)
