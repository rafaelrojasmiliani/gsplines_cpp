link_directories(${CMAKE_SOURCE_DIR}/build)
link_libraries(${GSL_LIBRARIES} ${Eigen3_LIBRARIES} gsplines)
include_directories(test_basislegendre PUBLIC ${CMAKE_SOURCE_DIR}/include
                    ${GSL_INCLUDE_DIRS} ${Eigen3_LIBRARIES})
add_executable(test_basislegendre ${CMAKE_SOURCE_DIR}/tests/basis_legendre.cpp)
add_executable(interpolator ${CMAKE_SOURCE_DIR}/tests/interpolator.cpp)
add_test(test_basislegendre test_basislegendre)
add_test(interpolator interpolator)
