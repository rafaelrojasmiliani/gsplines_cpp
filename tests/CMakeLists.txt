add_executable(test_basislegendre ${CMAKE_SOURCE_DIR}/tests/basis_legendre.cpp)
link_directories(${CMAKE_SOURCE_DIR}/build)
target_include_directories(
  test_basislegendre PUBLIC ${CMAKE_SOURCE_DIR}/include ${GSL_INCLUDE_DIRS}
                            ${Eigen3_LIBRARIES})
target_link_libraries(test_basislegendre PUBLIC ${GSL_LIBRARIES}
                                                ${Eigen3_LIBRARIES} gsplines)
add_test(test_basislegendre test_basislegendre)
