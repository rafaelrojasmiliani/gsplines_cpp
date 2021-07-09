link_directories(${CMAKE_SOURCE_DIR}/build)
link_libraries(${GSL_LIBRARIES} ${Eigen3_LIBRARIES} gsplines)
include_directories(PUBLIC ${CMAKE_SOURCE_DIR}/include ${GSL_INCLUDE_DIRS}
                    ${Eigen3_LIBRARIES})

file(GLOB_RECURSE test_list ${CMAKE_SOURCE_DIR}/tests/*.cpp)

foreach(file_path ${test_list})
  get_filename_component(barename ${file_path} NAME)
  string(REPLACE ".cpp" "" new_name ${barename})
  add_executable(${new_name} ${file_path})
  add_test(${new_name} ${new_name})
endforeach()
