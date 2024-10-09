include(GoogleTest)

file(GLOB_RECURSE test_list ${PROJECT_SOURCE_DIR}/tests/*.cpp)
file(GLOB_RECURSE test_py_list ${PROJECT_SOURCE_DIR}/tests/*.py)
list(
    REMOVE_ITEM test_py_list "${PROJECT_SOURCE_DIR}/tests/__init__.py"
)
list(
    REMOVE_ITEM test_py_list "${PROJECT_SOURCE_DIR}/tests/tools.py"
)

foreach(file_path ${test_list})
    get_filename_component(barename ${file_path} NAME)
    string(
        REPLACE ".cpp"
                ""
                new_name
                ${barename}
    )
    add_executable(${new_name} ${file_path})
    target_include_directories(${new_name} PUBLIC ${PROJECT_SOURCE_DIR}/include ${Eigen3_LIBRARIES})
    target_link_libraries(
        ${new_name}
        gsplines
        gtest
        pthread
    )
    if(GSPLINES_USE_RUCKIG)
        target_compile_definitions(pygsplines PUBLIC HAS_RUCKIG)
        target_link_libraries(pygsplines PUBLIC ${ruckig_LIBRARIES})
        # target_include_directories(pygsplines PUBLIC ${ruckig_INCLUDE_DIRS})
    endif()

    if(GSPLINES_USE_KDL)
        target_compile_definitions(pygsplines PUBLIC HAS_KDL)
        target_link_libraries(pygsplines PUBLIC ${orocos_kdl_LIBRARIES})
        # target_include_directories(pygsplines PUBLIC ${orocos_kdl_INCLUDE_DIRS})
    endif()
    gtest_discover_tests(${new_name})
endforeach()

foreach(file_path ${test_py_list})
    get_filename_component(barename ${file_path} NAME)
    string(
        REPLACE ".py"
                "_py"
                new_name
                ${barename}
    )
    string(
        REPLACE ".py"
                ""
                mod_name
                ${barename}
    )
    add_test(NAME ${new_name} COMMAND python3 -m unittest tests.${mod_name} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    set_tests_properties(
        ${new_name} PROPERTIES ENVIRONMENT
                               PYTHONPATH=${PROJECT_SOURCE_DIR}/python:$<TARGET_FILE_DIR:pygsplines>:$ENV{PYTHONPATH}
    )

endforeach()
