set(CTEST_PROJECT_NAME "gsplines")
# find_program(MEMORYCHECK_COMMAND dsfsdf)
set(MEMORYCHECK_COMMAND_OPTIONS
    "--trace-children=yes --leak-check=full --error-exitcode=1 --show-leak-kinds=all --track-origins=yes --dsymutil=yes"
)
set(CTEST_MEMORYCHECK_SUPPRESSIONS_FILE
    "${PROJECT_SOURCE_DIR}/valgrind-python.supp")
