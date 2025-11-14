# Python tests using pytest
find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
    # Check if pytest is available
    execute_process(
        COMMAND ${Python3_EXECUTABLE} -m pytest --version
        RESULT_VARIABLE PYTEST_result
        OUTPUT_QUIET
        ERROR_QUIET
    )

    if(PYTEST_result EQUAL 0)
        message(STATUS "pytest found: ${Python3_EXECUTABLE} -m pytest")

        # Add pytest tests
        file(GLOB PYTEST_FILES "${CMAKE_CURRENT_SOURCE_DIR}/python/test_*.py")
        foreach(test_file ${PYTEST_FILES})
            get_filename_component(test_name ${test_file} NAME_WE)
            add_test(
                NAME ${test_name}
                COMMAND ${Python3_EXECUTABLE} -m pytest ${test_file} -v
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python
            )
            # Ensure Python bindings are built before tests
            set_tests_properties(${test_name} PROPERTIES
                ENVIRONMENT "PYTHONPATH=${PROJECT_SOURCE_DIR}/python:$ENV{PYTHONPATH}"
            )
        endforeach()
    else()
        message(WARNING "pytest not found. Install with: pip install pytest")
    endif()
else()
    message(WARNING "Python3 not found. Tests disabled.")
endif()

add_subdirectory(integration)
