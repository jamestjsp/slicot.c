set(UNIT_TEST_SOURCES
    # Add unit test files here
)

if(UNIT_TEST_SOURCES)
    add_executable(unit_tests ${UNIT_TEST_SOURCES})
    target_link_libraries(unit_tests
        slicot
        gtest_main
        gmock
    )
    gtest_discover_tests(unit_tests)
endif()
