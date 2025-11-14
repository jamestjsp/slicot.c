set(INTEGRATION_TEST_SOURCES
    # Add integration test files here
)

if(INTEGRATION_TEST_SOURCES)
    add_executable(integration_tests ${INTEGRATION_TEST_SOURCES})
    target_link_libraries(integration_tests
        slicot
        gtest_main
        gmock
    )
    gtest_discover_tests(integration_tests)
endif()
