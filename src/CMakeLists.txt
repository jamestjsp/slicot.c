set(SLICOT_SOURCES
    # Add source files here as they are created
)

if(BUILD_SHARED_LIBS)
    add_library(slicot SHARED ${SLICOT_SOURCES})
else()
    add_library(slicot STATIC ${SLICOT_SOURCES})
endif()

target_include_directories(slicot PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

set_target_properties(slicot PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)
