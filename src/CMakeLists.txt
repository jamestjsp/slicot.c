set(SLICOT_SOURCES
    MB/mb01qd.c
    MB/mb03oy.c
    TG/tg01fd.c
)

if(BUILD_SHARED_LIBS)
    add_library(slicot SHARED ${SLICOT_SOURCES})
else()
    add_library(slicot STATIC ${SLICOT_SOURCES})
endif()

target_include_directories(slicot PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Link BLAS and LAPACK (already found in root CMakeLists.txt)
target_link_libraries(slicot PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

set_target_properties(slicot PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)
