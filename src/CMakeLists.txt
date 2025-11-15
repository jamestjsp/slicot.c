set(SLICOT_SOURCES
    MA/ma02ad.c
    MB/mb01qd.c
    MB/mb01uy.c
    MB/mb03oy.c
    SG/sg03bd.c
    SG/sg03br.c
    SG/sg03bu.c
    SG/sg03bv.c
    SG/sg03bw.c
    SG/sg03bx.c
    TB/tb01vy.c
    TB/tb01wd.c
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
