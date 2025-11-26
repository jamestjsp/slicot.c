set(SLICOT_SOURCES
    IB/ib01nd.c
    IB/ib01od.c
    IB/ib01oy.c
    MA/ma02ad.c
    MA/ma02ed.c
    SB/sb02mt.c
    SB/sb02nd.c
    MB/mb01qd.c
    MB/mb01rb.c
    MB/mb01rx.c
    MB/mb01td.c
    MB/mb01uy.c
    MB/mb02yd.c
    MB/mb03od.c
    MB/mb03oy.c
    MB/mb04id.c
    MB/mb04iy.c
    MB/mb03ud.c
    MB/mb04kd.c
    MB/mb04od.c
    MB/mb04ow.c
    MB/mb04oy.c
    MD/md03ba.c
    MD/md03bb.c
    MD/md03bd.c
    MD/md03bx.c
    MD/md03by.c
    NF/nf01ad.c
    NF/nf01ay.c
    NF/nf01bd.c
    NF/nf01be.c
    NF/nf01bf.c
    NF/nf01bp.c
    NF/nf01bq.c
    NF/nf01br.c
    NF/nf01bs.c
    NF/nf01by.c
    SG/sg03bd.c
    SG/sg03br.c
    SG/sg03bu.c
    SG/sg03bv.c
    SG/sg03bw.c
    SG/sg03bx.c
    TB/tb01vy.c
    TB/tb01wd.c
    TF/tf01mx.c
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
