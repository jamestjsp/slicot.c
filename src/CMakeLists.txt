set(SLICOT_SOURCES
    AB/ab04md.c
    AB/ab05md.c
    AB/ab05nd.c
    AB/ab05od.c
    AB/ab05pd.c
    AB/ab05qd.c
    AB/ab07md.c
    AB/ab07nd.c
    AB/ab13md.c
    DG/dg01md.c
    DK/dk01md.c
    MC/mc01td.c
    IB/ib01ad.c
    IB/ib01bd.c
    IB/ib01cd.c
    IB/ib01md.c
    IB/ib01my.c
    IB/ib01nd.c
    IB/ib01od.c
    IB/ib01oy.c
    IB/ib01pd.c
    IB/ib01px.c
    IB/ib01py.c
    IB/ib01qd.c
    IB/ib01rd.c
    IB/ib03bd.c
    MA/ma02ad.c
    MA/ma02ed.c
    MA/ma02fd.c
    SB/sb02mt.c
    SB/sb02nd.c
    SB/sb02rd.c
    SB/sb02ru.c
    SB/sb03od.c
    SB/sb03or.c
    SB/sb03ot.c
    SB/sb03ov.c
    SB/sb03oy.c
    SB/sb04px.c
    SB/sb10jd.c
    MB/mb01qd.c
    MB/mb01rb.c
    MB/mb02pd.c
    MB/mb01rx.c
    MB/mb01sd.c
    MB/mb01td.c
    MB/mb01uy.c
    MB/mb01vd.c
    MB/mb02qy.c
    MB/mb02ud.c
    MB/mb02yd.c
    MB/mb03od.c
    MB/mb03oy.c
    MB/mb04id.c
    MB/mb04iy.c
    MB/mb03ud.c
    MB/mb04kd.c
    MB/mb04nd.c
    MB/mb04ny.c
    MB/mb04od.c
    MB/mb04ow.c
    MB/mb04oy.c
    MB/mb04py.c
    MB/mb04zd.c
    MB/mb05nd.c
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
    TB/tb01id.c
    TB/tb01vd.c
    TB/tb01vy.c
    TB/tb01wd.c
    TC/tc01od.c
    TF/tf01md.c
    TF/tf01mx.c
    TF/tf01rd.c
    TG/tg01ad.c
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
