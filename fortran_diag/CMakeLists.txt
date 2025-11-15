# Fortran diagnostic tools for debugging C vs Fortran implementations
# Part of main SLICOT.C build - enable with: cmake -DBUILD_FORTRAN_DIAG=ON

message(STATUS "")
message(STATUS "=== Fortran Diagnostic Tools ===")

# Enable Fortran language
enable_language(Fortran)

# SLICOT-Reference via ExternalProject (complete isolation, avoids target name conflicts)
set(SLICOT_FORTRAN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../SLICOT-Reference" CACHE PATH "Path to SLICOT Fortran reference")
set(SLICOT_FORTRAN_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/slicot_fortran_external")
set(SLICOT_FORTRAN_LIB "${SLICOT_FORTRAN_BUILD_DIR}/lib/libslicot.a")

include(ExternalProject)
ExternalProject_Add(slicot_fortran_external
    SOURCE_DIR ${SLICOT_FORTRAN_DIR}
    BINARY_DIR ${SLICOT_FORTRAN_BUILD_DIR}
    CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
        -DBUILD_SHARED_LIBS=OFF
        -DSLICOT_TESTING=OFF
    BUILD_BYPRODUCTS ${SLICOT_FORTRAN_LIB}
    INSTALL_COMMAND ""  # Skip install
    STEP_TARGETS build
)

# Import as CMake library
add_library(slicot_fortran STATIC IMPORTED)
set_target_properties(slicot_fortran PROPERTIES
    IMPORTED_LOCATION ${SLICOT_FORTRAN_LIB}
)
add_dependencies(slicot_fortran slicot_fortran_external)

message(STATUS "SLICOT Fortran reference via ExternalProject (${SLICOT_FORTRAN_LIB})")

# C diagnostic program
add_executable(sg03bd_diag_c c/sg03bd_diag.c)
target_link_libraries(sg03bd_diag_c PRIVATE slicot ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} m)
target_include_directories(sg03bd_diag_c PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}/include
)

# Fortran diagnostic program
add_executable(sg03bd_diag_fortran fortran/sg03bd_diag.f)
target_link_libraries(sg03bd_diag_fortran PRIVATE
    slicot_fortran
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
)

# Test data paths
set(SG03BD_TEST_DATA "${SLICOT_FORTRAN_DIR}/examples/data/SG03BD.dat")
set(SG03BD_EXPECTED "${SLICOT_FORTRAN_DIR}/examples/results/SG03BD.res")

# Custom target: Run Fortran diagnostic
add_custom_target(run_fortran_diag
    COMMAND ${CMAKE_COMMAND} -E echo "Running Fortran diagnostic..."
    COMMAND sg03bd_diag_fortran < ${SG03BD_TEST_DATA} > fortran_trace.txt 2>&1 || true
    DEPENDS sg03bd_diag_fortran
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running Fortran diagnostic with SG03BD.dat"
)

# Custom target: Run C diagnostic
add_custom_target(run_c_diag
    COMMAND ${CMAKE_COMMAND} -E echo "Running C diagnostic..."
    COMMAND sg03bd_diag_c < ${SG03BD_TEST_DATA} > c_trace.txt 2>&1 || true
    DEPENDS sg03bd_diag_c
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running C diagnostic with SG03BD.dat"
)

# Custom target: Compare outputs
add_custom_target(compare_diag
    COMMAND ${CMAKE_COMMAND} -E make_directory output
    COMMAND ${CMAKE_COMMAND} -E copy_if_different fortran_trace.txt output/
    COMMAND ${CMAKE_COMMAND} -E copy_if_different c_trace.txt output/
    COMMAND ${CMAKE_COMMAND} -E echo "Comparing outputs..."
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare.py
            fortran_trace.txt c_trace.txt > output/diff.txt 2>&1 || true
    COMMAND ${CMAKE_COMMAND} -E echo "Comparison saved to: ${CMAKE_CURRENT_BINARY_DIR}/output/diff.txt"
    DEPENDS run_fortran_diag run_c_diag
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Comparing Fortran vs C diagnostic outputs"
)

# Custom target: Run full diagnostic workflow
add_custom_target(diag_all
    DEPENDS compare_diag
    COMMENT "Complete diagnostic workflow"
)

# Print usage instructions
message(STATUS "")
message(STATUS "Fortran Diagnostic Tools - Usage:")
message(STATUS "  cmake --build . --target sg03bd_diag_c        # Build C diagnostic")
message(STATUS "  cmake --build . --target sg03bd_diag_fortran  # Build Fortran diagnostic")
message(STATUS "  cmake --build . --target run_fortran_diag     # Run Fortran")
message(STATUS "  cmake --build . --target run_c_diag           # Run C")
message(STATUS "  cmake --build . --target compare_diag         # Compare results")
message(STATUS "  cmake --build . --target diag_all             # Run full workflow")
message(STATUS "")
message(STATUS "Output location: ${CMAKE_CURRENT_BINARY_DIR}/output/diff.txt")
message(STATUS "=================================")
message(STATUS "")
