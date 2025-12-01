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

# Include diagnostic target functions
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/DiagnosticTargets.cmake)

message(STATUS "")
message(STATUS "Adding diagnostic targets:")

# Add diagnostics for implemented routines
add_slicot_diagnostic(SG03BD)
add_slicot_diagnostic(AB13MD)

# Convenience targets (backward compatibility with old naming)
add_custom_target(run_fortran_diag DEPENDS run_sg03bd_fortran)
add_custom_target(run_c_diag DEPENDS run_sg03bd_c)
add_custom_target(compare_diag DEPENDS compare_sg03bd)
add_custom_target(diag_all DEPENDS sg03bd_diag_all)

# Print usage instructions
message(STATUS "")
message(STATUS "Fortran Diagnostic Tools - Usage:")
message(STATUS "  cmake --build . --target diag_all             # Run full workflow (recommended)")
message(STATUS "")
message(STATUS "Per-routine targets (example for SG03BD):")
message(STATUS "  cmake --build . --target sg03bd_diag_c        # Build C diagnostic")
message(STATUS "  cmake --build . --target sg03bd_diag_fortran  # Build Fortran diagnostic")
message(STATUS "  cmake --build . --target run_sg03bd_c         # Run C")
message(STATUS "  cmake --build . --target run_sg03bd_fortran   # Run Fortran")
message(STATUS "  cmake --build . --target compare_sg03bd       # Compare results")
message(STATUS "  cmake --build . --target sg03bd_diag_all      # Full workflow for SG03BD")
message(STATUS "")
message(STATUS "Backward compatibility aliases:")
message(STATUS "  run_fortran_diag, run_c_diag, compare_diag")
message(STATUS "")
message(STATUS "Output location: ${CMAKE_CURRENT_BINARY_DIR}/output/<routine>_diff.txt")
message(STATUS "=================================")
message(STATUS "")
