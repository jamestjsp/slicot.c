# CMake function for creating SLICOT diagnostic targets
# Reduces boilerplate from 76 lines to 1 line per routine
#
# Usage:
#   add_slicot_diagnostic(SG03BD)
#   add_slicot_diagnostic(AB01MD)

function(add_slicot_diagnostic ROUTINE)
    # Convert routine name to lowercase for file names
    string(TOLOWER ${ROUTINE} routine_lower)

    # Test data paths
    set(DATA_FILE "${SLICOT_FORTRAN_DIR}/examples/data/${ROUTINE}.dat")
    set(EXPECTED_FILE "${SLICOT_FORTRAN_DIR}/examples/results/${ROUTINE}.res")

    # C diagnostic executable
    add_executable(${routine_lower}_diag_c c/${routine_lower}_diag.c)
    target_link_libraries(${routine_lower}_diag_c PRIVATE
        slicot
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
        m
    )
    target_include_directories(${routine_lower}_diag_c PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_BINARY_DIR}/include
    )

    # Fortran diagnostic executable
    add_executable(${routine_lower}_diag_fortran fortran/${routine_lower}_diag.f)
    target_link_libraries(${routine_lower}_diag_fortran PRIVATE
        slicot_fortran
        ${BLAS_LIBRARIES}
        ${LAPACK_LIBRARIES}
    )

    # Custom target: Run Fortran diagnostic
    add_custom_target(run_${routine_lower}_fortran
        COMMAND ${CMAKE_COMMAND} -E echo "Running ${ROUTINE} Fortran diagnostic..."
        COMMAND ${routine_lower}_diag_fortran < ${DATA_FILE} > ${routine_lower}_fortran.txt 2>&1 || true
        DEPENDS ${routine_lower}_diag_fortran
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Running ${ROUTINE} Fortran diagnostic"
    )

    # Custom target: Run C diagnostic
    add_custom_target(run_${routine_lower}_c
        COMMAND ${CMAKE_COMMAND} -E echo "Running ${ROUTINE} C diagnostic..."
        COMMAND ${routine_lower}_diag_c < ${DATA_FILE} > ${routine_lower}_c.txt 2>&1 || true
        DEPENDS ${routine_lower}_diag_c
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Running ${ROUTINE} C diagnostic"
    )

    # Custom target: Compare outputs
    add_custom_target(compare_${routine_lower}
        COMMAND ${CMAKE_COMMAND} -E make_directory output
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${routine_lower}_fortran.txt output/
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${routine_lower}_c.txt output/
        COMMAND ${CMAKE_COMMAND} -E echo "Comparing ${ROUTINE} outputs..."
        COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/scripts/compare.py
                ${routine_lower}_fortran.txt ${routine_lower}_c.txt > output/${routine_lower}_diff.txt 2>&1 || true
        COMMAND ${CMAKE_COMMAND} -E echo "Comparison saved to: ${CMAKE_CURRENT_BINARY_DIR}/output/${routine_lower}_diff.txt"
        DEPENDS run_${routine_lower}_fortran run_${routine_lower}_c
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Comparing ${ROUTINE} Fortran vs C outputs"
    )

    # Custom target: Full workflow for this routine
    add_custom_target(${routine_lower}_diag_all
        DEPENDS compare_${routine_lower}
        COMMENT "Complete ${ROUTINE} diagnostic workflow"
    )

    # Print status
    message(STATUS "  Added diagnostic targets for ${ROUTINE}")
    message(STATUS "    - ${routine_lower}_diag_c, ${routine_lower}_diag_fortran")
    message(STATUS "    - run_${routine_lower}_c, run_${routine_lower}_fortran")
    message(STATUS "    - compare_${routine_lower}, ${routine_lower}_diag_all")
endfunction()
