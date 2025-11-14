cmake_minimum_required(VERSION 3.14)
project(slicot_c VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(BUILD_TESTING "Build tests" ON)
option(BUILD_EXAMPLES "Build examples" ON)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)
option(SLICOT_ILP64 "Use 64-bit integers for BLAS/LAPACK" OFF)

# Find BLAS/LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Detect BLAS/LAPACK symbol naming convention
include(CheckFunctionExists)
set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

# Try lowercase with underscore (most common)
check_function_exists(dgemm_ HAVE_DGEMM_UNDERSCORE)
if(HAVE_DGEMM_UNDERSCORE)
    set(SLC_FC_LOWER_US ON)
    message(STATUS "BLAS/LAPACK uses lowercase with underscore: dgemm_")
else()
    # Try lowercase only
    check_function_exists(dgemm HAVE_DGEMM_LOWER)
    if(HAVE_DGEMM_LOWER)
        set(SLC_FC_LOWER ON)
        message(STATUS "BLAS/LAPACK uses lowercase: dgemm")
    else()
        # Try uppercase
        check_function_exists(DGEMM HAVE_DGEMM_UPPER)
        if(HAVE_DGEMM_UPPER)
            set(SLC_FC_UPPER ON)
            message(STATUS "BLAS/LAPACK uses uppercase: DGEMM")
        else()
            # Default to lowercase with underscore
            set(SLC_FC_LOWER_US ON)
            message(WARNING "Could not detect BLAS/LAPACK naming convention, defaulting to dgemm_")
        endif()
    endif()
endif()

# Generate config header
configure_file(
    ${PROJECT_SOURCE_DIR}/include/slicot_config.h.in
    ${PROJECT_BINARY_DIR}/include/slicot_config.h
)

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_BINARY_DIR}/include)

add_subdirectory(src)

if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# Python bindings
if(BUILD_PYTHON_BINDINGS)
    find_package(Python3 COMPONENTS Interpreter Development NumPy)
    if(Python3_FOUND AND Python3_NumPy_FOUND)
        message(STATUS "Python3 found: ${Python3_EXECUTABLE}")
        message(STATUS "NumPy found: ${Python3_NumPy_INCLUDE_DIRS}")

        # Build Python extension
        add_custom_target(python_build ALL
            COMMAND ${Python3_EXECUTABLE} setup.py build_ext --inplace
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            DEPENDS slicot
            COMMENT "Building Python bindings"
        )

        # Install Python package (optional)
        install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install -e ${PROJECT_SOURCE_DIR})")
    else()
        message(WARNING "Python3 or NumPy not found. Python bindings disabled.")
    endif()
endif()
