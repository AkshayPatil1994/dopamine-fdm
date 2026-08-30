cmake_minimum_required(VERSION 3.20)
project(dopamine LANGUAGES Fortran)

#------------------------------------------------------------------------------
# GPU (OpenACC via nvfortran) — optional, off by default; CPU build is untouched
#------------------------------------------------------------------------------
option(ENABLE_GPU "Build with OpenACC GPU offload (requires nvfortran)" OFF)
if(ENABLE_GPU AND NOT CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
    message(FATAL_ERROR "ENABLE_GPU requires nvfortran: configure with -DCMAKE_Fortran_COMPILER=nvfortran")
endif()

#------------------------------------------------------------------------------
# MPI
#------------------------------------------------------------------------------
find_package(MPI REQUIRED COMPONENTS Fortran)

#------------------------------------------------------------------------------
# BLAS / LAPACK
#------------------------------------------------------------------------------
find_package(LAPACK REQUIRED)
find_package(BLAS   REQUIRED)

#------------------------------------------------------------------------------
# FFTW3 (local single-rank transforms/DCT only -- MPI-parallel transforms go
# through 2decomp&fft's pencil transposes instead, see src/projection.f90)
#------------------------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(FFTW REQUIRED COMPONENTS DOUBLE)

#------------------------------------------------------------------------------
# 2decomp&fft (2D pencil decomposition), vendored as a configure-time
# sub-build rather than FetchContent+add_subdirectory: 2decomp&fft's own
# CMakeLists.txt hardcodes ${CMAKE_SOURCE_DIR} (not
# ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_SOURCE_DIR}) for its own module path
# and config-file templates, which resolves to *our* source tree rather than
# theirs when embedded via add_subdirectory -- breaks its include(D2D_MPI)
# and configure_file(.../decomp2d-config.cmake.in) calls. Building it as an
# independent CMake invocation (its own CMAKE_SOURCE_DIR) sidesteps that;
# this mirrors the same synchronous configure+build pattern 2decomp&fft's
# own cmake/fft/fft.cmake uses to fetch+build its FindFFTW.cmake dependency.
#------------------------------------------------------------------------------
set(DECOMP2D_TAG        "v2.1.0")
set(DECOMP2D_SRC_DIR    "${CMAKE_BINARY_DIR}/_deps/2decomp_fft-src")
set(DECOMP2D_BIN_DIR    "${CMAKE_BINARY_DIR}/_deps/2decomp_fft-build")
set(DECOMP2D_INSTALL_DIR "${CMAKE_BINARY_DIR}/_deps/2decomp_fft-install")

if(NOT EXISTS "${DECOMP2D_SRC_DIR}/CMakeLists.txt")
    execute_process(
        COMMAND git clone --branch ${DECOMP2D_TAG} --depth 1
                https://github.com/2decomp-fft/2decomp-fft.git "${DECOMP2D_SRC_DIR}"
        RESULT_VARIABLE _d2d_clone_result
    )
    if(_d2d_clone_result)
        message(FATAL_ERROR "Failed to clone 2decomp&fft ${DECOMP2D_TAG}")
    endif()
endif()

if(NOT EXISTS "${DECOMP2D_INSTALL_DIR}/lib/decomp2d/decomp2d-config.cmake")
    execute_process(
        COMMAND ${CMAKE_COMMAND}
                -S "${DECOMP2D_SRC_DIR}" -B "${DECOMP2D_BIN_DIR}"
                -DCMAKE_INSTALL_PREFIX=${DECOMP2D_INSTALL_DIR}
                -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
                -DBUILD_TARGET=mpi -DFFT_Choice=fftw -DBUILD_TESTING=OFF
                -DCMAKE_MODULE_PATH=${CMAKE_CURRENT_SOURCE_DIR}/cmake
        RESULT_VARIABLE _d2d_cfg_result
    )
    if(_d2d_cfg_result)
        message(FATAL_ERROR "2decomp&fft configure failed")
    endif()
    execute_process(
        COMMAND ${CMAKE_COMMAND} --build "${DECOMP2D_BIN_DIR}" --target install -j${CMAKE_BUILD_PARALLEL_LEVEL}
        RESULT_VARIABLE _d2d_build_result
    )
    if(_d2d_build_result)
        message(FATAL_ERROR "2decomp&fft build/install failed")
    endif()
endif()

find_package(decomp2d REQUIRED CONFIG
    PATHS "${DECOMP2D_INSTALL_DIR}/lib/decomp2d"
    NO_DEFAULT_PATH
)

#------------------------------------------------------------------------------
# Compiler flags
#------------------------------------------------------------------------------
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -O3 -march=native -funroll-loops -std=gnu -fallow-argument-mismatch -fconvert=big-endian")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -xHost -ipo -convert big_endian")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -O3 -Mbyteswapio")
    if(ENABLE_GPU)
        # GPU_POISSON: single-GPU cuFFT periodic Poisson solve (nprocs==1 only), see
        # docs/Installation.md#building-single-gpu-openacc--cufft--cusparse
        # rdc (relocatable device code): required so nvlink can resolve cross-procedure
        # !$acc routine seq calls (e.g. wallmodel/ibm device-callable helpers)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -acc -gpu=cc86,rdc -Minfo=accel -DGPU_POISSON -cudalib=cufft,cusparse")
    endif()
endif()

#------------------------------------------------------------------------------
# Source files  (order must respect Fortran module USE dependencies)
#------------------------------------------------------------------------------
set(SRCS
    src/mpi.f90
    src/profiler.f90
    src/global.f90
    src/decomp.f90
    src/interpolation.f90
    src/sem.f90
    src/boundary_conditions.f90
    src/sgs_model.f90
    src/equations.f90
    src/readMask.f90
    src/ibm.f90
    src/scalar_transport.f90
    src/thermal_transport.f90
    src/genGridandIC.f90
    src/input_output.f90
    src/initialization.f90
    src/wallmodel.f90
    src/finalization.f90
)

if(ENABLE_GPU)
    list(APPEND SRCS src/poisson_gpu.f90)
endif()

list(APPEND SRCS
    src/projection.f90
    src/time_integration.f90
    src/monitor.f90
    src/reynolds_stress_budget.f90
    src/probe_output.f90
    src/main.f90
)

add_executable(dopamine ${SRCS})

# Fortran module output directory: dopamine and the verify_* targets below
# compile their own copies of src/mpi.f90, src/global.f90 and src/sem.f90, so
# without a per-target Fortran_MODULE_DIRECTORY they'd all write .mod files
# into the same shared build directory -- a parallel-build (-j) race ("Error
# copying Fortran module ...") since two targets' compiles of the same source
# can interleave writes to the same .mod path.
set_target_properties(dopamine PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/dopamine")

if(ENABLE_GPU)
    # projection.f90's `Use poisson_gpu` sits inside an #ifdef GPU_POISSON block;
    # CMake's Fortran dependency scanner doesn't run the preprocessor, so it never
    # sees this and won't order poisson_gpu.f90 before projection.f90 -- a -j race
    # ("Unable to open MODULE file poisson_gpu.mod") when they land in the same
    # parallel batch. State the dependency explicitly.
    set_source_files_properties(src/projection.f90 PROPERTIES
        OBJECT_DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/dopamine.dir/src/poisson_gpu.f90.o")
endif()

target_include_directories(dopamine PRIVATE
    ${FFTW_INCLUDE_DIRS}
)

# BEFORE: our own src/mpi.f90 module is named 'mpi', same as OpenMPI's real
# Fortran module -- gfortran searches -I paths in command-line order, so this
# must precede MPI::MPI_Fortran's include dir or 'Use mpi' resolves to the
# system module instead of ours (and none of our declared symbols exist there).
target_include_directories(dopamine BEFORE PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}/mod/dopamine"
)

target_link_libraries(dopamine
    MPI::MPI_Fortran
    ${FFTW_LIBRARIES}          # fftw3 (local single-rank transforms/DCT only)
    LAPACK::LAPACK
    BLAS::BLAS
    decomp2d
)

#------------------------------------------------------------------------------
# Standalone verification driver for the divergence-free SEM construction
# (src/sem.f90's raw_eddy_curl_sum, sem_divergence_free=1); not part of the dopamine
# executable, see docs/Numerics.md#111-inlet--constant-or-synthetic-eddy-method
#------------------------------------------------------------------------------
add_executable(verify_sem_divergence_free
    src/mpi.f90
    src/global.f90
    src/sem.f90
    tests/verify_sem_divergence_free.f90
)

set_target_properties(verify_sem_divergence_free PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/verify_sem_divergence_free")

target_include_directories(verify_sem_divergence_free PRIVATE
    ${FFTW_INCLUDE_DIRS}
)

target_include_directories(verify_sem_divergence_free BEFORE PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}/mod/verify_sem_divergence_free"
)

target_link_libraries(verify_sem_divergence_free
    MPI::MPI_Fortran
    ${FFTW_LIBRARIES}
    LAPACK::LAPACK
    BLAS::BLAS
)

#------------------------------------------------------------------------------
# Standalone verification driver for the wind-tunnel TI inflow profile parser
# (src/sem.f90's read_mean_profile_TI); not part of the dopamine executable, see
# docs/Numerics.md#111-inlet--constant-or-synthetic-eddy-method
#------------------------------------------------------------------------------
add_executable(verify_ti_profile
    src/mpi.f90
    src/global.f90
    src/sem.f90
    tests/verify_ti_profile.f90
)

set_target_properties(verify_ti_profile PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/verify_ti_profile")

target_include_directories(verify_ti_profile PRIVATE
    ${FFTW_INCLUDE_DIRS}
)

target_include_directories(verify_ti_profile BEFORE PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}/mod/verify_ti_profile"
)

target_link_libraries(verify_ti_profile
    MPI::MPI_Fortran
    ${FFTW_LIBRARIES}
    LAPACK::LAPACK
    BLAS::BLAS
)

#------------------------------------------------------------------------------
# Install
#------------------------------------------------------------------------------
install(TARGETS dopamine DESTINATION bin)
