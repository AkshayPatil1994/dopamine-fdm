cmake_minimum_required(VERSION 3.20)
project(dopamine LANGUAGES Fortran)

#------------------------------------------------------------------------------
# Testing (CTest) -- `ctest` after building runs the verify_* unit drivers
# plus the tgv_small MPI-rank-count parity regression (see tests/regression/),
# the same harness later reused for CPU-vs-GPU and 1-GPU-vs-N-GPU parity.
#------------------------------------------------------------------------------
enable_testing()

#------------------------------------------------------------------------------
# GPU (OpenACC via nvfortran) — optional, off by default; CPU build is untouched
#------------------------------------------------------------------------------
option(ENABLE_GPU "Build with OpenACC GPU offload (requires nvfortran)" OFF)
if(ENABLE_GPU AND NOT CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
    message(FATAL_ERROR "ENABLE_GPU requires nvfortran: configure with -DCMAKE_Fortran_COMPILER=nvfortran")
endif()
# cc86 = compute capability 8.6 (RTX A6000/3090/A40 etc.); override with
# -DCMAKE_CUDA_ARCHITECTURES=<cc> for other GPUs. Feeds both this project's
# own -gpu=cc<N> nvfortran flag below and 2decomp&fft's gpu-target sub-build.
# Optional CUDA toolkit pin for nvfortran (-gpu=cudaX.Y). Needed when the
# installed NVIDIA driver is older than the toolkit NVHPC bundles by default
# (e.g. driver 535 supports CUDA <=12.2 only; the 12.9 default segfaults in
# cuModuleLoad at startup). Empty = use nvfortran's default toolkit.
set(GPU_CUDA_VERSION "" CACHE STRING "CUDA toolkit version for nvfortran -gpu=cudaX.Y, e.g. 12.2 (empty = compiler default)")
option(DECOMP2D_GPU "With ENABLE_GPU, build 2decomp&fft in its GPU (CUDA-aware transpose) mode" ON)
if(ENABLE_GPU AND NOT CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES "86" CACHE STRING "CUDA compute capability (no dots, e.g. 86)")
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
    # ENABLE_GPU: build 2decomp&fft's own gpu target (BUILD_TARGET=gpu,
    # ENABLE_CUDA=ON) instead of mpi -- its transpose_x_to_y/y_to_z/etc.
    # (the same calls src/projection.f90 and, from here on, src/poisson_gpu.f90
    # use) compile to CUDA-aware-MPI device-to-device transposes under
    # `#if defined(_GPU)` (see its src/transpose_*.f90), instead of the host
    # staging path -- this is what makes a >nprocs=1 GPU build possible at
    # all without hand-writing our own GPU-aware pencil transpose code.
    set(_d2d_target_args -DBUILD_TARGET=mpi -DFFT_Choice=fftw)
    set(_d2d_fc "${CMAKE_Fortran_COMPILER}")
    if(ENABLE_GPU AND DECOMP2D_GPU)
        set(_d2d_target_args -DBUILD_TARGET=gpu -DENABLE_CUDA=ON -DFFT_Choice=fftw
                              -DCMAKE_CUDA_ARCHITECTURES=${CMAKE_CUDA_ARCHITECTURES})
        # 2decomp&fft replaces CMAKE_Fortran_FLAGS with its own D2D_FFLAGS, so extra
        # -gpu options can't be passed as flags; inject them via a compiler wrapper.
        # rdc must match this project's own -gpu=...,rdc: a non-rdc lib linked into an
        # rdc app crashes at startup, and the lib's CUDA version must not exceed ours.
        set(_d2d_gpu_opts "rdc")
        if(GPU_CUDA_VERSION)
            set(_d2d_gpu_opts "cuda${GPU_CUDA_VERSION},rdc")
        endif()
        set(_d2d_fc_wrapper "${CMAKE_BINARY_DIR}/nvfortran_d2d.sh")
        file(WRITE "${_d2d_fc_wrapper}" "#!/bin/sh\nexec ${CMAKE_Fortran_COMPILER} \"$@\" -gpu=${_d2d_gpu_opts}\n")
        file(CHMOD "${_d2d_fc_wrapper}" PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
        set(_d2d_fc "${_d2d_fc_wrapper}")
    endif()
    execute_process(
        COMMAND ${CMAKE_COMMAND}
                -S "${DECOMP2D_SRC_DIR}" -B "${DECOMP2D_BIN_DIR}"
                -DCMAKE_INSTALL_PREFIX=${DECOMP2D_INSTALL_DIR}
                -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                -DCMAKE_Fortran_COMPILER=${_d2d_fc}
                -DMPI_Fortran_COMPILER=${MPI_Fortran_COMPILER}
                ${_d2d_target_args} -DBUILD_TESTING=OFF
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
        set(_gpu_cuda_flag "")
        # 2decomp&fft's GPU lib is CUDA Fortran (-cuda); the app must link with the same language mode
        set(_gpu_cuda_lang_flag "")
        if(DECOMP2D_GPU)
            set(_gpu_cuda_lang_flag "-cuda ")
        endif()
        if(GPU_CUDA_VERSION)
            set(_gpu_cuda_flag "cuda${GPU_CUDA_VERSION},")
        endif()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -acc ${_gpu_cuda_lang_flag}-gpu=${_gpu_cuda_flag}cc${CMAKE_CUDA_ARCHITECTURES},rdc -Minfo=accel -DGPU_POISSON -cudalib=cufft,cusparse")
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
    src/halo_pad.f90
    src/waves.f90
    src/uav_actuator.f90
    src/interpolation.f90
    src/sem.f90
    src/particles.f90
    src/boundary_conditions.f90
    src/sgs_model.f90
    src/equations.f90
    src/readMask.f90
    src/ibm.f90
    src/scalar_transport.f90
    src/thermal_transport.f90
    src/vof_plic.f90
    src/vof_normals.f90
    src/vof_advect.f90
    src/genGridandIC.f90
    src/input_output.f90
    src/debug_trace.f90
)

if(ENABLE_GPU)
    # gpu_device.f90 before initialization.f90: initialization.f90's
    # `Use gpu_device` sits inside an #ifdef GPU_POISSON block (like
    # projection.f90's `Use poisson_gpu` below), so it must precede
    # initialization.f90 in SRCS for the (preprocessor-blind) module scanner
    # to order the build correctly; see the matching OBJECT_DEPENDS below.
    list(APPEND SRCS src/gpu_device.f90)
endif()

list(APPEND SRCS
    src/initialization.f90
    src/wallmodel.f90
    src/ibm_stress.f90
    src/finalization.f90
)

if(ENABLE_GPU)
    list(APPEND SRCS src/poisson_gpu.f90)
endif()

list(APPEND SRCS
    src/projection.f90
    src/mom_recon.f90
    src/vof_pressure.f90
    src/vof_state.f90
    src/vof_curv.f90
    src/monitor.f90
    src/vof_ibm.f90
    src/vof_twofluid.f90
    src/time_integration.f90
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
    set_source_files_properties(src/initialization.f90 PROPERTIES
        OBJECT_DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/dopamine.dir/src/gpu_device.f90.o")
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
# Standalone check of the wave generators (src/waves.f90): dispersion, stream-function surface conditions, JONSWAP Hs
#------------------------------------------------------------------------------
add_executable(verify_waves
    src/mpi.f90
    src/global.f90
    src/waves.f90
    tests/verify_waves.f90
)
set_target_properties(verify_waves PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/verify_waves")
target_include_directories(verify_waves PRIVATE ${FFTW_INCLUDE_DIRS})
target_include_directories(verify_waves BEFORE PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/mod/verify_waves")
target_link_libraries(verify_waves MPI::MPI_Fortran ${FFTW_LIBRARIES} LAPACK::LAPACK BLAS::BLAS)
add_test(NAME waves_generators COMMAND $<TARGET_FILE:verify_waves>)

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
# Standalone verification driver for the Bezier-parametrized inflow-optimization
# helpers (src/sem.f90's bezier_point/bezier_eval/init_inflow_opt/
# solve_inflow_opt_step2); not part of the dopamine executable, see
# docs/Numerics.md#111-inlet--constant-or-synthetic-eddy-method
#------------------------------------------------------------------------------
add_executable(verify_inflow_opt
    src/mpi.f90
    src/global.f90
    src/sem.f90
    tests/verify_inflow_opt.f90
)

set_target_properties(verify_inflow_opt PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/verify_inflow_opt")

target_include_directories(verify_inflow_opt PRIVATE
    ${FFTW_INCLUDE_DIRS}
)

target_include_directories(verify_inflow_opt BEFORE PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}/mod/verify_inflow_opt"
)

target_link_libraries(verify_inflow_opt
    MPI::MPI_Fortran
    ${FFTW_LIBRARIES}
    LAPACK::LAPACK
    BLAS::BLAS
)

#------------------------------------------------------------------------------
# dopamine-ESEM: standalone ESEM inflow precursor generator (precomputes an
# x-normal donor slice time series in probe_output.f90's format, consumed by
# a separate inflow_type=2 run -- see src/dopamine_esem_main.f90). Needs the
# real domain decomposition (decomp.f90/2decomp&fft) for MPI-parallel z-slab
# generation, and genGridAndIC/input_output for the same grid + namelist
# reading a live solver run uses; does not need equations/sgs/ibm/poisson.
#------------------------------------------------------------------------------
add_executable(dopamine-ESEM
    src/mpi.f90
    src/profiler.f90
    src/global.f90
    src/decomp.f90
    src/halo_pad.f90
    src/sem.f90
    src/particles.f90
    src/genGridandIC.f90
    src/input_output.f90
    src/dopamine_esem_main.f90
)

set_target_properties(dopamine-ESEM PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/dopamine-ESEM")

target_include_directories(dopamine-ESEM PRIVATE
    ${FFTW_INCLUDE_DIRS}
)

target_include_directories(dopamine-ESEM BEFORE PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}/mod/dopamine-ESEM"
)

target_link_libraries(dopamine-ESEM
    MPI::MPI_Fortran
    ${FFTW_LIBRARIES}
    LAPACK::LAPACK
    BLAS::BLAS
    decomp2d
)

#------------------------------------------------------------------------------
# Standalone verification drivers for the VOF/PLIC geometry and advection kernels
# (src/vof_plic.f90, vof_normals.f90, vof_advect.f90); serial, no solver dependencies
#------------------------------------------------------------------------------
set(VOF_SRCS src/vof_plic.f90 src/vof_normals.f90 src/vof_advect.f90)
add_executable(test_plic_geometry src/vof_plic.f90 tests/vof/test_plic_geometry.f90)
set_target_properties(test_plic_geometry PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/test_plic_geometry")
add_executable(test_vof_advect ${VOF_SRCS} tests/vof/test_vof_advect.f90)
set_target_properties(test_vof_advect PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/test_vof_advect")

#------------------------------------------------------------------------------
# Tests
#------------------------------------------------------------------------------
add_test(NAME verify_sem_divergence_free
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:verify_sem_divergence_free>)
add_test(NAME verify_ti_profile
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:verify_ti_profile>)
add_test(NAME verify_inflow_opt
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:verify_inflow_opt>)
add_test(NAME vof_plic_geometry COMMAND $<TARGET_FILE:test_plic_geometry>)
add_test(NAME vof_advect_prescribed COMMAND $<TARGET_FILE:test_vof_advect>)

# Parity regressions (tests/regression/<case>): deterministic small cases whose output_monitor
# diagnostics (src/monitor.f90:61) must agree within tolerance across rank counts and, in
# ENABLE_GPU builds, between CPU and GPU. See tests/regression/compare_monitor.py.
#   tgv_small      fully periodic Taylor-Green vortex (analytic IC)
#   chan_small     periodic-x/z channel, no-slip y walls, DNS (wall-normal tridiagonal Poisson path)
#   chan_les_small same with Vreman LES + equilibrium wall model on a stretched y grid
#   duct_small     4-wall duct (z walls, alpha_grid_z stretched z), DNS+Vreman+EQWM (eigenmode Poisson path)
#   sem_small      inflow/outflow x (DCT-IV Poisson path), uniform inflow
#   sed_small      chan_les_small + passive settling scalar (&SEDIMENT)
#   bouss_small    chan_les_small + Boussinesq temperature (&BOUSSINESQ)
#   uav_small      chan_les_small + static UAV actuator disk on the np=2 rank seam (radial loading + swirl)
#   uavpath_small  same, disk following a path across the seam with the tilt filter
#   ibm_small      LES channel + ghost-cell IBM sphere (synthetic SDF_in, 34x17x34) straddling the rank seam
#   ibm_wm_small   ibm_small with the IBM wall model on (ibm_wall_model_flag=1)
#   rsb_small      chan_les_small + Reynolds-stress-budget statistics (stats/rsb_*.bin compared too)
# The tests also compare the final snapshot field-by-field, ghost layers excluded (compare_monitor.py --fields),
# and restart == uninterrupted (restart_check.py), not just the monitor columns. nsave=nsteps in every case.
# Every build: np=1 vs np=TEST_PARITY_NPROCS of this same executable. GPU build additionally
# compares against a CPU build if CPU_REFERENCE_EXE (+ CPU_REFERENCE_MPIRUN) point at one.
set(TEST_PARITY_NPROCS "2" CACHE STRING
    "Second rank count for the MPI parity regressions (compared against np=1)")
set(CPU_REFERENCE_EXE "" CACHE FILEPATH "CPU-build dopamine executable to compare GPU runs against (optional)")
set(CPU_REFERENCE_MPIRUN "" CACHE FILEPATH "mpirun matching CPU_REFERENCE_EXE (default: this build's mpiexec)")
set(TEST_LD_LIBRARY_PATH "" CACHE STRING "LD_LIBRARY_PATH prepended for parity tests (e.g. a matching user-space libcuda)")
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    foreach(case tgv_small chan_small chan_les_small duct_small sem_small sed_small bouss_small uav_small uavpath_small ibm_small ibm_wm_small ibm_stair_small rsb_small)
        set(_case_dir ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/${case})
        set(_cmp ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/compare_monitor.py
                 --case-dir ${_case_dir} --mpirun ${MPIEXEC_EXECUTABLE} --fields --field-tol 1e-8)
        set(_kind cpu)
        if(ENABLE_GPU)
            set(_kind gpu)
        endif()
        set(_ftol 1e-8)
        set(_rtol 1e-12)
        set(_extra "")
        if(case STREQUAL "rsb_small")
            set(_extra --files stats/rsb_*.bin)
        endif()
        add_test(NAME ${case}_${_kind}_np1_vs_np${TEST_PARITY_NPROCS}
            COMMAND ${_cmp} --exe-a $<TARGET_FILE:dopamine> --exe-b $<TARGET_FILE:dopamine>
                    --np-a 1 --np-b ${TEST_PARITY_NPROCS} --field-tol ${_ftol} ${_extra}
                    --label-a ${_kind}-np1 --label-b ${_kind}-np${TEST_PARITY_NPROCS})
        if(TEST_LD_LIBRARY_PATH)
            set_tests_properties(${case}_${_kind}_np1_vs_np${TEST_PARITY_NPROCS} PROPERTIES
                ENVIRONMENT "LD_LIBRARY_PATH=${TEST_LD_LIBRARY_PATH}:$ENV{LD_LIBRARY_PATH}")
        endif()
        # Explicit 4-rank pencil layouts, so the x and z halos interact (2x2: corner planes) and x is split without z (4x1);
        # duct needs p_col=1.
        set(_layouts "2,2" "4,1")
        if(case STREQUAL "duct_small")
            set(_layouts "4,1")
        endif()
        foreach(_lay ${_layouts})
            string(REPLACE "," "x" _laytag ${_lay})
            if(TRUE)
                add_test(NAME ${case}_${_kind}_np1_vs_np4_${_laytag}
                    COMMAND ${_cmp} --exe-a $<TARGET_FILE:dopamine> --exe-b $<TARGET_FILE:dopamine>
                            --np-a 1 --np-b 4 --p-grid ${_lay} --p-grid-b-only --field-tol ${_ftol} ${_extra}
                            --label-a ${_kind}-np1 --label-b ${_kind}-${_laytag})
                if(TEST_LD_LIBRARY_PATH)
                    set_tests_properties(${case}_${_kind}_np1_vs_np4_${_laytag} PROPERTIES
                        ENVIRONMENT "LD_LIBRARY_PATH=${TEST_LD_LIBRARY_PATH}:$ENV{LD_LIBRARY_PATH}")
                endif()
            endif()
            if(ENABLE_GPU AND CPU_REFERENCE_EXE)
                set(_ref_mpirun_l ${MPIEXEC_EXECUTABLE})
                if(CPU_REFERENCE_MPIRUN)
                    set(_ref_mpirun_l ${CPU_REFERENCE_MPIRUN})
                endif()
                add_test(NAME ${case}_cpu_vs_gpu_np4_${_laytag}
                    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/compare_monitor.py
                            --case-dir ${_case_dir} --fields --field-tol 1e-8 ${_extra} --p-grid ${_lay}
                            --mpirun ${_ref_mpirun_l} --exe-a ${CPU_REFERENCE_EXE} --np-a 4
                            --mpirun-b ${MPIEXEC_EXECUTABLE} --exe-b $<TARGET_FILE:dopamine> --np-b 4
                            --label-a cpu-${_laytag} --label-b gpu-${_laytag})
                if(TEST_LD_LIBRARY_PATH)
                    set_tests_properties(${case}_cpu_vs_gpu_np4_${_laytag} PROPERTIES
                        ENVIRONMENT "LD_LIBRARY_PATH=${TEST_LD_LIBRARY_PATH}:$ENV{LD_LIBRARY_PATH}")
                endif()
            endif()
        endforeach()

        # hot-start restart reproduces an uninterrupted run (all field blocks incl. scalars), at np=1 and np=N
        foreach(_np 1 ${TEST_PARITY_NPROCS})
            set(_rtol 1e-12)
            add_test(NAME ${case}_${_kind}_restart_np${_np}
                COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/restart_check.py
                        --case-dir ${_case_dir} --mpirun ${MPIEXEC_EXECUTABLE} --exe $<TARGET_FILE:dopamine> --np ${_np}
                        --tol ${_rtol})
            if(TEST_LD_LIBRARY_PATH)
                set_tests_properties(${case}_${_kind}_restart_np${_np} PROPERTIES
                    ENVIRONMENT "LD_LIBRARY_PATH=${TEST_LD_LIBRARY_PATH}:$ENV{LD_LIBRARY_PATH}")
            endif()
        endforeach()
        if(ENABLE_GPU AND CPU_REFERENCE_EXE)
            set(_ref_mpirun ${MPIEXEC_EXECUTABLE})
            if(CPU_REFERENCE_MPIRUN)
                set(_ref_mpirun ${CPU_REFERENCE_MPIRUN})
            endif()
            foreach(_np 1 ${TEST_PARITY_NPROCS})
                add_test(NAME ${case}_cpu_vs_gpu_np${_np}
                    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/compare_monitor.py
                            --case-dir ${_case_dir}
                            --mpirun ${_ref_mpirun} --exe-a ${CPU_REFERENCE_EXE} --np-a ${_np}
                            --mpirun-b ${MPIEXEC_EXECUTABLE} --exe-b $<TARGET_FILE:dopamine> --np-b ${_np}
                            --label-a cpu-np${_np} --label-b gpu-np${_np})
                if(TEST_LD_LIBRARY_PATH)
                    set_tests_properties(${case}_cpu_vs_gpu_np${_np} PROPERTIES
                        ENVIRONMENT "LD_LIBRARY_PATH=${TEST_LD_LIBRARY_PATH}:$ENV{LD_LIBRARY_PATH}")
                endif()
            endforeach()
        endif()
    endforeach()

    # Point-particle parity (tests/regression/compare_particles.py): final particle state by ID, np=1 vs other layouts and, in a
    # GPU build with CPU_REFERENCE_EXE, CPU vs GPU at the same layout. Deterministic cases only: Brownian motion and the SGS
    # Langevin closure draw from per-rank random streams, so their results are statistically, not bitwise, layout independent.
    set(_particle_cases part_tracer_small part_inertial_small part_bouss_small part_outflow_small part_ibm_small)
    set(_pkind cpu)
    if(ENABLE_GPU)
        set(_pkind gpu)
    endif()
    foreach(_pc ${_particle_cases})
        set(_pdir ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/${_pc})
        set(_pcmp ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/compare_particles.py
                  --case-dir ${_pdir} --mpirun ${MPIEXEC_EXECUTABLE} --exe-a $<TARGET_FILE:dopamine> --exe-b $<TARGET_FILE:dopamine>)
        foreach(_lay "2:" "3:" "4:2,2" "4:4,1")
            string(REPLACE ":" ";" _lparts "${_lay}")
            list(GET _lparts 0 _lnp)
            list(LENGTH _lparts _llen)
            set(_lgrid "")
            set(_ltag np${_lnp})
            if(_llen GREATER 1)
                list(GET _lparts 1 _lgrid_raw)
                if(_lgrid_raw)
                    set(_lgrid --p-grid ${_lgrid_raw})
                    string(REPLACE "," "x" _lgt ${_lgrid_raw})
                    set(_ltag np${_lnp}_${_lgt})
                endif()
            endif()
            add_test(NAME ${_pc}_${_pkind}_np1_vs_${_ltag}
                COMMAND ${_pcmp} --np-a 1 --np-b ${_lnp} ${_lgrid})
            if(TEST_LD_LIBRARY_PATH)
                set_tests_properties(${_pc}_${_pkind}_np1_vs_${_ltag} PROPERTIES
                    ENVIRONMENT "LD_LIBRARY_PATH=${TEST_LD_LIBRARY_PATH}:$ENV{LD_LIBRARY_PATH}")
            endif()
        endforeach()
        if(ENABLE_GPU AND CPU_REFERENCE_EXE)
            set(_pref_mpirun ${MPIEXEC_EXECUTABLE})
            if(CPU_REFERENCE_MPIRUN)
                set(_pref_mpirun ${CPU_REFERENCE_MPIRUN})
            endif()
            foreach(_pn 1 2)
                add_test(NAME ${_pc}_cpu_vs_gpu_np${_pn}
                    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/compare_particles.py
                            --case-dir ${_pdir} --mpirun ${_pref_mpirun} --exe-a ${CPU_REFERENCE_EXE} --np-a ${_pn}
                            --mpirun-b ${MPIEXEC_EXECUTABLE} --exe-b $<TARGET_FILE:dopamine> --np-b ${_pn}
                            --label-a cpu --label-b gpu)
                if(TEST_LD_LIBRARY_PATH)
                    set_tests_properties(${_pc}_cpu_vs_gpu_np${_pn} PROPERTIES
                        ENVIRONMENT "LD_LIBRARY_PATH=${TEST_LD_LIBRARY_PATH}:$ENV{LD_LIBRARY_PATH}")
                endif()
            endforeach()
        endif()
    endforeach()

    # IBM ghost-cell bookkeeping must not depend on the decomposition (tests/regression/ibm_ghost_counts.py; ibm_ghost_diff.py lists
    # the cells of a DOPAMINE_TRACE_DIR run that differ). The image-point stencils use an extended halo (ibm.f90, halo_pad.f90).
    foreach(_lay 2 1x4 2x2 4x1)
        add_test(NAME ibm_small_ghost_counts_np1_vs_${_lay}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/ibm_ghost_counts.py
                    --case-dir ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/ibm_small
                    --exe $<TARGET_FILE:dopamine> --mpirun ${MPIEXEC_EXECUTABLE} --layout-a 1 --layout-b ${_lay})
        if(TEST_LD_LIBRARY_PATH)
            set_tests_properties(ibm_small_ghost_counts_np1_vs_${_lay} PROPERTIES
                ENVIRONMENT "LD_LIBRARY_PATH=${TEST_LD_LIBRARY_PATH}:$ENV{LD_LIBRARY_PATH}")
        endif()
    endforeach()

    # Two-fluid / VOF cases (tests/regression/vof_check.py): volume conservation, boundedness, shape error against the stored
    # bound and decomposition independence of vof_diag.dat. The reversal cases run the prescribed LeVeque/Enright fields through the
    # full solver; vof_translate_small must be identical to round-off for every layout.

    # Single-run physics checks of the two-fluid solver (limits on vof_diag.dat, see tests/regression/vof_check.py): Laplace pressure
    # jump and parasitic currents of a static cylinder, capillary-wave and gravity-wave frequencies against linear theory.
    add_test(NAME vof_drop_small_laplace
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_check.py
                --case-dir ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_drop_small --exe $<TARGET_FILE:dopamine>
                --mpirun ${MPIEXEC_EXECUTABLE} --limit 29:3.96:4.04 --limit 16:0:3e-5)
    add_test(NAME vof_capwave_small_frequency
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_check.py
                --case-dir ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_capwave_small --exe $<TARGET_FILE:dopamine>
                --mpirun ${MPIEXEC_EXECUTABLE} --freq 15:2.22206:0.01)
    add_test(NAME vof_slosh_small_frequency
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_check.py
                --case-dir ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_slosh_small --exe $<TARGET_FILE:dopamine>
                --mpirun ${MPIEXEC_EXECUTABLE} --freq 15:7.8290:0.01 --rowmax 16:0.4)
    # Variable-density steep wave: x rank seams cut the interface obliquely (np1 vs the x-split layout, tight PCG tolerance so only round-off remains)
    add_test(NAME vof_break_small_np1_vs_np4_4x1
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_check.py
                --case-dir ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_break_small --exe $<TARGET_FILE:dopamine>
                --mpirun ${MPIEXEC_EXECUTABLE} --np-a 1 --np-b 4 --p-grid 4,1 --same-tol 1e-9
                --same-cols 4,5,6,7,11,12,13,16,26,27,29,32)
    # Snapshot pressure identical for np1 and a 2x2 layout, seam planes included
    add_test(NAME vof_pressure_seams_np1_vs_2x2
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_p_seam/check_pressure.py
                --exe $<TARGET_FILE:dopamine> --mpirun ${MPIEXEC_EXECUTABLE})
    # Immersed body with a free surface: still water must stay at rest with the liquid of the fluid region conserved (1 vs 2x2 ranks),
    # and the IBM load output of a submerged sphere must be the buoyancy (1 and 4 ranks)
    add_test(NAME vof_ibm_still_np1_vs_np4_2x2
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_check.py
                --case-dir ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_ibm_still --exe $<TARGET_FILE:dopamine>
                --mpirun ${MPIEXEC_EXECUTABLE} --np-a 1 --np-b 4 --p-grid 2,2 --max-dvol 1e-12
                --rowmax 16:1e-8 --rowmax 17:1e-8 --same-cols 4,5,6,7)
    add_test(NAME vof_ibm_slosh_np1_vs_np4_2x2
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_check.py
                --case-dir ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_ibm_slosh --exe $<TARGET_FILE:dopamine>
                --mpirun ${MPIEXEC_EXECUTABLE} --np-a 1 --np-b 4 --p-grid 2,2 --max-dvol 1e-12 --same-tol 1e-7
                --same-cols 4,5,6,7,11,12,13,14,16,17,18,28,32)
    foreach(_np 1 4)
        add_test(NAME vof_ibm_buoyancy_np${_np}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_ibm_buoy/check_buoyancy.py
                    --exe $<TARGET_FILE:dopamine> --mpirun ${MPIEXEC_EXECUTABLE} --np ${_np})
    endforeach()
    # Hot-start restart of the two-fluid solver reproduces an uninterrupted run (pressure history, liquid fraction, time-dependent wave
    # inlet and relaxation zones, LES): the static drop to round-off, the wave flume to the PCG tolerance
    foreach(_np 1 4)
        add_test(NAME vof_restart_small_np${_np}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/restart_check.py
                    --case-dir ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_restart_small --mpirun ${MPIEXEC_EXECUTABLE}
                    --exe $<TARGET_FILE:dopamine> --np ${_np} --tol 1e-12)
        add_test(NAME vof_restart_flume_np${_np}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/restart_check.py
                    --case-dir ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_restart_flume --mpirun ${MPIEXEC_EXECUTABLE}
                    --exe $<TARGET_FILE:dopamine> --np ${_np} --tol 1e-7)
    endforeach()
    # Unsupported feature combinations must abort at input
    add_test(NAME vof_input_aborts
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_input_checks/check_aborts.py
                --exe $<TARGET_FILE:dopamine> --mpirun ${MPIEXEC_EXECUTABLE})
    # Single-phase staircase IBM (ibm_method = 1) on a flat plate: stream decay by the log-law stress, round-off divergence, same on 1 and 4 ranks
    add_test(NAME ibm_stair_plate
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/ibm_stair_plate/check_stair_plate.py
                --exe $<TARGET_FILE:dopamine> --mpirun ${MPIEXEC_EXECUTABLE})
    # Staircase IBM: walls one cell thick close their faces like thicker ones (fences of 1, 2, 3 cells)
    add_test(NAME ibm_stair_fence
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/ibm_stair_fence/check_fence.py
                --exe $<TARGET_FILE:dopamine> --mpirun ${MPIEXEC_EXECUTABLE})
    # Staircase IBM wall stress: the per-object roughness of a closed face does not depend on the orientation of the stream
    add_test(NAME ibm_stair_ids
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/ibm_stair_ids/check_ids.py
                --exe $<TARGET_FILE:dopamine> --mpirun ${MPIEXEC_EXECUTABLE})
    # Single-phase flat-wall EQWM: a uniform stream between two modelled walls decays by the log-law wall stress, same on 1 and 4 ranks
    add_test(NAME wm_stream
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/wm_stream/check_wm_stream.py
                --exe $<TARGET_FILE:dopamine> --mpirun ${MPIEXEC_EXECUTABLE})
    add_test(NAME wm_stream_rough
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/wm_stream/check_wm_stream.py
                --exe $<TARGET_FILE:dopamine> --mpirun ${MPIEXEC_EXECUTABLE} --rough 1e-4 --tol 0.08)
    # Wall-modelled bed (smooth log law) under a uniform stream: wall stress rho u_tau^2, density independent, same on 1 and 4 ranks
    add_test(NAME vof_wm_stream
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_wm_stream/check_wm_stream.py
                --exe $<TARGET_FILE:dopamine> --mpirun ${MPIEXEC_EXECUTABLE})
    add_test(NAME vof_wm_stream_rough
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_wm_stream/check_wm_stream.py
                --exe $<TARGET_FILE:dopamine> --mpirun ${MPIEXEC_EXECUTABLE} --rough 1e-4 --tol 0.08)
    # IBM wall model on the staircase faces of a flat plate: stream decay by the log-law stress, density independent, same on 1 and 4 ranks
    add_test(NAME vof_ibm_wm_plate
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_ibm_wm_plate/check_wm_plate.py
                --exe $<TARGET_FILE:dopamine> --mpirun ${MPIEXEC_EXECUTABLE})
    # Recycled-precursor turbulence at the wave inlet (synthetic donor slice), masked above the free surface
    add_test(NAME vof_wave_turb_inlet
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_wave_turb_inlet/check_turb_inlet.py
                --exe $<TARGET_FILE:dopamine> --mpirun ${MPIEXEC_EXECUTABLE})
    # Mean log-law current plus a small wave through the wave inlet and the relaxation zones
    add_test(NAME vof_current_flume
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_current_flume/check_current.py
                --exe $<TARGET_FILE:dopamine> --mpirun ${MPIEXEC_EXECUTABLE})
    # Mirror symmetry of the wall padding of the momentum stencils (free-slip Taylor-Green vortex, ratio 1)
    add_test(NAME vof_tgv_mirror_symmetry
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_tgv_mirror/check_mirror.py
                --exe $<TARGET_FILE:dopamine> --mpirun ${MPIEXEC_EXECUTABLE})
    # LES (Vreman) + no-slip bed + wave inlet and relaxation zones: decomposition independence including the pressure checksums
    add_test(NAME vof_les_flume_small_np1_vs_np4_4x1
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_check.py
                --case-dir ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_les_flume_small --exe $<TARGET_FILE:dopamine>
                --mpirun ${MPIEXEC_EXECUTABLE} --np-a 1 --np-b 4 --p-grid 4,1 --max-dvol 1e-3 --same-tol 1e-9
                --same-cols 4,5,6,7,11,12,13,16,29,31,32)
    set(_vof_cases "vof_vortex_small|1.7e-4" "vof_deform_small|1.5e-2" "vof_translate_small|-1" "vof_ibm_small|-1" "vof_ibm_wm_small|-1")
    foreach(_vc ${_vof_cases})
        string(REPLACE "|" ";" _vcl ${_vc})
        list(GET _vcl 0 _vname)
        list(GET _vcl 1 _vl1)
        foreach(_lay "2,2" "4,1")
            string(REPLACE "," "x" _laytag ${_lay})
            add_test(NAME ${_vname}_np1_vs_np4_${_laytag}
                COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/vof_check.py
                        --case-dir ${CMAKE_CURRENT_SOURCE_DIR}/tests/regression/${_vname} --exe $<TARGET_FILE:dopamine>
                        --mpirun ${MPIEXEC_EXECUTABLE} --np-a 1 --np-b 4 --p-grid ${_lay} --max-l1 ${_vl1}
                        --same-cols 4,5,6,7,11,12,13,14,16,17,18,28,32)
        endforeach()
    endforeach()
endif()

#------------------------------------------------------------------------------
# Install
#------------------------------------------------------------------------------
install(TARGETS dopamine DESTINATION bin)
install(TARGETS dopamine-ESEM DESTINATION bin)
