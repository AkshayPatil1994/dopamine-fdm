# FindFFTW.cmake
# Locates the FFTW3 library (double precision, serial).
#
# Imported targets created:
#   FFTW::Double        - serial double-precision library + headers
#
# Cache variables:
#   FFTW_ROOT           - root of FFTW installation (set via -DFFTW_ROOT=...)
#   FFTW_INCLUDE_DIRS   - include directory
#   FFTW_LIBRARIES      - combined list: fftw3 m
#
# Supports components: DOUBLE

include(FindPackageHandleStandardArgs)

# Allow the user to point us at an installation
set(_fftw_root_hints
    ${FFTW_ROOT}
    $ENV{FFTW_ROOT}
    $ENV{FFTW_HOME}
    $ENV{EBROOTFFTW}          # EasyBuild
    $ENV{TACC_FFTW3_DIR}      # TACC
)

# --- header (the solver includes the Fortran interface fftw3.f03, not only fftw3.h) ---
find_path(FFTW_INCLUDE_DIR
    NAMES fftw3.f03
    HINTS ${_fftw_root_hints} $ENV{FFTW_INC}     # FFTW_INC: Cray
    PATH_SUFFIXES include include/fftw
)

# --- serial library ---
find_library(FFTW_DOUBLE_LIB
    NAMES fftw3
    HINTS ${_fftw_root_hints} $ENV{FFTW_DIR}     # FFTW_DIR: Cray (the lib directory)
    PATH_SUFFIXES lib lib64
)

find_library(FFTW_MATH_LIB NAMES m)

set(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})
set(FFTW_LIBRARIES ${FFTW_DOUBLE_LIB})
if(FFTW_MATH_LIB)
    list(APPEND FFTW_LIBRARIES ${FFTW_MATH_LIB})
endif()

# Set per-component found variables required by HANDLE_COMPONENTS
if(FFTW_DOUBLE_LIB)
    set(FFTW_DOUBLE_FOUND TRUE)
    set(FFTW_DOUBLE_LIB_FOUND TRUE)   # decomp2d checks this name
endif()

find_package_handle_standard_args(FFTW
    REQUIRED_VARS FFTW_INCLUDE_DIR FFTW_DOUBLE_LIB
    HANDLE_COMPONENTS
)
mark_as_advanced(FFTW_INCLUDE_DIR FFTW_DOUBLE_LIB)

# Create imported target so consumers can do target_link_libraries(... FFTW::Double)
if(FFTW_FOUND AND NOT TARGET FFTW::Double)
    add_library(FFTW::Double UNKNOWN IMPORTED)
    set_target_properties(FFTW::Double PROPERTIES
        IMPORTED_LOCATION "${FFTW_DOUBLE_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIR}"
    )
endif()
