cmake_minimum_required(VERSION 3.15)

project(RRFS_Utl
	VERSION 1.0.0
	LANGUAGES C Fortran)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_DIRECTORY_LABELS ${PROJECT_NAME})

include(compiler_flags)
include(GNUInstallDirs)

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel|Production)$")
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo" "Production")
endif()

if(NOT CMAKE_C_COMPILER_ID MATCHES "^(GNU|Intel|Clang|AppleClang)$")
  message(WARNING "${CMAKE_C_COMPILER_ID} is not supported.")
endif()

if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Intel)$")
  message(WARNING "${CMAKE_Fortran_COMPILER_ID} is not supported.")
endif()

# User options
option(OPENMP "Enable OpenMP Threading" ON)

# Echo user options
message(STATUS "OPENMP ................. ${OPENMP}")

# Dependencies
find_package(MPI REQUIRED)
find_package(NetCDF REQUIRED Fortran)
if(OPENMP)
  find_package(OpenMP REQUIRED)
endif()

# NCEPLibs dependencies
find_package(bacio REQUIRED)
find_package(g2 REQUIRED)
find_package(g2tmpl REQUIRED)
find_package(sp REQUIRED)
find_package(ip REQUIRED)
find_package(w3emc REQUIRED)
find_package(w3nco REQUIRED)
find_package(bufr REQUIRED)
find_package(wrf_io REQUIRED)

# Other dependencies
find_package(ncdiag REQUIRED)
find_package(gsi REQUIRED)

add_subdirectory(baselib/regional_esg_grid.fd)
add_subdirectory(baselib)
add_subdirectory(mosaic)
add_subdirectory(lightning)
add_subdirectory(ref2tten)
add_subdirectory(cloudanalysis)
add_subdirectory(larccld)
add_subdirectory(metarcld)
add_subdirectory(update_sst)
add_subdirectory(use_raphrrr_sfc)
add_subdirectory(update_gvf)
add_subdirectory(update_snow)
add_subdirectory(update_ice)
add_subdirectory(lakesurgery)
add_subdirectory(gen_cs)
add_subdirectory(adjust_soiltq)
add_subdirectory(update_bc)
add_subdirectory(ens_mean_recenter)
add_subdirectory(FV3_ensmean_recenter)
