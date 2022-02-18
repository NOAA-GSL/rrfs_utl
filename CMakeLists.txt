cmake_minimum_required(VERSION 3.15)

project(RRFS_Utl
	VERSION 1.0.0
	LANGUAGES C Fortran)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(CMAKE_DIRECTORY_LABELS ${PROJECT_NAME})

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

# Dependencies
find_package(MPI REQUIRED)
find_package(NetCDF REQUIRED Fortran)

# NCEPLibs dependencies
find_package(bacio REQUIRED)
find_package(g2 REQUIRED)
find_package(g2tmpl REQUIRED)
find_package(sp REQUIRED)
find_package(ip REQUIRED)
find_package(w3emc REQUIRED)
find_package(w3nco REQUIRED)
find_package(bufr REQUIRED)

# Other dependencies
find_package(wrf_io REQUIRED)

SET(GSIBUILD  "${PROJECT_BINARY_DIR}/../../../build/src/GSI/src/GSI-build")


#----------------------------Normally no changes needed after this line ------------------------

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_INCLUDE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/include")
set(Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include")
set(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include")

SET(GSIINC  "${GSIBUILD}/include")

SET(GSILIB  "${GSIBUILD}/lib/libgsilib_shrd.a")
SET(wrflib  "${GSIBUILD}/lib/libgsilib_wrf.a")

SET(GSILIBS "${GSIBUILD}/lib")
SET(CONSTANTS "${GSIBUILD}/src/gsi/CMakeFiles/gsilib_wrf.dir/constants.f90.o")

SET(Fortran_FLAGS "-O3 -fp-model precise -assume byterecl -fpe0 -DLINUX -g -traceback")
SET(Fortran_LINK_FLAGS "-qopenmp")
SET(Fortran_FLAGS_bigEndian "-convert big_endian")
SET(C_FLAGS "-D_UNDERSCORE -DLINUX -Dfunder -DFortranByte=char -DFortranInt=int -DFortranLlong='long long'")

SET(LIBBASE "-L${PROJECT_BINARY_DIR}/lib -lmap_util -lesggrid_util -lncio  -lpesglib  -ltimelib")

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
add_subdirectory(update_snowice)
add_subdirectory(cycle_ice)
add_subdirectory(lakesurgery)
add_subdirectory(gen_cs)
