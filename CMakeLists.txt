cmake_minimum_required(VERSION 2.8)

set(ENV{CC} "icc")
set(ENV{CXX} "icpc")
set(ENV{FC} "mpif90 -f90=ifort")

project(data_process)
enable_language (Fortran)
enable_language (C)

site_name(HOSTNAME)
message("The hostname is  ${HOSTNAME}" )
string(REGEX MATCH "cheyenne" HOST-Cheyenne ${HOSTNAME}  )
string(REGEX MATCH "tfe[0-9]" HOST-Theia ${HOSTNAME} )
string(REGEX MATCH "hfe[0-9]" HOST-Hera ${HOSTNAME} )
if(EXISTS /jetmon)
  set(HOST-Jet "True" )
endif()

### Do the following to load required compiling modules on Jet
#     source /home/rtrr/PARM_EXEC/modulefiles/modulefile.jet.GSI_UPP_WRF 
###

if(HOST-Cheyenne)
  ###SET(GSIBUILD  "/glade/work/geguo/MPMC/b_release_v3.7/ci06.intel_17.0.1.mpt_2.15f")
  SET(GSIBUILD  "/glade/p/ral/jntp/gge/precompiled/GSI/build")
  set(COMUPP "/glade/p/ral/jntp/gge/precompiled/GSD_UPP/comupp")
  set(wrfframe "/glade/p/ral/jntp/gge/precompiled/WRFV3.9/frame/pack_utils.o")
  set(WPSBUILD "/glade/p/ral/jntp/gge/precompiled/WPSV3.9.0.01")
  ### MPI_Fortran_LIBRARIES
  SET(MPILIB  "/usr/lib64/libpthread.so;/glade/u/apps/ch/opt/mpt/2.19/lib/libmpi.so")
  ### BLAS_iomp5_LIBRARY
  SET(LIBIOMP5  "/glade/u/apps/opt/intel/2017u1/compilers_and_libraries/linux/lib/intel64/libiomp5.so")
  SET(LIBPNG    "/usr/lib64/libpng.so")
elseif(HOST-Jet)
  SET(GSIBUILD  "/lfs4/BMC/wrfruc/gge/precompiled/GSI")
  set(COMUPP "/lfs4/BMC/wrfruc/gge/precompiled/comupp")
  set(wrfframe "/lfs4/BMC/wrfruc/gge/precompiled/WRFV3.9.1.1/frame/pack_utils.o")
  set(WPSBUILD "/lfs4/BMC/wrfruc/gge/precompiled/WPSV3.9")
  SET(LIBIOMP5  "-qopenmp")
  SET(LIBPNG    "/usr/lib64/libpng12.so.0")
elseif(HOST-Hera)
  SET(GSIBUILD  "/scratch1/BMC/comgsi/precompiled/ProdGSI")
  set(COMUPP "/scratch1/BMC/comgsi/precompiled/comupp")
  set(wrfframe "/scratch1/BMC/comgsi/precompiled/WRFV3.9/frame/pack_utils.o")
  set(WPSBUILD "/scratch1/BMC/comgsi/precompiled/WPSV3.9")
  SET(MPILIB  "/apps/intel/compilers_and_libraries_2018/linux/mpi/intel64/lib/libmpifort.so;/apps/intel/compilers_and_libraries_2018/linux/mpi/intel64/lib/release_mt/libmpi.so;/apps/intel/compilers_and_libraries_2018/linux/mpi/intel64/lib/libmpigi.a;/usr/lib64/libdl.so;/usr/lib64/librt.so;/usr/lib64/libpthread.so")
  SET(LIBIOMP5  "-qopenmp")
  SET(LIBPNG    "/usr/lib64/libpng.so")
else()
  message(FATAL_ERROR "
       Unknow host ${HOSTNAME}

        ")
endif()

#----------------------------Normally no changes needed after this line ------------------------
SET(NETCDF  $ENV{NETCDF})

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_INCLUDE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/include")
set(Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include")
set(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include")

SET(GSIINC  "${GSIBUILD}/include")
SET(NETCDFINC  "${NETCDF}/include")
SET(G2INC  "${COMUPP}/include")

SET(GSILIB  "${GSIBUILD}/lib/libgsilib_shrd.a")
SET(wrflib  "${GSIBUILD}/lib/libWRFLIB.a")
SET(NETCDFLIB  "${NETCDF}/lib/libnetcdf.so;${NETCDF}/lib/libnetcdff.so")

SET(LIBZ      "/usr/lib64/libz.so")
SET(LIBJASPER "/usr/lib64/libjasper.so")
SET(GSILIBS "${GSIBUILD}/lib")
SET(LIBBUFR "${GSILIBS}/libbufr_v.a")
SET(LIBW3EMC "${GSILIBS}/libw3emc_v.a")
SET(LIBW3NCO "${GSILIBS}/libw3nco_v.a")
SET(LIBBACIO "${GSILIBS}/libbacio_v.a")

SET(LIBW3EMC_UPP "${COMUPP}/lib/libw3emc.a")
SET(LIBW3NCO_UPP "${COMUPP}/lib/libw3nco.a")
SET(LIBBACIO_UPP "${COMUPP}/lib/libbacio.a")
SET(LIBG2     "${COMUPP}/lib/libg2.a")
SET(LIBG2TMPL "${COMUPP}/lib/libg2tmpl.a")
SET(LIBSP     "${COMUPP}/lib/libsp.a")
SET(LIBIP     "${COMUPP}/lib/libip.a")

set(LIBG2_WPS "${WPSBUILD}/ungrib/src/ngl/g2/libg2_4.a")
SET(CONSTANTS "${GSIBUILD}/src/gsi/CMakeFiles/gsilib_shrd.dir/constants.f90.o")

SET(Fortran_FLAGS "-O3 -fp-model precise -assume byterecl -fpe0 -DLINUX -g -traceback")
SET(Fortran_FLAGS_bigEndian "-convert big_endian")
SET(C_FLAGS "-D_UNDERSCORE -DLINUX -Dfunder -DFortranByte=char -DFortranInt=int -DFortranLlong='long long'")

add_subdirectory(mosaic)
add_subdirectory(lightning/bufr)
add_subdirectory(lightning/netcdf)
add_subdirectory(ref2tten)
add_subdirectory(cloudanalysis)
add_subdirectory(larccld)
add_subdirectory(metarcld)
