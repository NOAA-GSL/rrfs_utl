cmake_minimum_required(VERSION 2.8)

file(GLOB src_all ${CMAKE_CURRENT_SOURCE_DIR}/*.f90)
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

include_directories( ${NETCDFINC} )

option(BUILD_GSDCLOUD_ARW "Build the GSD cloud analysis " ON)
set(GSDCLOUD_Fortran_FLAGS "-O3 -convert big_endian -g -traceback" CACHE INTERNAL "")
if(BUILD_GSDCLOUD_ARW)
   set(GSDCLOUD_DIR "${CMAKE_SOURCE_DIR}/cloudanalysis/NoVarCldLib")
   set(gsdcloud "gsdcloud_arw")
   add_subdirectory( ./NoVarCldLib )
   set(GSDCLOUD_LIBRARY ${gsdcloud} )
endif(BUILD_GSDCLOUD_ARW)

set(EXE "fv3sar_novarcldana.exe")
add_executable(${EXE} ${src_all})
set_target_properties( ${EXE} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
target_link_libraries( ${EXE} ${wrflib} ${NETCDFLIB} ${MPILIB} ${GSDCLOUD_LIBRARY} )
install(TARGETS ${EXE}
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
