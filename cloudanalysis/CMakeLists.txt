file(GLOB src_all ${CMAKE_CURRENT_SOURCE_DIR}/*.f90)
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

option(BUILD_GSDCLOUD_ARW "Build the GSD cloud analysis " ON)
#set(GSDCLOUD_Fortran_FLAGS "-O3 -convert big_endian -g -traceback" CACHE INTERNAL "")
set(GSDCLOUD_Fortran_FLAGS "-O3 -g -traceback" CACHE INTERNAL "")
if(BUILD_GSDCLOUD_ARW)
   set(GSDCLOUD_DIR "${CMAKE_SOURCE_DIR}/cloudanalysis/NonVarCldLib")
   set(gsdcloud "gsdcloud_arw")
   add_subdirectory( ./NonVarCldLib )
   set(GSDCLOUD_LIBRARY ${gsdcloud} )
endif(BUILD_GSDCLOUD_ARW)

set(EXE "fv3lam_nonvarcldana.exe")
add_executable(${EXE} ${src_all})
set_target_properties( ${EXE} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
target_link_libraries( ${EXE} ${LIBBASE}
                              NetCDF::NetCDF_Fortran
                              ${wrflib}
                              ${GSDCLOUD_LIBRARY}
                              MPI::MPI_Fortran )
install(TARGETS ${EXE}
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
