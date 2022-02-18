file(GLOB src_all ${CMAKE_CURRENT_SOURCE_DIR}/*90)
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
file(GLOB src_c ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
set_source_files_properties(${src_c} PROPERTIES COMPILE_FLAGS ${C_FLAGS})

include_directories( ${GSIINC} )

set(EXE "process_larccld.exe")
add_executable(${EXE} ${src_c} ${src_all} )
#set_target_properties( ${EXE} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
target_link_libraries( ${EXE} ${PROJECT_NAME}::ncio
                              ${PROJECT_NAME}::esggrid_util
                              ${PROJECT_NAME}::pesglib
                              NetCDF::NetCDF_Fortran
                              ${wrflib}
                              ${wrfframe}
                              ${GSILIB}
                              bufr::bufr_4
                              MPI::MPI_Fortran
                              MPI::MPI_C )
install(TARGETS ${EXE}
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
