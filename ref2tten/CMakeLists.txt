file(GLOB src_all ${CMAKE_CURRENT_SOURCE_DIR}/*.f90)
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

include_directories( ${GSIINC} ${NETCDFINC} )

set(EXE "ref2tten.exe")
add_executable(${EXE} ${src_all})
set_target_properties( ${EXE} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
target_link_libraries( ${EXE} NetCDF::NetCDF_Fortran ${wrflib} ${LIBZ} ${GSILIB} ${MPI_Fortran_LIBRARIES} )
install(TARGETS ${EXE}
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

