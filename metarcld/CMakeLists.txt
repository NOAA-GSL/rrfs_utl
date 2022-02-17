file(GLOB src_all ${CMAKE_CURRENT_SOURCE_DIR}/*90)
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

include_directories( ${GSIINC} ${NETCDFINC} )

set(EXE "process_metarcld.exe")
add_executable(${EXE} ${src_all} )
target_link_libraries( ${EXE} ${LIBBASE} NetCDF::NetCDF_Fortran ${wrflib} ${wrfframe} ${LIBZ} ${GSILIB} ${LIBW3EMC} ${LIBW3NCO} ${LIBBACIO} ${LIBBUFR} ${MPI_Fortran_LIBRARIES} )
install(TARGETS ${EXE}
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
