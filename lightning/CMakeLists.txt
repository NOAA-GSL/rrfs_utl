file(GLOB src_all ${CMAKE_CURRENT_SOURCE_DIR}/*90)
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
file(GLOB src_c ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
set_source_files_properties(${src_c} PROPERTIES COMPILE_FLAGS ${C_FLAGS})

include_directories( ${GSIINC} )

set(EXE process_Lightning.exe)
ADD_EXECUTABLE(${EXE} process_Lightning.f90 ${src_c} ${src_all} )
TARGET_LINK_LIBRARIES(${EXE} PUBLIC ${PROJECT_NAME}::ncio
                             PUBLIC ${PROJECT_NAME}::esggrid_util
                             PUBLIC ${PROJECT_NAME}::pesglib
                             PUBLIC NetCDF::NetCDF_Fortran
                             PUBLIC bufr::bufr_4
                             PUBLIC w3emc::w3emc_4
                             PUBLIC w3nco::w3nco_4
                             PUBLIC ${GSILIB}
                             PUBLIC MPI::MPI_Fortran
                             PUBLIC MPI::MPI_C )

install(
  TARGETS ${EXE}
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
