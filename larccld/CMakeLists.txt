file(GLOB src_all ${CMAKE_CURRENT_SOURCE_DIR}/*90)

include_directories( ${GSIINC} )

set(EXE "process_larccld.exe")
add_executable(${EXE} ${src_all} )
target_link_libraries( ${EXE} ${PROJECT_NAME}::ncio
                              ${PROJECT_NAME}::esggrid_util
                              ${PROJECT_NAME}::pesglib
                              NetCDF::NetCDF_Fortran
                              gsi::gsi
                              bufr::bufr_4
                              MPI::MPI_Fortran
                              MPI::MPI_C )

install(
  TARGETS ${EXE}
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
