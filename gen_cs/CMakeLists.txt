file(GLOB src_all)

include_directories( ${GSIINC} )

set(EXE "gen_cs.exe")
add_executable(${EXE} ${src_all} gen_cs.f90  )
target_link_libraries( ${EXE} ${PROJECT_NAME}::ncio
                              NetCDF::NetCDF_Fortran
                              ${GSILIB}
                              MPI::MPI_Fortran )

install(
  TARGETS ${EXE}
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
