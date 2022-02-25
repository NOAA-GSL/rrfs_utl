file(GLOB src_all ${CMAKE_CURRENT_SOURCE_DIR}/*.f90)

include_directories( ${GSIINC} )

set(EXE "ref2tten.exe")
add_executable(${EXE} ${src_all})
target_link_libraries( ${EXE} NetCDF::NetCDF_Fortran
                              ${wrflib}
                              ${GSILIB}
                              MPI::MPI_Fortran )

install(
  TARGETS ${EXE}
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

