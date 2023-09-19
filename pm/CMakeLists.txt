list(APPEND src_pm
  process_pm.f90)

add_executable(process_pm.exe ${src_pm})
target_link_libraries(process_pm.exe PRIVATE ${PROJECT_NAME}::ncio
                                          PRIVATE ${PROJECT_NAME}::esggrid_util
                                          PRIVATE ${PROJECT_NAME}::pesglib
                                          PRIVATE NetCDF::NetCDF_Fortran
                                          PRIVATE gsi::gsi
                                          PRIVATE bufr::bufr_4
                                          PRIVATE MPI::MPI_Fortran
                                          PRIVATE MPI::MPI_C)

install(
  TARGETS process_pm.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
