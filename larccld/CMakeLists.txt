list(APPEND src_larccld
  netCDFsub.f90
  process_NASALaRC_cloud.f90
  read_NASALaRC_cloud.f90
  write_bufr_NASALaRC.f90)

add_executable(process_larccld.exe ${src_larccld})
target_link_libraries(process_larccld.exe PRIVATE ${PROJECT_NAME}::ncio
                                          PRIVATE ${PROJECT_NAME}::esggrid_util
                                          PRIVATE ${PROJECT_NAME}::pesglib
                                          PRIVATE NetCDF::NetCDF_Fortran
                                          PRIVATE gsi::gsi
                                          PRIVATE bufr::bufr_4
                                          PRIVATE MPI::MPI_Fortran
                                          PRIVATE MPI::MPI_C)

install(
  TARGETS process_larccld.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
