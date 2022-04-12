list(APPEND src_lightning
  Check_Lightning_QC.f90
  Check_NLDN.f90
  netCDFsub_lightning.f90
  process_Lightning.f90
  read_lightning_bufr.f90
  write_bufr_lght.f90)

add_executable(process_Lightning.exe ${src_lightning})
target_link_libraries(process_Lightning.exe PRIVATE ${PROJECT_NAME}::ncio
                                            PRIVATE ${PROJECT_NAME}::esggrid_util
                                            PRIVATE ${PROJECT_NAME}::pesglib
                                            PRIVATE NetCDF::NetCDF_Fortran
                                            PRIVATE bufr::bufr_4
                                            PRIVATE w3emc::w3emc_4
                                            PRIVATE w3nco::w3nco_4
                                            PRIVATE gsi::gsi
                                            PRIVATE MPI::MPI_Fortran
                                            PRIVATE MPI::MPI_C)

install(
  TARGETS process_Lightning.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
