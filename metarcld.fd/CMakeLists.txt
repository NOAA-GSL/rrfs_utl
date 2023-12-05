list(APPEND src_metarcld
  cld_parm_array_module.f90
  process_metar_cloud.f90
  read_prepbufr_metarcld.f90
  reorg_metar_cloud_regular.f90)

add_executable(process_metarcld.exe ${src_metarcld})
target_link_libraries(process_metarcld.exe PRIVATE ${PROJECT_NAME}::ncio
                                           PRIVATE ${PROJECT_NAME}::esggrid_util
                                           PRIVATE ${PROJECT_NAME}::pesglib
                                           PRIVATE NetCDF::NetCDF_Fortran
                                           PRIVATE gsi::gsi
                                           PRIVATE w3emc::w3emc_4
                                           PRIVATE w3nco::w3nco_4
                                           PRIVATE bufr::bufr_4
                                           PRIVATE MPI::MPI_Fortran)

install(
  TARGETS process_metarcld.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
