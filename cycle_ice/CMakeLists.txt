list(APPEND src_cycle_ice
  cycle_ice.f90)

add_executable(cycle_ice.exe ${src_cycle_ice})
target_link_libraries(cycle_ice.exe PRIVATE ${PROJECT_NAME}::ncio
                                    PRIVATE NetCDF::NetCDF_Fortran
                                    PRIVATE MPI::MPI_Fortran)

install(
  TARGETS cycle_ice.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
