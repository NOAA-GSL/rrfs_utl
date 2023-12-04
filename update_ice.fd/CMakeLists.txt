list(APPEND src_update_ice
  update_ice.f90)

add_executable(update_ice.exe ${src_update_ice})
target_link_libraries(update_ice.exe PRIVATE ${PROJECT_NAME}::ncio
                                     PRIVATE NetCDF::NetCDF_Fortran
                                     PRIVATE MPI::MPI_Fortran)

install(
  TARGETS update_ice.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
