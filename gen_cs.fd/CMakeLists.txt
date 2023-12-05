list(APPEND src_gen_cs
  gen_cs.f90)

add_executable(gen_cs.exe ${src_gen_cs})
target_link_libraries(gen_cs.exe PRIVATE ${PROJECT_NAME}::ncio
                                 PRIVATE NetCDF::NetCDF_Fortran
                                 PRIVATE MPI::MPI_Fortran)

install(
  TARGETS gen_cs.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
