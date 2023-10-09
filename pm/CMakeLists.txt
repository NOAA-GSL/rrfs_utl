list(APPEND src_pm
  process_pm.f90)

add_executable(process_pm.exe ${src_pm})
target_link_libraries(process_pm.exe PRIVATE bufr::bufr_4
                                          PRIVATE MPI::MPI_Fortran)

install(
  TARGETS process_pm.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
