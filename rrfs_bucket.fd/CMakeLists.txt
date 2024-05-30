list(APPEND src_bucket
  bucket_diff.f)

add_executable(rrfs_bucket.exe ${src_bucket})
target_link_libraries(rrfs_bucket.exe PRIVATE ${PROJECT_NAME}::ncio
                                            PRIVATE w3emc::w3emc_4
                                            PRIVATE w3nco::w3nco_4
                                            PRIVATE g2tmpl::g2tmpl
                                            PRIVATE g2::g2_4
                                            PRIVATE bacio::bacio_4
                                            PRIVATE MPI::MPI_Fortran
                                            PRIVATE MPI::MPI_C)

install(
  TARGETS rrfs_bucket.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
