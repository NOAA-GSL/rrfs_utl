list(APPEND src_gen_be_ensmean
  gen_be_ensmean.f90 ncio_ensmean_recenter.f90)

add_executable(gen_be_ensmean.exe ${src_gen_be_ensmean})
target_link_libraries(gen_be_ensmean.exe PRIVATE ${PROJECT_NAME}::ncio
                                    PRIVATE gsi::gsi
                                    PRIVATE MPI::MPI_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(gen_be_ensmean.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS gen_be_ensmean.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
