list(APPEND src_gen_ensmean_recenter
  gen_ensmean_recenter.f90 ncio_ensmean_recenter.f90)

add_executable(gen_ensmean_recenter.exe ${src_gen_ensmean_recenter})
target_link_libraries(gen_ensmean_recenter.exe PRIVATE ${PROJECT_NAME}::ncio
                                    PRIVATE gsi::gsi
                                    PRIVATE MPI::MPI_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(gen_ensmean_recenter.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS gen_ensmean_recenter.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
