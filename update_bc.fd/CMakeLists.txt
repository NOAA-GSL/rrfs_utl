list(APPEND src_update_bc
  update_bc.f90
  module_update_bc.f90
  module_io_fv3lam_bdy.f90
  module_io_fv3lam_bk.f90
  mod_fv3lam_wind.f90)

add_executable(update_bc.exe ${src_update_bc})
target_link_libraries(update_bc.exe PRIVATE ${PROJECT_NAME}::ncio
                                    PRIVATE gsi::gsi
                                    PRIVATE MPI::MPI_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(update_bc.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS update_bc.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
