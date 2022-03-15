list(APPEND src_adjust_soiltq
  adjust_soiltq.f90
  unfill_fv3lam_grid2.f90
  gsl_adjust_soiltq_mod.f90
  module_bkio_fv3lam.f90
  genqsat.f90
  constants.f90)

set(EXE "adjust_soiltq.exe")
add_executable(adjust_soiltq.exe ${src_adjust_soiltq})
target_link_libraries(adjust_soiltq.exe PRIVATE ${PROJECT_NAME}::ncio
                                        PRIVATE gsi::gsi
                                        PRIVATE MPI::MPI_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(adjust_soiltq.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS adjust_soiltq.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
