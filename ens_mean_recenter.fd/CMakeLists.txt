list(APPEND src_ens_mean_recenter
  ens_mean_recenter.f90
  module_fv3lam_stats.f90
  module_mpi_arrange.f90)

add_executable(ens_mean_recenter_P2DIO.exe ${src_ens_mean_recenter})
target_link_libraries(ens_mean_recenter_P2DIO.exe PRIVATE ${PROJECT_NAME}::ncio
                                    PRIVATE gsi::gsi
                                    PRIVATE MPI::MPI_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(ens_mean_recenter_P2DIO.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS ens_mean_recenter_P2DIO.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
