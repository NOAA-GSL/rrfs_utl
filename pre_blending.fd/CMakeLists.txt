list(APPEND src_pre_blending
  pre_blending.f90
  module_fv3lam_stats.f90
  module_mpi_arrange.f90
  chgres_winds.f90
  remap_dwinds.f90
  remap_scalar.f90)

add_executable(fv3lam_pre_blending.exe ${src_pre_blending})
target_link_libraries(fv3lam_pre_blending.exe PRIVATE ${PROJECT_NAME}::ncio
                                              PRIVATE NetCDF::NetCDF_Fortran
                                              PRIVATE MPI::MPI_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(fv3lam_pre_blending.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS fv3lam_pre_blending.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
