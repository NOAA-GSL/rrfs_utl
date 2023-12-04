list(APPEND src_cloudanalysis
  cloudanalysis_fv3driver.f90
  module_gsi_rfv3io_tten.f90
  namelist_mod.f90
  pbl_height.f90
  module_fv3lam_stats.f90
  module_mpi_arrange.f90
  get_fv3sar_bk_parallel_mod.f90
  general_sub2grid_simple_mod.f90
  rapidrefresh_cldsurf_mod.f90)

add_executable(fv3lam_nonvarcldana.exe ${src_cloudanalysis})
target_link_libraries(fv3lam_nonvarcldana.exe PRIVATE ${PROJECT_NAME}::ncio
                                              PRIVATE ${PROJECT_NAME}::gsdcloud_arw
                                              PRIVATE gsi::gsi
                                              PRIVATE NetCDF::NetCDF_Fortran
                                              PRIVATE MPI::MPI_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(fv3lam_nonvarcldana.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS fv3lam_nonvarcldana.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_subdirectory(NonVarCldLib)
