list(APPEND src_use_raphrrr_sfc
  netCDFsub_geo.f90
  module_surface.f90
  use_raphrrr_sfc.f90)

add_executable(use_raphrrr_sfc.exe ${src_use_raphrrr_sfc})
target_link_libraries(use_raphrrr_sfc.exe PRIVATE ${PROJECT_NAME}::map_util
                                          PRIVATE gsi::gsi
                                          PRIVATE ${PROJECT_NAME}::ncio
                                          PRIVATE NetCDF::NetCDF_Fortran
                                          PRIVATE MPI::MPI_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(use_raphrrr_sfc.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS use_raphrrr_sfc.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
