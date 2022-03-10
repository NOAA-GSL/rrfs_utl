list(APPEND src_lakesurgery
  netCDFsub_geo.f90
  surface_lakesurge_netcdf_mass.f90 )

add_executable(lakesurgery.exe ${src_lakesurgery})
target_link_libraries(lakesurgery.exe PRIVATE wrf_io::wrf_io
                                      PRIVATE NetCDF::NetCDF_Fortran
                                      PRIVATE gsi::gsi
                                      PRIVATE MPI::MPI_Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(lakesurgery.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS lakesurgery.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
