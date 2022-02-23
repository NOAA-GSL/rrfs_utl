file(GLOB src_all  netCDFsub_geo.f90 surface_lakesurge_netcdf_mass.f90 )
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

include_directories( ${GSIINC} )

set(EXE "lakesurgery.exe")
add_executable(${EXE} ${src_all})
set_target_properties( ${EXE} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
target_link_libraries( ${EXE} PUBLIC wrf_io::wrf_io
                              PUBLIC NetCDF::NetCDF_Fortran
                              PUBLIC ${GSILIB}
                              PUBLIC MPI::MPI_Fortran )
if(OpenMP_Fortran_FOUND)
  target_link_libraries( ${EXE} PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS ${EXE}
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
