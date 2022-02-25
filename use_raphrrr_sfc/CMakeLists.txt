file(GLOB src_all netCDFsub_geo.f90 module_surface.f90 use_raphrrr_sfc.f90)

set(EXE "use_raphrrr_sfc.exe")
add_executable(${EXE} ${src_all})
target_link_libraries( ${EXE} PUBLIC ${PROJECT_NAME}::ncio
                              PUBLIC ${PROJECT_NAME}::map_util
                              PUBLIC ${PROJECT_NAME}::kinds
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
