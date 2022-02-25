file(GLOB src_all ${CMAKE_CURRENT_SOURCE_DIR}/*.f90)

set(EXE "fv3lam_nonvarcldana.exe")
add_executable(${EXE} ${src_all})
target_link_libraries( ${EXE} ${PROJECT_NAME}::ncio
                              ${PROJECT_NAME}::gsdcloud_arw
                              NetCDF::NetCDF_Fortran
                              ${wrflib}
                              MPI::MPI_Fortran )

install(
  TARGETS ${EXE}
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})


add_subdirectory( ./NonVarCldLib )
