list(APPEND src_rrfs_maxt
writemaxbg.f90
param.incl
maxmin.f
maxmin_ak.f
domain_dims.f
calcmaxbg.f90
makefile)

add_executable(rrfs_maxt.exe ${src_rrfs_maxt})
target_link_libraries(rrfs_maxt.exe PRIVATE w3nco::w3nco_4
                                    PRIVATE w3emc::w3emc_4
                                    PRIVATE bufr::bufr_4
                                    PRIVATE NetCDF::NetCDF_Fortran
                                    PRIVATE MPI::MPI_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(rrfs_maxt.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS rrfs_maxt.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
