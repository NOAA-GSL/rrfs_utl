list(APPEND src_rrfs_maxrh
pullandcalc.f90
param.incl
main.f90
kinds.f90
domain_dims.f
constants.f90
calcrh.f90)

add_executable(rrfs_maxrh.exe ${src_rrfs_maxrh})
target_link_libraries(rrfs_maxrh.exe PRIVATE w3nco::w3nco_4
                                    PRIVATE w3emc::w3emc_4
                                    PRIVATE bufr::bufr_4
                                    PRIVATE NetCDF::NetCDF_Fortran
                                    PRIVATE MPI::MPI_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(rrfs_maxrh.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS rrfs_maxrh.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
