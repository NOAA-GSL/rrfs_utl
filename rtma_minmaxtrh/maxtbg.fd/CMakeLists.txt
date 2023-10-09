list(APPEND src_rrfs_mint
writeminbg.f90
param.incl
maxmin.f
maxmin_ak.f
domain_dims.f
calcminbg.f90
makefile)

add_executable(rrfs_mint.exe ${src_rrfs_mint})
target_link_libraries(rrfs_mint.exe PRIVATE w3nco::w3nco_4
                                    PRIVATE w3emc::w3emc_4
                                    PRIVATE bufr::bufr_4
                                    PRIVATE NetCDF::NetCDF_Fortran
                                    PRIVATE MPI::MPI_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(rrfs_mint.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS rrfs_mint.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
