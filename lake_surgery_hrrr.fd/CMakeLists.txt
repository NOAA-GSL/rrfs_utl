list(APPEND src_lake_surgery_hrrr
  lake_surgery_hrrr.f90 )

add_executable(lake_surgery_hrrr.exe ${src_lake_surgery_hrrr})
target_link_libraries(lake_surgery_hrrr.exe PRIVATE ${PROJECT_NAME}::map_util
                                          PRIVATE gsi::gsi
                                          PRIVATE ${PROJECT_NAME}::ncio
                                          PRIVATE NetCDF::NetCDF_Fortran
                                          PRIVATE MPI::MPI_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(lake_surgery_hrrr.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS lake_surgery_hrrr.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

list(APPEND fix_soil_vegetation_pct
  fix_soil_vegetation_pct.f90 )

add_executable(fix_soil_vegetation_pct.exe ${fix_soil_vegetation_pct})
target_link_libraries(fix_soil_vegetation_pct.exe PRIVATE ${PROJECT_NAME}::map_util
                                          PRIVATE gsi::gsi
                                          PRIVATE ${PROJECT_NAME}::ncio
                                          PRIVATE NetCDF::NetCDF_Fortran
                                          PRIVATE MPI::MPI_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(fix_soil_vegetation_pct.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS fix_soil_vegetation_pct.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
