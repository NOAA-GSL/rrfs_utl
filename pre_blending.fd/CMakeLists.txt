list(APPEND src_pre_blending
  pre_blending.f90
  module_fv3lam_stats.f90
  module_mpi_arrange.f90
  general_sub2grid_simple_mod.f90
  chgres_winds.f90
  remap_dwinds.f90
  remap_scalar.f90)

find_library(LIBIOMP5_PATH NAMES libiomp5.so)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
   set(CMAKE_Fortran_FLAGS "-fp-model precise -L${LIBIOMP5_PATH} -liomp5")
   set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -L${LIBIOMP5_PATH} -liomp5")
   set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -traceback -fp-model precise -L${LIBIOMP5_PATH} -liomp5")
endif()

add_executable(fv3lam_pre_blending.exe ${src_pre_blending})
target_link_libraries(fv3lam_pre_blending.exe PRIVATE ${PROJECT_NAME}::ncio
	                                      PRIVATE gsi::gsi
                                              PRIVATE NetCDF::NetCDF_Fortran
                                              PRIVATE MPI::MPI_Fortran)

install(
  TARGETS fv3lam_pre_blending.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
