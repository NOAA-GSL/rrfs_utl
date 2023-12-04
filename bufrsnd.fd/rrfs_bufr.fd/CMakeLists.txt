list(APPEND src_rrfs_bufr
WRFBUFR.f
VRBLS2D_mod.f
VRBLS3D_mod.f
MASKS_mod.f
SOIL_mod.f
PROF_FV3SAR_NET.f
iw3jdn.f 
w3fs26.f 
w3movdat.f 
PARA_RANGE.f90
w3reddat.f)

add_executable(rrfs_bufr.exe ${src_rrfs_bufr})
target_link_libraries(rrfs_bufr.exe PRIVATE w3nco::w3nco_4
                                    PRIVATE bufr::bufr_4
                                    PRIVATE NetCDF::NetCDF_Fortran
                                    PRIVATE MPI::MPI_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(rrfs_bufr.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS rrfs_bufr.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS
  "-g -traceback -fp-model precise -assume noold_ldout_format")
   set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
   set(CMAKE_Fortran_FLAGS_DEBUG "-O0")

#elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
#   set(CMAKE_Fortran_FLAGS
#   "-g -fbacktrace -ffree-form -ffree-line-length-none -fconvert=big-endian")
#   set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
#   set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -ggdb -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check")
#  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
#     set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
#  endif()

endif()
