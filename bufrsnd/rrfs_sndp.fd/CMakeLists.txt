list(APPEND src_rrfs_sndp
  SNDPST.f
  BFRHDR.f
  BFRIZE.f
  CALHEL.f
  CALVIS.f 
  CALWXT_BOURG.f
  CALWXT_DOMINANT.f
  CALWXT_EXPLICIT.f 
  CALWXT.f
  CALWXT_RAMER.f
  CALWXT_REVISED.f 
  SMIX.f
  SNDPST_FV3S.f
  VAP.f
  WETBLB.f)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-line-length-120")
add_executable(rrfs_sndp.exe ${src_rrfs_sndp})
target_link_libraries(rrfs_sndp.exe PRIVATE bufr::bufr_4)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(rrfs_sndp.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS rrfs_sndp.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

#if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
#  set(CMAKE_Fortran_FLAGS
#  "-g -traceback -fp-model precise -convert big_endian -assume noold_ldout_format")
#   set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
#   set(CMAKE_Fortran_FLAGS_DEBUG "-O0")
#elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
#   set(CMAKE_Fortran_FLAGS
#   "-g -fbacktrace -ffree-form -ffree-line-length-none -fconvert=big-endian")
#   set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
#   set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -ggdb -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check")
#  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
#     set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
#  endif()
#endif()
