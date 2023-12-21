list(APPEND src_firewx_gridspecs firewx_gridspecs.f)

add_executable(firewx_gridspecs.exe ${src_firewx_gridspecs})
target_link_libraries(firewx_gridspecs.exe PRIVATE w3emc::w3emc_4
                                           PRIVATE w3nco::w3nco_4)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(firewx_gridspecs.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS firewx_gridspecs.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS
  "-g -traceback -fp-model precise -assume noold_ldout_format")
   set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
   set(CMAKE_Fortran_FLAGS_DEBUG "-O0")
endif()
