list(APPEND src_stnmlist stnmlist.f)

add_executable(rrfs_stnmlist.exe ${src_stnmlist})
target_link_libraries(rrfs_stnmlist.exe PRIVATE w3nco::w3nco_4
                                    PRIVATE bufr::bufr_4)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(rrfs_stnmlist.exe PRIVATE OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS rrfs_stnmlist.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
