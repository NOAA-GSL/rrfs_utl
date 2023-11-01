list(APPEND SRCS
  bad_chan.f90
  bad_penalty.f90
  create_ctl_time.f90
  low_count.f90
  time.f90
  valid.f90
)

add_executable(radmon_time.x ${SRCS})

target_link_libraries(radmon_time.x PRIVATE ${PROJECT_NAME}::radmonlib)
target_link_libraries(radmon_time.x PRIVATE NetCDF::NetCDF_Fortran)
target_link_libraries(radmon_time.x PRIVATE ncdiag::ncdiag_serial)
target_link_libraries(radmon_time.x PRIVATE w3emc::w3emc_4)

# Install executable targets
install(TARGETS radmon_time.x RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

