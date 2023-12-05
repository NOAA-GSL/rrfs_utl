
list(APPEND SRCS
  avgsdv.f90
  bcor.f90
  create_ctl_bcor.f90
)

add_executable(radmon_bcor.x ${SRCS})

target_link_libraries(radmon_bcor.x PRIVATE ${PROJECT_NAME}::radmonlib)
target_link_libraries(radmon_bcor.x PRIVATE NetCDF::NetCDF_Fortran)
target_link_libraries(radmon_bcor.x PRIVATE ncdiag::ncdiag_serial)
target_link_libraries(radmon_bcor.x PRIVATE w3emc::w3emc_4)

# Install executable targets
install(TARGETS radmon_bcor.x RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

