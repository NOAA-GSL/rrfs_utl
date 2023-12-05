
list(APPEND SRCS
  angle_bias.f90
  create_ctl_angle.f90
)

add_executable(radmon_angle.x ${SRCS})

target_link_libraries(radmon_angle.x PRIVATE ${PROJECT_NAME}::radmonlib)
target_link_libraries(radmon_angle.x PRIVATE NetCDF::NetCDF_Fortran)
target_link_libraries(radmon_angle.x PRIVATE ncdiag::ncdiag_serial)
target_link_libraries(radmon_angle.x PRIVATE w3emc::w3emc_4)

# Install executable targets
install(TARGETS radmon_angle.x RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

