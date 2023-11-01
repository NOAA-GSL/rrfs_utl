
list(APPEND SRCS
  bcoef.f90
  create_ctl_bcoef.f90
)

add_executable(radmon_bcoef.x ${SRCS})

#include_directories( /scratch1/BMC/wrfruc/mhu/rrfs/v0.7.1/ufs-srweather-app/include )
target_link_libraries(radmon_bcoef.x PRIVATE ${PROJECT_NAME}::radmonlib)
target_link_libraries(radmon_bcoef.x PRIVATE NetCDF::NetCDF_Fortran)
target_link_libraries(radmon_bcoef.x PRIVATE ncdiag::ncdiag_serial)
target_link_libraries(radmon_bcoef.x PRIVATE w3emc::w3emc_4)

# Install executable targets
install(TARGETS radmon_bcoef.x RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

