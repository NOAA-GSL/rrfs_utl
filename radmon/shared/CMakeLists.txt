list(APPEND SRCS
  kinds_radmon.F90
  read_diag.f90
  sparsearr.f90
)

set(module_dir "${CMAKE_CURRENT_BINARY_DIR}/include")
add_library(radmonlib STATIC ${SRCS})
add_library(${PROJECT_NAME}::radmonlib ALIAS radmonlib)

target_compile_definitions(radmonlib PRIVATE "_REAL8_")
target_link_libraries(radmonlib PRIVATE NetCDF::NetCDF_Fortran)
target_link_libraries(radmonlib PRIVATE ncdiag::ncdiag_serial)
set_target_properties(radmonlib PROPERTIES Fortran_MODULE_DIRECTORY "${module_dir}")
target_include_directories(radmonlib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                                          $<INSTALL_INTERFACE:include>)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})

install(
  TARGETS radmonlib
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
