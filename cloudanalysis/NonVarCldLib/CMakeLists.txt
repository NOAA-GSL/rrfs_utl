file(GLOB src_all ${CMAKE_CURRENT_SOURCE_DIR}/*.f90)
set(GSDCLOUD_Fortran_FLAGS "-O3 -g -traceback" CACHE INTERNAL "")
set_source_files_properties( ${src_all} COMPILE_FLAGS ${GSDCLOUD_Fortran_FLAGS})

set(module_dir "${CMAKE_CURRENT_BINARY_DIR}/include")
add_library( gsdcloud_arw STATIC ${src_all} )
add_library(${PROJECT_NAME}::gsdcloud_arw ALIAS gsdcloud_arw)
target_link_libraries( gsdcloud_arw PUBLIC ${PROJECT_NAME}::kinds)
set_target_properties( gsdcloud_arw PROPERTIES Fortran_MODULE_DIRECTORY "${module_dir}")
target_include_directories( gsdcloud_arw PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                                                $<INSTALL_INTERFACE:include>)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})

install(
  TARGETS gsdcloud_arw
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
