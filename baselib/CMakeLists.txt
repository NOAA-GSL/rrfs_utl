set(module_dir "${CMAKE_CURRENT_BINARY_DIR}/include")

add_library(ncio STATIC module_ncio.f90)
add_library(${PROJECT_NAME}::ncio ALIAS ncio)
set_target_properties(ncio PROPERTIES Fortran_MODULE_DIRECTORY "${module_dir}")
target_include_directories(ncio PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                                       $<INSTALL_INTERFACE:include>)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})
list(APPEND BASE_LIB_Targets ncio)

add_library (map_util STATIC module_map_utils.f90)
add_library(${PROJECT_NAME}::map_util ALIAS map_util)
target_link_libraries(map_util PUBLIC gsi::gsi)
set_target_properties(map_util PROPERTIES Fortran_MODULE_DIRECTORY "${module_dir}")
target_include_directories(map_util PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                                           $<INSTALL_INTERFACE:include>)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})
list(APPEND BASE_LIB_Targets map_util)

add_library (esggrid_util STATIC module_esggrid_util.f90)
add_library(${PROJECT_NAME}::esggrid_util ALIAS esggrid_util)
target_link_libraries( esggrid_util PUBLIC ${PROJECT_NAME}::pesglib)
set_target_properties(esggrid_util PROPERTIES Fortran_MODULE_DIRECTORY "${module_dir}")
target_include_directories(esggrid_util PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                                               $<INSTALL_INTERFACE:include>)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})
list(APPEND BASE_LIB_Targets esggrid_util)

install(
  TARGETS ${BASE_LIB_Targets}
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
