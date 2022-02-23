file(GLOB src_all pesg.f90 pfun.f90 pietc.f90 pietc_s.f90 pkind.f90 pmat.f90 pmat2.f90 pmat4.f90 pmat5.f90 psym2.f90)
set_source_files_properties(${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

set(module_dir "${CMAKE_CURRENT_BINARY_DIR}/include")
add_library(pesglib STATIC ${src_all})
add_library(${PROJECT_NAME}::pesglib ALIAS pesglib)
set_target_properties(pesglib PROPERTIES Fortran_MODULE_DIRECTORY "${module_dir}")
target_include_directories(pesglib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                                          $<INSTALL_INTERFACE:include>)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})

install(
  TARGETS pesglib
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
