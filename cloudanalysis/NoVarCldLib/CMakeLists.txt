cmake_minimum_required(VERSION 2.6)

#option(BUILD_GSDCLOUD_ARW "Build the GSD cloud analysis " ON)
#set(GSDCLOUD_DIR "${CMAKE_SOURCE_DIR}/cloudanalysis/NoVarCldLib")
#set(gsdcloud gsdcloud_arw)
#set(GSDCLOUD_Fortran_FLAGS "-O3 -convert big_endian" CACHE INTERNAL "")

if(BUILD_GSDCLOUD_ARW)
  file(GLOB GSDCLOUD_SRC ${GSDCLOUD_DIR}/*.f90)
  set_source_files_properties( ${GSDCLOUD_SRC} COMPILE_FLAGS ${GSDCLOUD_Fortran_FLAGS})
  add_library( ${gsdcloud} STATIC ${GSDCLOUD_SRC} )
  set_target_properties( ${gsdcloud} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_INCLUDE_OUTPUT_DIRECTORY} )
  message("HEY GSD Cloud analysis lib is called ${gsdcloud}")
endif()
