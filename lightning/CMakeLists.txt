file(GLOB src_all ${CMAKE_CURRENT_SOURCE_DIR}/*90)
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
file(GLOB src_c ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
set_source_files_properties(${src_c} PROPERTIES COMPILE_FLAGS ${C_FLAGS})

include_directories( ${GSIINC} )

set(EXE process_Lightning.exe)
ADD_EXECUTABLE(${EXE} process_Lightning.f90 ${src_c} ${src_all} )
TARGET_LINK_LIBRARIES(${EXE} ${PROJECT_NAME}::ncio
                             ${PROJECT_NAME}::esggrid_util
                             ${PROJECT_NAME}::pesglib
                             NetCDF::NetCDF_Fortran
                             bufr::bufr_4
                             w3emc::w3emc_4
                             w3nco::w3nco_4
                             ${CONSTANTS}
                             MPI::MPI_Fortran
                             MPI::MPI_C )

install(
  TARGETS ${EXE}
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
