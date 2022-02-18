file(GLOB src_all)
#set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS_bigEndian} )

include_directories( ${GSIINC} )

set(EXE "cycle_ice.exe")
add_executable(${EXE} ${src_all} cycle_ice.f90  )
target_link_libraries( ${EXE} ${LIBBASE}
                              NetCDF::NetCDF_Fortran
                              ${GSILIB}
                              ${MPI_Fortran_LIBRARIES} )
