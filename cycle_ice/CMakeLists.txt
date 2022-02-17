file(GLOB src_all)
#set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS_bigEndian} )

include_directories( ${GSIINC} ${NETCDFINC} ${G2INC} )

set(EXE "cycle_ice.exe")
add_executable(${EXE} ${src_all} cycle_ice.f90  )
target_link_libraries( ${EXE} ${LIBBASE} ${wrflibio} NetCDF::NetCDF_Fortran ${GSILIB} ${LIBW3NCO} ${LIBZ} ${LIBPNG} ${LIBJASPER} ${LIBBACIO} ${MPI_Fortran_LIBRARIES} )
