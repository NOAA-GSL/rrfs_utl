file(GLOB src_all)

include_directories( ${GSIINC} ${NETCDFINC} )

set(EXE "gen_cs.exe")
add_executable(${EXE} ${src_all} gen_cs.f90  )
target_link_libraries( ${EXE} ${LIBBASE} ${wrflibio} NetCDF::NetCDF_Fortran ${GSILIB} ${LIBW3NCO} ${LIBZ} ${LIBPNG} ${LIBJASPER} ${LIBBACIO} ${MPI_Fortran_LIBRARIES} )
