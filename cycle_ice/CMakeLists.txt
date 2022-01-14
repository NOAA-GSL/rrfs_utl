cmake_minimum_required(VERSION 2.8)

file(GLOB src_all)
#set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS_bigEndian} )

include_directories( ${GSIINC} ${NETCDFINC} ${G2INC} )

set(EXE "cycle_ice.exe")
add_executable(${EXE} ${src_all} cycle_ice.f90  )
target_link_libraries( ${EXE} ${LIBBASE} ${wrflibio} ${NETCDFLIB} ${GSILIB} ${LIBW3NCO} ${LIBZ} ${LIBPNG} ${LIBJASPER} ${MPILIB} ${LIBBACIO} ${LIBIOMP5})
