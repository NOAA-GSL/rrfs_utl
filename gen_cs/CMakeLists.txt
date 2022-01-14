cmake_minimum_required(VERSION 2.8)

file(GLOB src_all)

include_directories( ${GSIINC} ${NETCDFINC} )

set(EXE "gen_cs.exe")
add_executable(${EXE} ${src_all} gen_cs.f90  )
target_link_libraries( ${EXE} ${LIBBASE} ${wrflibio} ${NETCDFLIB} ${GSILIB} ${LIBW3NCO} ${LIBZ} ${LIBPNG} ${LIBJASPER} ${MPILIB} ${LIBBACIO} ${LIBIOMP5})
