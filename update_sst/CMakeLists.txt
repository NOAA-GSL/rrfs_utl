cmake_minimum_required(VERSION 2.8)

file(GLOB src_all grib2_read_mod.f90 netCDFsub_geo.f90  process_SST.f90  sstGlobal2RR.f90  update_SST_netcdf_mass.f90)
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

include_directories( ${GSIINC} ${G2INC} ${NETCDFINC} )

set(EXE "process_updatesst")
add_executable(${EXE} ${src_all})
set_target_properties( ${EXE} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
set_target_properties( ${EXE} PROPERTIES LINK_FLAGS ${Fortran_LINK_FLAGS} )
target_link_libraries( ${EXE} ${wrflibio} ${NETCDFLIB} ${GSILIB} ${LIBG2} ${LIBG2TMPL} ${LIBW3NCO} ${LIBZ} ${LIBPNG} ${LIBJASPER} ${MPILIB} ${LIBBACIO} ${LIBIOMP5})

