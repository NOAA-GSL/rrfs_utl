cmake_minimum_required(VERSION 2.8)

file(GLOB src_all  netCDFsub_geo.f90 grib2_read_snow_mod.f90 process_NESDIS_imssnow.f90  update_SNOWICE_netcdf_mass.f90 module_wrf_error.f90 )
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

include_directories( ${GSIINC} ${IPINC} ${G2INC} ${NETCDFINC} )

set(EXE "process_NESDIS_imssnow")
add_executable(${EXE} ${src_all})
set_target_properties( ${EXE} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
set_target_properties( ${EXE} PROPERTIES LINK_FLAGS ${Fortran_LINK_FLAGS} )
target_link_libraries( ${EXE} ${wrflibio} ${LIBIP} ${LIBSP} ${NETCDFLIB} ${GSILIB} ${LIBG2} ${LIBG2TMPL} ${LIBW3NCO} ${LIBZ} ${LIBPNG} ${LIBJASPER} ${MPILIB} ${LIBBACIO} ${LIBIOMP5})

