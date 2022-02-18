file(GLOB src_all  module_read_NSSL_mosaic.f90 netCDFsub.f90 netCDFsub_geo.f90 process_NSSL_mosaic.f90 read_ncep_binary.f90 read_nssl_binary.f90 write_bufr_ref.f90 read_grib2_mod.f90)
file(GLOB src_all_enkf  DART.f90 netCDFsub.f90 netCDFsub_geo.f90 read_ncep_binary.f90 read_nssl_binary.f90 write_bufr_ref.f90 read_grib2_mod.f90 write_netcdf_ref.f90)
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

include_directories( ${GSIINC} )

set(EXE "process_NSSL_mosaic.exe")
add_executable(${EXE} ${src_all})
set_target_properties( ${EXE} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
set_target_properties( ${EXE} PROPERTIES LINK_FLAGS ${Fortran_LINK_FLAGS} )
target_link_libraries( ${EXE} NetCDF::NetCDF_Fortran ${LIBBASE}
                              ${wrflib}
                              ${GSILIB}
                              g2::g2_4
                              g2tmpl::g2tmpl
                              bufr::bufr_4
                              bacio::bacio_4
                              ${MPI_Fortran_LIBRARIES} )

install(TARGETS ${EXE}
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

#set(EXE "process_NSSL_mosaic_enkf.exe")
#add_executable(${EXE} ${src_all_enkf} process_NSSL_mosaic_to_dart_and_netcdf.f90)
#set_target_properties( ${EXE} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
#set_target_properties( ${EXE} PROPERTIES LINK_FLAGS ${Fortran_LINK_FLAGS} )
#target_link_libraries( ${EXE} NetCDF::NetCDF_Fortran ${wrflib} ${GSILIB} g2::g2_4 g2tmpl::g2tmpl bufr::bufr_4 bacio::bacio_4 ${MPI_Fortran_LIBRARIES} )
#install(TARGETS ${EXE}
#        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
