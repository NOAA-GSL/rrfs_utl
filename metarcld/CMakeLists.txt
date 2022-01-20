cmake_minimum_required(VERSION 2.8)

file(GLOB src_all ${CMAKE_CURRENT_SOURCE_DIR}/*90)
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

include_directories( ${GSIINC} ${NETCDFINC} )

set(EXE "process_metarcld.exe")
add_executable(${EXE} ${src_all} )
target_link_libraries( ${EXE} ${LIBBASE} ${wrflib} ${wrfframe} ${NETCDFLIB} ${LIBZ} ${GSILIB} ${LIBW3EMC} ${LIBW3NCO} ${LIBBACIO} ${LIBBUFR} ${MPILIB} )
install(TARGETS ${EXE}
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
