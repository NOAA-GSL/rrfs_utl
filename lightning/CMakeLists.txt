cmake_minimum_required(VERSION 2.8)
file(GLOB src_all ${CMAKE_CURRENT_SOURCE_DIR}/*90)
set_source_files_properties( ${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )
file(GLOB src_c ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
set_source_files_properties(${src_c} PROPERTIES COMPILE_FLAGS ${C_FLAGS})

include_directories( ${GSIINC} ${NETCDFINC} )

ADD_EXECUTABLE(process_Lightning.exe process_Lightning.f90 ${src_c} ${src_all} )
TARGET_LINK_LIBRARIES(process_Lightning.exe  ${LIBBASE} ${NETCDFLIB} ${LIBZ} ${LIBBUFR} ${LIBW3EMC} ${LIBW3NCO} ${LIBBACIO} ${MPILIB} ${CONSTANTS})
install(TARGETS process_Lightning.exe
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
