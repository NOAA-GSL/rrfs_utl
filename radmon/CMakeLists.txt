set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fp-model precise -assume byterecl -convert big_endian")

add_subdirectory(shared)
add_subdirectory(verf_radang.fd)
add_subdirectory(verf_radbcoef.fd)
add_subdirectory(verf_radbcor.fd)
add_subdirectory(verf_radtime.fd)
