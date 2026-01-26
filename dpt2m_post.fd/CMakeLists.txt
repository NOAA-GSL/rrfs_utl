set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl")
set(CMAKE_EXE_LINKER_FLAGS "-mkl")

#set(CMAKE_Fortran_FLAGS "-O3 -fp-model precise -fpe0 -mkl" CACHE STRING "" FORCE)
#set(CMAKE_EXE_LINKER_FLAGS "-mkl" CACHE STRING "" FORCE)

message(STATUS "CHECK: Fortran Flags: ${CMAKE_Fortran_FLAGS}")
message(STATUS "CHECK: Release Flags: ${CMAKE_Fortran_FLAGS_RELEASE}")
message(STATUS "CHECK: Production Flags: ${CMAKE_Fortran_FLAGS_PRODUCTION}")
message(STATUS "CHECK: Linker Flags: ${CMAKE_EXE_LINKER_FLAGS}")

list(APPEND src_dpt2m_post
  dpt2m_post.f90)

include_directories(
    $ENV{WGRIB2_INC}
    $ENV{NETCDF_INCLUDES}
    $ENV{IP_INC4}
)

add_executable(dpt2m_post.exe ${src_dpt2m_post})

target_link_libraries(dpt2m_post.exe
    PRIVATE
    $ENV{WGRIB2_LIBDIR}/libwgrib2_ftn_api.a
    $ENV{WGRIB2_LIBDIR}/libwgrib2.a
    $ENV{IP_LIB4}
    $ENV{G2C_LIB}
    $ENV{SP_LIB4}
    $ENV{W3EMC_LIB4}
    -L$ENV{NETCDF_LIBRARIES}
    netcdff
    netcdf
    $ENV{LIBAEC_LIB}
    jasper
    png
    z
)

install(
  TARGETS dpt2m_post.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
