# This is the CMake file for the regional_esf_grid utility in
# UFS_UTILS.
#
# George Gayno
set(fortran_src
    pesg.f90
    pfun.f90
    pietc.f90
    pietc_s.f90
    pkind.f90
    pmat.f90
    pmat2.f90
    pmat4.f90
    pmat5.f90
    psym2.f90
    )

file(GLOB src_all pesg.f90 pfun.f90 pietc.f90 pietc_s.f90 pkind.f90 pmat.f90 pmat2.f90 pmat4.f90 pmat5.f90 psym2.f90)
set_source_files_properties(${src_all} PROPERTIES COMPILE_FLAGS ${Fortran_FLAGS} )

add_library(pesglib STATIC ${src_all})

#add_executable(regional_esg_grid ${fortran_src} regional_esg_grid.f90)
#target_link_libraries(regional_esg_grid ${NETCDFLIB} ${LIBZ} ${LIBPNG} ${LIBJASPER})
