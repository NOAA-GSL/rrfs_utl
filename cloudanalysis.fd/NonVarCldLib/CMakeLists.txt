list(APPEND src_gsdcloud_arw
  adaslib.f90
  ARPS_cldLib.f90
  BackgroundCld.f90
  BckgrndCC.f90
  build_missing_REFcone.f90
  cloudCover_NESDIS.f90
  cloudCover_radar.f90
  cloudCover_Surface.f90
  cloudLayers.f90
  cloudLWC.f90
  cloud_saturation.f90
  cloudType.f90
  constants.f90
  convert_lghtn2ref.f90
  get_sfm_1d_gnl.f90
  hydro_mxr_thompson.f90
  map_ctp.f90
  map_ctp_lar.f90
  mthermo.f90
  pbl_height.f90
  pcp_mxr_ARPSlib.f90
  PrecipMxr_radar.f90
  PrecipType.f90
  radar_ref2tten.f90
  read_Lightning_cld.f90
  read_nasalarc_cld.f90
  read_NESDIS.f90
  read_radar_ref.f90
  read_Surface.f90
  smooth.f90
  TempAdjust.f90
  vinterp_radar_ref.f90)

add_library(gsdcloud_arw STATIC ${src_gsdcloud_arw})
add_library(${PROJECT_NAME}::gsdcloud_arw ALIAS gsdcloud_arw)
target_link_libraries(gsdcloud_arw PUBLIC ${PROJECT_NAME}::ncio
                                   PUBLIC gsi::gsi)
set_target_properties(gsdcloud_arw PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
target_include_directories(gsdcloud_arw PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                                               $<INSTALL_INTERFACE:include>)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})

install(
  TARGETS gsdcloud_arw
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
