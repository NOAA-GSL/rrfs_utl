list(APPEND src_dpt2m_post
  dpt2m_post.f90)

add_executable(dpt2m_post.exe ${src_dpt2m_post})
target_link_libraries(dpt2m_post.exe PRIVATE wgrib2::wgrib2)

install(
  TARGETS dpt2m_post.exe
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

