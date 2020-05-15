include_directories(
  ${CMAKE_SOURCE_DIR}
  ${DelphesExternals_INCLUDE_DIR} 
)



# build all executables and put them into bin/
add_executable(DelphesEDM4HEP DelphesEDM4HEP.cpp Delphes.cc)
target_link_libraries(DelphesEDM4HEP Delphes)
install(TARGETS DelphesEDM4HEP DESTINATION bin)
