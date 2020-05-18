include_directories(
  ${CMAKE_SOURCE_DIR}
  ${DelphesExternals_INCLUDE_DIR} 
)

find_package(podio)
find_package(EDM4HEP)



# build all executables and put them into bin/
add_executable(DelphesEDM4HEP DelphesEDM4HEP.cpp Delphes.cc)
target_link_libraries(DelphesEDM4HEP Delphes)
install(TARGETS DelphesEDM4HEP DESTINATION bin)

add_executable(DelphesROOT_EDM4HEP DelphesROOT_EDM4HEP.cpp Delphes.cc)
target_link_libraries(DelphesROOT_EDM4HEP Delphes EDM4HEP::edm4hep podio::podioRootIO)
install(TARGETS DelphesROOT_EDM4HEP DESTINATION bin)
