##### Library ######
add_library(burstFit STATIC src/BurstModel.cxx src/NegativeStat.cxx)

target_link_libraries(
  burstFit
  PUBLIC optimizers
  PRIVATE evtbin st_app st_stream st_facilities st_graph tip CLHEP::CLHEP
)

target_include_directories(
  burstFit
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<INSTALL_INTERFACE:>
)

###### Executables ######
add_executable(gtburstfit src/gtburstfit/gtburstfit.cxx)
target_link_libraries(gtburstfit PRIVATE burstFit evtbin st_graph st_app st_graph st_stream tip)

###### Tests ######
add_executable(test_burstFit src/test/test_burstFit.cxx)
target_link_libraries(test_burstFit PRIVATE burstFit evtbin st_graph st_app st_graph st_stream tip)

###############################################################
# Installation
###############################################################
install(DIRECTORY burstFit DESTINATION ${FERMI_INSTALL_INCLUDEDIR})
install(DIRECTORY data/ DESTINATION ${FERMI_INSTALL_REFDATADIR}/burstFit)
install(DIRECTORY pfiles/ DESTINATION ${FERMI_INSTALL_PFILESDIR})

install(
  TARGETS burstFit gtburstfit test_burstFit
  EXPORT fermiTargets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
