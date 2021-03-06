#include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)
link_libraries( -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB})

art_make(
          LIBRARY_NAME GeneratorLevelAnalysis
          LIB_LIBRARIES
                        larpandora_LArPandoraInterface
                        lardataobj_RecoBase
                        nusimdata_SimulationBase
                        lardataobj_MCBase
			lardataobj_Simulation
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        art_Framework_Core
                        art_Framework_Services_Registry
                        art_Framework_Principal
			art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_COMMON}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${ART_UTILITIES}
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_BASIC_LIB_LIST}
                        ${SIMULATIONBASE}
                        larsim_MCCheater_BackTracker_service

          MODULE_LIBRARIES
                        larpandora_LArPandoraInterface
			larsim_MCCheater_BackTracker_service
                        lardataobj_RecoBase
                        nusimdata_SimulationBase
                        lardataobj_MCBase
			lardataobj_Simulation
                        larcore_Geometry
                        larcore_Geometry_Geometry_service
                        art_Framework_Core
                        art_Framework_Services_Registry
                        art_Framework_Principal
                        art_Persistency_Common canvas_Persistency_Common
			art_Persistency_Provenance canvas_Persistency_Provenance
			art_Utilities canvas_Utilities
			${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_COMMON}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${ART_UTILITIES}
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_BASIC_LIB_LIST}
                        ${SIMULATIONBASE}
                        GeneratorLevelAnalysis
)

add_subdirectory(job)

install_headers()
install_source()
install_fhicl()
