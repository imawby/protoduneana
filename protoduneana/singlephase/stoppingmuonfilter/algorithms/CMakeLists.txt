include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

cet_enable_asserts()

# define ROOTLIB
set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} )

art_make( BASENAME_ONLY
  LIBRARY_NAME protodune_singlephase_stopmufilter_algorithms
  
  LIB_LIBRARIES 
  lardataobj_RecoBase
  larana_OpticalDetector_OpHitFinder
  larsim_MCSTReco
  lardataobj_MCBase
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  larsim_PhotonPropagation_PhotonVisibilityService_service
  larana_OpticalDetector_OpDigiProperties_service
  lardata_Utilities
  nusimdata_SimulationBase
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_BASIC}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  art_Persistency_Common
  art_Utilities
  canvas
  ${MF_MESSAGELOGGER}
  ${MF_UTILITIES}
  ${FHICLCPP}
  cetlib cetlib_except
  ${ROOT_BASIC_LIB_LIST}
  
  MODULE_LIBRARIES 
  larana_OpticalDetector_OpDigiProperties_service
  larana_OpticalDetector_OpHitFinder
  larana_OpticalDetector
  larsim_PhotonPropagation_PhotonVisibilityService_service
  larsim_MCCheater_BackTrackerService_service
  larsim_MCCheater_ParticleInventoryService_service
  larevt_Filters
  lardata_RecoObjects
  lardataobj_RecoBase
  lardataobj_MCBase
  larreco_RecoAlg
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  larsim_Simulation 
  lardataobj_Simulation
  lardata_Utilities
  lardataobj_RawData
  larreco_Calorimetry
  lardataobj_AnalysisBase
  nusimdata_SimulationBase
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_BASIC}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}

  ${ART_PERSISTENCY_PROVENANCE}
  art_Persistency_Common
  art_Utilities
  canvas
  ${MF_MESSAGELOGGER}
  ${MF_UTILITIES}
  ${FHICLCPP}
  cetlib cetlib_except
  ${ROOT_BASIC_LIB_LIST}
  
  SERVICE_LIBRARIES
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_BASIC}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  art_Persistency_Common
  art_Utilities
  canvas
  ${MF_MESSAGELOGGER}
  ${MF_UTILITIES}
  ${FHICLCPP}
  cetlib cetlib_except
  ${CLHEP}
  ${ROOT_BASIC_LIB_LIST}
  )

install_headers()
install_source()
