### CMakeLists.txt ---

project(Core
VERSION
  ${${CMAKE_PROJECT_NAME}_VERSION}
LANGUAGES
  CXX)

## #################################################################
## Core Target Definition (add_library, export header, version, alias)
## #################################################################

add_library(${PROJECT_NAME} SHARED)

generate_export_header(${PROJECT_NAME}
  EXPORT_FILE_NAME "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}${PROJECT_NAME}Export.h")

set_target_properties(${PROJECT_NAME}
  PROPERTIES
    OUTPUT_NAME ${CMAKE_PROJECT_NAME}${PROJECT_NAME}
    VERSION   ${${PROJECT_NAME}_VERSION}
    SOVERSION ${${PROJECT_NAME}_VERSION_MAJOR})

add_library(${CMAKE_PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

## #################################################################
## Sources, Headers, Includes, Linking, Install
## #################################################################

target_sources(${PROJECT_NAME}
  PRIVATE
    Aligner.cpp
    Bestscenarioscounter.cpp
    CDR3SeqData.cpp
    Counter.cpp
    Coverageerrcounter.cpp
    Deletion.cpp
    Dinuclmarkov.cpp
    Errorrate.cpp
    Errorscounter.cpp
    EventUtils.cpp
    ExtractFeatures.cpp
    FastGenerator.cpp
    FastSampling.cpp
    GenModel.cpp
    Genechoice.cpp
    HypermutationfullNmererrorrate.cpp
    Hypermutationglobalerrorrate.cpp
    Insertion.cpp
    IntStr.cpp
    Model_Parms.cpp
    Model_marginals.cpp
    Pgencounter.cpp
    Rec_Event.cpp
    Singleerrorrate.cpp
    Utils.cpp

  PUBLIC FILE_SET HEADERS
    BASE_DIRS
      ${CMAKE_CURRENT_SOURCE_DIR}
      ${CMAKE_CURRENT_BINARY_DIR}
    FILES
      AccumulationContext.h
      Aligner.h
      AlignmentPreset.h
      SeqAn2Aligner.h
      Bestscenarioscounter.h
      CDR3SeqData.h
      Counter.h
      Coverageerrcounter.h
      Deletion.h
      Dinuclmarkov.h
      Errorrate.h
      Errorscounter.h
      EventUtils.h
      ExplorationContext.h
      ExtractFeatures.h
      FastGenerator.h
      FastSampling.h
      GenModel.h
      Genechoice.h
      HypermutationfullNmererrorrate.h
      Hypermutationglobalerrorrate.h
      Insertion.h
      IntStr.h
      ModelContext.h
      Model_Parms.h
      Model_marginals.h
      Pgencounter.h
      QuerySequenceContext.h
      Rec_Event.h
      Scenario.h
      ScenarioContext.h
      Singleerrorrate.h
      Utils.h
      ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_PROJECT_NAME}${PROJECT_NAME}Export.h
      ${CMAKE_CURRENT_BINARY_DIR}/Config
      ${CMAKE_CURRENT_BINARY_DIR}/Config.h
)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Config.h")
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Config")

target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src>
)

target_link_libraries(${PROJECT_NAME} PUBLIC GSL::gsl)
if(APPLE)
  target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP)
else()
  target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
endif()

if(WIN32)
	target_link_libraries(${PROJECT_NAME} PRIVATE ws2_32)
endif()

if(IGOR_WITH_SEQAN2)
  target_sources(${PROJECT_NAME} PRIVATE SeqAn2Aligner.cpp)
  target_include_directories(${PROJECT_NAME} PRIVATE ${SEQAN_INCLUDE_DIRS})
  target_compile_definitions(${PROJECT_NAME} PUBLIC IGOR_WITH_SEQAN2 ${SEQAN_DEFINITIONS})
  target_link_libraries(${PROJECT_NAME} PRIVATE ${SEQAN_LIBRARIES})
  if(IGOR_SEQAN2_SIMD)
    target_compile_options(${PROJECT_NAME} PRIVATE -msse4)
  endif()
endif()

## ###################################################################
## Install rules
## ###################################################################

install(
  TARGETS
    ${PROJECT_NAME}
  EXPORT
    ${CMAKE_PROJECT_NAME}${PROJECT_NAME}-targets
  LIBRARY
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE
    DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME
    DESTINATION ${CMAKE_INSTALL_BINDIR}
  FILE_SET HEADERS
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}/${PROJECT_NAME}
  INCLUDES
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
  EXPORT      ${CMAKE_PROJECT_NAME}${PROJECT_NAME}-targets
  FILE        ${CMAKE_PROJECT_NAME}${PROJECT_NAME}Targets.cmake
  NAMESPACE   ${CMAKE_PROJECT_NAME}::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${CMAKE_PROJECT_NAME}
)

export(
  EXPORT    ${CMAKE_PROJECT_NAME}${PROJECT_NAME}-targets
  NAMESPACE ${CMAKE_PROJECT_NAME}::
  FILE      ${CMAKE_BINARY_DIR}/${CMAKE_PROJECT_NAME}${PROJECT_NAME}Targets.cmake
)

######################################################################
### CMakeLists.txt ends here
