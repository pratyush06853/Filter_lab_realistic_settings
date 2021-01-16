# $Id: CMakeLists.txt 86065 2014-11-07 08:51:15Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(IronFilter)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Find ROOT
IF(NOT DEFINED ROOT_FOUND)
  message("Root not found... Let's try Find_package ;)")
  find_package(ROOT 6.18 CONFIG REQUIRED)
  #find_package(ROOT COMPONENTS MathMore Matrix Hist Tree Graf Gpad)
  IF(NOT ${ROOT_FOUND})
      message("Root not found... have to try an other way ;)")
      include($ENV{ROOTSYS}/etc/root/cmake/FindROOT.cmake)
    include_directories( PUBLIC ${ROOT_INCLUDE_DIRS})
    link_libraries( ${ROOT_LIBRARIES})
  ELSE(NOT ${ROOT_FOUND})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
      include(${ROOT_USE_FILE})
    include_directories(PUBLIC ${ROOT_INCLUDE_DIRS})
    link_libraries( ${ROOT_LIBRARIES})
  ENDIF()
ENDIF()
#-----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(IronFilter IronFilter.cc ${sources} ${headers})
target_link_libraries(IronFilter ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build IronFilter. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(IronFilter_SCRIPTS
  IronFilter.out
  IronFilter.mac
  IronFilterBasicTest.mac
  icons.mac
  gui.mac
  run.png
  init_vis.mac
  vis.mac
  )

foreach(_script ${IronFilter_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS IronFilter DESTINATION bin)
