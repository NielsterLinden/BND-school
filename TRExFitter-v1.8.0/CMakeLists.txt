# Set the minimum CMake version required to build the project.
cmake_minimum_required( VERSION 3.10 )

# Silence some warnings on macOS with new CMake versions.
if( NOT ${CMAKE_VERSION} VERSION_LESS 3.9 )
   cmake_policy( SET CMP0068 NEW )
endif()

# Set the project's name and version.
project( TRExFitter )

# Enable using CTest in the project.
include( CTest )

# Set up the "C++ version" to use.
set( CMAKE_CXX_STANDARD_REQUIRED 11 CACHE STRING
   "Minimum C++ standard required for the build" )
set( CMAKE_CXX_STANDARD 17 CACHE STRING
   "C++ standard to use for the build" )
set( CMAKE_CXX_EXTENSIONS FALSE CACHE BOOL
   "(Dis)allow using compiler extensions" )

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Specify the install locations for libraries and binaries.
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_INSTALL_LIBDIR ${CMAKE_BINARY_DIR}/lib )  # Needed by ROOT_GENERATE_DICTIONARY()

# Set the warning flag(s) to use.
set( CMAKE_CXX_FLAGS "-Wall -Wextra -Wshadow -pedantic -O3 -g -fPIC" )

set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

set(USE_LOCAL_XROOFIT TRUE CACHE BOOL "If true will use xroofit submodule")

set(TREXFITTER_AVX_OPTIMISATION TRUE CACHE BOOL "If true use AVX optimisation")

# ------------------------------------------------
# Dependencies and sub-libraries
# ------------------------------------------------


# use local xroofit if it is selected if the one from StatAnalysis does not exist
if( NOT USE_LOCAL_XROOFIT )
  # for some reason if this is done after find_package(ROOT ... ) the linking is incorrect
  find_package( xRooFit QUIET )
endif()
if( NOT xRooFit_FOUND )
  include_directories( xroofit )
  add_subdirectory( xroofit )
endif()


# Add ROOT system directory and require ROOT.
find_package( ROOT 6.24.00 REQUIRED COMPONENTS Core MathCore HistFactory Graf Hist RIO Tree Gpad )
find_package( AnalysisBase QUIET)

set( YAML_CPP_BUILD_TESTS OFF )
set( YAML_CPP_BUILD_TOOLS OFF )

# include yaml-cpp library
include_directories( yaml-cpp )
if(NOT ${AnalysisBase_FOUND})
  add_subdirectory( yaml-cpp )
endif()

# Handle the build of the SystematicSmoothingTool library. That
# library comes with a cmake configuration, but introduces
# dependencies on ATLAS cmake syntax – which we are trying to
# avoid here. So ignore the shipped CMakeLists.txt, and build our
# own version of SystematicSmoothingTool.
include_directories( SystematicSmoothingTool )
set( _syst_smoothing_tool_sources
  SystematicSmoothingTool/Root/PlotHist.cxx
  SystematicSmoothingTool/Root/SmoothHist.cxx
  SystematicSmoothingTool/Root/SmoothingTool.cxx
  SystematicSmoothingTool/Root/TRExTools.cxx
)
add_library( SystematicSmoothingTool SHARED ${_syst_smoothing_tool_sources} )
target_include_directories( SystematicSmoothingTool
   PUBLIC ${ROOT_INCLUDE_DIRS}
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:> )
target_link_libraries( SystematicSmoothingTool ${ROOT_LIBRARIES} )

# Add Unfolding code
include_directories( UnfoldingCode )
add_subdirectory( UnfoldingCode )

# Add the blinder code
include_directories( Blinder )
add_subdirectory( Blinder )

# Handle the build of the three auxiliary classes for TRExFitter style
# labels etc. We don't want to throw this into the TRExFitter
# library object.
set( _trex_utils_headers StyleUtils/TRExLabels.h StyleUtils/TRExStyle.h StyleUtils/TRExUtils.h )
set( _trex_utils_sources StyleUtils/TRExLabels.C StyleUtils/TRExStyle.C StyleUtils/TRExUtils.C )
add_library( TRExUtils SHARED ${_trex_utils_headers} ${_trex_utils_sources} )
target_include_directories( TRExUtils
   PUBLIC ${ROOT_INCLUDE_DIRS}
   $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:> )
target_link_libraries( TRExUtils ${ROOT_LIBRARIES} )


# ------------------------------------------------
# The actual TRExFitter library
# ------------------------------------------------

# Public header files for the shared/static library.
set( lib_headers
  TRExFitter/Common.h
  TRExFitter/ConfigParser.h
  TRExFitter/ConfigReader.h
  TRExFitter/ConfigReaderMulti.h
  TRExFitter/CorrelationMatrix.h
  TRExFitter/EFTConfig.h
  TRExFitter/EFTProcessor.h
  TRExFitter/FitResults.h
  TRExFitter/FittingTool.h
  TRExFitter/FitToys.h
  TRExFitter/FitUtils.h
  TRExFitter/HistoReader.h
  TRExFitter/HistoTools.h
  TRExFitter/LikelihoodScanManager.h
  TRExFitter/LimitEvent.h
  TRExFitter/LimitToys.h
  TRExFitter/Logger.h
  TRExFitter/Morphing.h
  TRExFitter/MultiFit.h
  TRExFitter/NormFactor.h
  TRExFitter/NtupleReader.h
  TRExFitter/NuisParameter.h
  TRExFitter/PruningUtil.h
  TRExFitter/RankingManager.h
  TRExFitter/ReparametrizationManager.h
  TRExFitter/Region.h
  TRExFitter/Sample.h
  TRExFitter/SampleHist.h
  TRExFitter/SignificanceToys.h
  TRExFitter/ShapeFactor.h
  TRExFitter/Systematic.h
  TRExFitter/SystematicHist.h
  TRExFitter/TemplateMorpher.h
  TRExFitter/TRExFit.h
  TRExFitter/TRExPlot.h
  TRExFitter/TruthSample.h
  TRExFitter/UncertaintyPlotter.h
  TRExFitter/Unfolding.h
  TRExFitter/UnfoldingSample.h
  TRExFitter/UnfoldingSystematic.h
  TRExFitter/YamlConverter.h)

# Source files for the shared/static library.
set( lib_sources
  Root/Common.cc
  Root/ConfigParser.cc
  Root/ConfigReader.cc
  Root/ConfigReaderMulti.cc
  Root/CorrelationMatrix.cc
  Root/EFTConfig.cc
  Root/EFTProcessor.cc
  Root/FitResults.cc
  Root/FittingTool.cc
  Root/FitToys.cc
  Root/FitUtils.cc
  Root/HistoReader.cc
  Root/HistoTools.cc
  Root/LikelihoodScanManager.cc
  Root/LimitEvent.cc
  Root/LimitToys.cc
  Root/Morphing.cc
  Root/MultiFit.cc
  Root/NormFactor.cc
  Root/NtupleReader.cc
  Root/NuisParameter.cc
  Root/PruningUtil.cc
  Root/RankingManager.cc
  Root/ReparametrizationManager.cc
  Root/Region.cc
  Root/Sample.cc
  Root/SampleHist.cc
  Root/SignificanceToys.cc
  Root/ShapeFactor.cc
  Root/Systematic.cc
  Root/SystematicHist.cc
  Root/TemplateMorpher.cc
  Root/TRExFit.cc
  Root/TRExPlot.cc
  Root/TruthSample.cc
  Root/UncertaintyPlotter.cc
  Root/Unfolding.cc
  Root/UnfoldingSample.cc
  Root/UnfoldingSystematic.cc
  Root/YamlConverter.cc )

# Build the shared library.
add_library( TRExFitter SHARED ${lib_headers} ${lib_sources} )
target_include_directories( TRExFitter PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:> ${ROOT_INCLUDE_DIRS} )
target_link_libraries( TRExFitter SystematicSmoothingTool TRExUtils UnfoldingCode Blinder xRooFit yaml-cpp ${ROOT_LIBRARIES} )
set_property( TARGET TRExFitter
   PROPERTY PUBLIC_HEADER ${lib_headers} )
target_include_directories(TRExFitter PUBLIC ${CMAKE_CURRENT_LIST_DIR} )

# Install the libraries.
install( TARGETS TRExFitter SystematicSmoothingTool TRExUtils
   EXPORT TRExFitterTargets
   ARCHIVE DESTINATION lib
   LIBRARY DESTINATION lib
   PUBLIC_HEADER DESTINATION include/TRExFitter )


# ------------------------------------------------
# TRExFitter executables
# ------------------------------------------------

# Helper macro for building the project's executables.
macro( TRExFitter_add_executable name )
  add_executable( ${name} ${ARGN} )
  target_include_directories( ${name} PUBLIC ${ROOT_INCLUDE_DIRS} )
  target_link_libraries( ${name} TRExFitter SystematicSmoothingTool TRExUtils UnfoldingCode Blinder xRooFit yaml-cpp ${ROOT_LIBRARIES} )
  install( TARGETS ${name}
    EXPORT TRExFitterTargets
    RUNTIME DESTINATION bin )
endmacro( TRExFitter_add_executable )

TRExFitter_add_executable( trex-fitter util/trex-fitter.cc )

# Compile the util/hupdate.cc executable.
add_executable( hupdate.exe util/hupdate.cc )
target_include_directories( hupdate.exe PUBLIC ${ROOT_INCLUDE_DIRS} )
target_link_libraries( hupdate.exe ${ROOT_LIBRARIES} )
install( TARGETS hupdate.exe
  EXPORT TRExFitterTargets
  RUNTIME DESTINATION bin )


# ------------------------------------------------
# Install target prescription
# ------------------------------------------------


# install( EXPORT TRExFitterTargets
#    FILE TRExFitterConfig-targets.cmake
#    DESTINATION cmake
#    NAMESPACE "TRExFitter::" )

# add compiler flag called TREXFITTER_HOME to not need to set environment variable any more
if(NOT TREXFITTER_NOHOME)
  add_definitions( -DTREXFITTER_HOME=${CMAKE_CURRENT_SOURCE_DIR} )
endif()

if(NOT TREXFITTER_AVX_OPTIMISATION)
  add_definitions( -DTREXFITTER_AVX_OPTIMISATION=0 )
else()
  add_definitions( -DTREXFITTER_AVX_OPTIMISATION=1 )
endif()
