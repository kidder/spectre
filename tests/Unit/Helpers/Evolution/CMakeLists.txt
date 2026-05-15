# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "EvolutionHelpers")

set(LIBRARY_SOURCES
  DomainCreator.cpp
  TimeDependence.cpp
  Wave.cpp
  )

add_spectre_library(${LIBRARY} ${SPECTRE_TEST_LIBS_TYPE} ${LIBRARY_SOURCES})

target_link_libraries(
  ${LIBRARY}
  PRIVATE
  DataStructures
  DomainCreators
  )

add_subdirectory(DgSubcell)
add_subdirectory(DiscontinuousGalerkin)
add_subdirectory(Imex)
add_subdirectory(Systems)
