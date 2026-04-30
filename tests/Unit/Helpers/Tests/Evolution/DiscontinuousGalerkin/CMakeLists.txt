# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "Test_EvolutionDgHelpers")

set(LIBRARY_SOURCES
  Test_ConservativeSystem.cpp
  Test_Variables.cpp
  )

add_test_library(${LIBRARY} "${LIBRARY_SOURCES}")

target_link_libraries(
  ${LIBRARY}
  PRIVATE
  CoordinateMaps
  DataStructures
  Domain
  DomainCreators
  DomainStructure
  EvolutionDgActionsHelpers
  LinearOperators
  Spectral
  )
