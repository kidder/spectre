# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "Test_Elasticity")

set(LIBRARY_SOURCES
  Test_Equations.cpp
  Test_Tags.cpp
  )

add_test_library(
  ${LIBRARY}
  "Elliptic/Systems/Elasticity/"
  "${LIBRARY_SOURCES}"
  "ConstitutiveRelations;DataStructures;Domain;Elasticity;LinearOperators;Utilities"
  )

add_subdirectory(BoundaryConditions)
