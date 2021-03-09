# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBS_TO_LINK
  Amr
  DomainCreators
  Informer
  Options
  Utilities
  Amr
  )

add_spectre_parallel_executable(
  ExploreAmr1D
  ExploreAmr
  Executables/Examples/ExploreAmr
  ExploreAmrMetavars<1>
  "${LIBS_TO_LINK}"
  )

add_spectre_parallel_executable(
  ExploreAmr2D
  ExploreAmr
  Executables/Examples/ExploreAmr
  ExploreAmrMetavars<2>
  "${LIBS_TO_LINK}"
  )

add_spectre_parallel_executable(
  ExploreAmr3D
  ExploreAmr
  Executables/Examples/ExploreAmr
  ExploreAmrMetavars<3>
  "${LIBS_TO_LINK}"
  )
