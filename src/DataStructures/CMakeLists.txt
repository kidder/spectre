# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY DataStructures)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  ApplyMatrices.cpp
  DynamicBuffer.cpp
  FloatingPointType.cpp
  Index.cpp
  IndexIterator.cpp
  LeviCivitaIterator.cpp
  SliceIterator.cpp
  StripeIterator.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  ApplyMatrices.hpp
  BoostMultiArray.hpp
  CachedTempBuffer.hpp
  ComplexDataVector.hpp
  ComplexDiagonalModalOperator.hpp
  ComplexModalVector.hpp
  DataVector.hpp
  DenseMatrix.hpp
  DenseVector.hpp
  DiagonalModalOperator.hpp
  DynamicBuffer.hpp
  FixedHashMap.hpp
  FloatingPointType.hpp
  GeneralIndexIterator.hpp
  IdPair.hpp
  Index.hpp
  IndexIterator.hpp
  LeviCivitaIterator.hpp
  Matrix.hpp
  ModalVector.hpp
  SliceIterator.hpp
  SliceTensorToVariables.hpp
  SliceVariables.hpp
  SpinWeighted.hpp
  StripeIterator.hpp
  Tags.hpp
  TempBuffer.hpp
  Transpose.hpp
  Variables.hpp
  VariablesTag.hpp
  VectorImpl.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  Blaze
  Boost::boost
  ErrorHandling
  Options
  Utilities
  )

add_subdirectory(DataBox)
add_subdirectory(Python)
add_subdirectory(Tensor)
