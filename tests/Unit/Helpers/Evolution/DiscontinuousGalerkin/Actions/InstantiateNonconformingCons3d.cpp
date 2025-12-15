// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/NonconformingComputeTimeDerivativeImpl.tpp"

namespace TestHelpers::evolution::dg::Actions {
  template void test_nonconforming<SystemType::Conservative, false, 3>();
  template void test_nonconforming<SystemType::Conservative, false, 2>();
  template void test_nonconforming<SystemType::Conservative, false, 1>();
}  // namespace TestHelpers::evolution::dg::Actions
