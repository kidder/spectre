// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/NonconformingComputeTimeDerivativeImpl.tpp"

namespace TestHelpers::evolution::dg::Actions {
template void test_nonconforming<SystemType::Conservative, false, 3>();
}  // namespace TestHelpers::evolution::dg::Actions
