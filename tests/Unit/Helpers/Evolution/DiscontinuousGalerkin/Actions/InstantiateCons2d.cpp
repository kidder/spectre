// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/ComputeTimeDerivative.tpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/SystemType.hpp"

namespace TestHelpers::evolution::dg::Actions {
template void test<SystemType::Conservative, false, 2>();
}  // namespace TestHelpers::evolution::dg::Actions
