// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/ComputeTimeDerivative.tpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/SystemType.hpp"

namespace TestHelpers::evolution::dg::Actions {
template void test<3, SystemType::Conservative, false>();
}  // namespace TestHelpers::evolution::dg::Actions
