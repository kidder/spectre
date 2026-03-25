// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/ComputeTimeDerivative.tpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/MixedSystemWithPrims.hpp"

namespace TestHelpers::evolution::dg::Actions {
template void test<3, MixedSystemWithPrims<3, false>>();
template void test<3, MixedSystemWithPrims<3, true>>();
}  // namespace TestHelpers::evolution::dg::Actions
