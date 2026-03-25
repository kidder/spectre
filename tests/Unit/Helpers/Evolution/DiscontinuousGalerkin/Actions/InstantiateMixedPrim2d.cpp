// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/ComputeTimeDerivative.tpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/MixedSystemWithPrims.hpp"

namespace TestHelpers::evolution::dg::Actions {
template void test<2, MixedSystemWithPrims<2, false>>();
template void test<2, MixedSystemWithPrims<2, true>>();
}  // namespace TestHelpers::evolution::dg::Actions
