// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/ComputeTimeDerivative.tpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/NonconservativeSystem.hpp"

namespace TestHelpers::evolution::dg::Actions {
template void test<1, NonconservativeSystem<1, false>>();
template void test<1, NonconservativeSystem<1, true>>();
}  // namespace TestHelpers::evolution::dg::Actions
