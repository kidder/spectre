// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/ComputeTimeDerivative.tpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/ConservativeSystem.hpp"

namespace TestHelpers::evolution::dg::Actions {
template void test<1, ConservativeSystem<1, false>>();
template void test<1, ConservativeSystem<1, true>>();
template void test<2, ConservativeSystem<2, false>>();
template void test<2, ConservativeSystem<2, true>>();
template void test<3, ConservativeSystem<3, false>>();
template void test<3, ConservativeSystem<3, true>>();
}  // namespace TestHelpers::evolution::dg::Actions
