// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/ComputeTimeDerivative.tpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/ConservativeSystem.hpp"

namespace TestHelpers::evolution::dg::Actions {
template void test<1, ConservativeSystem<1, false>>();
template void test<1, ConservativeSystem<1, true>>();
}  // namespace TestHelpers::evolution::dg::Actions
