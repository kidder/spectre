// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/ComputeTimeDerivative.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/SystemType.hpp"

namespace {
template <TestHelpers::evolution::dg::Actions::SystemType system_type,
          bool UsePrims>
void test() {
  TestHelpers::evolution::dg::Actions::test<system_type, UsePrims, 1>();
  // TestHelpers::evolution::dg::Actions::test<system_type, UsePrims, 2>();
  // TestHelpers::evolution::dg::Actions::test<system_type, UsePrims, 3>();
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.DG.ComputeTimeDerivative",
                  "[Unit][Evolution][Actions]") {
  test<TestHelpers::evolution::dg::Actions::SystemType::Conservative, false>();
}
