// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>

#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/ComputeTimeDerivative.hpp"

namespace TestHelpers::evolution::dg::Actions {
template <size_t Dim, bool PassVariables>
struct ConservativeSystem;

namespace {
template <bool PassVariables, size_t Dim>
void test_dim() {
  TestHelpers::evolution::dg::Actions::test<
      Dim, TestHelpers::evolution::dg::Actions::ConservativeSystem<
               Dim, PassVariables>>();
  // TestHelpers::evolution::dg::Actions::test<
  //     Dim, TestHelpers::evolution::dg::Actions::ConservativeSystemWithPrims<
  //              Dim, PassVariables>>();
  // TestHelpers::evolution::dg::Actions::test<
  //     Dim, TestHelpers::evolution::dg::Actions::NonconservativeSystem<
  //              Dim, PassVariables>>();
  // TestHelpers::evolution::dg::Actions::test<
  //     Dim, TestHelpers::evolution::dg::Actions::MixedSystem<
  //              Dim, PassVariables>>();
  // TestHelpers::evolution::dg::Actions::test<
  //     Dim, TestHelpers::evolution::dg::Actions::MixedSystemWithPrims<
  //              Dim, PassVariables>>();
}

template <bool PassVariables>
void test_pass_variables() {
  test_dim<PassVariables, 1>();
  // test_dim<PassVariables, 2>();
  // test_dim<PassVariables, 3>();
}
}  // namespace
}  // namespace TestHelpers::evolution::dg::Actions

SPECTRE_TEST_CASE("Unit.Evolution.DG.NonconformingBlocks",
                  "[Unit][Evolution][Actions]") {
  TestHelpers::evolution::dg::Actions::test_pass_variables<false>();
  //  TestHelpers::evolution::dg::Actions::test_pass_variables<true>();
}
