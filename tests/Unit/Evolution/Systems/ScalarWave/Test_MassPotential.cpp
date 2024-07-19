// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <random>
#include <string>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/Systems/ScalarWave/MassPotential.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Framework/CheckWithRandomValues.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Helpers/DataStructures/DataBox/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeWithRandomValues.hpp"

namespace {
template <typename DataType>
void test_compute_item_in_databox(const DataType& used_for_size) {
  TestHelpers::db::test_compute_tag<ScalarWave::Tags::MassPotentialCompute>(
      "MassPotential");

  MAKE_GENERATOR(generator);
  std::uniform_real_distribution<> dist(-1., 1.);

  const Scalar<DataVector> psi = make_with_random_values<Scalar<DataVector>>(
      make_not_null(&generator), dist, used_for_size);
  const double mass = make_with_random_values<double>(make_not_null(&generator),
                                                      dist, used_for_size);

  const auto box = db::create<
      db::AddSimpleTags<ScalarWave::Tags::Psi, ScalarWave::Tags::MassValue>,
      db::AddComputeTags<ScalarWave::Tags::MassPotentialCompute>>(psi, mass);

  const auto expected = ScalarWave::mass_potential(psi, mass);

  CHECK_ITERABLE_APPROX((db::get<ScalarWave::Tags::MassPotential>(box)),
                        expected);
}

void test_mass_potential(const DataVector& used_for_size) {
  void (*f)(const gsl::not_null<Scalar<DataVector>*>, const Scalar<DataVector>&,
            const double&) = &ScalarWave::mass_potential;
  pypp::check_with_random_values<1>(f, "MassPotential", {"mass_potential"},
                                    {{{-1., 1.}}}, used_for_size);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.Systems.ScalarWave.MassPotential",
                  "[Unit][Evolution]") {
  pypp::SetupLocalPythonEnvironment local_python_env(
      "Evolution/Systems/ScalarWave/");

  const DataVector used_for_size(5);
  test_mass_potential(used_for_size);
}
