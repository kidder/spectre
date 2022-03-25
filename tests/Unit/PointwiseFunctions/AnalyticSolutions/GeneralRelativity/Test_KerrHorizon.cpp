// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cmath>
#include <limits>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Framework/CheckWithRandomValues.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrHorizon.hpp"
#include "Utilities/ConstantExpressions.hpp"

namespace gr {
namespace Solutions {
namespace {

// This wraps the kerr_horizon_radius function with
// types that CheckWithRandomValues can understand.
// We use magnitude and angle of dimensionless spin so that
// we know their bounds.
template <typename DataType>
Scalar<DataType> wrap_kerr_horizon_radius(const Scalar<DataType>& theta,
                                          const Scalar<DataType>& phi,
                                          const double mass,
                                          const double dimless_spin_magnitude,
                                          const double dimless_spin_theta,
                                          const double dimless_spin_phi) {
  return kerr_horizon_radius<DataType>(
      {{get(theta), get(phi)}}, mass,
      {{dimless_spin_magnitude * sin(dimless_spin_theta) *
            cos(dimless_spin_phi),
        dimless_spin_magnitude * sin(dimless_spin_theta) *
            sin(dimless_spin_phi),
        dimless_spin_magnitude * cos(dimless_spin_theta)}});
}

template <typename DataType>
void test_kerr_horizon(const DataType& used_for_size) {
  pypp::check_with_random_values<6>(&wrap_kerr_horizon_radius<DataType>,
                                    "TestFunctions", "kerr_horizon_radius",
                                    {{{0.0, M_PI},
                                      {0.0, 2.0 * M_PI},
                                      {1.0, 2.0},
                                      {0.0, 1.0},
                                      {0.0, M_PI},
                                      {0.0, 2.0 * M_PI}}},
                                    used_for_size);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.AnalyticSolutions.Gr.KerrHorizon",
                  "[PointwiseFunctions][Unit]") {
  pypp::SetupLocalPythonEnvironment local_python_env(
      "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/");
  const DataVector dv(5);
  test_kerr_horizon(dv);
  test_kerr_horizon(0.0);

  // Test for Schwarzschild (mass=2) at randomly-chosen point
  CHECK(4.0 == get(kerr_horizon_radius<double>({{1.12345, 2.2222}}, 2.0,
                                               {{0.0, 0.0, 0.0}})));

  // Test for Kerr (mass=2) along pole, for spin not in z direction.
  CHECK(approx(2.0 * (1.0 + sqrt(0.86))) ==
        get(kerr_horizon_radius<double>(
            {{acos(0.3 / sqrt(0.14)), atan2(0.2, 0.1)}}, 2.0,
            {{0.1, 0.2, 0.3}})));

  // Test for Kerr (mass=2) along equator, for spin not in z direction.
  // Angles and radius worked out by hand.
  CHECK(approx(sqrt(square(2.0 * (1.0 + sqrt(0.86))) + 4.0 * 0.14)) ==
        get(kerr_horizon_radius<double>(
            {{acos(0.3 / sqrt(0.14)) + M_PI_2, atan2(0.2, 0.1)}}, 2.0,
            {{0.1, 0.2, 0.3}})));

  // Test for Kerr (mass=2) along pole, for extremal spin not in z direction.
  CHECK(approx(2.0) ==
        get(kerr_horizon_radius<double>(
            {{acos(0.3 / sqrt(0.14)), atan2(0.2, 0.1)}}, 2.0,
            {{0.1 / sqrt(0.14), 0.2 / sqrt(0.14), 0.3 / sqrt(0.14)}})));

  // Test for Kerr (mass=2) along equator, for extremal spin not in z direction.
  CHECK(approx(sqrt(8.0)) ==
        get(kerr_horizon_radius<double>(
            {{acos(0.3 / sqrt(0.14)) + M_PI_2, atan2(0.2, 0.1)}}, 2.0,
            {{0.1 / sqrt(0.14), 0.2 / sqrt(0.14), 0.3 / sqrt(0.14)}})));

  KerrHorizon kh_extremal{4.0, {{3.0 / 13.0, 4.0 / 13.0, 12.0 / 13.0}}};
  CHECK(kh_extremal.mass == 4.0);
  CHECK(kh_extremal.dimensionless_spin_magnitude == 1.0);
  CHECK(kh_extremal.polar_radius == 4.0);
  CHECK(kh_extremal.equatorial_radius == 4.0 * sqrt(2.0));
  CHECK(approx(kh_extremal.surface_area()) == 128.0 * M_PI);

  KerrHorizon kh{2.0, {{2.4 / 13.0, 3.2 / 13.0, 9.6 / 13.0}}};
  CHECK(kh.mass == 2.0);
  CHECK(approx(kh.dimensionless_spin_magnitude) == 0.8);
  CHECK(kh.polar_radius == 3.2);
  CHECK(approx(kh.equatorial_radius) == 2.0 * sqrt(3.2));
  CHECK(approx(kh.surface_area()) == 51.2 * M_PI);
}
}  // namespace Solutions
}  // namespace gr
