// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/KerrHorizon.hpp"

#include <array>
#include <cmath>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/ConstantExpressions.hpp"

namespace gr {
namespace Solutions {

template <typename DataType>
Scalar<DataType> kerr_horizon_radius(
    const std::array<DataType, 2>& theta_phi, const double mass,
    const std::array<double, 3>& dimensionless_spin) {

  const auto& theta = theta_phi[0];
  const auto& phi = theta_phi[1];

  const DataType ax_times_nx =
      (mass * dimensionless_spin[0]) * sin(theta) * cos(phi);
  const DataType ay_times_ny =
      (mass * dimensionless_spin[1]) * sin(theta) * sin(phi);
  const DataType az_times_nz = (mass * dimensionless_spin[2]) * cos(theta);

  KerrHorizon kh{mass, dimensionless_spin};

  const DataType denominator =
      sqrt(square(ax_times_nx + ay_times_ny + az_times_nz) +
           square(kh.polar_radius));

  return Scalar<DataType>{kh.equatorial_radius * kh.polar_radius / denominator};
}

template Scalar<DataVector> kerr_horizon_radius(
    const std::array<DataVector, 2>& theta_phi, const double mass,
    const std::array<double, 3>& dimensionless_spin);

template Scalar<double> kerr_horizon_radius(
    const std::array<double, 2>& theta_phi, const double mass,
    const std::array<double, 3>& dimensionless_spin);

KerrHorizon::KerrHorizon(double mass_in,
                         const std::array<double, 3>& dimensionless_spin_in)
    : mass(mass_in),
      dimensionless_spin(dimensionless_spin_in),
      dimensionless_spin_magnitude(std::hypot(
          dimensionless_spin[0], dimensionless_spin[1], dimensionless_spin[2])),
      polar_radius(mass * (1.0 + sqrt((1.0 + dimensionless_spin_magnitude) *
                                      (1.0 - dimensionless_spin_magnitude)))),
      equatorial_radius(sqrt(2.0 * mass * polar_radius)) {}
}  // namespace Solutions
}  // namespace gr

