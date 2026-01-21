// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/Variables.hpp"

namespace TestHelpers::evolution::dg::Actions {
namespace wave {
template <>
constexpr std::array<double, 1> k<1> = std::array{1.0};

template <>
constexpr std::array<double, 2> k<2> = std::array{0.6, 0.8};

template <>
constexpr std::array<double, 3> k<3> =
    std::array{3.0 / 13.0, 4.0 / 13.0, 12.0 / 13.0};

template <size_t Dim>
DataVector u(const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
             const double t) {
  if constexpr (Dim == 1) {
    return DataVector{k<1>[0] * get<0>(x) - c_s * t};
  } else if constexpr (Dim == 2) {
    return DataVector{k<2>[0] * get<0>(x) + k<2>[1] * get<1>(x) - c_s * t};
  } else {
    return DataVector{k<3>[0] * get<0>(x) + k<3>[1] * get<1>(x) +
                      k<3>[2] * get<2>(x) - c_s * t};
  }
}

template <size_t Dim>
DataVector f(const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
             const double t) {
  return DataVector{0.01 * square(u(x, t)) + 1.0};
}

template <size_t Dim>
DataVector df(const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
              const double t) {
  return DataVector{0.02 * u(x, t)};
}

template DataVector f(const tnsr::I<DataVector, 1, Frame::Inertial>& x,
                      const double t);
template DataVector f(const tnsr::I<DataVector, 2, Frame::Inertial>& x,
                      const double t);
template DataVector f(const tnsr::I<DataVector, 3, Frame::Inertial>& x,
                      const double t);

template DataVector df(const tnsr::I<DataVector, 1, Frame::Inertial>& x,
                      const double t);
template DataVector df(const tnsr::I<DataVector, 2, Frame::Inertial>& x,
                      const double t);
template DataVector df(const tnsr::I<DataVector, 3, Frame::Inertial>& x,
                      const double t);
}  // namespace wave
}  // namespace TestHelpers::evolution::dg::Actions
