// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Helpers/Evolution/Wave.hpp"

#include <array>
#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"

namespace TestHelpers::evolution::wave {
template <size_t Dim>
std::array<double, Dim> k() {
  if constexpr (Dim == 1) {
    return std::array{1.0};
  } else if constexpr (Dim == 2) {
    return std::array{0.6, 0.8};
  } else {
    return std::array{3.0 / 13.0, 4.0 / 13.0, 12.0 / 13.0};
  }
}

template <size_t Dim>
std::array<double, Dim> comoving_v() {
  return c_s * k<Dim>();
}

template <size_t Dim>
DataVector u(const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
             const double t) {
  if constexpr (Dim == 1) {
    return DataVector{k<1>()[0] * get<0>(x) - c_s * t};
  } else if constexpr (Dim == 2) {
    return DataVector{k<2>()[0] * get<0>(x) + k<2>()[1] * get<1>(x) - c_s * t};
  } else {
    return DataVector{k<3>()[0] * get<0>(x) + k<3>()[1] * get<1>(x) +
                      k<3>()[2] * get<2>(x) - c_s * t};
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

template std::array<double, 1> k<1>();
template std::array<double, 2> k<2>();
template std::array<double, 3> k<3>();

template std::array<double, 1> comoving_v<1>();
template std::array<double, 2> comoving_v<2>();
template std::array<double, 3> comoving_v<3>();

template DataVector u(const tnsr::I<DataVector, 1, Frame::Inertial>& x,
                      const double t);
template DataVector u(const tnsr::I<DataVector, 2, Frame::Inertial>& x,
                      const double t);
template DataVector u(const tnsr::I<DataVector, 3, Frame::Inertial>& x,
                      const double t);

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
}  // namespace TestHelpers::evolution::wave
