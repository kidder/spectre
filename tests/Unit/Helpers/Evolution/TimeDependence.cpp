// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Helpers/Evolution/TimeDependence.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <optional>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Creators/TimeDependence/CubicScale.hpp"
#include "Domain/Creators/TimeDependence/None.hpp"
#include "Domain/Creators/TimeDependence/TimeDependence.hpp"
#include "Domain/Creators/TimeDependence/UniformTranslation.hpp"
#include "Helpers/Evolution/Wave.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

namespace TestHelpers::evolution::grid {
template <Is grid_is, size_t Dim>
std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>> v(
    const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t) {
  if constexpr (grid_is == Is::Stationary) {
    return std::nullopt;
  } else if constexpr (grid_is == Is::Comoving) {
    auto result = make_with_value<tnsr::I<DataVector, Dim, Frame::Inertial>>(
        x, wave::c_s);
    for (size_t d = 0; d < Dim; ++d) {
      result.get(d) *= gsl::at(wave::k<Dim>(), d);
    }
    return result;
  } else {
    tnsr::I<DataVector, Dim, Frame::Inertial> result = x;
    const double a = a_0 + a_dot * t;
    for (size_t d = 0; d < Dim; ++d) {
      result.get(d) *= a_dot / a;
    }
    return result;
  }
}

template <Is grid_is, size_t Dim>
std::optional<Scalar<DataVector>> div_v(
    const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t) {
  if constexpr (grid_is == Is::Stationary) {
    return std::nullopt;
  } else if constexpr (grid_is == Is::Comoving) {
    return make_with_value<Scalar<DataVector>>(x, 0.0);
  } else {
    return make_with_value<Scalar<DataVector>>(
        x, a_dot * static_cast<double>(Dim) / (a_0 + a_dot * t));
  }
}

template <Is grid_is, size_t Dim>
std::unique_ptr<domain::creators::time_dependence::TimeDependence<Dim>>
create_time_dependence() {
  if constexpr (grid_is == Is::Stationary) {
    return std::make_unique<domain::creators::time_dependence::None<Dim>>();
  }
  if constexpr (grid_is == Is::Expanding) {
    return std::make_unique<domain::creators::time_dependence::CubicScale<Dim>>(
        0.0, std::numeric_limits<double>::infinity(), true,
        std::array{a_0, a_0}, std::array{a_dot, a_dot}, std::array{0.0, 0.0});
  }
  if constexpr (grid_is == Is::Comoving) {
    return std::make_unique<
        domain::creators::time_dependence::UniformTranslation<Dim>>(
        0.0, wave::comoving_v<Dim>());
  }
}

template std::optional<tnsr::I<DataVector, 1, Frame::Inertial>>
v<Is::Stationary>(const tnsr::I<DataVector, 1, Frame::Inertial>& x,
                  const double t);
template std::optional<tnsr::I<DataVector, 2, Frame::Inertial>>
v<Is::Stationary>(const tnsr::I<DataVector, 2, Frame::Inertial>& x,
                  const double t);
template std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>
v<Is::Stationary>(const tnsr::I<DataVector, 3, Frame::Inertial>& x,
                  const double t);
template std::optional<tnsr::I<DataVector, 1, Frame::Inertial>> v<Is::Comoving>(
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t);
template std::optional<tnsr::I<DataVector, 2, Frame::Inertial>> v<Is::Comoving>(
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t);
template std::optional<tnsr::I<DataVector, 3, Frame::Inertial>> v<Is::Comoving>(
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t);
template std::optional<tnsr::I<DataVector, 1, Frame::Inertial>>
v<Is::Expanding>(const tnsr::I<DataVector, 1, Frame::Inertial>& x,
                 const double t);
template std::optional<tnsr::I<DataVector, 2, Frame::Inertial>>
v<Is::Expanding>(const tnsr::I<DataVector, 2, Frame::Inertial>& x,
                 const double t);
template std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>
v<Is::Expanding>(const tnsr::I<DataVector, 3, Frame::Inertial>& x,
                 const double t);

template std::optional<Scalar<DataVector>> div_v<Is::Stationary>(
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t);
template std::optional<Scalar<DataVector>> div_v<Is::Stationary>(
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t);
template std::optional<Scalar<DataVector>> div_v<Is::Stationary>(
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t);
template std::optional<Scalar<DataVector>> div_v<Is::Comoving>(
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t);
template std::optional<Scalar<DataVector>> div_v<Is::Comoving>(
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t);
template std::optional<Scalar<DataVector>> div_v<Is::Comoving>(
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t);
template std::optional<Scalar<DataVector>> div_v<Is::Expanding>(
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t);
template std::optional<Scalar<DataVector>> div_v<Is::Expanding>(
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t);
template std::optional<Scalar<DataVector>> div_v<Is::Expanding>(
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t);

template std::unique_ptr<domain::creators::time_dependence::TimeDependence<1>>
create_time_dependence<Is::Stationary, 1>();
template std::unique_ptr<domain::creators::time_dependence::TimeDependence<2>>
create_time_dependence<Is::Stationary, 2>();
template std::unique_ptr<domain::creators::time_dependence::TimeDependence<3>>
create_time_dependence<Is::Stationary, 3>();
template std::unique_ptr<domain::creators::time_dependence::TimeDependence<1>>
create_time_dependence<Is::Comoving, 1>();
template std::unique_ptr<domain::creators::time_dependence::TimeDependence<2>>
create_time_dependence<Is::Comoving, 2>();
template std::unique_ptr<domain::creators::time_dependence::TimeDependence<3>>
create_time_dependence<Is::Comoving, 3>();
template std::unique_ptr<domain::creators::time_dependence::TimeDependence<1>>
create_time_dependence<Is::Expanding, 1>();
template std::unique_ptr<domain::creators::time_dependence::TimeDependence<2>>
create_time_dependence<Is::Expanding, 2>();
template std::unique_ptr<domain::creators::time_dependence::TimeDependence<3>>
create_time_dependence<Is::Expanding, 3>();
}  // namespace TestHelpers::evolution::grid
