// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>

#include "DataStructures/Tensor/TypeAliases.hpp"

/// \cond
class DataVector;
namespace domain::creators::time_dependence {
template <size_t Dim>
class TimeDependence;
}  // namespace domain::creators::time_dependence
/// \endcond

namespace TestHelpers::evolution::grid {
/// Initial expansion factor \f$a_0 = 1\f$
static constexpr double a_0 = 1.0;
/// Expansion rate \f$\frac{da}{dt} = \frac{1}{2}\f$
static constexpr double a_dot = 0.5;

/// \brief Motion of the grid in the inertial frame
///
/// \details For Comoving, the grid velocity corresponds to that of
/// TestHelpers::evolution::wave::comoving_v()
enum class Is : uint8_t {
  Stationary, /**< \f$\vec{v_g} = 0 */
  Comoving,   /**< \f$\vec{v_g} = c_s \vec{k}\f$ */
  Expanding   /**< \f$\vec{v_g} = \frac{1}{a} \frac{da}{dt} \vec{x_i}\f$ */
};

/// The grid velocity
template <Is grid_is, size_t Dim>
std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>> v(
    const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t);

/// The inertial divergence of the grid velocity
template <Is grid_is, size_t Dim>
std::optional<Scalar<DataVector>> div_v(
    const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t);

/// Create a concrete TimeDependence
template <Is grid_is, size_t Dim>
std::unique_ptr<domain::creators::time_dependence::TimeDependence<Dim>>
create_time_dependence();
}  // namespace TestHelpers::evolution::grid
