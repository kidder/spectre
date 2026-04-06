// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <optional>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Tags.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/Gsl.hpp"

/// \cond
namespace Parallel {
template <typename>
class GlobalCache;
}  // namespace Parallel
namespace tuples {
template <typename...>
class TaggedTuple;
}  // namespace tuples
/// \endcond

namespace TestHelpers::evolution::dg::Actions {

namespace wave {
/// The wave speed \f$c_s = \frac{1}{2}\f$
static constexpr double c_s = 0.5;

/// \brief The wave vector \f$\vec{k}\f$
///
/// \details If \f$d = 1\f$, \f$\vec{k} = (1.0)\f$.
/// If \f$d = 2\f$, \f$\vec{k} = (\frac{3}{5}, \frac{4}{5})\f$.
/// If \f$d = 3\f$, \f$\vec{k} = (\frac{3}{13}, \frac{4}{13}, \frac{12}{13})\f$.
template <size_t Dim>
std::array<double, Dim> k();

/// The velocity comoving with the wave \f$\vec{v} = c_s \vec{k}\f$
template <size_t Dim>
std::array<double, Dim> comoving_v();

/// The variables are given as functions of \f$u = k_i x^i - t\f$
template <size_t Dim>
DataVector u(const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
             const double t);

/// The profile function \f$f = 1 + \frac{1}{100} u^2\f$
template <size_t Dim>
DataVector f(const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
             const double t);

/// The derivative of the profile function \f$\frac{df}{du} = \frac{1}{50} u\f$
template <size_t Dim>
DataVector df(const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
              const double t);
}  // namespace wave

namespace grid {
/// Initial expansion factor \f$a_0 = 1\f$
static constexpr double a_0 = 1.0;
/// Expansion rate \f$\frac{da}{dt} = \frac{1}{2}\f$
static constexpr double a_dot = 0.5;

/// \brief Motion of the grid in the inertial frame
enum class Is : uint8_t {
  Stationary, /**< \f$\vec{v_g} = 0 */
  Comoving,   /**< \f$\vec{v_g} = c_s \vec{k}\f$ */
  Expanding   /**< \f$\vec{v_g} = \frac{1}{a} \frac{da}{dt} \vec{x_i}\f$ */
};

template <Is grid_is, size_t Dim>
std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>> v(
    const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t);

template <Is grid_is, size_t Dim>
std::optional<Scalar<DataVector>> div_v(
    const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t);
}  // namespace grid

struct Var1 : db::SimpleTag {
  using type = Scalar<DataVector>;

  template <size_t Dim>
  static void value(const gsl::not_null<type*> var1,
                    const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                    const double t);

  template <size_t Dim>
  static void flux(
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> flux1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& var2);

  template <grid::Is grid_is, size_t Dim>
  static void normal_dot_flux(
      const gsl::not_null<type*> normal_dot_flux,
      const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t,
      const double sign);

  template <grid::Is grid_is, size_t Dim>
  static void dt_value(const gsl::not_null<type*> dt_var1,
                       const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                       const double t);
};

template <size_t Dim>
struct Var2 : db::SimpleTag {
  using type = tnsr::I<DataVector, Dim, Frame::Inertial>;

  static void value(const gsl::not_null<type*> var2,
                    const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                    const double t);

  static void flux(
      const gsl::not_null<tnsr::IJ<DataVector, Dim, Frame::Inertial>*> flux2,
      const Scalar<DataVector>& var1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,
      const Scalar<DataVector>& prim_var);

  template <grid::Is grid_is>
  static void normal_dot_flux(
      const gsl::not_null<type*> normal_dot_flux,
      const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t,
      const double sign);

  template <grid::Is grid_is>
  static void dt_value(const gsl::not_null<type*> dt_var2,
                       const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                       const double t);
};

struct TempVar : db::SimpleTag {
  using type = Scalar<DataVector>;
  static void value(const gsl::not_null<type*> temp_var,
                    const Scalar<DataVector>& var1);
  template <size_t Dim>
  static void value(const gsl::not_null<type*> temp_var,
                    const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                    const double t);
};

struct PrimVar : db::SimpleTag {
  using type = Scalar<DataVector>;
  static void value(const gsl::not_null<type*> prim_var,
                    const Scalar<DataVector>& var1);
  template <size_t Dim>
  static void value(const gsl::not_null<type*> prim_var,
                    const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                    const double t);
};

struct Source1 : db::SimpleTag {
  using type = Scalar<DataVector>;
  template <size_t Dim>
  static void value(const gsl::not_null<type*> source1,
                    const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                    const double t);
};

template <size_t Dim>
struct Source2 : db::SimpleTag {
  using type = tnsr::I<DataVector, Dim, Frame::Inertial>;
  static void value(const gsl::not_null<type*> source2,
                    const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                    const double t);
};

template <typename System>
struct InitializeVars {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    static constexpr size_t volume_dim = System::volume_dim;
    using variables_tag = System::variables_tag;
    const auto& x =
        db::get<domain::Tags::Coordinates<volume_dim, Frame::Inertial>>(box);
    const double t = db::get<Tags::Time>(box);
    db::mutate<variables_tag, Source1, Source2<volume_dim>>(
        [&x, &t](
            const gsl::not_null<typename variables_tag::type*> evolved_vars,
            const gsl::not_null<Scalar<DataVector>*> source1,
            const gsl::not_null<
                tnsr::I<DataVector, volume_dim, Frame::Inertial>*>
                source2) {
          Var1::value(make_not_null(&get<Var1>(*evolved_vars)), x, t);
          Var2<volume_dim>::value(
              make_not_null(&get<Var2<volume_dim>>(*evolved_vars)), x, t);
          Source1::value(source1, x, t);
          Source2<volume_dim>::value(source2, x, t);
        },
        make_not_null(&box));
    if constexpr (System::has_primitive_and_conservative_vars) {
      db::mutate<PrimVar>(
          [&x, &t](const gsl::not_null<Scalar<DataVector>*> prim_var) {
            PrimVar::value(prim_var, x, t);
          },
          make_not_null(&box));
    }
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace TestHelpers::evolution::dg::Actions
