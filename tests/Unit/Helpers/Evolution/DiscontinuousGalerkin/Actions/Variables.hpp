// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Tags.hpp"
#include "Parallel/AlgorithmExecution.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/Gsl.hpp"

#include "Parallel/Printf/Printf.hpp"

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
static constexpr double c_s = 0.5;

template <size_t Dim>
constexpr std::array<double, Dim> k{};

template <size_t Dim>
DataVector u(const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
             const double t);

template <size_t Dim>
DataVector f(const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
             const double t);

template <size_t Dim>
DataVector df(const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
              const double t);
}  // namespace wave

struct Var1 : db::SimpleTag {
  using type = Scalar<DataVector>;

  template <size_t Dim>
  static void value(const gsl::not_null<type*> var1,
                    const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                    const double t);

  template <size_t Dim>
  static void flux(
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> flux,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t);

  template <size_t Dim>
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
      const gsl::not_null<tnsr::IJ<DataVector, Dim, Frame::Inertial>*> flux,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t);

  static void dt_value(const gsl::not_null<type*> dt_var2,
                       const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                       const double t);
};

template <size_t Dim>
struct Var3 : db::SimpleTag {
  using type = tnsr::I<DataVector, Dim, Frame::Inertial>;
  static void value(const gsl::not_null<type*> var3,
                    const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                    const double t);
};

struct Var4 : db::SimpleTag {
  using type = Scalar<DataVector>;
  static void value(const gsl::not_null<type*> var4,
                    const Scalar<DataVector>& var1);
  template <size_t Dim>
  static void value(const gsl::not_null<type*> var4,
                    const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                    const double t);
};

struct Source1: db::SimpleTag {
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
    db::mutate<variables_tag, Var3<volume_dim>, Source1, Source2<volume_dim>>(
        [&x, &t](
            const gsl::not_null<typename variables_tag::type*> evolved_vars,
            const gsl::not_null<
                tnsr::I<DataVector, volume_dim, Frame::Inertial>*>
                var3,
            const gsl::not_null<
                Scalar<DataVector>*>
                source1,
            const gsl::not_null<
                tnsr::I<DataVector, volume_dim, Frame::Inertial>*>
                source2) {
          Var1::value(make_not_null(&get<Var1>(*evolved_vars)), x, t);
          Var2<volume_dim>::value(
              make_not_null(&get<Var2<volume_dim>>(*evolved_vars)), x, t);
          Var3<volume_dim>::value(var3, x, t);
          Source1::value(source1, x, t);
          Source2<volume_dim>::value(source2, x, t);
        },
        make_not_null(&box));
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace TestHelpers::evolution::dg::Actions
