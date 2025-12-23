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
/// endcond

namespace TestHelpers::evolution::dg::Actions {
struct Var1 : db::SimpleTag {
  using type = Scalar<DataVector>;

  static void value(const gsl::not_null<type*> var1,
                    const DataVector& r_squared) {
    var1->get() = r_squared;
  }
  template <size_t Dim>
  static void dt_value(const gsl::not_null<type*> dt_var1,
                       const DataVector& r_squared,
                       const tnsr::I<DataVector, Dim, Frame::Inertial>& x) {
    dt_var1->get() = -6.0 * get<0>(x) * square(r_squared) -
                     12.0 * cube(get<0>(x)) * r_squared +
                     square(3.0 * r_squared - 1.0);
    if constexpr (Dim > 1) {
      dt_var1->get() -= 4.0 * get<1>(x) * square(r_squared) +
                        8.0 * cube(get<1>(x)) * r_squared;
    }
    if constexpr (Dim > 2) {
      dt_var1->get() -= 2.0 * get<2>(x) * square(r_squared) +
                        4.0 * cube(get<2>(x)) * r_squared;
    }
  }
};

template <size_t Dim>
struct Var2 : db::SimpleTag {
  using type = tnsr::I<DataVector, Dim, Frame::Inertial>;

  static void value(const gsl::not_null<type*> var2,
                    const tnsr::I<DataVector, Dim, Frame::Inertial>& x) {
    get<0>(*var2) = 3.0 * square(get<0>(x));
    if constexpr (Dim > 1) {
      get<1>(*var2) = 2.0 * square(get<1>(x));
    }
    if constexpr (Dim > 2) {
      get<2>(*var2) = square(get<2>(x));
    }
  }

  static void dt_value(const gsl::not_null<type*> dt_var2,
                       const DataVector& r_squared,
                       const tnsr::I<DataVector, Dim, Frame::Inertial>& x) {
    dt_var2->get(0) = -36.0 * r_squared * cube(get<0>(x)) -
                      18.0 * pow<5>(get<0>(x)) -
                      6.0 * square(r_squared) * get<0>(x);
    // if constexpr (Dim > 1) {
    //   dt_var2->get() -= 4.0 * get<1>(x) * square(r_squared) +
    //                     8.0 * cube(get<1>(x)) * r_squared;
    // }
    // if constexpr (Dim > 2) {
    //   dt_var2->get() -= 2.0 * get<2>(x) * square(r_squared) +
    //                     4.0 * cube(get<2>(x)) * r_squared;
    // }
  }
};

struct Var3 : db::SimpleTag {
  using type = Scalar<DataVector>;

  static void value(const gsl::not_null<type*> var3,
                    const DataVector& r_squared) {
    var3->get() = 3.0 * r_squared - 1.0;
  }
};

struct Var3Squared : db::SimpleTag {
  using type = Scalar<DataVector>;
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
    db::mutate<variables_tag, Var3>(
        [&x](const gsl::not_null<typename variables_tag::type*> evolved_vars,
             const gsl::not_null<Scalar<DataVector>*> var3) {
          const DataVector r_squared = get(dot_product(x, x));
          Var1::value(make_not_null(&get<Var1>(*evolved_vars)), r_squared);
          Var2<volume_dim>::value(
              make_not_null(&get<Var2<volume_dim>>(*evolved_vars)), x);
          Var3::value(var3, r_squared);
        },
        make_not_null(&box));
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};
}  // namespace TestHelpers::evolution::dg::Actions
