// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Utilities/TMPL.hpp"

/// \cond
namespace domain::BoundaryConditons {
class BoundaryCondition;
}  // namespace domain::BoundaryConditons
namespace Tags {
template <typename TagsList>
struct Variables;
}  // namespace Tags
/// \endcond

namespace TestHelpers::evolution::dg::Actions {
/// \cond
struct Var1;
template <size_t Dim>
struct Var2;
struct Var3;
struct Var3Squared;
/// \endcond

template <size_t Dim>
struct ConservativeTimeDerivativeTerms {
  using argument_tags = tmpl::list<Var1, Var2<Dim>, Var3>;
  using temporary_tags = tmpl::list<Var3Squared>;

  /// [dt_con]
  static ::evolution::dg::TimeDerivativeDecisions<Dim> apply(
      // Time derivatives returned by reference. All the tags in the
      // variables_tag in the system struct.
      const gsl::not_null<Scalar<DataVector>*> dt_var1,
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> dt_var2,

      // Fluxes returned by reference. Listed in the system struct as
      // flux_variables.
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> flux_var1,
      const gsl::not_null<tnsr::IJ<DataVector, Dim, Frame::Inertial>*>
          flux_var2,

      // Temporaries returned by reference. Listed in temporary_tags above.
      const gsl::not_null<Scalar<DataVector>*> square_var3,

      // Arguments listed in argument_tags above
      const Scalar<DataVector>& var1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,
      const Scalar<DataVector>& var3) {
    get(*square_var3) = square(get(var3));

    // Set source terms
    get(*dt_var1) = get(*square_var3);
    for (size_t d = 0; d < Dim; ++d) {
      dt_var2->get(d) = get(var3) * d;
    }

    // Set fluxes
    for (size_t i = 0; i < Dim; ++i) {
      flux_var1->get(i) = square(get(var1)) * var2.get(i);
      for (size_t j = 0; j < Dim; ++j) {
        flux_var2->get(i, j) = var2.get(i) * var2.get(j) * get(var1);
        if (i == j) {
          flux_var2->get(i, j) += cube(get(var1));
        }
      }
    }
    return {true};
  }
  /// [dt_con]
};

template <size_t Dim>
struct ConservativeSystem {
  static constexpr size_t volume_dim = Dim;
  static constexpr bool has_primitive_and_conservative_vars = false;
  using boundary_conditions_base =
      domain::BoundaryConditions::BoundaryCondition;
  using compute_volume_time_derivative_terms =
      ConservativeTimeDerivativeTerms<Dim>;
  using flux_variables = tmpl::list<Var1, Var2<Dim>>;
  using gradient_variables = tmpl::list<>;
  using variables_tag = Tags::Variables<tmpl::list<Var1, Var2<Dim>>>;
};
}  // namespace TestHelpers::evolution::dg::Actions
