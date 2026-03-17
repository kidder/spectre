// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <limits>
#include <memory>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/TaggedContainers.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Evolution/BoundaryCorrection.hpp"
#include "Evolution/PassVariables.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/Variables.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace domain::BoundaryConditons {
class BoundaryCondition;
}  // namespace domain::BoundaryConditons
/// \endcond

namespace TestHelpers::evolution::dg::Actions {
template <size_t Dim>
struct ConservativeBoundaryCorrection final
    : public ::evolution::BoundaryCorrection {
  /// \cond
  explicit ConservativeBoundaryCorrection(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(ConservativeBoundaryCorrection);  // NOLINT
  /// \endcond

  ConservativeBoundaryCorrection() = default;

  std::unique_ptr<::evolution::BoundaryCorrection> get_clone() const override {
    return std::make_unique<ConservativeBoundaryCorrection>(*this);
  }

  using dg_package_field_tags = tmpl::list<
      Var1, Var2<Dim>, ::Tags::Flux<Var1, tmpl::size_t<Dim>, Frame::Inertial>,
      ::Tags::Flux<Var2<Dim>, tmpl::size_t<Dim>, Frame::Inertial>, Var4>;
  using dg_package_data_temporary_tags = tmpl::list<Var4>;
  using dg_package_data_volume_tags = tmpl::list<::Tags::Time>;
  using dg_boundary_terms_volume_tags = tmpl::list<>;

  static Variables<dg_package_field_tags> expected_mortar_data(
      const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t) {
    Variables<dg_package_field_tags> result(get<0>(x).size());
    Var1::value(make_not_null(&get<Var1>(result)), x, t);
    Var2<Dim>::value(make_not_null(&get<Var2<Dim>>(result)), x, t);
    Var1::flux(make_not_null(
                   &get<::Tags::Flux<Var1, tmpl::size_t<Dim>, Frame::Inertial>>(
                       result)),
               x, t);
    Var2<Dim>::flux(
        make_not_null(
            &get<::Tags::Flux<Var2<Dim>, tmpl::size_t<Dim>, Frame::Inertial>>(
                result)),
        x, t);
    Var4::value(make_not_null(&get<Var4>(result)), x, t);
    return result;
  }

  static double dg_package_data(
      const gsl::not_null<Scalar<DataVector>*> out_var1,
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> out_var2,
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
          out_flux_var1,
      const gsl::not_null<tnsr::IJ<DataVector, Dim, Frame::Inertial>*>
          out_flux_var2,
      const gsl::not_null<Scalar<DataVector>*> out_var4,

      const Scalar<DataVector>& var1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,

      const tnsr::I<DataVector, Dim, Frame::Inertial>& flux_var1,
      const tnsr::IJ<DataVector, Dim, Frame::Inertial>& flux_var2,

      const Scalar<DataVector>& var4,

      const tnsr::i<DataVector, Dim, Frame::Inertial>& /*normal_covector*/,
      const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
      /*mesh_velocity*/,
      const std::optional<Scalar<DataVector>>& /*normal_dot_mesh_velocity*/,
      const double /*time*/) {
    *out_var1 = var1;
    *out_var2 = var2;
    *out_flux_var1 = flux_var1;
    *out_flux_var2 = flux_var2;
    *out_var4 = var4;
    return std::numeric_limits<double>::signaling_NaN();
  }
};

template <size_t Dim>
PUP::able::PUP_ID ConservativeBoundaryCorrection<Dim>::my_PUP_ID = 0;

template <size_t Dim>
struct ConservativeTimeDerivativeTerms {
  using argument_tags = tmpl::list<Var1, Var2<Dim>, Source1, Source2<Dim>>;
  using temporary_tags = tmpl::list<Var4>;

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
      const gsl::not_null<Scalar<DataVector>*> var4,

      // Arguments listed in argument_tags above
      const Scalar<DataVector>& var1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,
      const Scalar<DataVector>& source1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& source2) {
    Var4::value(var4, var1);

    // Set source terms
    get(*dt_var1) = get(source1);
    for (size_t d = 0; d < Dim; ++d) {
      dt_var2->get(d) = source2.get(d);
    }

    // Set fluxes
    for (size_t i = 0; i < Dim; ++i) {
      flux_var1->get(i) = var2.get(i);
      for (size_t j = 0; j < Dim; ++j) {
        flux_var2->get(i, j) = var2.get(i) * var2.get(j) / get(var1);
        if (i == j) {
          flux_var2->get(i, j) += get(*var4);
        }
      }
    }
    return {true};
  }
  /// [dt_con]
};

template <size_t Dim>
struct ConservativeTimeDerivativeTermsWithVariables
    : public ::evolution::PassVariables,
      private ConservativeTimeDerivativeTerms<Dim> {
  using base = ConservativeTimeDerivativeTerms<Dim>;

  using temporary_tags = typename base::temporary_tags;
  using argument_tags = typename base::argument_tags;

  using dt_var1 = ::Tags::dt<Var1>;
  using dt_var2 = ::Tags::dt<Var2<Dim>>;
  using flux_var1 = ::Tags::Flux<Var1, tmpl::size_t<Dim>, Frame::Inertial>;
  using flux_var2 = ::Tags::Flux<Var2<Dim>, tmpl::size_t<Dim>, Frame::Inertial>;

  /// [dt_con_variables]
  static ::evolution::dg::TimeDerivativeDecisions<Dim> apply(
      // Time derivatives returned by reference. All the tags in the
      // variables_tag in the system struct.
      const gsl::not_null<Variables<tmpl::list<dt_var1, dt_var2>>*> dt_vars,

      // Fluxes returned by reference. Listed in the system struct as
      // flux_variables.
      const gsl::not_null<Variables<tmpl::list<flux_var1, flux_var2>>*>
          flux_vars,

      // Temporaries returned by reference. Listed in temporary_tags above.
      const gsl::not_null<Variables<tmpl::list<Var4>>*> temporaries,

      // Arguments listed in argument_tags above
      const Scalar<DataVector>& var1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,
      const Scalar<DataVector>& source1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& source2) {
    // just forward to other implementation to reduce code duplication
    base::apply(get<dt_var1>(dt_vars), get<dt_var2>(dt_vars),
                get<flux_var1>(flux_vars), get<flux_var2>(flux_vars),
                get<Var4>(temporaries), var1, var2, source1, source2);
    return {true};
  }
};

template <size_t Dim, bool PassVariables>
struct ConservativeSystem {
  static constexpr size_t volume_dim = Dim;
  static constexpr bool has_primitive_and_conservative_vars = false;
  using boundary_conditions_base =
      domain::BoundaryConditions::BoundaryCondition;
  using boundary_correction = ConservativeBoundaryCorrection<Dim>;
  using compute_volume_time_derivative_terms = tmpl::conditional_t<
      PassVariables,
      ConservativeTimeDerivativeTermsWithVariables<Dim>,
      ConservativeTimeDerivativeTerms<Dim>>;
  using flux_variables = tmpl::list<Var1, Var2<Dim>>;
  using gradient_variables = tmpl::list<>;
  using variables_tag = Tags::Variables<tmpl::list<Var1, Var2<Dim>>>;
};
}  // namespace TestHelpers::evolution::dg::Actions
