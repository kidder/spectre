// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <limits>
#include <memory>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/TaggedContainers.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Evolution/BoundaryCorrection.hpp"
#include "Evolution/PassVariables.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/Variables.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
namespace dg {
enum class Formulation;
}  // namespace dg
namespace domain::BoundaryConditions {
class BoundaryCondition;
}  // namespace domain::BoundaryConditions
namespace evolution::dg {
template <size_t>
class TimeDerivativeDecisions;
}  // namespace evolution::dg
namespace gsl {
template <typename T>
class not_null;
}  // namespace gsl
template <typename>
class Variables;
/// \endcond

namespace TestHelpers::evolution::dg::Actions {
template <size_t Dim>
struct ConservativeBoundaryCorrection final
    : public ::evolution::BoundaryCorrection {
  /// \cond
  explicit ConservativeBoundaryCorrection(CkMigrateMessage* /*unused*/);
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(ConservativeBoundaryCorrection);  // NOLINT
  /// \endcond

  ConservativeBoundaryCorrection() = default;

  std::unique_ptr<::evolution::BoundaryCorrection> get_clone() const override;

  using variables_tags = tmpl::list<Var1, Var2<Dim>>;
  using dg_package_field_tags = tmpl::push_back<
      tmpl::append<db::wrap_tags_in<::Tags::NormalDotFlux, variables_tags>,
                   variables_tags>,
      TempVar>;
  using dg_package_data_temporary_tags = tmpl::list<TempVar>;
  using dg_package_data_volume_tags = tmpl::list<::Tags::Time>;
  using dg_boundary_terms_volume_tags = tmpl::list<>;

  template <grid::Is grid_is>
  static Variables<dg_package_field_tags> expected_mortar_data(
      const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t,
      const double sign);

  static double dg_package_data(
      const gsl::not_null<Scalar<DataVector>*> out_normal_dot_flux_var1,
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
          out_normal_dot_flux_var2,
      const gsl::not_null<Scalar<DataVector>*> out_var1,
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> out_var2,
      const gsl::not_null<Scalar<DataVector>*> out_temp_var,

      const Scalar<DataVector>& var1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,

      const tnsr::I<DataVector, Dim, Frame::Inertial>& flux_var1,
      const tnsr::IJ<DataVector, Dim, Frame::Inertial>& flux_var2,

      const Scalar<DataVector>& temp_var,

      const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
      const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
      /*mesh_velocity*/,
      const std::optional<Scalar<DataVector>>& normal_dot_mesh_velocity,
      const double /*time*/);

  static void dg_boundary_terms(
      const gsl::not_null<Scalar<DataVector>*> boundary_correction_var1,
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
          boundary_correction_var2,
      const Scalar<DataVector>& interior_normal_dot_flux_var1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>&
          interior_normal_dot_flux_var2,
      const Scalar<DataVector>& interior_var1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& interior_var2,
      const Scalar<DataVector>& interior_max_abs_char_speed,
      const Scalar<DataVector>& exterior_normal_dot_flux_var1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>&
          exterior_normal_dot_flux_var2,
      const Scalar<DataVector>& exterior_var1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& exterior_var2,
      const Scalar<DataVector>& exterior_max_abs_char_speed,
      const ::dg::Formulation dg_formulation);
};

template <size_t Dim>
PUP::able::PUP_ID ConservativeBoundaryCorrection<Dim>::my_PUP_ID = 0;

template <size_t Dim>
struct ConservativeTimeDerivativeTerms {
  using argument_tags = tmpl::list<Var1, Var2<Dim>, Source1, Source2<Dim>>;
  using temporary_tags = tmpl::list<TempVar>;

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
      const gsl::not_null<Scalar<DataVector>*> temp_var,

      // Arguments listed in argument_tags above
      const Scalar<DataVector>& var1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,
      const Scalar<DataVector>& source1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& source2);
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

  static ::evolution::dg::TimeDerivativeDecisions<Dim> apply(
      // Time derivatives returned by reference. All the tags in the
      // variables_tag in the system struct.
      const gsl::not_null<Variables<tmpl::list<dt_var1, dt_var2>>*> dt_vars,

      // Fluxes returned by reference. Listed in the system struct as
      // flux_variables.
      const gsl::not_null<Variables<tmpl::list<flux_var1, flux_var2>>*>
          flux_vars,

      // Temporaries returned by reference. Listed in temporary_tags above.
      const gsl::not_null<Variables<tmpl::list<TempVar>>*> temporaries,

      // Arguments listed in argument_tags above
      const Scalar<DataVector>& var1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,
      const Scalar<DataVector>& source1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& source2);
};

template <size_t Dim, bool PassVariables>
struct ConservativeSystem {
  static constexpr size_t volume_dim = Dim;
  static constexpr bool has_primitive_and_conservative_vars = false;
  using boundary_conditions_base =
      domain::BoundaryConditions::BoundaryCondition;
  using boundary_correction = ConservativeBoundaryCorrection<Dim>;
  using compute_volume_time_derivative_terms =
      tmpl::conditional_t<PassVariables,
                          ConservativeTimeDerivativeTermsWithVariables<Dim>,
                          ConservativeTimeDerivativeTerms<Dim>>;
  using flux_variables = tmpl::list<Var1, Var2<Dim>>;
  using gradient_variables = tmpl::list<>;
  using variables_tag = Tags::Variables<tmpl::list<Var1, Var2<Dim>>>;
};
}  // namespace TestHelpers::evolution::dg::Actions
