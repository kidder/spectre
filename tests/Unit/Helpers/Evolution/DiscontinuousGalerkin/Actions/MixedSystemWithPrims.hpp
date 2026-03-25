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
struct MixedWithPrimsBoundaryCorrection final
    : public ::evolution::BoundaryCorrection {
  /// \cond
  explicit MixedWithPrimsBoundaryCorrection(CkMigrateMessage* /*unused*/);
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(MixedWithPrimsBoundaryCorrection);  // NOLINT
  /// \endcond

  MixedWithPrimsBoundaryCorrection() = default;

  std::unique_ptr<::evolution::BoundaryCorrection> get_clone() const override;

  using dg_package_field_tags =
      tmpl::list<Var1, Var2<Dim>,
                 ::Tags::Flux<Var2<Dim>, tmpl::size_t<Dim>, Frame::Inertial>,
                 Var3>;
  using dg_package_data_primitive_tags = tmpl::list<Var4>;
  using dg_package_data_temporary_tags = tmpl::list<>;
  using dg_package_data_volume_tags = tmpl::list<::Tags::Time>;
  using dg_boundary_terms_volume_tags = tmpl::list<>;

  template <bool UseMovingMesh>
  static Variables<dg_package_field_tags> expected_mortar_data(
      const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t);

  static double dg_package_data(
      const gsl::not_null<Scalar<DataVector>*> out_var1,
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> out_var2,
      const gsl::not_null<tnsr::IJ<DataVector, Dim, Frame::Inertial>*>
          out_flux_var2,
      const gsl::not_null<Scalar<DataVector>*> out_var3,

      const Scalar<DataVector>& var1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,

      const tnsr::IJ<DataVector, Dim, Frame::Inertial>& flux_var2,

      const Scalar<DataVector>& var4,

      const tnsr::i<DataVector, Dim, Frame::Inertial>& /*normal_covector*/,
      const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
      /*mesh_velocity*/,
      const std::optional<Scalar<DataVector>>& /*normal_dot_mesh_velocity*/,
      const double /*time*/);
};

template <size_t Dim>
PUP::able::PUP_ID MixedWithPrimsBoundaryCorrection<Dim>::my_PUP_ID = 0;

template <size_t Dim>
struct MixedWithPrimsTimeDerivativeTerms {
  using argument_tags =
      tmpl::list<Var1, Var2<Dim>, Var4, Source1, Source2<Dim>>;
  using temporary_tags = tmpl::list<Var3>;

  static ::evolution::dg::TimeDerivativeDecisions<Dim> apply(
      // Time derivatives returned by reference. All the tags in the
      // variables_tag in the system struct.
      const gsl::not_null<Scalar<DataVector>*> dt_var1,
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> dt_var2,

      // Fluxes returned by reference. Listed in the system struct as
      // flux_variables.
      const gsl::not_null<tnsr::IJ<DataVector, Dim, Frame::Inertial>*>
          flux_var2,

      // Temporaries returned by reference. Listed in temporary_tags above.
      const gsl::not_null<Scalar<DataVector>*> var3,

      // Partial derivative arguments. Listed in the system struct as
      // gradient_variables.
      const tnsr::i<DataVector, Dim, Frame::Inertial>& d_var1,

      // Arguments listed in argument_tags above
      const Scalar<DataVector>& var1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,
      const Scalar<DataVector>& var4, const Scalar<DataVector>& source1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& source2);
};

template <size_t Dim>
struct MixedWithPrimsTimeDerivativeTermsWithVariables
    : public ::evolution::PassVariables,
      private MixedWithPrimsTimeDerivativeTerms<Dim> {
  using base = MixedWithPrimsTimeDerivativeTerms<Dim>;

  using temporary_tags = typename base::temporary_tags;
  using argument_tags = typename base::argument_tags;

  using dt_var1 = ::Tags::dt<Var1>;
  using dt_var2 = ::Tags::dt<Var2<Dim>>;
  using flux_var2 = ::Tags::Flux<Var2<Dim>, tmpl::size_t<Dim>, Frame::Inertial>;

  static ::evolution::dg::TimeDerivativeDecisions<Dim> apply(
      // Time derivatives returned by reference. All the tags in the
      // variables_tag in the system struct.
      const gsl::not_null<Variables<tmpl::list<dt_var1, dt_var2>>*> dt_vars,

      // Fluxes returned by reference. Listed in the system struct as
      // flux_variables.
      const gsl::not_null<Variables<tmpl::list<flux_var2>>*> flux_vars,

      // Temporaries returned by reference. Listed in temporary_tags above.
      const gsl::not_null<Variables<tmpl::list<Var3>>*> temporaries,

      // Partial derivative arguments. Listed in the system struct as
      // gradient_variables.
      const tnsr::i<DataVector, Dim, Frame::Inertial>& d_var1,

      // Arguments listed in argument_tags above
      const Scalar<DataVector>& var1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,
      const Scalar<DataVector>& var4, const Scalar<DataVector>& source1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& source2);
};

template <size_t Dim, bool PassVariables>
struct MixedSystemWithPrims {
  static constexpr size_t volume_dim = Dim;
  static constexpr bool has_primitive_and_conservative_vars = true;
  using boundary_conditions_base =
      domain::BoundaryConditions::BoundaryCondition;
  using boundary_correction = MixedWithPrimsBoundaryCorrection<Dim>;
  using compute_volume_time_derivative_terms =
      tmpl::conditional_t<PassVariables,
                          MixedWithPrimsTimeDerivativeTermsWithVariables<Dim>,
                          MixedWithPrimsTimeDerivativeTerms<Dim>>;
  using flux_variables = tmpl::list<Var2<Dim>>;
  using gradient_variables = tmpl::list<Var1>;
  using primitive_variables_tag = Tags::Variables<tmpl::list<Var4>>;
  using variables_tag = Tags::Variables<tmpl::list<Var1, Var2<Dim>>>;
};
}  // namespace TestHelpers::evolution::dg::Actions
