// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/ConservativeSystem.hpp"

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/DiscontinuousGalerkin/TimeDerivativeDecisions.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "Utilities/Gsl.hpp"

#include "Parallel/Printf/Printf.hpp"

namespace TestHelpers::evolution::dg::Actions {
template <size_t Dim>
ConservativeBoundaryCorrection<Dim>::ConservativeBoundaryCorrection(
    CkMigrateMessage* /*unused*/) {}

template <size_t Dim>
std::unique_ptr<::evolution::BoundaryCorrection>
ConservativeBoundaryCorrection<Dim>::get_clone() const {
  return std::make_unique<ConservativeBoundaryCorrection>(*this);
}

template <size_t Dim>
template <grid::Is grid_is>
Variables<typename ConservativeBoundaryCorrection<Dim>::dg_package_field_tags>
ConservativeBoundaryCorrection<Dim>::expected_mortar_data(
    const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t,
    const double sign) {
  Variables<dg_package_field_tags> result(get<0>(x).size());
  Var1::value(make_not_null(&get<Var1>(result)), x, t);
  Var2<Dim>::value(make_not_null(&get<Var2<Dim>>(result)), x, t);
  Var1::normal_dot_flux<grid_is>(
      make_not_null(&get<::Tags::NormalDotFlux<Var1>>(result)), normal_covector,
      x, t, sign);
  Var2<Dim>::template normal_dot_flux<grid_is>(
      make_not_null(&get<::Tags::NormalDotFlux<Var2<Dim>>>(result)),
      normal_covector, x, t, sign);
  TempVar::value(make_not_null(&get<TempVar>(result)), x, t);
  // Parallel::printf("---- Expected--------------------------------------\n");
  // Parallel::printf("Result = %s\n", result);
  // Parallel::printf("n = %s\n", normal_covector);
  // Parallel::printf("x = %s\n", x);
  // Parallel::printf("t = %f\n", t);
  // Parallel::printf("--------------------------------------------------\n");
  return result;
}

template <size_t Dim>
double ConservativeBoundaryCorrection<Dim>::dg_package_data(
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
    const std::optional<Scalar<DataVector>>& /*normal_dot_mesh_velocity*/,
    const double /*time*/) {
  *out_normal_dot_flux_var1 = dot_product(flux_var1, normal_covector);
  for (size_t i = 0; i < Dim; ++i) {
    out_normal_dot_flux_var2->get(i) =
        flux_var2.get(0, i) * normal_covector.get(0);
    for (size_t j = 1; j < Dim; ++j) {
      out_normal_dot_flux_var2->get(i) +=
          flux_var2.get(j, i) * normal_covector.get(j);
    }
  }
  *out_var1 = var1;
  *out_var2 = var2;
  *out_temp_var = temp_var;
  // Parallel::printf("---- Packaged--------------------------------------\n");
  // Parallel::printf("flux1 = %s\n", flux_var1);
  // Parallel::printf("flux2 = %s\n", flux_var2);
  // Parallel::printf("out_n_dot_flux1 = %s\n", *out_normal_dot_flux_var1);
  // Parallel::printf("out_n_dot_flux2 = %s\n", *out_normal_dot_flux_var2);
  // Parallel::printf("var1 = %s\n", var1);
  // Parallel::printf("var2 = %s\n", var2);
  // Parallel::printf("temp = %s\n", temp_var);
  // Parallel::printf("n_dot_v = %s\n", normal_dot_mesh_velocity);
  // Parallel::printf("v = %s\n", mesh_velocity);
  // Parallel::printf("t = %f\n", time);
  // Parallel::printf("---------------------------------------------------\n");
  return std::numeric_limits<double>::signaling_NaN();
}

template <size_t Dim>
void ConservativeBoundaryCorrection<Dim>::dg_boundary_terms(
    const gsl::not_null<Scalar<DataVector>*> boundary_correction_var1,
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
        boundary_correction_var2,
    const Scalar<DataVector>& interior_normal_dot_flux_var1,
    const tnsr::I<DataVector, Dim, Frame::Inertial>&
        interior_normal_dot_flux_var2,
    const Scalar<DataVector>& /* interior_var1 */,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& /* interior_var2 */,
    const Scalar<DataVector>& /* interior_max_abs_char_speed */,
    const Scalar<DataVector>& exterior_normal_dot_flux_var1,
    const tnsr::I<DataVector, Dim, Frame::Inertial>&
        exterior_normal_dot_flux_var2,
    const Scalar<DataVector>& /* exterior_var1 */,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& /* exterior_var2 */,
    const Scalar<DataVector>& /* exterior_max_abs_char_speed */,
    const ::dg::Formulation dg_formulation) {
  get(*boundary_correction_var1) =
      -0.5 *
      ((dg_formulation == ::dg::Formulation::StrongInertial ? 1.0 : -1.0) *
           get(interior_normal_dot_flux_var1) +
       get(exterior_normal_dot_flux_var1));
  for (size_t i = 0; i < Dim; ++i) {
    boundary_correction_var2->get(i) =
        -0.5 *
        ((dg_formulation == ::dg::Formulation::StrongInertial ? 1.0 : -1.0) *
             interior_normal_dot_flux_var2.get(i) +
         exterior_normal_dot_flux_var2.get(i));
  }
}

template <size_t Dim>
::evolution::dg::TimeDerivativeDecisions<Dim>
ConservativeTimeDerivativeTerms<Dim>::apply(
    // Time derivatives returned by reference. All the tags in the
    // variables_tag in the system struct.
    const gsl::not_null<Scalar<DataVector>*> dt_var1,
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> dt_var2,

    // Fluxes returned by reference. Listed in the system struct as
    // flux_variables.
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> flux_var1,
    const gsl::not_null<tnsr::IJ<DataVector, Dim, Frame::Inertial>*> flux_var2,

    // Temporaries returned by reference. Listed in temporary_tags above.
    const gsl::not_null<Scalar<DataVector>*> temp_var,

    // Arguments listed in argument_tags above
    const Scalar<DataVector>& var1,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,
    const Scalar<DataVector>& source1,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& source2) {
  TempVar::value(temp_var, var1);

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
        flux_var2->get(i, j) += get(*temp_var);
      }
    }
  }
  return {true};
}

template <size_t Dim>
::evolution::dg::TimeDerivativeDecisions<Dim>
ConservativeTimeDerivativeTermsWithVariables<Dim>::apply(
    // Time derivatives returned by reference. All the tags in the
    // variables_tag in the system struct.
    const gsl::not_null<Variables<tmpl::list<dt_var1, dt_var2>>*> dt_vars,

    // Fluxes returned by reference. Listed in the system struct as
    // flux_variables.
    const gsl::not_null<Variables<tmpl::list<flux_var1, flux_var2>>*> flux_vars,

    // Temporaries returned by reference. Listed in temporary_tags above.
    const gsl::not_null<Variables<tmpl::list<TempVar>>*> temporaries,

    // Arguments listed in argument_tags above
    const Scalar<DataVector>& var1,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,
    const Scalar<DataVector>& source1,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& source2) {
  // just forward to other implementation to reduce code duplication
  base::apply(get<dt_var1>(dt_vars), get<dt_var2>(dt_vars),
              get<flux_var1>(flux_vars), get<flux_var2>(flux_vars),
              get<TempVar>(temporaries), var1, var2, source1, source2);
  return {true};
}

template struct ConservativeBoundaryCorrection<1>;
template struct ConservativeBoundaryCorrection<2>;
template struct ConservativeBoundaryCorrection<3>;

template Variables<ConservativeBoundaryCorrection<1>::dg_package_field_tags>
    ConservativeBoundaryCorrection<1>::template expected_mortar_data<
        grid::Is::Stationary>(
        const tnsr::i<DataVector, 1, Frame::Inertial>& normal_covector,
        const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t,
        const double sign);
template Variables<ConservativeBoundaryCorrection<2>::dg_package_field_tags>
    ConservativeBoundaryCorrection<2>::template expected_mortar_data<
        grid::Is::Stationary>(
        const tnsr::i<DataVector, 2, Frame::Inertial>& normal_covector,
        const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t,
        const double sign);
template Variables<ConservativeBoundaryCorrection<3>::dg_package_field_tags>
    ConservativeBoundaryCorrection<3>::template expected_mortar_data<
        grid::Is::Stationary>(
        const tnsr::i<DataVector, 3, Frame::Inertial>& normal_covector,
        const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t,
        const double sign);
template Variables<ConservativeBoundaryCorrection<1>::dg_package_field_tags>
    ConservativeBoundaryCorrection<1>::template expected_mortar_data<
        grid::Is::Comoving>(
        const tnsr::i<DataVector, 1, Frame::Inertial>& normal_covector,
        const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t,
        const double sign);
template Variables<ConservativeBoundaryCorrection<2>::dg_package_field_tags>
    ConservativeBoundaryCorrection<2>::template expected_mortar_data<
        grid::Is::Comoving>(
        const tnsr::i<DataVector, 2, Frame::Inertial>& normal_covector,
        const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t,
        const double sign);
template Variables<ConservativeBoundaryCorrection<3>::dg_package_field_tags>
    ConservativeBoundaryCorrection<3>::template expected_mortar_data<
        grid::Is::Comoving>(
        const tnsr::i<DataVector, 3, Frame::Inertial>& normal_covector,
        const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t,
        const double sign);
template Variables<ConservativeBoundaryCorrection<1>::dg_package_field_tags>
    ConservativeBoundaryCorrection<1>::template expected_mortar_data<
        grid::Is::Expanding>(
        const tnsr::i<DataVector, 1, Frame::Inertial>& normal_covector,
        const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t,
        const double sign);
template Variables<ConservativeBoundaryCorrection<2>::dg_package_field_tags>
    ConservativeBoundaryCorrection<2>::template expected_mortar_data<
        grid::Is::Expanding>(
        const tnsr::i<DataVector, 2, Frame::Inertial>& normal_covector,
        const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t,
        const double sign);
template Variables<ConservativeBoundaryCorrection<3>::dg_package_field_tags>
    ConservativeBoundaryCorrection<3>::template expected_mortar_data<
        grid::Is::Expanding>(
        const tnsr::i<DataVector, 3, Frame::Inertial>& normal_covector,
        const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t,
        const double sign);

template struct ConservativeTimeDerivativeTerms<1>;
template struct ConservativeTimeDerivativeTerms<2>;
template struct ConservativeTimeDerivativeTerms<3>;

template struct ConservativeTimeDerivativeTermsWithVariables<1>;
template struct ConservativeTimeDerivativeTermsWithVariables<2>;
template struct ConservativeTimeDerivativeTermsWithVariables<3>;

template struct ConservativeSystem<1, true>;
template struct ConservativeSystem<2, true>;
template struct ConservativeSystem<3, true>;
template struct ConservativeSystem<1, false>;
template struct ConservativeSystem<2, false>;
template struct ConservativeSystem<3, false>;
}  // namespace TestHelpers::evolution::dg::Actions
