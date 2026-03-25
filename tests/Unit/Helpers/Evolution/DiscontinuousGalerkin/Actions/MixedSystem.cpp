// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/MixedSystem.hpp"

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/DiscontinuousGalerkin/TimeDerivativeDecisions.hpp"
#include "Utilities/Gsl.hpp"

namespace TestHelpers::evolution::dg::Actions {
template <size_t Dim>
MixedBoundaryCorrection<Dim>::MixedBoundaryCorrection(
    CkMigrateMessage* /*unused*/) {}

template <size_t Dim>
std::unique_ptr<::evolution::BoundaryCorrection>
MixedBoundaryCorrection<Dim>::get_clone() const {
  return std::make_unique<MixedBoundaryCorrection>(*this);
}

template <size_t Dim>
template <bool UseMovingMesh>
Variables<typename MixedBoundaryCorrection<Dim>::dg_package_field_tags>
MixedBoundaryCorrection<Dim>::expected_mortar_data(
    const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t) {
  Variables<dg_package_field_tags> result(get<0>(x).size());
  Var1::value(make_not_null(&get<Var1>(result)), x, t);
  Var2<Dim>::value(make_not_null(&get<Var2<Dim>>(result)), x, t);
  Var2<Dim>::template flux<UseMovingMesh>(
      make_not_null(
          &get<::Tags::Flux<Var2<Dim>, tmpl::size_t<Dim>, Frame::Inertial>>(
              result)),
      x, t);
  Var4::value(make_not_null(&get<Var4>(result)), x, t);
  return result;
}

template <size_t Dim>
double MixedBoundaryCorrection<Dim>::dg_package_data(
    const gsl::not_null<Scalar<DataVector>*> out_var1,
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> out_var2,
    const gsl::not_null<tnsr::IJ<DataVector, Dim, Frame::Inertial>*>
        out_flux_var2,
    const gsl::not_null<Scalar<DataVector>*> out_var4,

    const Scalar<DataVector>& var1,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,

    const tnsr::IJ<DataVector, Dim, Frame::Inertial>& flux_var2,

    const Scalar<DataVector>& var4,

    const tnsr::i<DataVector, Dim, Frame::Inertial>& /*normal_covector*/,
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
    /*mesh_velocity*/,
    const std::optional<Scalar<DataVector>>& /*normal_dot_mesh_velocity*/,
    const double /*time*/) {
  *out_var1 = var1;
  *out_var2 = var2;
  *out_flux_var2 = flux_var2;
  *out_var4 = var4;
  return std::numeric_limits<double>::signaling_NaN();
}

template <size_t Dim>
::evolution::dg::TimeDerivativeDecisions<Dim>
MixedTimeDerivativeTerms<Dim>::apply(
    // Time derivatives returned by reference. All the tags in the
    // variables_tag in the system struct.
    const gsl::not_null<Scalar<DataVector>*> dt_var1,
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> dt_var2,

    // Fluxes returned by reference. Listed in the system struct as
    // flux_variables.
    const gsl::not_null<tnsr::IJ<DataVector, Dim, Frame::Inertial>*> flux_var2,

    // Temporaries returned by reference. Listed in temporary_tags above.
    const gsl::not_null<Scalar<DataVector>*> var4,

    // Partial derivative arguments. Listed in the system struct as
    // gradient_variables.
    const tnsr::i<DataVector, Dim, Frame::Inertial>& d_var1,

    // Arguments listed in argument_tags above
    const Scalar<DataVector>& var1,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,
    const Scalar<DataVector>& source1,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& source2) {
  Var4::value(var4, var1);

  // Set source terms
  get(*dt_var1) = get(source1);
  for (size_t i = 0; i < Dim; ++i) {
    get(*dt_var1) -= 2.0 * var2.get(i) * d_var1.get(i) / get(var1);
    dt_var2->get(i) = source2.get(i);
  }

  // Set fluxes
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = 0; j < Dim; ++j) {
      flux_var2->get(i, j) = var2.get(i) * var2.get(j) / get(var1);
      if (i == j) {
        flux_var2->get(i, j) += get(*var4);
      }
    }
  }
  return {true};
}

template <size_t Dim>
::evolution::dg::TimeDerivativeDecisions<Dim>
MixedTimeDerivativeTermsWithVariables<Dim>::apply(
    // Time derivatives returned by reference. All the tags in the
    // variables_tag in the system struct.
    const gsl::not_null<Variables<tmpl::list<dt_var1, dt_var2>>*> dt_vars,

    // Fluxes returned by reference. Listed in the system struct as
    // flux_variables.
    const gsl::not_null<Variables<tmpl::list<flux_var2>>*> flux_vars,

    // Temporaries returned by reference. Listed in temporary_tags above.
    const gsl::not_null<Variables<tmpl::list<Var4>>*> temporaries,

    // Partial derivative arguments. Listed in the system struct as
    // gradient_variables.
    const tnsr::i<DataVector, Dim, Frame::Inertial>& d_var1,

    // Arguments listed in argument_tags above
    const Scalar<DataVector>& var1,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,
    const Scalar<DataVector>& source1,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& source2) {
  // just forward to other implementation to reduce code duplication
  base::apply(get<dt_var1>(dt_vars), get<dt_var2>(dt_vars),
              get<flux_var2>(flux_vars), get<Var4>(temporaries), d_var1, var1,
              var2, source1, source2);
  return {true};
}

template struct MixedBoundaryCorrection<1>;
template struct MixedBoundaryCorrection<2>;
template struct MixedBoundaryCorrection<3>;
template Variables<MixedBoundaryCorrection<1>::dg_package_field_tags>
    MixedBoundaryCorrection<1>::template expected_mortar_data<true>(
        const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t);
template Variables<MixedBoundaryCorrection<2>::dg_package_field_tags>
    MixedBoundaryCorrection<2>::template expected_mortar_data<true>(
        const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t);
template Variables<MixedBoundaryCorrection<3>::dg_package_field_tags>
    MixedBoundaryCorrection<3>::template expected_mortar_data<true>(
        const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t);
template Variables<MixedBoundaryCorrection<1>::dg_package_field_tags>
    MixedBoundaryCorrection<1>::template expected_mortar_data<false>(
        const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t);
template Variables<MixedBoundaryCorrection<2>::dg_package_field_tags>
    MixedBoundaryCorrection<2>::template expected_mortar_data<false>(
        const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t);
template Variables<MixedBoundaryCorrection<3>::dg_package_field_tags>
    MixedBoundaryCorrection<3>::template expected_mortar_data<false>(
        const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t);

template struct MixedTimeDerivativeTerms<1>;
template struct MixedTimeDerivativeTerms<2>;
template struct MixedTimeDerivativeTerms<3>;

template struct MixedTimeDerivativeTermsWithVariables<1>;
template struct MixedTimeDerivativeTermsWithVariables<2>;
template struct MixedTimeDerivativeTermsWithVariables<3>;

template struct MixedSystem<1, true>;
template struct MixedSystem<2, true>;
template struct MixedSystem<3, true>;
template struct MixedSystem<1, false>;
template struct MixedSystem<2, false>;
template struct MixedSystem<3, false>;
}  // namespace TestHelpers::evolution::dg::Actions
