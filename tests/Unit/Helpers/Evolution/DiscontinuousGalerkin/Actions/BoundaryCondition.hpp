// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/BoundaryConditions/Type.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/Variables.hpp"

namespace TestHelpers::evolution::dg::Actions {
template <size_t Dim>
struct BoundaryCondition final
    : public domain::BoundaryConditions::BoundaryCondition {
  // This determines the function called in BoundaryConditionsImpl.hpp
  // namely dg_demand_outgoing_char_speeds
  static constexpr ::evolution::BoundaryConditions::Type bc_type =
      ::evolution::BoundaryConditions::Type::Ghost;

  /// \cond
  explicit BoundaryCondition(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(BoundaryCondition);  // NOLINT
  /// \endcond

  BoundaryCondition() = default;

  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition> get_clone()
      const override {
    return std::make_unique<BoundaryCondition>(*this);
  }

  using dg_gridless_tags = tmpl::list<>;
  using dg_interior_evolved_variables_tags = tmpl::list<Var1, Var2<Dim>>;
  using dg_interior_primitive_variables_tags = tmpl::list<>;
  using dg_interior_temporary_tags = tmpl::list<TempVar>;

  // Just count the number of times the boundary condition is called
  static std::optional<std::string> dg_ghost(
      const gsl::not_null<Scalar<DataVector>*> out_var1,
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> out_var2,
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
          out_flux_var1,
      const gsl::not_null<tnsr::IJ<DataVector, Dim, Frame::Inertial>*>
          out_flux_var2,
      const gsl::not_null<Scalar<DataVector>*> out_temp_var,
      const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
          face_mesh_velocity,
      const tnsr::i<DataVector, Dim, Frame::Inertial>&
      /*outward_directed_normal_covector*/,
      const Scalar<DataVector>& var1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,
      const Scalar<DataVector>& temp_var) {
    ++number_of_times_called;
    *out_var1 = var1;
    *out_var2 = var2;
    *out_temp_var = temp_var;
    for (size_t i = 0; i < Dim; ++i) {
      out_flux_var1->get(i) = var2.get(i);
      // if(face_mesh_velocity.has_value()) {
      //   out_flux_var1->get(i) -= face_mesh_velocity.value().get(i) *
      //   get(var1);
      // }
      for (size_t j = 0; j < Dim; ++j) {
        out_flux_var2->get(i, j) = var2.get(i) * var2.get(j) / get(var1);
        // if (face_mesh_velocity.has_value()) {
        //   out_flux_var2->get(i, j) -=
        //       face_mesh_velocity.value().get(i) * var2.get(j);
        // }
        if (i == j) {
          out_flux_var2->get(i, j) += get(temp_var);
        }
      }
    }
    return std::nullopt;
  }

  static size_t number_of_times_called;
};

template <size_t Dim>
PUP::able::PUP_ID BoundaryCondition<Dim>::my_PUP_ID = 0;

template <size_t Dim>
size_t BoundaryCondition<Dim>::number_of_times_called = 0;
}  // namespace TestHelpers::evolution::dg::Actions
