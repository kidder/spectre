// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Domain/BoundaryConditions/BoundaryCondition.hpp"
#include "Evolution/BoundaryConditions/Type.hpp"

namespace TestHelpers::evolution::dg::Actions {
template <size_t Dim>
struct BoundaryCondition final
    : public domain::BoundaryConditions::BoundaryCondition {
  // This determines the function called in BoundaryConditionsImpl.hpp
  // namely dg_demand_outgoing_char_speeds
  static constexpr ::evolution::BoundaryConditions::Type bc_type =
      ::evolution::BoundaryConditions::Type::DemandOutgoingCharSpeeds;

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
  using dg_interior_evolved_variables_tags = tmpl::list<>;
  using dg_interior_temporary_tags = tmpl::list<>;

  // Just count the number of times the boundary condition is called
  static std::optional<std::string> dg_demand_outgoing_char_speeds(
      const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
      /*face_mesh_velocity*/,
      const tnsr::i<DataVector, Dim, Frame::Inertial>&
      /*outward_directed_normal_covector*/) {
    ++number_of_times_called;
    return std::nullopt;
  }

  static size_t number_of_times_called;
};

template <size_t Dim>
PUP::able::PUP_ID BoundaryCondition<Dim>::my_PUP_ID = 0;

template <size_t Dim>
size_t BoundaryCondition<Dim>::number_of_times_called = 0;
}  // namespace TestHelpers::evolution::dg::Actions
