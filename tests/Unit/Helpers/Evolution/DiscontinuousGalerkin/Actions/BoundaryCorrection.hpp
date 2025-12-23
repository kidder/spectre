// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>

#include "Evolution/BoundaryCorrection.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
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
/// \endcond

template <size_t Dim>
struct BoundaryCorrection final : public ::evolution::BoundaryCorrection {
  struct MaxAbsCharSpeed : db::SimpleTag {
    using type = Scalar<DataVector>;
  };

  /// \cond
  explicit BoundaryCorrection(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(BoundaryCorrection);  // NOLINT
  /// \endcond

  BoundaryCorrection() = default;

  std::unique_ptr<::evolution::BoundaryCorrection> get_clone() const override {
    return std::make_unique<BoundaryCorrection>(*this);
  }

  using variables_tags = tmpl::list<Var1, Var2<Dim>>;
  using variables_tag = Tags::Variables<variables_tags>;
  using dg_package_field_tags = tmpl::push_back<
      tmpl::append<db::wrap_tags_in<::Tags::NormalDotFlux, variables_tags>,
                   variables_tags>,
      MaxAbsCharSpeed>;
  using dg_package_data_temporary_tags = tmpl::list<>;
  using dg_package_data_volume_tags = tmpl::list<>;
  using dg_boundary_terms_volume_tags = tmpl::list<>;

  double dg_package_data(
      const gsl::not_null<Scalar<DataVector>*> out_normal_dot_flux_var1,
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
          out_normal_dot_flux_var2,
      const gsl::not_null<Scalar<DataVector>*> out_var1,
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> out_var2,
      const gsl::not_null<Scalar<DataVector>*> max_abs_char_speed,

      const Scalar<DataVector>& var1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,

      const tnsr::I<DataVector, Dim, Frame::Inertial>& flux_var1,
      const tnsr::IJ<DataVector, Dim, Frame::Inertial>& flux_var2,

      const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
      const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
          mesh_velocity,
      const std::optional<Scalar<DataVector>>& normal_dot_mesh_velocity) const {
    return -999.0;
  }
};

template <size_t Dim>
PUP::able::PUP_ID BoundaryCorrection<Dim>::my_PUP_ID = 0;
}  // namespace TestHelpers::evolution::dg::Actions
