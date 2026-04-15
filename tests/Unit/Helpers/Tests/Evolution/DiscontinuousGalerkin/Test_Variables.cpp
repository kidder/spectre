// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <memory>
#include <random>
#include <tuple>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/Creators/TimeDependence/CubicScale.hpp"
#include "Domain/Creators/TimeDependence/None.hpp"
#include "Domain/Creators/TimeDependence/TimeDependence.hpp"
#include "Domain/Creators/TimeDependence/UniformTranslation.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeRandomVectorInMagnitudeRange.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/DomainCreator.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/Variables.hpp"
#include "NumericalAlgorithms/LinearOperators/Divergence.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Utilities/Gsl.hpp"

namespace TestHelpers::evolution::dg::Actions {
namespace {
template <grid::Is grid_is, size_t Dim>
void check_vars(const gsl::not_null<std::mt19937*> generator,
                const Mesh<Dim>& mesh,
                const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                                      Frame::Inertial>& inv_jac,
                const tnsr::I<DataVector, Dim, Frame::Inertial>& x_i,
                const double t,
                const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
                    grid_velocity) {
  const size_t n_pts = get<0>(x_i).size();

  Variables<tmpl::list<Var1, Var2<Dim>>> vars{n_pts};

  auto& var1 = get<Var1>(vars);
  Var1::value(make_not_null(&var1), x_i, t);

  auto& var2 = get<Var2<Dim>>(vars);
  Var2<Dim>::value(make_not_null(&var2), x_i, t);

  Scalar<DataVector> temp_var{n_pts};
  TempVar::value(make_not_null(&temp_var), x_i, t);

  Scalar<DataVector> prim_var{n_pts};
  PrimVar::value(make_not_null(&prim_var), x_i, t);

  using flux_var1 = ::Tags::Flux<Var1, tmpl::size_t<Dim>, Frame::Inertial>;
  using flux_var2 = ::Tags::Flux<Var2<Dim>, tmpl::size_t<Dim>, Frame::Inertial>;
  using flux_tags = tmpl::list<flux_var1, flux_var2>;
  Variables<flux_tags> fluxes(n_pts);

  auto& flux1 = get<flux_var1>(fluxes);
  Var1::flux(make_not_null(&flux1), var2);

  auto& flux2 = get<flux_var2>(fluxes);
  Var2<Dim>::flux(&flux2, var1, var2, prim_var);

  Scalar<DataVector> source1{n_pts};
  Source1::value(make_not_null(&source1), x_i, t);

  tnsr::I<DataVector, Dim, Frame::Inertial> source2{n_pts};
  Source2<Dim>::value(make_not_null(&source2), x_i, t);

  const auto div_fluxes = divergence(fluxes, mesh, inv_jac);

  Scalar<DataVector> dt_var1{n_pts};
  Var1::dt_value<grid_is>(make_not_null(&dt_var1), x_i, t);

  tnsr::I<DataVector, Dim, Frame::Inertial> dt_var2{n_pts};
  Var2<Dim>::template dt_value<grid_is>(make_not_null(&dt_var2), x_i, t);

  const auto n =
      make_random_vector_in_magnitude_range_flat<DataVector, Dim, UpLo::Lo,
                                                 Frame::Inertial>(
          generator, x_i, 1.0, 1.0);

  Scalar<DataVector> n_dot_flux1{n_pts};
  Var1::normal_dot_flux<grid_is>(make_not_null(&n_dot_flux1), n, x_i, t, 1.0);

  tnsr::I<DataVector, Dim, Frame::Inertial> n_dot_flux2{n_pts};
  Var2<Dim>::template normal_dot_flux<grid_is>(make_not_null(&n_dot_flux2), n,
                                               x_i, t, 1.0);

  if constexpr (grid_is == grid::Is::Stationary) {
    DataVector u{n_pts, -wave::c_s * t};
    for (size_t d = 0; d < Dim; ++d) {
      u += gsl::at(wave::k<Dim>(), d) * x_i.get(d);
    }
    CHECK_ITERABLE_APPROX(u, wave::u(x_i, t));

    const DataVector f = 0.01 * square(u) + 1.0;
    CHECK_ITERABLE_APPROX(f, wave::f(x_i, t));

    CHECK_ITERABLE_APPROX(get(var1), f);

    tnsr::I<DataVector, Dim, Frame::Inertial> expected_var2{n_pts};
    for (size_t i = 0; i < Dim; ++i) {
      expected_var2.get(i) = gsl::at(wave::k<Dim>(), i) * square(f);
    }
    CHECK_ITERABLE_APPROX(var2, expected_var2);

    const DataVector expected_temp_var = 0.5 * square(wave::c_s) * get(var1);
    CHECK_ITERABLE_APPROX(get(temp_var), expected_temp_var);
    TempVar::value(make_not_null(&temp_var), var1);
    CHECK_ITERABLE_APPROX(get(temp_var), expected_temp_var);

    const DataVector expected_prim_var = 0.5 * square(wave::c_s) * get(var1);
    CHECK_ITERABLE_APPROX(get(prim_var), expected_prim_var);
    PrimVar::value(make_not_null(&prim_var), var1);
    CHECK_ITERABLE_APPROX(get(prim_var), expected_prim_var);

    tnsr::I<DataVector, Dim, Frame::Inertial> expected_flux1{n_pts};
    for (size_t i = 0; i < Dim; ++i) {
      expected_flux1.get(i) = gsl::at(wave::k<Dim>(), i) * square(f);
    }
    CHECK_ITERABLE_APPROX(flux1, expected_flux1);

    tnsr::IJ<DataVector, Dim, Frame::Inertial> expected_flux2{n_pts};
    for (size_t i = 0; i < Dim; ++i) {
      for (size_t j = 0; j < Dim; ++j) {
        expected_flux2.get(i, j) =
            gsl::at(wave::k<Dim>(), i) * gsl::at(wave::k<Dim>(), j) * cube(f);
        if (i == j) {
          expected_flux2.get(i, j) += 0.5 * square(wave::c_s) * f;
        }
      }
    }
    CHECK_ITERABLE_APPROX(flux2, expected_flux2);

    const DataVector expected_source1 =
        wave::df(x_i, t) * (2.0 * f - wave::c_s);
    CHECK_ITERABLE_APPROX(get(source1), expected_source1);

    tnsr::I<DataVector, Dim, Frame::Inertial> expected_source2{n_pts};
    for (size_t i = 0; i < Dim; ++i) {
      expected_source2.get(i) =
          gsl::at(wave::k<Dim>(), i) * wave::df(x_i, t) *
          (3.0 * square(f) - 2.0 * wave::c_s * f + 0.5 * square(wave::c_s));
    }
    CHECK_ITERABLE_APPROX(source2, expected_source2);

    DataVector expected_dt_var1 = get(source1);
    expected_dt_var1 -= get(get<Tags::div<flux_var1>>(div_fluxes));
    const auto approx = Approx::custom().scale(1.0).epsilon(1.e-11);
    CHECK_ITERABLE_CUSTOM_APPROX(get(dt_var1), expected_dt_var1, approx);

    tnsr::I<DataVector, Dim, Frame::Inertial> expected_dt_var2{n_pts};
    const auto& div_flux2 = get<Tags::div<flux_var2>>(div_fluxes);
    for (size_t i = 0; i < Dim; ++i) {
      expected_dt_var2.get(i) = source2.get(i) - div_flux2.get(i);
    }
    CHECK_ITERABLE_CUSTOM_APPROX(dt_var2, expected_dt_var2, approx);

    DataVector expected_n_dot_flux1{n_pts, 0.0};
    for (size_t d = 0; d < Dim; ++d) {
      expected_n_dot_flux1 += flux1.get(d) * n.get(d);
    }
    CHECK_ITERABLE_APPROX(get(n_dot_flux1), expected_n_dot_flux1);

    tnsr::I<DataVector, Dim, Frame::Inertial> expected_n_dot_flux2{n_pts, 0.0};
    for (size_t i = 0; i < Dim; ++i) {
      for (size_t d = 0; d < Dim; ++d) {
        expected_n_dot_flux2.get(i) += flux2.get(d, i) * n.get(d);
      }
    }
    CHECK_ITERABLE_CUSTOM_APPROX(n_dot_flux2, expected_n_dot_flux2, approx);
  } else {
    using var_tags = tmpl::list<Var1, Var2<Dim>>;
    using d_var1 = ::Tags::deriv<Var1, tmpl::size_t<Dim>, Frame::Inertial>;
    using d_var2 = ::Tags::deriv<Var2<Dim>, tmpl::size_t<Dim>, Frame::Inertial>;
    const auto deriv_vars = partial_derivatives<var_tags>(vars, mesh, inv_jac);
    DataVector expected_dt_var1 = get(source1);
    expected_dt_var1 -= get(get<Tags::div<flux_var1>>(div_fluxes));
    for (size_t d = 0; d < Dim; ++d) {
      expected_dt_var1 +=
          grid_velocity.value().get(d) * get<d_var1>(deriv_vars).get(d);
    }
    const auto approx = Approx::custom().scale(1.0).epsilon(1.e-11);
    CHECK_ITERABLE_CUSTOM_APPROX(get(dt_var1), expected_dt_var1, approx);

    tnsr::I<DataVector, Dim, Frame::Inertial> expected_dt_var2{n_pts};
    const auto& div_flux2 = get<Tags::div<flux_var2>>(div_fluxes);
    for (size_t i = 0; i < Dim; ++i) {
      expected_dt_var2.get(i) = source2.get(i) - div_flux2.get(i);
      for (size_t d = 0; d < Dim; ++d) {
        expected_dt_var2.get(i) +=
            grid_velocity.value().get(d) * get<d_var2>(deriv_vars).get(d, i);
      }
    }
    CHECK_ITERABLE_CUSTOM_APPROX(dt_var2, expected_dt_var2, approx);

    const auto n_dot_v = dot_product(n, grid_velocity.value());
    DataVector expected_n_dot_flux1{n_pts, 0.0};
    for (size_t d = 0; d < Dim; ++d) {
      expected_n_dot_flux1 += flux1.get(d) * n.get(d);
    }
    expected_n_dot_flux1 -= get(n_dot_v) * get(var1);
    CHECK_ITERABLE_APPROX(get(n_dot_flux1), expected_n_dot_flux1);

    tnsr::I<DataVector, Dim, Frame::Inertial> expected_n_dot_flux2{n_pts, 0.0};
    for (size_t i = 0; i < Dim; ++i) {
      for (size_t d = 0; d < Dim; ++d) {
        expected_n_dot_flux2.get(i) += flux2.get(d, i) * n.get(d);
      }
      expected_n_dot_flux2.get(i) -= get(n_dot_v) * var2.get(i);
    }
    CHECK_ITERABLE_CUSTOM_APPROX(n_dot_flux2, expected_n_dot_flux2, approx);
  }
  if constexpr (grid_is == grid::Is::Comoving) {
    Scalar<DataVector> expected_dt_var1{n_pts, 0.0};
    tnsr::I<DataVector, Dim, Frame::Inertial> expected_dt_var2{n_pts, 0.0};
    const auto approx = Approx::custom().scale(1.0).epsilon(1.e-11);
    CHECK_ITERABLE_CUSTOM_APPROX(dt_var1, expected_dt_var1, approx);
    CHECK_ITERABLE_CUSTOM_APPROX(dt_var2, expected_dt_var2, approx);
  }
}

template <grid::Is grid_is, size_t Dim>
void test(const gsl::not_null<std::mt19937*> generator) {
  CAPTURE(grid_is);
  CAPTURE(Dim);
  const double t = 1.5;
  std::unique_ptr<domain::creators::time_dependence::TimeDependence<Dim>>
      time_dependence{nullptr};
  if constexpr (grid_is == grid::Is::Stationary) {
    time_dependence =
        std::make_unique<domain::creators::time_dependence::None<Dim>>();
  }
  if constexpr (grid_is == grid::Is::Expanding) {
    time_dependence =
        std::make_unique<domain::creators::time_dependence::CubicScale<Dim>>(
            0.0, 1000.0, true, std::array{grid::a_0, grid::a_0},
            std::array{grid::a_dot, grid::a_dot}, std::array{0.0, 0.0});
  }
  if constexpr (grid_is == grid::Is::Comoving) {
    time_dependence = std::make_unique<
        domain::creators::time_dependence::UniformTranslation<Dim>>(
        0.0, wave::comoving_v<Dim>());
  }
  const auto creator = domain_creator<true>(*time_dependence);
  const auto domain = creator->create_domain();
  const size_t num_blocks = domain.blocks().size();
  const auto& initial_extents = creator->initial_extents();
  const auto functions_of_time = creator->functions_of_time();
  for (size_t b = 0; b < num_blocks; ++b) {
    const ElementId<Dim> element_id{b};
    const auto& block = domain.blocks()[b];
    const Mesh<Dim> mesh{initial_extents[b], Spectral::Basis::Legendre,
                         Spectral::Quadrature::GaussLobatto};
    const auto x_l = logical_coordinates(mesh);
    if constexpr (grid_is == grid::Is::Stationary) {
      const ElementMap<Dim, Frame::Inertial> element_map{element_id, block};
      const auto x_i = element_map(x_l);
      const auto grid_velocity = grid::v<grid_is>(x_i, t);
      const auto div_v = grid::div_v<grid_is, Dim>(x_i, t);
      CHECK_FALSE(grid_velocity.has_value());
      CHECK_FALSE(div_v.has_value());
      const auto inv_jac = element_map.inv_jacobian(x_l);

      check_vars<grid_is>(generator, mesh, inv_jac, x_i, t, grid_velocity);

    } else {
      const ElementMap<Dim, Frame::Grid> element_map_g{element_id, block};
      const auto x_g = element_map_g(x_l);
      const auto& grid_to_inertial_map =
          block.moving_mesh_grid_to_inertial_map();
      const auto moving_mesh_quantities =
          grid_to_inertial_map.coords_frame_velocity_jacobians(
              x_g, t, functions_of_time);
      const auto& x_i = std::get<0>(moving_mesh_quantities);
      const auto grid_velocity = grid::v<grid_is>(x_i, t);
      const auto div_v = grid::div_v<grid_is, Dim>(x_i, t);
      const auto& expected_grid_velocity = std::get<3>(moving_mesh_quantities);
      CHECK_ITERABLE_APPROX(expected_grid_velocity, grid_velocity.value());
      const ElementMap<Dim, Frame::Inertial> element_map_i{element_id, block};
      const auto inv_jac =
          element_map_i.inv_jacobian(x_l, t, functions_of_time);
      const auto expected_div_v =
          divergence(expected_grid_velocity, mesh, inv_jac);
      CHECK_ITERABLE_APPROX(expected_div_v, div_v.value());

      check_vars<grid_is>(generator, mesh, inv_jac, x_i, t, grid_velocity);
    }
  }
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Helpers.Evolution.DG.Variables",
                  "[Unit][Evolution][Actions]") {
  MAKE_GENERATOR(generator);

  test<grid::Is::Stationary, 1>(&generator);
  test<grid::Is::Stationary, 2>(&generator);
  test<grid::Is::Stationary, 3>(&generator);
  test<grid::Is::Comoving, 1>(&generator);
  test<grid::Is::Comoving, 2>(&generator);
  test<grid::Is::Comoving, 3>(&generator);
  test<grid::Is::Expanding, 1>(&generator);
  test<grid::Is::Expanding, 2>(&generator);
  test<grid::Is::Expanding, 3>(&generator);
}
}  // namespace TestHelpers::evolution::dg::Actions
