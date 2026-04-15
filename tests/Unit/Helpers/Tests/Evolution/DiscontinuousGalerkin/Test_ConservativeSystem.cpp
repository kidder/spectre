// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>

#include <array>
#include <cstddef>
#include <memory>
#include <random>
#include <tuple>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/Creators/TimeDependence/CubicScale.hpp"
#include "Domain/Creators/TimeDependence/None.hpp"
#include "Domain/Creators/TimeDependence/TimeDependence.hpp"
#include "Domain/Creators/TimeDependence/UniformTranslation.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Evolution/DiscontinuousGalerkin/TimeDerivativeDecisions.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/DataStructures/MakeRandomVectorInMagnitudeRange.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/ConservativeSystem.hpp"
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
void test_dg_package_data(
    const gsl::not_null<std::mt19937*> generator,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& x_i, const double t,
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
        grid_velocity) {
  const size_t n_pts = get<0>(x_i).size();

  const auto n =
      make_random_vector_in_magnitude_range_flat<DataVector, Dim, UpLo::Lo,
                                                 Frame::Inertial>(
          generator, x_i, 1.0, 1.0);

  std::optional<Scalar<DataVector>> n_dot_grid_velocity = std::nullopt;
  if (grid_velocity.has_value()) {
    n_dot_grid_velocity = dot_product(n, grid_velocity.value());
  }

  Scalar<DataVector> var1{n_pts};
  Var1::value(make_not_null(&var1), x_i, t);

  tnsr::I<DataVector, Dim, Frame::Inertial> var2{n_pts};
  Var2<Dim>::value(make_not_null(&var2), x_i, t);

  Scalar<DataVector> temp_var{n_pts};
  TempVar::value(make_not_null(&temp_var), x_i, t);

  tnsr::I<DataVector, Dim, Frame::Inertial> flux1{n_pts};
  Var1::flux(make_not_null(&flux1), var2);

  tnsr::IJ<DataVector, Dim, Frame::Inertial> flux2{n_pts};
  Var2<Dim>::flux(&flux2, var1, var2, temp_var);

  Scalar<DataVector> n_dot_flux1{n_pts};
  tnsr::I<DataVector, Dim, Frame::Inertial> n_dot_flux2{n_pts};
  Scalar<DataVector> out_var1{n_pts};
  tnsr::I<DataVector, Dim, Frame::Inertial> out_var2{n_pts};
  Scalar<DataVector> out_temp_var{n_pts};

  ConservativeBoundaryCorrection<Dim>::dg_package_data(
      &n_dot_flux1, &n_dot_flux2, &out_var1, &out_var2, &out_temp_var, var1,
      var2, flux1, flux2, temp_var, n, grid_velocity, n_dot_grid_velocity, t);
  Scalar<DataVector> expected_n_dot_flux1{n_pts};
  Var1::normal_dot_flux<grid_is>(&expected_n_dot_flux1, n, x_i, t, 1.0);

  tnsr::I<DataVector, Dim, Frame::Inertial> expected_n_dot_flux2{n_pts};
  Var2<Dim>::template normal_dot_flux<grid_is>(&expected_n_dot_flux2, n, x_i, t,
                                               1.0);

  CHECK_ITERABLE_APPROX(expected_n_dot_flux1, n_dot_flux1);
  CHECK_ITERABLE_APPROX(expected_n_dot_flux2, n_dot_flux2);
  CHECK_ITERABLE_APPROX(var1, out_var1);
  CHECK_ITERABLE_APPROX(var2, out_var2);
  CHECK_ITERABLE_APPROX(temp_var, out_temp_var);
}

template <size_t Dim>
void test_time_derivative_terms(
    const tnsr::I<DataVector, Dim, Frame::Inertial>& x_i, const double t) {
  const size_t n_pts = get<0>(x_i).size();

  using dt_var1 = ::Tags::dt<Var1>;
  using dt_var2 = ::Tags::dt<Var2<Dim>>;
  using dt_tags = tmpl::list<dt_var1, dt_var2>;
  Variables<dt_tags> dt_vars{n_pts};

  using flux_var1 = ::Tags::Flux<Var1, tmpl::size_t<Dim>, Frame::Inertial>;
  using flux_var2 = ::Tags::Flux<Var2<Dim>, tmpl::size_t<Dim>, Frame::Inertial>;
  using flux_tags = tmpl::list<flux_var1, flux_var2>;
  Variables<flux_tags> flux_vars(n_pts);

  Variables<tmpl::list<TempVar>> temp_vars{n_pts};

  Variables<tmpl::list<Var1, Var2<Dim>>> vars{n_pts};

  auto& var1 = get<Var1>(vars);
  Var1::value(make_not_null(&var1), x_i, t);

  auto& var2 = get<Var2<Dim>>(vars);
  Var2<Dim>::value(make_not_null(&var2), x_i, t);

  Variables<dt_tags> source_vars{n_pts};
  auto& source1 = get<dt_var1>(source_vars);
  Source1::value(make_not_null(&source1), x_i, t);

  auto& source2 = get<dt_var2>(source_vars);
  Source2<Dim>::value(make_not_null(&source2), x_i, t);

  ConservativeTimeDerivativeTermsWithVariables<Dim>::apply(
      &dt_vars, &flux_vars, &temp_vars, var1, var2, source1, source2);

  const auto& expected_dt_vars = source_vars;

  Variables<tmpl::list<TempVar>> expected_temp_vars{n_pts};
  auto& temp_var = get<TempVar>(expected_temp_vars);
  TempVar::value(make_not_null(&temp_var), x_i, t);

  Variables<flux_tags> expected_flux_vars{n_pts};
  auto& flux1 = get<flux_var1>(expected_flux_vars);
  auto& flux2 = get<flux_var2>(expected_flux_vars);
  Var1::flux(make_not_null(&flux1), var2);
  Var2<Dim>::flux(&flux2, var1, var2, temp_var);

  CHECK_VARIABLES_APPROX(expected_dt_vars, dt_vars);
  CHECK_VARIABLES_APPROX(expected_flux_vars, flux_vars);
  CHECK_VARIABLES_APPROX(expected_temp_vars, temp_vars);
}

template <grid::Is grid_is, size_t Dim>
void test(const gsl::not_null<std::mt19937*> generator) {
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
    const ElementMap<Dim, Frame::Inertial> element_map{element_id, block};
    const auto x_i = element_map(x_l, t, functions_of_time);
    const auto grid_velocity = grid::v<grid_is>(x_i, t);
    test_time_derivative_terms(x_i, t);
    test_dg_package_data<grid_is>(generator, x_i, t, grid_velocity);
  }
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Helpers.Evolution.DG.ConservativeSystem",
                  "[Unit][Evolution][Actions]") {
  MAKE_GENERATOR(generator);
  test<grid::Is::Stationary, 1>(&generator);
  test<grid::Is::Comoving, 1>(&generator);
  test<grid::Is::Expanding, 1>(&generator);
  test<grid::Is::Stationary, 2>(&generator);
  test<grid::Is::Comoving, 2>(&generator);
  test<grid::Is::Expanding, 2>(&generator);
  test<grid::Is::Stationary, 3>(&generator);
  test<grid::Is::Comoving, 3>(&generator);
  test<grid::Is::Expanding, 3>(&generator);
}
}  // namespace TestHelpers::evolution::dg::Actions
