// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <optional>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Block.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/Creators/TimeDependence/TimeDependence.hpp"
#include "Domain/Domain.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/FunctionsOfTime/FunctionOfTime.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Helpers/Evolution/DomainCreator.hpp"
#include "Helpers/Evolution/TimeDependence.hpp"
#include "NumericalAlgorithms/LinearOperators/Divergence.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/LogicalCoordinates.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"

namespace TestHelpers::evolution::grid {
namespace {
template <Is grid_is, size_t Dim>
void test() {
  CAPTURE(grid_is);
  CAPTURE(Dim);
  const double t = 1.5;
  const auto time_dependence = create_time_dependence<grid_is, Dim>();
  const auto creator =
      TestHelpers::evolution::domain_creator<true>(*time_dependence);
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
    if constexpr (grid_is == Is::Stationary) {
      const ElementMap<Dim, Frame::Inertial> element_map{element_id, block};
      const auto x_i = element_map(x_l);
      const auto grid_velocity = v<grid_is>(x_i, t);
      const auto div_grid_velocity = div_v<grid_is, Dim>(x_i, t);
      CHECK_FALSE(grid_velocity.has_value());
      CHECK_FALSE(div_grid_velocity.has_value());
      const auto inv_jac = element_map.inv_jacobian(x_l);
    } else {
      const ElementMap<Dim, Frame::Grid> element_map_g{element_id, block};
      const auto x_g = element_map_g(x_l);
      const auto& grid_to_inertial_map =
          block.moving_mesh_grid_to_inertial_map();
      const auto moving_mesh_quantities =
          grid_to_inertial_map.coords_frame_velocity_jacobians(
              x_g, t, functions_of_time);
      const auto& x_i = std::get<0>(moving_mesh_quantities);
      const auto grid_velocity = v<grid_is>(x_i, t);
      const auto div_grid_velocity = div_v<grid_is, Dim>(x_i, t);
      const auto& expected_grid_velocity = std::get<3>(moving_mesh_quantities);
      CHECK_ITERABLE_APPROX(expected_grid_velocity, grid_velocity.value());
      const ElementMap<Dim, Frame::Inertial> element_map_i{element_id, block};
      const auto inv_jac =
          element_map_i.inv_jacobian(x_l, t, functions_of_time);
      const auto expected_div_v =
          divergence(expected_grid_velocity, mesh, inv_jac);
      CHECK_ITERABLE_APPROX(expected_div_v, div_grid_velocity.value());
    }
  }
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Helpers.Evolution.TimeDependence",
                  "[Unit][Evolution][Actions]") {
  test<Is::Stationary, 1>();
  test<Is::Stationary, 2>();
  test<Is::Stationary, 3>();
  test<Is::Comoving, 1>();
  test<Is::Comoving, 2>();
  test<Is::Comoving, 3>();
  test<Is::Expanding, 1>();
  test<Is::Expanding, 2>();
  test<Is::Expanding, 3>();
}
}  // namespace TestHelpers::evolution::grid
