// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <unordered_set>

#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Framework/ActionTesting.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Amr/Actions/PingNeighbors.hpp"
#include "ParallelAlgorithms/Amr/Tags.hpp"
#include "Utilities/TMPL.hpp"

namespace {
template <typename Metavariables>
struct Component {
  using metavariables = Metavariables;
  static constexpr size_t volume_dim = Metavariables::volume_dim;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = ElementId<volume_dim>;
  using const_global_cache_tags = tmpl::list<>;
  using simple_tags = tmpl::list<amr::Tags::ReceivedNeighborPings<volume_dim>>;
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      Parallel::Phase::Initialization,
      tmpl::list<ActionTesting::InitializeDataBox<simple_tags>>>>;
};

template <size_t VolumeDim>
struct Metavariables {
  static constexpr size_t volume_dim = VolumeDim;
  using component_list = tmpl::list<Component<Metavariables>>;
};

template <size_t VolumeDim>
void test() {
  using metavariables = Metavariables<VolumeDim>;
  using my_component = Component<metavariables>;
  DirectionMap<VolumeDim, std::unordered_set<ElementId<VolumeDim>>> neighbors{};
  ActionTesting::MockRuntimeSystem<metavariables> runner{{}};
  const ElementId<VolumeDim> self_id{0};
  ActionTesting::emplace_component_and_initialize<my_component>(
      &runner, self_id, {neighbors});
  size_t block_id = 0;
  for (const auto& direction : Direction<VolumeDim>::all_directions()) {
    const ElementId<VolumeDim> first_neighbor_id{++block_id};
    const ElementId<VolumeDim> second_neighbor_id{++block_id};

    neighbors.emplace(
        direction, std::unordered_set{first_neighbor_id, second_neighbor_id});
  }
  CHECK(ActionTesting::get_databox_tag<
            my_component, amr::Tags::ReceivedNeighborPings<VolumeDim>>(runner,
                                                                       self_id)
            .empty());
  for (const auto& [direction, neighbors_in_direction] : neighbors) {
    for (const auto& neighbor : neighbors_in_direction) {
      ActionTesting::simple_action<my_component, amr::Actions::PingNeighbors>(
          make_not_null(&runner), self_id, neighbor, direction);
    }
  }
  const auto& pinging_neighbors = ActionTesting::get_databox_tag<
      my_component, amr::Tags::ReceivedNeighborPings<VolumeDim>>(runner,
                                                                 self_id);
  CHECK(pinging_neighbors.size() == 2 * VolumeDim);
  block_id = 0;
  for (const auto& direction : Direction<VolumeDim>::all_directions()) {
    REQUIRE(pinging_neighbors.contains(direction));
    const auto& neighbors_in_direction = neighbors.at(direction);
    CHECK(neighbors_in_direction.size() == 2);
    CHECK(neighbors_in_direction.contains(ElementId<VolumeDim>{++block_id}));
    CHECK(neighbors_in_direction.contains(ElementId<VolumeDim>{++block_id}));
  }
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Amr.Actions.PingNeighbors",
                  "[Unit][ParallelAlgorithms]") {
  test<1>();
  test<2>();
  test<3>();
}
