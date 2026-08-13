// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <unordered_set>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Amr/Tags.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"

namespace amr::Actions {
struct PingNeighbors {
  template <typename ParallelComponent, typename DbTagList,
            typename Metavariables>
  static void apply(
      db::DataBox<DbTagList>& box,
      Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ElementId<Metavariables::volume_dim>& /*element_id*/,
      const ElementId<Metavariables::volume_dim>& neighbor_id,
      const Direction<Metavariables::volume_dim>& direction_to_neighbor) {
    constexpr size_t volume_dim = Metavariables::volume_dim;
    db::mutate<amr::Tags::ReceivedNeighborPings<volume_dim>>(
        [&neighbor_id, direction_to_neighbor](
            const gsl::not_null<DirectionMap<
                volume_dim, std::unordered_set<ElementId<volume_dim>>>*>
                neighbor_pings) {
          if (neighbor_pings->contains(direction_to_neighbor)) {
            auto& neighbors = neighbor_pings->at(direction_to_neighbor);
            if (neighbors.contains(neighbor_id)) {
              ERROR("Attempted to insert " << neighbor_id
                                           << " a second time into "
                                           << *neighbor_pings);
            } else {
              neighbors.emplace(neighbor_id);
            }
          } else {
            neighbor_pings->emplace(direction_to_neighbor,
                                    std::unordered_set{neighbor_id});
          }
        },
        make_not_null(&box));
  }
};
}  // namespace amr::Actions
