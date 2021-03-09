// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <unordered_map>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Amr/Actions/AdjustDomain.hpp"
#include "Domain/Amr/Flag.hpp"
#include "Domain/Amr/Tags.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/Neighbors.hpp"
#include "Domain/Structure/Side.hpp"
#include "Domain/Tags.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Printf.hpp"
#include "Utilities/Gsl.hpp"

namespace amr::Actions {

struct AdjustDomain;

struct UpdateParentFromUpperChild {
  template <typename ParallelComponent, typename DbTagList,
            typename Metavariables, typename ArrayIndex>
  static void apply(
      db::DataBox<DbTagList>& box,
      const Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index,
      typename Element<Metavariables::volume_dim>::Neighbors_t
          potential_neighbors_from_upper_child,
      std::unordered_map<ElementId<Metavariables::volume_dim>,
                         std::array<amr::Flag, Metavariables::volume_dim>>
          amr_flags_of_potential_neighbors_from_upper_child) noexcept {
    constexpr size_t volume_dim = Metavariables::volume_dim;
    if constexpr (tmpl::list_contains_v<
                      DbTagList, amr::Tags::NeighborFlags<volume_dim>> and
                  tmpl::list_contains_v<DbTagList,
                                        ::domain::Tags::Element<volume_dim>>) {
      Parallel::printf("Amr flags = %s\n",
                       get<amr::Tags::Flags<volume_dim>>(box));
      const ElementId<volume_dim> element_id{array_index};
      Parallel::printf("Updating parent %s\n", element_id);
      db::mutate<::domain::Tags::Element<volume_dim>,
                 amr::Tags::NeighborFlags<volume_dim>>(
          make_not_null(&box),
          [
            &element_id, &potential_neighbors_from_upper_child, &
            amr_flags_of_potential_neighbors_from_upper_child
          ](const gsl::not_null<Element<volume_dim>*> element,
            const gsl::not_null<std::unordered_map<
                ElementId<volume_dim>, std::array<amr::Flag, volume_dim>>*>
                amr_flags_of_potential_neighbors) noexcept {
            auto potential_neighbors = element->neighbors();
            for (const auto& direction_neighbors :
                 potential_neighbors_from_upper_child) {
              const Direction<volume_dim>& direction =
                  direction_neighbors.first;
              if (0 == potential_neighbors.count(direction)) {
                potential_neighbors.insert(direction_neighbors);
              } else {
                potential_neighbors.at(direction).add_ids(
                    direction_neighbors.second.ids());
              }
            }
            *element = Element<volume_dim>(element_id, potential_neighbors);
            for (const auto& neighbor_flags :
                 amr_flags_of_potential_neighbors_from_upper_child) {
              if (0 == amr_flags_of_potential_neighbors->count(
                           neighbor_flags.first)) {
                amr_flags_of_potential_neighbors->insert(neighbor_flags);
              }
            }
          });
      auto& element_array =
          Parallel::get_parallel_component<ParallelComponent>(cache);
      Parallel::simple_action<AdjustDomain>(element_array[element_id]);
    } else {
      ERROR("Could not find the tag "
            << pretty_type::get_name<amr::Tags::NeighborFlags<volume_dim>>()
            << " or "
            << pretty_type::get_name<::domain::Tags::Element<volume_dim>>()
            << " in the DataBox.");
    }
  }
};
}  // namespace amr::Actions
