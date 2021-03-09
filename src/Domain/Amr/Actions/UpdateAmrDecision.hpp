// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <unordered_map>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Amr/Flag.hpp"
#include "Domain/Amr/Tags.hpp"
#include "Domain/Amr/UpdateAmrDecision.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Tags.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Printf.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace amr::Actions {
struct UpdateAmrDecision {
  template <typename ParallelComponent, typename DbTagList,
            typename Metavariables, typename ArrayIndex>
  static void apply(db::DataBox<DbTagList>& box,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/,
                    const ElementId<Metavariables::volume_dim>& neighbor_id,
                    const std::array<amr::Flag, Metavariables::volume_dim>&
                        neighbor_amr_flags) noexcept {
    if constexpr (tmpl::list_contains_v<
                      DbTagList,
                      amr::Tags::Flags<Metavariables::volume_dim>> and
                  tmpl::list_contains_v<
                      DbTagList,
                      ::domain::Tags::Element<Metavariables::volume_dim>>) {
      constexpr size_t volume_dim = Metavariables::volume_dim;
      std::array<amr::Flag, Metavariables::volume_dim> my_current_amr_flags =
          get<amr::Tags::Flags<volume_dim>>(box);
      const auto& element = get<::domain::Tags::Element<volume_dim>>(box);
      const bool my_amr_decision_changed =
          amr::update_amr_decision(make_not_null(&my_current_amr_flags),
                                   element, neighbor_id, neighbor_amr_flags);
      db::mutate<amr::Tags::Flags<volume_dim>,
                 amr::Tags::NeighborFlags<volume_dim>>(
          make_not_null(&box),
          [&my_current_amr_flags, &neighbor_id, &neighbor_amr_flags ](
              const gsl::not_null<std::array<amr::Flag, volume_dim>*> amr_flags,
              const gsl::not_null<std::unordered_map<
                  ElementId<volume_dim>, std::array<amr::Flag, volume_dim>>*>
                  amr_flags_of_neighbors) noexcept {
            auto& my_neighbor_amr_flags =
                (*amr_flags_of_neighbors)[neighbor_id];
            for (size_t d = 0; d < volume_dim; ++d) {
              (*amr_flags)[d] = my_current_amr_flags[d];
              my_neighbor_amr_flags[d] = neighbor_amr_flags[d];
            }
          });

      if (my_amr_decision_changed) {
        auto& amr_element_array =
            Parallel::get_parallel_component<ParallelComponent>(cache);
        for (const auto direction_neighbors : element.neighbors()) {
          for (const auto id : direction_neighbors.second.ids()) {
            Parallel::simple_action<UpdateAmrDecision>(
                amr_element_array[id], element.id(), my_current_amr_flags);
          }
        }

        Parallel::printf("Updated element %s refinement flags to : %s\n",
                         element.id(), my_current_amr_flags);
      }
    } else {
      ERROR("Could not find the tag "
            << pretty_type::get_name<
                   amr::Tags::Flags<Metavariables::volume_dim>>()
            << " or "
            << pretty_type::get_name<
                   ::domain::Tags::Element<Metavariables::volume_dim>>()
            << " in the DataBox.");
    }
  }
};
}  // namespace amr::Actions
