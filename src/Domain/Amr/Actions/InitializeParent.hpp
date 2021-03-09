// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <tuple>
#include <unordered_map>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Amr/Actions/SendUpperSiblingNeighborsToParent.hpp"
#include "Domain/Amr/Flag.hpp"
#include "Domain/Amr/Tags.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Printf.hpp"
#include "ParallelAlgorithms/Initialization/MutateAssign.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/TMPL.hpp"

namespace amr::Actions {
template <size_t Dim>
struct InitializeParent {
  using simple_tags =
      db::AddSimpleTags<::domain::Tags::Element<Dim>, amr::Tags::Flags<Dim>,
                        amr::Tags::NeighborFlags<Dim>>;

  template <typename ParallelComponent, typename DbTagList,
            typename Metavariables, typename ArrayIndex>
  static void apply(
      db::DataBox<DbTagList>& box,
      const Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, size_t dim_to_join,
      std::array<amr::Flag, Dim> amr_flags,
      typename Element<Dim>::Neighbors_t potential_neighbors_from_child,
      std::unordered_map<ElementId<Dim>, std::array<amr::Flag, Dim>>
          amr_flags_of_potential_neighbors_from_child) noexcept {
    if constexpr (tmpl::list_contains_v<DbTagList,
                                        amr::Tags::NeighborFlags<Dim>> and
                  tmpl::list_contains_v<DbTagList, amr::Tags::Flags<Dim>> and
                  tmpl::list_contains_v<DbTagList,
                                        ::domain::Tags::Element<Dim>>) {
      // This action will be called twice, once by each child (joining is done
      // one dimension at a time), but which call is executed first is
      // arbitrary.  Determine whether this is the first call by examining the
      // amr flags in the DataBox which are default initialized to Undefined
      // when the chare is created.
      const bool is_first_call = (make_array<Dim>(amr::Flag::Undefined) ==
                                  get<amr::Tags::Flags<Dim>>(box));
      const ElementId<Dim> element_id{array_index};

      if (is_first_call) {
        Parallel::printf("Initializing parent %s\n", element_id);
        Element<Dim> element(element_id,
                             std::move(potential_neighbors_from_child));

        ::Initialization::mutate_assign<simple_tags>(
            make_not_null(&box), std::move(element), std::move(amr_flags),
            std::move(amr_flags_of_potential_neighbors_from_child));
      } else {
        Parallel::printf("Updating parent %s\n", element_id);
        db::mutate<::domain::Tags::Element<Dim>, amr::Tags::NeighborFlags<Dim>>(
            make_not_null(&box),
            [
              &element_id, &potential_neighbors_from_child, &
              amr_flags_of_potential_neighbors_from_child
            ](const gsl::not_null<Element<Dim>*> element,
              const gsl::not_null<std::unordered_map<
                  ElementId<Dim>, std::array<amr::Flag, Dim>>*>
                  amr_flags_of_potential_neighbors) noexcept {
              auto potential_neighbors = element->neighbors();
              for (const auto& direction_neighbors :
                   potential_neighbors_from_child) {
                const Direction<Dim>& direction = direction_neighbors.first;
                if (0 == potential_neighbors.count(direction)) {
                  potential_neighbors.insert(direction_neighbors);
                } else {
                  potential_neighbors.at(direction).add_ids(
                      direction_neighbors.second.ids());
                }
              }
              *element = Element<Dim>(element_id, potential_neighbors);
              for (const auto& neighbor_flags :
                   amr_flags_of_potential_neighbors_from_child) {
                if (0 == amr_flags_of_potential_neighbors->count(
                             neighbor_flags.first)) {
                  amr_flags_of_potential_neighbors->insert(neighbor_flags);
                }
              }
            });
        auto& element_array =
            Parallel::get_parallel_component<ParallelComponent>(cache);
        Parallel::simple_action<AdjustDomain>(element_array[element_id]);
        element_array(element_id.id_of_child(dim_to_join, Side::Lower))
            .ckDestroy();
        element_array(element_id.id_of_child(dim_to_join, Side::Upper))
            .ckDestroy();
      }
    } else {
      ERROR("Could not find the tag "
            << pretty_type::get_name<amr::Tags::Flags<Dim>>() << " or "
            << pretty_type::get_name<amr::Tags::NeighborFlags<Dim>>() << " or "
            << pretty_type::get_name<::domain::Tags::Element<Dim>>()
            << " in the DataBox.");
    }
  }
};
}  // namespace amr::Actions
