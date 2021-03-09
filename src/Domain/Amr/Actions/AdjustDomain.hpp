// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <unordered_map>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Amr/Actions/InitializeChild.hpp"
#include "Domain/Amr/Actions/InitializeParent.hpp"
#include "Domain/Amr/Actions/UpdateParentFromUpperChild.hpp"
#include "Domain/Amr/Flag.hpp"
#include "Domain/Amr/Helpers.hpp"
#include "Domain/Amr/Tags.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/Neighbors.hpp"
#include "Domain/Structure/OrientationMap.hpp"
#include "Domain/Structure/Side.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Printf.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"

namespace amr::Actions {

struct AdjustDomain {
  template <typename DbTagList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagList>&&> apply(
      db::DataBox<DbTagList>& box,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    apply<ParallelComponent>(box, cache, array_index);
    return {std::move(box)};
  }

  template <typename ParallelComponent, typename DbTagList,
            typename Metavariables, typename ArrayIndex>
  static void apply(db::DataBox<DbTagList>& box,
                    const Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& array_index) noexcept {
    constexpr size_t volume_dim = Metavariables::volume_dim;
    if constexpr (tmpl::list_contains_v<
                      DbTagList, amr::Tags::NeighborFlags<volume_dim>> and
                  tmpl::list_contains_v<DbTagList,
                                        amr::Tags::Flags<volume_dim>> and
                  tmpl::list_contains_v<DbTagList,
                                        ::domain::Tags::Element<volume_dim>>) {
      const ElementId<volume_dim> element_id{array_index};

      const auto& my_amr_flags = get<amr::Tags::Flags<volume_dim>>(box);
      const auto& amr_flags_of_neighbors =
          get<amr::Tags::NeighborFlags<volume_dim>>(box);
      const auto& element = get<::domain::Tags::Element<volume_dim>>(box);

      Parallel::printf("Adjusting %s with flags %s and neighbor flags %s\n",
                       element, my_amr_flags, amr_flags_of_neighbors);

      auto& element_array =
          Parallel::get_parallel_component<ParallelComponent>(cache);
      tuples::tagged_tuple_from_typelist<
          typename ParallelComponent::initialization_tags>
          initialization_items{};

      // splitting is done first, in case it enables joining in another
      // dimension
      const auto first_split_flag = std::find(
          std::begin(my_amr_flags), std::end(my_amr_flags), amr::Flag::Split);
      if (first_split_flag != std::end(my_amr_flags)) {
        const size_t dim_to_split = static_cast<size_t>(
            std::distance(std::begin(my_amr_flags), first_split_flag));
        ASSERT(dim_to_split < volume_dim, "dim_to_split = " << dim_to_split);
        std::array<amr::Flag, volume_dim> child_amr_decision = my_amr_flags;
        child_amr_decision[dim_to_split] = amr::Flag::DoNothing;

        typename Element<volume_dim>::Neighbors_t
            lower_child_potential_neighbors = element.neighbors();
        const Direction<volume_dim> upper_direction =
            Direction<volume_dim>(dim_to_split, Side::Upper);
        const ElementId<volume_dim> upper_child_id =
            element_id.id_of_child(dim_to_split, Side::Upper);
        lower_child_potential_neighbors[upper_direction] =
            Neighbors<volume_dim>({{upper_child_id}},
                                  OrientationMap<volume_dim>{});

        std::unordered_map<ElementId<volume_dim>,
                           std::array<amr::Flag, volume_dim>>
            lower_child_potential_neighbor_amr_flags = amr_flags_of_neighbors;
        lower_child_potential_neighbor_amr_flags[upper_child_id] =
            child_amr_decision;

        ElementId<volume_dim> lower_child_id =
            element_id.id_of_child(dim_to_split, Side::Lower);
        Parallel::printf("Creating Element %s\n", lower_child_id);
        element_array(lower_child_id)
            .insert(cache.thisProxy, Metavariables::Phase::AdjustDomain);
        Parallel::simple_action<amr::Actions::InitializeChild<volume_dim>>(
            element_array[lower_child_id], child_amr_decision,
            lower_child_potential_neighbors,
            lower_child_potential_neighbor_amr_flags);
        Parallel::simple_action<amr::Actions::AdjustDomain>(
            element_array[lower_child_id]);

        typename Element<volume_dim>::Neighbors_t
            upper_child_potential_neighbors = element.neighbors();
        const Direction<volume_dim> lower_direction =
            Direction<volume_dim>(dim_to_split, Side::Lower);
        upper_child_potential_neighbors[lower_direction] =
            Neighbors<volume_dim>({{lower_child_id}},
                                  OrientationMap<volume_dim>{});

        std::unordered_map<ElementId<volume_dim>,
                           std::array<amr::Flag, volume_dim>>
            upper_child_potential_neighbor_amr_flags = amr_flags_of_neighbors;
        upper_child_potential_neighbor_amr_flags[lower_child_id] =
            child_amr_decision;

        Parallel::printf("Creating Element %s\n", upper_child_id);
        element_array(upper_child_id)
            .insert(cache.thisProxy, Metavariables::Phase::AdjustDomain);
        Parallel::simple_action<amr::Actions::InitializeChild<volume_dim>>(
            element_array[upper_child_id], child_amr_decision,
            upper_child_potential_neighbors,
            upper_child_potential_neighbor_amr_flags);
        Parallel::simple_action<amr::Actions::AdjustDomain>(
            element_array[upper_child_id]);

        // delete myself
        element_array(array_index).ckDestroy();
      } else {  // not splitting, so check for joining
        const auto first_join_flag = std::find(
            std::begin(my_amr_flags), std::end(my_amr_flags), amr::Flag::Join);
        if (first_join_flag != std::end(my_amr_flags)) {
          const size_t dim_to_join = amr::dimension_to_join(
              element, my_amr_flags, amr_flags_of_neighbors);
          ASSERT(dim_to_join < volume_dim, "dim_to_join = " << dim_to_join);
          ElementId<volume_dim> parent_id =
              element_id.id_of_parent(dim_to_join);
          typename Element<volume_dim>::Neighbors_t parent_potential_neighbors =
              element.neighbors();
          std::array<amr::Flag, volume_dim> parent_amr_decision = my_amr_flags;
          parent_amr_decision[dim_to_join] = amr::Flag::DoNothing;
          // Of the two Elements to be joined, the lower sibling will create the
          // parent chare.
          if (element_id == parent_id.id_of_child(dim_to_join, Side::Lower)) {
            parent_potential_neighbors.erase(
                Direction<volume_dim>(dim_to_join, Side::Upper));
            Parallel::printf("Creating Parent %s\n", parent_id);
            element_array(parent_id).insert(cache.thisProxy,
                                            Metavariables::Phase::AdjustDomain);
          } else {
            parent_potential_neighbors.erase(
                Direction<volume_dim>(dim_to_join, Side::Lower));
            Parallel::printf("Joining Parent %s\n", parent_id);
          }
          Parallel::simple_action<amr::Actions::InitializeParent<volume_dim>>(
              element_array[parent_id], dim_to_join, parent_amr_decision,
              parent_potential_neighbors, amr_flags_of_neighbors);
          //         element_array(array_index).ckDestroy();
          Parallel::printf("\n");
        } else {  // neither splitting nor joining, determine neighbors
          Parallel::printf("Element %s is adjusting neighbors\n", element_id);

          db::mutate<::domain::Tags::Element<volume_dim>,
                     amr::Tags::Flags<volume_dim>,
                     amr::Tags::NeighborFlags<volume_dim>>(
              make_not_null(&box),
              [&element_id](
                  const gsl::not_null<Element<volume_dim>*> element,
                  const gsl::not_null<std::array<amr::Flag, volume_dim>*>
                      amr_flags,
                  const gsl::not_null<
                      std::unordered_map<ElementId<volume_dim>,
                                         std::array<amr::Flag, volume_dim>>*>
                      amr_flags_of_neighbors) noexcept {
                auto new_neighbors = element->neighbors();
                for (auto& direction_neighbors : new_neighbors) {
                  const Direction<volume_dim>& direction =
                      direction_neighbors.first;
                  auto& neighbors = direction_neighbors.second;
                  neighbors.set_ids_to(
                      amr::new_neighbor_ids(element_id, direction, neighbors,
                                            *amr_flags_of_neighbors));
                }
                *element = Element<volume_dim>(element_id, new_neighbors);
                amr_flags_of_neighbors->clear();
                for (size_t d = 0; d < volume_dim; ++d) {
                  (*amr_flags)[d] = amr::Flag::Undefined;
                }
              });
          Parallel::printf("\nFinished %s\n",
                           get<::domain::Tags::Element<volume_dim>>(box));
          ASSERT(amr::neighbors_are_within_one_refinement_level(element),
                 "Element = " << element);
        }
      }
    } else {
      ERROR("Could not find the tag "
            << pretty_type::get_name<amr::Tags::Flags<volume_dim>>() << " or "
            << pretty_type::get_name<amr::Tags::NeighborFlags<volume_dim>>()
            << " or "
            << pretty_type::get_name<::domain::Tags::Element<volume_dim>>()
            << " in the DataBox.");
    }
  }
};
}  // namespace amr::Actions
