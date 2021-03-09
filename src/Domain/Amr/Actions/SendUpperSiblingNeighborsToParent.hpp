// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <unordered_map>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Amr/Actions/AdjustDomain.hpp"
#include "Domain/Amr/Actions/UpdateParentFromUpperChild.hpp"
#include "Domain/Amr/Flag.hpp"
#include "Domain/Amr/Tags.hpp"
#include "Domain/Structure//ElementId.hpp"
#include "Domain/Structure//Neighbors.hpp"
#include "Domain/Structure//Side.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Tags.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Printf.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"

namespace amr::Actions {

struct AdjustDomain;

struct SendUpperSiblingNeighborsToParent {
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

      const auto first_join_flag = std::find(
          std::begin(my_amr_flags), std::end(my_amr_flags), amr::Flag::Join);
      ASSERT(first_join_flag != std::end(my_amr_flags),
             "Upper sibling does not want to join?");
      size_t dim_to_join = static_cast<size_t>(
          std::distance(std::begin(my_amr_flags), first_join_flag));
      ASSERT(dim_to_join < volume_dim, "dim_to_join = " << dim_to_join);
      ElementId<volume_dim> parent_id = element_id.id_of_parent(dim_to_join);
      typename Element<volume_dim>::Neighbors_t parent_potential_neighbors =
          element.neighbors();
      parent_potential_neighbors.erase(
          Direction<volume_dim>(dim_to_join, Side::Lower));
      auto& element_array =
          Parallel::get_parallel_component<ParallelComponent>(cache);
      Parallel::printf("Sending upper sibling %s neighbors to parent %s\n",
                       element_id, parent_id);
      Parallel::simple_action<UpdateParentFromUpperChild>(
          element_array[parent_id], parent_potential_neighbors,
          amr_flags_of_neighbors);
      Parallel::simple_action<AdjustDomain>(element_array[parent_id]);
      element_array(array_index).ckDestroy();
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
