// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <tuple>
#include <unordered_map>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Amr/Flag.hpp"
#include "Domain/Amr/Tags.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Tags.hpp"
#include "ParallelAlgorithms/Initialization/MutateAssign.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel
/// \endcond

namespace amr::Actions {
template <size_t Dim>
struct InitializeChild {
  using simple_tags =
      db::AddSimpleTags<::domain::Tags::Element<Dim>, amr::Tags::Flags<Dim>,
                        amr::Tags::NeighborFlags<Dim>>;
  template <typename ParallelComponent, typename DbTagList,
            typename Metavariables, typename ArrayIndex>
  static void apply(
      db::DataBox<DbTagList>& box,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& array_index, std::array<amr::Flag, Dim> amr_flags,
      typename Element<Dim>::Neighbors_t potential_neighbors,
      std::unordered_map<ElementId<Dim>, std::array<amr::Flag, Dim>>
          amr_flags_of_potential_neighbors) noexcept {
    if constexpr (tmpl::list_contains_v<DbTagList,
                                        amr::Tags::NeighborFlags<Dim>> and
                  tmpl::list_contains_v<DbTagList, amr::Tags::Flags<Dim>> and
                  tmpl::list_contains_v<DbTagList,
                                        ::domain::Tags::Element<Dim>>) {
      const ElementId<Dim> element_id{array_index};
      Element<Dim> element(element_id, std::move(potential_neighbors));
      ::Initialization::mutate_assign<simple_tags>(
          make_not_null(&box), std::move(element), std::move(amr_flags),
          std::move(amr_flags_of_potential_neighbors));
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
