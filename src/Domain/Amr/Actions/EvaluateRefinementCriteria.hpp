// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Amr/Actions/UpdateAmrDecision.hpp"
#include "Domain/Amr/Criteria/Random.hpp"
#include "Domain/Amr/Flag.hpp"
#include "Domain/Amr/Tags.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/Neighbors.hpp"
#include "Domain/Tags.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Printf.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace tuples {
template <class... Tags>
class TaggedTuple;
}  // namespace tuples
/// \endcond

namespace amr::Actions {
/// \brief Evaluates the refinement criteria in order to set the amr::Flag%s of
/// an Element and sends this information to the neighbors of the Element.
///
/// DataBox changes:
/// - Modifies:
///   * `amr::Tags::Flags<Dim>
///
/// Invokes:
/// - UpdateAmrDecision on all neighboring Element%s
struct EvaluateRefinementCriteria {
  using const_global_cache_tags = tmpl::list<amr::Tags::MaximumRefinementLevel>;

  template <typename DbTagList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagList>&&> apply(
      db::DataBox<DbTagList>& box,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    if constexpr (tmpl::list_contains_v<
                      DbTagList, Parallel::Tags::FromGlobalCache<
                                     amr::Tags::MaximumRefinementLevel>> and
                  tmpl::list_contains_v<
                      DbTagList,
                      amr::Tags::Flags<Metavariables::volume_dim>> and
                  tmpl::list_contains_v<
                      DbTagList,
                      ::domain::Tags::Element<Metavariables::volume_dim>>) {
      constexpr size_t volume_dim = Metavariables::volume_dim;
      const ElementId<volume_dim> element_id{array_index};
      const auto& maximum_refinement_level =
          db::get<amr::Tags::MaximumRefinementLevel>(box);

      db::mutate<amr::Tags::Flags<Metavariables::volume_dim>>(
          make_not_null(&box),
          [&element_id, &maximum_refinement_level ](
              const gsl::not_null<std::array<amr::Flag, volume_dim>*>
                  amr_flags) noexcept {
            for (size_t d = 0; d < volume_dim; ++d) {
              (*amr_flags)[d] = amr::RefinementCriteria::random_flag(
                  element_id.segment_ids()[d].refinement_level(),
                  maximum_refinement_level);
            }
          });

      auto& amr_element_array =
          Parallel::get_parallel_component<ParallelComponent>(cache);

      const auto& my_element = get<::domain::Tags::Element<volume_dim>>(box);
      const std::array<amr::Flag, Metavariables::volume_dim>& my_flags =
          get<amr::Tags::Flags<volume_dim>>(box);
      for (const auto [direction, neighbors] : my_element.neighbors()) {
        for (const auto neighbor_id : neighbors.ids()) {
          //          Parallel::printf("Sending to %s\n", neighbor_id);
          Parallel::simple_action<UpdateAmrDecision>(
              amr_element_array[neighbor_id], element_id, my_flags);
        }
      }

      Parallel::printf("Element %s initial refinement flag: %s\n", element_id,
                       get<amr::Tags::Flags<volume_dim>>(box));
    } else {
      ERROR("Could not find the tag "
            << pretty_type::get_name<
                   amr::Tags::Flags<Metavariables::volume_dim>>()
            << " or "
            << pretty_type::get_name<
                   ::domain::Tags::Element<Metavariables::volume_dim>>()
            << " in the DataBox.");
    }

    return {std::move(box)};
  }
};
}  // namespace amr::Actions
