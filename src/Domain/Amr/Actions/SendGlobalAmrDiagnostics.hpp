// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <boost/rational.hpp>
#include <cstddef>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Amr/Actions/CheckDomain.hpp"
#include "Domain/Amr/Components/AmrMonitor.hpp"
#include "Domain/Amr/Helpers.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Reduction.hpp"
#include "Utilities/BoostHelpers.hpp"

namespace amr::Actions {

/// \brief Send global AMR diagnostics about Element%s to the AmrMonitor
///
/// Sends the following:
/// - The fraction of a Block volume (in the logical coordinate frame) covered
///   by the Element
///
/// The information is sent to AmrMonitor which runs the action
/// amr::Actions::CheckDomain after all Element%s have contributed to the
/// reduction
struct SendGlobalAmrDiagnostics {
  template <typename DbTagList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagList>&&> apply(
      db::DataBox<DbTagList>& box,
      tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, const ActionList /*meta*/,
      const ParallelComponent* /*meta*/) noexcept {
    const ElementId<Metavariables::volume_dim> element_id{array_index};
    const auto& my_proxy =
        Parallel::get_parallel_component<ParallelComponent>(cache)[array_index];
    const auto& target_proxy = Parallel::get_parallel_component<
        Parallel::Components::AmrMonitor<Metavariables>>(cache);
    Parallel::contribute_to_reduction<amr::Actions::CheckDomain>(
        Parallel::ReductionData<
            Parallel::ReductionDatum<boost::rational<size_t>, funcl::Plus<>>>{
            amr::fraction_of_block_volume(element_id)},
        my_proxy, target_proxy);
    return {std::move(box)};
  }
};
}  // namespace amr::Actions
