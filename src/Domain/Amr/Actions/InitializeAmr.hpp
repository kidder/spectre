// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <tuple>
#include <utility>

#include "Domain/Amr/Flag.hpp"
#include "Domain/Amr/Tags.hpp"
#include "ParallelAlgorithms/Initialization/MutateAssign.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
template <size_t VolumeDim>
class ElementId;
namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel
namespace tuples {
template <class... Tags>
class TaggedTuple;
}  // namespace tuples
/// \endcond

namespace amr::Actions {
/*!
 * \ingroup InitializationGroup
 * \brief Initialize items related to adaptive mesh refinement
 *
 * DataBox:
 * - Adds:
 *   - `amr::Tags::Flags<Dim>`
 *   - `amr::Tags::NeighborFlags<Dim>`
 * - Removes: nothing
 * - Modifies: nothing
 *
 * \note This action relies on the `SetupDataBox` aggregated initialization
 * mechanism, so `Actions::SetupDataBox` must be present in the `Initialization`
 * phase action list prior to this action.
 */
template <size_t Dim>
struct InitializeAmr {
  using simple_tags =
      tmpl::list<amr::Tags::Flags<Dim>, amr::Tags::NeighborFlags<Dim>>;

  template <typename DataBox, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static auto apply(DataBox& box,
                    const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ElementId<Dim>& /*array_index*/,
                    const ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) noexcept {
    std::array<amr::Flag, Dim> amr_flags =
        make_array<Dim>(amr::Flag::Undefined);

    // default initialization of NeighborFlags is okay
    ::Initialization::mutate_assign<tmpl::list<amr::Tags::Flags<Dim>>>(
        make_not_null(&box), std::move(amr_flags));
    return std::make_tuple(std::move(box));
  }
};
}  // namespace amr::Actions
