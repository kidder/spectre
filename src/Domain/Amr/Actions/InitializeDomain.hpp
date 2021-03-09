// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <tuple>
#include <utility>

#include "Domain/Block.hpp"
#include "Domain/CreateInitialElement.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Tags.hpp"
#include "ParallelAlgorithms/Initialization/MutateAssign.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
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
 * \brief Initialize items related to the basic structure of the element
 *
 * GlobalCache:
 * - Uses:
 *   - `domain::Tags::Domain<Dim, Frame::Inertial>`
 * DataBox:
 * - Uses:
 *   - `domain::Tags::InitialRefinementLevels<Dim>`
 * - Adds:
 *   - `Tags::Element<Dim>`
 * - Removes: nothing
 * - Modifies: nothing
 *
 * \note This action relies on the `SetupDataBox` aggregated initialization
 * mechanism, so `Actions::SetupDataBox` must be present in the `Initialization`
 * phase action list prior to this action.
 */
template <size_t Dim>
struct InitializeDomain {
  using const_global_cache_tags = tmpl::list<domain::Tags::Domain<Dim>>;
  using initialization_tags =
      tmpl::list<domain::Tags::InitialRefinementLevels<Dim>>;

  using simple_tags = tmpl::list<domain::Tags::Element<Dim>>;

  template <typename DataBox, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static auto apply(DataBox& box,
                    const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ElementId<Dim>& array_index,
                    const ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) noexcept {
    const auto& initial_refinement =
        db::get<domain::Tags::InitialRefinementLevels<Dim>>(box);
    const auto& domain = db::get<domain::Tags::Domain<Dim>>(box);
    const ElementId<Dim> element_id{array_index};
    const auto& my_block = domain.blocks()[element_id.block_id()];
    Element<Dim> element = domain::Initialization::create_initial_element(
        element_id, my_block, initial_refinement);
    ::Initialization::mutate_assign<simple_tags>(make_not_null(&box),
                                                 std::move(element));
    return std::make_tuple(std::move(box));
  }
};
}  // namespace amr::Actions
