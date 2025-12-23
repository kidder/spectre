// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Domain/Structure/ElementId.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Utilities/TMPL.hpp"

namespace TestHelpers::evolution::dg::Actions {
template <typename Metavariables, typename InitializationActions,
          typename TestActions>
struct Component {
  using metavariables = Metavariables;
  using initialization_actions = InitializationActions;
  using test_actions = TestActions;

  static constexpr size_t volume_dim = metavariables::volume_dim;
  using array_index = ElementId<volume_dim>;
  using chare_type = ActionTesting::MockArrayChare;

  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<Parallel::Phase::Initialization,
                             typename metavariables::initialization_actions>,
      Parallel::PhaseActions<Parallel::Phase::Testing,
                             typename metavariables::test_actions>>;
};
}  // namespace TestHelpers::evolution::dg::Actions
