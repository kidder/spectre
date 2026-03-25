// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Domain/Structure/ElementId.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Utilities/TMPL.hpp"

namespace TestHelpers::evolution::dg::Actions {
template <typename Metavariables>
struct Component {
  static constexpr size_t volume_dim = Metavariables::volume_dim;

  using metavariables = Metavariables;

  using array_index = ElementId<volume_dim>;
  using chare_type = ActionTesting::MockArrayChare;

  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<Parallel::Phase::Initialization,
                             typename Metavariables::initialization_actions>,
      Parallel::PhaseActions<Parallel::Phase::Testing,
                             typename Metavariables::test_actions>>;
};
}  // namespace TestHelpers::evolution::dg::Actions
