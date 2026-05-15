// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>

#include "Domain/BoundaryConditions/BoundaryCondition.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"

namespace TestHelpers::evolution {
template <size_t Dim>
struct BoundaryCondition final
    : public domain::BoundaryConditions::BoundaryCondition {
  /// \cond
  explicit BoundaryCondition(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(BoundaryCondition);  // NOLINT
  /// \endcond

  BoundaryCondition() = default;

  std::unique_ptr<domain::BoundaryConditions::BoundaryCondition> get_clone()
      const override {
    return std::make_unique<BoundaryCondition>(*this);
  }
};

template <size_t Dim>
PUP::able::PUP_ID BoundaryCondition<Dim>::my_PUP_ID = 0;
}  // namespace TestHelpers::evolution
