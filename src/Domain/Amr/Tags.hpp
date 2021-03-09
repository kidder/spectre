// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines tags related to adaptive mesh refinement

#pragma once

#include <array>
#include <cstddef>
#include <unordered_map>

#include "DataStructures/DataBox/Tag.hpp"
#include "Domain/Amr/Flag.hpp"
#include "Domain/Structure/SegmentId.hpp"
#include "Options/Options.hpp"

/// \cond
template <size_t VolumeDim>
class ElementId;
/// \encdoncd

namespace amr {
namespace OptionTags {
struct MaximumRefinementLevel {
  using type = size_t;
  static constexpr Options::String help = {"Maximum refinement level"};
  static type upper_bound() noexcept { return SegmentId::max_refinement_level; }
};

struct NumberOfIterations {
  using type = size_t;
  static constexpr Options::String help = {"Number of AMR iterations"};
};
}  // namespace OptionTags

namespace Tags {
/// \ingroup DataBoxTagsGroup
/// \ingroup ComputationalDomainGroup
/// \brief amr::Flag%s for an Element.
template <size_t VolumeDim>
struct Flags : db::SimpleTag {
  using type = std::array<amr::Flag, VolumeDim>;
};

/// \ingroup DataBoxTagsGroup
/// \ingroup ComputationalDomainGroup
/// \brief amr::Flag%s for the neighbors of an Element.
template <size_t VolumeDim>
struct NeighborFlags : db::SimpleTag {
  using type = std::unordered_map<ElementId<VolumeDim>,
                                  std::array<amr::Flag, VolumeDim>>;
};

struct MaximumRefinementLevel : db::SimpleTag {
  using type = size_t;
  static constexpr bool pass_metavariables = false;
  using option_tags = tmpl::list<OptionTags::MaximumRefinementLevel>;
  static type create_from_options(const type& option) { return option; }
};

struct NumberOfIterations : db::SimpleTag {
  using type = size_t;
  static constexpr bool pass_metavariables = false;
  using option_tags = tmpl::list<OptionTags::NumberOfIterations>;
  static type create_from_options(const type& option) { return option; }
};
}  // namespace Tags
}  // namespace amr
