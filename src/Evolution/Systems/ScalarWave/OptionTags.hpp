// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "Evolution/Tags.hpp"
#include "Options/String.hpp"

namespace ScalarWave::OptionTags {
/// \ingroup OptionTagsGroup
/// \ingroup TimeGroup
/// \brief One-index constraint parameter gamma_2
struct Gamma2 {
  using type = double;
  static constexpr Options::String help = {
      "One-index constraint parameter gamma_2."};
  static type lower_bound() { return 0.0; }
  using group = evolution::OptionTags::Group;
};

struct Mass {
  using type = double;
  static constexpr Options::String help = {"Mass of the scalar field"};
  static type lower_bound() { return 0.0; }
  using group = evolution::OptionTags::Group;
};
}  // namespace ScalarWave::OptionTags
