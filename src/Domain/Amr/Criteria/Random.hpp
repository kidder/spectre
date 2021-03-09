// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <random>

#include "Domain/Amr/Flag.hpp"

namespace amr {
namespace RefinementCriteria {
amr::Flag random_flag(const size_t current_refinement_level,
                      const size_t maximum_refinement_level) {
  static constexpr double fraction_do_something = 0.85;
  static std::random_device r;
  static const auto seed = r();
  static std::mt19937 generator(seed);
  ;
  static std::uniform_real_distribution<> distribution(0.0, 1.0);

  const double random_number = distribution(generator);
  if (random_number > fraction_do_something) {
    return amr::Flag::DoNothing;
  }
  const double join_fraction =
      current_refinement_level / static_cast<double>(maximum_refinement_level);
  if (random_number < join_fraction * fraction_do_something) {
    return amr::Flag::Join;
  }
  return amr::Flag::Split;
}
}  // namespace RefinementCriteria
}  // namespace amr
