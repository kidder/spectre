// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/DomainCreator.hpp"

#include <array>
#include <cstddef>
#include <memory>

#include "Domain/Creators/RotatedBricks.hpp"
#include "Domain/Creators/RotatedIntervals.hpp"
#include "Domain/Creators/RotatedRectangles.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/BoundaryCondition.hpp"

namespace TestHelpers::evolution::dg::Actions {
template <size_t Dim>
std::unique_ptr<DomainCreator<Dim>> domain_creator() {
  if constexpr (Dim == 1) {
    return std::make_unique<domain::creators::RotatedIntervals>(
        std::array{-3.0}, std::array{1.5}, std::array{2.0}, std::array{0_st},
        std::array<std::array<size_t, 2>, 1>{std::array{12_st, 12_st}},
        std::make_unique<BoundaryCondition<1>>(),
        std::make_unique<BoundaryCondition<1>>(), nullptr);
  }
  if constexpr (Dim == 2) {
    return std::make_unique<domain::creators::RotatedRectangles>(
        std::array{-3.0, -2.0}, std::array{1.5, 4.5}, std::array{2.0, 7.0},
        std::array{0_st, 0_st},
        std::array<std::array<size_t, 2>, 2>{std::array{6_st, 6_st},
                                             std::array{6_st, 6_st}},
        std::make_unique<BoundaryCondition<2>>());
  }
  if constexpr (Dim == 3) {
    return std::make_unique<domain::creators::RotatedBricks>(
        std::array{-3.0, -2.0, -1.0}, std::array{1.5, 4.5, 2.5},
        std::array{2.0, 7.0, 3.0}, std::array{0_st, 0_st, 0_st},
        std::array<std::array<size_t, 2>, 3>{std::array{6_st, 6_st},
                                             std::array{6_st, 6_st},
                                             std::array{6_st, 6_st}},
        std::make_unique<BoundaryCondition<Dim>>());
  }
}

template std::unique_ptr<DomainCreator<1>> domain_creator();
template std::unique_ptr<DomainCreator<2>> domain_creator();
template std::unique_ptr<DomainCreator<3>> domain_creator();
}  // namespace TestHelpers::evolution::dg::Actions
