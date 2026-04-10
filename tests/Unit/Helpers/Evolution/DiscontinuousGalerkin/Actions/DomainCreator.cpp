// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/DomainCreator.hpp"

#include <array>
#include <cstddef>
#include <memory>

#include "Domain/Creators/NonconformingSphericalShells.hpp"
#include "Domain/Creators/RotatedBricks.hpp"
#include "Domain/Creators/RotatedIntervals.hpp"
#include "Domain/Creators/RotatedRectangles.hpp"
#include "Domain/Creators/TimeDependence/TimeDependence.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/BoundaryCondition.hpp"

namespace TestHelpers::evolution::dg::Actions {
template <bool BlocksAreConforming, size_t Dim>
std::unique_ptr<DomainCreator<Dim>> domain_creator(
    const domain::creators::time_dependence::TimeDependence<Dim>&
        time_dependence) {
  if constexpr (BlocksAreConforming) {
    if constexpr (Dim == 1) {
      return std::make_unique<domain::creators::RotatedIntervals>(
          std::array{-3.0}, std::array{1.5}, std::array{2.0}, std::array{0_st},
          std::array<std::array<size_t, 2>, 1>{std::array{7_st, 7_st}},
          std::make_unique<BoundaryCondition<1>>(),
          std::make_unique<BoundaryCondition<1>>(),
          time_dependence.get_clone());
    }
    if constexpr (Dim == 2) {
      return std::make_unique<domain::creators::RotatedRectangles>(
          std::array{-3.0, -2.0}, std::array{1.5, 4.5}, std::array{2.0, 7.0},
          std::array{0_st, 0_st},
          std::array<std::array<size_t, 2>, 2>{std::array{7_st, 8_st},
                                               std::array{7_st, 8_st}},
          std::make_unique<BoundaryCondition<2>>(),
          time_dependence.get_clone());
    }
    if constexpr (Dim == 3) {
      return std::make_unique<domain::creators::RotatedBricks>(
          std::array{-3.0, -2.0, -1.0}, std::array{1.5, 4.5, 2.5},
          std::array{2.0, 7.0, 3.0}, std::array{0_st, 0_st, 0_st},
          std::array<std::array<size_t, 2>, 3>{std::array{7_st, 8_st},
                                               std::array{7_st, 8_st},
                                               std::array{7_st, 8_st}},
          std::make_unique<BoundaryCondition<3>>(),
          time_dependence.get_clone());
    }
  } else {
    if constexpr (Dim == 3) {
      return std::make_unique<domain::creators::NonconformingSphericalShells>(
          1.9, 2.4, 2.9, 0, 0, 8, 7, 20,
          std::make_unique<BoundaryCondition<3>>(),
          std::make_unique<BoundaryCondition<3>>());
    } else {
      return nullptr;
    }
  }
}

template std::unique_ptr<DomainCreator<1>> domain_creator<true>(
    const domain::creators::time_dependence::TimeDependence<1>&);
template std::unique_ptr<DomainCreator<2>> domain_creator<true>(
    const domain::creators::time_dependence::TimeDependence<2>&);
template std::unique_ptr<DomainCreator<3>> domain_creator<true>(
    const domain::creators::time_dependence::TimeDependence<3>&);

template std::unique_ptr<DomainCreator<1>> domain_creator<false>(
    const domain::creators::time_dependence::TimeDependence<1>&);
template std::unique_ptr<DomainCreator<2>> domain_creator<false>(
    const domain::creators::time_dependence::TimeDependence<2>&);
template std::unique_ptr<DomainCreator<3>> domain_creator<false>(
    const domain::creators::time_dependence::TimeDependence<3>&);
}  // namespace TestHelpers::evolution::dg::Actions
