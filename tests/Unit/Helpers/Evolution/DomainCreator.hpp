// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>

/// \cond
template <size_t Dim>
class DomainCreator;
namespace domain::creators::time_dependence {
template <size_t Dim>
struct TimeDependence;
}
/// \endcond

namespace TestHelpers::evolution {
template <bool BlocksAreConforming, size_t Dim>
std::unique_ptr<DomainCreator<Dim>> domain_creator(
    const domain::creators::time_dependence::TimeDependence<Dim>&
        time_dependence);
}  // namespace TestHelpers::evolution
