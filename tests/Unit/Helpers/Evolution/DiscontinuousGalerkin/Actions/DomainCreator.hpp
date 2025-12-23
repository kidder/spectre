// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>

/// \cond
template <size_t Dim>
class DomainCreator;
/// \endcond

namespace TestHelpers::evolution::dg::Actions {
template <size_t Dim>
std::unique_ptr<DomainCreator<Dim>> domain_creator();
}  // namespace TestHelpers::evolution::dg::Actions
