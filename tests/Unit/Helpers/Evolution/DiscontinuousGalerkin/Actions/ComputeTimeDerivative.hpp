// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

namespace TestHelpers::evolution::dg::Actions {
/// \cond
enum class SystemType;
/// \endcond

template <size_t Dim, SystemType system_type, bool UsePrims>
void test();
}  // namespace TestHelpers::evolution::dg::Actions
