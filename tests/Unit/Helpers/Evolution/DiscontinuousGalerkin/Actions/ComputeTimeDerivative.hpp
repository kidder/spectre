// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

namespace TestHelpers::evolution::dg::Actions {
/// \cond
enum class SystemType;
/// \endcond

template <SystemType system_type, bool UsePrims, size_t Dim>
void test();
}  // namespace TestHelpers::evolution::dg::Actions
