// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <bit>
#include <concepts>
#include <limits>

#include "Utilities/Requires.hpp"

// NOLINTNEXTLINE(cert-dcl58-cpp)
namespace std {
template <
    std::unsigned_integral T,
    Requires<not std::same_as<T, bool> and not std::same_as<T, char> and
             not std::same_as<T, char8_t> and not std::same_as<T, char16_t> and
             not std::same_as<T, char32_t> and not std::same_as<T, wchar_t>> =
        nullptr>
constexpr T bit_floor(T x) noexcept {
  if (x != 0) {
    return T(1) << (std::numeric_limits<T>::digits - std::countl_zero(x) - 1);
  }
  return 0;
}

template <
    std::unsigned_integral T,
    Requires<not std::same_as<T, bool> and not std::same_as<T, char> and
             not std::same_as<T, char8_t> and not std::same_as<T, char16_t> and
             not std::same_as<T, char32_t> and not std::same_as<T, wchar_t>> =
        nullptr>
constexpr bool has_single_bit(T x) noexcept {
  return std::popcount(x) == 1;
}
}  // namespace std
