// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "Utilities/System/Abort.hpp"

// [[OutputRegex, I failed]]
[[noreturn]] SPECTRE_TEST_CASE("Unit.TestingFramework.Abort", "[Unit]") {
  OUTPUT_TEST();
  sys::abort("I failed");
}
