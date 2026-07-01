// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include "Parallel/Info.hpp"
#include "Utilities/Serialization/Serialize.hpp"
#include "Utilities/System/Info.hpp"

namespace {
void test(const sys::Info& info) {
  CHECK(1 == info.number_of_procs());
  CHECK(0 == info.my_proc());
  CHECK(1 == info.number_of_nodes());
  CHECK(0 == info.my_node());
  CHECK(1 == info.procs_on_node(info.my_node()));
  CHECK(0 == info.my_local_rank());
  CHECK(0 == info.first_proc_on_node(info.my_node()));
  CHECK(0 == info.local_rank_of(info.my_proc()));
  CHECK(0 == info.node_of(info.my_proc()));
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Parallel.Info", "[Unit][Parallel]") {
  const Parallel::Info info{};
  test(info);
  const auto pupped_info = serialize_and_deserialize(info);
  test(pupped_info);
  const auto cloned_info = info.get_clone();
  test(*cloned_info);
}
