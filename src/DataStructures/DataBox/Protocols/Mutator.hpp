// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <type_traits>

#include "Utilities/ProtocolHelpers.hpp"

namespace db::protocols {

/// \brief A DataBox mutator
///
/// A class conforming to this protocol can be used as the template argument for
/// a call to db::mutate_apply(const gsl::not_null<DataBox<BoxTags>*> box,
/// Args&&... args).  The conforming class must provide the following:
///
/// - `return_tags`: A type list of tags corresponding to mutable items in the
///   DataBox passed to `db::mutate_apply` that may be modified.
/// - `argument_tags`: A type list of tags corresponding to items in the DataBox
///   passed to `db::mutate_apply` that may not be modified.
/// - `apply`: A static function whose return value is returned by
///   `db::mutate_apply`, and that takes as arguments:
///      - A `const gsl::not_null<Tag::type*>` for each `Tag` in `return_tags`
///      - A `const db::const_item_type<Tag, BoxTags>` for each `Tag` in
///        `argument_tags`
///      - The additional arguments passed to `db::mutate_apply`
struct Mutator {
  template <typename ConformingType>
  struct test {
    using argument_tags = typename ConformingType::argument_tags;
    using return_tags = typename ConformingType::return_tags;
    static_assert(
        std::is_same_v<void, std::void_t<decltype(&ConformingType::apply)>>);
  };
};
}  // namespace db::protocols
