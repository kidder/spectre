// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <boost/rational.hpp>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Tags.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Reduction.hpp"
#include "Utilities/MakeString.hpp"
#include "Utilities/System/Abort.hpp"

namespace amr::Actions {

/// \brief Check the validity of the computational domain after adaptive mesh
/// refinement
///
/// Checks the following:
/// - That the fraction of Block volume (in the logical coordinate frame)
/// covered by all Element%s is equal to the number of Block%s in the Domain
struct CheckDomain {
  template <typename ParallelComponent, typename DbTagList,
            typename Metavariables, typename ArrayIndex>
  static void apply(db::DataBox<DbTagList>& box,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/,
                    const boost::rational<size_t>& volume) noexcept {
    if constexpr (tmpl::list_contains_v<
                      DbTagList,
                      Parallel::Tags::FromGlobalCache<
                          ::domain::Tags::Domain<Metavariables::volume_dim>>>) {
      const boost::rational<size_t> number_of_blocks{
          db::get<::domain::Tags::Domain<Metavariables::volume_dim>>(box)
              .blocks()
              .size()};
      if (number_of_blocks != volume) {
        sys::abort(MakeString{} << "Check Domain failed!  Expected volume "
                                << number_of_blocks << ", not " << volume
                                << "\n");
      }
    } else {
      ERROR("Could not find the tag "
            << pretty_type::get_name<
                   ::domain::Tags::Domain<Metavariables::volume_dim>>()
            << " in the DataBox.");
    }
  }
};
}  // namespace amr::Actions
