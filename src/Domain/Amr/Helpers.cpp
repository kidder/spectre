// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Domain/Amr/Helpers.hpp"

#include "DataStructures/Index.hpp"
#include "DataStructures/IndexIterator.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/Neighbors.hpp"
#include "Domain/Structure/OrientationMap.hpp"
#include "Domain/Structure/SegmentId.hpp"
#include "Parallel/Printf.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/StlStreamDeclarations.hpp"

namespace {
bool overlapping_within_one_level(const SegmentId& segment1,
                                  const SegmentId& segment2) {
  if (0 == segment1.refinement_level())
    return (2 > segment2.refinement_level());
  if (0 == segment2.refinement_level())
    return (2 > segment1.refinement_level());
  return (segment1 == segment2 || segment1.id_of_parent() == segment2 ||
          segment1 == segment2.id_of_parent());
}
}  // namespace

namespace amr {
template <size_t VolumeDim>
std::array<size_t, VolumeDim> desired_refinement_levels(
    const ElementId<VolumeDim>& id,
    const std::array<amr::Flag, VolumeDim>& flags) noexcept {
  std::array<size_t, VolumeDim> result{};

  for (size_t d = 0; d < VolumeDim; ++d) {
    ASSERT(amr::Flag::Undefined != gsl::at(flags, d),
           "Undefined amr::Flag in dimension " << d);
    gsl::at(result, d) = gsl::at(id.segment_ids(), d).refinement_level();
    if (amr::Flag::Join == gsl::at(flags, d)) {
      --gsl::at(result, d);
    } else if (amr::Flag::Split == gsl::at(flags, d)) {
      ++gsl::at(result, d);
    }
  }
  return result;
}

template <size_t VolumeDim>
std::array<size_t, VolumeDim> desired_refinement_levels_of_neighbor(
    const ElementId<VolumeDim>& neighbor_id,
    const std::array<amr::Flag, VolumeDim>& neighbor_flags,
    const OrientationMap<VolumeDim>& orientation) noexcept {
  if (orientation.is_aligned()) {
    return desired_refinement_levels(neighbor_id, neighbor_flags);
  }
  std::array<size_t, VolumeDim> result{};
  for (size_t d = 0; d < VolumeDim; ++d) {
    ASSERT(amr::Flag::Undefined != gsl::at(neighbor_flags, d),
           "Undefined amr::Flag in dimension " << d);
    const size_t mapped_dim = orientation(d);
    gsl::at(result, d) =
        gsl::at(neighbor_id.segment_ids(), mapped_dim).refinement_level();
    if (amr::Flag::Join == gsl::at(neighbor_flags, mapped_dim)) {
      --gsl::at(result, d);
    } else if (amr::Flag::Split == gsl::at(neighbor_flags, mapped_dim)) {
      ++gsl::at(result, d);
    }
  }
  return result;
}

template <size_t VolumeDim>
boost::rational<size_t> fraction_of_block_volume(
    const ElementId<VolumeDim>& element_id) noexcept {
  const auto& segment_ids = element_id.segment_ids();
  size_t sum_of_refinement_levels = 0;
  for (const auto& segment_id : segment_ids) {
    sum_of_refinement_levels += segment_id.refinement_level();
  }
  return {1, two_to_the(sum_of_refinement_levels)};
}

template <size_t VolumeDim>
bool has_potential_sibling(const ElementId<VolumeDim>& element_id,
                           const Direction<VolumeDim>& direction) noexcept {
  const auto& segment_id_to_check =
      gsl::at(element_id.segment_ids(), direction.dimension());
  if (0 == segment_id_to_check.refinement_level()) {
    return false;
  }
  return direction.side() == segment_id_to_check.side_of_sibling();
}

template <size_t VolumeDim>
bool neighbors_are_within_one_refinement_level(
    const Element<VolumeDim>& element) noexcept {
  for (const auto& kv : element.neighbors()) {
    const auto& neighbors_in_dir = kv.second;
    const OrientationMap<VolumeDim>& orientation =
        neighbors_in_dir.orientation();
    for (size_t d = 0; d < VolumeDim; ++d) {
      const size_t my_dim_in_neighbor = orientation(d);
      const size_t my_level = element.id().segment_ids()[d].refinement_level();
      for (const auto& neighbor_id : neighbors_in_dir.ids()) {
        const size_t neighbor_level =
            neighbor_id.segment_ids()[my_dim_in_neighbor].refinement_level();
        if (my_level > neighbor_level + 1 or neighbor_level > my_level + 1) {
          return false;
        }
      }
    }
  }
  return true;
}

template <size_t VolumeDim>
std::unordered_set<ElementId<VolumeDim>> new_neighbor_ids(
    const ElementId<VolumeDim>& my_id, const Direction<VolumeDim>& direction,
    const Neighbors<VolumeDim>& previous_neighbors,
    const std::unordered_map<ElementId<VolumeDim>,
                             std::array<amr::Flag, VolumeDim>>&
        previous_neighbor_amr_decisions) noexcept {
  std::unordered_set<ElementId<VolumeDim>> new_neighbors;
  const OrientationMap<VolumeDim>& orientation =
      previous_neighbors.orientation();
  const Direction<VolumeDim> neighbor_direction =
      orientation(direction.opposite());
  const size_t neighbor_dim = neighbor_direction.dimension();
  const std::array<SegmentId, VolumeDim> my_mapped_segment_ids =
      orientation(my_id.segment_ids());
  for (const auto& previous_neighbor_id : previous_neighbors.ids()) {
    std::array<std::vector<SegmentId>, VolumeDim> valid_segments;
    bool there_is_no_neighbor = false;
    for (size_t d = 0; d < VolumeDim; ++d) {
      const amr::Flag flag =
          previous_neighbor_amr_decisions.at(previous_neighbor_id)[d];
      const SegmentId& segment_id = previous_neighbor_id.segment_ids()[d];
      if (neighbor_dim == d) {
        valid_segments[d].push_back(
            amr::Flag::Join == flag
                ? segment_id.id_of_parent()
                : (amr::Flag::Split == flag
                       ? segment_id.id_of_child(neighbor_direction.side())
                       : segment_id));
      } else {
        if (amr::Flag::Join == flag) {
          if (overlapping_within_one_level(segment_id.id_of_parent(),
                                           my_mapped_segment_ids[d])) {
            valid_segments[d].push_back(segment_id.id_of_parent());
          }
        } else if (amr::Flag::Split == flag) {
          if (overlapping_within_one_level(segment_id.id_of_child(Side::Lower),
                                           my_mapped_segment_ids[d])) {
            valid_segments[d].push_back(segment_id.id_of_child(Side::Lower));
          }
          if (overlapping_within_one_level(segment_id.id_of_child(Side::Upper),
                                           my_mapped_segment_ids[d])) {
            valid_segments[d].push_back(segment_id.id_of_child(Side::Upper));
          }
        } else {
          if (overlapping_within_one_level(segment_id,
                                           my_mapped_segment_ids[d])) {
            valid_segments[d].push_back(segment_id);
          }
        }  // if-elseif-else on flag
        // There are either 0, 1, or 2 valid_segments
        if (0 == valid_segments[d].size()) {
          there_is_no_neighbor = true;
        }
      }  // if-else on neighbor_dim
      if (there_is_no_neighbor) {
        break;
      }
    }  // loop over d

    if (there_is_no_neighbor) {
      continue;
    }
    std::array<SegmentId, VolumeDim> new_neighbor_segment_ids;
    Index<VolumeDim> extents(0);
    for (size_t d = 0; d < VolumeDim; ++d) {
      extents[d] = valid_segments[d].size();
    }
    for (IndexIterator<VolumeDim> index(extents); index; ++index) {
      for (size_t d = 0; d < VolumeDim; ++d) {
        new_neighbor_segment_ids[d] = valid_segments[d][index()[d]];
      }
      new_neighbors.insert(ElementId<VolumeDim>(previous_neighbor_id.block_id(),
                                                new_neighbor_segment_ids));
    }
  }  // loop over previous_neighbors.ids()
  return new_neighbors;
}

template std::unordered_set<ElementId<2>> new_neighbor_ids(
    const ElementId<2>& my_id, const Direction<2>& direction,
    const Neighbors<2>& previous_neighbors,
    const std::unordered_map<ElementId<2>, std::array<amr::Flag, 2>>&
        previous_neighbor_amr_flags) noexcept;

template std::unordered_set<ElementId<3>> new_neighbor_ids(
    const ElementId<3>& my_id, const Direction<3>& direction,
    const Neighbors<3>& previous_neighbors,
    const std::unordered_map<ElementId<3>, std::array<amr::Flag, 3>>&
        previous_neighbor_amr_flags) noexcept;

template <>
std::unordered_set<ElementId<1>> new_neighbor_ids(
    const ElementId<1>& my_id, const Direction<1>& direction,
    const Neighbors<1>& previous_neighbors,
    const std::unordered_map<ElementId<1>, std::array<amr::Flag, 1>>&
        previous_neighbor_amr_flags) noexcept {
  std::unordered_set<ElementId<1>> ids_of_neighbors;
  const OrientationMap<1>& orientation = previous_neighbors.orientation();
  const Direction<1> neighbor_direction = orientation(direction.opposite());
  // Only 1 neighbor in 1-D in a given direction
  const ElementId<1>& previous_neighbor_id =
      *(previous_neighbors.ids().begin());
  const amr::Flag flag =
      previous_neighbor_amr_flags.at(previous_neighbor_id)[0];
  SegmentId previous_segment_id = previous_neighbor_id.segment_ids()[0];
  SegmentId new_segment_id =
      (amr::Flag::Join == flag
           ? previous_segment_id.id_of_parent()
           : (amr::Flag::Split == flag
                  ? previous_segment_id.id_of_child(neighbor_direction.side())
                  : previous_segment_id));
  ids_of_neighbors.emplace(previous_neighbor_id.block_id(),
                           std::array<SegmentId, 1>({{new_segment_id}}));
  Parallel::printf("Element %s now has neighbors %s in direction %s\n", my_id,
                   ids_of_neighbors, direction);
  return ids_of_neighbors;
}

template <size_t VolumeDim>
size_t dimension_to_join(
    const Element<VolumeDim>& element,
    const std::array<amr::Flag, VolumeDim>& flags,
    const std::unordered_map<ElementId<VolumeDim>,
                             std::array<amr::Flag, VolumeDim>>&
        neighbor_flags) noexcept {
  const auto number_of_join_flags = alg::count(flags, amr::Flag::Join);
  if (1 == number_of_join_flags) {
    return static_cast<size_t>(
        std::distance(std::begin(flags), alg::find(flags, amr::Flag::Join)));
  }
  const auto& element_id = element.id();
  for (size_t d = 0; d < VolumeDim; ++d) {
    const auto parent_id = element_id.id_of_parent(d);
    const auto sibling_id = (element_id == parent_id.id_of_child(d, Side::Lower)
                                 ? parent_id.id_of_child(d, Side::Upper)
                                 : parent_id.id_of_child(d, Side::Lower));
    if (neighbor_flags.find(sibling_id) != std::end(neighbor_flags)) {
      return d;
    }
  }
  Parallel::printf("Element id: %s\n", element_id);
  ERROR("Could not determine dimension.\n");
}

/// \cond

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                                   \
  template std::array<size_t, DIM(data)> desired_refinement_levels<DIM(data)>( \
      const ElementId<DIM(data)>&,                                             \
      const std::array<amr::Flag, DIM(data)>&) noexcept;                       \
  template std::array<size_t, DIM(data)>                                       \
  desired_refinement_levels_of_neighbor<DIM(data)>(                            \
      const ElementId<DIM(data)>&, const std::array<amr::Flag, DIM(data)>&,    \
      const OrientationMap<DIM(data)>&) noexcept;                              \
  template boost::rational<size_t> fraction_of_block_volume<DIM(data)>(        \
      const ElementId<DIM(data)>& element_id) noexcept;                        \
  template bool has_potential_sibling(                                         \
      const ElementId<DIM(data)>& element_id,                                  \
      const Direction<DIM(data)>& direction) noexcept;                         \
  template bool neighbors_are_within_one_refinement_level(                     \
      const Element<DIM(data)>& element) noexcept;                             \
  template size_t dimension_to_join(                                           \
      const Element<DIM(data)>&, const std::array<amr::Flag, DIM(data)>&,      \
      const std::unordered_map<ElementId<DIM(data)>,                           \
                               std::array<amr::Flag, DIM(data)>>&) noexcept;

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3))

#undef DIM
#undef INSTANTIATE
/// \endcond
}  // namespace amr
