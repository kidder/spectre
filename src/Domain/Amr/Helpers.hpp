// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Functions used for adaptive mesh refinement decisions.

#pragma once

#include <array>
#include <boost/rational.hpp>
#include <cstddef>
#include <unordered_map>
#include <unordered_set>

#include "Domain/Amr/Flag.hpp"

/// \cond
template <size_t VolumeDim>
class Direction;

template <size_t VolumeDim>
class Element;

template <size_t VolumeDim>
class ElementId;

template <size_t VolumeDim>
class OrientationMap;

template <size_t VolumeDim>
class Neighbors;
/// \endcond

namespace amr {
/// \ingroup ComputationalDomainGroup
/// \brief Computes the desired refinement level of the Element with ElementId
/// `id` given the desired amr::Flag%s `flags`
template <size_t VolumeDim>
std::array<size_t, VolumeDim> desired_refinement_levels(
    const ElementId<VolumeDim>& id,
    const std::array<amr::Flag, VolumeDim>& flags) noexcept;

/// \ingroup ComputationalDomainGroup
/// \brief Computes the desired refinement level of a neighboring Element with
/// ElementId `neighbor_id` given its desired amr::Flag%s `neighbor_flags`
/// taking into account the OrientationMap `orientation` of the neighbor
///
/// \details The OrientationMap `orientation` is that from the Element that has
/// a neighbor with ElementId `neighbor_id`
template <size_t VolumeDim>
std::array<size_t, VolumeDim> desired_refinement_levels_of_neighbor(
    const ElementId<VolumeDim>& neighbor_id,
    const std::array<amr::Flag, VolumeDim>& neighbor_flags,
    const OrientationMap<VolumeDim>& orientation) noexcept;

/// \ingroup ComputationalDomainGroup
/// Fraction of the logical volume of a block covered by an element
///
/// \note The sum of this over all the elements of a block should be one
template <size_t VolumeDim>
boost::rational<size_t> fraction_of_block_volume(
    const ElementId<VolumeDim>& element_id) noexcept;

/// \ingroup ComputationalDomainGroup
/// \brief Whether or not the Element with `element_id` can have a sibling
/// in the given `direction`
template <size_t VolumeDim>
bool has_potential_sibling(const ElementId<VolumeDim>& element_id,
                           const Direction<VolumeDim>& direction) noexcept;

/// \ingroup ComputationalDomainGroup
/// \brief returns true if all the neighbors of `element` are within one
/// refinement level.
template <size_t VolumeDim>
bool neighbors_are_within_one_refinement_level(
    const Element<VolumeDim>& element) noexcept;

/// \ingroup ComputationalDomainGroup
/// \brief returns the ElementId%s of the neighbors in the given `direction` of
/// the Element whose ElementId is `my_id` given the `previous_neighbors` and
/// their amr::Flag%s.
template <size_t VolumeDim>
std::unordered_set<ElementId<VolumeDim>> new_neighbor_ids(
    const ElementId<VolumeDim>& my_id, const Direction<VolumeDim>& direction,
    const Neighbors<VolumeDim>& previous_neighbors,
    const std::unordered_map<ElementId<VolumeDim>,
                             std::array<amr::Flag, VolumeDim>>&
        previous_neighbor_amr_flags) noexcept;

/// \ingroup ComputationalDomainGroup
/// \brief returns the dimension in which an Element will join
template <size_t VolumeDim>
size_t dimension_to_join(
    const Element<VolumeDim>& element,
    const std::array<amr::Flag, VolumeDim>& flags,
    const std::unordered_map<ElementId<VolumeDim>,
                             std::array<amr::Flag, VolumeDim>>&
        neighbor_flags) noexcept;
}  // namespace amr
