// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>

/// \cond
template <size_t Dim>
class Direction;
template <size_t Dim, typename Frame>
class ElementMap;
template <size_t Dim>
class Mesh;
template <typename T>
class Variables;
namespace domain {
template <typename SourceFrame, typename TargetFrame, size_t Dim>
class CoordinateMapBase;
namespace FunctionsOfTime {
class FunctionOfTime;
}  // namespace FunctionsOfTime
}  // namespace domain
namespace Frame {
struct Grid;
struct Inertial;
}  // namespace Frame
namespace evolution::dg::Tags {
struct MagnitudeOfNormal;
template <size_t Dim>
struct NormalCovector;
}  // namespace evolution::dg::Tags
   /// \endcond

/// Compute the inertial frame Euclidean normal covector and its magnitude on
/// the face of an element in a given direction
template <size_t Dim>
Variables<tmpl::list<::evolution::dg::Tags::MagnitudeOfNormal,
                     ::evolution::dg::Tags::NormalCovector<Dim>>>
normal_covector_and_magnitude(
    const Mesh<Dim - 1>& mesh, const Direction<Dim>& direction,
    const ElementMap<Dim, Frame::Grid>& element_map,
    const domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, Dim>&
        coordinate_map,
    const double t,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time);
