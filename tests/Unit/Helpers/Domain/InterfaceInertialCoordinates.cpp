// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Helpers/Domain/InterfaceInertialCoordinates.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>

#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/ElementMap.hpp"
#include "Domain/InterfaceLogicalCoordinates.hpp"

template <size_t Dim>
tnsr::I<DataVector, Dim, Frame::Inertial> interface_inertial_coordinates(
    const Mesh<Dim - 1>& mesh, const Direction<Dim>& direction,
    const ElementMap<Dim, Frame::Grid>& element_map,
    const domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, Dim>&
        coordinate_map,
    const double t,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) {
  return coordinate_map(
      element_map(interface_logical_coordinates(mesh, direction)), t,
      functions_of_time);
}

template tnsr::I<DataVector, 1, Frame::Inertial> interface_inertial_coordinates(
    const Mesh<0>& mesh, const Direction<1>& direction,
    const ElementMap<1, Frame::Grid>& element_map,
    const domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, 1>&
        coordinate_map,
    const double t,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time);

template tnsr::I<DataVector, 2, Frame::Inertial> interface_inertial_coordinates(
    const Mesh<1>& mesh, const Direction<2>& direction,
    const ElementMap<2, Frame::Grid>& element_map,
    const domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, 2>&
        coordinate_map,
    const double t,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time);

template tnsr::I<DataVector, 3, Frame::Inertial> interface_inertial_coordinates(
    const Mesh<2>& mesh, const Direction<3>& direction,
    const ElementMap<3, Frame::Grid>& element_map,
    const domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, 3>&
        coordinate_map,
    const double t,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time);
